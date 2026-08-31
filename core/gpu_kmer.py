"""GPU-accelerated canonical k-mer counting (Stage 1 of GPU acceleration).

Replicates the exact semantics of the Jellyfish path used by the probe
pipeline: for a set of query (probe) k-mers, count how many positions in a
reference (genome/transcript) have the same *canonical* k-mer, where the
canonical form is min(kmer, reverse_complement(kmer)) under the 2-bit
ordering A<C<G<T. This matches `jellyfish count -C` + `jellyfish query`.

Only exact canonical matching is done here (fully sensitive, no heuristics),
so results are bit-for-bit comparable to Jellyfish. k-mers containing any
non-ACGT base are skipped, exactly as Jellyfish does.

A pure-NumPy CPU fallback is provided so the pipeline keeps working when no
GPU / CuPy is available. Only k <= 32 is handled on the GPU (a k-mer fits in
one uint64); for k > 32 the caller should fall back to Jellyfish.
"""

import gzip
import os
import sys

import numpy as np

MAX_GPU_K = 32

_BASE_LUT = np.full(256, 255, dtype=np.uint8)
for _b, _code in ((b'A', 0), (b'C', 1), (b'G', 2), (b'T', 3),
                  (b'a', 0), (b'c', 1), (b'g', 2), (b't', 3)):
    _BASE_LUT[_b[0]] = _code


def _ensure_cuda_path():
    """CuPy JIT needs CUDA_PATH. The toolkit lives inside the conda env on the
    host, or at /usr/local/cuda in the CUDA Docker image."""
    if os.environ.get('CUDA_PATH'):
        return
    for cand in (os.path.join(sys.prefix, 'targets', 'x86_64-linux'),
                 '/usr/local/cuda'):
        if os.path.exists(os.path.join(cand, 'include', 'cuda_runtime.h')):
            os.environ['CUDA_PATH'] = cand
            return


_CP = None
_GPU_CHECKED = False


def _cupy():
    """Import CuPy lazily; return the module or None if unusable."""
    global _CP, _GPU_CHECKED
    if _GPU_CHECKED:
        return _CP
    _GPU_CHECKED = True
    try:
        _ensure_cuda_path()
        import cupy as cp
        if cp.cuda.runtime.getDeviceCount() < 1:
            _CP = None
        else:
            _CP = cp
    except Exception:
        _CP = None
    return _CP


def gpu_available():
    return _cupy() is not None


_COUNT_KERNEL_SRC = r'''
extern "C" __global__
void count_canonical(const unsigned char* codes, const long long n_codes,
                     const int k, const unsigned long long* sorted_keys,
                     const int n_keys, unsigned long long* counts,
                     const long long n_starts) {
    long long i = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n_starts) return;
    unsigned long long fwd = 0ULL, rc = 0ULL;
    const int shift = 2 * (k - 1);
    for (int j = 0; j < k; ++j) {
        unsigned char c = codes[i + j];
        if (c > 3) return;
        fwd = (fwd << 2) | (unsigned long long)c;
        rc = (rc >> 2) | ((unsigned long long)(3 - c) << shift);
    }
    unsigned long long canon = fwd < rc ? fwd : rc;
    int lo = 0, hi = n_keys - 1;
    while (lo <= hi) {
        int mid = (lo + hi) >> 1;
        unsigned long long v = sorted_keys[mid];
        if (v == canon) { atomicAdd(&counts[mid], 1ULL); return; }
        if (v < canon) lo = mid + 1; else hi = mid - 1;
    }
}
'''

_COUNT_KERNEL = None


def _get_kernel(cp):
    global _COUNT_KERNEL
    if _COUNT_KERNEL is None:
        _COUNT_KERNEL = cp.RawKernel(_COUNT_KERNEL_SRC, 'count_canonical')
    return _COUNT_KERNEL


def _encode_kmer(kmer):
    """2-bit encode a k-mer string to (code, valid). code is a Python int."""
    val = 0
    for ch in kmer:
        c = _BASE_LUT[ord(ch)]
        if c > 3:
            return 0, False
        val = (val << 2) | int(c)
    return val, True


def _canonical(code, k):
    """Canonical 2-bit code = min(code, reverse_complement(code))."""
    rc = 0
    x = code
    for _ in range(k):
        rc = (rc << 2) | (3 - (x & 3))
        x >>= 2
    return min(code, rc)


_NL = ord('\n')
_GT = ord('>')


def _encode_file(path):
    """Vectorized FASTA -> code array (bases 0-3, 255 = separator/invalid).

    Sequence lines are concatenated; header lines become runs of 255 so no
    counted k-mer window spans a header or two records. Pure NumPy, no per-line
    or per-record Python loop over sequence data — already fast, so the k-mer
    engine encodes directly rather than reading the on-disk cache.
    """
    if path.endswith('.gz'):
        with gzip.open(path, 'rb') as fh:
            data = np.frombuffer(fh.read(), dtype=np.uint8)
    else:
        data = np.fromfile(path, dtype=np.uint8)
    if data.size == 0:
        return np.empty(0, dtype=np.uint8)

    codes = _BASE_LUT[data]

    nl = np.flatnonzero(data == _NL)
    starts = np.empty(nl.size + 1, dtype=np.int64)
    starts[0] = 0
    starts[1:] = nl + 1
    starts = starts[starts < data.size]
    hdr_starts = starts[data[starts] == _GT]
    if hdr_starts.size:
        pos = np.searchsorted(nl, hdr_starts)
        hdr_ends = np.where(pos < nl.size, nl[pos], data.size - 1)
        for a, b in zip(hdr_starts.tolist(), hdr_ends.tolist()):
            codes[a:b + 1] = 255

    return codes[data != _NL]


def _chunk_array(codes, k, chunk_codes):
    """Yield (codes segment, n_starts) covering every k-mer start in `codes`.

    Segments carry a (k-1) tail so windows near a boundary have all k bases,
    and n_starts partitions the start positions so each is counted once.
    """
    n = codes.size
    if n < k:
        return
    if n <= chunk_codes + k - 1:
        yield codes, n - k + 1
        return
    off = 0
    while n - off >= k:
        if n - off <= chunk_codes + k - 1:
            yield codes[off:], (n - off) - k + 1
            return
        yield codes[off:off + chunk_codes + k - 1], chunk_codes
        off += chunk_codes


_SEP = np.array([255], dtype=np.uint8)


def _encode_batch(paths):
    """Encode several FASTA files into one array, 255-separated between files."""
    out = []
    for p in paths:
        codes = _encode_file(p)
        if codes.size:
            out.append(codes)
            out.append(_SEP)
    if not out:
        return np.empty(0, dtype=np.uint8)
    return np.concatenate(out)


def _batch_by_bytes(files, target_bytes):
    """Group files so each batch is ~target_bytes of raw FASTA (balances load)."""
    batches, cur, cur_bytes = [], [], 0
    for f in files:
        try:
            sz = os.path.getsize(f)
        except OSError:
            sz = 0
        cur.append(f)
        cur_bytes += sz
        if cur_bytes >= target_bytes:
            batches.append(cur)
            cur, cur_bytes = [], 0
    if cur:
        batches.append(cur)
    return batches


def _untrack_shm(name):
    """Detach a shared-memory segment from this process's resource_tracker.

    Worker processes create the segments but the GPU consumer owns their
    lifecycle (it unlinks them). Without this, each worker's resource_tracker
    reports the already-unlinked segments as leaks at shutdown (bpo-39959).
    """
    try:
        from multiprocessing import resource_tracker
        resource_tracker.unregister('/' + name, 'shared_memory')
    except Exception:
        pass


def _encode_to_shm(paths):
    """Worker task: encode a batch into shared memory; return (name, n_codes).

    Only the small name/size tuple crosses the process boundary, so the
    encoded gigabytes never go through pickle — the GPU consumer maps the
    shared block with zero copy.
    """
    from multiprocessing import shared_memory
    codes = _encode_batch(paths)
    if codes.size == 0:
        return None, 0
    shm = shared_memory.SharedMemory(create=True, size=int(codes.nbytes))
    buf = np.ndarray(codes.shape, dtype=np.uint8, buffer=shm.buf)
    buf[:] = codes
    name = shm.name
    _untrack_shm(name)
    shm.close()
    return name, int(codes.size)


def _shm_available(probe_bytes):
    """True if a shared-memory block of ~probe_bytes can be created (e.g.
    /dev/shm large enough). Docker defaults /dev/shm to 64 MB, which would
    break the zero-copy feed; in that case the caller uses the pickle path.
    """
    from multiprocessing import shared_memory
    try:
        shm = shared_memory.SharedMemory(create=True, size=max(1, probe_bytes))
        shm.close()
        shm.unlink()
        return True
    except Exception:
        return False


def _safe_shm_window(workers, batch_bytes):
    """Number of outstanding shared-memory batches that fit in /dev/shm.

    Each in-flight batch holds ~batch_bytes of codes; too many at once exhausts
    /dev/shm (small by default in Docker) and the worker processes get killed.
    Budget ~60% of the free space and cap at the desired 2*workers window.
    """
    want = max(2, (workers or 1) * 2)
    try:
        st = os.statvfs('/dev/shm')
        free = st.f_bavail * st.f_frsize
    except OSError:
        return want
    fit = int(free * 0.6) // max(1, batch_bytes)
    return max(2, min(want, fit))


def _iter_shm_arrays(files, workers, batch_bytes):
    """Yield code arrays backed by shared memory, encoded in parallel.

    Keeps a bounded window of outstanding tasks so transient shared memory
    stays limited, and unlinks each block after the consumer has used it.
    """
    from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
    from multiprocessing import shared_memory

    batches = _batch_by_bytes(files, batch_bytes)
    if not batches:
        return
    window = _safe_shm_window(workers, batch_bytes)
    with ProcessPoolExecutor(max_workers=workers) as ex:
        it = iter(batches)
        futures = set()
        for _ in range(window):
            try:
                futures.add(ex.submit(_encode_to_shm, next(it)))
            except StopIteration:
                break
        while futures:
            done, futures = wait(futures, return_when=FIRST_COMPLETED)
            for fut in done:
                try:
                    futures.add(ex.submit(_encode_to_shm, next(it)))
                except StopIteration:
                    pass
                name, size = fut.result()
                if not name or size == 0:
                    continue
                shm = shared_memory.SharedMemory(name=name)
                try:
                    yield np.ndarray((size,), dtype=np.uint8, buffer=shm.buf)
                finally:
                    shm.close()
                    shm.unlink()


def _iter_encoded(files, workers, batch=64):
    """Yield encoded code arrays (pickle IPC path; used by the CPU fallback)."""
    batches = _batch_by_bytes(files, batch * 2_000_000)
    if workers and workers > 1 and len(batches) > 1:
        import concurrent.futures
        with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as ex:
            for codes in ex.map(_encode_batch, batches):
                if codes.size:
                    yield codes
    else:
        for b in batches:
            codes = _encode_batch(b)
            if codes.size:
                yield codes


def build_canonical_keys(wanted_kmers, k):
    """Canonical 2-bit keys for the wanted k-mers, plus mappings.

    Returns (sorted_keys uint64, canon_to_originals, key_index). Shared with
    the fused aligner so it can count these k-mers on the GPU during the
    alignment scan.
    """
    canon_to_originals = {}
    for w in wanted_kmers:
        code, ok = _encode_kmer(w)
        if not ok:
            continue
        canon = _canonical(code, k)
        canon_to_originals.setdefault(canon, []).append(w)
    sorted_keys = np.array(sorted(canon_to_originals.keys()), dtype=np.uint64)
    key_index = {int(v): i for i, v in enumerate(sorted_keys)}
    return sorted_keys, canon_to_originals, key_index


def counts_array_to_dict(counts_np, canon_to_originals, key_index):
    """Map a per-key counts array back to {original_kmer: count}."""
    result = {}
    for canon, originals in canon_to_originals.items():
        c = int(counts_np[key_index[int(canon)]])
        for w in originals:
            result[w] = c
    return result


def count_kernel(cp):
    """The canonical-count CUDA kernel (shared with the fused aligner)."""
    return _get_kernel(cp)


def canonical_counts(files, k, wanted_kmers, device=None, chunk_mb=256,
                     workers=None):
    """Canonical-k-mer counts for `wanted_kmers` across FASTA `files`.

    Returns {wanted_kmer: count}, equivalent to querying a `jellyfish -C`
    database built from the same files. Uses the GPU when available and
    k <= 32; otherwise a NumPy CPU fallback with identical semantics.
    Files are encoded in parallel across `workers` processes (default: all
    CPUs) to keep the GPU fed.
    """
    wanted = [w for w in wanted_kmers]
    if not wanted:
        return {}
    if workers is None:
        workers = os.cpu_count() or 1

    canon_to_originals = {}
    for w in wanted:
        code, ok = _encode_kmer(w)
        if not ok:
            continue
        canon = _canonical(code, k)
        canon_to_originals.setdefault(canon, []).append(w)

    if not canon_to_originals:
        return {w: 0 for w in wanted}

    sorted_keys = np.array(sorted(canon_to_originals.keys()), dtype=np.uint64)
    key_index = {int(v): i for i, v in enumerate(sorted_keys)}

    cp = _cupy() if k <= MAX_GPU_K else None
    chunk_codes = max(1 << 20, chunk_mb * 1024 * 1024)

    if cp is not None:
        counts = _count_gpu(cp, files, k, sorted_keys, chunk_codes, device, workers)
    else:
        counts = _count_cpu(files, k, sorted_keys, key_index, chunk_codes, workers)

    result = {}
    for canon, originals in canon_to_originals.items():
        c = int(counts[key_index[int(canon)]])
        for w in originals:
            result[w] = c
    for w in wanted:
        result.setdefault(w, 0)
    return result


def _count_gpu(cp, files, k, sorted_keys, chunk_codes, device, workers,
               batch_bytes=64 * 1024 * 1024):
    from concurrent.futures.process import BrokenProcessPool
    dev = 0 if device is None else device
    with cp.cuda.Device(dev):
        keys_d = cp.asarray(sorted_keys)
        counts_d = cp.zeros(sorted_keys.shape[0], dtype=cp.uint64)
        kernel = _get_kernel(cp)
        threads = 256
        use_shm = workers and workers > 1 and _shm_available(batch_bytes)

        def _run(source):
            for codes in source:
                for seg, n_starts in _chunk_array(codes, k, chunk_codes):
                    if n_starts <= 0:
                        continue
                    codes_d = cp.asarray(seg)
                    blocks = (n_starts + threads - 1) // threads
                    kernel((blocks,), (threads,),
                           (codes_d, np.int64(codes_d.size), np.int32(k), keys_d,
                            np.int32(keys_d.size), counts_d, np.int64(n_starts)))

        if use_shm:
            try:
                _run(_iter_shm_arrays(files, workers, batch_bytes))
            except BrokenProcessPool:
                counts_d = cp.zeros(sorted_keys.shape[0], dtype=cp.uint64)
                _run(_iter_encoded(files, workers))
        else:
            _run(_iter_encoded(files, workers))
        cp.cuda.Stream.null.synchronize()
        return cp.asnumpy(counts_d)


def _count_cpu(files, k, sorted_keys, key_index, chunk_codes, workers):
    """Vectorized NumPy fallback with identical canonical semantics."""
    counts = np.zeros(sorted_keys.shape[0], dtype=np.uint64)
    key_set = sorted_keys
    segments = ((seg, ns) for codes in _iter_encoded(files, workers)
                for seg, ns in _chunk_array(codes, k, chunk_codes))
    for seg, n_starts in segments:
        if n_starts <= 0:
            continue
        codes = seg.astype(np.int64)
        fwd = np.zeros(n_starts, dtype=np.uint64)
        rc = np.zeros(n_starts, dtype=np.uint64)
        valid = np.ones(n_starts, dtype=bool)
        shift = 2 * (k - 1)
        for j in range(k):
            c = codes[j:j + n_starts]
            valid &= (c <= 3)
            cc = np.where(c <= 3, c, 0).astype(np.uint64)
            fwd = (fwd << np.uint64(2)) | cc
            rc = (rc >> np.uint64(2)) | ((np.uint64(3) - cc) << np.uint64(shift))
        canon = np.minimum(fwd, rc)
        canon = canon[valid]
        if canon.size == 0:
            continue
        idx = np.searchsorted(key_set, canon)
        idx = np.clip(idx, 0, key_set.size - 1)
        hit = key_set[idx] == canon
        np.add.at(counts, idx[hit], 1)
    return counts
