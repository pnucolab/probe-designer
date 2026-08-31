"""GPU-accelerated fully-sensitive short-probe aligner (Stage 2).

Replaces RazerS3 for the ungapped (Hamming) case: for each probe it reports
EVERY reference position, on both strands, where the probe matches with
<= max_mismatches substitutions. Full sensitivity is guaranteed by
pigeonhole seeding: a window with <= e mismatches leaves at least one of
e+1 disjoint seeds mismatch-free, so an exact seed lookup cannot miss it.

Non-ACGT reference bases (N) count as mismatches, matching RazerS3. The
gapped case (max_bulges > 0) is not handled here; the caller falls back to
RazerS3 for that.

Output is a single SAM file compatible with what the pipeline parses from
RazerS3: QNAME FLAG RNAME POS 255 {L}M * 0 0 SEQ QUAL NM:i:n MD:Z:md.
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import gpu_kmer  # noqa: E402  (reuses _cupy, shared-memory ingestion helpers)

MAX_SEED_LEN = 32
MAX_E = 8
MIN_SEED_LEN = 4
MAX_B = 4            # max indels the gapped kernel's band supports
MAX_GAP_LEN = 64     # max probe length for the gapped kernel's local DP arrays

_NL = ord('\n')
_GT = ord('>')

_ALIGN_LUT = np.full(256, 4, dtype=np.uint8)
for _b, _c in ((b'A', 0), (b'C', 1), (b'G', 2), (b'T', 3),
               (b'a', 0), (b'c', 1), (b'g', 2), (b't', 3)):
    _ALIGN_LUT[_b[0]] = _c
_SEPARATOR = 255


def available_devices():
    """Indices of all usable CUDA devices (empty if no GPU/CuPy)."""
    cp = gpu_kmer._cupy()
    if cp is None:
        return []
    try:
        return list(range(cp.cuda.runtime.getDeviceCount()))
    except Exception:
        return []


def _split_files_by_bytes(files, n):
    """Partition files into n shards balanced by total byte size (greedy)."""
    order = sorted(files, key=lambda f: -(os.path.getsize(f) if os.path.exists(f) else 0))
    shards = [[] for _ in range(n)]
    sizes = [0] * n
    for f in order:
        i = min(range(n), key=lambda j: sizes[j])
        shards[i].append(f)
        try:
            sizes[i] += os.path.getsize(f)
        except OSError:
            pass
    return shards


def _merge_sam_parts(parts, sam_out):
    """Concatenate per-device SAM shards into one file (single header)."""
    with open(sam_out, 'w', encoding='utf-8') as out:
        out.write('@HD\tVN:1.4\tSO:unsorted\n')
        for p in parts:
            if not os.path.exists(p):
                continue
            with open(p, 'r', encoding='utf-8') as f:
                for line in f:
                    if not line.startswith('@'):
                        out.write(line)
            try:
                os.remove(p)
            except OSError:
                pass


def gpu_align_multi(probes_file, chrom_files, sam_out, max_mismatches, max_bulges=0,
                    devices=None, workers=None, batch_mb=256, capacity=8_000_000,
                    per_probe_cap=100):
    """gpu_align across all GPUs: shard the reference files by size, align each
    shard on its own device concurrently, then concatenate the SAM shards.

    Contigs are never split across shards, so the merged SAM needs no cross-shard
    dedup. Returns True on success, False if any shard fails."""
    import threading
    devs = devices if devices is not None else available_devices()
    if len(devs) <= 1:
        return gpu_align(probes_file, chrom_files, sam_out, max_mismatches,
                         max_bulges=max_bulges, device=(devs[0] if devs else 0),
                         workers=workers, batch_mb=batch_mb, capacity=capacity,
                         per_probe_cap=per_probe_cap)
    n = len(devs)
    shards = _split_files_by_bytes(list(chrom_files), n)
    dev_workers = max(1, (workers or os.cpu_count() or 1) // n)
    parts = [f'{sam_out}.dev{i}' for i in range(n)]
    ok = [False] * n

    def _run(i):
        if not shards[i]:
            ok[i] = True
            open(parts[i], 'w').close()
            return
        ok[i] = gpu_align(probes_file, shards[i], parts[i], max_mismatches,
                          max_bulges=max_bulges, device=devs[i], workers=dev_workers,
                          batch_mb=batch_mb, capacity=capacity, per_probe_cap=per_probe_cap)

    threads = [threading.Thread(target=_run, args=(i,)) for i in range(n)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    if not all(ok):
        for p in parts:
            try:
                os.remove(p)
            except OSError:
                pass
        return False
    _merge_sam_parts(parts, sam_out)
    return True


def gpu_align_and_count_multi(probes_file, groups, sam_out, max_mismatches, kmer_length,
                              max_bulges=0, devices=None, workers=None, batch_mb=256,
                              capacity=8_000_000, per_probe_cap=None):
    """Fused align+count across all GPUs: shard each group's files by size, run
    the fused scan per device concurrently, then merge SAM shards and sum the
    per-group k-mer counts. Returns (True, {label: {kmer: count}}) or (False, None)."""
    import threading
    devs = devices if devices is not None else available_devices()
    if len(devs) <= 1:
        return gpu_align_and_count(probes_file, groups, sam_out, max_mismatches,
                                   kmer_length, max_bulges=max_bulges,
                                   device=(devs[0] if devs else 0), workers=workers,
                                   batch_mb=batch_mb, capacity=capacity,
                                   per_probe_cap=per_probe_cap)
    n = len(devs)
    groups = list(groups)
    per_dev_groups = [[] for _ in range(n)]
    for label, files in groups:
        shards = _split_files_by_bytes(list(files), n)
        for i in range(n):
            per_dev_groups[i].append((label, shards[i]))
    dev_workers = max(1, (workers or os.cpu_count() or 1) // n)
    parts = [f'{sam_out}.dev{i}' for i in range(n)]
    results = [None] * n

    def _run(i):
        results[i] = gpu_align_and_count(
            probes_file, per_dev_groups[i], parts[i], max_mismatches, kmer_length,
            max_bulges=max_bulges, device=devs[i], workers=dev_workers,
            batch_mb=batch_mb, capacity=capacity, per_probe_cap=per_probe_cap)

    threads = [threading.Thread(target=_run, args=(i,)) for i in range(n)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    if any(r is None or not r[0] for r in results):
        for p in parts:
            try:
                os.remove(p)
            except OSError:
                pass
        return False, None
    _merge_sam_parts(parts, sam_out)
    merged = {}
    for _ok, gc in results:
        for label, counts in (gc or {}).items():
            m = merged.setdefault(label, {})
            for kmer, c in counts.items():
                m[kmer] = m.get(kmer, 0) + c
    return True, merged


def gpu_available():
    return gpu_kmer.gpu_available()


_KERNEL_SRC = r'''
extern "C" __global__
void seed_align(const unsigned char* ref, const long long n_ref,
                const int seed_len, const int e,
                const unsigned long long* sorted_seeds, const int n_seeds,
                const int* post_start, const int* post_q, const int* post_soff,
                const unsigned char* qcodes, const int* q_off, const int* q_len,
                int* out_count, const int capacity,
                int* out_q, long long* out_pos, int* out_nmm,
                short* out_mmpos, signed char* out_mmref) {
    long long i = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (i + seed_len > n_ref) return;

    unsigned long long seed = 0ULL;
    for (int j = 0; j < seed_len; ++j) {
        unsigned char c = ref[i + j];
        if (c > 3) return;
        seed = (seed << 2) | (unsigned long long)c;
    }

    int lo = 0, hi = n_seeds - 1, found = -1;
    while (lo <= hi) {
        int mid = (lo + hi) >> 1;
        unsigned long long v = sorted_seeds[mid];
        if (v == seed) { found = mid; break; }
        if (v < seed) lo = mid + 1; else hi = mid - 1;
    }
    if (found < 0) return;

    for (int p = post_start[found]; p < post_start[found + 1]; ++p) {
        int q = post_q[p];
        int soff = post_soff[p];
        long long cand = i - soff;
        if (cand < 0) continue;
        int L = q_len[q];
        if (cand + L > n_ref) continue;
        int qoff = q_off[q];

        int mm = 0, boundary = 0;
        short mmpos[MAX_E_PLACEHOLDER];
        signed char mmref[MAX_E_PLACEHOLDER];
        for (int j = 0; j < L; ++j) {
            unsigned char rc = ref[cand + j];
            if (rc == 255) { boundary = 1; break; }
            unsigned char qc = qcodes[qoff + j];
            if (rc != qc) {
                if (mm < e) { mmpos[mm] = (short)j; mmref[mm] = (signed char)rc; }
                mm++;
                if (mm > e) break;
            }
        }
        if (boundary || mm > e) continue;

        int jcur = soff / seed_len;
        int jmin = -1;
        for (int j0 = 0; j0 <= e; ++j0) {
            int a = j0 * seed_len, b = a + seed_len;
            int clean = 1;
            for (int t = 0; t < mm; ++t) {
                if (mmpos[t] >= a && mmpos[t] < b) { clean = 0; break; }
            }
            if (clean) { jmin = j0; break; }
        }
        if (jmin != jcur) continue;

        int outi = atomicAdd(out_count, 1);
        if (outi < capacity) {
            out_q[outi] = q;
            out_pos[outi] = cand;
            out_nmm[outi] = mm;
            for (int t = 0; t < e; ++t) {
                if (t < mm) { out_mmpos[(long long)outi * e + t] = mmpos[t];
                              out_mmref[(long long)outi * e + t] = mmref[t]; }
                else { out_mmpos[(long long)outi * e + t] = -1; }
            }
        }
    }
}
'''

_KERNEL = None


def _get_kernel(cp, e):
    global _KERNEL
    if _KERNEL is None:
        src = _KERNEL_SRC.replace('MAX_E_PLACEHOLDER', str(MAX_E))
        _KERNEL = cp.RawKernel(src, 'seed_align')
    return _KERNEL


# Gapped (edit-distance) kernel: seed exactly, then a banded DP verifies every
# candidate diagonal, allowing up to `e` total edits of which up to `B` are
# single-base indels. It emits a compact edit script (substitutions + indels);
# the host turns that into a left-aligned, RazerS3-compatible CIGAR + MD.
#
# Pigeonhole still holds for edit distance: an alignment with <= e edits leaves
# at least one of the e+1 disjoint probe seeds free of substitutions AND indels,
# so an exact seed lookup cannot miss it. The seed anchors probe[soff] to
# ref[i]; the true probe start is i - soff + net, net in [-B, B] (net indels
# before the seed), so each seed hit probes a 2B+1 band of diagonals.
_GAPPED_KERNEL_SRC = r'''
#define MAXLEN LEN_PLACEHOLDER
#define BB B_PLACEHOLDER
#define BW (2*BB+1)
#define IDX(a,d) ((a)*BW + ((d)+BB))
#define INF 1000000

extern "C" __global__
void seed_align_gapped(const unsigned char* ref, const long long n_ref,
        const int seed_len, const int e,
        const unsigned long long* sorted_seeds, const int n_seeds,
        const int* post_start, const int* post_q, const int* post_soff,
        const unsigned char* qcodes, const int* q_off, const int* q_len,
        int* out_count, const int capacity,
        int* out_q, long long* out_pos, int* out_edits,
        short* out_sub_pos, signed char* out_sub_ref, int* out_nsub,
        short* out_ind_pos, signed char* out_ind_type, signed char* out_ind_ref,
        int* out_nind) {
    long long gi = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (gi + seed_len > n_ref) return;

    unsigned long long seed = 0ULL;
    for (int j = 0; j < seed_len; ++j) {
        unsigned char c = ref[gi + j];
        if (c > 3) return;
        seed = (seed << 2) | (unsigned long long)c;
    }
    int lo = 0, hi = n_seeds - 1, found = -1;
    while (lo <= hi) {
        int mid = (lo + hi) >> 1;
        unsigned long long v = sorted_seeds[mid];
        if (v == seed) { found = mid; break; }
        if (v < seed) lo = mid + 1; else hi = mid - 1;
    }
    if (found < 0) return;

    int dp[(MAXLEN + 1) * BW];
    signed char bt[(MAXLEN + 1) * BW];

    for (int p = post_start[found]; p < post_start[found + 1]; ++p) {
        int q = post_q[p];
        int soff = post_soff[p];
        int L = q_len[q];
        if (L > MAXLEN) continue;
        int qoff = q_off[q];
        int blk = soff / seed_len;

        for (int net = -BB; net <= BB; ++net) {
            long long cand = gi - (long long)soff + net;
            if (cand < 0) continue;
            if (cand + L + BB > n_ref) continue;

            // Banded DP: dp[a][d] = min edits aligning probe[0..a) to
            // ref[cand .. cand+a+d), with |a-t| <= B. Start-anchored: only d=0
            // has 0 cost at a=0 (no leading indels).
            for (int d = -BB; d <= BB; ++d) dp[IDX(0, d)] = (d == 0) ? 0 : INF;
            for (int a = 1; a <= L; ++a) {
                unsigned char pc = qcodes[qoff + a - 1];
                for (int d = -BB; d <= BB; ++d) {
                    int t = a + d;
                    int best = INF; signed char bb = 3;
                    if (t >= 0) {
                        long long rj = cand + t - 1;   // ref index for diag/del
                        unsigned char rc = (t - 1 >= 0) ? ref[rj] : 255;
                        // diagonal (match / substitution)
                        if (t - 1 >= 0 && rc != 255 && dp[IDX(a - 1, d)] < INF) {
                            int cost = (rc == pc && rc <= 3) ? 0 : 1;
                            int v = dp[IDX(a - 1, d)] + cost;
                            if (v < best) { best = v; bb = 0; }
                        }
                        // deletion of ref[rj] (ref-only), from dp[a][d-1]
                        if (d - 1 >= -BB && t - 1 >= 0 && rc != 255 &&
                            dp[IDX(a, d - 1)] < INF) {
                            int v = dp[IDX(a, d - 1)] + 1;
                            if (v < best) { best = v; bb = 1; }
                        }
                    }
                    // insertion of probe[a-1] (read-only), from dp[a-1][d+1]
                    if (d + 1 <= BB && dp[IDX(a - 1, d + 1)] < INF) {
                        int v = dp[IDX(a - 1, d + 1)] + 1;
                        if (v < best) { best = v; bb = 2; }
                    }
                    dp[IDX(a, d)] = best; bt[IDX(a, d)] = bb;
                }
            }

            // pick end-anchored (last op diagonal) min-edit d
            int bestd = 99, beste = INF;
            for (int d = -BB; d <= BB; ++d) {
                int val = dp[IDX(L, d)];
                if (val <= e && bt[IDX(L, d)] == 0) {
                    if (val < beste || (val == beste && (d < 0 ? -d : d) <
                                        (bestd < 0 ? -bestd : bestd))) {
                        beste = val; bestd = d;
                    }
                }
            }
            if (bestd == 99) continue;

            // traceback: collect subs + indels, count indels, find blocks hit
            int nsub = 0, nind = 0;
            short sub_pos[MAX_E_PLACEHOLDER]; signed char sub_ref[MAX_E_PLACEHOLDER];
            short ind_pos[MAX_E_PLACEHOLDER]; signed char ind_type[MAX_E_PLACEHOLDER];
            signed char ind_ref[MAX_E_PLACEHOLDER];
            int a = L, d = bestd;
            int edited_block = 0;   // bitmask of blocks touched by an edit
            while (a > 0) {
                signed char op = bt[IDX(a, d)];
                int t = a + d;
                if (op == 0) {
                    long long rj = cand + t - 1;
                    unsigned char rc = ref[rj];
                    unsigned char pc = qcodes[qoff + a - 1];
                    if (!(rc == pc && rc <= 3)) {
                        if (nsub < e) { sub_pos[nsub] = (short)(a - 1);
                                        sub_ref[nsub] = (signed char)rc; }
                        nsub++;
                        edited_block |= (1 << ((a - 1) / seed_len));
                    }
                    a -= 1; /* d unchanged */
                } else if (op == 1) {
                    long long rj = cand + t - 1;
                    if (nind < e) { ind_pos[nind] = (short)a;
                                    ind_type[nind] = 0; /* del */
                                    ind_ref[nind] = (signed char)ref[rj]; }
                    nind++;
                    edited_block |= (1 << ((a < L ? a : L - 1) / seed_len));
                    d -= 1; /* a unchanged */
                } else if (op == 2) {
                    if (nind < e) { ind_pos[nind] = (short)(a - 1);
                                    ind_type[nind] = 1; /* ins */
                                    ind_ref[nind] = -1; }
                    nind++;
                    edited_block |= (1 << ((a - 1) / seed_len));
                    a -= 1; d += 1;
                } else break;
            }
            if (nind > BB || nsub + nind != beste) continue;

            // canonical emitter: only the first edit-free block reports the hit
            int first_free = -1;
            for (int j = 0; j <= e; ++j) {
                if (!(edited_block & (1 << j))) { first_free = j; break; }
            }
            if (first_free != blk) continue;

            int outi = atomicAdd(out_count, 1);
            if (outi < capacity) {
                out_q[outi] = q;
                out_pos[outi] = cand;
                out_edits[outi] = beste;
                out_nsub[outi] = nsub;
                out_nind[outi] = nind;
                for (int t2 = 0; t2 < e; ++t2) {
                    if (t2 < nsub) { out_sub_pos[(long long)outi * e + t2] = sub_pos[t2];
                                     out_sub_ref[(long long)outi * e + t2] = sub_ref[t2]; }
                    else out_sub_pos[(long long)outi * e + t2] = -1;
                    if (t2 < nind) { out_ind_pos[(long long)outi * e + t2] = ind_pos[t2];
                                     out_ind_type[(long long)outi * e + t2] = ind_type[t2];
                                     out_ind_ref[(long long)outi * e + t2] = ind_ref[t2]; }
                    else out_ind_pos[(long long)outi * e + t2] = -1;
                }
            }
        }
    }
}
'''

_GAPPED_KERNEL = None
_GAPPED_KEY = None


def _get_gapped_kernel(cp, max_len, B):
    global _GAPPED_KERNEL, _GAPPED_KEY
    key = (max_len, B)
    if _GAPPED_KERNEL is None or _GAPPED_KEY != key:
        src = (_GAPPED_KERNEL_SRC
               .replace('LEN_PLACEHOLDER', str(max_len))
               .replace('B_PLACEHOLDER', str(B))
               .replace('MAX_E_PLACEHOLDER', str(MAX_E)))
        _GAPPED_KERNEL = cp.RawKernel(src, 'seed_align_gapped')
        _GAPPED_KEY = key
    return _GAPPED_KERNEL


def _revcomp_codes(codes):
    return (3 - codes[::-1]).astype(np.uint8)


def _encode_seq(seq):
    arr = np.frombuffer(seq.encode('ascii', 'ignore'), dtype=np.uint8)
    return _ALIGN_LUT[arr]


def _build_probe_index(probes, e, seed_len):
    """Build query codes + seed postings for forward and RC of every probe.

    Query 2*p = forward probe p ('+' strand), 2*p+1 = reverse complement
    ('-' strand). Returns arrays ready to upload to the GPU.
    """
    q_codes_parts, q_off, q_len, q_probe, q_strand = [], [], [], [], []
    off = 0
    seed_map = {}
    for pi, (_pid, seq) in enumerate(probes):
        fwd = _encode_seq(seq)
        rc = _revcomp_codes(fwd)
        for strand, codes in ((0, fwd), (1, rc)):
            q = len(q_len)
            q_codes_parts.append(codes)
            q_off.append(off)
            q_len.append(codes.size)
            q_probe.append(pi)
            q_strand.append(strand)
            off += codes.size
            for j in range(e + 1):
                a = j * seed_len
                seg = codes[a:a + seed_len]
                if seg.size < seed_len or (seg > 3).any():
                    continue
                code = 0
                for c in seg.tolist():
                    code = (code << 2) | int(c)
                seed_map.setdefault(code, []).append((q, a))

    qcodes = (np.concatenate(q_codes_parts) if q_codes_parts
              else np.empty(0, dtype=np.uint8))
    sorted_seeds = np.array(sorted(seed_map.keys()), dtype=np.uint64)
    post_start = np.zeros(sorted_seeds.size + 1, dtype=np.int32)
    post_q, post_soff = [], []
    for idx, code in enumerate(sorted_seeds.tolist()):
        postings = seed_map[code]
        post_start[idx + 1] = post_start[idx] + len(postings)
        for (q, a) in postings:
            post_q.append(q)
            post_soff.append(a)
    return {
        'qcodes': qcodes,
        'q_off': np.array(q_off, dtype=np.int32),
        'q_len': np.array(q_len, dtype=np.int32),
        'q_probe': np.array(q_probe, dtype=np.int32),
        'q_strand': np.array(q_strand, dtype=np.int32),
        'sorted_seeds': sorted_seeds,
        'post_start': post_start,
        'post_q': np.array(post_q, dtype=np.int32),
        'post_soff': np.array(post_soff, dtype=np.int32),
    }


_REF_EMPTY = (np.empty(0, dtype=np.uint8), [],
              np.empty(0, np.int64), np.empty(0, np.int64))


def _encode_file_records(path):
    """FASTA -> (codes, names, starts, lens).

    codes: uint8, bases 0-3, N/other 4, 255 header runs separating contigs.
    Contig spans are derived vectorially (kept-bytes-before = position minus
    newlines-before via searchsorted) rather than a per-contig cumulative
    count, which keeps encoding roughly twice as fast as a naive loop.
    """
    if path.endswith('.gz'):
        import gzip
        with gzip.open(path, 'rb') as fh:
            data = np.frombuffer(fh.read(), dtype=np.uint8)
    else:
        data = np.fromfile(path, dtype=np.uint8)
    if data.size == 0:
        return _REF_EMPTY

    keep = data != _NL
    codes_full = _ALIGN_LUT[data]

    nl = np.flatnonzero(~keep)
    line_starts = np.empty(nl.size + 1, dtype=np.int64)
    line_starts[0] = 0
    line_starts[1:] = nl + 1
    line_starts = line_starts[line_starts < data.size]
    hdr = line_starts[data[line_starts] == _GT]
    if hdr.size == 0:
        return _REF_EMPTY

    nl_after = nl[np.searchsorted(nl, hdr)]
    for a, b in zip(hdr.tolist(), nl_after.tolist()):
        codes_full[a:b + 1] = 255
    codes = codes_full[keep]

    seq_start = nl_after + 1
    seq_end = np.empty(hdr.size, dtype=np.int64)
    seq_end[:-1] = hdr[1:]
    seq_end[-1] = data.size
    starts = seq_start - np.searchsorted(nl, seq_start)
    ends = seq_end - np.searchsorted(nl, seq_end)
    lens = ends - starts

    valid = lens > 0
    if not valid.any():
        return _REF_EMPTY
    starts = starts[valid].astype(np.int64)
    lens = lens[valid].astype(np.int64)
    hdr_v = hdr[valid]
    nla_v = nl_after[valid]
    names = []
    for h, he in zip(hdr_v.tolist(), nla_v.tolist()):
        tok = data[h + 1:he].tobytes().split(None, 1)
        names.append(tok[0].decode('ascii', 'ignore') if tok else '')
    return codes, names, starts, lens


def _encode_batch_records(paths):
    all_codes, all_names, all_starts, all_lens = [], [], [], []
    off = 0
    for p in paths:
        codes, names, starts, lens = _encode_file_records(p)
        if codes.size == 0:
            continue
        all_codes.append(codes)
        all_names.extend(names)
        all_starts.append(starts + off)
        all_lens.append(lens)
        off += codes.size
    if not all_codes:
        return _REF_EMPTY
    return (np.concatenate(all_codes), all_names,
            np.concatenate(all_starts), np.concatenate(all_lens))


def _batch_task(paths):
    from multiprocessing import shared_memory
    codes, names, starts, lens = _encode_batch_records(paths)
    names_blob = '\n'.join(names)
    if codes.size == 0:
        return None, 0, names_blob, starts, lens
    shm = shared_memory.SharedMemory(create=True, size=int(codes.nbytes))
    np.ndarray(codes.shape, dtype=np.uint8, buffer=shm.buf)[:] = codes
    name = shm.name
    gpu_kmer._untrack_shm(name)
    shm.close()
    return name, int(codes.size), names_blob, starts, lens


def _iter_batches(files, workers, batch_bytes):
    """Yield (codes, names, starts, lens) batches; encode in parallel, ship
    codes via shared memory and metadata via pickle."""
    batches = gpu_kmer._batch_by_bytes(files, batch_bytes)
    if not batches:
        return
    if not workers or workers <= 1 or len(batches) == 1 or not gpu_kmer._shm_available(batch_bytes):
        for b in batches:
            res = _encode_batch_records(b)
            if res[0].size:
                yield res
        return

    from concurrent.futures import (ProcessPoolExecutor, wait, FIRST_COMPLETED,
                                    process as _fut_process)
    from multiprocessing import shared_memory
    # Cap outstanding shared-memory batches so their combined size stays within
    # /dev/shm (small by default in Docker); otherwise the worker processes are
    # killed writing their blocks and the pool breaks.
    window = gpu_kmer._safe_shm_window(workers, batch_bytes)

    def _free(shm):
        if shm is not None:
            shm.close()
            shm.unlink()

    try:
        with ProcessPoolExecutor(max_workers=workers) as ex:
            it = iter(batches)
            futures = set()
            for _ in range(window):
                try:
                    futures.add(ex.submit(_batch_task, next(it)))
                except StopIteration:
                    break
            # Hand the consumer a view straight into the shared-memory block so
            # the GPU upload reads from it directly (no host-side copy). The block
            # is only freed once the consumer comes back for the next batch, by
            # which point it has finished uploading the previous one.
            prev = None
            while futures:
                done, futures = wait(futures, return_when=FIRST_COMPLETED)
                for fut in done:
                    try:
                        futures.add(ex.submit(_batch_task, next(it)))
                    except StopIteration:
                        pass
                    name, size, names_blob, starts, lens = fut.result()
                    if not name or size == 0:
                        continue
                    shm = shared_memory.SharedMemory(name=name)
                    codes = np.ndarray((size,), dtype=np.uint8, buffer=shm.buf)
                    names = names_blob.split('\n') if names_blob else []
                    _free(prev)
                    prev = shm
                    yield codes, names, starts, lens
            _free(prev)
    except _fut_process.BrokenProcessPool:
        # Shared-memory ingestion failed (e.g. /dev/shm exhausted); fall back to
        # single-process encoding so the scan still completes.
        for b in batches:
            res = _encode_batch_records(b)
            if res[0].size:
                yield res


def _md_and_seq(query_seq, strand, mmpos, mmref):
    """Build MD tag and SEQ (reference-oriented) for one hit."""
    L = len(query_seq)
    order = sorted(range(len(mmpos)), key=lambda t: mmpos[t])
    md, last = [], 0
    letters = 'ACGTN'
    for t in order:
        p = mmpos[t]
        md.append(str(p - last))
        md.append(letters[mmref[t]] if 0 <= mmref[t] < 5 else 'N')
        last = p + 1
    md.append(str(L - last))
    return ''.join(md), query_seq


def _probe_query_strings(probes):
    from Bio.Seq import Seq
    out = []
    for _pid, seq in probes:
        out.append((seq.upper(), str(Seq(seq.upper()).reverse_complement())))
    return out


def _align_batch(cp, out, ref_d, n_ref, seed_len, e, d, n_seeds, kernel,
                 threads, capacity, q_probe, q_strand, probes, qstrings,
                 names, starts, lens, seen, probe_hit_counts, per_probe_cap):
    """Run the seed-align kernel on one device batch and write its hits.

    Returns True on success, False on hit-buffer overflow (caller aborts)."""
    out_count = cp.zeros(1, dtype=cp.int32)
    out_q = cp.empty(capacity, dtype=cp.int32)
    out_pos = cp.empty(capacity, dtype=cp.int64)
    out_nmm = cp.empty(capacity, dtype=cp.int32)
    out_mmpos = cp.empty(capacity * max(1, e), dtype=cp.int16)
    out_mmref = cp.empty(capacity * max(1, e), dtype=cp.int8)

    n_pos = n_ref - seed_len + 1
    blocks = (n_pos + threads - 1) // threads
    kernel((blocks,), (threads,), (
        ref_d, np.int64(n_ref), np.int32(seed_len), np.int32(e),
        d['sorted_seeds'], np.int32(n_seeds),
        d['post_start'], d['post_q'], d['post_soff'],
        d['qcodes'], d['q_off'], d['q_len'],
        out_count, np.int32(capacity),
        out_q, out_pos, out_nmm, out_mmpos, out_mmref))
    cp.cuda.Stream.null.synchronize()

    nhits = int(out_count.get()[0])
    if nhits > capacity:
        return False
    if nhits == 0:
        return True

    h_q = cp.asnumpy(out_q[:nhits])
    h_pos = cp.asnumpy(out_pos[:nhits])
    h_nmm = cp.asnumpy(out_nmm[:nhits])
    h_mmpos = cp.asnumpy(out_mmpos[:nhits * max(1, e)]).reshape(nhits, max(1, e))
    h_mmref = cp.asnumpy(out_mmref[:nhits * max(1, e)]).reshape(nhits, max(1, e))
    _write_hits(out, h_q, h_pos, h_nmm, h_mmpos, h_mmref, e,
                q_probe, q_strand, probes, qstrings,
                names, starts, lens, seen, probe_hit_counts, per_probe_cap)
    return True


def _align_batch_gapped(cp, out, ref_d, n_ref, seed_len, e, d, n_seeds, kernel,
                        threads, capacity, q_probe, q_strand, probes, qstrings,
                        names, starts, lens, seen, probe_hit_counts,
                        per_probe_cap):
    """Gapped variant of _align_batch: runs the banded-DP kernel and writes
    hits whose CIGAR/MD carry left-aligned indels."""
    ew = max(1, e)
    out_count = cp.zeros(1, dtype=cp.int32)
    out_q = cp.empty(capacity, dtype=cp.int32)
    out_pos = cp.empty(capacity, dtype=cp.int64)
    out_edits = cp.empty(capacity, dtype=cp.int32)
    out_nsub = cp.empty(capacity, dtype=cp.int32)
    out_nind = cp.empty(capacity, dtype=cp.int32)
    out_sub_pos = cp.empty(capacity * ew, dtype=cp.int16)
    out_sub_ref = cp.empty(capacity * ew, dtype=cp.int8)
    out_ind_pos = cp.empty(capacity * ew, dtype=cp.int16)
    out_ind_type = cp.empty(capacity * ew, dtype=cp.int8)
    out_ind_ref = cp.empty(capacity * ew, dtype=cp.int8)

    n_pos = n_ref - seed_len + 1
    blocks = (n_pos + threads - 1) // threads
    kernel((blocks,), (threads,), (
        ref_d, np.int64(n_ref), np.int32(seed_len), np.int32(e),
        d['sorted_seeds'], np.int32(n_seeds),
        d['post_start'], d['post_q'], d['post_soff'],
        d['qcodes'], d['q_off'], d['q_len'],
        out_count, np.int32(capacity),
        out_q, out_pos, out_edits,
        out_sub_pos, out_sub_ref, out_nsub,
        out_ind_pos, out_ind_type, out_ind_ref, out_nind))
    cp.cuda.Stream.null.synchronize()

    nhits = int(out_count.get()[0])
    if nhits > capacity:
        return False
    if nhits == 0:
        return True

    h_q = cp.asnumpy(out_q[:nhits])
    h_pos = cp.asnumpy(out_pos[:nhits])
    h_edits = cp.asnumpy(out_edits[:nhits])
    h_nsub = cp.asnumpy(out_nsub[:nhits])
    h_nind = cp.asnumpy(out_nind[:nhits])
    h_sub_pos = cp.asnumpy(out_sub_pos[:nhits * ew]).reshape(nhits, ew)
    h_sub_ref = cp.asnumpy(out_sub_ref[:nhits * ew]).reshape(nhits, ew)
    h_ind_pos = cp.asnumpy(out_ind_pos[:nhits * ew]).reshape(nhits, ew)
    h_ind_type = cp.asnumpy(out_ind_type[:nhits * ew]).reshape(nhits, ew)
    h_ind_ref = cp.asnumpy(out_ind_ref[:nhits * ew]).reshape(nhits, ew)
    _write_hits_gapped(out, h_q, h_pos, h_edits, h_nsub, h_nind,
                       h_sub_pos, h_sub_ref, h_ind_pos, h_ind_type, h_ind_ref,
                       q_probe, q_strand, probes, qstrings,
                       names, starts, lens, seen, probe_hit_counts,
                       per_probe_cap)
    return True


_BASES = 'ACGTN'


def _left_align_columns(cols):
    """Left-align single-base indels in a column list of (op, q_base, r_base).

    op is 'M' (q vs r), 'I' (q base, r '-') or 'D' (r base, q '-'). Shifting a
    gap left through equal bases keeps the alignment and edit count identical,
    matching RazerS3's leftmost placement (e.g. a deletion inside a homopolymer
    reported at its first position)."""
    cols = [list(c) for c in cols]
    changed = True
    while changed:
        changed = False
        for i in range(1, len(cols)):
            op, q, r = cols[i]
            pop, pq, pr = cols[i - 1]
            if pop != 'M' or pq != pr:
                continue
            if op == 'D' and q == '-' and pq == r:      # deletion: swap q only
                cols[i][1] = pq
                cols[i - 1][1] = '-'
                cols[i - 1][0] = 'D'
                cols[i][0] = 'M'
                changed = True
            elif op == 'I' and r == '-' and pr == q:    # insertion: swap r only
                cols[i][2] = pr
                cols[i - 1][2] = '-'
                cols[i - 1][0] = 'I'
                cols[i][0] = 'M'
                changed = True
    return cols


def _cigar_md_from_columns(cols):
    """Build (CIGAR, MD) from an aligned column list (RazerS3 conventions)."""
    ops = []
    for op, _q, _r in cols:
        if ops and ops[-1][1] == op:
            ops[-1][0] += 1
        else:
            ops.append([1, op])
    cigar = ''.join(f'{n}{op}' for n, op in ops)

    md = []
    run = 0
    i = 0
    n = len(cols)
    while i < n:
        op, q, r = cols[i]
        if op == 'I':
            i += 1
            continue
        if op == 'M':
            if q == r:
                run += 1
            else:
                md.append(str(run)); md.append(r); run = 0
            i += 1
        else:  # 'D' — group consecutive deletions under one caret
            md.append(str(run)); md.append('^')
            while i < n and cols[i][0] == 'D':
                md.append(cols[i][2]); i += 1
            run = 0
    if run > 0 or not md:                 # RazerS3 omits a trailing zero run
        md.append(str(run))
    return cigar, ''.join(md)


def _write_hits_gapped(out, h_q, h_pos, h_edits, h_nsub, h_nind,
                       h_sub_pos, h_sub_ref, h_ind_pos, h_ind_type, h_ind_ref,
                       q_probe, q_strand, probes, qstrings,
                       names, contig_starts, contig_lens,
                       seen, probe_hit_counts, per_probe_cap):
    # Collect this batch's candidate hits, then apply RazerS3's duplicate
    # removal (drop a hit whose reference span overlaps a strictly-better hit of
    # the same probe/strand/contig). All hits to a contig live in one batch, so
    # overlap suppression is complete here.
    records = []
    for r in range(h_q.shape[0]):
        q = int(h_q[r])
        pi = int(q_probe[q])
        strand = int(q_strand[q])
        gpos = int(h_pos[r])
        edits = int(h_edits[r])
        nsub = int(h_nsub[r])
        nind = int(h_nind[r])

        query_seq = qstrings[pi][0] if strand == 0 else qstrings[pi][1]
        L = len(query_seq)

        subs = {int(h_sub_pos[r][t]): _BASES[int(h_sub_ref[r][t])]
                for t in range(nsub) if 0 <= int(h_sub_ref[r][t]) < 5}
        dels = {}
        inss = set()
        for t in range(nind):
            p = int(h_ind_pos[r][t])
            if int(h_ind_type[r][t]) == 0:
                dels.setdefault(p, []).append(_BASES[int(h_ind_ref[r][t])]
                                              if 0 <= int(h_ind_ref[r][t]) < 5
                                              else 'N')
            else:
                inss.add(p)

        cols = []
        for p in range(L + 1):
            if p in dels:
                for rb in dels[p]:
                    cols.append(('D', '-', rb))
            if p == L:
                break
            qb = query_seq[p]
            if p in inss:
                cols.append(('I', qb, '-'))
            elif p in subs:
                cols.append(('M', qb, subs[p]))
            else:
                cols.append(('M', qb, qb))

        cols = _left_align_columns(cols)
        if cols[0][0] != 'M' or cols[-1][0] != 'M':
            continue                       # leading/trailing indel: wrong locus
        cigar, md = _cigar_md_from_columns(cols)

        ci = int(np.searchsorted(contig_starts, gpos, side='right') - 1)
        if ci < 0 or ci >= len(names):
            continue
        cstart = int(contig_starts[ci])
        clen = int(contig_lens[ci])
        local = gpos - cstart
        ref_consumed = sum(1 for op, _q, _r in cols if op in ('M', 'D'))
        if local < 0 or local + ref_consumed > clen:
            continue

        records.append({
            'pi': pi, 'probe_id': probes[pi][0], 'strand': strand,
            'rname': names[ci], 'pos1': local + 1,
            'rstart': local, 'rend': local + ref_consumed,
            'edits': edits, 'cigar': cigar, 'md': md, 'seq': query_seq, 'L': L,
        })

    # Suppress hits dominated by a strictly-better overlapping hit (same
    # probe/strand/contig); process fewest-errors first so keepers win.
    kept = []
    by_group = {}
    records.sort(key=lambda h: (h['edits'], h['pos1']))
    for h in records:
        g = (h['probe_id'], h['strand'], h['rname'])
        dominated = False
        for k in by_group.get(g, ()):
            if k['edits'] < h['edits'] and k['rstart'] < h['rend'] and h['rstart'] < k['rend']:
                dominated = True
                break
        if dominated:
            continue
        by_group.setdefault(g, []).append(h)
        kept.append(h)

    for h in kept:
        probe_id = h['probe_id']
        if per_probe_cap and probe_hit_counts.get(probe_id, 0) >= per_probe_cap:
            continue
        dedup = (probe_id, h['rname'], h['pos1'], h['strand'])
        if dedup in seen:
            continue
        seen.add(dedup)
        flag = 16 if h['strand'] == 1 else 0
        qual = 'I' * h['L']
        out.write(f"{probe_id}\t{flag}\t{h['rname']}\t{h['pos1']}\t255\t"
                  f"{h['cigar']}\t*\t0\t0\t{h['seq']}\t{qual}\t"
                  f"NM:i:{h['edits']}\tMD:Z:{h['md']}\n")
        probe_hit_counts[probe_id] = probe_hit_counts.get(probe_id, 0) + 1


def _gpu_scan(cp, probes, groups, sam_out, e, B, seed_len, max_len, idx,
              qstrings, kmer_length, keys_info, device, workers, batch_mb,
              capacity, per_probe_cap):
    """Scan each group's files once: align (write SAM) and, if kmer_length is
    set, canonical-count the probe k-mers per group in the same pass.

    Returns {group_label: {kmer: count}} (empty when not counting), or None on
    hit-buffer overflow."""
    q_probe = idx['q_probe']
    q_strand = idx['q_strand']
    batch_bytes = max(1 << 20, batch_mb * 1024 * 1024)
    threads = 256
    n_seeds = int(idx['sorted_seeds'].size)
    do_count = kmer_length is not None
    gapped = B > 0
    group_counts = {}

    with cp.cuda.Device(device):
        d = {k: cp.asarray(v) for k, v in (
            ('sorted_seeds', idx['sorted_seeds']),
            ('post_start', idx['post_start']),
            ('post_q', idx['post_q']),
            ('post_soff', idx['post_soff']),
            ('qcodes', idx['qcodes']),
            ('q_off', idx['q_off']),
            ('q_len', idx['q_len']),
        )}
        kernel = (_get_gapped_kernel(cp, max_len, B) if gapped
                  else _get_kernel(cp, e))
        if do_count:
            sorted_keys, canon_to_orig, key_index = keys_info
            keys_d = cp.asarray(sorted_keys)
            ckernel = gpu_kmer.count_kernel(cp)

        def _ref_source(files):
            """Yield (ref_d on GPU, n_ref, names, starts, lens), uploading each
            batch's codes straight to the GPU from its shared-memory buffer."""
            for codes, names, starts, lens in _iter_batches(
                    list(files), workers, batch_bytes):
                yield cp.asarray(codes), int(codes.size), names, starts, lens

        seen = set()
        probe_hit_counts = {}
        with open(sam_out, 'w', encoding='utf-8') as out:
            out.write('@HD\tVN:1.4\tSO:unsorted\n')
            for label, files in groups:
                counts_d = (cp.zeros(sorted_keys.shape[0], dtype=cp.uint64)
                            if do_count else None)
                for ref_d, n_ref, names, starts, lens in _ref_source(list(files)):
                    if n_ref >= seed_len:
                        batch_fn = _align_batch_gapped if gapped else _align_batch
                        ok = batch_fn(
                            cp, out, ref_d, n_ref, seed_len, e, d, n_seeds,
                            kernel, threads, capacity, q_probe, q_strand,
                            probes, qstrings, names, starts, lens, seen,
                            probe_hit_counts, per_probe_cap)
                        if not ok:
                            return None
                    if do_count and n_ref >= kmer_length:
                        n_starts = n_ref - kmer_length + 1
                        cblocks = (n_starts + threads - 1) // threads
                        ckernel((cblocks,), (threads,), (
                            ref_d, np.int64(n_ref), np.int32(kmer_length),
                            keys_d, np.int32(keys_d.size), counts_d,
                            np.int64(n_starts)))
                if do_count:
                    cp.cuda.Stream.null.synchronize()
                    group_counts[label] = gpu_kmer.counts_array_to_dict(
                        cp.asnumpy(counts_d), canon_to_orig, key_index)
    return group_counts


def _prepare_scan(probes, max_mismatches, max_bulges=0):
    """Validate params and build the alignment seed index. Returns
    (e, B, seed_len, max_len, idx, qstrings) or None if unsupported (caller
    raises, since razers3 is no longer available as a fallback)."""
    e = int(max_mismatches)
    B = int(max_bulges)
    if e > MAX_E or B > MAX_B or B > e:
        return None
    min_len = min(len(s) for _, s in probes)
    max_len = max(len(s) for _, s in probes)
    if B > 0 and max_len > MAX_GAP_LEN:
        return None
    seed_len = min(min_len // (e + 1), MAX_SEED_LEN)
    if seed_len < MIN_SEED_LEN:
        return None
    idx = _build_probe_index(probes, e, seed_len)
    if idx['sorted_seeds'].size == 0:
        return None
    return e, B, seed_len, max_len, idx, _probe_query_strings(probes)


def gpu_align(probes_file, chrom_files, sam_out, max_mismatches, max_bulges=0,
              device=0, workers=None, batch_mb=256, capacity=8_000_000,
              per_probe_cap=100):
    """Align probes to reference FASTA files, writing a RazerS3-style SAM.

    Handles both the ungapped (max_bulges == 0) and gapped cases. Returns True
    on success, False if the parameters are unsupported or a hit-buffer overflow
    could drop hits (the caller treats False as a hard error).
    """
    from Bio import SeqIO
    cp = gpu_kmer._cupy()
    if cp is None:
        return False
    probes = [(r.id, str(r.seq).upper()) for r in SeqIO.parse(probes_file, 'fasta')]
    if not probes:
        with open(sam_out, 'w', encoding='utf-8') as f:
            f.write('@HD\tVN:1.4\tSO:unsorted\n')
        return True
    prep = _prepare_scan(probes, max_mismatches, max_bulges)
    if prep is None:
        return False
    e, B, seed_len, max_len, idx, qstrings = prep
    if workers is None:
        workers = os.cpu_count() or 1
    res = _gpu_scan(cp, probes, [('all', chrom_files)], sam_out, e, B, seed_len,
                    max_len, idx, qstrings, None, None, device, workers,
                    batch_mb, capacity, per_probe_cap)
    if res is None:
        try:
            os.remove(sam_out)
        except OSError:
            pass
        return False
    return True


def gpu_align_and_count(probes_file, groups, sam_out, max_mismatches,
                        kmer_length, max_bulges=0, device=0, workers=None,
                        batch_mb=256, capacity=8_000_000, per_probe_cap=None):
    """Fused single scan: align probes (write SAM) AND canonical-count every
    candidate probe's `kmer_length`-mers per group, in one pass over the
    reference.

    `groups` is a list of (label, files). Returns (True, {label: {kmer: count}})
    on success, or (False, None) if the parameters are unsupported (the caller
    treats that as a hard error since razers3 is gone).
    """
    from Bio import SeqIO
    cp = gpu_kmer._cupy()
    if cp is None:
        return False, None
    probes = [(r.id, str(r.seq).upper()) for r in SeqIO.parse(probes_file, 'fasta')]
    if not probes:
        with open(sam_out, 'w', encoding='utf-8') as f:
            f.write('@HD\tVN:1.4\tSO:unsorted\n')
        return True, {label: {} for label, _ in groups}
    prep = _prepare_scan(probes, max_mismatches, max_bulges)
    if prep is None:
        return False, None
    e, B, seed_len, max_len, idx, qstrings = prep
    if kmer_length > gpu_kmer.MAX_GPU_K:
        return False, None
    if workers is None:
        workers = os.cpu_count() or 1

    unique_kmers = set()
    for _pid, seq in probes:
        for i in range(len(seq) - kmer_length + 1):
            unique_kmers.add(seq[i:i + kmer_length])
    keys_info = gpu_kmer.build_canonical_keys(unique_kmers, kmer_length)
    if keys_info[0].size == 0:
        return False, None

    res = _gpu_scan(cp, probes, list(groups), sam_out, e, B, seed_len, max_len,
                    idx, qstrings, kmer_length, keys_info, device, workers,
                    batch_mb, capacity, per_probe_cap)
    if res is None:
        try:
            os.remove(sam_out)
        except OSError:
            pass
        return False, None
    return True, res


def _write_hits(out, h_q, h_pos, h_nmm, h_mmpos, h_mmref, e,
                q_probe, q_strand, probes, qstrings,
                names, contig_starts, contig_lens,
                seen, probe_hit_counts, per_probe_cap):
    for r in range(h_q.shape[0]):
        q = int(h_q[r])
        pi = int(q_probe[q])
        strand = int(q_strand[q])
        gpos = int(h_pos[r])
        nmm = int(h_nmm[r])

        ci = int(np.searchsorted(contig_starts, gpos, side='right') - 1)
        if ci < 0 or ci >= len(names):
            continue
        cstart = int(contig_starts[ci])
        clen = int(contig_lens[ci])
        local = gpos - cstart
        L = len(qstrings[pi][0])
        if local < 0 or local + L > clen:
            continue

        rname = names[ci]
        probe_id = probes[pi][0]
        if per_probe_cap:
            c = probe_hit_counts.get(probe_id, 0)
            if c >= per_probe_cap:
                continue
        pos1 = local + 1
        dedup = (probe_id, rname, pos1, strand)
        if dedup in seen:
            continue
        seen.add(dedup)

        mmpos = [int(h_mmpos[r][t]) for t in range(e) if int(h_mmpos[r][t]) >= 0]
        mmref = [int(h_mmref[r][t]) for t in range(len(mmpos))]
        query_seq = qstrings[pi][0] if strand == 0 else qstrings[pi][1]
        md, seq = _md_and_seq(query_seq, strand, mmpos, mmref)
        flag = 16 if strand == 1 else 0

        qual = 'I' * L
        out.write(f'{probe_id}\t{flag}\t{rname}\t{pos1}\t255\t{L}M\t*\t0\t0\t'
                  f'{seq}\t{qual}\tNM:i:{nmm}\tMD:Z:{md}\n')
        probe_hit_counts[probe_id] = probe_hit_counts.get(probe_id, 0) + 1
