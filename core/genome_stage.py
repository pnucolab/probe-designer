"""Optional local-NVMe staging for reference genomes.

The genome data lives on NFS, where cold reads are ~0.15 GB/s (network
round-trips per file) versus ~2 GB/s on the local NVMe. Because the
reference is static, staging it once on local disk makes every job's scan
fast regardless of page-cache state, removing the cold-NFS variance that
makes wall-clock time swing (e.g. 52s cold vs 26s warm for the same job).

Enable by setting PROBESET_LOCAL_GENOME_ROOT to a local directory (e.g.
/scratch/probeset_genomes). Reference files are mirrored there under their
absolute path; the pipeline transparently reads the local copy when present
and falls back to the original path otherwise. Staging is done once (or when
genomes change) via `stage()` or the CLI — never inside a job's hot path.
"""

import os
import shutil

ENV_VAR = 'PROBESET_LOCAL_GENOME_ROOT'


def local_root():
    root = os.environ.get(ENV_VAR)
    return root or None


def _mirror_path(path, root):
    return os.path.join(root, os.path.abspath(path).lstrip(os.sep))


def local_path(path, root=None):
    """Return the local mirror of `path` if it exists, else `path` unchanged."""
    root = root or local_root()
    if not root:
        return path
    lp = _mirror_path(path, root)
    try:
        if os.path.exists(lp):
            return lp
    except OSError:
        pass
    return path


def map_files(files, root=None):
    """Redirect a list of reference files to their local mirror where present."""
    root = root or local_root()
    if not root:
        return list(files)
    return [local_path(f, root) for f in files]


def _needs_copy(src, dst):
    try:
        s, d = os.stat(src), os.stat(dst)
    except OSError:
        return True
    return s.st_size != d.st_size or s.st_mtime > d.st_mtime


def _copy_one(args):
    src, root = args
    dst = _mirror_path(src, root)
    if not _needs_copy(src, dst):
        return 0
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    tmp = f"{dst}.tmp.{os.getpid()}.{id(args)}"
    shutil.copy2(src, tmp)
    os.replace(tmp, dst)
    return 1


def stage(files, root=None, verbose=True, workers=16):
    """Copy `files` into the local mirror (size/mtime-skip), in parallel.

    Copying is I/O-bound (cold NFS read + local write), so threads overlap the
    network latency. Returns the number of files actually copied.
    """
    root = root or local_root()
    if not root:
        raise RuntimeError(f"{ENV_VAR} is not set; nowhere to stage")
    from concurrent.futures import ThreadPoolExecutor
    args = [(f, root) for f in files]
    with ThreadPoolExecutor(max_workers=workers) as ex:
        copied = sum(ex.map(_copy_one, args))
    if verbose:
        print(f"Staged {copied}/{len(files)} files to {root}")
    return copied


def _iter_fasta(paths):
    for p in paths:
        if os.path.isdir(p):
            for fn in sorted(os.listdir(p)):
                if fn.endswith(('.fa', '.fna', '.fasta')):
                    yield os.path.join(p, fn)
        elif p.endswith(('.fa', '.fna', '.fasta')):
            yield p


if __name__ == '__main__':
    import sys
    import time
    if len(sys.argv) < 2:
        print(f"usage: {ENV_VAR}=/scratch/probeset_genomes "
              f"python core/genome_stage.py <fasta-or-dir> [...]")
        sys.exit(1)
    if not local_root():
        print(f"Set {ENV_VAR} first (e.g. export {ENV_VAR}=/scratch/probeset_genomes)")
        sys.exit(1)
    targets = list(_iter_fasta(sys.argv[1:]))
    print(f"Staging {len(targets)} FASTA files to {local_root()} ...")
    t = time.time()
    n = stage(targets)
    print(f"Done in {time.time()-t:.1f}s ({n} copied, {len(targets)-n} up-to-date)")
