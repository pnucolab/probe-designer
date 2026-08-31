#!/usr/bin/env python
"""Pre-build Jellyfish k-mer databases for all supported species and k-mer lengths.

Usage:
    python build_kmer_dbs.py                          # Build all (14-19mer) for all species
    python build_kmer_dbs.py --species gut-microbe    # Build for one species
    python build_kmer_dbs.py --kmer-lengths 16 17 18  # Build specific k-mer lengths
    python build_kmer_dbs.py --threads 40             # Set parallelism
"""

import argparse
import os
import shlex
import subprocess
import time


SUPPORTED_SPECIES = [
    "gut-microbe",
    "human-oral-microbiome",
    "human-skin-microbiome",
    "human-vaginal-microbiome",
    "mouse-gut-microbiome",
    "human",
    "mouse",
]

KMER_RANGE = range(14, 20)  # 14 to 19


def get_genome_files(species):
    """Get genome FASTA files for a species."""
    microbiome_species = [
        "gut-microbe", "human-oral-microbiome", "human-skin-microbiome",
        "human-vaginal-microbiome", "mouse-gut-microbiome",
    ]

    if species in microbiome_species:
        genome_dir = os.path.join('data', species)
    else:
        genome_dir = os.path.join('data', 'gencode_data', species, 'transcript_chunks')

    if not os.path.exists(genome_dir):
        return [], genome_dir

    files = [
        os.path.join(genome_dir, f) for f in os.listdir(genome_dir)
        if f.endswith('.fa') or f.endswith('.fna')
    ]
    return sorted(files), genome_dir


def build_db(genome_files, k, db_path, threads):
    """Build a single Jellyfish k-mer database."""
    import tempfile
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as flist:
        for gf in genome_files:
            flist.write(gf + '\n')
        file_list_path = flist.name
    try:
        # shell=True is required: this is a real shell pipeline feeding the
        # concatenated genomes to Jellyfish through /dev/fd/0. Paths are passed
        # through shlex.quote() instead of hand-written single quotes so that a
        # path containing a quote cannot break out into arbitrary commands
        # (db_path is derived from the --species CLI value).
        cmd = (
            f"xargs cat < {shlex.quote(file_list_path)} | jellyfish count "
            f"-m {int(k)} -s 5G -t {int(threads)} -C -o {shlex.quote(db_path)} /dev/fd/0"
        )
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=False)
        return proc.returncode == 0
    finally:
        os.remove(file_list_path)


def main():
    parser = argparse.ArgumentParser(description="Pre-build Jellyfish k-mer databases")
    parser.add_argument("--species", nargs='+', default=None,
                       help=f"Species to build for (default: all). Options: {', '.join(SUPPORTED_SPECIES)}")
    parser.add_argument("--kmer-lengths", type=int, nargs='+', default=None,
                       help=f"K-mer lengths to build (default: {min(KMER_RANGE)}-{max(KMER_RANGE)})")
    parser.add_argument("--threads", type=int, default=40, help="Number of threads (default: 40)")
    args = parser.parse_args()

    species_list = args.species or SUPPORTED_SPECIES
    kmer_lengths = args.kmer_lengths or list(KMER_RANGE)

    for k in kmer_lengths:
        if k < 14 or k > 19:
            print(f"Error: k-mer length {k} out of supported range (14-19)")
            return

    total_builds = 0
    skipped = 0

    for species in species_list:
        genome_files, genome_dir = get_genome_files(species)
        if not genome_files:
            print(f"\n[{species}] No genome files found in {genome_dir}, skipping")
            continue

        total_size = sum(os.path.getsize(f) for f in genome_files) / (1024**3)
        print(f"\n{'='*60}")
        print(f"[{species}] {len(genome_files)} genome files ({total_size:.1f} GB)")
        print(f"{'='*60}")

        for k in kmer_lengths:
            db_path = os.path.join(genome_dir, f'kmer_{k}mer_counts.jf')

            if os.path.exists(db_path):
                size = os.path.getsize(db_path) / (1024**3)
                print(f"  {k}-mer DB already exists ({size:.1f} GB): {db_path}")
                skipped += 1
                continue

            print(f"  Building {k}-mer DB...", end='', flush=True)
            start = time.time()
            success = build_db(genome_files, k, db_path, args.threads)
            elapsed = time.time() - start

            if success:
                size = os.path.getsize(db_path) / (1024**3)
                print(f" done in {elapsed:.0f}s ({size:.1f} GB)")
                total_builds += 1
            else:
                print(f" FAILED after {elapsed:.0f}s")

    print(f"\nSummary: {total_builds} databases built, {skipped} already existed")


if __name__ == '__main__':
    main()
