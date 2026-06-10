#!/usr/bin/env python3
"""Split the human transcriptome FASTA into ~7 MB chunks for RazerS3 alignment.

Default: read /home/abyot/projects/test/oligominer/OligoMiner/GRCh38.p14.genome.release_49_transcripts.fa
         (the unfiltered version, matching OligoMiner's bowtie2 index content)
         and write 94 chunks into data/gencode_data/human/transcript_chunks/.
"""

import os
import shutil
import sys
import math

SRC = "/home/abyot/projects/test/oligominer/OligoMiner/GRCh38.p14.genome.release_49_transcripts.fa"
OUT_DIR = "data/gencode_data/human/transcript_chunks"
N_CHUNKS = 94
CHUNK_PREFIX = "GRCh38.p14.genome.release_49_transcripts.part"


def main():
    if not os.path.exists(SRC):
        sys.exit(f"ERROR: source not found: {SRC}")

    src_size = os.path.getsize(SRC)
    target_chunk_size = math.ceil(src_size / N_CHUNKS)
    print(f"Source: {SRC}")
    print(f"  size: {src_size / 1024 / 1024:.1f} MB")
    print(f"  splitting into ~{N_CHUNKS} chunks of ~{target_chunk_size / 1024 / 1024:.1f} MB each")

    # Move existing chunks aside in case the user wants to revert.
    existing_kmer_dbs = [f for f in os.listdir(OUT_DIR) if f.endswith('.jf')] if os.path.isdir(OUT_DIR) else []
    if os.path.isdir(OUT_DIR):
        old_fa = [f for f in os.listdir(OUT_DIR) if f.endswith('.fa') or f.endswith('.fna')]
        for f in old_fa:
            os.remove(os.path.join(OUT_DIR, f))
        print(f"  removed {len(old_fa)} existing chunk file(s); kept {len(existing_kmer_dbs)} Jellyfish DB(s)")
    else:
        os.makedirs(OUT_DIR, exist_ok=True)

    chunk_idx = 0
    chunk_path = lambda i: os.path.join(OUT_DIR, f"{CHUNK_PREFIX}{i:03d}.fa")
    out = open(chunk_path(chunk_idx), 'w', encoding='utf-8')
    bytes_in_chunk = 0
    seqs_in_chunk = 0
    total_seqs = 0
    chunks_written = 1

    with open(SRC, 'r', encoding='utf-8') as f:
        current_seq_lines = []
        for line in f:
            if line.startswith('>'):
                # Flush previous sequence to current chunk.
                if current_seq_lines:
                    seq_bytes = sum(len(l) for l in current_seq_lines)
                    if (bytes_in_chunk + seq_bytes > target_chunk_size
                            and chunk_idx + 1 < N_CHUNKS
                            and seqs_in_chunk > 0):
                        # Roll over to next chunk.
                        out.close()
                        chunk_idx += 1
                        out = open(chunk_path(chunk_idx), 'w', encoding='utf-8')
                        chunks_written += 1
                        bytes_in_chunk = 0
                        seqs_in_chunk = 0
                    out.writelines(current_seq_lines)
                    bytes_in_chunk += seq_bytes
                    seqs_in_chunk += 1
                    total_seqs += 1
                current_seq_lines = [line]
            else:
                current_seq_lines.append(line)
        # Final sequence.
        if current_seq_lines:
            out.writelines(current_seq_lines)
            total_seqs += 1
            seqs_in_chunk += 1

    out.close()
    print(f"\nWrote {chunks_written} chunks containing {total_seqs:,} sequences")

    # Cross-check.
    chunk_seq_total = 0
    for i in range(chunks_written):
        p = chunk_path(i)
        with open(p, 'r') as f:
            n = sum(1 for line in f if line.startswith('>'))
        chunk_seq_total += n
    print(f"  cross-check seq count: {chunk_seq_total:,} (expected {total_seqs:,})")

    if existing_kmer_dbs:
        print(f"\nNote: {len(existing_kmer_dbs)} Jellyfish DB(s) preserved in {OUT_DIR}/")
        print(f"  they were built from the .filtered.fa version (95 fewer transcripts)")
        print(f"  rebuild via build_kmer_dbs.py if you want them aligned with the new chunks")


if __name__ == "__main__":
    main()
