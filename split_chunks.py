#!/usr/bin/env python3
"""Split each gut-microbe .fna chunk into parts targeting ~8MB per file, split by FASTA headers."""

import os
import glob
import math

SRC_DIR = "data/gut-microbe"
OUT_DIR = os.path.join(SRC_DIR, "split_chunks")
TARGET_SIZE = 8 * 1024 * 1024  # 8 MB

# Clean previous split
if os.path.exists(OUT_DIR):
    for f in os.listdir(OUT_DIR):
        os.remove(os.path.join(OUT_DIR, f))
os.makedirs(OUT_DIR, exist_ok=True)

chunk_files = sorted(glob.glob(os.path.join(SRC_DIR, "chunk_*.fna")))
print(f"Found {len(chunk_files)} chunk files to split (target: {TARGET_SIZE // (1024*1024)} MB per part)\n")

total_parts = 0

for chunk_path in chunk_files:
    basename = os.path.splitext(os.path.basename(chunk_path))[0]
    file_size = os.path.getsize(chunk_path)
    num_parts = max(1, math.ceil(file_size / TARGET_SIZE))

    # Read sequences by FASTA header
    sequences = []
    current = []
    with open(chunk_path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                if current:
                    sequences.append(current)
                current = [line]
            else:
                current.append(line)
    if current:
        sequences.append(current)

    total_seqs = len(sequences)
    seqs_per_part = math.ceil(total_seqs / num_parts)

    parts_written = 0
    for i in range(num_parts):
        start = i * seqs_per_part
        end = min(start + seqs_per_part, total_seqs)
        part_seqs = sequences[start:end]
        if not part_seqs:
            continue
        part_path = os.path.join(OUT_DIR, f"{basename}_part{i:03d}.fna")
        with open(part_path, 'w') as f:
            for seq in part_seqs:
                f.writelines(seq)
        parts_written += 1

    total_parts += parts_written
    print(f"{basename}: {file_size / (1024*1024):.1f} MB, {total_seqs} seqs -> {parts_written} parts")

print(f"\nDone. {total_parts} split chunks saved to {OUT_DIR}/")
