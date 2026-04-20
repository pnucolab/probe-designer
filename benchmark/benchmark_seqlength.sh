#!/usr/bin/env bash
#
# Benchmark probe_designer.py execution time across sequence lengths.
# Generates truncated/tiled sequences from the input to reach target lengths.
#
# Usage:
#   ./benchmark_seqlength.sh <input_fasta> <species> [extra_args...]
#
# Example:
#   ./benchmark_seqlength.sh gene_sequences/052837ce_pasted_gene.fasta human

set -euo pipefail

if [ $# -lt 2 ]; then
    echo "Usage: $0 <input_fasta> <species> [extra_args...]"
    exit 1
fi

INPUT_FASTA="$1"
SPECIES="$2"
shift 2
EXTRA_ARGS=("$@")

if [ ! -f "$INPUT_FASTA" ]; then
    echo "Error: Input file not found: $INPUT_FASTA"
    exit 1
fi

export PATH="/home/abyot/miniconda3/envs/probeset/bin:$PATH"
SEQ_LENGTHS_KB=(1 5 10 15 20 25 30 40 50 60 80 100)
RESULTS_FILE="benchmark_seqlength_results.txt"

# Extract header and raw sequence from input
HEADER=$(head -1 "$INPUT_FASTA")
RAW_SEQ=$(grep -v '^>' "$INPUT_FASTA" | tr -d '\n')
RAW_LEN=${#RAW_SEQ}

echo "========================================"
echo "Sequence Length Benchmark"
echo "========================================"
echo "Input:   $INPUT_FASTA (${RAW_LEN} bp)"
echo "Species: $SPECIES"
echo "Target lengths (kb): ${SEQ_LENGTHS_KB[*]}"
echo "Extra args: ${EXTRA_ARGS[*]:-none}"
echo "========================================"
echo ""

if [ ! -f "$RESULTS_FILE" ]; then
    printf "%-20s %-20s\n" "Seq Length (kb)" "Time (seconds)" | tee "$RESULTS_FILE"
    printf "%-20s %-20s\n" "---------------" "--------------" | tee -a "$RESULTS_FILE"
else
    printf "%-20s %-20s\n" "Seq Length (kb)" "Time (seconds)"
    printf "%-20s %-20s\n" "---------------" "--------------"
fi

for KB in "${SEQ_LENGTHS_KB[@]}"; do
    TARGET_LEN=$((KB * 1000))
    TASK_ID="benchmark_seq${KB}kb"
    mkdir -p "outputs/alignments/${TASK_ID}"

    # Build sequence of target length by tiling the input
    TILED_SEQ=""
    while [ ${#TILED_SEQ} -lt "$TARGET_LEN" ]; do
        TILED_SEQ="${TILED_SEQ}${RAW_SEQ}"
    done
    TILED_SEQ="${TILED_SEQ:0:$TARGET_LEN}"

    # Write temp FASTA
    TEMP_FASTA=$(mktemp /tmp/bench_seq_XXXXXX.fa)
    echo "${HEADER}" > "$TEMP_FASTA"
    echo "$TILED_SEQ" >> "$TEMP_FASTA"

    SEQUENCE=$(cat "$TEMP_FASTA")

    echo ""
    echo ">>> Running probe_designer.py with sequence length ${KB} kb (${TARGET_LEN} bp) ..."

    START=$(date +%s%N)

    /home/abyot/miniconda3/envs/probeset/bin/python probe_designer.py \
        --gene-sequence "$SEQUENCE" \
        --species "$SPECIES" \
        --task-id "$TASK_ID" \
        "${EXTRA_ARGS[@]}" \
        > "outputs/alignments/${TASK_ID}/run.log" 2>&1 || true

    END=$(date +%s%N)
    ELAPSED=$(( (END - START) / 1000000000 )).$(( ((END - START) % 1000000000) / 1000000 ))

    printf "%-20s %-20s\n" "$KB" "$ELAPSED" | tee -a "$RESULTS_FILE"

    rm -f "$TEMP_FASTA"
done

echo ""
echo "========================================"
echo "Summary"
echo "========================================"
cat "$RESULTS_FILE"
echo ""
echo "Results saved to $RESULTS_FILE"
