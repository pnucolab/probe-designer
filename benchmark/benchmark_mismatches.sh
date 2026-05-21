#!/usr/bin/env bash
#
# Benchmark probe_designer.py execution time across max mismatch values.
#
# Usage:
#   ./benchmark_mismatches.sh <input_fasta> <species> [extra_args...]
#
# Example:
#   ./benchmark_mismatches.sh gene_sequences/052837ce_pasted_gene.fasta human
#   ./benchmark_mismatches.sh gene_sequences/my_gene.fasta gut-microbe --align-microbiome

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
MISMATCHES=(0 1 2)
RESULTS_FILE="benchmark_mismatches_results.txt"
SEQUENCE=$(cat "$INPUT_FASTA")

echo "========================================"
echo "Max Mismatches Benchmark"
echo "========================================"
echo "Input:   $INPUT_FASTA"
echo "Species: $SPECIES"
echo "Max mismatches: ${MISMATCHES[*]}"
echo "Extra args: ${EXTRA_ARGS[*]:-none}"
echo "========================================"
echo ""

if [ ! -f "$RESULTS_FILE" ]; then
    printf "%-20s %-20s\n" "Max Mismatches" "Time (seconds)" | tee "$RESULTS_FILE"
    printf "%-20s %-20s\n" "--------------" "--------------" | tee -a "$RESULTS_FILE"
else
    printf "%-20s %-20s\n" "Max Mismatches" "Time (seconds)"
    printf "%-20s %-20s\n" "--------------" "--------------"
fi

for MM in "${MISMATCHES[@]}"; do
    TASK_ID="benchmark_mm${MM}"
    mkdir -p "output/alignments/${TASK_ID}"
    echo ""
    echo ">>> Running probe_designer.py with --max-mismatches $MM ..."

    START=$(date +%s%N)

    /home/abyot/miniconda3/envs/probeset/bin/python probe_designer.py \
        --gene-sequence "$SEQUENCE" \
        --species "$SPECIES" \
        --max-mismatches "$MM" \
        --task-id "$TASK_ID" \
        "${EXTRA_ARGS[@]}" \
        > "output/alignments/${TASK_ID}/run.log" 2>&1 || true

    END=$(date +%s%N)
    ELAPSED=$(( (END - START) / 1000000000 )).$(( ((END - START) % 1000000000) / 1000000 ))

    printf "%-20s %-20s\n" "$MM" "$ELAPSED" | tee -a "$RESULTS_FILE"
done

echo ""
echo "========================================"
echo "Summary"
echo "========================================"
cat "$RESULTS_FILE"
echo ""
echo "Results saved to $RESULTS_FILE"
