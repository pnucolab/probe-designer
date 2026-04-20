#!/usr/bin/env bash
#
# Benchmark probe_designer.py execution time across probe lengths.
#
# Usage:
#   ./benchmark_probe_length.sh <input_fasta> <species> [extra_args...]
#
# Example:
#   ./benchmark_probe_length.sh gene_sequences/052837ce_pasted_gene.fasta human
#   ./benchmark_probe_length.sh gene_sequences/my_gene.fasta gut-microbe --align-microbiome

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
PROBE_LENGTHS=(25 30 35)
RESULTS_FILE="benchmark_results.txt"
SEQUENCE=$(cat "$INPUT_FASTA")

echo "========================================"
echo "Probe Length Benchmark"
echo "========================================"
echo "Input:   $INPUT_FASTA"
echo "Species: $SPECIES"
echo "Probe lengths: ${PROBE_LENGTHS[*]}"
echo "Extra args: ${EXTRA_ARGS[*]:-none}"
echo "========================================"
echo ""

if [ ! -f "$RESULTS_FILE" ]; then
    printf "%-20s %-20s\n" "Probe Length (bp)" "Time (seconds)" | tee "$RESULTS_FILE"
    printf "%-20s %-20s\n" "-----------------" "--------------" | tee -a "$RESULTS_FILE"
else
    printf "%-20s %-20s\n" "Probe Length (bp)" "Time (seconds)"
    printf "%-20s %-20s\n" "-----------------" "--------------"
fi

for PL in "${PROBE_LENGTHS[@]}"; do
    TASK_ID="benchmark_pl${PL}"
    mkdir -p "outputs/alignments/${TASK_ID}"
    echo ""
    echo ">>> Running probe_designer.py with --probe-length $PL ..."

    START=$(date +%s%N)

    /home/abyot/miniconda3/envs/probeset/bin/python probe_designer.py \
        --gene-sequence "$SEQUENCE" \
        --species "$SPECIES" \
        --probe-length "$PL" \
        --task-id "$TASK_ID" \
        "${EXTRA_ARGS[@]}" \
        > "outputs/alignments/${TASK_ID}/run.log" 2>&1 || true

    END=$(date +%s%N)
    ELAPSED=$(( (END - START) / 1000000000 )).$(( ((END - START) % 1000000000) / 1000000 ))

    printf "%-20s %-20s\n" "$PL" "$ELAPSED" | tee -a "$RESULTS_FILE"
done

echo ""
echo "========================================"
echo "Summary"
echo "========================================"
cat "$RESULTS_FILE"
echo ""
echo "Results saved to $RESULTS_FILE"
