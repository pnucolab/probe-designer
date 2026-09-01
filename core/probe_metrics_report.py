"""
Writes the per-probe metrics table consumed by the frontend.

Measures each probe with ProbeMetricsCalculator and renders a fixed-width
report. No probe is ranked or filtered here.
"""

from Bio import SeqIO
from probe_metrics import ProbeMetricsCalculator


def _reverse_complement(sequence: str) -> str:
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(complement.get(base, 'N') for base in reversed(sequence.upper()))


def write_probe_metrics_report(fasta_file: str, output_csv: str) -> int:
    """
    Measure all probes in a FASTA file and save results to a formatted table.

    Args:
        fasta_file: Path to input FASTA file with probe sequences
        output_csv: Path to output file for metrics (will be .txt format)

    Returns:
        Number of probes measured
    """

    calculator = ProbeMetricsCalculator(
        temperature_celsius=37.0,
        formamide_percent=50.0,
        na_concentration_mM=390.0,
        dnac1_nM=25.0,
        dnac2_nM=25.0,
        max_homopolymer=5,
        enable_hard_filters=True
    )

    sequences = []
    probe_ids = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences.append(str(record.seq))
        probe_ids.append(record.id)

    if not sequences:
        print(f"No sequences found in {fasta_file}")
        return 0


    results = calculator.calculate_metrics_for_set(sequences)

    seq_to_id = {seq: probe_id for seq, probe_id in zip(sequences, probe_ids)}
    for result in results:
        result['probe_id'] = seq_to_id[result['sequence']]

    output_file = output_csv.replace('.csv', '.txt')

    with open(output_file, 'w', encoding='utf-8') as f:

        f.write("=" * 245 + "\n")
        f.write("NON-ALIGNED PROBE METRICS\n")
        f.write("=" * 245 + "\n\n")

        header_parts = [
            "Probe ID",
            "Sequence",
            "Tm(°C)",
            "GC%",
            "Len",
            "Complexity",
            "SecStruct",
            "Homopoly",
            "OrderReadySeq"
        ]
        f.write("  ".join(header_parts) + "\n")
        f.write("-" * 245 + "\n")
        for result in results:
            row_parts = [
                result['probe_id'],
                result['sequence'],
                f"{result['tm']:.2f}",
                f"{result['gc_content']:.1f}",
                str(result['probe_length']),
                f"{result['complexity']:.3f}",
                f"{result['secondary_structure_penalty']:.3f}",
                'Yes' if result['has_homopolymer'] else 'No',
                _reverse_complement(result['sequence'])
            ]
            f.write("  ".join(row_parts) + "\n")


    print(f"\nMeasured {len(results)} probes")
    print(f"Results saved to: {output_file}")
