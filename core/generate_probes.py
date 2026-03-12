from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import argparse
import glob

def calculate_gc_content(sequence):
    """Calculate GC content as percentage"""
    sequence = sequence.upper()
    g_count = sequence.count('G')
    c_count = sequence.count('C')
    total_length = len(sequence)
    if total_length == 0:
        return 0.0
    gc_percentage = (g_count + c_count) / total_length * 100
    return gc_percentage

def generate_candidate_probes(fasta_file, probe_length=36, step_size=3, min_gc=40, max_gc=80):
    candidate_probes = []
    total_checked = 0
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = str(record.seq).upper()
        for i in range(0, len(seq) - probe_length + 1, step_size):
            probe = seq[i:i + probe_length]
            total_checked += 1
            gc_content = calculate_gc_content(probe)
            if min_gc <= gc_content <= max_gc:
                candidate_probes.append({
                    'isoform_id': record.id,
                    'start': i+1,
                    'end': i + probe_length,
                    'probe_seq': probe
                })
    return candidate_probes, total_checked

def save_probes_to_fasta(probes, output_file):
    fasta_records = []
    for idx, p in enumerate(probes):
        probe_id = f"probe_{idx}|start={p['start']}|end={p['end']}|transcript:{p['isoform_id']}"
        fasta_records.append(
            SeqRecord(
                Seq(p['probe_seq']),
                id=probe_id,
                description=""
            )
        )
    SeqIO.write(fasta_records, output_file, "fasta")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate candidate probes from FASTA sequences")
    parser.add_argument("--probe-length", type=int, default=36,
                        help="Length of probes to generate (default: 36)")
    parser.add_argument("--step-size", type=int, default=1,
                        help="Step size for sliding window (default: 1)")
    parser.add_argument("--min-gc", type=float, default=40,
                        help="Minimum GC content percentage (default: 40)")
    parser.add_argument("--max-gc", type=float, default=80,
                        help="Maximum GC content percentage (default: 80)")
    parser.add_argument("--max-mismatches", type=int, default=2,
                        help="Maximum number of mismatches allowed in alignments (default: 2)")
    parser.add_argument("--species", default="human",
                        help="Species name (default: human)")
    parser.add_argument("--input", default=None,
                        help="Input FASTA file (default: (required) provide --input or use pipeline to pass custom FASTA)")
    parser.add_argument("--output", default="outputs/candidate_probes.fa",
                        help="Output FASTA file (default: outputs/candidate_probes.fa)")
    
    args = parser.parse_args()
    
    fasta_path = args.input
    probes_fasta = args.output
    
    os.makedirs(os.path.dirname(probes_fasta), exist_ok=True)
    
    probes, total_checked = generate_candidate_probes(
        fasta_path,
        probe_length=args.probe_length,
        step_size=args.step_size,
        min_gc=args.min_gc,
        max_gc=args.max_gc
    )
    
    print(f"\nChecked {total_checked} candidate positions.")
    print(f"Found {len(probes)} candidate probes passing GC filter ({args.min_gc}-{args.max_gc}%, length={args.probe_length}).")
    
    save_probes_to_fasta(probes, probes_fasta)
    print(f"Saved {len(probes)} probes to {probes_fasta}")
