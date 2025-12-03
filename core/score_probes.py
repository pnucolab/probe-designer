"""
Utility module for scoring probe sequences using thermodynamic metrics.
Integrates with probe_designer.py pipeline.
"""

import csv
import os
import sys
from typing import List, Dict
from Bio import SeqIO
from scorer import ThermodynamicProbeScorer


def get_quality_label(score, rejected):
    """Get quality label based on score."""
    if rejected:
        return 'Rejected'
    if score >= 90:
        return 'Excellent'
    elif score >= 80:
        return 'Very Good'
    elif score >= 70:
        return 'Good'
    elif score >= 60:
        return 'Acceptable'
    else:
        return 'Marginal'

def score_and_save_probes(fasta_file: str, output_csv: str) -> int:
    """
    Score all probes in a FASTA file and save results to formatted table.
    
    Args:
        fasta_file: Path to input FASTA file with probe sequences
        output_csv: Path to output file for scores (will be .txt format)
        
    Returns:
        Number of probes scored
    """
   
    scorer = ThermodynamicProbeScorer(
        temperature_celsius=37.0,
        formamide_percent=50.0,
        na_concentration_mM=390.0,
        dnac1_nM=25.0,
        dnac2_nM=25.0,
        target_tm=47.0,
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
    

    results = scorer.score_probe_set(sequences)
    
    for i, result in enumerate(results):
        result['probe_id'] = probe_ids[i]
 
    output_file = output_csv.replace('.csv', '.txt')
  
    with open(output_file, 'w', encoding='utf-8') as f:
       
        f.write("=" * 245 + "\n")
        f.write("NON-ALIGNED PROBE SCORING RESULTS\n")
        f.write("=" * 245 + "\n\n")
        
        header_parts = [
            "Rank",
            "Probe ID",
            "Sequence",
            "Score",
            "Tm(°C)",
            "GC%",
            "Len",
            "Complexity",
            "SecStruct",
            "Homopoly",
            "Status"
        ]
        f.write("  ".join(header_parts) + "\n")
        f.write("-" * 245 + "\n")
        for rank, result in enumerate(results, 1):
            quality = get_quality_label(result['on_target_score'], result['rejected'])
            row_parts = [
                str(rank),
                result['probe_id'],
                result['sequence'],
                f"{result['on_target_score']:.2f}",
                f"{result['tm']:.2f}",
                f"{result['gc_content']:.1f}",
                str(result['probe_length']),
                f"{result['complexity']:.3f}",
                f"{result['secondary_structure_penalty']:.3f}",
                'Yes' if result['has_homopolymer'] else 'No',
                quality
            ]
            f.write("  ".join(row_parts) + "\n")
        

    print(f"\nScored {len(results)} probes")
    print(f"Results saved to: {output_file}")