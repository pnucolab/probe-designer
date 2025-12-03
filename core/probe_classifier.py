"""
Determines if probes are safe based on alignment patterns.
"""

import re
from typing import Dict


def classify_probes_from_sam(sam_file: str, is_microbiome: bool = False) -> Dict:
    """
    Classify probes as safe or risky based on SAM alignments.
    """
    probe_data = {}
    sp_pattern = re.compile(r'SP:Z:(\S+)')
    nm_pattern = re.compile(r'NM:i:(\d+)')
    
    with open(sam_file, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('@'):
                continue
                
            parts = line.strip().split('\t')
            if len(parts) < 11:
                continue
            
            probe_id_full = parts[0]
            probe_id_match = re.search(r'probe_\d+', probe_id_full)
            if not probe_id_match:
                continue
            probe_id = probe_id_match.group(0)  
            nm_match = nm_pattern.search(line)
            if not nm_match:
                continue
            mismatches = int(nm_match.group(1))
            species_value = None
            if is_microbiome:
                sp_match = sp_pattern.search(line)
                if sp_match:
                    species_value = sp_match.group(1)
            if probe_id not in probe_data:
                probe_data[probe_id] = {
                    'total': 0,
                    'min_mm': float('inf'),
                    'is_microbiome_species': False
                }
            
            probe_data[probe_id]['total'] += 1
            probe_data[probe_id]['min_mm'] = min(probe_data[probe_id]['min_mm'], mismatches)
            
            if species_value and species_value != 'Unknown':
                probe_data[probe_id]['is_microbiome_species'] = True
    
    results = {}
    for probe_id, data in probe_data.items():
        if data['is_microbiome_species'] and data['total'] == 1 and data['min_mm'] <= 2:
            status = 'safe'
        elif data['min_mm'] <= 1:
            status = 'high_risk'
        else:
            status = 'medium_risk'
        
        results[probe_id] = {
            'status': status,
            'total_alignments': data['total'],
            'min_mismatches': data['min_mm'] if data['min_mm'] != float('inf') else 0
        }
    
    return results