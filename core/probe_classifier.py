"""
Determines if probes are safe based on alignment patterns.
"""

import re
from typing import Dict

def classify_probes_from_sam(sam_file: str, is_microbiome: bool = False) -> Dict:
    """
    Classify probes as safe or risky based on SAM alignments.
    Now considers gene-level specificity: multiple alignments to the same gene are safe,
    but alignments to different genes are risky.
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
            
            probe_id = parts[0]
            target_info = parts[2]  

            gene_name = None
            if len(parts) > 3:
                # Try to get gene name from the 4th column if it exists
                try:
                    gene_name = parts[3]
                except IndexError:
                    gene_name = target_info
            else:
                gene_name = target_info
                
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
                    'unique_genes': set(),
                    'is_microbiome_species': False,
                    'targets': []
                }
            probe_data[probe_id]['targets'].append(target_info) 
            probe_data[probe_id]['total'] += 1
            probe_data[probe_id]['min_mm'] = min(probe_data[probe_id]['min_mm'], mismatches)
            probe_data[probe_id]['unique_genes'].add(gene_name)
            
            if species_value and species_value != 'Unknown':
                probe_data[probe_id]['is_microbiome_species'] = True
    
    results = {}
    for probe_id, data in probe_data.items():
        unique_gene_count = len(data['unique_genes'])
        
        # For microbiome: check if ANY alignment is to a different species
        if is_microbiome and data['is_microbiome_species']:
            # Extract source transcript from probe_id
            probe_source = probe_id.split('|')[-1].replace('transcript:', '')
            
            # Check if ANY target is different from source
            is_offtarget = False
            for target in data['targets']:
                if probe_source not in target:
                    is_offtarget = True
                    break
            
            if is_offtarget:
                # Off-target to microbiome genome
                if data['min_mm'] <= 1:
                    status = 'high_risk'
                else:
                    status = 'medium_risk'
            else:
                # All alignments are self-alignments
                status = 'safe'
        
        # For regular species: safe if aligns to only one gene
        elif unique_gene_count == 1:
            probe_source = probe_id.split('|')[-1].replace('transcript:', '')
            is_self = any(probe_source in target for target in data['targets'])
            
            if is_self:
                status = 'safe'
            elif data['min_mm'] <= 1:
                status = 'high_risk'
            else:
                status = 'medium_risk'
        # High risk if perfect/near-perfect matches to multiple genes
        elif data['min_mm'] <= 1:
            status = 'high_risk'
        # Medium risk for other multi-gene alignments
        else:
            status = 'medium_risk'
        
        results[probe_id] = {
            'status': status,
            'total_alignments': data['total'],
            'unique_genes': unique_gene_count,
            'min_mismatches': data['min_mm'] if data['min_mm'] != float('inf') else 0
        }
    
    return results