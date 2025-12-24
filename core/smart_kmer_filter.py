"""
Smart k-mer filtering with coding awareness.
Simple logic: Discard probes if k-mer matches coding region, keep if matches non-coding.
"""

import os
import re
import time
from pathlib import Path


def parse_gff_for_coding(gff_file, target_species):
    """
    Parse GFF/GTF file to extract coding regions by transcript.
    Returns: dict mapping transcript_id -> list of coding regions
    """
    coding_regions = {}
    
    if not os.path.exists(gff_file):
        print(f"Warning: GFF file not found: {gff_file}")
        return coding_regions
    
    print(f"Parsing GFF file for coding regions: {gff_file}")
    
    with open(gff_file, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            
            transcript_id = parts[0]
            feature_type = parts[2]
            start = int(parts[3])
            end = int(parts[4])
            
            if feature_type != 'coding':
                continue
            
            if transcript_id not in coding_regions:
                coding_regions[transcript_id] = []
            
            coding_regions[transcript_id].append({
                'start': start,
                'end': end
            })
    
    print(f"Loaded coding regions for {len(coding_regions)} transcripts")
    return coding_regions


def check_overlap(start1, end1, start2, end2):
    """Check if two genomic regions overlap."""
    return start1 <= end2 and end1 >= start2


def smart_kmer_filter(match_details, mer_to_probe, probes, k, species, output_base, chrom_files, razers, args):
    """
    Simple coding-aware filtering for microbiome species only.
    For human/mouse: all k-mer matches mark probe as unsafe (transcriptome is all coding)
    For microbiome: discard if k-mer matches coding, keep if matches non-coding.
    """
    
    microbiome_species_list = ["gut-microbe", "human-oral-microbiome", "human-skin-microbiome", 
                                "human-vaginal-microbiome", "mouse-gut-microbiome"]
    
    if species not in microbiome_species_list:
        print("Using simple k-mer filtering (transcriptome = all coding)")
        organizing_start = time.time()
        unsafe_probes = set()
        probe_match_details = {}
        
        for mer_id, _, _ in match_details:
            if mer_id in mer_to_probe:
                probe_list = mer_to_probe[mer_id]
                for probe_id, _ in probe_list:
                    unsafe_probes.add(probe_id)
                    if probe_id not in probe_match_details:
                        probe_match_details[probe_id] = {
                            'total_matches': 0,
                            'coding_matches': 0,
                            'non_coding_matches': 0
                        }
                    probe_match_details[probe_id]['total_matches'] += 1
                    probe_match_details[probe_id]['coding_matches'] += 1  
        organizing_elapsed = time.time() - organizing_start
        print(f"Simple filtering took {organizing_elapsed:.2f}s")
        print(f"  Unsafe probes (any match in transcriptome): {len(unsafe_probes)}")
        
        return unsafe_probes, probe_match_details, set()
    
    print("Using coding-aware filtering for microbiome species")
    
    gff_file = None
    gff_dir = os.path.join('data', species, 'gff')
    if os.path.exists(gff_dir):
        gff_files = [f for f in os.listdir(gff_dir) if f.endswith('.gff') or f.endswith('.gtf')]
        if gff_files:
            gff_file = os.path.join(gff_dir, gff_files[0])
            print(f"Found microbiome GFF file: {gff_file}")
    
    coding_regions = {}
    if gff_file:
        coding_regions = parse_gff_for_coding(gff_file, species)
    else:
        print("Warning: No GFF file found - using simple matching")
        
        unsafe_probes = set()
        probe_match_details = {}
        for mer_id, _, _ in match_details:
            if mer_id in mer_to_probe:
                probe_list = mer_to_probe[mer_id]
                for probe_id, _ in probe_list:
                    unsafe_probes.add(probe_id)
                    if probe_id not in probe_match_details:
                        probe_match_details[probe_id] = {'total_matches': 0, 'coding_matches': 0, 'non_coding_matches': 0}
                    probe_match_details[probe_id]['total_matches'] += 1
        return unsafe_probes, probe_match_details, set()
    
    organizing_start = time.time()
    unsafe_probes = set()
    probe_match_details = {}
    
    print("Processing matches with coding checking...")
    for mer_id, chromosome, position in match_details:
        if mer_id not in mer_to_probe:
            continue
        
        probe_list = mer_to_probe[mer_id]
        match_pos = int(position)
        match_end = match_pos + k
        
        is_in_coding = False
        
        if chromosome in coding_regions:
            for coding in coding_regions[chromosome]:
                if check_overlap(match_pos, match_end, coding['start'], coding['end']):
                    is_in_coding = True
                    break
        
        for probe_id, _ in probe_list:
            if probe_id not in probe_match_details:
                probe_match_details[probe_id] = {
                    'total_matches': 0,
                    'coding_matches': 0,
                    'non_coding_matches': 0
                }
            
            probe_match_details[probe_id]['total_matches'] += 1
            
            if is_in_coding:
                probe_match_details[probe_id]['coding_matches'] += 1
                unsafe_probes.add(probe_id)
            else:
                probe_match_details[probe_id]['non_coding_matches'] += 1
    
    organizing_elapsed = time.time() - organizing_start
    print(f"coding filtering took {organizing_elapsed:.2f}s")
    print(f"  Unsafe probes (coding matches): {len(unsafe_probes)}")
    print(f"  Safe probes (non-coding only): {len([p for p in probe_match_details if p not in unsafe_probes])}")
    
    return unsafe_probes, probe_match_details, set()