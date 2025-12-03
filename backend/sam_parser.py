"""
SAM file parser for alignment visualization
Parses filtered_probe_alignments.sam and returns formatted data
"""

import re
from typing import List, Dict, Any


def parse_cigar(cigar: str) -> List[tuple]:
    """
    Parse CIGAR string into operations.
    Returns list of (count, operation) tuples.
    Example: "2M1I17M" -> [(2, 'M'), (1, 'I'), (17, 'M')]
    """
    pattern = re.compile(r'(\d+)([MIDNSHPX=])')
    return [(int(count), op) for count, op in pattern.findall(cigar)]


def parse_md_tag(md: str) -> List[tuple]:
    """
    Parse MD tag to find mismatch positions.
    Returns list of (position, type, value) tuples.
    Example: "5T10A3" -> positions where mismatches occur
    """
    mismatches = []
    pattern = re.compile(r'(\d+)|([A-Z]+)|\^([A-Z]+)')
    pos = 0
    
    for match in pattern.finditer(md):
        if match.group(1):  
            pos += int(match.group(1))
        elif match.group(2): 
            ref_base = match.group(2)
            for base in ref_base:
                mismatches.append((pos, 'substitution', base))
                pos += 1
        elif match.group(3): 
            deleted = match.group(3)
            mismatches.append((pos, 'deletion', deleted))
    
    return mismatches


def apply_mismatches_to_sequence(sequence: str, cigar: str, md: str, flag: int) -> str:
    """
    Apply mismatch highlighting to sequence.
    - Substitutions: lowercase
    - Insertions: lowercase
    - Deletions: insert '-' character
    """
    cigar_ops = parse_cigar(cigar)
    md_mismatches = parse_md_tag(md)
    
    result = list(sequence)
    seq_pos = 0
    ref_pos = 0
    
    insertion_positions = set()
    temp_seq_pos = 0
    for count, op in cigar_ops:
        if op == 'M': 
            temp_seq_pos += count
        elif op == 'I': 
            for i in range(count):
                insertion_positions.add(temp_seq_pos + i)
            temp_seq_pos += count
        elif op == 'D': 
            pass
    
    for pos in insertion_positions:
        if pos < len(result):
            result[pos] = result[pos].lower()
    
    seq_pos = 0
    for count, op in cigar_ops:
        if op == 'M':
            for i in range(count):
                for md_pos, md_type, _ in md_mismatches:
                    if md_type == 'substitution' and md_pos == ref_pos:
                        if seq_pos < len(result):
                            result[seq_pos] = result[seq_pos].lower()
                seq_pos += 1
                ref_pos += 1
        elif op == 'I':
            seq_pos += count
        elif op == 'D':
            for md_pos, md_type, deleted_bases in md_mismatches:
                if md_type == 'deletion' and md_pos == ref_pos:
                    for _ in deleted_bases:
                        result.insert(seq_pos, '-')
                        seq_pos += 1
            ref_pos += count
    
    return ''.join(result)


def count_sam_alignments(sam_path: str) -> int:
    """
    Count total number of valid alignments in SAM file (excluding headers).
    """
    count = 0
    try:
        with open(sam_path, 'r') as f:
            for line in f:
                if not line.startswith('@'):
                    parts = line.strip().split('\t')
                    if len(parts) >= 11:
                        count += 1
        return count
    except Exception as e:
        print(f"Error counting alignments: {e}")
        return 0


def parse_sam_file(sam_path: str, offset: int = 0, limit: int = None, mismatch_filter: int = None, probe_id_filter: str = None) -> List[Dict[str, Any]]:
    """
    Parse SAM file and extract alignment information with pagination support.
    Handles * sequences by tracking the last valid sequence per probe.
    """
    alignments = []
    skipped = 0
    sequence_cache = {}  
    if probe_id_filter:
        base_probe_id = probe_id_filter.split('|')[0].strip()
        probe_numeric = re.search(r'probe_(\d+)', probe_id_filter)
        probe_num = probe_numeric.group(1) if probe_numeric else None
    try:
        with open(sam_path, 'r') as f:
            for line in f:
                if line.startswith('@'):
                    continue
                
                parts = line.strip().split('\t')
                if len(parts) < 11:
                    continue

                qname = parts[0]
                flag = int(parts[1])
                rname = parts[2]

                # Detect format
                try:
                    pos = int(parts[3])
                    # Standard SAM: no gene_id
                    gene_id = None
                    cigar = parts[5]
                    seq = parts[9]
                except ValueError:
                    # Annotated SAM: has gene_id in column 3
                    gene_id = parts[3]
                    pos = int(parts[4])
                    cigar = parts[6]  # Shifted by 1
                    seq = parts[10]  
                
                if seq == '*':
                    if qname in sequence_cache:
                        seq = sequence_cache[qname]
                    else:
                        continue
                else:
                    sequence_cache[qname] = seq
                
                if skipped < offset:
                    skipped += 1
                    continue
                
                nm = None
                md = None
                species = None
                for field in parts[11:]:
                    if field.startswith('NM:i:'):
                        nm = int(field.split(':')[2])
                    elif field.startswith('MD:Z:'):
                        md = field.split(':')[2]
                    elif field.startswith('SP:Z:'):  
                        species = field.split(':')[2]
                
                if nm is None or md is None:
                    continue
                
                if mismatch_filter is not None and nm != mismatch_filter:
                    continue

                is_reverse = bool(flag & 0x10)
                strand = '-' if is_reverse else '+'
                
                formatted_seq = apply_mismatches_to_sequence(seq, cigar, md, flag)
                
                alignments.append({
                    'probe_id': qname,
                    'sequence': formatted_seq,
                    'target_transcript': rname,
                    'gene_id': gene_id,
                    'mismatches': nm,
                    'position': pos,
                    'strand': strand,
                    'species': species
                })
                
                if limit and len(alignments) >= limit:
                    break
        
        return alignments
        
    except Exception as e:
        print(f"Error parsing SAM file: {e}")
        return []


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        sam_file = sys.argv[1]
        alignments = parse_sam_file(sam_file, limit=10)
        
        print("\n" + "="*100)
        print(f"{'Probe ID':<15} {'Sequence':<25} {'Target':<20} {'Gene':<15} {'MM':<5} {'Pos':<8} {'Strand':<6}")  
        print("="*100)
        
        for aln in alignments:
            print(f"{aln['probe_id']:<15} {aln['sequence']:<25} {aln['target_transcript']:<20} "
            f"{(aln.get('gene_id') or '-'):<15} {aln['mismatches']:<5} {aln['position']:<8} {aln['strand']:<6}")  

        
        print("="*100)