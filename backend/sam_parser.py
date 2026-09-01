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


def compute_alignment_diffs(cigar: str, md: str) -> tuple:
    """
    Split an alignment's differences into (substitutions, indel_bases).

    Substitutions come from the MD tag; indel bases (insertions/deletions,
    i.e. bulges) come from the CIGAR. razers3's NM:i tag drops some indels
    (e.g. trailing deletions), so both are recomputed from CIGAR + MD.
    """
    md_subs = sum(1 for _, typ, _ in parse_md_tag(md) if typ == 'substitution')
    indel_bases = sum(count for count, op in parse_cigar(cigar) if op in ('I', 'D'))
    return md_subs, indel_bases


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


def reconstruct_reference_from_sam(sequence: str, cigar: str, md: str) -> str:
    """
    Reconstruct the reference (target) sequence for an alignment using SEQ + CIGAR + MD.
    Mismatch bases (from MD) are lowercased to highlight divergence from the probe.
    Deletions (bases absent in probe) are inserted from MD. Insertions (extra probe
    bases) emit '-' so target columns stay aligned with the probe display.
    """
    cigar_ops = parse_cigar(cigar)
    md_mismatches = parse_md_tag(md)

    md_subs = {pos: base for pos, typ, base in md_mismatches if typ == 'substitution'}
    md_dels = {pos: base for pos, typ, base in md_mismatches if typ == 'deletion'}

    result = []
    seq_pos = 0
    ref_pos = 0

    for count, op in cigar_ops:
        if op == 'M':
            for _ in range(count):
                if ref_pos in md_subs:
                    result.append(md_subs[ref_pos].lower())
                elif seq_pos < len(sequence):
                    result.append(sequence[seq_pos].upper())
                seq_pos += 1
                ref_pos += 1
        elif op == 'I':
            for _ in range(count):
                result.append('-')
            seq_pos += count
        elif op == 'D':
            deleted = md_dels.get(ref_pos, '')
            for base in deleted:
                result.append(base.lower())
            ref_pos += count

    return ''.join(result)


def same_genome(probe_id, target):
    """True if reference `target` belongs to the genome named by `probe_id`.

    Mirror of core.probe_classifier.same_genome (kept local to avoid a
    core<->backend import). Used for microbiome probe-input where the probe id
    is the genome id used to name that genome's contigs, so self-alignments can
    be excluded from the off-target views regardless of contig-header style.
    """
    if not probe_id or not target:
        return False
    return target == probe_id or target.startswith(probe_id + '_')


def _is_self_alignment(rname, gene_label, source_gene, source_transcripts, source_transcripts_versionless):
    """Return True if this alignment row points at the source gene/transcript."""
    if source_gene and gene_label and gene_label == source_gene:
        return True
    if rname:
        if rname in source_transcripts:
            return True
        base = rname.split('.', 1)[0]
        if base in source_transcripts_versionless:
            return True
    return False


def count_sam_alignments(sam_path: str, source_gene=None, source_transcripts=None,
                         genome_self_match=False) -> int:
    """
    Count number of unique alignments in SAM file (excluding headers).
    Dedupes by (qname, rname, pos, strand) so primary+secondary records of
    the same hit count once. When `source_gene`/`source_transcripts` is
    provided, alignments that resolve to the source are excluded so the
    count matches the off-target view.
    """
    source_transcripts = set(source_transcripts or [])
    source_transcripts_versionless = {t.split('.', 1)[0] for t in source_transcripts}
    seen = set()
    try:
        with open(sam_path, 'r', encoding='utf-8') as f:
            for line in f:
                if line.startswith('@'):
                    continue
                parts = line.split('\t')
                if len(parts) < 11:
                    continue
                qname = parts[0]
                try:
                    flag = int(parts[1])
                except ValueError:
                    continue
                rname = parts[2]
                # Annotated SAM has gene_id in col 4 (non-numeric); position
                # is then col 5. Standard SAM has position in col 4.
                col4 = parts[3]
                gene_label = None
                try:
                    pos = int(col4)
                except ValueError:
                    gene_label = col4
                    try:
                        pos = int(parts[4])
                    except (ValueError, IndexError):
                        continue
                if genome_self_match and same_genome(qname, rname):
                    continue
                if (source_gene or source_transcripts) and _is_self_alignment(
                    rname, gene_label, source_gene, source_transcripts, source_transcripts_versionless
                ):
                    continue
                strand = '-' if (flag & 0x10) else '+'
                seen.add((qname, rname, pos, strand))
        return len(seen)
    except Exception as e:
        print(f"Error counting alignments: {e}")
        return 0


def parse_sam_file(sam_path: str, offset: int = 0, limit: int = None, mismatch_filter: int = None, probe_id_filter: str = None, source_gene=None, source_transcripts=None, genome_self_match=False) -> List[Dict[str, Any]]:
    """
    Parse SAM file and extract alignment information with pagination support.
    Handles * sequences by tracking the last valid sequence per probe.
    Dedupes records that share (qname, rname, pos, strand) — some alignment
    runs emit both primary and secondary records for the same hit. When
    `source_gene`/`source_transcripts` is provided, alignments to the source
    are skipped so callers see off-target hits only.
    """
    alignments = []
    skipped = 0
    sequence_cache = {}
    seen_keys = set()
    source_transcripts_set = set(source_transcripts or [])
    source_transcripts_versionless = {t.split('.', 1)[0] for t in source_transcripts_set}
    self_filter_active = bool(source_gene or source_transcripts_set)
    if probe_id_filter:
        base_probe_id = probe_id_filter.split('|')[0].strip()
        # probe_numeric = re.search(r'probe_(\d+)', probe_id_filter)
        # probe_num = probe_numeric.group(1) if probe_numeric else None
    try:
        with open(sam_path, 'r', encoding='utf-8') as f:
            for line in f:
                if line.startswith('@'):
                    continue

                parts = line.strip().split('\t')
                if len(parts) < 11:
                    continue

                qname = parts[0]
                if probe_id_filter:
                    qname_base = qname.split('|')[0].strip()
                    if qname.strip() != probe_id_filter and qname_base != base_probe_id:
                        continue

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

                if genome_self_match and same_genome(qname, rname):
                    continue
                if self_filter_active and _is_self_alignment(
                    rname, gene_id, source_gene, source_transcripts_set, source_transcripts_versionless
                ):
                    continue

                strand_key = '-' if (flag & 0x10) else '+'
                dedup_key = (qname, rname, pos, strand_key)
                if dedup_key in seen_keys:
                    continue
                seen_keys.add(dedup_key)

                if skipped < offset:
                    skipped += 1
                    continue
                
                md = None
                species = None
                for field in parts[11:]:
                    if field.startswith('MD:Z:'):
                        md = field.split(':')[2]
                    elif field.startswith('SP:Z:'):
                        species = field.split(':')[2]

                if md is None:
                    continue

                substitutions, bulges = compute_alignment_diffs(cigar, md)
                nm = substitutions + bulges

                if mismatch_filter is not None and nm != mismatch_filter:
                    continue

                is_reverse = bool(flag & 0x10)
                strand = '-' if is_reverse else '+'
                
                formatted_seq = apply_mismatches_to_sequence(seq, cigar, md, flag)
                target_seq = reconstruct_reference_from_sam(seq, cigar, md)

                alignments.append({
                    'probe_id': qname,
                    'sequence': formatted_seq,
                    'target_sequence': target_seq,
                    'target_transcript': rname,
                    'gene_id': gene_id,
                    'mismatches': nm,
                    'substitutions': substitutions,
                    'bulges': bulges,
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
        sample_alignments = parse_sam_file(sam_file, limit=10)
        
        print("\n" + "="*100)
        print(f"{'Probe ID':<15} {'Sequence':<25} {'Target':<20} {'Gene':<15} {'MM':<5} {'Pos':<8} {'Strand':<6}")  
        print("="*100)
        
        for aln in sample_alignments:
            print(f"{aln['probe_id']:<15} {aln['sequence']:<25} {aln['target_transcript']:<20} "
            f"{(aln.get('gene_id') or '-'):<15} {aln['mismatches']:<5} {aln['position']:<8} {aln['strand']:<6}")  

        
        print("="*100)