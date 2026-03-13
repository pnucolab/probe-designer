"""
Determines if probes are safe based on alignment patterns.
"""

import re
from collections import defaultdict
from typing import Dict, Optional, Set


def infer_source_transcripts(sam_file, gene_mappings=None, is_microbiome=False):
    """
    Infer source transcript(s) and gene name from probe alignment patterns.

    Uses unique probes with perfect matches (NM:i:0) to identify the source.
    With mismatches allowed, probes hit many transcripts — but only the
    source transcript will have many probes mapping perfectly (NM:i:0).

    Args:
        sam_file: Filtered SAM file with probe alignments
        gene_mappings: dict of transcript_id -> gene_name (from samannotator)

    Returns:
        dict with:
            'source_transcripts': set of transcript IDs identified as source
            'source_gene': gene name (resolved via mappings) or None
    """
    nm_pattern = re.compile(r'NM:i:(\d+)')

    perfect_probes = defaultdict(set)
    all_probes = defaultdict(set)

    with open(sam_file, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('@'):
                continue
            parts = line.split('\t')
            if len(parts) < 11:
                continue
            probe_id = parts[0]
            transcript_id = parts[2]
            if transcript_id == '*':
                continue

            all_probes[transcript_id].add(probe_id)

            m = nm_pattern.search(line)
            if m and int(m.group(1)) == 0:
                perfect_probes[transcript_id].add(probe_id)

    if not all_probes:
        return {'source_transcripts': set(), 'source_gene': None}

    if perfect_probes:
        ranked = sorted(perfect_probes.items(), key=lambda x: len(x[1]), reverse=True)
        top_transcript = ranked[0][0]
        top_count = len(ranked[0][1])
        print(f"Source inference: {top_transcript} has {top_count} unique probes with perfect matches")
    else:
        ranked = sorted(all_probes.items(), key=lambda x: len(x[1]), reverse=True)
        top_transcript = ranked[0][0]
        top_count = len(ranked[0][1])
        print(f"Source inference (no perfect matches): {top_transcript} has {top_count} unique probes")

    source_gene = None
    source_transcripts = set()

    if gene_mappings:
        source_gene = gene_mappings.get(top_transcript)
        if source_gene:
            for tid, gname in gene_mappings.items():
                if gname == source_gene:
                    source_transcripts.add(tid)
        else:
            source_transcripts.add(top_transcript)
    else:
        # No gene mappings — detect isoforms by probe overlap.
        # Isoforms of the same gene share most perfectly-mapping probes.
        source_transcripts.add(top_transcript)
        if not is_microbiome:
            probe_set = perfect_probes if perfect_probes else all_probes
            top_probes = probe_set[top_transcript]
            if top_probes:
                for tid, probes in probe_set.items():
                    if tid == top_transcript:
                        continue
                    overlap = len(probes & top_probes)
                    # If >50% of this transcript's probes also map to the top transcript,
                    # it's likely an isoform of the same gene
                    if overlap > 0 and overlap / len(probes) > 0.5:
                        source_transcripts.add(tid)

    if source_gene:
        print(f"Source gene: {source_gene} ({len(source_transcripts)} transcript(s))")
    elif len(source_transcripts) > 1:
        print(f"Source transcript: {top_transcript} + {len(source_transcripts) - 1} isoform(s) detected by probe overlap")
    else:
        print(f"Source transcript: {top_transcript} (no gene mapping available)")

    return {'source_transcripts': source_transcripts, 'source_gene': source_gene}


def classify_probes_from_sam(sam_file: str, is_microbiome: bool = False,
                             source_transcripts: Optional[Set[str]] = None,
                             source_gene: Optional[str] = None) -> Dict:
    """
    Classify probes as safe or risky based on SAM alignments.

    Self-alignment is determined by comparing alignment targets against
    known source transcripts/gene. Source identity is inferred from
    probe alignment patterns (see infer_source_transcripts in probe_designer.py).

    Args:
        sam_file: Path to SAM file (annotated SAMs have gene name in column 4)
        is_microbiome: Whether this is a microbiome species
        source_transcripts: Set of transcript IDs identified as the source
        source_gene: Gene name of the source (resolved from mappings)
    """
    probe_data = {}
    sp_pattern = re.compile(r'SP:Z:(\S+)')
    nm_pattern = re.compile(r'NM:i:(\d+)')

    if source_transcripts is None:
        source_transcripts = set()

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
                    'targets': [],
                    'target_transcripts': []
                }
            probe_data[probe_id]['targets'].append(target_info)
            probe_data[probe_id]['target_transcripts'].append(target_info)
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
            is_offtarget = _check_offtarget(
                probe_id, data['targets'], data['target_transcripts'],
                source_transcripts, source_gene
            )

            if is_offtarget:
                if data['min_mm'] <= 1:
                    status = 'high_risk'
                else:
                    status = 'medium_risk'
            else:
                status = 'safe'

        # Primary check: if source_transcripts known, use them directly
        elif source_transcripts and all(t in source_transcripts for t in data['target_transcripts']):
            status = 'safe'

        # Check if aligns to only one gene and it's self
        elif unique_gene_count == 1:
            is_self = _check_is_self(
                probe_id, data['targets'], data['target_transcripts'],
                data['unique_genes'], source_transcripts, source_gene
            )

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


def _check_is_self(probe_id, targets, target_transcripts, unique_genes,
                   source_transcripts, source_gene):
    """
    Check if a probe's alignments are self-alignments.

    Uses three strategies in order of reliability:
    1. Source transcript match (from infer_source_transcripts)
    2. Source gene name match (from annotation)
    3. Legacy probe ID substring match (fallback)
    """
    # Strategy 1: check if all target transcripts are in the inferred source set
    if source_transcripts:
        return all(t in source_transcripts for t in target_transcripts)

    # Strategy 2: check if the only gene hit matches the source gene
    if source_gene and unique_genes:
        only_gene = next(iter(unique_genes))
        return only_gene == source_gene

    # Strategy 3: legacy fallback — substring match from probe ID
    probe_source = probe_id.split('|')[-1].replace('transcript:', '')
    return any(probe_source in target for target in targets)


def _check_offtarget(probe_id, targets, target_transcripts,
                     source_transcripts, source_gene):
    """
    Check if ANY target is off-target (for microbiome mode).

    Uses source transcripts/gene when available, falls back to probe ID substring.
    """
    # Strategy 1: check against inferred source transcripts
    if source_transcripts:
        return any(t not in source_transcripts for t in target_transcripts)

    # Strategy 2: legacy fallback
    probe_source = probe_id.split('|')[-1].replace('transcript:', '')
    for target in targets:
        if probe_source not in target:
            return True
    return False
