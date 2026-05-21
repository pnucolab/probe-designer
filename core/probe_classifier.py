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


_HEADER_TOKEN_SPLIT = re.compile(r'[\s|,;:()\[\]]+')

_REFSEQ_LOOKUP: dict = {}
_REFSEQ_LOADED = False


def _load_refseq_lookup() -> dict:
    """Load NCBI RefSeq accession -> gene_symbol map (versionless keys).

    Read once from data/ncbi/refseq_to_gene.tsv (built by
    data/ncbi/build_refseq_lookup.py). Returns {} if the file is missing.
    """
    global _REFSEQ_LOADED
    if _REFSEQ_LOADED:
        return _REFSEQ_LOOKUP
    import os
    here = os.path.dirname(os.path.abspath(__file__))
    root = os.path.dirname(here)
    path = os.path.join(root, "data", "ncbi", "refseq_to_gene.tsv")
    if os.path.exists(path):
        try:
            with open(path, "r", encoding="utf-8") as fh:
                next(fh, None)
                for line in fh:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) >= 2 and parts[0] and parts[1]:
                        _REFSEQ_LOOKUP.setdefault(parts[0], parts[1])
        except OSError:
            pass
    _REFSEQ_LOADED = True
    return _REFSEQ_LOOKUP


def resolve_source_from_header(header, gene_mappings):
    """
    Resolve the source gene/transcripts from a user-supplied FASTA header.

    Scans every token in the header (split on whitespace and the punctuation
    used in typical FASTA descriptions: `|,;:()[]`) and accepts the first
    token that resolves to either:
      - a transcript ID present as a key in `gene_mappings`, or
      - a gene symbol present as a value in `gene_mappings`.

    Each token is also tried with a trailing `.\\d+` version stripped, and
    underscore-separated parts are tried individually. This catches identifiers
    anywhere in the description, e.g. `EGFR` inside `(EGFR)`.

    Returns the same dict shape as `infer_source_transcripts`. If nothing
    resolves, returns empty source -- caller should treat every alignment as
    off-target.
    """
    empty = {'source_transcripts': set(), 'source_gene': None}
    if not header or not gene_mappings:
        return empty

    raw = header.strip()
    if not raw:
        return empty

    gene_symbols = set(gene_mappings.values())
    refseq_lookup = _load_refseq_lookup()

    source_gene = None
    matched_token = None
    for token in _HEADER_TOKEN_SPLIT.split(raw):
        if not token:
            continue
        candidates = [token, re.sub(r'\.\d+$', '', token)]
        if '_' in token:
            candidates.extend(p for p in token.split('_') if p)
        for cand in candidates:
            if cand in gene_mappings:
                source_gene = gene_mappings[cand]
                matched_token = cand
                break
            if cand in gene_symbols:
                source_gene = cand
                matched_token = cand
                break
            if cand in refseq_lookup:
                # NCBI RefSeq accession -> gene symbol; only accept if that
                # symbol is in the host transcriptome (i.e. has Ensembl
                # transcripts) so source_transcripts can be populated.
                sym = refseq_lookup[cand]
                if sym in gene_symbols:
                    source_gene = sym
                    matched_token = cand
                    break
        if source_gene:
            break

    if not source_gene:
        first_token = raw.split()[0]
        print(f"Source header token '{first_token}' not found in transcriptome mappings - "
              f"skipping source inference (all alignments treated as off-target)")
        return empty

    source_transcripts = {tid for tid, gname in gene_mappings.items() if gname == source_gene}
    print(f"Source from header: gene '{source_gene}' "
          f"(token '{matched_token}', {len(source_transcripts)} transcript(s))")
    return {'source_transcripts': source_transcripts, 'source_gene': source_gene}


def classify_probes_from_sam(sam_file: str, is_microbiome: bool = False,
                             source_transcripts: Optional[Set[str]] = None,
                             source_gene: Optional[str] = None,
                             is_probe_input: bool = False) -> Dict:
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
                    'target_transcripts': [],
                    'alignments': [],  # (transcript, gene_name, mismatches) per hit
                }
            probe_data[probe_id]['targets'].append(target_info)
            probe_data[probe_id]['target_transcripts'].append(target_info)
            probe_data[probe_id]['alignments'].append((target_info, gene_name, mismatches))
            probe_data[probe_id]['total'] += 1
            probe_data[probe_id]['min_mm'] = min(probe_data[probe_id]['min_mm'], mismatches)
            probe_data[probe_id]['unique_genes'].add(gene_name)

            if species_value and species_value != 'Unknown':
                probe_data[probe_id]['is_microbiome_species'] = True

    results = {}
    for probe_id, data in probe_data.items():
        unique_gene_count = len(data['unique_genes'])

        # Minimum mismatch among OFF-TARGET alignments only (i.e. hits that are
        # NOT to the source gene / its isoforms). The risk level and the
        # reported mismatch count must reflect the off-target, not a perfect
        # self-match. Without this, a probe that matches its own gene at 0 mm
        # but a paralog at 2 mm would be mislabeled "high risk: 0 mismatches".
        off_mm = _offtarget_min_mm(
            data['alignments'], source_transcripts, source_gene,
            treat_all_as_offtarget=is_probe_input,
        )

        # Probe-input mode: source is unknown, so any alignment is off-target.
        # Skip self-detection entirely to avoid the probe-ID substring fallback.
        if is_probe_input:
            status = 'high_risk' if (off_mm is not None and off_mm <= 1) else 'medium_risk'

        # For microbiome: check if ANY alignment is to a different species
        elif is_microbiome and data['is_microbiome_species']:
            is_offtarget = _check_offtarget(
                probe_id, data['targets'], data['target_transcripts'],
                source_transcripts, source_gene
            )

            if is_offtarget:
                status = 'high_risk' if (off_mm is not None and off_mm <= 1) else 'medium_risk'
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
            elif off_mm is not None and off_mm <= 1:
                status = 'high_risk'
            else:
                status = 'medium_risk'
        # Multiple genes hit, at least one off-target.
        elif off_mm is not None and off_mm <= 1:
            status = 'high_risk'
        else:
            status = 'medium_risk'

        # For non-safe probes report the off-target mismatch count; for safe
        # probes the mismatch count isn't shown so the global min is fine.
        if status == 'safe':
            reported_mm = data['min_mm'] if data['min_mm'] != float('inf') else 0
        else:
            reported_mm = off_mm if off_mm is not None else (
                data['min_mm'] if data['min_mm'] != float('inf') else 0
            )

        results[probe_id] = {
            'status': status,
            'total_alignments': data['total'],
            'unique_genes': unique_gene_count,
            'min_mismatches': reported_mm,
        }

    return results


def _offtarget_min_mm(alignments, source_transcripts, source_gene,
                      treat_all_as_offtarget=False):
    """Minimum mismatch count among off-target alignments.

    An alignment is off-target if its transcript is not in
    `source_transcripts` (version-tolerant) and its gene != `source_gene`.
    Returns None when there are no off-target alignments.
    `treat_all_as_offtarget` (probe-input mode) ignores the source entirely.
    """
    src = set(source_transcripts or [])
    src_versionless = {t.split('.', 1)[0] for t in src}
    offs = []
    for transcript, gene, mm in alignments:
        if treat_all_as_offtarget:
            offs.append(mm)
            continue
        base = transcript.split('.', 1)[0]
        in_source = (transcript in src) or (base in src_versionless)
        if not in_source and source_gene:
            in_source = (gene == source_gene)
        if not in_source:
            offs.append(mm)
    return min(offs) if offs else None


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
