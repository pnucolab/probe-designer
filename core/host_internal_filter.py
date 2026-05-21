"""
Host-internal probe safety filter.

Used when the input is a host (human/mouse) gene/transcript and the alignment
target is the host transcriptome. A probe is kept only if every alignment
points to the same gene — multiple isoforms of one gene count as a single
target, any cross-gene hit makes the probe unsafe.

Non-aligned probes are absent from the result; the caller drops them.
"""

from collections import defaultdict
from typing import Dict, Optional, Set, Tuple


def _strip_version(transcript_id: str) -> str:
    return transcript_id.split('.', 1)[0]


def filter_probes_by_self_alignment(
    sam_file: str,
    source_transcripts: Optional[Set[str]] = None,
    source_gene: Optional[str] = None,
    gene_mappings: Optional[Dict[str, str]] = None,
) -> Tuple[Set[str], Dict[str, Dict]]:
    """
    Return (safe_probe_ids, per_probe_details).

    A probe is safe iff every aligned RNAME resolves (via `gene_mappings`) to
    the same gene. Different isoforms of the same gene collapse to one group
    and are treated as self. Transcripts missing from `gene_mappings` fall
    back to their version-stripped RNAME as their own group, so unmapped hits
    that don't share an ID are treated as cross-gene.

    `source_transcripts` and `source_gene` are kept in the signature for
    caller compatibility but are not required for the safety decision.
    """
    if gene_mappings is None:
        gene_mappings = {}

    probe_alignments: Dict[str, Set[str]] = defaultdict(set)
    with open(sam_file, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('@'):
                continue
            parts = line.split('\t')
            if len(parts) < 11:
                continue
            rname = parts[2]
            if rname == '*':
                continue
            probe_alignments[parts[0]].add(rname)

    safe: Set[str] = set()
    details: Dict[str, Dict] = {}
    for probe_id, rnames in probe_alignments.items():
        groups: Set[str] = set()
        for rname in rnames:
            base = _strip_version(rname)
            group = gene_mappings.get(rname) or gene_mappings.get(base) or base
            groups.add(group)
        is_safe = len(groups) == 1
        details[probe_id] = {
            'aligned_to': sorted(rnames),
            'genes': sorted(groups),
            'safe': is_safe,
        }
        if is_safe:
            safe.add(probe_id)

    return safe, details
