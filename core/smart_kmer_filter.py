"""
K-mer safety filtering.
Any k-mer match to a gene other than the source marks the probe as unsafe.
"""

import time


def smart_kmer_filter(match_details, mer_to_probe, probes, k, species, output_base, chrom_files, razers, args, source_transcripts=None):
    """
    K-mer filtering: matches to other genes mark the probe as unsafe.
    Matches to the source gene's own transcripts are excluded.
    """

    start_time = time.time()
    unsafe_probes = set()
    probe_match_details = {}

    if source_transcripts is None:
        source_transcripts = set()

    self_matches = 0
    other_matches = 0

    for mer_id, transcript_id, _ in match_details:
        # Skip matches to the source gene's own transcripts
        base_transcript = transcript_id.split()[0] if transcript_id else transcript_id
        if base_transcript in source_transcripts:
            self_matches += 1
            continue

        other_matches += 1
        if mer_id in mer_to_probe:
            probe_list = mer_to_probe[mer_id]
            for probe_id, _ in probe_list:
                unsafe_probes.add(probe_id)
                if probe_id not in probe_match_details:
                    probe_match_details[probe_id] = {'total_matches': 0}
                probe_match_details[probe_id]['total_matches'] += 1

    elapsed = time.time() - start_time
    print(f"K-mer filtering took {elapsed:.2f}s")
    if source_transcripts:
        print(f"  Self-gene matches (excluded): {self_matches:,}")
        print(f"  Other-gene matches: {other_matches:,}")
    print(f"  Unsafe probes (off-target k-mer matches): {len(unsafe_probes)}")

    return unsafe_probes, probe_match_details, set()
