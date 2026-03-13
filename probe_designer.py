"""
Simple Probe Alignment Pipeline - Modified

This script handles two input scenarios:
1. Gene sequence (pasted) -> generate probes from gene -> align to species genome
2. Probe sequence (pasted) -> skip generation -> align to species genome

Usage:
  # Scenario 1: Gene sequence (pasted)
  python simple_pipeline.py --gene-sequence ">gene1\nATCG..." --species human --probe-length 36 --max-mismatches 2
  
  # Scenario 2: Probe sequence (pasted)
  python simple_pipeline.py --probe-sequence ">probe1\nATCG..." --species human --max-mismatches 2

Notes:
- Only ONE input type allowed per run (mutually exclusive)
- Gene/probe sequences must be in proper FASTA format
- Probe length auto-detected when probes provided; otherwise uses --probe-length
"""

import argparse
import os
import subprocess
import time
import re
import concurrent.futures
from multiprocessing import cpu_count
import sys
from datetime import datetime
import tempfile
import traceback
from pathlib import Path
import shutil
from Bio import SeqIO
from Bio.Seq import Seq
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'core'))
from score_probes import score_and_save_probes
from scorer import ThermodynamicProbeScorer
from probe_classifier import classify_probes_from_sam, infer_source_transcripts
from smart_kmer_filter import smart_kmer_filter


def get_local_release_number(filename_pattern, local_dir):
    """Extract highest release number from local files."""
    os.makedirs(local_dir, exist_ok=True)
    pattern_gz = re.compile(filename_pattern)
    pattern_no_gz = re.compile(filename_pattern.replace(r'\.gz', ''))
    max_release = None
    for filename in os.listdir(local_dir):
        file_path = os.path.join(local_dir, filename)
        if filename.endswith('.fai'):
            os.remove(file_path)
            continue
        if match := pattern_gz.match(filename):
            release = int(match.group(1))
            max_release = release if max_release is None else max(max_release, release)
        elif match := pattern_no_gz.match(filename):
            release = int(match.group(1))
            max_release = release if max_release is None else max(max_release, release)
    return max_release


def validate_local_gencode(species, base_dir=None):
    """Locate and validate pre-existing local Gencode data (genome, chunks)."""
    if base_dir is None:
        base_dir = os.path.dirname(os.path.abspath(__file__))
    species_data_dir = os.path.join(base_dir, 'data', 'gencode_data', species)
    if species == "human":
        genome_pattern = r'GRCh38\.p14\.genome\.release_(\d+)\.fa'
    else:
        genome_pattern = r'GRCm39\.genome\.release_(\d+)\.fa'
    local_release = get_local_release_number(genome_pattern, species_data_dir)
    if not local_release:
        print(f"Error: No local genome file found for {species} in {species_data_dir}")
        return False
    if species == "human":
        genome_file = f"{species_data_dir}/GRCh38.p14.genome.release_{local_release}.fa"
    else:
        genome_file = f"{species_data_dir}/GRCm39.genome.release_{local_release}.fa"
    chunk_dir = f"{species_data_dir}/transcript_chunks"
    if not os.path.exists(genome_file):
        print(f"Error: Genome file not found: {genome_file}")
        return False
    if not os.path.exists(chunk_dir):
        print(f"Error: Transcript chunks directory not found: {chunk_dir}")
        return False
    existing_chunks = [f for f in os.listdir(chunk_dir)
                      if f.endswith('.fa') and f'release_{local_release}' in f]
    if not existing_chunks:
        print(f"Error: No chunk files found for release_{local_release} in {chunk_dir}")
        return False
    print(f"Using local {species} data (release_{local_release}): {len(existing_chunks)} chunk files")
    return True

DEFAULT_DIRS = {
    'gene': 'gene_sequences',
    'probe': 'probe_sequences'
}
class samannotator:
    """SAM annotator """
    
    def __init__(self, gtf_dir='data/gencode_data'):
        self.gtf_dir = Path(gtf_dir)
        self.mappings = {}
        
    def _find_gtf_file(self, species):
        """Find GTF file for species and extract release number."""
        species_dir = self.gtf_dir / species
        if not species_dir.exists():
            return None, None
        
        for file in species_dir.glob('*.gtf'):
            match = re.search(r'release_(\d+)', file.name)
            release = int(match.group(1)) if match else None
            return file, release
        return None, None
    
    def _get_cache_file(self, species):
        """Get cache file path for species - stored in species directory."""
        species_cache_dir = self.gtf_dir / species / 'annotation_data'
        species_cache_dir.mkdir(parents=True, exist_ok=True)
        return species_cache_dir / f'{species}_transcript_gene_mappings.txt'
    
    def _get_cache_release(self, cache_file):
        """Get release number from cache file."""
        if not cache_file.exists():
            return None
        
        try:
            with open(cache_file, 'r', encoding='utf-8') as f:
                first_line = f.readline().strip()
                if first_line.startswith('#') and 'release_' in first_line:
                    match = re.search(r'release_(\d+)', first_line)
                    return int(match.group(1)) if match else None
        except Exception:
            pass
        return None
    
    def _load_cache(self, cache_file):
        """Load mappings from cache file (lightning fast)."""
        if not cache_file.exists():
            return False
        
        start = time.time()
        with open(cache_file, 'r', encoding='utf-8') as f:
            next(f)  
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) == 2:
                    self.mappings[parts[0]] = parts[1]
        
        print(f"Loaded {len(self.mappings):,} transcript→gene mappings ({time.time()-start:.3f}s)")
        return True
    
    def _extract_from_gtf(self, gtf_file, release, cache_file):
        """Extract mappings from GTF file (comprehensive but slower)."""
        print(f"Extracting transcript→gene mappings from GTF (release_{release})...")
        start = time.time()
        
        transcript_pattern = re.compile(r'transcript_id "([^"]+)"')
        gene_pattern = re.compile(r'gene_name "([^"]+)"')
        
        with open(gtf_file, 'r', encoding='utf-8') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                
                if 'transcript_id' in line and 'gene_name' in line:
                    t_match = transcript_pattern.search(line)
                    g_match = gene_pattern.search(line)
                    
                    if t_match and g_match:
                        self.mappings[t_match.group(1)] = g_match.group(1)
        
        elapsed = time.time() - start
        print(f"Extracted {len(self.mappings):,} mappings in {elapsed:.2f}s")
        
        with open(cache_file, 'w', encoding='utf-8') as f:
            f.write(f"# Transcript to gene mappings from release_{release}\n")
            for tid, gname in self.mappings.items():
                f.write(f"{tid}\t{gname}\n")
        print(f"Cached mappings to {cache_file}")
    
    def _ensure_mappings_loaded(self, species):
        """Smart loading: use cache if valid, otherwise extract from GTF."""
        gtf_file, gtf_release = self._find_gtf_file(species)
        
        if not gtf_file:
            print(f"No GTF file found for {species} - annotation will be skipped")
            return False
        
        cache_file = self._get_cache_file(species)
        cache_release = self._get_cache_release(cache_file)
    
        if cache_release == gtf_release and self._load_cache(cache_file):
            print(f"Using cached mappings (release_{cache_release})")
            return True
    
        if cache_release and cache_release != gtf_release:
            print(f"GTF updated: release_{cache_release} → release_{gtf_release}")
        
        self._extract_from_gtf(gtf_file, gtf_release, cache_file)
        return True
    
    def annotate_sam(self, input_sam, output_sam, species):
        """Annotate SAM file with gene names."""
        print(f"\n{'='*80}")
        print("ANNOTATING SAM WITH GENE NAMES")
        print(f"{'='*80}")
    
        if not self.mappings:
            if not self._ensure_mappings_loaded(species):
                print("Skipping annotation - using original SAM")
                return input_sam
        
        print(f"Annotating: {os.path.basename(input_sam)} → {os.path.basename(output_sam)}")
        start = time.time()
        
        processed = 0
        annotated = 0
        headers = 0
        
        with open(input_sam, 'r', encoding='utf-8') as infile, open(output_sam, 'w', encoding='utf-8') as outfile:
            for line in infile:
                processed += 1
                
                if line.startswith('@'):
                    outfile.write(line)
                    headers += 1
                else:
                    fields = line.rstrip('\n').split('\t')
                    if len(fields) >= 3:
                        transcript_id = fields[2]
                        gene_name = self.mappings.get(transcript_id, transcript_id)
                        fields.insert(3, gene_name)
                        annotated += 1
                    
                    outfile.write('\t'.join(fields) + '\n')
        
        elapsed = time.time() - start
        print("Annotation complete!")
        print(f"Lines processed: {processed:,}")
        print(f"Headers: {headers:,}")
        print(f"Alignments annotated: {annotated:,}")
        print(f"Time: {elapsed:.3f}s ({processed/elapsed:,.0f} lines/sec)")
        print(f"Output: {output_sam}")
        
        return output_sam


def count_and_filter_gc(fasta_file, min_gc=40, max_gc=80):
    """Count probes and filter by GC content."""
    count = 0
    passing = []
    
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = str(record.seq).upper()
        if 'N' in seq:
            count += 1
            continue
        gc = (seq.count('G') + seq.count('C')) / len(seq) * 100
        if min_gc <= gc <= max_gc:
            passing.append(record)
        count += 1
    
    return count, passing

def validate_fasta_format(filepath):
    """Validate FASTA format and DNA alphabet. Returns (valid, error_message, sequences)"""
    sequences = []
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            current_header = None
            current_seq = []
            line_num = 0
            
            for line in f:
                line_num += 1
                line = line.strip()
                if not line:
                    continue
                    
                if line.startswith('>'):
                    if current_header is not None:
                        seq_str = ''.join(current_seq)
                        if not seq_str:
                            return False, f"Empty sequence for header: {current_header}", []
                        sequences.append((current_header, seq_str))
                    current_header = line[1:].strip()
                    if not current_header:
                        return False, f"Empty header at line {line_num}", []
                    current_seq = []
                else:
                    if not re.match(r'^[ATCGNatcgn]+$', line):
                        return False, f"Invalid DNA characters at line {line_num}. Only ATCGN allowed.", []
                    current_seq.append(line.upper())
            
            if current_header is not None:
                seq_str = ''.join(current_seq)
                if not seq_str:
                    return False, f"Empty sequence for header: {current_header}", []
                sequences.append((current_header, seq_str))
            
            if not sequences:
                return False, "No sequences found in FASTA file", []
                
            return True, None, sequences
            
    except (FileNotFoundError, PermissionError, UnicodeDecodeError, OSError) as e:
        return False, f"File error: {e}", []


def resolve_file_path(filename, file_type):
    """Resolve file path: check if full path exists, otherwise search in default directory"""
    if os.path.exists(filename):
        return filename
    
    default_dir = DEFAULT_DIRS.get(file_type, '.')
    default_path = os.path.join(default_dir, filename)
    
    if os.path.exists(default_path):
        return default_path
    
    return None


def detect_probe_length(sequences):
    """Detect probe length from provided sequences. Returns length or None if inconsistent."""
    lengths = set()
    for _, seq in sequences:
        lengths.add(len(seq))
    
    if len(lengths) == 1:
        return lengths.pop()
    else:
        return None


def run(cmd, description, check_output_file=None):
    """Run a command with logging and timing."""
    print(f"\n{description}")
    print("-" * 60)
    print(f"Command: {' '.join(cmd) if isinstance(cmd, list) else cmd}")
    start = time.time()
    result = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if result.stdout:
        print(result.stdout.strip())
    if result.returncode != 0:
        print(f"{description} failed (exit {result.returncode})")
        if result.stderr:
            print(result.stderr)
        return False
    if check_output_file and not os.path.exists(check_output_file):
        print(f"Warning: expected output file not found: {check_output_file}")
        return False
    dur = time.time() - start
    print(f"{description} completed in {dur:.1f}s")
    return True

def annotate_microbiome_sam(sam_file, species, _, align_host=False):
    """Add species annotations to SAM file for microbiome data."""
    metadata_file = os.path.join('data', species, 'genomes-all_metadata.tsv')
    
    if not os.path.exists(metadata_file):
        print(f"Warning: Metadata file not found: {metadata_file}")
        return False
    
    print(f"\nAnnotating SAM file with species information from {metadata_file}...")
    
    genome_to_species = {}
    try:
        with open(metadata_file, 'r', encoding='utf-8') as f:
            for i, line in enumerate(f):
                if i == 0:
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 15:
                    genome = parts[0]
                    lineage = parts[14]
                    match = re.search(r's__([^;]+)', lineage)
                    if match:
                        species_name = match.group(1)
                        genome_to_species[genome] = species_name
    except Exception as e:
        print(f"Error reading metadata: {e}")
        return False
    
    print(f"Loaded {len(genome_to_species)} genome-to-species mappings")
    
    annotated_sam = sam_file.replace('.sam', '_annotated.sam')
    
    try:
        transcript_to_gene = {}
        if align_host:
            host_species = 'human' if species in ['gut-microbe', 'human-oral-microbiome', 'human-skin-microbiome', 'human-vaginal-microbiome'] else 'mouse'
            ann = samannotator()
            ann._ensure_mappings_loaded(host_species)
            transcript_to_gene = ann.mappings

        with open(sam_file, 'r', encoding='utf-8') as inf, \
            open(annotated_sam, 'w', encoding='utf-8') as outf:
            
            for line in inf:
                if line.startswith('@'):
                    outf.write(line)
                else:
                    parts = line.strip().split('\t')
                    if len(parts) >= 11:
                        rname = parts[2]
                        
                        if rname.startswith('MGYG'):
                            genome = rname.split('_')[0] if '_' in rname else rname
                            species_name = genome_to_species.get(genome, 'Unknown')
                        elif rname.startswith('ENST') or rname.startswith('ENSMUST'):
                            transcript_base = rname.split('.')[0]
                            species_name = transcript_to_gene.get(transcript_base, transcript_to_gene.get(rname, 'Unknown'))
                        else:
                            species_name = 'Unknown'
                        
                        outf.write(line.strip() + f'\tSP:Z:{species_name}\n')
                    else:
                        outf.write(line)
        
        print(f"Created annotated SAM: {annotated_sam}")
        return True
        
    except Exception as e:
        print(f"Error annotating SAM file: {e}")
        return False
def generate_14mers(sequence, k=14):
    """Generate all possible k-mers from a sequence."""
    return [sequence[i:i+k] for i in range(len(sequence) - k + 1)]


def get_jellyfish_db(genome_files, k, species, threads=1):
    """Get or build a Jellyfish k-mer count database for a species.
    Database is cached in the species data directory for reuse across runs.
    Uses -L 2 to only store k-mers appearing 2+ times, reducing DB size.
    """
    jellyfish = 'jellyfish'
    # Cache DB alongside genome data for reuse
    species_dir = os.path.dirname(genome_files[0])
    db_path = os.path.join(species_dir, f'kmer_{k}mer_counts.jf')

    if os.path.exists(db_path):
        print(f"  Using pre-built Jellyfish database: {db_path}")
        return db_path

    print(f"  Building Jellyfish database (one-time, will be cached for future runs)...")
    file_list = ' '.join(f"'{f}'" for f in genome_files)
    cmd = f"cat {file_list} | {jellyfish} count -m {k} -s 5G -t {threads} -L 2 -C -o '{db_path}' /dev/fd/0"
    proc = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=False)
    if proc.returncode != 0:
        print(f"Warning: Jellyfish count failed: {proc.stderr}")
        return None
    return db_path


def query_jellyfish_db(db_path, kmer_fasta):
    """Query k-mer counts from a Jellyfish database. Returns dict of {original_kmer: count}.
    Jellyfish outputs canonical forms, so we read the input FASTA to map back to original k-mers.
    """
    jellyfish = 'jellyfish'
    # Read original k-mer sequences in order
    original_kmers = []
    with open(kmer_fasta, 'r', encoding='utf-8') as f:
        for line in f:
            if not line.startswith('>'):
                original_kmers.append(line.strip().upper())

    cmd = [jellyfish, 'query', '-s', kmer_fasta, db_path]
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    counts = {}
    if proc.returncode == 0:
        lines = [l for l in proc.stdout.strip().split('\n') if l]
        for i, line in enumerate(lines):
            parts = line.split()
            if len(parts) == 2 and i < len(original_kmers):
                counts[original_kmers[i]] = int(parts[1])
    return counts


def count_source_kmers(source_sequences, k):
    """Count canonical k-mer occurrences in source sequences (matching Jellyfish -C behavior).
    Each k-mer on either strand increments the count of its canonical form.
    Returns counts keyed by BOTH the k-mer and its reverse complement for easy lookup.
    """
    canonical_counts = {}
    for seq in source_sequences:
        seq_upper = seq.upper()
        for i in range(len(seq_upper) - k + 1):
            kmer = seq_upper[i:i+k]
            rc = str(Seq(kmer).reverse_complement())
            canonical = min(kmer, rc)
            canonical_counts[canonical] = canonical_counts.get(canonical, 0) + 1

    # Build lookup that maps both forward and RC to the canonical count
    counts = {}
    for canonical, count in canonical_counts.items():
        rc = str(Seq(canonical).reverse_complement())
        counts[canonical] = count
        counts[rc] = count
    return counts


def find_source_genome_sequences(chrom_files, source_transcripts):
    """Find all sequences from the source genome in the database chunk files.
    E.g., if source transcript is MGYG000000001_1, find all MGYG000000001_* contigs."""
    source_prefixes = set()
    for tid in source_transcripts:
        parts = tid.rsplit('_', 1)
        if len(parts) == 2 and parts[1].isdigit():
            source_prefixes.add(parts[0] + '_')

    source_seqs = []
    for chrom_file in chrom_files:
        for rec in SeqIO.parse(chrom_file, "fasta"):
            if any(rec.id.startswith(prefix) for prefix in source_prefixes):
                source_seqs.append(str(rec.seq))
    return source_seqs


def check_14mer_safety(probes_file, species, output_base, args, source_transcripts=None, input_fasta=None):
    """Check if k-mers from on-target probes have off-target genome matches using Jellyfish."""
    k = args.kmer_length
    print(f"\nChecking {k}-mer safety for on-target probes...")
    probes = list(SeqIO.parse(probes_file, "fasta"))
    if not probes:
        return probes_file

    # Generate k-mers per probe
    probe_to_kmers = {}
    unique_kmers = set()
    for record in probes:
        kmers = generate_14mers(str(record.seq), k)
        probe_to_kmers[record.id] = kmers
        unique_kmers.update(kmers)

    total_mers = sum(len(v) for v in probe_to_kmers.values())
    print(f"Generated {len(unique_kmers)} unique k-mers from {len(probes)} probes")
    print(f"  Total k-mers: {total_mers}, Unique: {len(unique_kmers)}, Duplicates removed: {total_mers - len(unique_kmers)} ({(total_mers-len(unique_kmers))/total_mers*100:.1f}%)")

    microbiome_species_list = ["gut-microbe", "human-oral-microbiome", "human-skin-microbiome",
                            "human-vaginal-microbiome", "mouse-gut-microbiome"]

    chrom_files = []

    if species in microbiome_species_list:
        if args.align_microbiome:
            microbiome_dir = os.path.join('data', species)
            if os.path.exists(microbiome_dir):
                chrom_files = [os.path.join(microbiome_dir, f) for f in os.listdir(microbiome_dir)
                            if f.endswith('.fa') or f.endswith('.fna')]
        elif args.align_host:
            if species == "mouse-gut-microbiome":
                host_dir = os.path.join('data', 'gencode_data', 'mouse', 'transcript_chunks')
            else:
                host_dir = os.path.join('data', 'gencode_data', 'human', 'transcript_chunks')

            if os.path.exists(host_dir):
                chrom_files = [os.path.join(host_dir, f) for f in os.listdir(host_dir) if f.endswith('.fa')]
    else:
        chrom_dir = os.path.join('data', 'gencode_data', species, 'transcript_chunks')
        if os.path.exists(chrom_dir):
            chrom_files = [os.path.join(chrom_dir, f) for f in os.listdir(chrom_dir) if f.endswith('.fa')]

    if not chrom_files:
        print(f"Warning: No genome files found for {k}-mer safety check")
        return probes_file

    threads = args.parallelism or cpu_count()

    # Get or build Jellyfish database (cached per species)
    print(f"Loading Jellyfish {k}-mer database for {species} ({len(chrom_files)} genome files)...")
    db_start = time.time()
    db_path = get_jellyfish_db(chrom_files, k, species, threads)
    if not db_path:
        print("Error: Failed to get Jellyfish database, falling back to no k-mer filter")
        return probes_file
    db_elapsed = time.time() - db_start
    print(f"  Jellyfish database ready in {db_elapsed:.2f}s")

    # Write unique k-mers to FASTA for querying
    temp_kmer_file = os.path.join(output_base, 'temp_kmers.fa')
    kmer_list = sorted(unique_kmers)
    with open(temp_kmer_file, 'w', encoding='utf-8') as f:
        for i, kmer in enumerate(kmer_list):
            f.write(f">kmer_{i}\n{kmer}\n")

    # Query all k-mers against genome database
    print("Querying k-mer counts from genome database...")
    query_start = time.time()
    genome_counts = query_jellyfish_db(db_path, temp_kmer_file)
    query_elapsed = time.time() - query_start
    print(f"  Queried {len(genome_counts)} k-mers in {query_elapsed:.2f}s")

    # Count k-mer occurrences in source genome for self-match exclusion
    # Use full source genome from database (all contigs), not just the input fragment
    source_counts = {}
    source_seqs = []
    if source_transcripts and chrom_files:
        source_seqs = find_source_genome_sequences(chrom_files, source_transcripts)

    if not source_seqs and input_fasta and os.path.exists(input_fasta):
        source_seqs = [str(rec.seq) for rec in SeqIO.parse(input_fasta, "fasta")]

    if source_seqs:
        source_counts = count_source_kmers(source_seqs, k)
        print(f"  Counted {len(source_counts)} k-mers in source genome for self-match exclusion")

    # Filter probes: a probe is unsafe if any of its k-mers has off-target matches
    unsafe_probes = set()
    probe_match_details = {}
    total_self_matches = 0
    total_other_matches = 0

    for probe_id, kmers in probe_to_kmers.items():
        for kmer in kmers:
            genome_count = genome_counts.get(kmer, 0)
            source_count = source_counts.get(kmer, 0)
            off_target = genome_count - source_count
            if off_target > 0:
                total_other_matches += 1
                unsafe_probes.add(probe_id)
                if probe_id not in probe_match_details:
                    probe_match_details[probe_id] = {'total_matches': 0}
                probe_match_details[probe_id]['total_matches'] += 1
            elif genome_count > 0:
                total_self_matches += 1

    print("K-mer filtering complete:")
    print(f"  Self-gene matches (excluded): {total_self_matches:,}")
    print(f"  Other-gene matches: {total_other_matches:,}")
    print(f"  Unsafe probes (off-target k-mer matches): {len(unsafe_probes)}")

    # Write report
    report_start = time.time()
    match_report_file = f'{output_base}/kmer_matches_report.txt'

    all_probe_ids = {p.id for p in probes}
    safe_probes_ids = all_probe_ids - unsafe_probes

    with open(match_report_file, 'w', encoding='utf-8') as f:
        f.write(f"{k}-mer Match Report\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"Total probes checked: {len(probes)}\n")
        f.write(f"Unsafe probes (k-mer matches): {len(unsafe_probes)}\n")
        f.write(f"Safe probes (no matches): {len(safe_probes_ids)}\n\n")

        f.write("=" * 80 + "\n")
        f.write("UNSAFE PROBES\n")
        f.write("=" * 80 + "\n\n")
        for probe_id in sorted(unsafe_probes):
            if probe_id in probe_match_details:
                d = probe_match_details[probe_id]
                f.write(f"{probe_id}:\n")
                f.write(f"  Total off-target k-mer matches: {d['total_matches']}\n\n")

        f.write("\n" + "=" * 80 + "\n")
        f.write("SAFE PROBES\n")
        f.write("=" * 80 + "\n\n")

        for probe_id in sorted(safe_probes_ids):
            f.write(f"{probe_id}:\n")
            f.write(f"  Off-target k-mer matches: 0\n\n")

    report_elapsed = time.time() - report_start
    print(f"Report generation took {report_elapsed:.2f}s")

    safe_probes_file = f'{output_base}/safe_probes.fa'
    safe_probes = [p for p in probes if p.id not in unsafe_probes]

    # Secondary host check for microbiome species
    host_species = None
    if species in microbiome_species_list and args.align_microbiome and args.align_host:
        if species == "mouse-gut-microbiome":
            host_species = "mouse"
        else:
            host_species = "human"

        print(f"\nSecondary safety check: Checking safe probe k-mers against {host_species} transcripts...")

    if host_species:
        host_chrom_dir = os.path.join('data', 'gencode_data', host_species, 'transcript_chunks')
        if os.path.exists(host_chrom_dir):
            host_chrom_files = [os.path.join(host_chrom_dir, f) for f in os.listdir(host_chrom_dir) if f.endswith('.fa')]

            # Get or build host Jellyfish database (cached)
            print(f"  Loading {host_species} Jellyfish database...")
            host_db_path = get_jellyfish_db(host_chrom_files, k, host_species, threads)
            if not host_db_path:
                print(f"  Warning: Failed to build {host_species} Jellyfish database, skipping host check")
                host_species = None

            # Collect k-mers from safe probes only
            safe_probe_ids = {p.id for p in safe_probes}
            safe_kmers = set()
            safe_kmer_to_probes = {}
            for probe_id in safe_probe_ids:
                for kmer in probe_to_kmers[probe_id]:
                    safe_kmers.add(kmer)
                    if kmer not in safe_kmer_to_probes:
                        safe_kmer_to_probes[kmer] = []
                    safe_kmer_to_probes[kmer].append(probe_id)

            temp_safe_kmer_file = os.path.join(output_base, 'temp_safe_kmers.fa')
            with open(temp_safe_kmer_file, 'w', encoding='utf-8') as f:
                for i, kmer in enumerate(sorted(safe_kmers)):
                    f.write(f">kmer_{i}\n{kmer}\n")

            host_counts = query_jellyfish_db(host_db_path, temp_safe_kmer_file)

            # Any k-mer with count > 0 in host means probe matches host
            host_unsafe = set()
            for kmer, count in host_counts.items():
                if count > 0 and kmer in safe_kmer_to_probes:
                    for probe_id in safe_kmer_to_probes[kmer]:
                        host_unsafe.add(probe_id)

            safe_probes = [p for p in safe_probes if p.id not in host_unsafe]

            print(f"  Probes matching {host_species} transcripts: {len(host_unsafe)}")
            print(f"  Final safe probes after {host_species} check: {len(safe_probes)}")

            os.remove(temp_safe_kmer_file)
        else:
            print(f"  Warning: {host_species.capitalize()} transcript directory not found, skipping host check")

    SeqIO.write(safe_probes, safe_probes_file, "fasta")

    print("Safety check results:")
    print(f"  Total non-aligned probes: {len(probes)}")
    print(f"  Unsafe probes (have {k}-mer matches): {len(unsafe_probes)}")
    print(f"  Safe probes: {len(safe_probes)}")
    print(f"  Safe probes written to: {safe_probes_file}")
    print(f"  Match details report: {match_report_file}")

    # Cleanup temp file (DB is cached for reuse)
    os.remove(temp_kmer_file)

    return safe_probes_file
    

def main():
    """Main pipeline logic."""
    parser = argparse.ArgumentParser(description="Probe design")

    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument("--gene-sequence", "-gs", help="Gene sequence in FASTA format (pasted text)")
    input_group.add_argument("--probe-sequence", "-ps", help="Probe sequence in FASTA format (pasted text)")
    input_group.add_argument("--gene-sequence-stdin", action="store_true", help="Read gene sequence from stdin")
    input_group.add_argument("--probe-sequence-stdin", action="store_true", help="Read probe sequence from stdin")

    parser.add_argument("--species", "-s", default="human", help="Species (human or mouse)")
    parser.add_argument("--probe-length", type=int, default=36, help="Probe length (auto-detected if probes provided)")
    parser.add_argument("--max-mismatches", type=int, default=2, help="Max mismatches allowed in alignments")
    parser.add_argument("--parallelism", "-p", type=int, default=None, help="Number of parallel worker threads")
    parser.add_argument("--task-id", help="Task ID for organizing outputs")
    parser.add_argument("--kmer-length", type=int, default=18,
                       help="K-mer length for safety check (minimum: 14, default: 18)")
    parser.add_argument("--skip-annotation", action="store_true",
                       help="Skip gene name annotation of SAM files")
    parser.add_argument("--align-microbiome", action="store_true", default=False,
                   help="Align against microbiome genomes")
    parser.add_argument("--align-host", action="store_true", default=False,
                    help="Align against host transcriptome")
    parser.add_argument("--tm-range", type=str, default="42-47",
                    help="Tm filter range in °C, e.g. 42-47 (default: 42-47)")

    args = parser.parse_args()
    if args.kmer_length < 14:
        print(f"Error: K-mer length {args.kmer_length} is too small. Minimum k-mer length is 14.")
        sys.exit(1)
    if args.kmer_length > args.probe_length:
        print(f"Error: K-mer length ({args.kmer_length}) exceeds probe length ({args.probe_length}). Please adjust k-mer length and try again.")
        sys.exit(1)
    if args.gene_sequence_stdin:
        args.gene_sequence = sys.stdin.read()
    elif args.probe_sequence_stdin:
        args.probe_sequence = sys.stdin.read()
    print("PROBE DESIGN")
    print("=" * 80)
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    annotator = samannotator() if not args.skip_annotation else None
    if annotator:
        print("Gene annotation: Enabled")
    else:
        print("Gene annotation: Skipped")

    input_fasta = None
    skip_probe_generation = False
    scenario = None
    
    if args.gene_sequence:
        scenario = "gene_sequence"
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as tmp:
            tmp.write(args.gene_sequence)
            input_fasta = tmp.name
        print(f"Scenario 1: Gene sequence (pasted) -> {input_fasta}")
        
    elif args.probe_sequence:
        scenario = "probe_sequence"
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as tmp:
            tmp.write(args.probe_sequence)
            input_fasta = tmp.name
        skip_probe_generation = True
        print(f"Scenario 2: Probe sequence (pasted) -> {input_fasta}")

    print("\nValidating FASTA format...")
    valid, error_msg, sequences = validate_fasta_format(input_fasta)
    if not valid:
        print(f"FASTA validation failed: {error_msg}")
        return False
    print(f"FASTA valid: {len(sequences)} sequence(s) found")

    if not skip_probe_generation and len(sequences) > 1:
        print(f"Error: Only a single FASTA sequence is allowed per job. Found {len(sequences)} sequences.")
        return False

    if skip_probe_generation:
        detected_length = detect_probe_length(sequences)
        if detected_length:
            args.probe_length = detected_length
            print(f"Detected probe length: {args.probe_length} bp")
        else:
            print(f"Warning: Inconsistent probe lengths detected, using default: {args.probe_length} bp")
        
        for header, seq in sequences:
            if len(seq) != args.probe_length:
                print(f"Probe length mismatch: {header} has {len(seq)} bp, expected {args.probe_length} bp")
                return False
    if args.species in ["gut-microbe", "human-oral-microbiome", "human-skin-microbiome", "human-vaginal-microbiome", "mouse-gut-microbiome"]:
        genome = None
        print(f"Using pre-existing {args.species} data")
    else:
        print(f"\nValidating genome for {args.species}...")
        species_data_dir = os.path.join('data', 'gencode_data', args.species)
        try:
            success = validate_local_gencode(args.species, base_dir='.')

            if not success:
                print("Genome validation failed. Ensure Gencode files are in: " + species_data_dir)
                return False

            local_genome_release = get_local_release_number(
                r'GRCh38\.p14\.genome\.release_(\d+)\.fa' if args.species == "human" else r'GRCm39\.genome\.release_(\d+)\.fa', species_data_dir
            )

            if args.species == "human":
                genome_file = os.path.join(species_data_dir, f'GRCh38.p14.genome.release_{local_genome_release}.fa')
            else:  # mouse
                genome_file = os.path.join(species_data_dir, f'GRCm39.genome.release_{local_genome_release}.fa')

            if not os.path.exists(genome_file):
                print(f"Genome file not found: {genome_file}")
                return False

            genome = genome_file
            print(f"Genome ready: {genome}")
        except (FileNotFoundError, PermissionError, OSError) as e:
            print(f"Error processing genome: {e}")
            traceback.print_exc()
            return False

    if annotator:
        print("\n📋 PREPARING GENE ANNOTATION MAPPINGS")
        print("-" * 60)
        annotator._ensure_mappings_loaded(args.species)

    if args.task_id:
        output_base = f'outputs/alignments/{args.task_id}'
    else:
        output_base = 'outputs'
    os.makedirs(output_base, exist_ok=True)

    probes_out = f'{output_base}/candidate_probes.fa'
    sam_out = f'{output_base}/probe_alignments.sam'
    filtered_sam = f'{output_base}/filtered_probe_alignments.sam'
    annotated_sam = f'{output_base}/filtered_probe_alignments_annotated.sam'
    align_dir = f'{output_base}/chroms'

    if args.species in ["gut-microbe", "human-oral-microbiome", "human-skin-microbiome", "human-vaginal-microbiome", "mouse-gut-microbiome"]:
        microbiome_dir = os.path.join('data', args.species)
        if args.species == "mouse-gut-microbiome":
            host_dir = os.path.join('data', 'gencode_data', 'mouse', 'transcript_chunks')
        else:
            host_dir = os.path.join('data', 'gencode_data', 'human', 'transcript_chunks')
        dirs = []
        if args.align_microbiome:
            dirs.append(microbiome_dir)
        if args.align_host and os.path.exists(host_dir):
            dirs.append(host_dir)
        
        if not dirs:
            print("Error: At least one alignment target must be selected for microbiome species")
            return False
        
        chrom_dir = ':'.join(dirs)
    else:
        chrom_dir = os.path.join('data', 'gencode_data', args.species, 'transcript_chunks')

    os.makedirs(align_dir, exist_ok=True)

    if skip_probe_generation:
        print("\nSkipping probe generation - using provided probes")

        total, passing = count_and_filter_gc(input_fasta, min_gc=40, max_gc=80)
        
        SeqIO.write(passing, probes_out, "fasta")
        
        print(f"Found {len(passing)} candidate probes passing GC filter (40-80%, length={args.probe_length}).")
        
        if len(passing) == 0:
            print("No probes passed GC filter")
            return False
    else:
        print("\nGenerating probes from input sequence...")
        gen_cmd = [
            sys.executable, os.path.join('core', 'generate_probes.py'),
            '--probe-length', str(args.probe_length),
            '--max-mismatches', str(args.max_mismatches),
            '--input', input_fasta,
            '--output', probes_out
        ]
        if not run(gen_cmd, 'Generate candidate probes', check_output_file=probes_out):
            return False

    # Tm and homopolymer filter: remove probes outside Tm range or with homopolymer runs
    tm_min, tm_max = [float(x) for x in args.tm_range.split('-')]
    scorer = ThermodynamicProbeScorer()
    candidates = list(SeqIO.parse(probes_out, 'fasta'))
    passed = []
    tm_rejected = 0
    homopolymer_rejected = 0
    for record in candidates:
        seq = str(record.seq)
        tm = scorer.calculate_tm(seq)
        if not (tm_min <= tm <= tm_max):
            tm_rejected += 1
            continue
        if scorer.check_homopolymer_runs(seq):
            homopolymer_rejected += 1
            continue
        passed.append(record)
    print(f"\nTm filter ({tm_min}-{tm_max}°C): rejected {tm_rejected}/{len(candidates)} probes")
    print(f"Homopolymer filter (≥{scorer.max_homopolymer}bp runs): rejected {homopolymer_rejected}/{len(candidates)} probes")
    print(f"Passed both filters: {len(passed)}/{len(candidates)} probes")
    if len(passed) == 0:
        print("No probes passed Tm and homopolymer filters")
        return False
    with open(probes_out, 'w') as f:
        SeqIO.write(passed, f, 'fasta')

    razers = os.path.join('bin', 'razers3')
    if not os.path.exists(razers):
        print(f"razers3 binary not found at {razers}")
        return False
    identity_percent = int(((args.probe_length - args.max_mismatches) / args.probe_length) * 100)
    print("\nAlignment parameters:")
    print(f"Probe length: {args.probe_length} bp")
    print(f"Max mismatches: {args.max_mismatches}")
    print(f"Calculated identity: {identity_percent}%")

    parallelism = args.parallelism
    if parallelism is None or parallelism <= 0:
        parallelism = cpu_count()

    chrom_dirs = chrom_dir.split(':')
    existing_chroms = []
    for dir_path in chrom_dirs:
        if not os.path.exists(dir_path):
            print(f"Warning: Directory not found: {dir_path}")
            continue
        existing_chroms.extend([os.path.join(dir_path, f) for f in os.listdir(dir_path) if f.endswith('.fa') or f.endswith('.fna')])

    if not existing_chroms:
        print(f"No chromosome files found in {chrom_dir}")
        return False

    def align_chrom(chrom_fasta):
        base = os.path.splitext(os.path.basename(chrom_fasta))[0]
        out_sam = os.path.join(align_dir, f'{base}.sam')
        cmd = [razers, '-ng', '-i', str(identity_percent), '-rr', '100', '-m', '100','-tc', '1', '-o', out_sam, chrom_fasta, probes_out]
        start = time.time()
        proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
        elapsed = time.time() - start
        success = proc.returncode == 0 and os.path.exists(out_sam) and os.path.getsize(out_sam) > 0
        return (out_sam, success, elapsed) if success else (None, False, elapsed)

  
    print(f"\nUsing {len(existing_chroms)} pre-chunked transcript files from {chrom_dir}/")
    chrom_files = sorted(existing_chroms)

    print(f"\nRunning razers3 on {len(chrom_files)} files with {parallelism} workers")
    alignment_start = time.time()
    sam_files = []
    failed_count = 0
    
    with concurrent.futures.ThreadPoolExecutor(max_workers=parallelism) as ex:
        futures = {ex.submit(align_chrom, cf): cf for cf in chrom_files}
        for fut in concurrent.futures.as_completed(futures):
            chrom = futures[fut]
            try:
                out_sam, success, _ = fut.result()
                if success:
                    sam_files.append(out_sam)
                else:
                    failed_count += 1
            except (subprocess.SubprocessError, FileNotFoundError, PermissionError, OSError, TimeoutError) as e:
                print(f"Exception aligning {chrom}: {e}")
                failed_count += 1
                
    alignment_elapsed = time.time() - alignment_start

    successful = len(sam_files)
    total_chroms = len(chrom_files)
    
    print(f"\nAlignment summary: {successful} successful, {failed_count} failed (Total: {total_chroms}) in {alignment_elapsed:.2f}s")
    if not sam_files:
        print("No SAM files were produced by the aligner. Exiting.")
        return False

    print(f"\nMerging {len(sam_files)} SAM files into {sam_out}")
    seen = set()
    header_lines = []
    with open(sam_out, 'w', encoding='utf-8') as outf:
        for idx, sf in enumerate(sam_files):
            with open(sf, 'r', encoding='utf-8') as inf:
                for line in inf:
                    if line.startswith('@'):
                        if idx == 0 and line not in header_lines:
                            header_lines.append(line)
                    else:
                        break

        for sf in sam_files:
            with open(sf, 'r', encoding='utf-8') as inf:
                for line in inf:
                    if line.startswith('@'):
                        continue
                    key = line.strip()
                    if key not in seen:
                        seen.add(key)
                        outf.write(line)

    print(f"Merged SAM written to {sam_out} ({len(seen):,} alignments)")

    print(f"\nFiltering merged SAM -> keeping alignments with NM <= {args.max_mismatches}")
    nm_pattern = re.compile(r'NM:i:(\d+)')
    total = 0
    kept = 0
    skipped_no_nm = 0
    aligned_probes = set()
    non_aligned_probes_fa = f'{output_base}/non_aligned_probes.fa'

    with open(sam_out, 'r', encoding='utf-8') as inf, \
        open(filtered_sam, 'w', encoding='utf-8') as outf:

        for line in inf:
            if line.startswith('@'):
                continue
            total += 1
            parts = line.split('\t')
            if len(parts) < 11:
                continue

            qname = parts[0]
            aligned_probes.add(qname)

            m = nm_pattern.search(line)
            if not m:
                skipped_no_nm += 1
                continue
            nm = int(m.group(1))

            if nm <= args.max_mismatches:
                outf.write(line)
                kept += 1
    
    microbiome_species = ["gut-microbe", "human-oral-microbiome", "human-skin-microbiome",
                      "human-vaginal-microbiome", "mouse-gut-microbiome"]
    
    if args.species in microbiome_species:
        annotate_microbiome_sam(filtered_sam, args.species, output_base, args.align_host)

    print("\nExtracting non-aligned probes...")
    non_aligned_count = 0
    candidate_probes = list(SeqIO.parse(probes_out, "fasta"))

    is_microbiome = args.species in microbiome_species
    if is_microbiome and args.align_host:
        host_species = 'human' if args.species in ['gut-microbe', 'human-oral-microbiome', 'human-skin-microbiome', 'human-vaginal-microbiome'] else 'mouse'
        
        temp_sam = filtered_sam.replace('filtered_probe_alignments.sam', 'filtered_probe_alignments_temp.sam')
        annotator.annotate_sam(annotated_sam, temp_sam, host_species)
        shutil.move(temp_sam, annotated_sam) 

    elif annotator and not is_microbiome:
        annotator.annotate_sam(filtered_sam, annotated_sam, args.species)

    sam_to_analyze = annotated_sam if os.path.exists(annotated_sam) else filtered_sam

    # Infer source transcripts from alignment patterns
    gene_mappings = annotator.mappings if annotator and annotator.mappings else None
    source_info = infer_source_transcripts(filtered_sam, gene_mappings=gene_mappings, is_microbiome=is_microbiome)
    source_transcripts = source_info['source_transcripts']
    source_gene = source_info['source_gene']

    probe_classifications = classify_probes_from_sam(
        sam_to_analyze, is_microbiome=is_microbiome,
        source_transcripts=source_transcripts, source_gene=source_gene
    )
    # Count self-aligning vs off-target probes
    self_aligning_probes = sum(1 for p_id, data in probe_classifications.items() if data['status'] == 'safe')
    off_target_probes = len(aligned_probes) - self_aligning_probes

    # Extract all on-target probes (self-aligned safe + non-aligned) for k-mer check.
    # K-mer check is mandatory for all on-target probes: if any k-mer matches
    # a gene other than the source, the probe is not safe.
    print("\nExtracting probes...")
    self_aligned_safe_count = 0

    with open(non_aligned_probes_fa, 'w', encoding='utf-8') as out:
        for record in candidate_probes:
            probe_id = record.id
            if probe_id not in probe_classifications:
                SeqIO.write(record, out, "fasta")
                non_aligned_count += 1
            elif probe_classifications[probe_id]['status'] == 'safe':
                SeqIO.write(record, out, "fasta")
                self_aligned_safe_count += 1

    print("\nFiltering summary:")
    print(f"Total alignments scanned: {total:,}")
    print(f"Kept in filtered (NM<={args.max_mismatches}): {kept:,}")
    print(f"Total candidate probes: {len(candidate_probes):,}")
    print(f"Aligned probes found: {len(aligned_probes):,}")
    print(f"  - Self-alignments only (safe): {self_aligning_probes:,}")
    print(f"  - Off-target alignments: {off_target_probes:,}")
    if source_gene:
        print(f"  - Inferred source gene: {source_gene}")
    print(f"Self-aligned safe probes: {self_aligned_safe_count:,}")
    print(f"Non-aligned probes: {non_aligned_count:,}")
    print(f"Skipped (no NM tag): {skipped_no_nm:,}")

    total_on_target = self_aligned_safe_count + non_aligned_count
    safe_probes_file = f'{output_base}/safe_probes.fa'
    safe_scores_csv = f'{output_base}/safe_probes_scores.csv'

    if total_on_target > 0:
        print(f"\nPerforming {args.kmer_length}-mer safety check on {total_on_target:,} on-target probes...")
        print(f"  ({self_aligned_safe_count:,} self-aligned + {non_aligned_count:,} non-aligned)")
        kmer_result_file = check_14mer_safety(
            non_aligned_probes_fa, args.species, output_base, args,
            source_transcripts=source_transcripts,
            input_fasta=input_fasta
        )
        total_safe_count = len(list(SeqIO.parse(kmer_result_file, 'fasta')))
    else:
        print(f"\nNo on-target probes to check for {args.kmer_length}-mer safety")
        total_safe_count = 0
        kmer_result_file = None

    # Copy k-mer safe probes to final safe_probes.fa
    safe_records = list(SeqIO.parse(kmer_result_file, 'fasta')) if total_safe_count > 0 and kmer_result_file else []
    with open(safe_probes_file, 'w', encoding='utf-8') as out:
        for record in safe_records:
            SeqIO.write(record, out, "fasta")

    print(f"\nTotal safe probes: {total_safe_count:,}")

    if total_safe_count > 0:
        print(f"\nScoring {total_safe_count} safe probes...")
        score_and_save_probes(safe_probes_file, safe_scores_csv)
    else:
        print("\nNo safe probes to score - creating empty score file...")
        with open(safe_scores_csv, 'w', encoding='utf-8') as f:
            f.write("No safe probes found\n")
        safe_scores_txt = safe_scores_csv.replace('.csv', '.txt')
        with open(safe_scores_txt, 'w', encoding='utf-8') as f:
            f.write("No safe probes found\n")

    print('\n' + '='*80)
    print('PIPELINE COMPLETED SUCCESSFULLY')
    print('='*80)
    print(f'Filtered SAM: {filtered_sam}')
    if annotator and os.path.exists(annotated_sam):
        size = os.path.getsize(annotated_sam)
        if size > 1024*1024:
            size_str = f"({size/(1024*1024):.1f} MB)"
        else:
            size_str = f"({size/1024:.1f} KB)"
        print(f'Annotated SAM: {annotated_sam} {size_str}')
    print(f'Non-aligned probes FASTA: {non_aligned_probes_fa}')
    if total_safe_count > 0:
        print(f'Safe probes FASTA: {safe_probes_file}')
        print(f'Safe probes scores: {safe_scores_csv}')
    if total_on_target > 0:
        print(f'k-mer match report: {output_base}/kmer_matches_report.txt')
    
    cleanup_dirs = [
        os.path.join(output_base, 'chroms'),
    ]
    
    for cleanup_dir in cleanup_dirs:
        if os.path.exists(cleanup_dir):
            shutil.rmtree(cleanup_dir)

    print('\n' + '='*80)
    print('PIPELINE COMPLETED SUCCESSFULLY')
    print('='*80)
    return True


if __name__ == '__main__':
    os.makedirs('core', exist_ok=True)
    os.makedirs('bin', exist_ok=True)
    os.makedirs(DEFAULT_DIRS['gene'], exist_ok=True)
    
    temp_file = None
    if len(sys.argv) > 1:
        parser = argparse.ArgumentParser()
        parser.add_argument("--gene-sequence", "-gs")
        parser.add_argument("--probe-sequence", "-ps")
        temp_args, unknown = parser.parse_known_args()
        if temp_args.gene_sequence or temp_args.probe_sequence:
            pass

    main()