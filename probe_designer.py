"""
This script handles two input scenarios:
1. Gene sequence (pasted) -> generate probes from gene -> align to species genome
2. Probe sequence (pasted) -> skip generation -> align to species genome

Usage:
  # Scenario 1: Gene sequence (pasted)
  python probe_designer.py --gene-sequence ">gene1\nATCG..." --species human --probe-length 36 --max-mismatches 2
  
  # Scenario 2: Probe sequence (pasted)
  python probe_designer.py --probe-sequence ">probe1\nATCG..." --species human --max-mismatches 2

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
import shlex
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
from probe_metrics_report import write_probe_metrics_report
from probe_metrics import ProbeMetricsCalculator
from probe_classifier import classify_probes_from_sam, infer_source_transcripts, resolve_source_from_header
from host_internal_filter import filter_probes_by_self_alignment
from backend.sam_parser import compute_alignment_diffs
import genome_stage


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
    """Locate and validate pre-existing local host data (genome + chunks).

    All host-specific paths and genome file patterns are read from
    config/organisms.yml via the organism registry — no per-host
    hard-coded paths or regex live here.
    """
    from organism_registry import get_host
    host = get_host(species)
    if host is None:
        print(f"Error: Unknown host '{species}'. Add an entry to config/organisms.yml.")
        return False
    if base_dir is None:
        base_dir = os.path.dirname(os.path.abspath(__file__))
    species_data_dir = (host.genome_dir if os.path.isabs(host.genome_dir)
                        else os.path.join(base_dir, host.genome_dir))
    chunk_dir = (host.transcript_chunks_dir if os.path.isabs(host.transcript_chunks_dir)
                 else os.path.join(base_dir, host.transcript_chunks_dir))
    local_release = get_local_release_number(host.genome_file_regex, species_data_dir)
    if not local_release:
        print(f"Error: No genome file matching '{host.genome_file_regex}' in {species_data_dir}")
        return False
    pat = re.compile(host.genome_file_regex)
    genome_file = None
    if os.path.isdir(species_data_dir):
        for fn in os.listdir(species_data_dir):
            m = pat.match(fn)
            if m and int(m.group(1)) == local_release:
                genome_file = os.path.join(species_data_dir, fn)
                break
    if not genome_file or not os.path.exists(genome_file):
        print(f"Error: Genome file for release {local_release} not found in {species_data_dir}")
        return False
    if not os.path.exists(chunk_dir):
        print(f"Error: Transcript chunks directory not found: {chunk_dir}")
        return False
    existing_chunks = [f for f in os.listdir(chunk_dir) if f.endswith('.fa')]
    if not existing_chunks:
        print(f"Error: No chunk .fa files in {chunk_dir}")
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
        """Find GTF file for species and extract release number.

        The GTF directory is taken from the organism registry
        (config/organisms.yml → host.genome_dir). Falls back to the
        legacy `data/gencode_data/<species>/` layout if the host
        isn't in the registry (keeps microbiome catalogs working).
        """
        try:
            from organism_registry import get_host
            host = get_host(species)
            species_dir = Path(host.genome_dir) if host else (self.gtf_dir / species)
        except Exception:
            species_dir = self.gtf_dir / species
        if not species_dir.exists():
            return None, None

        
        for pattern in ('*.gtf', '*.gff3', '*.gff'):
            for file in species_dir.glob(pattern):
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
        """Extract transcript→gene mappings from GTF or GFF3.

        GTF format (Ensembl vertebrates):
            ... transcript_id "X"; gene_name "Y"; ...
        GFF3 format (Ensembl Plants, some non-vertebrate genomes):
            ... ID=transcript:X;Parent=gene:Y;Name=Z;...
        Both are handled here so a single attribute file (whatever Ensembl
        ships for the organism) is enough.
        """
        is_gff3 = str(gtf_file).lower().endswith(('.gff3', '.gff'))
        fmt = 'GFF3' if is_gff3 else 'GTF'
        print(f"Extracting transcript→gene mappings from {fmt} (release_{release})...")
        start = time.time()

        if not is_gff3:
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
        else:
    
            t_id_re = re.compile(r'(?:^|;)ID=(?:transcript:)?([^;]+)')
            t_parent_re = re.compile(r'(?:^|;)Parent=(?:gene:)?([^;]+)')
            t_name_re = re.compile(r'(?:^|;)Name=([^;]+)')
            transcript_types = {
                'mRNA', 'transcript', 'ncRNA', 'rRNA', 'tRNA', 'snRNA',
                'snoRNA', 'miRNA', 'lnc_RNA', 'pseudogenic_transcript',
            }
            gene_names = {}
            with open(gtf_file, 'r', encoding='utf-8') as f:
                for line in f:
                    if line.startswith('#') or not line.strip():
                        continue
                    parts = line.split('\t')
                    if len(parts) < 9 or parts[2] != 'gene':
                        continue
                    attrs = parts[8]
                    g_match = re.search(r'(?:^|;)ID=(?:gene:)?([^;]+)', attrs)
                    n_match = t_name_re.search(attrs)
                    if g_match:
                        gid = g_match.group(1)
                        gene_names[gid] = n_match.group(1) if n_match else gid
            with open(gtf_file, 'r', encoding='utf-8') as f:
                for line in f:
                    if line.startswith('#') or not line.strip():
                        continue
                    parts = line.split('\t')
                    if len(parts) < 9 or parts[2] not in transcript_types:
                        continue
                    attrs = parts[8]
                    t_match = t_id_re.search(attrs)
                    p_match = t_parent_re.search(attrs)
                    if t_match and p_match:
                        tid = t_match.group(1)
                        gid = p_match.group(1)
                        self.mappings[tid] = gene_names.get(gid, gid)
        
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

def _selected_microbiome_ids(args):
    """Microbiome IDs the run targets: the --microbiomes list, or [--species] as fallback."""
    if getattr(args, 'microbiomes', ''):
        return [m for m in (m.strip() for m in args.microbiomes.split(',')) if m]
    return [args.species]


def annotate_microbiome_sam(sam_file, microbiome_ids, _, align_host=False):
    """Add species annotations to SAM file for microbiome data."""
    metadata_files = [os.path.join('data', mb_id, 'genomes-all_metadata.tsv') for mb_id in microbiome_ids]
    metadata_files = [f for f in metadata_files if os.path.exists(f)]

    if not metadata_files:
        print(f"Warning: No metadata files found for: {', '.join(microbiome_ids)}")
        return False

    print(f"\nAnnotating SAM file with species information from {', '.join(metadata_files)}...")

    genome_to_species = {}
    try:
        for metadata_file in metadata_files:
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
            from organism_registry import host_for_microbiome
            host_obj = host_for_microbiome(microbiome_ids[0])
            host_species = host_obj.id if host_obj else None
            if host_species:
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
    """
    jellyfish = 'jellyfish'
    species_dir = os.path.dirname(genome_files[0])
    db_path = os.path.join(species_dir, f'kmer_{k}mer_counts.jf')

    if os.path.exists(db_path):
        print(f"  Using pre-built Jellyfish database: {db_path}")
        return db_path

    print(f"  Building Jellyfish database (one-time, will be cached for future runs)...")
    # shell=True is required here: the command is a genuine shell pipeline that
    # streams the concatenated genomes into Jellyfish via /dev/fd/0. Every
    # interpolated path therefore goes through shlex.quote() rather than the
    # hand-written single quotes used previously, which a path containing a
    # quote character would have escaped straight into arbitrary command
    # execution.
    file_list = ' '.join(shlex.quote(f) for f in genome_files)
    cmd = (
        f"cat {file_list} | {shlex.quote(jellyfish)} count "
        f"-m {int(k)} -s 5G -t {int(threads)} -C -o {shlex.quote(db_path)} /dev/fd/0"
    )
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


def _gpu_kmer_enabled(args, k):
    """Whether GPU k-mer counting should be used for this step."""
    mode = getattr(args, 'gpu', 'auto')
    if mode == 'off':
        return False
    try:
        from gpu_kmer import gpu_available, MAX_GPU_K
    except Exception as e:
        if mode == 'on':
            print(f"  Warning: --gpu on but GPU engine unavailable ({e}); using Jellyfish")
        return False
    if k > MAX_GPU_K:
        if mode == 'on':
            print(f"  Warning: k={k} exceeds GPU limit ({MAX_GPU_K}); using Jellyfish for this step")
        return False
    if not gpu_available():
        if mode == 'on':
            print("  Warning: --gpu on but no usable GPU/CuPy detected; using Jellyfish")
        return False
    return True


def genome_kmer_counts(files, k, wanted_kmers, args, kmer_fasta, device=0,
                       threads=1, label='genome'):
    """Canonical genome k-mer counts for `wanted_kmers`.

    Uses the GPU engine when enabled/available (equivalent to `jellyfish -C`
    query, validated for parity); otherwise builds/queries a Jellyfish DB.
    Returns {wanted_kmer: count} or None on failure.
    """
    if _gpu_kmer_enabled(args, k):
        from gpu_kmer import canonical_counts
        workers = threads if threads and threads > 0 else None
        return canonical_counts(files, k, wanted_kmers,
                                device=device, workers=workers)
    db_path = get_jellyfish_db(files, k, label, threads)
    if not db_path:
        return None
    return query_jellyfish_db(db_path, kmer_fasta)


def count_self_genome_kmer_hits(sequences, k, probe_kmers):
    """Canonical occurrences of each probe k-mer within its own-genome sequences.

    Matches Jellyfish `-C` canonical counting (each genome position contributes
    to the canonical form of its k-mer) but only tracks the given probe k-mers,
    so there is no per-position reverse-complement over whole genomes. Returns
    {probe_kmer: count}. Used to subtract a probe's own-genome hits from the
    combined-genome count so only OTHER genomes count as off-target.
    """
    rc_of = {}
    wanted = set()
    for pk in probe_kmers:
        rc = str(Seq(pk).reverse_complement())
        rc_of[pk] = rc
        wanted.add(pk)
        wanted.add(rc)

    literal_counts = {}
    for seq in sequences:
        s = seq.upper()
        for i in range(len(s) - k + 1):
            kmer = s[i:i+k]
            if kmer in wanted:
                literal_counts[kmer] = literal_counts.get(kmer, 0) + 1

    counts = {}
    for pk in probe_kmers:
        rc = rc_of[pk]
        if rc == pk:
            counts[pk] = literal_counts.get(pk, 0)
        else:
            counts[pk] = literal_counts.get(pk, 0) + literal_counts.get(rc, 0)
    return counts


def find_source_genome_sequences(chrom_files, source_transcripts):
    """Find all sequences from the source in the database chunk files.

    Microbiome (MGYG-style, e.g. `MGYG000000001_1`): match every contig whose
    ID starts with the genome prefix (`MGYG000000001_*`).
    Host (ENST/ENSMUST, e.g. `ENST00000275493.7`): exact-match the transcript
    ID, version-tolerant — so all isoforms of the source gene contribute their
    k-mers to the self-exclusion set.
    """
    source_prefixes = set()
    source_ids = set()
    source_ids_versionless = set()
    for tid in source_transcripts:
        parts = tid.rsplit('_', 1)
        if len(parts) == 2 and parts[1].isdigit():
            source_prefixes.add(parts[0] + '_')
        else:
            source_ids.add(tid)
            source_ids_versionless.add(tid.split('.', 1)[0])

    # Fast path: a microbiome genome's contigs (MGYG000000001_* etc.) live in a
    # single file named by the genome id (MGYG000000001.fa), so read only that
    # file instead of scanning the whole database. Host transcript ids (ENST…)
    # are not filename-addressable, so if any are present fall back to a full
    # scan. The prefix/id matching below is unchanged, so the result is identical.
    candidates = None
    if source_prefixes and not source_ids:
        genome_ids = {p[:-1] for p in source_prefixes}
        cand = [f for f in chrom_files
                if os.path.splitext(os.path.basename(f))[0] in genome_ids]
        if cand:
            candidates = cand
    if candidates is None:
        candidates = chrom_files

    source_seqs = []
    for chrom_file in candidates:
        for rec in SeqIO.parse(chrom_file, "fasta"):
            rid = rec.id
            if source_prefixes and any(rid.startswith(p) for p in source_prefixes):
                source_seqs.append(str(rec.seq))
                continue
            if rid in source_ids or rid.split('.', 1)[0] in source_ids_versionless:
                source_seqs.append(str(rec.seq))
    return source_seqs


def check_14mer_safety(probes_file, species, output_base, args, source_transcripts=None,
                       genome_self_match=False, precomputed_group_counts=None):
    """Check if k-mers from on-target probes have off-target genome matches.

    When `precomputed_group_counts` ({group_label: {kmer: count}}) is provided
    — computed during the fused GPU alignment scan — the genome database is not
    scanned again; per-group counts are looked up instead.
    """
    k = args.kmer_length
    print(f"\nChecking {k}-mer safety for on-target probes...")
    probes = list(SeqIO.parse(probes_file, "fasta"))
    if not probes:
        return probes_file

    probe_to_kmers = {}
    unique_kmers = set()
    for record in probes:
        kmers = generate_14mers(str(record.seq), k)
        probe_to_kmers[record.id] = kmers
        unique_kmers.update(kmers)

    total_mers = sum(len(v) for v in probe_to_kmers.values())
    print(f"Generated {len(unique_kmers)} unique k-mers from {len(probes)} probes")
    print(f"  Total k-mers: {total_mers}, Unique: {len(unique_kmers)}, Duplicates removed: {total_mers - len(unique_kmers)} ({(total_mers-len(unique_kmers))/total_mers*100:.1f}%)")

    from organism_registry import load_registry as _load_reg_kmer
    _reg = _load_reg_kmer()

    chrom_files = []
    db_groups = [] 

    if _reg.is_microbiome_species(species):
        host_obj, _ = _reg.get_microbiome(species)
        if args.align_microbiome:
            for mb_id in _selected_microbiome_ids(args):
                _, mb_obj = _reg.get_microbiome(mb_id)
                if mb_obj and os.path.exists(mb_obj.data_dir):
                    files = [os.path.join(mb_obj.data_dir, f) for f in os.listdir(mb_obj.data_dir)
                             if f.endswith('.fa') or f.endswith('.fna')]
                    if files:
                        db_groups.append((mb_id, files))
                        chrom_files.extend(files)
        elif args.align_host:
            host_dir = host_obj.transcript_chunks_dir
            if os.path.exists(host_dir):
                files = [os.path.join(host_dir, f) for f in os.listdir(host_dir) if f.endswith('.fa')]
                if files:
                    db_groups.append((host_obj.id, files))
                    chrom_files.extend(files)
    else:
        host_obj = _reg.get_host(species)
        chrom_dir = host_obj.transcript_chunks_dir if host_obj else None
        if chrom_dir and os.path.exists(chrom_dir):
            files = [os.path.join(chrom_dir, f) for f in os.listdir(chrom_dir) if f.endswith('.fa')]
            if files:
                db_groups.append((species, files))
                chrom_files.extend(files)

    if not db_groups:
        print(f"Warning: No genome files found for {k}-mer safety check")
        return probes_file

    threads = args.parallelism or cpu_count()

    temp_kmer_file = os.path.join(output_base, 'temp_kmers.fa')
    kmer_list = sorted(unique_kmers)
    with open(temp_kmer_file, 'w', encoding='utf-8') as f:
        for i, kmer in enumerate(kmer_list):
            f.write(f">kmer_{i}\n{kmer}\n")

    have_precomputed = (precomputed_group_counts is not None and
                        all(lbl in precomputed_group_counts for lbl, _ in db_groups))
    if have_precomputed:
        backend = "fused GPU scan (precomputed)"
    else:
        backend = f"GPU (device {args.gpu_device})" if _gpu_kmer_enabled(args, k) else "Jellyfish"
    print(f"Querying k-mer counts from genome database(s) [{backend}]...")
    query_start = time.time()
    genome_counts = {}
    for label, files in db_groups:
        if have_precomputed:
            counts = precomputed_group_counts[label]
        else:
            print(f"  Counting {k}-mers for {label} ({len(files)} genome files)...")
            counts = genome_kmer_counts(files, k, kmer_list, args, temp_kmer_file,
                                        device=args.gpu_device, threads=threads,
                                        label=label)
        if counts is None:
            print(f"  Error: Failed to count k-mers for {label}, skipping")
            continue
        for kmer, c in counts.items():
            genome_counts[kmer] = genome_counts.get(kmer, 0) + c
    query_elapsed = time.time() - query_start
    print(f"  Queried {len(genome_counts)} k-mers in {query_elapsed:.2f}s")

    if not genome_counts:
        print("Error: Failed to query any genome database, falling back to no k-mer filter")
        return probes_file

    
    source_counts = {}
    source_seqs = []
    if source_transcripts and chrom_files:
        source_seqs = find_source_genome_sequences(chrom_files, source_transcripts)

    if source_seqs:
        # Only the probe k-mers are ever looked up in source_counts, so count
        # just those (canonical, both strands) instead of every k-mer in the
        # source genome. Identical lookups, without the per-position Seq() calls.
        source_counts = count_self_genome_kmer_hits(source_seqs, k, unique_kmers)
        print(f"  Counted source-genome self-matches for {len(source_counts)} probe k-mers (self-match exclusion)")

    # Microbiome probe-input: each probe's own genome (matched by probe id ->
    # contig prefix) is excluded per-probe, so only OTHER genomes count as
    # off-target. A single global source set can't express this because a
    # k-mer that is self for one probe is off-target for another.
    per_probe_source = {}
    if genome_self_match and chrom_files:
        from probe_classifier import same_genome
        own_seqs = {p.id: [] for p in probes}
        for chrom_file in chrom_files:
            for rec in SeqIO.parse(chrom_file, "fasta"):
                for pid in own_seqs:
                    if same_genome(pid, rec.id):
                        own_seqs[pid].append(str(rec.seq))
                        break
        for p in probes:
            per_probe_source[p.id] = count_self_genome_kmer_hits(
                own_seqs.get(p.id, []), k, probe_to_kmers[p.id]
            )
        matched = sum(1 for p in probes if per_probe_source[p.id])
        print(f"  Per-probe self-genome k-mer exclusion enabled ({matched}/{len(probes)} probes matched a genome)")

    unsafe_probes = set()
    probe_match_details = {}
    total_self_matches = 0
    total_other_matches = 0

    for probe_id, kmers in probe_to_kmers.items():
        src = per_probe_source.get(probe_id, source_counts)
        for kmer in kmers:
            genome_count = genome_counts.get(kmer, 0)
            source_count = src.get(kmer, 0)
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

    report_elapsed = time.time() - report_start
    print(f"Report generation took {report_elapsed:.2f}s")

    safe_probes_file = f'{output_base}/safe_probes.fa'
    safe_probes = [p for p in probes if p.id not in unsafe_probes]

    host_species = None
    host_chrom_dir = None
    if _reg.is_microbiome_species(species) and args.align_microbiome and args.align_host:
        _host_obj, _ = _reg.get_microbiome(species)
        if _host_obj:
            host_species = _host_obj.id
            host_chrom_dir = _host_obj.transcript_chunks_dir
            print(f"\nSecondary safety check: Checking safe probe k-mers against {host_species} transcripts...")

    if host_species:
        if os.path.exists(host_chrom_dir):
            host_chrom_files = [os.path.join(host_chrom_dir, f) for f in os.listdir(host_chrom_dir) if f.endswith('.fa')]

            _host_precomputed = (precomputed_group_counts is not None and
                                 host_species in precomputed_group_counts)
            if _host_precomputed:
                _hbackend = "fused GPU scan (precomputed)"
            else:
                _hbackend = f"GPU (device {args.gpu_device})" if _gpu_kmer_enabled(args, k) else "Jellyfish"
            print(f"  Counting safe-probe {k}-mers against {host_species} [{_hbackend}]...")

            safe_probe_ids = {p.id for p in safe_probes}
            safe_kmers = set()
            safe_kmer_to_probes = {}
            for probe_id in safe_probe_ids:
                for kmer in probe_to_kmers[probe_id]:
                    safe_kmers.add(kmer)
                    if kmer not in safe_kmer_to_probes:
                        safe_kmer_to_probes[kmer] = []
                    safe_kmer_to_probes[kmer].append(probe_id)

            safe_kmer_list = sorted(safe_kmers)
            temp_safe_kmer_file = os.path.join(output_base, 'temp_safe_kmers.fa')
            with open(temp_safe_kmer_file, 'w', encoding='utf-8') as f:
                for i, kmer in enumerate(safe_kmer_list):
                    f.write(f">kmer_{i}\n{kmer}\n")

            if _host_precomputed:
                host_counts = precomputed_group_counts[host_species]
            else:
                host_counts = genome_kmer_counts(
                    host_chrom_files, k, safe_kmer_list, args, temp_safe_kmer_file,
                    device=args.gpu_device, threads=threads, label=host_species)
            if host_counts is None:
                print(f"  Warning: Failed to count {host_species} k-mers, skipping host check")
                host_counts = {}
                host_species = None

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
    parser.add_argument("--microbiomes", default="", help="Comma-separated microbiome IDs to align against (defaults to --species)")
    parser.add_argument("--probe-length", type=int, default=36, help="Probe length (auto-detected if probes provided)")
    parser.add_argument("--max-mismatches", type=int, default=2, help="Max mismatches allowed in alignments")
    parser.add_argument("--max-bulges", type=int, default=0,
                        help="Max bulges (indels) allowed in alignments (default: 0, no bulges)")
    parser.add_argument("--parallelism", "-p", type=int, default=None, help="Number of parallel worker threads")
    parser.add_argument("--task-id", help="Task ID for organizing outputs")
    parser.add_argument("--kmer-length", type=int, default=16,
                       help="K-mer length for safety check (minimum: 14, default: 16)")
    parser.add_argument("--skip-annotation", action="store_true",
                       help="Skip gene name annotation of SAM files")
    parser.add_argument("--align-microbiome", action="store_true", default=False,
                   help="Align against microbiome genomes")
    parser.add_argument("--align-host", action="store_true", default=False,
                    help="Align against host transcriptome")
    parser.add_argument("--tm-range", type=str, default=None,
                    help="Tm filter range in °C, e.g. 42-47. Default depends on "
                         "input: 42-47 for gene input, 20-90 for probe input "
                         "(wide enough to keep every provided probe).")
    parser.add_argument("--gc-range", type=str, default=None,
                    help="GC content filter range in percent, e.g. 40-80. Default "
                         "depends on input: 40-80 for gene input, 20-90 for probe "
                         "input (wide enough to keep every provided probe).")
    parser.add_argument("--host-internal-mode", action="store_true", default=False,
                    help="Host-internal design: keep only probes whose every alignment is to the source transcript (rejects non-aligned and cross-transcript-aligned probes)")
    parser.add_argument("--microbe-mode", action="store_true", default=False,
                    help="Microbial probe design: treat input as microbial regardless of header. Skip source inference for host targets so every host alignment counts as off-target.")
    parser.add_argument("--gpu", choices=["auto", "on", "off"], default="auto",
                    help="GPU-accelerated k-mer counting: auto (use GPU if available), on (require GPU), off (force Jellyfish). Default: auto")
    parser.add_argument("--gpu-device", type=int, default=0,
                    help="CUDA device index for single-GPU k-mer counting (default: 0)")
    parser.add_argument("--gpu-devices", default="all",
                    help="CUDA devices for alignment: 'all' to shard across every GPU "
                         "(default), or a comma-separated list e.g. '0,1' or '0' to pin.")

    args = parser.parse_args()
    if args.probe_length < 14 or args.probe_length > 50:
        print(f"Error: Probe length must be between 14 and 50 bp, got {args.probe_length}.")
        sys.exit(1)
    if args.kmer_length < 14:
        print(f"Error: K-mer length {args.kmer_length} is too small. Minimum k-mer length is 14.")
        sys.exit(1)
    if args.max_bulges < 0:
        print(f"Error: Max bulges must be >= 0, got {args.max_bulges}.")
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
        probe_lengths = [len(seq) for _, seq in sequences]
        min_len, max_len = min(probe_lengths), max(probe_lengths)

        if min_len < 14 or max_len > 50:
            print(f"Error: Probe lengths must be between 14 and 50 bp, got {min_len}-{max_len} bp.")
            return False
        if args.kmer_length > min_len:
            print(f"Error: K-mer length ({args.kmer_length}) exceeds the shortest probe ({min_len} bp). Lower --kmer-length to at most {min_len} bp and try again.")
            return False

        args.probe_length = min_len
        if min_len == max_len:
            print(f"Detected probe length: {min_len} bp")
        else:
            print(f"Mixed probe lengths detected: {min_len}-{max_len} bp (aligner identity threshold set from the shortest probe; each probe is still filtered to <= {args.max_mismatches} mismatches)")
        
    from organism_registry import load_registry as _load_reg
    _registry = _load_reg()
    if _registry.is_microbiome_species(args.species):
        genome = None
        print(f"Using pre-existing {args.species} data")
    else:
        _host = _registry.get_host(args.species)
        if _host is None:
            print(f"Unknown species '{args.species}'. Add it to config/organisms.yml.")
            return False
        print(f"\nValidating genome for {args.species}...")
        species_data_dir = _host.genome_dir
        try:
            success = validate_local_gencode(args.species, base_dir='.')

            if not success:
                print(f"Genome validation failed. Ensure files are in: {species_data_dir}")
                return False

            local_genome_release = get_local_release_number(
                _host.genome_file_regex, species_data_dir
            )
            pat = re.compile(_host.genome_file_regex)
            genome_file = None
            for fn in os.listdir(species_data_dir):
                m = pat.match(fn)
                if m and int(m.group(1)) == local_genome_release:
                    genome_file = os.path.join(species_data_dir, fn)
                    break
            if not genome_file or not os.path.exists(genome_file):
                print(f"Genome file for release {local_genome_release} not found in {species_data_dir}")
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
        output_base = f'output/alignments/{args.task_id}'
    else:
        output_base = 'output'
    os.makedirs(output_base, exist_ok=True)

    probes_out = f'{output_base}/candidate_probes.fa'
    sam_out = f'{output_base}/probe_alignments.sam'
    filtered_sam = f'{output_base}/filtered_probe_alignments.sam'
    annotated_sam = f'{output_base}/filtered_probe_alignments_annotated.sam'
    align_dir = f'{output_base}/chroms'

    if _registry.is_microbiome_species(args.species):
        _mb_host, _ = _registry.get_microbiome(args.species)
        host_dir = _mb_host.transcript_chunks_dir
        dirs = []
        if args.align_microbiome:
            for _mb_id in _selected_microbiome_ids(args):
                _, _mb = _registry.get_microbiome(_mb_id)
                if _mb and os.path.exists(_mb.data_dir):
                    dirs.append(_mb.data_dir)
        if args.align_host and os.path.exists(host_dir):
            dirs.append(host_dir)

        if not dirs:
            print("Error: At least one alignment target must be selected for microbiome species")
            return False

        chrom_dir = ':'.join(dirs)
    else:
        chrom_dir = _host.transcript_chunks_dir

    os.makedirs(align_dir, exist_ok=True)

    # Groups for the fused single-scan aligner: same (label -> files) partition
    # the k-mer safety check uses as db_groups (+ host), so its per-group counts
    # can be computed during the alignment pass instead of a second scan.
    scan_groups = []
    if _registry.is_microbiome_species(args.species):
        _grp_host, _ = _registry.get_microbiome(args.species)
        if args.align_microbiome:
            for _mb_id in _selected_microbiome_ids(args):
                _, _mb = _registry.get_microbiome(_mb_id)
                if _mb and os.path.exists(_mb.data_dir):
                    _gf = [os.path.join(_mb.data_dir, f) for f in os.listdir(_mb.data_dir)
                           if f.endswith('.fa') or f.endswith('.fna')]
                    if _gf:
                        scan_groups.append((_mb_id, _gf))
        if args.align_host and _grp_host and os.path.exists(_grp_host.transcript_chunks_dir):
            _gf = [os.path.join(_grp_host.transcript_chunks_dir, f)
                   for f in os.listdir(_grp_host.transcript_chunks_dir) if f.endswith('.fa')]
            if _gf:
                scan_groups.append((_grp_host.id, _gf))
    else:
        _gf = [os.path.join(chrom_dir, f) for f in os.listdir(chrom_dir)
               if f.endswith('.fa') or f.endswith('.fna')]
        if _gf:
            scan_groups.append((args.species, _gf))

    # Read reference from the local-NVMe mirror when staged (PROBESET_LOCAL_GENOME_ROOT),
    # avoiding slow cold NFS reads. Falls back to the original path per-file.
    if genome_stage.local_root():
        scan_groups = [(label, genome_stage.map_files(files)) for label, files in scan_groups]
        _staged = sum(1 for _, fs in scan_groups for f in fs if f.startswith(genome_stage.local_root()))
        _total = sum(len(fs) for _, fs in scan_groups)
        print(f"Local genome staging active ({genome_stage.local_root()}): "
              f"{_staged}/{_total} reference files served from local NVMe")

    # Range defaults depend on the input mode. Probes supplied by the user are
    # kept as given, so an unspecified range widens to 20-90 rather than the
    # design defaults, which would silently discard probes the user chose. An
    # explicitly supplied range is always honoured, in either mode.
    if args.tm_range is None:
        args.tm_range = "20-90" if skip_probe_generation else "42-47"
    if args.gc_range is None:
        args.gc_range = "20-90" if skip_probe_generation else "40-80"
    if skip_probe_generation:
        print(f"Filter ranges - Tm {args.tm_range}°C, GC {args.gc_range}%")

    gc_parts = args.gc_range.split('-')
    if len(gc_parts) == 1:
        gc_min = gc_max = float(gc_parts[0])
    else:
        gc_min, gc_max = float(gc_parts[0]), float(gc_parts[1])
    gc_min = max(0.0, gc_min)
    gc_max = min(100.0, gc_max)

    if skip_probe_generation:
        print("\nSkipping probe generation - using provided probes")

        total, passing = count_and_filter_gc(input_fasta, min_gc=gc_min, max_gc=gc_max)

        SeqIO.write(passing, probes_out, "fasta")

        length_desc = f"{min_len}" if min_len == max_len else f"{min_len}-{max_len}"
        print(f"Input probes: {total}")
        print(f"Found {len(passing)} candidate probes passing GC filter ({gc_min:g}-{gc_max:g}%, length={length_desc} bp).")

        if len(passing) == 0:
            print("No probes passed GC filter")
            return False
    else:
        print("\nGenerating probes from input sequence...")
        gen_cmd = [
            sys.executable, os.path.join('core', 'generate_probes.py'),
            '--probe-length', str(args.probe_length),
            '--min-gc', str(gc_min),
            '--max-gc', str(gc_max),
            '--input', input_fasta,
            '--output', probes_out
        ]
        if not run(gen_cmd, 'Generate candidate probes', check_output_file=probes_out):
            return False

    tm_parts = args.tm_range.split('-')
    if len(tm_parts) == 1:
        tm_min = tm_max = float(tm_parts[0])
    else:
        tm_min, tm_max = float(tm_parts[0]), float(tm_parts[1])
    tm_min = max(10.0, tm_min)
    tm_max = min(100.0, tm_max)
    metrics = ProbeMetricsCalculator()
    candidates = list(SeqIO.parse(probes_out, 'fasta'))
    passed = []
    tm_rejected = 0
    homopolymer_rejected = 0
    for record in candidates:
        seq = str(record.seq)
        tm = metrics.calculate_tm(seq)
        if not (tm_min <= tm <= tm_max):
            tm_rejected += 1
            continue
        if metrics.check_homopolymer_runs(seq):
            homopolymer_rejected += 1
            # There is no user-facing setting for homopolymer runs, so a
            # provided probe is never dropped for one — it is only reported.
            if not skip_probe_generation:
                continue
        passed.append(record)
    print(f"\nTm filter ({tm_min}-{tm_max}°C): rejected {tm_rejected}/{len(candidates)} probes")
    if skip_probe_generation:
        print(f"Homopolymer (≥{metrics.max_homopolymer}bp runs): {homopolymer_rejected}/{len(candidates)} probes "
              f"(reported only, not filtered)")
        print(f"Retained {len(passed)}/{len(candidates)} probes")
    else:
        print(f"Homopolymer filter (≥{metrics.max_homopolymer}bp runs): rejected {homopolymer_rejected}/{len(candidates)} probes")
        print(f"Passed both filters: {len(passed)}/{len(candidates)} probes")
    if len(passed) == 0:
        print("No probes passed Tm and homopolymer filters")
        return False
    with open(probes_out, 'w') as f:
        SeqIO.write(passed, f, 'fasta')

    identity_percent = int(((args.probe_length - args.max_mismatches) / args.probe_length) * 100)
    no_gaps = args.max_bulges == 0
    print("\nAlignment parameters:")
    if skip_probe_generation and min_len != max_len:
        print(f"Probe length: {min_len}-{max_len} bp (identity from shortest {min_len} bp)")
    else:
        print(f"Probe length: {args.probe_length} bp")
    print(f"Max mismatches (total budget, substitutions + bulges): {args.max_mismatches}")
    print(f"Max bulges: {args.max_bulges}")
    print(f"Bulges (indels): {'disabled (Hamming distance)' if no_gaps else 'enabled'}")
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

    if genome_stage.local_root():
        existing_chroms = genome_stage.map_files(existing_chroms)

    print(f"\nUsing {len(existing_chroms)} pre-chunked transcript files from {chrom_dir}/")
    chrom_files = sorted(existing_chroms)

    # GPU-only alignment (razers3 removed). The GPU aligner handles both the
    # ungapped (Hamming) and gapped (bulges) cases; there is no CPU fallback, so
    # any failure aborts the pipeline rather than silently changing the method.
    if getattr(args, 'gpu', 'auto') == 'off':
        print("Error: --gpu off, but this pipeline is GPU-only (razers3 removed). "
              "Enable the GPU to run alignment.")
        return False

    from gpu_aligner import (gpu_align_multi, gpu_align_and_count_multi,
                             available_devices,
                             gpu_available as _gpu_aln_available)
    from gpu_kmer import MAX_GPU_K as _MAX_GPU_K
    if not _gpu_aln_available():
        print("Error: no usable GPU/CuPy detected. This pipeline is GPU-only "
              "(razers3 removed) and cannot run alignment without a GPU.")
        return False

    if str(getattr(args, 'gpu_devices', 'all')).lower() == 'all':
        gpu_devices = available_devices() or [args.gpu_device]
    else:
        gpu_devices = [int(x) for x in str(args.gpu_devices).split(',') if x.strip() != '']
    _dev_label = f"device{'s' if len(gpu_devices) > 1 else ''} {','.join(map(str, gpu_devices))}"

    gpu_aligned = False
    precomputed_group_counts = None
    alignment_start = time.time()
    can_fuse = scan_groups and args.kmer_length <= _MAX_GPU_K
    if can_fuse:
        print(f"\nAligning {len(chrom_files)} files on GPU ({_dev_label}), "
              f"fused with {args.kmer_length}-mer counting...")
        ok, _gc = gpu_align_and_count_multi(
            probes_out, scan_groups, sam_out, args.max_mismatches,
            args.kmer_length, max_bulges=args.max_bulges, devices=gpu_devices,
            workers=parallelism, per_probe_cap=None)
        if ok:
            gpu_aligned = True
            precomputed_group_counts = _gc
        else:
            print("Fused GPU scan declined; trying plain GPU alignment")
    if not gpu_aligned:
        print(f"\nAligning {len(chrom_files)} files on GPU ({_dev_label})...")
        if gpu_align_multi(probes_out, chrom_files, sam_out, args.max_mismatches,
                           max_bulges=args.max_bulges, devices=gpu_devices,
                           workers=parallelism, per_probe_cap=None):
            gpu_aligned = True

    if not gpu_aligned:
        print("Error: GPU aligner could not align with the given parameters "
              f"(probe length {args.probe_length}, max_mismatches {args.max_mismatches}, "
              f"max_bulges {args.max_bulges}). Aborting — no razers3 fallback.")
        return False

    with open(sam_out, 'r', encoding='utf-8') as _sf:
        n_aln = sum(1 for ln in _sf if not ln.startswith('@'))
    _fused = " (fused k-mer counts precomputed)" if precomputed_group_counts is not None else ""
    print(f"GPU alignment complete in {time.time()-alignment_start:.2f}s "
          f"({n_aln:,} alignments){_fused} -> {sam_out}")

    print(f"\nFiltering merged SAM -> keeping alignments with <= {args.max_mismatches} total error(s) (substitutions + bulges), of which <= {args.max_bulges} may be bulge(s)")
    md_pattern = re.compile(r'MD:Z:(\S+)')
    total = 0
    kept = 0
    skipped_no_nm = 0
    aligned_probes = set()
    filter_seen = set()
    non_aligned_probes_fa = f'{output_base}/non_aligned_probes.fa'

    
    strict_host_reject = _registry.is_microbiome_species(args.species) and args.align_host
    
    _all_host_prefixes = tuple(
        p for h in _registry.hosts for p in h.transcript_id_prefixes
    )
    host_aligned_probes = set()

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
            if strict_host_reject:
                rname = parts[2]
                if _all_host_prefixes and rname.split('.', 1)[0].startswith(_all_host_prefixes):
                    host_aligned_probes.add(qname)

            md_match = md_pattern.search(line)
            if not md_match:
                skipped_no_nm += 1
                continue
            cigar = parts[5]
            subs, bulges = compute_alignment_diffs(cigar, md_match.group(1))

            if (subs + bulges) <= args.max_mismatches and bulges <= args.max_bulges:

                try:
                    flag = int(parts[1])
                except ValueError:
                    flag = 0
                strand = '-' if (flag & 0x10) else '+'
                dedup_key = (qname, parts[2], parts[3], strand)
                if dedup_key in filter_seen:
                    continue
                filter_seen.add(dedup_key)
                outf.write(line)
                kept += 1
    
    if _registry.is_microbiome_species(args.species):
        annotate_microbiome_sam(filtered_sam, _selected_microbiome_ids(args), output_base, args.align_host)

    print("\nExtracting non-aligned probes...")
    non_aligned_count = 0
    candidate_probes = list(SeqIO.parse(probes_out, "fasta"))

    is_microbiome = _registry.is_microbiome_species(args.species)
    if annotator and is_microbiome and args.align_host:
        _h = _registry.host_for_microbiome(args.species)
        host_species = _h.id if _h else None

        temp_sam = filtered_sam.replace('filtered_probe_alignments.sam', 'filtered_probe_alignments_temp.sam')
        annotator.annotate_sam(annotated_sam, temp_sam, host_species)
        shutil.move(temp_sam, annotated_sam)

    elif annotator and not is_microbiome:
        annotator.annotate_sam(filtered_sam, annotated_sam, args.species)

    sam_to_analyze = annotated_sam if os.path.exists(annotated_sam) else filtered_sam

    
    gene_mappings = annotator.mappings if annotator and annotator.mappings else None
    genome_self_match = False
    if scenario == "probe_sequence" and not args.host_internal_mode:
        source_transcripts = set()
        source_gene = None
        if is_microbiome:
            genome_self_match = True
            print("\nProbe input (microbiome): self-alignments identified by probe id <-> genome id match; "
                  "cross-genome alignments treated as off-target")
        else:
            print("\nSkipping source inference (probe input mode): all alignments treated as off-target")
    elif args.microbe_mode and not is_microbiome:
        
        print("\nMicrobe mode: skipping source inference; all host alignments treated as off-target")
        source_transcripts = set()
        source_gene = None
    elif not is_microbiome:
        
        header = sequences[0][0] if sequences else ''
        source_info = resolve_source_from_header(header, gene_mappings)
        source_transcripts = source_info['source_transcripts']
        source_gene = source_info['source_gene']
    else:
        source_info = infer_source_transcripts(filtered_sam, gene_mappings=gene_mappings, is_microbiome=is_microbiome)
        source_transcripts = source_info['source_transcripts']
        source_gene = source_info['source_gene']

    if args.host_internal_mode:
        safe_ids, _details = filter_probes_by_self_alignment(
            sam_to_analyze, source_transcripts, source_gene,
            gene_mappings=gene_mappings,
        )
        probe_classifications = {pid: {'status': 'safe'} for pid in safe_ids}
        print(f"Host-internal filter: {len(safe_ids)} probes align only to source")

        
        if not source_transcripts and safe_ids:
            inferred_transcripts = set()
            inferred_genes = set()
            for pid in safe_ids:
                d = _details.get(pid) or {}
                inferred_transcripts.update(d.get('aligned_to', []))
                inferred_genes.update(d.get('genes', []))
            source_transcripts = inferred_transcripts
            if not source_gene and len(inferred_genes) == 1:
                source_gene = next(iter(inferred_genes))
            print(f"Host-internal filter: inferred source from safe probes: "
                  f"{source_gene or '(no single gene)'} "
                  f"({len(source_transcripts)} transcript(s))")
    else:
        probe_classifications = classify_probes_from_sam(
            sam_to_analyze, is_microbiome=is_microbiome,
            source_transcripts=source_transcripts, source_gene=source_gene,
            is_probe_input=(scenario == "probe_sequence"),
            genome_self_match=genome_self_match,
        )

    
    try:
        import json as _json
        with open(f'{output_base}/source_info.json', 'w', encoding='utf-8') as _sf:
            _json.dump({
                'source_gene': source_gene,
                'source_transcripts': sorted(source_transcripts),
                'is_microbiome': bool(is_microbiome),
                'scenario': scenario,
                'host_internal_mode': bool(args.host_internal_mode),
                'genome_self_match': bool(genome_self_match),
            }, _sf)
    except OSError:
        pass

    # Count self-aligning vs off-target probes
    self_aligning_probes = sum(1 for p_id, data in probe_classifications.items() if data['status'] == 'safe')
    off_target_probes = len(aligned_probes) - self_aligning_probes

    # Extract all on-target probes (self-aligned safe + non-aligned) for k-mer check.
    # K-mer check is mandatory for all on-target probes: if any k-mer matches
    # a gene other than the source, the probe is not safe.
    print("\nExtracting probes...")
    self_aligned_safe_count = 0

    host_rejected_count = 0
    with open(non_aligned_probes_fa, 'w', encoding='utf-8') as out:
        for record in candidate_probes:
            probe_id = record.id
            
            if strict_host_reject and probe_id in host_aligned_probes:
                host_rejected_count += 1
                continue
            if probe_id not in probe_classifications:
                if args.host_internal_mode:
                    continue
                SeqIO.write(record, out, "fasta")
                non_aligned_count += 1
            elif probe_classifications[probe_id]['status'] == 'safe':
                
                if strict_host_reject and not args.align_microbiome:
                    continue
                SeqIO.write(record, out, "fasta")
                self_aligned_safe_count += 1

    print("\nFiltering summary:")
    print(f"Total alignments scanned: {total:,}")
    print(f"Kept in filtered (edit distance <= {args.max_mismatches}): {kept:,}")
    print(f"Total candidate probes: {len(candidate_probes):,}")
    print(f"Aligned probes found: {len(aligned_probes):,}")
    print(f"  - Self-alignments only (safe): {self_aligning_probes:,}")
    print(f"  - Off-target alignments: {off_target_probes:,}")
    if source_gene:
        print(f"  - Inferred source gene: {source_gene}")
    print(f"Self-aligned safe probes: {self_aligned_safe_count:,}")
    print(f"Non-aligned probes: {non_aligned_count:,}")
    if strict_host_reject:
        print(f"Rejected for host cross-reactivity: {host_rejected_count:,}")
    print(f"Skipped (no MD tag): {skipped_no_nm:,}")

    total_on_target = self_aligned_safe_count + non_aligned_count
    safe_probes_file = f'{output_base}/safe_probes.fa'
    safe_scores_csv = f'{output_base}/safe_probes_scores.csv'

    if total_on_target > 0:
        print(f"\nPerforming {args.kmer_length}-mer safety check on {total_on_target:,} on-target probes...")
        print(f"  ({self_aligned_safe_count:,} self-aligned + {non_aligned_count:,} non-aligned)")
        kmer_result_file = check_14mer_safety(
            non_aligned_probes_fa, args.species, output_base, args,
            source_transcripts=source_transcripts,
            genome_self_match=genome_self_match,
            precomputed_group_counts=precomputed_group_counts,
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
        print(f"\nMeasuring {total_safe_count} safe probes...")
        write_probe_metrics_report(safe_probes_file, safe_scores_csv)
    else:
        print("\nNo safe probes to measure - creating empty metrics file...")
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
    os.makedirs(DEFAULT_DIRS['gene'], exist_ok=True)

    main()