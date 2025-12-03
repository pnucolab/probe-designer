"""
Simple Probe Alignment Pipeline - Modified

This script handles two input scenarios:
1. Gene sequence (pasted) -> generate probes from gene -> align to species genome
2. Probe sequence (pasted) -> skip generation -> align to species genome

Usage:
  # Scenario 1: Gene sequence (pasted)
  python simple_pipeline.py --gene-sequence ">gene1\nATCG..." --species human --probe-length 30 --max-mismatches 2
  
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
from Bio import SeqIO
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'core'))
from score_probes import score_and_save_probes
from probe_classifier import classify_probes_from_sam

try:
    import genome_downloader
except Exception:
    genome_downloader = None

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
            with open(cache_file, 'r') as f:
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
        with open(cache_file, 'r') as f:
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
        
        with open(gtf_file, 'r') as f:
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
        
        with open(cache_file, 'w') as f:
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
        
        with open(input_sam, 'r') as infile, open(output_sam, 'w') as outfile:
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

def count_and_filter_gc(fasta_file, min_gc=40, max_gc=60):
    """Count probes and filter by GC content."""
    count = 0
    passing = []
    
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = str(record.seq).upper()
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

def annotate_microbiome_sam(sam_file, species, output_base):
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
        with open(sam_file, 'r', encoding='utf-8') as inf, \
             open(annotated_sam, 'w', encoding='utf-8') as outf:
            
            for line in inf:
                if line.startswith('@'):
                    outf.write(line)
                else:
                    parts = line.strip().split('\t')
                    if len(parts) >= 11:
                        rname = parts[2]
                        genome = rname.split('_')[0] if '_' in rname else rname
                        species_name = genome_to_species.get(genome, 'Unknown')
                        outf.write(line.strip() + f'\tSP:Z:{species_name}\n')
                    else:
                        outf.write(line)
        
        print(f"Created annotated SAM: {annotated_sam}")
        return True
        
    except Exception as e:
        print(f"Error annotating SAM file: {e}")
        return False
    
def main():
    """Main pipeline logic."""
    parser = argparse.ArgumentParser(description="Simple probe -> genome alignment pipeline")

    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument("--gene-sequence", "-gs", help="Gene sequence in FASTA format (pasted text)")
    input_group.add_argument("--probe-sequence", "-ps", help="Probe sequence in FASTA format (pasted text)")
    input_group.add_argument("--gene-sequence-stdin", action="store_true", help="Read gene sequence from stdin")
    input_group.add_argument("--probe-sequence-stdin", action="store_true", help="Read probe sequence from stdin")

    parser.add_argument("--species", "-s", default="human", help="Species (human or mouse)")
    parser.add_argument("--probe-length", type=int, default=30, help="Probe length (auto-detected if probes provided)")
    parser.add_argument("--max-mismatches", type=int, default=2, help="Max mismatches allowed in alignments")
    parser.add_argument("--parallelism", "-p", type=int, default=None, help="Number of parallel worker threads")
    parser.add_argument("--task-id", help="Task ID for organizing outputs")
    parser.add_argument("--skip-annotation", action="store_true",
                       help="Skip gene name annotation of SAM files")
    
    args = parser.parse_args()
    if args.gene_sequence_stdin:
        args.gene_sequence = sys.stdin.read()
    elif args.probe_sequence_stdin:
        args.probe_sequence = sys.stdin.read()
    print("SIMPLE PROBE -> GENOME ALIGNMENT PIPELINE")
    print("=" * 80)
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    # Initialize annotator
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
    
    if skip_probe_generation:
        detected_length = detect_probe_length(sequences)
        if detected_length is None:
            print("Inconsistent probe lengths detected. All probes must have the same length.")
            return False
        args.probe_length = detected_length
        print(f"Probe length auto-detected: {args.probe_length} bp")
        
        for header, seq in sequences:
            if len(seq) != args.probe_length:
                print(f"Probe length mismatch: {header} has {len(seq)} bp, expected {args.probe_length} bp")
                return False
    if args.species in ["gut-microbe", "human-oral-microbiome", "human-skin-microbiome", "human-vaginal-microbiome", "mouse-gut-microbiome"]:
        genome = None
        print(f"Using pre-existing {args.species} data")
    else:
        if not genome_downloader:
            print("genome_downloader module not available. Ensure core/genome_downloader.py is present.")
            return False

        print(f"\nChecking/downloading genome for {args.species}...")
        try:
            success = genome_downloader.download_and_process_gencode(args.species, base_dir='.')
            
            if not success:
                print("Genome download/processing failed")
                return False
            
            species_data_dir = os.path.join('data', 'gencode_data', args.species)
            
            latest_release = genome_downloader.get_latest_release_number(args.species)
            if not latest_release:
                print("Could not determine latest release")
                return False
            
            if args.species == "human":
                genome_file = os.path.join(species_data_dir, f'GRCh38.p14.genome.release_{latest_release}.fa')
            else:  # mouse
                genome_file = os.path.join(species_data_dir, f'GRCm39.genome.release_{latest_release}.fa')
            
            if not os.path.exists(genome_file):
                print(f"Genome file not found: {genome_file}")
                return False
            
            genome = genome_file
            print(f"Genome ready: {genome}")
        except (FileNotFoundError, PermissionError, ConnectionError, TimeoutError, OSError) as e:
            print(f"Error downloading genome: {e}")
            traceback.print_exc()
            return False

    # Pre-load annotation mappings if enabled
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
        chrom_dir = os.path.join('data', args.species)
    else:
        chrom_dir = os.path.join('data', 'gencode_data', args.species, 'transcript_chunks')

    os.makedirs(align_dir, exist_ok=True)

    if skip_probe_generation:
        print("\nSkipping probe generation - using provided probes")

        total, passing = count_and_filter_gc(input_fasta, min_gc=40, max_gc=60)
        
        SeqIO.write(passing, probes_out, "fasta")
        
        print(f"Found {len(passing)} candidate probes passing GC filter (40-60%, length={args.probe_length}).")
        
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

    if not os.path.exists(chrom_dir):
        print(f"Chromosome chunks directory not found: {chrom_dir}")
        print("Chromosomes should have been created during genome download")
        return False

    existing_chroms = [os.path.join(chrom_dir, f) for f in os.listdir(chrom_dir) if f.endswith('.fa') or f.endswith('.fna')]
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
        annotate_microbiome_sam(filtered_sam, args.species, output_base)

    print("\nExtracting non-aligned probes...")
    non_aligned_count = 0
    candidate_probes = list(SeqIO.parse(probes_out, "fasta"))

    is_microbiome = args.species in microbiome_species
    if annotator and not is_microbiome:
        annotator.annotate_sam(filtered_sam, annotated_sam, args.species)
    
    sam_to_analyze = annotated_sam if os.path.exists(annotated_sam) else filtered_sam
    probe_classifications = classify_probes_from_sam(sam_to_analyze, is_microbiome=is_microbiome)
    
    print("\nExtracting non-aligned probes...")
    
    with open(non_aligned_probes_fa, 'w', encoding='utf-8') as non_aligned_out:
        for record in candidate_probes:
            probe_id = record.id.split('|')[0]
            if probe_id not in probe_classifications:
                SeqIO.write(record, non_aligned_out, "fasta")
                non_aligned_count += 1
            elif probe_classifications[probe_id]['status'] == 'safe':
                SeqIO.write(record, non_aligned_out, "fasta")
                non_aligned_count += 1

    print("\nFiltering summary:")
    print(f"Total alignments scanned: {total:,}")
    print(f"Kept in filtered (NM<={args.max_mismatches}): {kept:,}")
    print(f"Total candidate probes: {len(candidate_probes):,}")
    print(f"Aligned probes found: {len(aligned_probes):,}")
    print(f"Non-aligned probes: {non_aligned_count:,}")
    print(f"Skipped (no NM tag): {skipped_no_nm:,}")
    
    if non_aligned_count > 0:
        print("\nScoring {non_aligned_count} non-aligned probes...")
        non_aligned_scores_csv = non_aligned_probes_fa.replace('.fa', '_scores.csv')
        score_and_save_probes(non_aligned_probes_fa, non_aligned_scores_csv)
    else:
        print("\nNo non-aligned probes to score")

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