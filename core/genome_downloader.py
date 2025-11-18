"""Download, process, filter and chunk Gencode genome data for human and mouse species."""
import os
import re
import subprocess
import requests
from lxml import html


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


def get_latest_release_number(species):
    """Scrape Gencode FTP directory for latest release number."""
    base_url = "https://ftp.ebi.ac.uk/pub/databases/gencode/"
    ftp_url = f"{base_url}Gencode_{species}/"
    try:
        response = requests.get(ftp_url, timeout=30)
        if response.status_code != 200:
            print(f"Error: Failed to access FTP directory for {species}. Status code: {response.status_code}")
            return None
        
        tree = html.fromstring(response.content)
        directories = tree.xpath('//a/text()')
        pattern = r'release_(\d+)' if species == "human" else r'release_M(\d+)'
        release_numbers = [int(match.group(1)) for dir_name in directories if (match := re.match(pattern, dir_name.strip('/')))]
        
        return max(release_numbers) if release_numbers else None
    
    except requests.exceptions.RequestException as e:
        print(f"Error accessing FTP directory for {species}: {e}")
        return None


def download_file(url, save_path):
    """Download file from URL and save locally."""
    try:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        response = requests.get(url, stream=True, timeout=30)
        if response.status_code != 200:
            print(f"Error: Failed to download file from {url}. Status code: {response.status_code}")
            return False
        
        with open(save_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)
        print(f"Downloaded: {save_path}")
        return True
    
    except requests.exceptions.RequestException as e:
        print(f"Error downloading {url}: {e}")
        return False


def decompress_file(gz_path):
    """Decompress .gz file using gunzip."""
    try:
        subprocess.run(f"gunzip -f {gz_path}", shell=True, check=True)
        print(f"Decompressed: {gz_path}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error decompressing {gz_path}: Command '{e.cmd}' returned non-zero exit status {e.returncode}.")
        if e.stderr:
            print(e.stderr)
        return False
    except FileNotFoundError:
        print("Error: 'gunzip' command not found. Ensure gzip is installed.")
        return False


def filter_transcripts(input_file, output_file, min_length=50):
    """Filter FASTA sequences by minimum length."""
    kept_count = 0
    filtered_count = 0
    
    current_header = None
    current_sequence = ""
    
    print(f"Filtering transcripts (minimum length: {min_length}nt)...")
    
    with open(input_file, 'r') as inf, open(output_file, 'w') as outf:
        for line in inf:
            line = line.strip()
            
            if line.startswith('>'):
                if current_header is not None:
                    if len(current_sequence) >= min_length:
                        outf.write(f"{current_header}\n")
                        for i in range(0, len(current_sequence), 70):
                            outf.write(f"{current_sequence[i:i+70]}\n")
                        kept_count += 1
                    else:
                        filtered_count += 1
                
                current_header = line
                current_sequence = ""
            else:
                current_sequence += line
        
        if current_header is not None:
            if len(current_sequence) >= min_length:
                outf.write(f"{current_header}\n")
                for i in range(0, len(current_sequence), 70):
                    outf.write(f"{current_sequence[i:i+70]}\n")
                kept_count += 1
            else:
                filtered_count += 1
    
    print(f"Filtering complete: kept {kept_count}, removed {filtered_count}")
    return kept_count, filtered_count


def chunk_filtered_transcripts(input_file, release, species, max_headers=3000):
    """Split filtered transcripts into chunks with max_headers per file."""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_dir = os.path.join(script_dir, '..', 'data', 'gencode_data', species, 'transcript_chunks')
    os.makedirs(output_dir, exist_ok=True)
    
    existing_chunks = []
    if os.path.exists(output_dir):
        existing_chunks = [f for f in os.listdir(output_dir) 
                          if f.endswith('.fa') and f'release_{release}' in f]
    
    if existing_chunks:
        print(f"Chunks already exist for release_{release} ({len(existing_chunks)} files)")
        return len(existing_chunks)
    
    print(f"Chunking filtered transcripts ({max_headers} headers per chunk)...")
    part_num = 1
    header_count = 0
    out_file = None
    chunks_created = 0
    
    base_name = os.path.basename(input_file).replace('.fa', '')
    
    with open(input_file, 'r') as f:
        for line in f:
            if line.startswith(">"):
                header_count += 1
                
                if out_file is None or header_count > max_headers:
                    if out_file:
                        out_file.close()
                        chunks_created += 1
                    
                    out_path = os.path.join(output_dir, f"{base_name}.part{part_num}.fa")
                    out_file = open(out_path, 'w')
                    part_num += 1
                    header_count = 1
            
            if out_file:
                out_file.write(line)
    
    if out_file:
        out_file.close()
        chunks_created += 1
    
    print(f"Created {chunks_created} chunk files in {output_dir}/")
    return chunks_created


def cleanup_all_outdated_files(species_data_dir, latest_release, genome_file, annotation_file, transcripts_file, filtered_file):
    """Remove ALL files from older releases including transcripts, filtered, and chunks."""
    print(f"Cleaning up ALL outdated files (keeping only release_{latest_release})")
    files_deleted = 0
    
    all_source_files_current = True
    
    for file_path in [genome_file, annotation_file, transcripts_file, filtered_file]:
        if os.path.exists(file_path):
            release_match = re.search(r'release_(\d+)', file_path)
            if release_match:
                file_release = int(release_match.group(1))
                if file_release < latest_release:
                    all_source_files_current = False
                    print(f"Found outdated: {os.path.basename(file_path)} (release_{file_release})")
            else:
                all_source_files_current = False
        else:
            all_source_files_current = False
            print(f"Missing: {os.path.basename(file_path)}")
    
    if os.path.exists(species_data_dir):
        for filename in os.listdir(species_data_dir):
            file_path = os.path.join(species_data_dir, filename)
            
            if os.path.isdir(file_path):
                continue
            
            should_delete = False
            
            if filename.endswith('.fai'):
                should_delete = True
            else:
                patterns = [
                    r'release_(\d+)',
                    r'\.v(\d+)\.',
                    r'\.vM(\d+)\.',
                ]
                
                for pattern in patterns:
                    match = re.search(pattern, filename)
                    if match:
                        file_release = int(match.group(1))
                        if file_release < latest_release:
                            should_delete = True
                            print(f"  Deleting outdated: {filename} (release_{file_release})")
                        break
            
            if should_delete:
                try:
                    os.remove(file_path)
                    files_deleted += 1
                except (OSError, FileNotFoundError) as e:
                    print(f"Warning: Could not delete {filename}: {e}")
        
        chunk_dir = os.path.join(species_data_dir, "transcript_chunks")
        if os.path.exists(chunk_dir):
            if not all_source_files_current:
                print("Source files are outdated/missing - deleting ALL chunks to force regeneration")
                for filename in os.listdir(chunk_dir):
                    if filename.endswith('.fa'):
                        file_path = os.path.join(chunk_dir, filename)
                        try:
                            os.remove(file_path)
                            files_deleted += 1
                        except (OSError, FileNotFoundError) as e:
                            print(f"  Warning: Could not delete {filename}: {e}")
            else:
                has_outdated_chunks = False
                for filename in os.listdir(chunk_dir):
                    if filename.endswith('.fa'):
                        chunk_release_match = re.search(r'release_(\d+)', filename)
                        chunk_release = int(chunk_release_match.group(1)) if chunk_release_match else None
                        
                        if not chunk_release or chunk_release < latest_release:
                            has_outdated_chunks = True
                            break
                
                if has_outdated_chunks:
                    print("Found outdated chunks - deleting ALL chunks to force complete regeneration")
                    for filename in os.listdir(chunk_dir):
                        if filename.endswith('.fa'):
                            file_path = os.path.join(chunk_dir, filename)
                            try:
                                os.remove(file_path)
                                files_deleted += 1
                            except (OSError, FileNotFoundError) as e:
                                print(f"  Warning: Could not delete {filename}: {e}")
    
    if files_deleted > 0:
        print(f"Cleaned up {files_deleted} outdated file(s)")
    else:
        print("No outdated files found")
    
    return all_source_files_current


def download_and_process_gencode(species, base_dir=None):
    """Main pipeline: download, generate transcripts, filter, and chunk.
    
    Args:
        species: 'human' or 'mouse'
        base_dir: Base directory for data storage. If None, uses '../data' relative to script location.
    """
    print(f"\n{'='*60}")
    print(f"GENCODE PIPELINE FOR {species.upper()}")
    print(f"{'='*60}\n")
    
    latest_release = get_latest_release_number(species)
    if not latest_release:
        print(f"Error: Could not determine latest release for {species}.")
        return False
    
    print(f"Latest release: {latest_release}\n")
    
    if base_dir is None:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        base_dir = os.path.join(script_dir, '..')
    
    species_data_dir = os.path.join(base_dir, 'data', 'gencode_data', species)
    os.makedirs(species_data_dir, exist_ok=True)
    
    if species == "human":
        genome_pattern = r'GRCh38\.p14\.genome\.release_(\d+)\.fa'
        annotation_pattern = r'gencode\.v(\d+)\.basic\.annotation\.release_(\d+)\.gtf'
        genome_url = f"https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{latest_release}/GRCh38.p14.genome.fa.gz"
        annotation_url = f"https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{latest_release}/gencode.v{latest_release}.basic.annotation.gtf.gz"
        genome_file = f"{species_data_dir}/GRCh38.p14.genome.release_{latest_release}.fa"
        annotation_file = f"{species_data_dir}/gencode.v{latest_release}.basic.annotation.release_{latest_release}.gtf"
    else:  # mouse
        genome_pattern = r'GRCm39\.genome\.release_(\d+)\.fa'
        annotation_pattern = r'gencode\.vM(\d+)\.basic\.annotation\.release_(\d+)\.gtf'
        genome_url = f"https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M{latest_release}/GRCm39.genome.fa.gz"
        annotation_url = f"https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M{latest_release}/gencode.vM{latest_release}.basic.annotation.gtf.gz"
        genome_file = f"{species_data_dir}/GRCm39.genome.release_{latest_release}.fa"
        annotation_file = f"{species_data_dir}/gencode.vM{latest_release}.basic.annotation.release_{latest_release}.gtf"
    
    transcripts_file = genome_file.replace('.fa', '_transcripts.fa')
    filtered_file = genome_file.replace('.fa', '_transcripts.filtered.fa')
    chunk_dir = f"{species_data_dir}/transcript_chunks"
    
    local_genome_release = get_local_release_number(genome_pattern, species_data_dir)
    local_annotation_release = get_local_release_number(annotation_pattern, species_data_dir)
    
    print("Checking local file versions:")
    print(f"Genome: release_{local_genome_release if local_genome_release else 'none'}")
    print(f"Annotation: release_{local_annotation_release if local_annotation_release else 'none'}")
    print(f" Latest available: release_{latest_release}\n")
    
    all_source_files_current = cleanup_all_outdated_files(
        species_data_dir, latest_release, genome_file, annotation_file, transcripts_file, filtered_file
    )
    print()
    
    genome_needs_update = (not local_genome_release or local_genome_release < latest_release or not os.path.exists(genome_file))
    annotation_needs_update = (not local_annotation_release or local_annotation_release < latest_release or not os.path.exists(annotation_file))
    needs_reprocess = genome_needs_update or annotation_needs_update
    
    if needs_reprocess:
        print("Update needed:")
        if genome_needs_update:
            if not os.path.exists(genome_file):
                print("- Genome is missing")
            else:
                print("- Genome is outdated")
        if annotation_needs_update:
            if not os.path.exists(annotation_file):
                print("- Annotation is missing")
            else:
                print("- Annotation is outdated")
        print()
    else:
        print(f"Both genome and annotation are up to date (release_{latest_release})")
        
        chunk_dir = f"{species_data_dir}/transcript_chunks"
        if os.path.exists(chunk_dir) and os.path.exists(genome_file) and os.path.exists(annotation_file):
            existing_chunks = [f for f in os.listdir(chunk_dir) 
                              if f.endswith('.fa') and f'release_{latest_release}' in f]
            if existing_chunks and os.path.exists(filtered_file):
                print(f"Pipeline already complete for release_{latest_release}")
                print(f"Found {len(existing_chunks)} chunk files")
                print(f" Filtered transcripts: {filtered_file}")
                return True
        
        print("Transcripts need to be regenerated")
    
    print(f"\n{'─'*60}")
    print(" STEP 1: Downloading Genome")
    print(f"{'─'*60}")
    if genome_needs_update or not os.path.exists(genome_file):
        genome_gz = genome_file + '.gz'
        
        if os.path.exists(genome_gz):
            print(f"  {genome_gz} exists, attempting decompression...")
            if decompress_file(genome_gz):
                print(f"Genome ready: {genome_file}")
            else:
                print(f"Decompression failed for {genome_gz}. Removing corrupted file to force re-download.")
                try:
                    os.remove(genome_gz)
                except OSError as e:
                    print(f"Warning: Could not remove corrupted file {genome_gz}: {e}")
                if not download_file(genome_url, genome_gz):
                    return False
                if not decompress_file(genome_gz):
                    print(f"Decompression failed again for {genome_gz}. Removing corrupted file.")
                    try:
                        os.remove(genome_gz)
                    except OSError as e:
                        print(f"Warning: Could not remove corrupted file {genome_gz}: {e}")
                    return False
        else:
            if not download_file(genome_url, genome_gz):
                return False
            if not decompress_file(genome_gz):
                print(f"Decompression failed for {genome_gz}. Removing corrupted file.")
                try:
                    os.remove(genome_gz)
                except OSError as e:
                    print(f"Warning: Could not remove corrupted file {genome_gz}: {e}")
                return False
    else:
        print(f"Genome already exists: {genome_file}")
    
    print(f"\n{'─'*60}")
    print("STEP 2: Downloading Annotation")
    print(f"{'─'*60}")
    if annotation_needs_update or not os.path.exists(annotation_file):
        annotation_gz = annotation_file + '.gz'
        
        if os.path.exists(annotation_gz):
            print(f"  {annotation_gz} exists, attempting decompression...")
            if decompress_file(annotation_gz):
                print(f"Annotation ready: {annotation_file}")
            else:
                print(f"Decompression failed for {annotation_gz}. Removing corrupted file to force re-download.")
                try:
                    os.remove(annotation_gz)
                except OSError as e:
                    print(f"Warning: Could not remove corrupted file {annotation_gz}: {e}")
                if not download_file(annotation_url, annotation_gz):
                    return False
                if not decompress_file(annotation_gz):
                    print(f"Decompression failed again for {annotation_gz}. Removing corrupted file.")
                    try:
                        os.remove(annotation_gz)
                    except OSError as e:
                        print(f"Warning: Could not remove corrupted file {annotation_gz}: {e}")
                    return False
        else:
            if not download_file(annotation_url, annotation_gz):
                return False
            if not decompress_file(annotation_gz):
                print(f"Decompression failed for {annotation_gz}. Removing corrupted file.")
                try:
                    os.remove(annotation_gz)
                except OSError as e:
                    print(f"Warning: Could not remove corrupted file {annotation_gz}: {e}")
                return False
    else:
        print(f"Annotation already exists: {annotation_file}")
    
    print(f"\n{'─'*60}")
    print("STEP 3: Generating Transcripts")
    print(f"{'─'*60}")
    
    transcripts_needs_generation = True
    if os.path.exists(transcripts_file):
        transcripts_release_match = re.search(r'release_(\d+)', transcripts_file)
        if transcripts_release_match:
            transcripts_release = int(transcripts_release_match.group(1))
            if transcripts_release == latest_release:
                print(f"Transcripts already exist: {transcripts_file}")
                transcripts_needs_generation = False
            else:
                print(f"Transcripts are outdated (release_{transcripts_release}), regenerating...")
        else:
            print("Cannot determine transcripts version, regenerating...")
    
    if transcripts_needs_generation:
        try:
            print("Running gffread...")
            subprocess.run(
                f"gffread {annotation_file} -g {genome_file} -w {transcripts_file}",
                shell=True,
                check=True
            )
            print(f"Generated: {transcripts_file}")
        except subprocess.CalledProcessError as e:
            print(f"Error running gffread: {e}")
            return False
        except FileNotFoundError:
            print("Error: 'gffread' command not found.")
            return False
    
    # Filter transcripts
    print(f"\n{'─'*60}")
    print("STEP 4: Filtering Transcripts")
    print(f"{'─'*60}")
    
    filtered_needs_generation = True
    if os.path.exists(filtered_file):
        filtered_release_match = re.search(r'release_(\d+)', filtered_file)
        if filtered_release_match:
            filtered_release = int(filtered_release_match.group(1))
            if filtered_release == latest_release:
                print(f"Filtered transcripts already exist: {filtered_file}")
                filtered_needs_generation = False
            else:
                print(f"Filtered transcripts are outdated (release_{filtered_release}), regenerating...")
        else:
            print("Cannot determine filtered version, regenerating...")
    
    if filtered_needs_generation:
        filter_transcripts(transcripts_file, filtered_file, min_length=50)
    
    print(f"\n{'─'*60}")
    print("STEP 5: Chunking Filtered Transcripts")
    print(f"{'─'*60}")
    chunk_filtered_transcripts(filtered_file, latest_release, species, max_headers=3000)
    
    print(f"\n{'='*60}")
    print("PIPELINE COMPLETE")
    print(f"{'='*60}")
    print(f"Species: {species}")
    print(f"Release: {latest_release}")
    print(f"Genome: {genome_file}")
    print(f"Annotation: {annotation_file}")
    print(f"Transcripts: {transcripts_file}")
    print(f"Filtered: {filtered_file}")
    print(f"Chunks: {chunk_dir}/")
    print(f"{'='*60}\n")
    
    return True


if __name__ == "__main__":
    import sys

    if len(sys.argv) != 2:
        print("Usage: python genome_downloader.py <species>")
        print("Species options: 'human' or 'mouse'")
        sys.exit(1)
    
    species = sys.argv[1].lower()
    if species not in ['human', 'mouse']:
        print("Error: Species must be 'human' or 'mouse'")
        sys.exit(1)
    
    success = download_and_process_gencode(species)
    sys.exit(0 if success else 1)