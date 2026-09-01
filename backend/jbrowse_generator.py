"""
JBrowse 2 File Generator for Probe Visualization

Generates GFF3 and indexed FASTA files for visualizing probe locations
on the input gene sequence, colored by off-target risk from transcriptome alignment.
"""

import re
from pathlib import Path
from typing import List, Dict
import logging
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'core'))
from probe_classifier import classify_probes_from_sam, infer_source_transcripts


logger = logging.getLogger(__name__)


class JBrowseFileGenerator:
    """Generate JBrowse-compatible files from probe analysis results"""
    
    def __init__(self, output_dir: str, job_id: str):
        self.output_dir = Path(output_dir)
        self.job_id = job_id
        self.output_dir.mkdir(parents=True, exist_ok=True)
    
    def extract_input_sequence(
        self, 
        gene_sequences_dir: str,
        probe_sequences_dir: str,
        input_type: str
    ) -> tuple[str, str]:
        
        if input_type == 'gene_sequence':
            search_dir = Path(gene_sequences_dir)
            prefix = "pasted_gene"
        elif input_type == 'probe_sequence':
            search_dir = Path(probe_sequences_dir)
            prefix = "pasted_probes"
        else:
            raise ValueError(f"Unknown input_type: {input_type}")
        pattern = f"{self.job_id}_{prefix}.fasta"
        input_file = search_dir / pattern
        
        if not input_file.exists():
            raise FileNotFoundError(f"Input file not found: {input_file}")
        
        sequence_id = None
        sequence_lines = []
        
        with open(input_file, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                
                if line.startswith('>'):
                    if sequence_id is not None:
                        break
                    sequence_id = line[1:].split()[0] 
                else:
                    sequence_lines.append(line.upper())
        
        if not sequence_id or not sequence_lines:
            raise ValueError("No valid sequence found in input file")
        
        sequence = ''.join(sequence_lines)
        
        return sequence_id, sequence
    
    def generate_reference_fasta(self, sequence: str, sequence_id: str = "reference"):
        """
        Generate reference FASTA file and index it.
        
        Args:
            sequence: The reference sequence string
            sequence_id: ID for the sequence (default: "reference")
            
        Returns:
            Path to the generated FASTA file
        """
        fasta_path = self.output_dir / "reference.fasta"
        with open(fasta_path, 'w', encoding='utf-8') as f:
            f.write(f">{sequence_id}\n")
            for i in range(0, len(sequence), 80):
                f.write(sequence[i:i+80] + "\n")
    
        self._create_fasta_index(fasta_path, sequence_id, len(sequence))
        
        logger.info(f"Generated: {fasta_path}")
        logger.info(f"Generated: {fasta_path}.fai")
        
        return fasta_path
    
    def _create_fasta_index(self, fasta_path: Path, seq_id: str, seq_length: int):
        """
        Create a FASTA index file (.fai) manually.
        
        Format: NAME LENGTH OFFSET LINEBASES LINEWIDTH
        """
        fai_path = Path(str(fasta_path) + '.fai')
        offset = len(f">{seq_id}\n")
        linebases = 80
        linewidth = 81  
        
        with open(fai_path, 'w', encoding='utf-8') as f:
            f.write(f"{seq_id}\t{seq_length}\t{offset}\t{linebases}\t{linewidth}\n")
    
    def parse_probe_regions_from_fasta(
        self, 
        candidate_probes_file: Path,
        aligned_sam_file: Path
    ) -> List[Dict]:
        """Parse probe regions and classify risk using shared classifier."""
        sys.path.insert(0, str(Path(__file__).parent))
        
        # Detect if this is microbiome data by checking for SP:Z: tags
        is_microbiome = False
        if aligned_sam_file.exists():
            with open(aligned_sam_file, 'r', encoding='utf-8') as f:
                for line in f:
                    if line.startswith('@'):
                        continue
                    if 'SP:Z:' in line and 'SP:Z:Unknown' not in line:
                        is_microbiome = True
                        break
        
        # Source-of-truth: prefer the pipeline's persisted decision (source_info.json
        # written by probe_designer.py). Fall back to local alignment-pattern inference
        # only when that file is absent (older jobs).
        probe_classifications = {}
        host_internal_mode = False
        genome_self_match = False
        if aligned_sam_file.exists():
            source_info_file = self.output_dir / "source_info.json"
            if source_info_file.exists():
                import json as _json
                with open(source_info_file, 'r', encoding='utf-8') as _sf:
                    _src = _json.load(_sf)
                source_transcripts = set(_src.get('source_transcripts') or [])
                source_gene = _src.get('source_gene')
                host_internal_mode = bool(_src.get('host_internal_mode'))
                genome_self_match = bool(_src.get('genome_self_match'))
            else:
                source_info = infer_source_transcripts(str(aligned_sam_file))
                source_transcripts = source_info['source_transcripts']
                source_gene = source_info['source_gene']
            probe_classifications = classify_probes_from_sam(
                str(aligned_sam_file),
                is_microbiome=is_microbiome,
                source_transcripts=source_transcripts,
                source_gene=source_gene,
                genome_self_match=genome_self_match,
            )
        
        # Load safe probes from safe_probes_scores.txt if it exists
        safe_probe_ids = set()
        safe_probes_file = self.output_dir / "safe_probes_scores.txt"
        if safe_probes_file.exists():
            with open(safe_probes_file, 'r', encoding='utf-8') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#') or line.startswith('=') or line.startswith('-') or line.startswith('Probe'):
                        continue
                    parts = line.split()
                    if len(parts) >= 1:
                        # Extract probe_id (e.g., "probe_5")
                        probe_id = parts[0]
                        probe_id_match = re.search(r'probe_\d+', probe_id)
                        if probe_id_match:
                            safe_probe_ids.add(probe_id_match.group(0))
        
        # Parse probe regions from FASTA
        probe_regions = []
    
        with open(candidate_probes_file, 'r', encoding='utf-8') as f:
            for line in f:
                if line.startswith('>'):
                    header = line.strip()[1:]  
                    probe_id_match = re.search(r'probe_\d+', header)
                    start_match = re.search(r'start=(\d+)', header)
                    end_match = re.search(r'end=(\d+)', header)
                    
                    if probe_id_match and start_match and end_match:
                        probe_id_display = probe_id_match.group(0)  
                        full_probe_id = header  
                        start = int(start_match.group(1))
                        end = int(end_match.group(1))
                        
                        if full_probe_id in probe_classifications:
                            classification = probe_classifications[full_probe_id]
                            status = classification['status']
                            mismatches = classification['min_mismatches']
                            # Self-aligned safe probes must also pass k-mer check
                            if status == 'safe' and probe_id_display not in safe_probe_ids:
                                status = 'medium_risk'
                                mismatches = None
                        else:
                            # No alignments. In Host Probe Design that means the
                            # probe doesn't bind the source gene at all — surface
                            # as its own `no_alignment` category so the UI can
                            # distinguish "no hit anywhere" from "k-mer failed".
                            if probe_id_display in safe_probe_ids:
                                status = 'safe'
                            elif host_internal_mode:
                                status = 'no_alignment'
                            else:
                                status = 'medium_risk'
                            mismatches = None
                        
                        probe_regions.append({
                            "probe_id": probe_id_display,
                            "start": start,
                            "end": end,
                            "status": status,
                            "mismatches": mismatches
                        })
        return probe_regions

    def generate_probes_gff3(
        self, 
        probe_regions: List[Dict], 
        reference_id: str = "reference"
    ):
        
        gff3_path = self.output_dir / "probes.gff3"
        
        with open(gff3_path, 'w', encoding='utf-8') as f:
            f.write("##gff-version 3\n")
            
            if probe_regions:
                max_pos = max(r['end'] for r in probe_regions)
                f.write(f"##sequence-region {reference_id} 1 {max_pos}\n")
            sorted_probes = sorted(probe_regions, key=lambda x: x['start'])
            
            for probe in sorted_probes:
                seqid = reference_id
                source = "probe_design"
                feature_type = "probe"
                start = probe['start'] + 1
                end = probe['end']
                score = "."
                strand = "+"
                phase = "."
                attributes = [
                    f"ID={probe['probe_id']}",
                    f"Name={probe['probe_id']}",
                    f"risk_level={probe['status']}"
                ]
                if probe['mismatches'] is not None:
                    attributes.append(f"mismatches={probe['mismatches']}")
                color_map = {
                    'high_risk': '#dc2626',
                    'medium_risk': '#ea580c',
                    'safe': '#10b981',
                    'no_alignment': '#6b7280',
                }
                if probe['status'] in color_map:
                    attributes.append(f"color={color_map[probe['status']]}")
                desc_map = {
                    'high_risk': f"High risk: {probe['mismatches']} mismatch(es)" if probe['mismatches'] is not None else "High risk: Off-target alignment",
                    'medium_risk': f"Medium risk: {probe['mismatches']} mismatches" if probe['mismatches'] is not None else "Medium risk: Failed k-mer safety analysis",
                    'safe': "Safe: No off-target alignments",
                    'no_alignment': "No alignment to source gene",
                }
                if probe['status'] in desc_map:
                    attributes.append(f"description={desc_map[probe['status']]}")
                
                attributes_str = ";".join(attributes)
                gff_line = f"{seqid}\t{source}\t{feature_type}\t{start}\t{end}\t{score}\t{strand}\t{phase}\t{attributes_str}\n"
                f.write(gff_line)
        
        logger.info(f"Generated: {gff3_path} with {len(probe_regions)} probes")
        return gff3_path
    def convert_sam_to_indexed_bam(self, sam_path: Path, skip_if_exists: bool = True):
        """
        Convert SAM to sorted, indexed BAM with improved error handling.
        Automatically fixes missing headers using reference.fasta if available.
        """
        import subprocess

        try:
            bam_path = sam_path.with_suffix('.bam')
            sorted_bam_path = sam_path.with_name(sam_path.stem + '.sorted.bam')
            if skip_if_exists and bam_path.exists():
                bai_path = Path(str(bam_path) + '.bai')
                if bai_path.exists():
                    logger.info(f"BAM already exists: {bam_path}")
                    return bam_path
            try:
                subprocess.run(['samtools', '--version'], check=True, capture_output=True)
            except FileNotFoundError:
                logger.error("samtools not found in PATH. Please install samtools.")
                raise FileNotFoundError("samtools is required but not installed")
            if not sam_path.exists():
                logger.error(f"SAM file does not exist: {sam_path}")
                raise FileNotFoundError(f"SAM file not found: {sam_path}")
            
            if sam_path.stat().st_size == 0:
                logger.warning(f"SAM file is empty: {sam_path}")
                return None
            
            has_header = False
            with open(sam_path, 'r', encoding='utf-8') as f:
                first_line = f.readline()
                has_header = first_line.startswith('@')
            
            if not has_header:
                logger.warning("SAM file missing header, attempting to rebuild")
                reference_fasta = self.output_dir / "reference.fasta"
                fixed_sam = sam_path.with_name(sam_path.stem + "_fixed.sam")

                if not reference_fasta.exists():
                    logger.error(f"Cannot rebuild SAM header: {reference_fasta} not found")
                    raise FileNotFoundError(f"Missing reference FASTA: {reference_fasta}")
                fai_path = reference_fasta.with_suffix('.fasta.fai')
                if not fai_path.exists():
                    fai_path = Path(str(reference_fasta) + '.fai')
                
                if not fai_path.exists():
                    logger.info(f"Creating FASTA index: {fai_path}")
                    try:
                        subprocess.run(
                            ['samtools', 'faidx', str(reference_fasta)], 
                            check=True,
                            capture_output=True,
                            text=True
                        )
                        fai_path = Path(str(reference_fasta) + '.fai')
                    except subprocess.CalledProcessError as e:
                        logger.error(f"Failed to create FASTA index: {e.stderr}")
                        raise

                header_path = sam_path.with_name("header.sam")
                logger.info(f"Building SAM header from {fai_path}")
                
                try:
                    with open(header_path, 'w', encoding='utf-8') as hdr, \
                         open(fai_path, 'r', encoding='utf-8') as fai:
                        hdr.write("@HD\tVN:1.6\tSO:unsorted\n")
                        for line in fai:
                            cols = line.strip().split('\t')
                            if len(cols) >= 2:
                                hdr.write(f"@SQ\tSN:{cols[0]}\tLN:{cols[1]}\n")
                except Exception as e:
                    logger.error(f"Failed to build SAM header: {e}")
                    raise

                logger.info("Merging header with SAM content")
                try:
                    with open(fixed_sam, 'w', encoding='utf-8') as out, \
                         open(header_path, 'r', encoding='utf-8') as hdr, \
                         open(sam_path, 'r', encoding='utf-8') as body:
                        out.writelines(hdr.readlines())
                        for line in body:
                            if not line.startswith('@'):
                                out.write(line)
                except Exception as e:
                    logger.error(f"Failed to merge SAM files: {e}")
                    raise

                sam_path = fixed_sam
                logger.info(f"Using fixed SAM file: {fixed_sam}")

            
            logger.info(f"Converting {sam_path.name} to BAM format")
            try:
                result = subprocess.run(
                    ['samtools', 'view', '-b', str(sam_path), '-o', str(bam_path)],
                    check=True,
                    capture_output=True,
                    text=True
                )
                if result.stderr:
                    logger.warning(f"samtools view stderr: {result.stderr}")
            except subprocess.CalledProcessError as e:
                logger.error(f"samtools view failed: {e.stderr.strip() if e.stderr else e}")
                raise

            if not bam_path.exists() or bam_path.stat().st_size == 0:
                logger.error(f"BAM file was not created or is empty: {bam_path}")
                raise ValueError("BAM conversion failed")

           
            logger.info("Sorting BAM file")
            try:
                result = subprocess.run(
                    ['samtools', 'sort', str(bam_path), '-o', str(sorted_bam_path)],
                    check=True,
                    capture_output=True,
                    text=True
                )
                if result.stderr:
                    logger.warning(f"samtools sort stderr: {result.stderr}")
            except subprocess.CalledProcessError as e:
                logger.error(f"samtools sort failed: {e.stderr.strip() if e.stderr else e}")
                raise

            
            bam_path.unlink(missing_ok=True)
            sorted_bam_path.rename(bam_path)

            if not bam_path.exists() or bam_path.stat().st_size == 0:
                logger.error(f"Sorted BAM file is missing or empty: {bam_path}")
                raise ValueError("BAM sorting failed")

            logger.info("Indexing BAM file")
            try:
                result = subprocess.run(
                    ['samtools', 'index', str(bam_path)],
                    check=True,
                    capture_output=True,
                    text=True
                )
                if result.stderr:
                    logger.warning(f"samtools index stderr: {result.stderr}")
            except subprocess.CalledProcessError as e:
                logger.error(f"samtools index failed: {e.stderr.strip() if e.stderr else e}")
                raise

            bai_path = Path(str(bam_path) + '.bai')
            if not bai_path.exists():
                logger.error(f"BAM index was not created: {bai_path}")
                raise ValueError("BAM indexing failed")

            logger.info(f"Successfully generated: {bam_path}")
            logger.info(f"Successfully generated: {bai_path}")
            if sam_path.name == "filtered_probe_alignments.sam":
                output_bam = sam_path.parent / "probe_alignments.bam"
                if not output_bam.exists():
                    import shutil
                    shutil.copy(bam_path, output_bam)
                    shutil.copy(str(bam_path) + '.bai', str(output_bam) + '.bai')
                    logger.info(f"Created copy: {output_bam}")
            return bam_path

        except subprocess.CalledProcessError as e:
            logger.error(f"samtools command failed: {e.stderr.strip() if e.stderr else str(e)}")
            raise 
        except FileNotFoundError as e:
            logger.error(f"File not found: {e}")
            raise  
        except Exception as e:
            logger.error(f"Failed to convert SAM to BAM: {e}")
            logger.exception("Detailed error traceback:")
            raise  


def generate_jbrowse_files(
    job_id: str,
    output_dir: str,
    gene_sequences_dir: str,
    probe_sequences_dir: str,
    input_type: str
) -> Dict[str, str]:
    
    try:
        generator = JBrowseFileGenerator(output_dir, job_id)
        
        
        logger.info(f"Extracting input sequence for job {job_id}")
        sequence_id, sequence = generator.extract_input_sequence(
            gene_sequences_dir,
            probe_sequences_dir,
            input_type
        )
        
        
        logger.info(f"Generating reference FASTA (length: {len(sequence)} bp)")
        reference_fasta = generator.generate_reference_fasta(
            sequence=sequence,
            sequence_id=sequence_id
        )
        
        candidate_probes_file = Path(output_dir) / "candidate_probes.fa"
        # aligned_sam_file = Path(output_dir) / "filtered_probe_alignments_annotated.sam"
        classification_sam_file = Path(output_dir) / "filtered_probe_alignments_annotated.sam"
        if not classification_sam_file.exists():
            classification_sam_file = Path(output_dir) / "filtered_probe_alignments.sam"
        # bam_conversion_sam_file = Path(output_dir) / "filtered_probe_alignments.sam"

        logger.info("Parsing probe regions and risk classifications")
        probe_regions = generator.parse_probe_regions_from_fasta(
            candidate_probes_file,
            classification_sam_file
        )
        
        logger.info(f"Generating probes GFF3 with {len(probe_regions)} probes")
        probes_gff3 = generator.generate_probes_gff3(
            probe_regions=probe_regions,
            reference_id=sequence_id
        )
        
        # bam_file = None

        # if bam_conversion_sam_file.exists():
        #     logger.info("Converting filtered_probe_alignments.sam to BAM")
        #     try:
        #         bam_file = generator.convert_sam_to_indexed_bam(bam_conversion_sam_file)
        #         logger.info(f"Successfully created filtered BAM: {bam_file}")
        #     except Exception as e:
        #         logger.warning(f"Failed to convert filtered_probe_alignments.sam: {e}")

        # main_sam_file = Path(output_dir) / "probe_alignments.sam"
        # if main_sam_file.exists():
        #     logger.info("Converting probe_alignments.sam to BAM")
        #     try:
        #         main_bam_file = generator.convert_sam_to_indexed_bam(main_sam_file)
        #         logger.info(f"Successfully created main BAM: {main_bam_file}")
        #         # Use this as the primary BAM if filtered one failed
        #         if bam_file is None:
        #             bam_file = main_bam_file
        #     except Exception as e:
        #         logger.warning(f"Failed to convert probe_alignments.sam: {e}")
        
        result = {
            'reference_fasta': str(reference_fasta),
            'reference_fasta_index': str(reference_fasta) + '.fai',
            'probes_gff3': str(probes_gff3),
            #'bam_file': str(bam_file) if bam_file else None,
            #'bam_index': str(bam_file) + '.bai' if bam_file else None,
            'sequence_id': sequence_id,
            'sequence_length': len(sequence),
            'total_probes': len(probe_regions)
        }
        
        logger.info(f"JBrowse files generated successfully for job {job_id}")
        return result
        
    except Exception as e:
        logger.exception(f"Failed to generate JBrowse files for job {job_id}")
        raise

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    
   
    result = generate_jbrowse_files(
        job_id="4f7850df-0503-48f8-9da7-6b2ca474695d",
        output_dir="../output/alignments/4f7850df-0503-48f8-9da7-6b2ca474695d",
        gene_sequences_dir="../gene_sequences",
        probe_sequences_dir="../probe_sequences",
        input_type="gene_sequence"
    )
    
    print("\nGenerated files:")
    for key, value in result.items():
        print(f"  {key}: {value}")