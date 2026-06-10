#!/usr/bin/env python3
"""Celery worker for running probe pipeline tasks."""

import os
import sys
import subprocess
import logging
from datetime import datetime
from pathlib import Path
from celery import Celery
from .jbrowse_generator import generate_jbrowse_files, JBrowseFileGenerator

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

CELERY_BROKER_URL = os.getenv("CELERY_BROKER_URL", "redis://localhost:6379/1")
CELERY_RESULT_BACKEND = os.getenv("CELERY_RESULT_BACKEND", "redis://localhost:6379/1")

celery_app = Celery(
    "probe_pipeline",
    broker=CELERY_BROKER_URL,
    backend=CELERY_RESULT_BACKEND
)

celery_app.conf.update(
    task_serializer='json',
    accept_content=['json'],
    result_serializer='json',
    timezone='UTC',
    enable_utc=True,
    task_track_started=True,
    task_time_limit=3600,
    worker_prefetch_multiplier=1,
)

ROOT = Path(__file__).resolve().parent.parent


@celery_app.task(bind=True, name="run_pipeline_task")
def run_pipeline_task(self, species, probe_length, max_mismatches, job_id, input_type, input_file, storage_dir, kmer_length=18, microbiomes="", align_microbiome=False, align_host=False, tm_range="42-47", gc_range="40-80", host_internal_mode=False, microbe_mode=False):
    """
    Run the probe design pipeline.
    
    Args:
        species: Target species (human/mouse)
        probe_length: Length of probes
        max_mismatches: Maximum allowed mismatches
        job_id: Unique job identifier
        input_type: Type of input (gene_sequence or probe_sequence)
        input_file: Filename of the saved input
        storage_dir: Directory where input file is stored
        align_microbiome: Whether to align against microbiome genomes (for microbiome species)
        align_host: Whether to align against host transcriptome (for microbiome species)
    """
    logger.info("Starting pipeline task %s", job_id)
    logger.info("Input type: %s, File: %s, Storage: %s", input_type, input_file, storage_dir)
    
    start_time = datetime.utcnow()
    
    try:
        self.update_state(
            state='PROGRESS',
            meta={'status': 'Starting pipeline...'}
        )
        
        input_path = Path(storage_dir) / input_file
        output_dir = ROOT / "output" / "alignments" / job_id
        output_dir.mkdir(parents=True, exist_ok=True)
        
        log_file = output_dir / "pipeline_log.txt"
        
        if not input_path.exists():
            raise FileNotFoundError(f"Input file not found: {input_path}")
        
        logger.info("Reading sequence from %s", input_path)
        with open(input_path, 'r', encoding='utf-8') as f:
            sequence_content = f.read()
        
        logger.info("Sequence content length: %d characters", len(sequence_content))
        
        pipeline_script = ROOT / "probe_designer.py"
        
        cmd = [
            sys.executable,
            str(pipeline_script),
            "--species", species,
            "--kmer-length", str(kmer_length),
            "--max-mismatches", str(max_mismatches),
            "--task-id", job_id,
        ]
        cmd.extend(["--tm-range", tm_range])
        cmd.extend(["--gc-range", gc_range])
        if microbiomes:
            cmd.extend(["--microbiomes", microbiomes])
        if align_microbiome:
            cmd.append("--align-microbiome")
        if align_host:
            cmd.append("--align-host")
        if host_internal_mode:
            cmd.append("--host-internal-mode")
        if microbe_mode:
            cmd.append("--microbe-mode")

        stdin_data = None
        if input_type == "gene_sequence":
            cmd.extend(["--probe-length", str(probe_length)])
            cmd.extend(["--gene-sequence-stdin"])
            stdin_data = sequence_content
        elif input_type == "probe_sequence":
            cmd.extend(["--probe-sequence-stdin"])
            stdin_data = sequence_content
        else:
            raise ValueError(f"Unknown input type: {input_type}")

        logger.info("Running command: %s (sequence via stdin)", ' '.join(cmd))

        with open(log_file, 'w', encoding='utf-8') as log:
            process = subprocess.Popen(
                cmd,
                stdin=subprocess.PIPE,
                stdout=log,
                stderr=subprocess.PIPE,
                text=True,
                cwd=str(ROOT)
            )
            
            _, stderr = process.communicate(input=stdin_data)
            exit_code = process.returncode

        
        if exit_code != 0:
            error_msg = f"Pipeline failed with exit code {exit_code}"
            if stderr:
                stderr_lines = stderr.strip().split('\n')[-10:]
                error_msg += "\n\nLast stderr lines:\n" + '\n'.join(stderr_lines)
            
            logger.error(error_msg)
            raise RuntimeError(error_msg)
        logger.info("Pipeline completed successfully for job %s", job_id)

        try:
            logger.info("Generating JBrowse files for job %s", job_id)
            
            jbrowse_result = generate_jbrowse_files(
                job_id=job_id,
                output_dir=str(output_dir),
                gene_sequences_dir=str(ROOT / "gene_sequences"),
                probe_sequences_dir=str(ROOT / "probe_sequences"),
                input_type=input_type
            )
            
            logger.info("JBrowse files generated successfully for job %s", job_id)
        except Exception as e:
            logger.error("Failed to generate JBrowse files for job %s: %s", job_id, str(e))
            logger.exception("JBrowse generation error details:")

        end_time = datetime.utcnow()
        execution_time = (end_time - start_time).total_seconds()
        logger.info("Job %s completed successfully in %.2fs", job_id, execution_time)

        

        return {
            'status': 'completed',
            'job_id': job_id,
            'execution_time': execution_time,
            'completed_at': end_time.isoformat() + 'Z',
            'output_dir': str(output_dir)
        }
        
    except (subprocess.SubprocessError, FileNotFoundError, PermissionError, OSError) as e:
        logger.exception("Pipeline task %s failed", job_id)

        end_time = datetime.utcnow()
        execution_time = (end_time - start_time).total_seconds()
        
        return {
            'status': 'failed',
            'job_id': job_id,
            'error': str(e),
            'execution_time': execution_time,
            'completed_at': end_time.isoformat() + 'Z'
        }


if __name__ == "__main__":
    celery_app.worker_main([
        'worker',
        '--loglevel=info',
        '--concurrency=2',
    ])