#!/usr/bin/env python3
"""Production FastAPI backend for probe design pipeline.

This module exposes job management endpoints with proper security,
validation, and error handling for production deployment.

Supports two input scenarios:
1. Gene sequence (pasted text)
2. Probe sequence (pasted text)
"""

import os
import shutil
import logging
from datetime import datetime
from pathlib import Path
import uuid
from typing import Optional, List
import re
from fastapi import FastAPI, HTTPException, Form, Request, Query
from fastapi.responses import FileResponse, JSONResponse
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel, Field

from .sam_parser import parse_sam_file, count_sam_alignments
from .celery_worker import celery_app, run_pipeline_task

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

app = FastAPI(
    title="Probe Pipeline API",
    version="1.0.0",
    description="API for probe design"
)


ALLOWED_ORIGINS = os.getenv("ALLOWED_ORIGINS", "http://localhost:4000").split(",")
app.add_middleware(
    CORSMiddleware,
    allow_origins=ALLOWED_ORIGINS,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


ROOT = Path(__file__).resolve().parent.parent
GENE_SEQUENCES_DIR = ROOT / "gene_sequences"
PROBE_SEQUENCES_DIR = ROOT / "probe_sequences"


for directory in [GENE_SEQUENCES_DIR, PROBE_SEQUENCES_DIR]:
    directory.mkdir(parents=True, exist_ok=True)


jobs = {}



class JobResponse(BaseModel):
    """Response model for job status."""
    job_id: str
    status: str
    submitted_at: str
    input_type: str
    info: Optional[dict] = None


class FileInfo(BaseModel):
    """Information about a job output file."""
    filename: str
    size: int


class FilesResponse(BaseModel):
    """Response model for listing job files."""
    job_id: str
    files: list[FileInfo]


class MessageResponse(BaseModel):
    """Generic message response."""
    message: str

class AlignmentRecord(BaseModel):
    """Parsed alignment record from SAM file."""
    probe_id: str
    sequence: str
    target_transcript: str
    gene_id: Optional[str] = None
    mismatches: int
    position: int
    strand: str
    species: Optional[str] = None
class AlignmentsResponse(BaseModel):
    """Response model for parsed alignments with pagination."""
    job_id: str
    total_alignments: int
    alignments: List[AlignmentRecord]
    page: int = 1
    page_size: int = 25
    total_pages: int = 0

class ProbeRegion(BaseModel):
    """A probe region with alignment status."""
    start: int
    end: int
    probe_id: str
    status: str  
    mismatches: Optional[int] = None

class RegionsResponse(BaseModel):
    """Response model for probe regions."""
    job_id: str
    total_regions: int
    gene_length: int
    regions: List[ProbeRegion]

@app.get("/jobs/{job_id}/probe-regions", response_model=RegionsResponse)
def get_probe_regions(job_id: str):
    """
    Get all probe regions with their alignment status.
    - Green (safe): Non-aligned probes
    - Red (high_risk): 0-1 mismatch
    - Orange (medium_risk): 2 mismatches
    """
    
 
    try:
        uuid.UUID(job_id)
    except ValueError as e:
        raise HTTPException(status_code=400, detail="Invalid job ID format") from e

    job_dir = ROOT / "outputs" / "alignments" / job_id
    candidate_probes_file = job_dir / "candidate_probes.fa"
    aligned_sam_file = job_dir / "filtered_probe_alignments.sam"
    
    if not job_dir.exists():
        raise HTTPException(status_code=404, detail="Job results not found")
    
    if not candidate_probes_file.exists():
        raise HTTPException(status_code=404, detail="Candidate probes file not found")
    
    try:
        regions = []
        gene_length = 0
        aligned_probes = {}  
        
        if aligned_sam_file.exists():
            nm_pattern = re.compile(r'NM:i:(\d+)')
            with open(aligned_sam_file, 'r', encoding='utf-8') as f:
                for line in f:
                    if line.startswith('@'):
                        continue
                    
                    parts = line.strip().split('\t')
                    if len(parts) < 11:
                        continue
                    
                    probe_id_full = parts[0]
                    probe_id_match = re.search(r'probe_\d+', probe_id_full)
                    if not probe_id_match:
                        continue
                    
                    probe_id = probe_id_match.group(0)
                    
                    nm_match = nm_pattern.search(line)
                    if nm_match:
                        nm = int(nm_match.group(1))
                       
                        if probe_id not in aligned_probes or nm < aligned_probes[probe_id]:
                            aligned_probes[probe_id] = nm
        
        with open(candidate_probes_file, 'r', encoding='utf-8') as f:
            for line in f:
                if line.startswith('>'):
                    header = line.strip()[1:]
                    probe_id_match = re.search(r'probe_\d+', header)
                    start_match = re.search(r'start=(\d+)', header)
                    end_match = re.search(r'end=(\d+)', header)
                    
                    if probe_id_match and start_match and end_match:
                        probe_id = probe_id_match.group(0)
                        start = int(start_match.group(1))
                        end = int(end_match.group(1))
                        
                        if probe_id in aligned_probes:
                            mismatches = aligned_probes[probe_id]
                            if mismatches <= 1:
                                status = 'high_risk' 
                            else: 
                                status = 'medium_risk'  
                        else:
                            mismatches = None
                            status = 'safe'  
                        
                        regions.append({
                            "probe_id": probe_id,
                            "start": start,
                            "end": end,
                            "status": status,
                            "mismatches": mismatches
                        })
                        
                        gene_length = max(gene_length, end)
        
        return {
            "job_id": job_id,
            "total_regions": len(regions),
            "gene_length": gene_length,
            "regions": regions
        }
        
    except Exception as e:
        logger.exception("Failed to parse probe regions for job %s", job_id)
        raise HTTPException(status_code=500, detail="Failed to parse probe regions") from e

def validate_fasta_content(content: bytes, content_type: str) -> tuple[bool, Optional[str]]:
    """Validate FASTA format and DNA alphabet.
    
    Returns: (is_valid, error_message)
    """
    try:
        text = content.decode('utf-8')
    except UnicodeDecodeError:
        return False, "File is not valid UTF-8 text"
    
    lines = text.strip().split('\n')
    if not lines:
        return False, "File is empty"

    if not lines[0].startswith('>'):
        return False, "Invalid FASTA format: must start with header line (>)"
    
    has_sequence = False
    for i, line in enumerate(lines, 1):
        line = line.strip()
        if not line:
            continue
        
        if line.startswith('>'):
            if len(line) == 1:
                return False, f"Empty header at line {i}"
            continue
        
        if not re.match(r'^[ATCGNatcgn]+$', line):
            return False, f"Invalid DNA characters at line {i}. Only ATCGN allowed."
        has_sequence = True
    
    if not has_sequence:
        return False, "No sequence data found in FASTA file"
    
    return True, None


@app.post("/jobs", response_model=JobResponse, summary="Create probe design job")
async def create_job(
    species: str = Form("human"),
    probe_length: int = Form(30),
    max_mismatches: int = Form(2),
    gene_sequence: str = Form(""),
    probe_sequence: str = Form(""),
):
    
    input_type = 'gene' if gene_sequence.strip() else ('probe' if probe_sequence.strip() else 'none')
    logger.info("Job input type: %s", input_type)
    
  
    provided = []
    
  
    if gene_sequence.strip():
        provided.append('gene_sequence')
        logger.info("Detected gene_sequence input")
    if probe_sequence.strip():
        provided.append('probe_sequence')
        logger.info("Detected probe_sequence input")

    logger.info("Provided inputs: %s", provided)

    if len(provided) == 0:
        raise HTTPException(
            status_code=400,
            detail="No input provided. Must provide one of: gene_sequence or probe_sequence"
        )
    
    if len(provided) > 1:
        raise HTTPException(
            status_code=400,
            detail=f"Multiple inputs provided: {', '.join(provided)}. Only one input type allowed per job."
        )
    
    input_type = provided[0]
    logger.info("Creating job with input_type: %s", input_type)
    content = None
    filename = None
    storage_dir = None
    
    if input_type == 'gene_sequence':
        content = gene_sequence.encode('utf-8')
        storage_dir = GENE_SEQUENCES_DIR
        storage_dir.mkdir(parents=True, exist_ok=True)
        filename = 'pasted_gene.fasta'
        
    elif input_type == 'probe_sequence':
        content = probe_sequence.encode('utf-8')
        storage_dir = PROBE_SEQUENCES_DIR
        storage_dir.mkdir(parents=True, exist_ok=True)
        filename = 'pasted_probes.fasta'
    
    is_valid, validation_error = validate_fasta_content(content, input_type)
    if not is_valid:
        raise HTTPException(status_code=400, detail=f"FASTA validation failed: {validation_error}")
    
    job_id = str(uuid.uuid4())
    safe_filename = Path(filename).name
    dest_name = f"{job_id}_{safe_filename}"
    dest = storage_dir / dest_name
    
    try:
        with open(dest, "wb") as out:
            out.write(content)
        logger.info("Saved %s file to %s", input_type, dest)
    except (PermissionError, OSError) as e:
        logger.exception("Failed to save file")
        raise HTTPException(status_code=500, detail=f"Failed to save file: {e}") from e
    
    pipeline_args = {
        'species': species,
        'probe_length': probe_length,
        'max_mismatches': max_mismatches,
        'job_id': job_id,
        'input_type': input_type,
        'input_file': dest_name,
        'storage_dir': str(storage_dir),
    }
    
    try:
        task = run_pipeline_task.apply_async(
            kwargs=pipeline_args,
            task_id=job_id
        )
        
        submitted_at = datetime.utcnow().isoformat() + "Z"
        
        jobs[job_id] = {
            "submitted_at": submitted_at,
            "input_type": input_type,
            "input_file": dest_name,
            "storage_dir": str(storage_dir),
            "species": species,
            "probe_length": probe_length,
            "max_mismatches": max_mismatches,
            "status": "PENDING",
        }
        
        logger.info("Submitted job %s with input_type=%s", job_id, input_type)
        
        return {
            "job_id": job_id,
            "status": "PENDING",
            "submitted_at": submitted_at,
            "input_type": input_type
        }
        
    except Exception as e:
        logger.exception("Failed to submit Celery task")
        if dest.exists():
            dest.unlink()
        raise HTTPException(status_code=500, detail=f"Failed to create job: {e}") from e


@app.get("/jobs/{job_id}", response_model=JobResponse)
def get_job(job_id: str):
    """Get job status and information."""
    try:
        uuid.UUID(job_id)
    except ValueError as e:
        raise HTTPException(status_code=400, detail=f"Invalid job ID format: {e}") from e
    try:
        async_result = celery_app.AsyncResult(job_id)
        state = async_result.state
        raw_info = async_result.info
        
        if isinstance(raw_info, dict):
            info = raw_info
        elif isinstance(raw_info, Exception):
            info = {
                "status": "failed",
                "error": str(raw_info),
                "error_type": type(raw_info).__name__
            }
            state = 'FAILURE'  
        elif raw_info is not None:
            info = {
                "status": "unknown",
                "info": str(raw_info)
            }
        else:
            info = {}
        if state in ['SUCCESS', 'FAILURE']:
            job_dir = ROOT / "outputs" / "alignments" / job_id
            log_file = job_dir / "pipeline_log.txt"
            
            if log_file.exists():
                stats = parse_pipeline_stats(log_file)
                info["stats"] = stats
        
        submitted_at = "unknown"
        input_type = "unknown"
        
        if job_id in jobs:
            jobs[job_id]["status"] = state
            submitted_at = jobs[job_id].get("submitted_at", "unknown")
            input_type = jobs[job_id].get("input_type", "unknown")
        else:
            job_dir = ROOT / "outputs" / "alignments" / job_id
            if job_dir.exists():
                try:
                    submitted_at = datetime.fromtimestamp(job_dir.stat().st_ctime).isoformat() + "Z"
                except (KeyError, ValueError, TypeError):
                    submitted_at = "unknown"
                    logger.exception("Failed to get submission time for job %s", job_id)

                if (job_dir / f"{job_id}_pasted_gene.fasta").exists() or any(f.startswith(job_id) and 'gene' in f for f in os.listdir(GENE_SEQUENCES_DIR) if os.path.isfile(GENE_SEQUENCES_DIR / f)):
                    input_type = "gene_sequence"
                elif (job_dir / f"{job_id}_pasted_probes.fasta").exists() or any(f.startswith(job_id) and 'probe' in f for f in os.listdir(PROBE_SEQUENCES_DIR) if os.path.isfile(PROBE_SEQUENCES_DIR / f)):
                    input_type = "probe_sequence"
            elif state == 'PENDING':
                pass
            else:
                raise HTTPException(status_code=404, detail="Job not found")
            
        info["species"] = jobs[job_id].get("species", "unknown")
        return {
            "job_id": job_id,
            "status": state,
            "submitted_at": submitted_at,
            "input_type": input_type,
            "info": info
        }
    except HTTPException:
        raise
    except (KeyError, ValueError, FileNotFoundError, PermissionError, OSError) as e:
        logger.exception("Failed to get job status for %s", job_id)
        raise HTTPException(status_code=500, detail=f"Failed to retrieve job status: {str(e)}") from e


def parse_pipeline_stats(log_path):
    """Parse pipeline_log.txt to extract key statistics."""
    stats = {
        "candidate_probes": None,
        "total_alignments": None,
        "non_aligned_probes": None,
        "filtered_alignments": None
    }
    
    try:
        with open(log_path, 'r', encoding='utf-8') as f:
            content = f.read()
        
        match = re.search(r'Found ([\d,]+) candidate probes passing GC filter', content)
        if match:
            stats["candidate_probes"] = int(match.group(1).replace(',', ''))
        
        match = re.search(r'Total alignments scanned:\s*([\d,]+)', content)
        if match:
            stats["total_alignments"] = int(match.group(1).replace(',', ''))
        
        match = re.search(r'Non-aligned probes:\s*([\d,]+)', content)
        if match:
            stats["non_aligned_probes"] = int(match.group(1).replace(',', ''))
        
        match = re.search(r'Kept in filtered \(NM<=\d+\):\s*([\d,]+)', content)
        print(f"DEBUG: Filtered alignments match: {match}")  
        if match:
            print(f"DEBUG: Captured value: {match.group(1)}")  
            stats["filtered_alignments"] = int(match.group(1).replace(',', ''))
    except (FileNotFoundError, PermissionError, OSError, ValueError):
        logger.exception("Failed to parse pipeline stats from %s", log_path)
    
    return stats


@app.get("/jobs/{job_id}/files", response_model=FilesResponse)
def list_job_files(job_id: str):
    """List output files for a completed job."""
    
    try:
        uuid.UUID(job_id)
    except ValueError as e:
        raise HTTPException(status_code=400, detail="Invalid job ID format") from e
    
    job_dir = ROOT / "outputs" / "alignments" / job_id
    
    if not job_dir.exists():
        raise HTTPException(status_code=404, detail="Job results not found")
    
    if not job_dir.is_dir():
        raise HTTPException(status_code=500, detail="Invalid job directory")
    
    files = []
    try:
        for p in job_dir.iterdir():
            if p.is_file():
                if p.name.lower() in ["transcripts.fa", "pipeline_log.txt", "job_result.json"]:
                    continue
                files.append({
                    "filename": p.name,
                    "size": p.stat().st_size
                })
        
        return {"job_id": job_id, "files": files}
        
    except Exception as e:
        logger.exception("Failed to list files for job %s", job_id)
        raise HTTPException(status_code=500, detail=f"Failed to list job files: {str(e)}") from e

@app.get("/jobs/{job_id}/alignments", response_model=AlignmentsResponse)
def get_job_alignments(
    job_id: str, 
    page: int = Query(1, ge=1, description="Page number (1-indexed)"),
    page_size: int = Query(25, ge=1, le=100, description="Items per page"),
    mismatch: Optional[int] = Query(None, ge=0, description="Filter by specific mismatch count"),
    probe_id: Optional[str] = Query(None, description="Get all alignments for specific probe")  # ADD THIS LINE
):
    """
    Get parsed alignment data from filtered_probe_alignments_annotated.sam with pagination
    
    Parameters:
    - job_id: The job UUID
    - page: Page number (default 1)
    - page_size: Items per page (default 25, max 100)
    - probe_id: Filter by specific probe ID (returns all alignments for that probe)  # ADD THIS
    """
    try:
        uuid.UUID(job_id)
    except ValueError as e:
        raise HTTPException(status_code=400, detail="Invalid job ID format") from e

    job_dir = ROOT / "outputs" / "alignments" / job_id
    annotated_sam = job_dir / "filtered_probe_alignments_annotated.sam"
    sam_file = annotated_sam if annotated_sam.exists() else job_dir / "filtered_probe_alignments.sam"

    if not job_dir.exists():
        raise HTTPException(status_code=404, detail="Job results not found")
    sam_file = job_dir / "filtered_probe_alignments_annotated.sam"

    if not sam_file.exists():
        sam_file = job_dir / "filtered_probe_alignments.sam"
        logger.warning("Annotated SAM not found for job %s, using non-annotated", job_id)
    
    if not sam_file.exists():
        logger.info("No alignment file found for job %s - all probes are safe", job_id)
        return {
            "job_id": job_id,
            "total_alignments": 0,
            "alignments": [],
            "page": page,
            "page_size": page_size,
            "total_pages": 0
        }
    
    if probe_id is not None:
        try:
            alignments = parse_sam_file(str(sam_file), probe_id_filter=probe_id)
            return {
                "job_id": job_id,
                "total_alignments": len(alignments),
                "alignments": alignments,
                "page": 1,
                "page_size": len(alignments),
                "total_pages": 1
            }
        except Exception as e:
            logger.exception("Failed to parse alignments for probe %s in job %s", probe_id, job_id)
            raise HTTPException(status_code=500, detail="Failed to parse alignment data") from e
    
    try:
        if mismatch is not None:
            all_alignments = parse_sam_file(str(sam_file), mismatch_filter=mismatch)
            total_count = len(all_alignments)
            total_pages = (total_count + page_size - 1) // page_size if total_count > 0 else 0
            
            if page > total_pages and total_pages > 0:
                raise HTTPException(status_code=400, detail=f"Page {page} exceeds total pages {total_pages}")
            
            start_idx = (page - 1) * page_size
            end_idx = start_idx + page_size
            alignments = all_alignments[start_idx:end_idx]
        else:
            total_count = count_sam_alignments(str(sam_file))
            total_pages = (total_count + page_size - 1) // page_size
            
            if page > total_pages and total_pages > 0:
                raise HTTPException(status_code=400, detail=f"Page {page} exceeds total pages {total_pages}")
            
            offset = (page - 1) * page_size
            alignments = parse_sam_file(str(sam_file), limit=page_size, offset=offset)
        
        return {
            "job_id": job_id,
            "total_alignments": total_count,
            "alignments": alignments,
            "page": page,
            "page_size": page_size,
            "total_pages": total_pages
        }
        
    except Exception as e:
        logger.exception("Failed to parse alignments for job %s", job_id)
        raise HTTPException(status_code=500, detail="Failed to parse alignment data") from e
    

@app.get("/jobs/{job_id}/download/{filename}")
def download_job_file(job_id: str, filename: str):
    """Download a specific output file from a completed job."""
    
    try:
        uuid.UUID(job_id)
    except ValueError as e:
        raise HTTPException(status_code=400, detail="Invalid job ID format") from e

    if "/" in filename or "\\" in filename or filename.startswith("."):
        raise HTTPException(status_code=400, detail="Invalid filename")
    
    job_dir = ROOT / "outputs" / "alignments" / job_id
    
    if not job_dir.exists():
        raise HTTPException(status_code=404, detail="Job results not found")
    
    file_path = job_dir / filename
    
    allowed_files = [
        "probe_alignments.sam",
        "filtered_probe_alignments.sam",
        "filtered_probe_alignments_annotated.sam", 
        "non_aligned_probes_scores.txt",
        "filtered_probe_alignments.bam",
        "filtered_probe_alignments.bam.bai",
        "probe_alignments.bam",
        "probe_alignments.bam.bai",
        "non_aligned_probes.fa",
        "non_aligned_probes_scores.txt",
        "candidate_probes.fa",
        "reference.fasta",
        "reference.fasta.fai",
        "probes.gff3"
    ]
    allowed_normalized = {f.lower().strip() for f in allowed_files}
    if filename.lower() not in allowed_normalized:
        raise HTTPException(status_code=403, detail="File not available for download")
    
    if not file_path.exists():
        raise HTTPException(status_code=404, detail="File not found")
    
    try:
        file_path.resolve().relative_to(job_dir.resolve())
    except ValueError as e:
        logger.warning("Path traversal attempt: %s for job %s", filename, job_id)
        raise HTTPException(status_code=403, detail="Access denied") from e

    return FileResponse(
        path=str(file_path),
        filename=filename,
        media_type="application/octet-stream"
    )
@app.delete("/jobs/{job_id}", response_model=MessageResponse)
def cancel_job(job_id: str):
    """Cancel a running job and delete its results."""
    
    try:
        uuid.UUID(job_id)
    except ValueError as e:
        raise HTTPException(status_code=400, detail="Invalid job ID format") from e

    if job_id not in jobs:
        raise HTTPException(status_code=404, detail="Job not found")
    


@app.exception_handler(Exception)
async def general_exception_handler(request: Request, exc: Exception):
    """Global exception handler for unhandled errors."""
    logger.exception("Unhandled exception on %s %s", request.method, request.url.path)
    return JSONResponse(
        status_code=500,
        content={
            "detail": "Internal server error",
            "path": str(request.url.path)
        }
    )


if __name__ == "__main__":
    import uvicorn
    port = int(os.getenv("API_PORT", "8002"))
    host = os.getenv("API_HOST", "0.0.0.0")
    uvicorn.run(
        "backend.main_app:app",
        host=host,
        port=port,
        reload=False,  
        access_log=True
    )