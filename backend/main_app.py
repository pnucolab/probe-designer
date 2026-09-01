#!/usr/bin/env python3
"""Production FastAPI backend for probe design pipeline.

This module exposes job management endpoints with proper security,
validation, and error handling for production deployment.

Supports two input scenarios:
1. Gene sequence (pasted text)
2. Probe sequence (pasted text)
"""

import os
import json
import logging
from datetime import datetime
from pathlib import Path
import uuid
from typing import Optional, List, Dict
import re
from fastapi import FastAPI, HTTPException, Form, Request, Query, Response
from fastapi.responses import FileResponse, JSONResponse, StreamingResponse
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel

from .sam_parser import parse_sam_file, count_sam_alignments, same_genome
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
    target_sequence: Optional[str] = None
    target_transcript: str
    gene_id: Optional[str] = None
    mismatches: int
    substitutions: Optional[int] = None
    bulges: Optional[int] = None
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

@app.get("/api/organisms", summary="Organism registry from config/organisms.yml")
def get_organism_registry():
    """Return the host + microbiome catalog config so the frontend builds
    its dropdowns dynamically from `config/organisms.yml` (no hard-coded
    species in the UI)."""
    import sys as _sys
    _root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    if _root not in _sys.path:
        _sys.path.insert(0, _root)
    from core.organism_registry import load_registry
    return load_registry().to_json_dict()


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

    job_dir = ROOT / "output" / "alignments" / job_id
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

_HOST_LOOKUP: Dict[str, object] = {"loaded": False, "transcript": {}, "symbol": {}}
_HOST_TOKEN_SPLIT = re.compile(r'[\s|,;:()\[\]]+')
_HOST_BINOMIAL_RE = re.compile(r'\b(Homo\s+sapiens|Mus\s+musculus)\b', re.IGNORECASE)


def _load_host_lookup() -> None:
    """Lazy-load human + mouse transcript→gene mappings into reverse-lookup dicts.

    Both maps are merged into a single transcript-id → species dict and a single
    gene-symbol → species dict. Ensembl version suffix is stripped so both
    'ENST00000372348' and 'ENST00000372348.5' resolve. Idempotent.
    """
    if _HOST_LOOKUP["loaded"]:
        return
    import sys as _sys
    _root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    if _root not in _sys.path:
        _sys.path.insert(0, _root)
    from probe_designer import samannotator  # reuses existing cache loader
    transcript_map: Dict[str, str] = {}
    symbol_map: Dict[str, str] = {}
    for sp in ("human", "mouse"):
        ann = samannotator()
        try:
            ann._ensure_mappings_loaded(sp)
        except Exception:
            logger.exception("Failed to load %s mappings for host-token check", sp)
            continue
        for tid, gname in ann.mappings.items():
            transcript_map.setdefault(tid, sp)
            transcript_map.setdefault(tid.split('.')[0], sp)
            if gname:
                symbol_map.setdefault(gname, sp)

    # Additional source: NCBI RefSeq RNA accessions (e.g. NM_005228, XM_…).
    # Built once via data/ncbi/build_refseq_lookup.py. Added on top of the
    # GENCODE/Ensembl lookup; existing entries are not overwritten.
    refseq_path = os.path.join(_root, "data", "ncbi", "refseq_to_gene.tsv")
    refseq_count = 0
    if os.path.exists(refseq_path):
        try:
            with open(refseq_path, "r", encoding="utf-8") as fh:
                next(fh, None)  # header
                for line in fh:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 3:
                        continue
                    rid, sym, sp = parts[0], parts[1], parts[2]
                    if not rid or not sp:
                        continue
                    transcript_map.setdefault(rid, sp)
                    if sym:
                        symbol_map.setdefault(sym, sp)
                    refseq_count += 1
        except OSError:
            logger.exception("Failed to load NCBI RefSeq lookup at %s", refseq_path)

    _HOST_LOOKUP["transcript"] = transcript_map
    _HOST_LOOKUP["symbol"] = symbol_map
    _HOST_LOOKUP["loaded"] = True
    logger.info("Host-token lookup loaded: %d transcript IDs, %d gene symbols (incl. %d RefSeq)",
                len(transcript_map), len(symbol_map), refseq_count)


def detect_host_tokens(fasta_text: str) -> List[Dict[str, str]]:
    """Scan FASTA headers for host (human/mouse) transcript IDs or gene symbols.

    Returns a deduped list of {token, kind, species}. Empty list means clean.
    """
    if not fasta_text:
        return []
    _load_host_lookup()
    tmap = _HOST_LOOKUP["transcript"]
    smap = _HOST_LOOKUP["symbol"]
    seen = set()
    matches: List[Dict[str, str]] = []
    for line in fasta_text.splitlines():
        line = line.strip()
        if not line.startswith('>'):
            continue
        # Latin binomial ("Homo sapiens" / "Mus musculus") in the header is a
        # definitive host signal — independent of any gene-name lookup.
        bm = _HOST_BINOMIAL_RE.search(line)
        if bm:
            binom = bm.group(1)
            sp = "human" if binom.lower().startswith("homo") else "mouse"
            key = ("species_name", binom.lower(), sp)
            if key not in seen:
                seen.add(key)
                matches.append({"token": binom, "kind": "species_name", "species": sp})
        for raw in _HOST_TOKEN_SPLIT.split(line[1:]):
            if not raw:
                continue
            # Candidates: full token, version-stripped form, underscore-split parts.
            # Trying the full form first preserves exact-match for the ~941 host
            # symbols that contain '_' (e.g. '5_8S_rRNA', 'Metazoa_SRP'); the
            # underscore-split parts catch user-supplied composites like
            # 'EGFR_transcript' or 'ABL1_var2'.
            candidates = [raw, re.sub(r'\.\d+$', '', raw)]
            if '_' in raw:
                candidates.extend(p for p in raw.split('_') if p)
            for cand in candidates:
                if cand in tmap:
                    sp = tmap[cand]
                    key = ('transcript', cand, sp)
                    if key not in seen:
                        seen.add(key)
                        matches.append({"token": cand, "kind": "transcript", "species": sp})
                    break
                if cand in smap:
                    sp = smap[cand]
                    key = ('symbol', cand, sp)
                    if key not in seen:
                        seen.add(key)
                        matches.append({"token": cand, "kind": "symbol", "species": sp})
                    break
    return matches


def _host_token_message(matches: List[Dict[str, str]]) -> str:
    if not matches:
        return ""
    # Show the most user-recognizable example: gene symbol > species name > transcript ID.
    priority = {"symbol": 0, "species_name": 1, "transcript": 2}
    best = min(matches, key=lambda m: priority.get(m["kind"], 9))
    return (
        f"Host sequence detected (e.g., {best['token']}). "
        "Please paste a microbial gene/transcript instead."
    )


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


@app.post("/validate-host-token", summary="Detect host (human/mouse) tokens in FASTA headers")
async def validate_host_token(fasta: str = Form(""), species: str = Form("")):
    """Detect host identifiers in the FASTA header and return a routing verdict.

    Response fields:
        is_host: bool — any host token detected.
        matches: list — raw match records.
        detected_species: list — unique host species inferred (e.g. ["human"]).
        should_block: bool — frontend should reject submission.
        host_internal_mode: bool — backend will run host-internal design.
        message: str — user-facing message explaining the verdict.
    """
    import sys as _sys
    _root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    if _root not in _sys.path:
        _sys.path.insert(0, _root)
    from core.organism_registry import load_registry as _load_reg
    _registry = _load_reg()

    matches = detect_host_tokens(fasta)
    if not matches:
        return {
            "is_host": False, "matches": [], "detected_species": [],
            "should_block": False, "host_internal_mode": False, "message": "",
        }

    detected = sorted({m['species'] for m in matches})
    sel = (species or "").strip()
    should_block = False
    host_internal_mode = False
    message = _host_token_message(matches)

    if _registry.is_host_species(sel):
        if len(detected) > 1:
            should_block = True
            message = "Pasted FASTA contains identifiers from multiple host species; cannot determine target."
        elif detected[0] != sel:
            should_block = True
            message = f"Pasted FASTA is {detected[0]} but selected target is {sel}. Select {detected[0]} or paste a {sel} sequence."
        else:
            host_internal_mode = True
            message = f"Host sequence detected — designing probes within the {sel} transcriptome (source-only filter)."
    elif _registry.is_microbiome_species(sel):
        message = f"Host sequence detected; running standard pipeline against {sel}."
    elif sel:
        should_block = True
    else:
        # No species supplied: keep legacy behavior so callers without context
        # still get the original blocking message.
        should_block = True

    return {
        "is_host": True,
        "matches": matches,
        "detected_species": detected,
        "should_block": should_block,
        "host_internal_mode": host_internal_mode,
        "message": message,
    }


@app.post("/jobs", response_model=JobResponse, summary="Create probe design job")
async def create_job(
    species: str = Form("human"),
    probe_length: int = Form(36),
    max_mismatches: int = Form(2),
    max_bulges: int = Form(0),
    kmer_length: int = Form(16),
    gene_sequence: str = Form(""),
    probe_sequence: str = Form(""),
    microbiomes: str = Form(""),
    align_microbiome: str = Form("false"),
    align_host: str = Form("false"),
    tm_range: str = Form("42-47"),
    gc_range: str = Form("40-80"),
    mode: str = Form("microbe"),
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
    if probe_length < 20 or probe_length > 50:
        raise HTTPException(
            status_code=400,
            detail="Probe length must be between 20 and 50 bp"
        )
    if max_bulges < 0 or max_bulges > 2:
        raise HTTPException(
            status_code=400,
            detail="Max bulges must be a whole number from 0 to 2 (inclusive)."
        )
    if max_mismatches < 0:
        raise HTTPException(
            status_code=400,
            detail="Max mismatches must be a whole number of 0 or greater."
        )
    if kmer_length < 0:
        raise HTTPException(
            status_code=400,
            detail="K-mer length must be a whole number of 0 or greater."
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
        try:
            lines = probe_sequence.strip().split('\n')
            sequences = [line.strip() for line in lines if line.strip() and not line.startswith('>')]
            if sequences:
                lengths = set(len(seq) for seq in sequences)
                if len(lengths) == 1:
                    probe_length = lengths.pop()
                    logger.info(f"Detected probe length from input: {probe_length} bp")
        except Exception:
            pass
    
    is_valid, validation_error = validate_fasta_content(content, input_type)
    if not is_valid:
        raise HTTPException(status_code=400, detail=f"FASTA validation failed: {validation_error}")

    # Defense-in-depth: handle host transcript/gene names. Frontend runs the
    # same routing on Submit; this also catches direct API consumers
    # (curl, /docs, scripted clients).
    #
    # Routing:
    #   - Host detected + host species (human/mouse) selected, species matches
    #       -> proceed with --host-internal-mode (design within source gene)
    #   - Host detected + host species selected, species mismatch
    #       -> block with explicit mismatch error
    #   - Host detected + microbiome species selected
    #       -> proceed with current pipeline unchanged
    # Pull host/microbiome membership from the organism registry instead of
    # a hard-coded set, so new hosts (e.g. zebrafish) added via YAML are
    # recognized here automatically.
    import sys as _sys
    _root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    if _root not in _sys.path:
        _sys.path.insert(0, _root)
    from core.organism_registry import load_registry as _load_reg
    _registry = _load_reg()
    is_host_species_id = _registry.is_host_species
    is_microbiome_species_id = _registry.is_microbiome_species

    host_matches = detect_host_tokens(content.decode('utf-8', errors='ignore'))
    host_internal_mode = False
    # Microbial Probe Design mode: the user explicitly chose to treat input as
    # microbial. Skip host-token routing entirely — any host alignment further
    # down the pipeline is handled as cross-reactivity, not as a redirect.
    if (mode or "").lower() == "microbe":
        host_matches = []
    if host_matches:
        detected_species = {m['species'] for m in host_matches}
        if is_host_species_id(species):
            if len(detected_species) > 1:
                raise HTTPException(
                    status_code=400,
                    detail="Pasted FASTA contains identifiers from multiple host species; cannot determine target.",
                )
            detected = next(iter(detected_species))
            if detected != species:
                raise HTTPException(
                    status_code=400,
                    detail=f"Pasted FASTA is {detected} but selected target is {species}. Select {detected} or paste a {species} sequence.",
                )
            host_internal_mode = True
        elif not is_microbiome_species_id(species):
            raise HTTPException(status_code=400, detail=_host_token_message(host_matches))

    # Honor the user's explicit Host Probe Design choice even when the header
    # carries no recognizable token (e.g. `>epidermal mRNA`). The host-internal
    # filter groups alignments by gene and keeps probes that hit only one gene
    # (all isoforms count as one), so it doesn't need a header-derived source.
    if (mode or "").lower() == "host" and is_host_species_id(species):
        host_internal_mode = True

    if input_type == 'gene_sequence':
        lines = gene_sequence.strip().split('\n')
        current_header = ''
        current_seq = ''
        short_seqs = []
        header_count = 0
        for line in lines:
            stripped = line.strip()
            if stripped.startswith('>'):
                header_count += 1
                if current_header and 0 < len(current_seq) < probe_length:
                    short_seqs.append(f'"{current_header}" ({len(current_seq)} bp)')
                current_header = stripped[1:].split()[0] if stripped[1:].strip() else 'unnamed'
                current_seq = ''
            else:
                current_seq += stripped.replace(' ', '')
        if current_header and 0 < len(current_seq) < probe_length:
            short_seqs.append(f'"{current_header}" ({len(current_seq)} bp)')
        if header_count > 1:
            raise HTTPException(
                status_code=400,
                detail=f"Only a single FASTA sequence is allowed per job. Found {header_count} sequences. Please submit one sequence at a time."
            )
        if short_seqs:
            raise HTTPException(
                status_code=400,
                detail=f"Sequence(s) shorter than probe length ({probe_length} bp): {', '.join(short_seqs)}. Remove short sequences or reduce the probe length."
            )

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
        'microbiomes': microbiomes,
        'probe_length': probe_length,
        'max_mismatches': max_mismatches,
        'max_bulges': max_bulges,
        'kmer_length': kmer_length,
        'job_id': job_id,
        'input_type': input_type,
        'input_file': dest_name,
        'storage_dir': str(storage_dir),
        'align_microbiome': align_microbiome.lower() == 'true',
        'align_host': align_host.lower() == 'true',
        'tm_range': tm_range,
        'gc_range': gc_range,
        'host_internal_mode': host_internal_mode,
        'microbe_mode': (mode or "").lower() == "microbe",
    }
    
    submitted_at = datetime.utcnow().isoformat() + "Z"
    pipeline_args['submitted_at'] = submitted_at

    try:
        task = run_pipeline_task.apply_async(
            kwargs=pipeline_args,
            task_id=job_id
        )

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
def get_job(job_id: str, response: Response):
    """Get job status and information."""
    response.headers["Cache-Control"] = "no-store"
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
        # The Celery/Redis result is ephemeral (expires, and is lost if Redis
        # restarts), which would leave a finished job stuck showing PENDING. If
        # the job completed on disk, trust that over the missing result.
        if state not in ['SUCCESS', 'FAILURE']:
            _log = ROOT / "output" / "alignments" / job_id / "pipeline_log.txt"
            if _log.exists():
                try:
                    _tail = _log.read_text(errors='ignore')[-4000:]
                except OSError:
                    _tail = ""
                if "PIPELINE COMPLETED SUCCESSFULLY" in _tail:
                    state = 'SUCCESS'
                    if not isinstance(info, dict):
                        info = {}
        # Durable timestamps written by the worker (survive Redis loss/expiry).
        job_meta = {}
        _meta_path = ROOT / "output" / "alignments" / job_id / "job_meta.json"
        if _meta_path.exists():
            try:
                with open(_meta_path, encoding="utf-8") as _mf:
                    job_meta = json.load(_mf) or {}
            except (OSError, ValueError):
                job_meta = {}

        if state in ['SUCCESS', 'FAILURE']:
            job_dir = ROOT / "output" / "alignments" / job_id
            log_file = job_dir / "pipeline_log.txt"

            if log_file.exists():
                stats = parse_pipeline_stats(log_file)
                info["stats"] = stats
            # completed_at lives only in the (ephemeral) Celery result; recover it
            # from job_meta.json, then the pipeline log's mtime, so the UI still
            # shows a completion time after a restart/expiry.
            if not info.get("completed_at"):
                info["completed_at"] = job_meta.get("completed_at") or (
                    datetime.fromtimestamp(log_file.stat().st_mtime).isoformat() + "Z"
                    if log_file.exists() else None)

        submitted_at = job_meta.get("submitted_at") or "unknown"
        input_type = "unknown"

        if job_id in jobs:
            jobs[job_id]["status"] = state
            if submitted_at == "unknown":
                submitted_at = jobs[job_id].get("submitted_at", "unknown")
            input_type = jobs[job_id].get("input_type", "unknown")
        else:
            job_dir = ROOT / "output" / "alignments" / job_id
            if job_dir.exists():
                if submitted_at == "unknown":
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
            
        info["species"] = jobs.get(job_id, {}).get("species", "unknown")
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
        "safe_probes": None,
        "filtered_alignments": None,
        "input_probes": None,
        "gc_passed": None,
        "gc_rejected": None,
        "tm_rejected": None,
        "homopolymer_rejected": None,
        "max_bulges": 0
    }

    try:
        with open(log_path, 'r', encoding='utf-8') as f:
            content = f.read()

        match = re.search(r'Max bulges:\s*(\d+)', content)
        if match:
            stats["max_bulges"] = int(match.group(1))

        match = re.search(r'Input probes:\s*([\d,]+)', content)
        if match:
            stats["input_probes"] = int(match.group(1).replace(',', ''))

        match = re.search(r'Found ([\d,]+) candidate probes passing GC filter', content)
        if match:
            stats["gc_passed"] = int(match.group(1).replace(',', ''))

        if stats["input_probes"] is not None and stats["gc_passed"] is not None:
            stats["gc_rejected"] = stats["input_probes"] - stats["gc_passed"]

        match = re.search(r'Tm filter.*rejected ([\d,]+)/', content)
        if match:
            stats["tm_rejected"] = int(match.group(1).replace(',', ''))

        match = re.search(r'Homopolymer filter.*rejected ([\d,]+)/', content)
        if match:
            stats["homopolymer_rejected"] = int(match.group(1).replace(',', ''))

        match = re.search(r'Passed both filters: ([\d,]+)/', content)
        if match:
            stats["candidate_probes"] = int(match.group(1).replace(',', ''))
        else:
            if stats["gc_passed"] is not None:
                stats["candidate_probes"] = stats["gc_passed"]
        
        match = re.search(r'Total alignments scanned:\s*([\d,]+)', content)
        if match:
            stats["total_alignments"] = int(match.group(1).replace(',', ''))
        
        # Potential on-target = self-aligned safe + non-aligned probes
        non_aligned = 0
        self_aligned_safe = 0
        match = re.search(r'Non-aligned probes:\s*([\d,]+)', content)
        if match:
            non_aligned = int(match.group(1).replace(',', ''))
        match = re.search(r'Self-aligned safe probes:\s*([\d,]+)', content)
        if match:
            self_aligned_safe = int(match.group(1).replace(',', ''))
        stats["non_aligned_probes"] = self_aligned_safe + non_aligned
        # Fall back to Initial safe probes for older logs
        if stats["non_aligned_probes"] == 0:
            match = re.search(r'Initial safe probes:\s*([\d,]+)', content)
            if match:
                stats["non_aligned_probes"] = int(match.group(1).replace(',', ''))
        
        match = re.search(r'Total safe probes:\s*([\d,]+)', content)
        if not match:
            match = re.search(r'Final safe probes:\s*([\d,]+)', content)
        if match:
            stats["safe_probes"] = int(match.group(1).replace(',', ''))
        
        match = re.search(r'Kept in filtered \(NM<=\d+\):\s*([\d,]+)', content)
        if match:
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
    
    job_dir = ROOT / "output" / "alignments" / job_id
    
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
    page_size: int = Query(25, ge=1, le=1000000, description="Items per page"),
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

    job_dir = ROOT / "output" / "alignments" / job_id
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
    
    src_gene, src_transcripts = _load_job_source(job_dir)
    genome_self_match = _load_job_genome_self_match(job_dir)

    if probe_id is not None:
        try:
            alignments = parse_sam_file(
                str(sam_file), probe_id_filter=probe_id,
                source_gene=src_gene, source_transcripts=src_transcripts,
                genome_self_match=genome_self_match,
            )
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
        # Streaming fast path: skip to the requested offset and parse only the
        # records for this page. Avoids the full-SAM parse + sort that
        # previously cost ~10s on 500K-row jobs.
        offset = (page - 1) * page_size
        alignments = parse_sam_file(
            str(sam_file), offset=offset, limit=page_size, mismatch_filter=mismatch,
            source_gene=src_gene, source_transcripts=src_transcripts,
            genome_self_match=genome_self_match,
        )
        total_count = count_sam_alignments(
            str(sam_file), source_gene=src_gene, source_transcripts=src_transcripts,
            genome_self_match=genome_self_match,
        )
        total_pages = (total_count + page_size - 1) // page_size if total_count > 0 else 0
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


def _load_job_source(job_dir) -> tuple:
    """Read source_info.json for a job. Returns (source_gene, source_transcripts).
    Both are empty when the file is missing or has no inferred source."""
    src_path = job_dir / "source_info.json"
    if not src_path.exists():
        return None, []
    try:
        import json as _json
        with open(src_path, 'r', encoding='utf-8') as fh:
            info = _json.load(fh)
        return (info.get('source_gene') or None, list(info.get('source_transcripts') or []))
    except (OSError, ValueError):
        logger.warning("Failed to read source_info.json at %s", src_path)
        return None, []


def _load_job_genome_self_match(job_dir) -> bool:
    """Whether the job used probe-id<->genome-id self matching (microbiome
    probe-input). When true, off-target views drop alignments to a probe's own
    genome (rname belongs to the probe's genome)."""
    src_path = job_dir / "source_info.json"
    if not src_path.exists():
        return False
    try:
        import json as _json
        with open(src_path, 'r', encoding='utf-8') as fh:
            return bool(_json.load(fh).get('genome_self_match'))
    except (OSError, ValueError):
        return False


def _gc_percent(seq: str) -> str:
    if not seq:
        return ''
    clean = seq.replace('-', '').upper()
    if not clean:
        return ''
    gc = sum(1 for c in clean if c in 'GC')
    return f"{(gc / len(clean)) * 100:.1f}%"


def _load_probe_risk_levels(job_dir):
    """Map probe_id -> risk label for probes that failed because of off-target
    alignments (high risk, or medium risk from 1-2 mismatches). Probes that
    failed k-mer analysis, safe probes and no-alignment probes are omitted."""
    gff3 = job_dir / "probes.gff3"
    labels = {}
    if not gff3.exists():
        return labels
    try:
        with open(gff3, 'r', encoding='utf-8') as fh:
            for line in fh:
                if not line or line.startswith('#'):
                    continue
                cols = line.rstrip('\n').split('\t')
                if len(cols) < 9:
                    continue
                attrs = {}
                for kv in cols[8].split(';'):
                    i = kv.find('=')
                    if i > 0:
                        attrs[kv[:i]] = kv[i + 1:]
                pid = attrs.get('ID')
                risk = attrs.get('risk_level')
                desc = (attrs.get('description') or '').lower()
                if not pid:
                    continue
                if risk == 'high_risk':
                    labels[pid] = 'High risk'
                elif risk == 'medium_risk' and 'k-mer' not in desc:
                    labels[pid] = 'Medium risk (mismatch)'
    except OSError:
        logger.warning("Failed to read probes.gff3 at %s", gff3)
    return labels


@app.get("/jobs/{job_id}/alignments/download")
def download_job_alignments_tsv(job_id: str, mismatch: Optional[int] = Query(None, ge=0),
                                risk_report: bool = Query(False)):
    """
    Stream all parsed alignments as a TSV. Used by the UI's 'download full
    alignments' button so the browser doesn't have to hold the full set in memory.

    When risk_report=true, restrict rows to medium-risk (mismatch) and high-risk
    probes (k-mer / safe / no-alignment probes excluded).
    """
    try:
        uuid.UUID(job_id)
    except ValueError as e:
        raise HTTPException(status_code=400, detail="Invalid job ID format") from e

    job_dir = ROOT / "output" / "alignments" / job_id
    if not job_dir.exists():
        raise HTTPException(status_code=404, detail="Job results not found")

    sam_file = job_dir / "filtered_probe_alignments_annotated.sam"
    if not sam_file.exists():
        sam_file = job_dir / "filtered_probe_alignments.sam"
    if not sam_file.exists():
        raise HTTPException(status_code=404, detail="Alignment file not found")

    headers = [
        'Probe ID', 'Probe Sequence', 'GC%', 'Target Sequence',
        'Target Transcript', 'Gene ID', 'Species',
        'Mismatches', 'Bulges', 'Position', 'Strand'
    ]

    src_gene, src_transcripts = _load_job_source(job_dir)
    genome_self_match = _load_job_genome_self_match(job_dir)
    risk_map = _load_probe_risk_levels(job_dir) if risk_report else None

    def row_iter():
        yield '\t'.join(headers) + '\n'
        try:
            alignments = parse_sam_file(
                str(sam_file), mismatch_filter=mismatch,
                source_gene=src_gene, source_transcripts=src_transcripts,
                genome_self_match=genome_self_match,
            )
        except Exception:
            logger.exception("Failed to parse alignments for download (job %s)", job_id)
            return
        for a in alignments:
            pid = a.get('probe_id') or ''
            if risk_map is not None:
                if pid not in risk_map and pid.split('|')[0].strip() not in risk_map:
                    continue
            probe_seq = (a.get('sequence') or '').upper().replace('-', '')
            target = a.get('target_transcript') or ''
            gene_id = a.get('gene_id') or ''
            species = a.get('species') or ''
            # When species/gene_id is unclassified, fall back to the target rname
            # so the user gets a meaningful label instead of literal 'Unknown'.
            if not gene_id or gene_id == 'Unknown':
                gene_id = target
            if not species or species == 'Unknown':
                species = target
            row = [
                a.get('probe_id') or '',
                probe_seq,
                _gc_percent(probe_seq),
                a.get('target_sequence') or '',
                target,
                gene_id,
                species,
                str(a.get('substitutions') if a.get('substitutions') is not None else (a.get('mismatches') if a.get('mismatches') is not None else '')),
                str(a.get('bulges') if a.get('bulges') is not None else ''),
                str(a.get('position') if a.get('position') is not None else ''),
                a.get('strand') or '',
            ]
            yield '\t'.join(row) + '\n'

    if risk_report:
        filename = 'medium_high_risk_report.tsv'
    else:
        suffix = 'all' if mismatch is None else f'mm{mismatch}'
        filename = f'probe_alignments_{suffix}.tsv'
    return StreamingResponse(
        row_iter(),
        media_type='text/tab-separated-values',
        headers={'Content-Disposition': f'attachment; filename="{filename}"'}
    )


OFFTARGET_TOP_N_DEFAULT = 12

@app.get("/jobs/{job_id}/offtarget-summary")
def get_offtarget_summary(job_id: str, top_n: int = Query(OFFTARGET_TOP_N_DEFAULT, ge=1, le=100)):
    """
    Lightweight aggregation of off-target alignments grouped by species (microbiome)
    or gene_id (host transcripts). When the species/gene label is 'Unknown', falls
    back to the SAM target reference name (rname). Streams the SAM once without
    parsing CIGAR/MD/sequences, so it's ~17x faster than /alignments for large jobs.
    """
    try:
        uuid.UUID(job_id)
    except ValueError as e:
        raise HTTPException(status_code=400, detail="Invalid job ID format") from e

    job_dir = ROOT / "output" / "alignments" / job_id
    if not job_dir.exists():
        raise HTTPException(status_code=404, detail="Job results not found")

    sam_file = job_dir / "filtered_probe_alignments_annotated.sam"
    if not sam_file.exists():
        sam_file = job_dir / "filtered_probe_alignments.sam"
    if not sam_file.exists():
        return {
            "job_id": job_id,
            "group_by": None,
            "total_alignments": 0,
            "total_groups": 0,
            "hidden_groups": 0,
            "excluded_unknown": 0,
            "groups": [],
        }

    # Self-alignment exclusion: when the pipeline identified a source gene/
    # transcripts (e.g. host-internal mode), drop those hits from the
    # off-target tally so the chart reflects truly off-target alignments only.
    source_gene, source_transcripts_list = _load_job_source(job_dir)
    source_transcripts = set(source_transcripts_list)
    source_transcripts_versionless = {t.split('.', 1)[0] for t in source_transcripts}
    genome_self_match = _load_job_genome_self_match(job_dir)

    sp_pattern = re.compile(r'\tSP:Z:([^\t\n]+)')
    counts: Dict[str, int] = {}
    excluded_unknown = 0
    excluded_self = 0
    total = 0
    group_by = None
    seen_keys: set = set()
    try:
        with open(sam_file, 'r', encoding='utf-8') as f:
            for line in f:
                if line.startswith('@') or not line.strip():
                    continue
                parts = line.split('\t', 11)
                if len(parts) < 5:
                    continue
                # Annotated host SAM has gene_id in column 4 (non-numeric); microbiome has SP:Z tag.
                col4 = parts[3]
                rname = parts[2]
                qname = parts[0]
                try:
                    flag = int(parts[1])
                except ValueError:
                    flag = 0
                key = None
                try:
                    pos = int(col4)
                except ValueError:
                    key = col4 if col4 and col4 != '*' else None
                    if key is not None and group_by is None:
                        group_by = 'gene'
                    try:
                        pos = int(parts[4])
                    except (ValueError, IndexError):
                        pos = -1
                if key is None:
                    m = sp_pattern.search(line)
                    if m:
                        key = m.group(1)
                        if group_by is None:
                            group_by = 'species'
                if key is None or key == 'Unknown':
                    # Fall back to the target reference name when species/gene
                    # is unclassified, so the user gets a meaningful label.
                    key = rname if rname and rname != '*' else None
                if key is None:
                    continue
                # Skip self-alignments to the source gene/transcript, or (for
                # microbiome probe-input) to the probe's own genome.
                rname_base = rname.split('.', 1)[0] if rname else ''
                is_self = (
                    (genome_self_match and same_genome(qname, rname)) or
                    (source_gene and key == source_gene) or
                    (rname and rname in source_transcripts) or
                    (rname_base and rname_base in source_transcripts_versionless)
                )
                if is_self:
                    excluded_self += 1
                    continue
                strand = '-' if (flag & 0x10) else '+'
                dedup_key = (qname, rname, pos, strand)
                if dedup_key in seen_keys:
                    continue
                seen_keys.add(dedup_key)
                total += 1
                counts[key] = counts.get(key, 0) + 1
    except OSError as e:
        logger.exception("Failed to read SAM for job %s", job_id)
        raise HTTPException(status_code=500, detail="Failed to read alignment file") from e

    sorted_groups = sorted(counts.items(), key=lambda kv: -kv[1])
    top = sorted_groups[:top_n]
    return {
        "job_id": job_id,
        "group_by": group_by or 'species',
        "total_alignments": total,
        "total_groups": len(sorted_groups),
        "hidden_groups": max(0, len(sorted_groups) - len(top)),
        "excluded_unknown": excluded_unknown,
        "groups": [{"label": label, "count": count} for label, count in top],
    }


@app.get("/jobs/{job_id}/download/{filename}")
def download_job_file(job_id: str, filename: str):
    """Download a specific output file from a completed job."""
    
    try:
        uuid.UUID(job_id)
    except ValueError as e:
        raise HTTPException(status_code=400, detail="Invalid job ID format") from e

    if "/" in filename or "\\" in filename or filename.startswith("."):
        raise HTTPException(status_code=400, detail="Invalid filename")
    
    job_dir = ROOT / "output" / "alignments" / job_id
    
    if not job_dir.exists():
        raise HTTPException(status_code=404, detail="Job results not found")
    
    file_path = job_dir / filename
    
    allowed_files = [
        "probe_alignments.sam",
        "filtered_probe_alignments.sam",
        "filtered_probe_alignments_annotated.sam", 
        "safe_probes_scores.txt",
        "filtered_probe_alignments.bam",
        "filtered_probe_alignments.bam.bai",
        "probe_alignments.bam",
        "probe_alignments.bam.bai",
        "non_aligned_probes.fa",
        "safe_probes.fa",
        "safe_probes_scores.txt",
        "candidate_probes.fa",
        "reference.fasta",
        "reference.fasta.fai",
        "probes.gff3",
        "kmer_matches_report.txt"
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