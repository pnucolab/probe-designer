# Microbial Probe Designer Pipeline

## Overview

The Microbial Probe Designer Pipeline is a web-based application for designing and analyzing oligonucleotide probes for microbial gene sequences while filtering out cross-hybridizing probes against host genomes. It supports both gene sequence input (for probe generation) and direct probe sequence input (for alignment analysis), with specialized support for microbiome datasets.

## Architecture

### Technology Stack

**Backend:**
- FastAPI (Python web framework)
- Celery (distributed task queue)
- Redis (message broker)
- CUDA GPU + CuPy (sequence alignment; the pipeline is GPU-only)
- Jellyfish (k-mer counting)
- samtools (BAM processing)

**Frontend:**
- Svelte (reactive UI framework)
- Vite (build tool)
- Caddy (web server and reverse proxy)
- JBrowse 2 (genome browser integration)

**Containerization:**
- Docker & Docker Compose
- Rootless Docker support

---

## Backend Functionality

### Main Components

#### 1. FastAPI Application ([`backend/main_app.py`](backend/main_app.py))

**Purpose:** REST API for job submission, status tracking, and result retrieval.

**Key Endpoints:**

- **POST `/jobs`** - Submit new probe design job
  - Parameters:
    - `species`: Host organism (human/mouse/microbiome datasets)
    - `probe_length`: Length of probes in base pairs (default: 30)
    - `max_mismatches`: Maximum allowed mismatches (0-2, default: 2)
    - `gene_sequence`: Gene sequence in FASTA format (optional)
    - `probe_sequence`: Probe sequences in FASTA format (optional)
  - Returns: Job ID and initial status

- **GET `/jobs/{job_id}`** - Get job status and results
  - Returns: Job status, execution time, output files, pipeline statistics

- **GET `/jobs/{job_id}/files`** - List available result files
  - Returns: Array of files with names, sizes, and paths

- **GET `/jobs/{job_id}/download/{filename}`** - Download specific result file
  - Returns: File content with appropriate headers
  - Allowed files: SAM, BAM, BAI, GFF3, FASTA, probe scores, logs

- **GET `/jobs/{job_id}/alignments`** - Get alignment statistics
  - Query params: `page`, `page_size`, `mismatch`, `probe_id`
  - Returns: Paginated alignment data with statistics

- **GET `/jobs/{job_id}/probe-regions`** - Get probe regions for visualization
  - Returns: Probe coordinates and risk classifications

**Directory Structure:**
```
microbe_genome/
├── gene_sequences/          # Stored gene input files
├── probe_sequences/         # Stored probe input files
├── data/                    # Reference genomes and annotations
│   ├── gencode_data/        # GENCODE human/mouse data
│   └── {species}/           # Microbiome genome collections
└── output/
    └── alignments/
        └── {job_id}/        # Job-specific output directory
            ├── chroms/      # Per-chromosome alignment files
            ├── jbrowse/     # JBrowse visualization files
            ├── *.sam        # Alignment results
            ├── *.bam        # Binary alignment files
            └── *.txt        # Statistics and logs
```

#### 2. Celery Worker ([`backend/celery_worker.py`](backend/celery_worker.py))

**Purpose:** Asynchronous task execution for pipeline jobs.

**Task: `run_pipeline_task`**

Workflow:
1. Read input sequence from saved file
2. Execute [`probe_designer.py`](probe_designer.py) with appropriate parameters
3. Monitor pipeline execution and capture logs
4. Generate JBrowse visualization files
5. Return job status and results

**Configuration:**
- Broker: Redis
- Result Backend: Redis
- Task time limit: 3600 seconds (1 hour)
- Concurrency: 2 workers

#### 3. Pipeline Script ([`probe_designer.py`](probe_designer.py))

**Purpose:** Core probe design and alignment pipeline.

**Workflow:**

For **Gene Sequence Input:**
1. Generate tiled probes from gene sequence using [`core/generate_probes.py`](core/generate_probes.py)
2. Apply GC content filtering (40-80%)
3. Align probes against reference genome/transcriptome using the GPU aligner
4. Parse alignment results (SAM format) 
5. Filter alignments by mismatch threshold
6. **14-mer Safety Check** on non-aligned probes
7. Score final safe probes using thermodynamic calculations

For **Probe Sequence Input:**
1. Read pre-designed probes
2. Apply GC content filtering (40-80%)
3. Align against reference genome/transcriptome
4. Parse and analyze alignments
5. Filter alignments by mismatch threshold
6. **14-mer Safety Check** on non-aligned probes
7. Score final safe probes

**Key Functions:**
- Genome download and indexing (GENCODE via [`core/genome_downloader.py`](core/genome_downloader.py))
- Probe generation with GC content filtering
- **14-mer generation and safety analysis**
- GPU alignment execution (RazerS3-compatible SAM output)
- SAM file parsing and statistics
- Gene annotation using cached transcript-to-gene mappings
- **Cross-reactivity detection via k-mer matching**

**Supported Species:**

Hosts (registered in [`config/organisms.yml`](config/organisms.yml)):
- `human`: Human (Homo sapiens, GRCh38)
- `mouse`: Mouse (Mus musculus, GRCm39)
- `zebrafish`: Zebrafish (Danio rerio, GRCz11)
- `drosophila`: Fruit Fly (Drosophila melanogaster, BDGP6)
- `celegans`: Nematode (Caenorhabditis elegans, WBcel235)
- `xtropicalis`: Tropical Clawed Frog (Xenopus tropicalis, UCB_Xtro_10.0)
- `chicken`: Chicken (Gallus gallus, GRCg7b)
- `pig`: Pig (Sus scrofa, Sscrofa11.1)
- `arabidopsis`: Thale Cress (Arabidopsis thaliana, TAIR10)
- `maize`: Maize (Zea mays, Zm-B73-REFERENCE-NAM-5.0)
- `soybean`: Soybean (Glycine max, Glycine_max_v2.1)

Microbiome catalogs:
- `gut-microbe`: Human gut microbiome
- `human-oral-microbiome`: Human oral microbiome
- `human-skin-microbiome`: Human skin microbiome
- `human-vaginal-microbiome`: Human vaginal microbiome
- `mouse-gut-microbiome`: Mouse gut microbiome

New hosts or catalogs are added by editing `config/organisms.yml` — no code changes required.

#### 4. JBrowse Generator ([`backend/jbrowse_generator.py`](backend/jbrowse_generator.py))

**Purpose:** Create genome browser visualization files.

**Generated Files:**
- `reference.fasta` - Reference sequence in FASTA format
- `reference.fasta.fai` - FASTA index file
- `probes.gff3` - Probe annotations in GFF3 format
- `alignments.bam` - Binary alignment file (if available)
- `alignments.bam.bai` - BAM index file

**Features:**
- Risk-based probe coloring (safe/medium_risk/high_risk)
- Automatic SAM to indexed BAM conversion
- Header reconstruction for headerless SAM files

#### 5. SAM Parser ([`backend/sam_parser.py`](backend/sam_parser.py))

**Purpose:** Parse and analyze SAM alignment files.

**Functions:**
- [`parse_sam_file()`](backend/sam_parser.py) - Extract alignment records with filtering
- [`count_sam_alignments()`](backend/sam_parser.py) - Count total alignments
- Mismatch detection and annotation
- Position and strand information extraction
- Support for both standard and annotated SAM formats

#### 6. Core Utilities

**[`core/genome_downloader.py`](core/genome_downloader.py):**
- Downloads GENCODE genome and annotation files
- Filters transcripts by length (≥50bp)
- Chunks filtered transcripts for parallel processing
- Automatic cleanup of outdated releases

**[`core/generate_probes.py`](core/generate_probes.py):**
- Generates tiled probes from input sequences
- GC content filtering (40-80% default)
- Configurable probe length and step size

**[`core/scorer.py`](core/scorer.py):**
- Thermodynamic probe scoring using nearest-neighbor parameters
- Melting temperature calculation
- Secondary structure prediction
- Complexity and homopolymer analysis

**[`core/probe_classifier.py`](core/probe_classifier.py):**
- Classifies probes as safe, medium_risk, or high_risk
- Based on alignment mismatch patterns
- Supports microbiome species annotation

#### 7. K-mer Safety Check ([`probe_designer.py`](probe_designer.py) - `check_14mer_safety()`)

**Purpose:** Final safety verification for non-aligned probes using 14-mer analysis.

**Workflow:**
1. **14-mer Generation:** Extract all possible 14-mers from each non-aligned probe
2. **Parallel Alignment:** Align all 14-mers against reference genome on the GPU
3. **Match Detection:** Identify probes containing 14-mers with genome matches
4. **Risk Classification:** 
   - **Unsafe Probes:** Have one or more 14-mers matching the host genome
   - **Safe Probes:** No 14-mer matches found
5. **Detailed Reporting:** Generate comprehensive match reports with genomic locations

**Features:**
- **High Sensitivity:** 14-mer length ensures detection of potential cross-hybridization sites
- **Parallel Processing:** Multi-threaded 14-mer alignment for performance
- **Gene Annotation:** Maps matches to gene names (human/mouse) or species (microbiome)
- **Deduplication:** Removes duplicate matches for cleaner reporting
- **Position Tracking:** Records exact 14-mer positions within probes and genome locations

**Output Files:**
- `safe_probes.fa` - Final safe probes passing all filters
- `14mer_matches_report.txt` - Detailed match analysis
- `safe_probes_scores.txt` - Thermodynamic metrics of safe probes

**Match Report Format:**
```
14-mer Match Report
==================

Total probes checked: 1,250
Unsafe probes (with 14-mer matches): 156  
Safe probes: 1,094

UNSAFE PROBES WITH MATCH DETAILS:
================================================================

Probe ID: probe_123
Sequence: ATGCGATCGATCGATCGATCGATCGATCGAT
Total matches: 3

14-mer Sequence  Position    Transcript/Location           Gene Name    Genome Position
-------------------------------------------------------------------------------------
ATGCGATCGATCGA   0-13       ENST00000123456.7            GAPDH        12345
CGATCGATCGATCG   8-21       ENST00000789012.3            ACTB         67890
GATCGATCGATCGA   12-25      ENST00000345678.2            TUBB         54321
```

---

## Frontend Functionality

### Main Components

#### 1. Main Application ([`frontend/src/App.svelte`](frontend/src/App.svelte))

**Purpose:** Root component managing application state and routing.

**Features:**
- Multi-view navigation (home, upload, results, about)
- Job submission and tracking
- Error handling and notifications
- Landing page with feature overview

#### 2. Landing Page ([`frontend/src/components/LandingPage.svelte`](frontend/src/components/LandingPage.svelte))

**Purpose:** Welcome page with feature highlights.

**Features:**
- Hero section with tool description
- Feature cards (High Specificity, Fast Processing, Interactive Results)
- Visual genome mixing illustration
- Get Started call-to-action

#### 3. Upload Form ([`frontend/src/components/UploadForm.svelte`](frontend/src/components/UploadForm.svelte))

**Purpose:** User interface for job submission.

**Features:**

**Input Types:**
- Gene Sequence: Paste FASTA sequence to generate probes
- Probe Sequence: Paste pre-designed probes for alignment

**Configuration Options:**
- Host Organism: Human, Mouse, or Microbiome datasets
- Probe Length: 20-50 bp (gene input only, default: 30)
- Max Mismatches: 0-2 allowed mismatches (default: 2)

**Default Sequences:**
- Pre-populated with example gene/probe sequences
- Users can paste custom sequences or use defaults

**Validation:**
- Real-time FASTA format validation
- DNA alphabet checking (ATCGN only)
- Required field validation
- Error feedback with specific line numbers

#### 4. Job Status ([`frontend/src/components/JobStatus.svelte`](frontend/src/components/JobStatus.svelte))

**Purpose:** Display job progress and results.

**Features:**

**Status Display:**
- PENDING: Job queued
- PROGRESS: Pipeline running with detailed steps
- SUCCESS: Job completed with full results
- FAILURE: Error occurred with details

**Results View:**
- Pipeline statistics (candidate probes, alignments, filtered results)
- Probe scoring table with thermodynamic metrics
- Per-chromosome alignment breakdown
- Downloadable result files

**Alignment Browser:**
- Paginated alignment table (25 per page)
- Mismatch filtering (0, 1, or 2 mismatches)
- Probe-specific alignment details
- Search by probe ID functionality
- Expandable probe details with sequence information

**Non-Aligned Probes Display:**
- Probe table with per-probe metrics
- Thermodynamic properties (Tm, GC%, complexity, secondary structure, homopolymer)
- Color-coded metric visualization

**File Downloads:**
- SAM/BAM files (per chromosome and merged)
- Probe scoring tables
- JBrowse configuration files
- Pipeline logs and statistics

#### 5. JBrowse Viewer ([`frontend/src/components/JBrowseViewer.svelte`](frontend/src/components/JBrowseViewer.svelte))

**Purpose:** Embedded genome browser for visual analysis.

**Features:**
- Interactive genome visualization of input sequence
- Probe track display with risk-based coloring
- Alignment track (if available)
- Feature click handling for probe details
- Pan and zoom navigation
- Responsive container with manual resizing

---

## Deployment

### Docker Compose Setup

**Services:**

1. **redis** - Message broker for Celery
   - Image: `redis:7`
   - Port: 6379 (internal)

2. **probe-designer-api** - FastAPI backend
   - Build: [`backend/Dockerfile`](backend/Dockerfile)
   - Port: 8002 (internal)
   - Volumes: Code, data, output
   - Dependencies: samtools, Jellyfish, CUDA/CuPy (GPU aligner)

3. **probe-designer-worker** - Celery worker
   - Build: [`backend/Dockerfile`](backend/Dockerfile)
   - Command: Celery worker with 2 concurrent tasks
   - Volumes: Shared with API

4. **probe-designer-frontend** - Caddy + Svelte
   - Build: [`frontend/Dockerfile`](frontend/Dockerfile)
   - Port: 4000 (exposed)
   - Reverse proxy to backend API

**Network:**
- Internal network for service communication
- Frontend exposed on port 4000

### Environment Variables

```bash
CELERY_BROKER_URL=redis://redis:6379/1
CELERY_RESULT_BACKEND=redis://redis:6379/1
ALLOWED_ORIGINS=http://localhost:4000
```

### Starting the Application

```bash
# Build and start all services
docker-compose up --build -d

# View logs
docker-compose logs -f

# Stop services
docker-compose down
```

**Access:**
- Frontend: http://localhost:4000
- API Docs: http://localhost:4000/docs (via reverse proxy)

---

## Data Flow

### Job Submission Flow

1. User submits form in frontend
2. Frontend sends POST to `/jobs` endpoint
3. Backend validates FASTA input and generates job ID
4. Input sequence saved to `gene_sequences/` or `probe_sequences/`
5. Celery task queued with job parameters
6. Frontend polls `/jobs/{job_id}` for status updates

### Pipeline Execution Flow

1. Celery worker picks up task
2. Worker executes [`probe_designer.py`](probe_designer.py) subprocess
3. Pipeline downloads/indexes genome (if needed) using [`core/genome_downloader.py`](core/genome_downloader.py)
4. **Probe Generation & Initial Filtering:**
   - Generates probes (gene input) or reads probes (probe input)
   - Applies GC content filtering (40-80%)
5. **Primary Alignment & Filtering:**
   - Aligns probes on the GPU (RazerS3-compatible SAM output)
   - Filters alignments by mismatch threshold
   - Annotates alignments with gene names (human/mouse) or species (microbiome)
6. **14-mer Safety Analysis:**
   - Extracts non-aligned probes for safety check
   - Generates all 14-mers from candidate probes
   - Performs parallel 14-mer alignment against genome
   - Identifies and removes probes with matching 14-mers
7. **Final Processing:**
   - Scores remaining safe probes using thermodynamic calculations
   - Classifies probes by risk level
   - Creates JBrowse visualization files
   - Updates job status to SUCCESS/FAILURE

### Results Retrieval Flow

1. Frontend requests `/jobs/{job_id}/files`
2. Backend lists available output files
3. User can download individual files via `/jobs/{job_id}/download/{filename}`
4. Alignment data accessible via `/jobs/{job_id}/alignments`
5. JBrowse viewer loads visualization files for interactive browsing

---

## File Formats

### Input Formats

**Gene Sequence (FASTA):**
```
>transcript:ENSB:ID CDS=start-end
ATGCGATCGATCGATCG...
```

**Probe Sequence (FASTA):**
```
>probe_0|start=1|end=30|transcript:ID
ATGCGATCGATCGATCGATCGATCGATCGAT
>probe_1|start=31|end=60|transcript:ID
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT
```

### Output Formats

**SAM Files:**
- Standard sequence alignment format
- One file per chromosome/transcript chunk
- Contains alignment positions, mismatches, quality scores
- Annotated versions include gene names or species information

**GFF3 Files:**
- Probe annotations for genome browser
- Features: probe positions, IDs, risk levels, mismatch counts
- Color coding for visualization

**BAM Files:**
- Binary compressed alignment format
- Sorted and indexed for fast random access
- Used by JBrowse for visualization

**Probe Scoring Files:**
- Whitespace-delimited text table
- Thermodynamic properties and quality metrics
- Listed in probe order (no composite score or ranking)

**14-mer Match Reports:**
- Text format with detailed match analysis
- Probe-by-probe breakdown of unsafe 14-mer matches
- Genomic coordinates and gene annotations
- Statistics summary of safety check results

---

## Error Handling

### Backend Errors

- **400 Bad Request**: Invalid FASTA format, unsupported characters, missing sequences
- **404 Not Found**: Job ID or file not found
- **403 Forbidden**: Unauthorized file access or path traversal attempts
- **500 Internal Server Error**: Pipeline execution failure, genome download errors

### Frontend Error Display

- Form validation errors shown inline with specific line numbers
- API errors displayed in styled error boxes
- Network errors caught and presented to user
- Loading states during asynchronous operations
- Real-time pipeline progress updates

---

## Performance Considerations

### Backend Optimizations

- Asynchronous task processing via Celery
- Parallel alignment processing using multiprocessing
- Pagination for large alignment datasets (25 records per page)
- File streaming for large downloads
- Connection pooling for Redis
- Efficient SAM parsing with sequence caching

### Frontend Optimizations

- Lazy loading of alignment data
- Virtual scrolling for large probe tables
- Debounced API calls for search functionality
- Efficient Svelte reactivity patterns
- Compressed static assets via Vite

### Pipeline Optimizations

- **Multi-stage Filtering:** Progressive filtering reduces computational load
  1. GC content filtering (fast, eliminates poor probes early)
  2. Primary alignment filtering (medium cost, removes obvious matches)
  3. 14-mer safety check (thorough, ensures no missed cross-reactivity)
- **Parallel 14-mer Processing:** Multi-threaded alignment of k-mers
- **Pre-chunked Transcript Files:** Enables parallel alignment processing
- **Cached Transcript-to-Gene Mappings:** Avoids repeated GTF parsing
- **Incremental Genome Updates:** Cleanup of outdated files
- **Deduplication:** Efficient match reporting without redundancy

---

## Troubleshooting

### Common Issues

**Directory Creation Errors:**
- Ensure `gene_sequences/` and `probe_sequences/` directories exist
- Backend creates them automatically with `parents=True, exist_ok=True`

**Celery Worker Not Processing:**
- Check Redis connection
- Verify worker is running: `docker-compose logs probe-designer-worker`
- Check broker URL configuration

**JBrowse Not Loading:**
- Verify JBrowse files generated in `output/alignments/{job_id}/`
- Check BAM file creation and indexing
- Review browser console for React/JBrowse errors

**Pipeline Failures:**
- Check genome download and decompression
- Verify a CUDA GPU is visible to the worker (`nvidia-smi`) — the aligner is GPU-only
- Review pipeline logs in `output/alignments/{job_id}/pipeline_log.txt`

**14-mer Safety Check Issues:**
- Check 14-mer generation and alignment steps
- Verify sufficient disk space for temporary 14-mer files
- Review 14-mer match reports for detailed error information

**Docker Issues:**
- Restart Docker daemon: `systemctl --user restart docker` (rootless)
- Remove old containers: `docker-compose down -v`
- Rebuild with cache clearing: `docker-compose build --no-cache`

**SAM Processing Errors:**
- Check for missing SAM headers (auto-reconstructed from reference.fasta.fai)
- Verify samtools installation in container
- Review BAM conversion logs

---

## Development

### Local Setup Without Docker

**Backend:**
```bash
pip install -r requirements.txt
redis-server &
celery -A backend.celery_worker worker --loglevel=info &
uvicorn backend.main_app:app --reload --port 8002
```

**Frontend:**
```bash
cd frontend
npm install --legacy-peer-deps
npm run dev
```

**Dependencies:**
- samtools (for BAM processing)
- Jellyfish (for k-mer counting)
- A CUDA-capable GPU with CuPy (the aligner; there is no CPU fallback)

### Testing

**Backend:**
```bash
pytest tests/
```

**Frontend:**
```bash
cd frontend
npm test
```

**Pipeline Testing:**
```bash
python probe_designer.py --gene-sequence ">test\nATCGATCG" --species human --probe-length 20
```

**14-mer Safety Testing:**
```bash
python probe_designer.py --gene-sequence ">test\nATCGATCGATCGATCGATCGATCGATCGAT" --species human --probe-length 30
# Check output/alignments/{job_id}/14mer_matches_report.txt for detailed results
```

---

## Technical Specifications

### Alignment Parameters

- **Identity Threshold:** Calculated as `((probe_length - max_mismatches) / probe_length) * 100`
- **Aligner Options** (retained from the RazerS3-compatible interface): 
  - `-ng`: No gaps allowed
  - `-rr 100`: Report up to 100 alignments per probe
  - `-m 100`: Consider up to 100 matches
  - `-tc 1`: Use 1 thread per alignment job

### K-mer Safety Parameters

- **K-mer Length:** 14 nucleotides (optimal for specificity vs. sensitivity)
- **Alignment Identity:** 100% (perfect matches only for safety)
- **Coverage:** All possible 14-mers extracted from each probe
- **Parallel Workers:** Configurable (default: CPU count)
- **Match Reporting:** Position-specific with gene annotation

### Filtering Cascade

1. **GC Content Filter:** 40-80% GC content (eliminates ~10-20% of probes)
2. **Primary Alignment Filter:** ≤2 mismatches allowed (species-specific filtering)
3. **14-mer Safety Filter:** Zero 14-mer matches required (final safety verification)
4. **Quality Scoring:** Thermodynamic and complexity metrics for ranking

### Probe Scoring Metrics

- **Melting Temperature (Tm):** Nearest-neighbor thermodynamics
- **GC Content:** Percentage of G and C bases
- **Complexity:** Sequence entropy calculation
- **Secondary Structure:** Hairpin formation energy
- **Homopolymer Runs:** Maximum consecutive identical bases

### Risk Classification

- **Safe:** No alignments found in host genome/transcriptome and passes 14-mer safety check
- **Medium Risk:** Alignments with mismatches within threshold but no perfect 14-mer matches
- **High Risk:** Perfect or near-perfect alignments to host or failed 14-mer safety check

---

## Future Enhancements

- Support for additional host organisms (plant, fungal genomes)
- Batch job submission with multiple sequences
- User authentication and job history management
- Advanced probe design algorithms (avoiding repetitive regions)
- Export results in multiple formats (CSV, Excel, BED)
- Real-time progress updates via WebSocket
- Integration with public sequence databases
- Primer design functionality
- Cross-reactivity prediction algorithms
- Configurable k-mer lengths for different sensitivity levels
- Machine learning-based probe quality prediction

---

## Support

For issues or questions, please refer to:
- Project repository documentation
- API interactive docs at `/docs` endpoint  
- Pipeline logs in `output/alignments/{job_id}/pipeline_log.txt`
- 14-mer match reports in `output/alignments/{job_id}/14mer_matches_report.txt`
- JBrowse integration logs in browser console
- Celery worker logs via `docker-compose logs probe-designer-worker`

**License:** GNU General Public License v3.0 (see [`LICENSE`](LICENSE))

**Development:** Computational Omics Laboratory at Pusan National University

---