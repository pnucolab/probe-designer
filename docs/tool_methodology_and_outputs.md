# Microbe Genome Probe Design Pipeline

## Overview

The Microbe Genome Probe Design Pipeline is a web-based application for designing and analyzing oligonucleotide probes for gene sequences. It supports both gene sequence input (for probe generation) and direct probe sequence input (for alignment analysis).

## Architecture

### Technology Stack

**Backend:**
- FastAPI (Python web framework)
- Celery (distributed task queue)
- Redis (message broker)
- RazerS3 (sequence alignment tool)

**Frontend:**
- Svelte (reactive UI framework)
- Caddy (web server and reverse proxy)
- JBrowse 2 (genome browser integration)

**Containerization:**
- Docker & Docker Compose
- Rootless Docker support

---

## Backend Functionality

### Main Components

#### 1. FastAPI Application (`backend/main_app.py`)

**Purpose:** REST API for job submission, status tracking, and result retrieval.

**Key Endpoints:**

- **POST `/jobs`** - Submit new probe design job
  - Parameters:
    - `species`: Host organism (human/mouse)
    - `probe_length`: Length of probes in base pairs
    - `max_mismatches`: Maximum allowed mismatches (0-2)
    - `gene_sequence`: Gene sequence in FASTA format (optional)
    - `probe_sequence`: Probe sequences in FASTA format (optional)
  - Returns: Job ID and initial status

- **GET `/jobs/{job_id}`** - Get job status and results
  - Returns: Job status, execution time, output files

- **GET `/jobs/{job_id}/files`** - List available result files
  - Returns: Array of files with names, sizes, and paths

- **GET `/jobs/{job_id}/files/{filename}`** - Download specific result file
  - Returns: File content with appropriate headers

- **GET `/jobs/{job_id}/alignments`** - Get alignment statistics
  - Query params: `page`, `page_size`, `mismatch_filter`
  - Returns: Paginated alignment data with statistics

- **GET `/jobs/{job_id}/probe/{probe_id}/alignments`** - Get alignments for specific probe
  - Returns: Detailed alignment data for a single probe

**Directory Structure:**
```
microbe_genome/
├── gene_sequences/          # Stored gene input files
├── probe_sequences/         # Stored probe input files
└── outputs/
    └── alignments/
        └── {job_id}/        # Job-specific output directory
            ├── chroms/      # Per-chromosome alignment files
            ├── jbrowse/     # JBrowse visualization files
            └── *.sam        # Alignment results
```

#### 2. Celery Worker (`backend/celery_worker.py`)

**Purpose:** Asynchronous task execution for pipeline jobs.

**Task: `run_pipeline_task`**

Workflow:
1. Read input sequence from saved file
2. Execute `simple_pipeline.py` with appropriate parameters
3. Monitor pipeline execution and capture logs
4. Generate JBrowse visualization files
5. Return job status and results

**Configuration:**
- Broker: Redis
- Result Backend: Redis
- Task time limit: 3600 seconds (1 hour)
- Concurrency: 2 workers

#### 3. Pipeline Script (`simple_pipeline.py`)

**Purpose:** Core probe design and alignment pipeline.

**Workflow:**

For **Gene Sequence Input:**
1. Generate tiled probes from gene sequence
2. Align probes against reference genome
3. Parse alignment results (SAM format)
4. Filter and annotate alignments

For **Probe Sequence Input:**
1. Read pre-designed probes
2. Align against reference genome
3. Parse and analyze alignments

**Key Functions:**
- Genome download and indexing (GENCODE)
- Probe generation with configurable parameters
- RazerS3 alignment execution
- SAM file parsing and statistics

#### 4. JBrowse Generator (`backend/jbrowse_generator.py`)

**Purpose:** Create genome browser visualization files.

**Generated Files:**
- `reference.fa` - Reference sequence in FASTA format
- `probes.gff3` - Probe annotations in GFF3 format
- `alignments.bam` - Binary alignment file (if available)
- `alignments.bam.bai` - BAM index file

#### 5. SAM Parser (`backend/sam_parser.py`)

**Purpose:** Parse and analyze SAM alignment files.

**Functions:**
- `parse_sam_file()` - Extract alignment records with filtering
- `count_sam_alignments()` - Count total alignments
- Mismatch detection and annotation
- Position and strand information extraction

---

## Frontend Functionality

### Main Components

#### 1. Main Application (`frontend/src/App.svelte`)

**Purpose:** Root component managing application state and routing.

**Features:**
- Job submission and tracking
- Navigation between upload and status views
- Error handling and notifications

#### 2. Upload Form (`frontend/src/components/UploadForm.svelte`)

**Purpose:** User interface for job submission.

**Features:**

**Input Types:**
- Gene Sequence: Paste FASTA sequence to generate probes
- Probe Sequence: Paste pre-designed probes for alignment

**Configuration Options:**
- Host Organism: Human or Mouse
- Probe Length: Configurable in base pairs (gene input only)
- Max Mismatches: 0-2 allowed mismatches

**Default Sequences:**
- Pre-populated with example gene/probe sequences
- Users can paste custom sequences or use defaults

**Validation:**
- FASTA format validation
- Required field checking
- Real-time error feedback

#### 3. Job Status (`frontend/src/components/JobStatus.svelte`)

**Purpose:** Display job progress and results.

**Features:**

**Status Display:**
- PENDING: Job queued
- PROGRESS: Pipeline running
- SUCCESS: Job completed
- FAILURE: Error occurred

**Results View:**
- Alignment statistics (total, unique locations)
- Per-chromosome breakdown
- Downloadable result files

**Alignment Browser:**
- Paginated alignment table (25 per page)
- Mismatch filtering
- Probe-specific alignment details
- Search by probe ID

**File Downloads:**
- SAM files (per chromosome)
- Summary statistics
- JBrowse configuration files

#### 4. JBrowse Viewer (`frontend/src/components/JBrowseViewer.svelte`)

**Purpose:** Embedded genome browser for visual analysis.

**Features:**
- Interactive genome visualization
- Probe track display
- Alignment track (if available)
- Pan and zoom navigation
- Track configuration

---

## Deployment

### Docker Compose Setup

**Services:**

1. **redis** - Message broker for Celery
   - Image: `redis:7`
   - Port: 6379 (internal)

2. **probe-designer-api** - FastAPI backend
   - Build: `./Dockerfile`
   - Port: 8000 (internal)
   - Volumes: Code, data, outputs

3. **probe-designer-worker** - Celery worker
   - Build: `./Dockerfile`
   - Command: Celery worker with 2 concurrent tasks
   - Volumes: Shared with API

4. **probe-designer-frontend** - Caddy + Svelte
   - Build: `./frontend/Dockerfile`
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
3. Backend validates input and generates job ID
4. Input sequence saved to `gene_sequences/` or `probe_sequences/`
5. Celery task queued with job parameters
6. Frontend polls `/jobs/{job_id}` for status updates

### Pipeline Execution Flow

1. Celery worker picks up task
2. Worker executes `simple_pipeline.py` subprocess
3. Pipeline downloads/indexes genome (if needed)
4. Generates probes (gene input) or reads probes (probe input)
5. Aligns probes using RazerS3
6. Parses SAM output and generates statistics
7. Creates JBrowse visualization files
8. Updates job status to SUCCESS/FAILURE

### Results Retrieval Flow

1. Frontend requests `/jobs/{job_id}/files`
2. Backend lists available output files
3. User can download individual files via `/jobs/{job_id}/files/{filename}`
4. Alignment data accessible via `/jobs/{job_id}/alignments`
5. JBrowse viewer loads visualization files

---

## File Formats

### Input Formats

**Gene Sequence (FASTA):**
```
>transcript:ENSB:ID CDS=start-end
ATGCGATCGATCG...
```

**Probe Sequence (FASTA):**
```
>probe_0|start=1|end=20|transcript:ID
ATGCGATCGATCG
>probe_1|start=21|end=40|transcript:ID
GCTAGCTAGCTAG
```

### Output Formats

**SAM Files:**
- Standard sequence alignment format
- One file per chromosome
- Contains alignment positions, mismatches, quality scores

**GFF3 Files:**
- Probe annotations for genome browser
- Features: probe positions, IDs, attributes

**BAM Files:**
- Binary compressed alignment format
- Indexed for fast random access
- Used by JBrowse for visualization

---

## Error Handling

### Backend Errors

- **400 Bad Request**: Invalid input format or parameters
- **404 Not Found**: Job ID or file not found
- **500 Internal Server Error**: Pipeline execution failure

### Frontend Error Display

- Form validation errors shown inline
- API errors displayed in error boxes
- Network errors caught and presented to user
- Loading states during asynchronous operations

---

## Performance Considerations

### Backend Optimizations

- Asynchronous task processing via Celery
- Pagination for large alignment datasets
- File streaming for downloads
- Connection pooling for Redis

### Frontend Optimizations

- Lazy loading of alignment data
- Virtual scrolling for large tables
- Debounced API calls
- Efficient Svelte reactivity

---

## Troubleshooting

### Common Issues

**Directory Creation Errors:**
- Ensure `gene_sequences/` and `probe_sequences/` directories exist
- Backend creates them automatically with `parents=True, exist_ok=True`

**Celery Worker Not Processing:**
- Check Redis connection
- Verify worker is running: `docker-compose logs worker`
- Check broker URL configuration

**JBrowse Not Loading:**
- Verify JBrowse files generated in `outputs/alignments/{job_id}/jbrowse/`
- Check file permissions
- Review browser console for errors

**Docker Shim Errors:**
- Restart Docker daemon: `systemctl --user restart docker` (rootless)
- Remove old containers: `docker-compose down -v`
- Rebuild: `docker-compose up --build -d`

---

## Development

### Local Setup Without Docker

**Backend:**
```bash
pip install -r requirements.txt
redis-server &
celery -A backend.celery_worker worker --loglevel=info &
uvicorn backend.main_app:app --reload --port 8000
```

**Frontend:**
```bash
cd frontend
npm install
npm run dev
```

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

---

## Future Enhancements

- Support for additional organisms
- Batch job submission
- Job history and management
- Advanced filtering and search
- Export results in multiple formats
- User authentication and authorization
- Pipeline parameter optimization
- Real-time progress updates via WebSocket

---

## Support

For issues or questions, please refer to:
- Project repository documentation
- API interactive docs at `/docs` endpoint
- Log files in `outputs/alignments/{job_id}/`

---