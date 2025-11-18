# Probe Designer Tool

A web-based tool for probe sequence alignment and visualization with JBrowse 2 browser.

## Features

- Submit gene sequences or probe sequences for analysis
- Align against reference genomes (Human GRCh38 or Mouse GRCm39)
- Interactive alignment visualization with JBrowse 2
- Filter alignments by mismatch count
- Real-time job status monitoring
- Download alignment files (BAM, BAI, SAM)

## Tech Stack

**Backend:**
- FastAPI
- Celery (task queue)
- Redis (message broker)

**Frontend:**
- Svelte
- JBrowse 2
- Vite
- TailwindCSS

## Quick Start with Docker

```bash
docker-compose up -d
```

Access the application:
- Frontend: http://localhost:4000
- Backend API: http://localhost:8000

## Manual Setup

### Backend

```bash
# Create and activate a new conda environment (recommended)
conda create -n probe-designer python=3.11 -y
conda activate probe-designer
# Install dependencies
pip install -r requirements.txt

# Start Celery worker
celery -A backend.celery_worker worker --loglevel=info

# Start FastAPI server
uvicorn backend.main_app:app --host 0.0.0.0 --port 8000 --reload
```

> **Note:** The first run may take a while because the reference FASTA file and annotation 
file for the host are downloaded automatically. On subsequent runs, the tool will use the previously downloaded files unless new releases are available.

### Frontend

```bash
cd frontend
conda activate probe-designer
npm install
npm run dev
```

## Usage

### Web Interface

1. Select input type: Gene Sequence or Probe Sequence
2. Paste your sequence
3. Choose species (Human or Mouse)
4. Submit and monitor job progress
5. View alignments in JBrowse 2
6. Filter by mismatches (0, 1, 2, 3, or All)
7. Search specific probes
8. Download alignment files

### API Endpoints

**Submit Job:**
```
POST /submit-job
Form data: input_type, input_data, species
```

**Job Status:**
```
GET /job-status/{job_id}
```

**Get Alignments:**
```
GET /api/probes/{job_id}?page=1&page_size=25&mismatch_filter=0
```

**Download Files:**
```
GET /download/{job_id}/{filename}
```

## Configuration

Environment variables:
- `CELERY_BROKER_URL` (default: redis://localhost:6379/1)
- `CELERY_RESULT_BACKEND` (default: redis://localhost:6379/1)
- `ALLOWED_ORIGINS` (default: http://localhost:4000)

## Project Structure

```
├── backend/
│   ├── main_app.py          # FastAPI application
│   ├── celery_worker.py     # Celery tasks
│   ├── jbrowse_generator.py # JBrowse file generation
│   └── sam_parser.py        # SAM/BAM parsing
├── frontend/
│   └── src/components/
│       ├── UploadForm.svelte
│       ├── JobStatus.svelte
│       └── JBrowseViewer.svelte
├── bin/                      # Binary executables
├── core/                     # Core utilities
├── docs/                     # Documentation
├── probe_designer.py         # Main pipeline script
├── docker-compose.yml
└── requirements.txt
```

## Output Files

Each job generates:
- `probe_alignments.bam`
- `probe_alignments.bam.bai`
- `filtered_probe_alignments.bam`
- `filtered_probe_alignments.bam.bai`
- `probes.gff3`
- `reference.fasta.fai`
- JBrowse 2 configuration files

## License

GNU General Public License v3.0
