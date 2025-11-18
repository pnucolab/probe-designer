# Probe Designer Backend Startup Instructions

To start the application, you must run both the Celery worker and the FastAPI server. These steps are mandatory:

1. **Start the Celery worker:**
   ```sh
   celery -A combined_tasks.celery_app worker --loglevel=info
   ```

2. **Start the FastAPI server:**
   ```sh
   uvicorn combined_api:app --host 0.0.0.0 --port 8003 --reload
   ```

Make sure you have all dependencies installed and Redis running before starting these services.
# Probe Designer - Clean Structure

This is a reorganized version of the probe designer pipeline with a clean, maintainable directory structure.

## Quick Start

1. **Setup Environment**:
   ```bash
   bash setup.sh
   ```

2. **Run Pipeline**:
   ```bash
   python simple_pipeline.py --transcript-id ENST00000634225.2 --species human
   ```

3. **Start Web API** (optional):
   ```bash
   cd backend && python start_api.py
   ```

## Directory Structure

```
probe_designer_clean/
├── 📄 README.md                    # This file
├── 🚀 simple_pipeline.py           # Main pipeline entry point
├── ⚙️  setup.sh                    # Environment setup script
├── 📋 requirements.txt            # Python dependencies
├── 🔨 CMakeLists.txt             # Build configuration
│
├── 🔬 core/                      # Core pipeline modules
│   ├── gencode_genome.py         # Gencode data download and management
│   ├── generate_probes.py        # Probe generation logic
│   ├── chunk_transcripts.py      # Transcript chunking for parallel processing
│   ├── razers_test_on_chunktranscripts.py # Alignment testing
│   ├── filter_unique_probes.py   # Probe uniqueness filtering
│   └── config.py                 # Configuration settings
│
├── 🌐 backend/                   # Web API and backend services
│   ├── app.py                    # FastAPI web application
│   ├── celery_worker.py          # Celery task worker
│   ├── simple_api.py             # Simplified API endpoint
│   ├── main_app.py              # Application entry point
│   ├── start_api.py             # API startup script
│   ├── start_celery.sh          # Celery worker startup script
│   ├── requirements.txt         # Backend-specific dependencies
│   ├── conda_environment.yml    # Conda environment specification
│   └── *.py                     # Other backend utilities
│
├── ⚡ bin/                      # Binary executables
│   └── razers3                  # Sequence alignment tool
│
├── 💾 data/                     # Data storage
│   └── gencode_data/            # Downloaded Gencode reference data
│       ├── human/               # Human genome data
│       └── mouse/               # Mouse genome data
│
├── 📚 docs/                     # Documentation
│   ├── CLI_USAGE_GUIDE.md
│   ├── INTEGRATED_PIPELINE_GUIDE.md
│   └── README.md
│
└── 📊 outputs/                  # Pipeline outputs and results
    ├── split_outputs/           # Chunked transcript files
    ├── split_outputs_probe_alignments/ # Alignment results
    └── files/                   # Job-specific output files
```

## Usage

### Command Line Interface
```bash
# Run from the probe_designer_clean directory
python simple_pipeline.py --transcript-id ENST00000634225.2 --species human

# With custom probe length
python simple_pipeline.py --transcript-id ENST00000634225.2 --species human --probe-length 25

# Multiple transcripts
python simple_pipeline.py --transcript-id ENST00000634225.2,ENST00000700202.2 --species human
```

### Web API
```bash
# Start the backend API server
cd backend
python start_api.py

# Start Celery workers (in separate terminal)
cd backend
bash start_celery.sh
```

## Installation

1. Install Python dependencies:
```bash
pip install -r requirements.txt
```

2. For the backend API, install additional dependencies:
```bash
cd backend
pip install -r requirements.txt
```

## Notes

- The `bin/razers3` executable must be in the PATH or the scripts will look for it in the `bin/` directory
- Gencode data is automatically downloaded to `data/gencode_data/` when needed
- Pipeline outputs are stored in the `outputs/` directory
- All import paths have been adjusted to work with the new directory structure
