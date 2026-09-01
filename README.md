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

> **Note:** Genome data is not downloaded automatically. You must obtain and
place the reference files yourself before running the tool — see
[Preparing Genome Data](#preparing-genome-data) below.

### Frontend

```bash
cd frontend
conda activate probe-designer
npm install
npm run dev
```

## Preparing Genome Data

All genome data lives under `data/` and is config-driven through
`config/organisms.yml`. Adding or updating an organism requires editing that
file plus placing the matching files on disk — no code changes. Genome files
are **not** downloaded automatically; you obtain and place them yourself.

### 1. The organism registry (`config/organisms.yml`)

Each host entry declares where its files live and how to find them:

| Field | Meaning |
|-------|---------|
| `genome_dir` | Directory holding the genome FASTA + GTF |
| `genome_file_regex` | Python regex matching the genome FASTA; capture group 1 = release number |
| `transcript_chunks_dir` | Directory of split transcript chunks aligned against |
| `microbiomes` | Microbiome catalogs hostable on this organism (each has `id`, `data_dir`, `source_url`) |

### 2. Host genomes

Place files in the host's `genome_dir` (GENCODE hosts under
`data/gencode_data/<id>/`, Ensembl hosts under `data/ensembl_data/<id>/`):

- Genome FASTA named to match `genome_file_regex`
  (e.g. `GRCh38.p14.genome.release_49.fa`)
- Matching annotation GTF (same release number)
- Transcript FASTA (e.g. `..._transcripts.fa`)

Sources: [GENCODE](https://www.gencodegenes.org/) (human, mouse),
[Ensembl](https://www.ensembl.org/) (other organisms).

Split the transcript FASTA into chunks of roughly 7–8 MB each and place them in
the organism's `transcript_chunks_dir`. The pipeline aligns chunks in parallel,
which is what keeps a whole-transcriptome scan tractable; for the human
transcriptome this works out to around 90–100 files.

Only two things are required of the split:

- **Every chunk file ends in `.fa`.** The pipeline collects chunks by listing
  the directory and taking that extension — the filenames themselves are never
  parsed, so any naming scheme works.
- **No FASTA record is split across two chunks.** Break only at record
  boundaries, i.e. immediately before a `>` header line.

Any splitter that respects record boundaries will do (`seqkit split2`, or a
short script that opens the next file once the current one passes the size
target and the next line starts with `>`).

### 3. Microbiome catalogs

Place one FASTA per genome (`.fa`/`.fna`) in the catalog's `data_dir`
(e.g. `data/gut-microbe/`). Catalog sources are linked from each
`source_url` in `organisms.yml` (e.g. EBI MGnify genome catalogues).
An optional `genomes-all_metadata.tsv` adds species annotation.

#### Custom catalog for a known community

A full public catalog holds thousands of genomes and treats every one of them as
potentially present, so probes are rejected for cross-reactivity with organisms
that may not occur in your samples at all. If you know which organisms are in
your system — from 16S or shotgun profiling, or because it is a defined
community — screen against just those genomes instead:

1. Create `data/<your-catalog-id>/` and place one FASTA per reference genome in
   it (same layout as any other catalog).
2. Register it under the relevant host's `microbiomes:` list in
   `config/organisms.yml`:

   ```yaml
   - id: my-community
     display_name: My Sample Community
     data_dir: data/my-community
   ```

3. Build its k-mer databases as in step 4 below
   (`python build_kmer_dbs.py --species my-community`).
4. Select it in the web UI, or pass `--microbiomes my-community
   --align-microbiome` on the command line.

Probe yield is typically much higher, because a probe is only rejected for
organisms that are actually there. The trade-off is transferability: probes
validated against a narrow catalog may cross-react in samples whose community
differs, so a catalog-wide screen remains the conservative choice when the
community is unknown or variable between subjects.

### 4. Build k-mer databases

Off-target screening uses pre-built Jellyfish k-mer databases. Build them once
per species (re-used on later runs):

```bash
# Build all k-mer lengths (14–19) for every species
python build_kmer_dbs.py

# One species only
python build_kmer_dbs.py --species gut-microbe

# Specific k-mer lengths / thread count
python build_kmer_dbs.py --species human mouse --kmer-lengths 16 17 18 --threads 40
```

Outputs `kmer_<k>mer_counts.jf` next to each species' genome files. Requires
[Jellyfish](https://github.com/gmarcais/Jellyfish) on your `PATH`.

### Expected layout

```
data/
├── gencode_data/
│   ├── human/
│   │   ├── GRCh38.p14.genome.release_49.fa
│   │   ├── gencode.v49.basic.annotation.release_49.gtf
│   │   ├── GRCh38.p14.genome.release_49_transcripts.fa
│   │   └── transcript_chunks/
│   │       ├── ...part000.fa ... part093.fa
│   │       └── kmer_14mer_counts.jf ... kmer_19mer_counts.jf
│   └── mouse/ ...
├── ensembl_data/            # zebrafish, drosophila, ... (same shape)
└── gut-microbe/             # one FASTA per genome + kmer_*.jf
    ├── MGYG000000001.fa
    └── genomes-all_metadata.tsv
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

**Submit job:**
```
POST /jobs
Form fields: gene_sequence | probe_sequence  (exactly one),
             species, microbiomes, mode, probe_length, max_mismatches,
             max_bulges, kmer_length, tm_range, gc_range,
             align_host, align_microbiome
```

**Job status:**
```
GET /jobs/{job_id}
```

**List result files:**
```
GET /jobs/{job_id}/files
```

**Alignments (paginated):**
```
GET /jobs/{job_id}/alignments?page=1&page_size=25&mismatch=0&probe_id=probe_1
GET /jobs/{job_id}/alignments/download          # TSV; add ?risk_report=true
```

**Off-target summary:**
```
GET /jobs/{job_id}/offtarget-summary?top_n=12
```

**Download a result file:**
```
GET /jobs/{job_id}/download/{filename}
```

**Organism registry / host-token check:**
```
GET  /api/organisms
POST /validate-host-token       Form fields: fasta, species
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
├── core/                     # Core utilities
├── docs/                     # Documentation
├── probe_designer.py         # Main pipeline script
├── docker-compose.yml
└── requirements.txt
```

## Output Files

Each job writes these to `output/alignments/{job_id}/`:

| File | Contents |
| --- | --- |
| `safe_probes_scores.txt` | Per-probe metrics table (Tm, GC%, length, complexity, secondary structure) |
| `kmer_matches_report.txt` | K-mer safety check results for non-aligned probes |
| `candidate_probes.fa` | Probes after GC, Tm and homopolymer filtering |
| `non_aligned_probes.fa` | Probes with no alignment to the reference |
| `probe_alignments.sam` / `.bam` / `.bam.bai` | All alignments |
| `filtered_probe_alignments.sam` / `.bam` / `.bam.bai` | Alignments within the mismatch budget |
| `filtered_probe_alignments_annotated.sam` | Filtered alignments with gene annotations |
| `probes.gff3` | Probe positions for the JBrowse 2 track |
| `reference.fasta` / `.fai` | Reference sequence and index used by JBrowse 2 |

## License

GNU General Public License v3.0
