
# Probe Pipeline  Frontend

This frontend is a small Svelte + Vite web UI for the Probe Pipeline API. It lets users upload (optionally gzipped) FASTA files, submit a probe-generation job to the backend, monitor job status, and download resulting files when the job completes. The UI is designed for lightweight local development and for quickly exercising the pipeline's upload and status APIs.

Minimal Svelte + Vite frontend for the Probe Pipeline API.

This document explains how to run the frontend locally using the Conda environment named `probeset` (the project standard), how the dev server proxies API calls to the backend, and some troubleshooting tips.

## Prerequisites

- Conda (mamba/conda) with an environment named `probeset` available.
- Node.js (v16 or newer recommended) and npm available in your `probeset` environment.

> The project expects the backend API (FastAPI) to be running at `http://localhost:8001` by default. The Vite dev server proxies `/jobs` requests to that address; update `vite.config.js` if your backend is hosted elsewhere.

## Activate Conda environment

Open a terminal and activate the `probeset` environment before running any commands:

```zsh
conda activate probeset
```

(If you use `mamba`: `mamba activate probeset`.)

## Install dependencies

From the repository root or the `frontend/` directory run:

```zsh
cd frontend
npm install --no-audit --no-fund
```

If you encounter peer-dependency resolution problems, you can retry with the legacy npm resolver:

```zsh
npm install --legacy-peer-deps --no-audit --no-fund
```

## Run the dev server

Start the dev server (Vite). The dev server proxies API requests for convenience during development.

```zsh
npm run dev
```

Open the app at:

- http://localhost:5173/

## Build / Preview

To build the production bundle:

```zsh
npm run build
```

Preview the production build locally:

```zsh
npm run preview
```

## How the upload flow works

- The upload form POSTs multipart/form-data to `POST /jobs` on the backend.
- After a successful job submission the frontend redirects to `/?job_id=<task_id>` so the page shows the job status and results. The job status panel polls `/jobs/{job_id}` and lists available files when the job completes.

## Customizing the API target

If your FastAPI backend runs on a different host/port, update the proxy in `vite.config.js`.

Example change in `frontend/vite.config.js`:

```js
proxy: {
	'/jobs': {
		target: 'http://MY_BACKEND_HOST:8001',
		changeOrigin: true
	}
}
```

## Troubleshooting

- "Failed to resolve \"@sveltejs/vite-plugin-svelte\"" or ESM errors: ensure `package.json` contains `"type": "module"` and your Node.js is recent (>=16).
- PostCSS/Tailwind messages: this scaffold runs without Tailwind. If you want Tailwind utilities, re-add `tailwindcss`, `postcss`, and `autoprefixer` and follow Tailwind setup.
- If a file upload returns `400 Bad Request` for gzipped FASTA files, ensure the backend is the updated version (the API accepts `.gz`) and that the file field in the form is named `file` (the form already uses `file`).

## What's in this folder

- `index.html` — app entry
- `src/` — Svelte source files
	- `components/UploadForm.svelte` — upload form
	- `components/JobStatus.svelte` — job status + files list
- `vite.config.js` — dev server + proxy config
- `package.json` — frontend deps and scripts

## Next steps (optional)

- Re-enable Tailwind/Flowbite for consistent, production-ready styling.
- Add client-side validation and upload progress indicator.
- Add E2E tests (curl or Playwright) that upload a gzipped FASTA and verify job submission.

---

If you want I can add a short curl example to this README that demonstrates a gzipped FASTA upload and the expected response. Would you like that?
