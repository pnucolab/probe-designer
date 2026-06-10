<script>
  import { onMount } from 'svelte';
  import UploadForm from './components/UploadForm.svelte';
  import JobStatus from './components/JobStatus.svelte';
  import PageHeader from './components/PageHeader.svelte';
  import Bottom from './components/bottom.svelte';
  import Navigation from './components/Navigation.svelte';
  import LandingPage from './components/LandingPage.svelte';
  import Button from './components/Button.svelte';
  import './styles/common.css';
  import './app.css';

  let currentView = 'home';
  let currentJob = null;
  let designMode = 'microbe';

  onMount(() => {
    try {
      const params = new URLSearchParams(window.location.search);
      const jid = params.get('job_id');
      const view = params.get('view');
      const mode = params.get('mode');
      if (mode === 'host' || mode === 'microbe') {
        designMode = mode;
      }

      if (jid) {
        currentJob = { job_id: jid };
        currentView = 'results';
      } else if (view === 'about') {
        currentView = 'about';
      } else if (view === 'upload') {
        currentView = 'upload';
      } else {
        currentView = 'home';
      }
    } catch (e) {
      currentView = 'home';
    }
  });

  function handleSubmitted(e) {
    const jobData = e.detail;
    if (jobData && jobData.job_id) {
      currentJob = jobData;
      currentView = 'results';
      window.history.pushState({}, '', `/?job_id=${encodeURIComponent(jobData.job_id)}`);
    }
  }

  function handleNavigate(e) {
    const { view, mode } = e.detail;
    currentView = view;
    if (mode === 'host' || mode === 'microbe') {
      designMode = mode;
    }

    if (view === 'home') {
      window.history.pushState({}, '', '/');
    } else if (view === 'upload') {
      window.history.pushState({}, '', `/?view=upload&mode=${designMode}`);
    } else if (view === 'results' && currentJob && currentJob.job_id) {
      window.history.pushState({}, '', `/?job_id=${encodeURIComponent(currentJob.job_id)}`);
    } else if (view === 'about') {
      window.history.pushState({}, '', '/?view=about');
    }
  }

  function handleGetStarted(mode) {
    if (mode === 'host' || mode === 'microbe') {
      designMode = mode;
    }
    handleNavigate({ detail: { view: 'upload' } });
    // Reflect the mode in the URL so reloads keep the same form.
    const params = new URLSearchParams(window.location.search);
    params.set('view', 'upload');
    params.set('mode', designMode);
    window.history.replaceState({}, '', `/?${params.toString()}`);
  }
</script>

<Navigation {currentView} hasJob={!!currentJob} activeMode={designMode} on:navigate={handleNavigate} />

<main class="main-container">
  {#if currentView === 'home'}
    <LandingPage
      onGetStarted={handleGetStarted}
      onAbout={() => handleNavigate({ detail: { view: 'about' } })}
    />

  {:else if currentView === 'upload'}
    <PageHeader />
    <div class="content-wrapper">
      <UploadForm mode={designMode} on:submitted={handleSubmitted} />
    </div>

  {:else if currentView === 'results'}
    <PageHeader showJobResults={true} />
    <div class="content-wrapper">
      <JobStatus jobId={currentJob.job_id} />
    </div>

  {:else if currentView === 'about'}
    <div class="section">
      <h1 class="heading-1">About SHARP-FISH</h1>

      <section>
        <h2 class="heading-2">Overview</h2>
        <p class="text-body">
          SHARP-FISH is an interactive web platform for designing oligonucleotide probes
          that target host or microbial genes with high specificity. Users can provide
          host transcripts, microbial sequences, or pre-existing probe sets, and the tool
          screens each candidate against the relevant host transcriptome and the co-residing
          microbes present in complex microbial communities, removing candidates with up to
          two mismatches and a stretch of greater than or equal to 14 consecutive matches to
          off-target sequences. The final probe sets can be applied to detect host or
          microbial transcripts within complex, host-associated tissues.
        </p>
      </section>

      <section>
        <h2 class="heading-2">How It Works</h2>
        <div class="info-box">
          <ol class="ordered-list">
            <li><strong>Sequence Input:</strong> Users paste either microbial gene sequences (for probe generation) or probe sequences (for direct alignment analysis) in FASTA format directly into the web interface</li>
            <li><strong>Probe Generation:</strong> The pipeline generates short oligonucleotide probes of specified length from the input genome</li>
            <li><strong>Alignment:</strong> Each probe is aligned against the selected host transcriptome and/or co-residing microbiome communities (e.g. gut, oral, skin, vaginal) using high-performance alignment algorithms</li>
            <li><strong>Specificity Filtering:</strong> Probes that match the host or microbiome transcriptomes within the specified mismatch tolerance are flagged as potential cross-hybridizers and filtered out</li>
            <li><strong>Results & Visualization:</strong> For gene sequence input, results are displayed in an interactive genome browser showing probe positions and alignment details.
              For probe sequence input, a detailed alignment table shows all matches with
              mismatch counts and filtering options</li>
          </ol>
        </div>
      </section>

      <section>
        <h2 class="heading-2">How to Use</h2>
        <div class="grid-auto">
          <a href="https://probe-designer-rtd.readthedocs.io/en/latest/index.html" target="_blank" rel="noopener" class="doc-link-card">
            <div class="doc-link-row">
              <div class="doc-link-icon">
                <svg width="22" height="22" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.5" stroke-linecap="round" stroke-linejoin="round">
                  <path d="M2 3h6a4 4 0 0 1 4 4v14a3 3 0 0 0-3-3H2z"/>
                  <path d="M22 3h-6a4 4 0 0 0-4 4v14a3 3 0 0 1 3-3h7z"/>
                </svg>
              </div>
              <div>
                <h3 class="doc-link-title">Documentation</h3>
                <p class="doc-link-desc">Usage guides, setup instructions, and best practices</p>
              </div>
            </div>
          </a>
          <a href="https://probe-designer-rtd.readthedocs.io/en/latest/userguide/07-output-interpretation.html" target="_blank" rel="noopener" class="doc-link-card">
            <div class="doc-link-row">
              <div class="doc-link-icon">
                <svg width="22" height="22" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.5" stroke-linecap="round" stroke-linejoin="round">
                  <path d="M4 4h16v16H4z"/>
                  <line x1="8" y1="8" x2="16" y2="8"/>
                  <line x1="8" y1="12" x2="16" y2="12"/>
                  <line x1="8" y1="16" x2="12" y2="16"/>
                </svg>
              </div>
              <div>
                <h3 class="doc-link-title">Input & Output Interpretation</h3>
                <p class="doc-link-desc">File formats, probe scores, and result analysis</p>
              </div>
            </div>
          </a>
        </div>
      </section>

      <section>
        <h2 class="heading-2">Version History</h2>
        <p class="text-body">
          See the <a href="https://probe-designer-rtd.readthedocs.io/en/latest/changelog.html" target="_blank" rel="noopener" style="color: #3b82f6; text-decoration: underline;">changelog</a> for a full list of releases, features, and updates.
        </p>
      </section>

      <section>
        <h2 class="heading-2">Use Cases</h2>
        <ul class="unordered-list">
          <li>Designing specific oligonucleotide probes for microbial spatial transcriptomics within host tissues</li>
          <li>Single-cell transcriptomics probe design for host-associated microbiome studies</li>
          <li>Screening existing probe sets for cross-hybridization against host and co-residing microbiome genomes</li>
          <li>Quality control of probe specificity before experimental use</li>
        </ul>
      </section>

      <section class="info-box-blue">
        <p class="blue-box-text"><strong>Citation:</strong> If you use this tool in your research, please cite</p>
      </section>

      <section>
        <h2 class="heading-2">Getting Started</h2>
        <p class="text-body" style="margin-bottom: 1rem;">
          Pick a probe design mode below. SHARP-FISH will validate every candidate
          against the relevant host transcriptome and/or co-residing microbiome
          catalogs in minutes.
        </p>
        <div style="display:flex; flex-wrap:wrap; gap:12px;">
          <Button
            variant="gradient"
            size="md"
            onClick={() => handleGetStarted('host')}
          >
            Within-Organism Probe Design →
          </Button>
          <Button
            variant="gradient"
            size="md"
            onClick={() => handleGetStarted('microbe')}
          >
            Microbial Probe Design →
          </Button>
        </div>
      </section>
    </div>
  {/if}
  
  <Bottom />
</main>