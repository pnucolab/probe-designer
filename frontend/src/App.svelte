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

  let currentView = 'home';
  let currentJob = null;

  onMount(() => {
    try {
      const params = new URLSearchParams(window.location.search);
      const jid = params.get('job_id');
      const view = params.get('view');
      
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
    const { view } = e.detail;
    currentView = view;
    
    if (view === 'home') {
      window.history.pushState({}, '', '/');
    } else if (view === 'upload') {
      window.history.pushState({}, '', '/?view=upload');
    } else if (view === 'results' && currentJob && currentJob.job_id) {
      window.history.pushState({}, '', `/?job_id=${encodeURIComponent(currentJob.job_id)}`);
    } else if (view === 'about') {
      window.history.pushState({}, '', '/?view=about');
    }
  }

  function handleGetStarted() {
    handleNavigate({ detail: { view: 'upload' } });
  }
</script>

<Navigation {currentView} hasJob={!!currentJob} on:navigate={handleNavigate} />

<main class="main-container">
  {#if currentView === 'home'}
    <LandingPage onGetStarted={handleGetStarted} />

  {:else if currentView === 'upload'}
    <PageHeader />
    <div class="content-wrapper">
      <UploadForm on:submitted={handleSubmitted} />
    </div>

  {:else if currentView === 'results'}
    <PageHeader showJobResults={true} />
    <div class="content-wrapper">
      <JobStatus jobId={currentJob.job_id} />
    </div>

  {:else if currentView === 'about'}
    <div class="section">
      <h1 class="heading-1">About SoloMicrobe</h1>
      
      <section>
        <h2 class="heading-2">Overview</h2>
        <p class="text-body">
          SoloMicrobe is an interactive web platform for designing oligonucleotide probes 
          that target microbial genes with high specificity. Users can provide microbial sequences or 
          pre-existing probe sets, and the tool screens each candidate against the host genome and 
          the co-residing microbes present in complex microbial communities, removing candidates with 
          up to two mismatches and a stretch of greater than or equal to 14 consecutive matchesto 
          off-target sequences.The final probe sets can be applied to detect microbial transcripts 
          within complex, host-associated tissues.
        </p>
      </section>

      <section>
        <h2 class="heading-2">How It Works</h2>
        <div class="info-box">
          <ol class="ordered-list">
            <li><strong>Sequence Input:</strong> Users paste either microbial gene sequences (for probe generation) or probe sequences (for direct alignment analysis) in FASTA format directly into the web interface</li>
            <li><strong>Probe Generation:</strong> The pipeline generates short oligonucleotide probes of specified length from the input genome</li>
            <li><strong>Host Genome Alignment:</strong> Each probe is aligned against the selected host genome (e.g. human) using high-performance alignment algorithms</li>
            <li><strong>Specificity Filtering:</strong> Probes that match the host genome within the specified mismatch tolerance have a potential to hybridize with host genomes.</li>
            <li><strong>Results & Visualization:</strong> For gene sequence input, results are displayed in an interactive genome browser showing probe positions and alignment details. 
              For probe sequence input, a detailed alignment table shows all matches against the host genome with 
              mismatch counts and filtering options</li>
          </ol>
        </div>
      </section>

      <section>
        <h2 class="heading-2">Key Features</h2>
        <div class="grid-auto">
          <div class="card">
            <h3 class="heading-3">Customizable Parameters</h3>
            <p class="text-small">Adjust probe length and mismatch tolerance to balance specificity and sensitivity</p>
          </div>
          <div class="card">
            <h3 class="heading-3">High Performance</h3>
            <p class="text-small">Optimized algorithms handle large genomes efficiently with parallel processing</p>
          </div>
          <div class="card">
            <h3 class="heading-3">Real-time Monitoring</h3>
            <p class="text-small">Track job progress in real-time with detailed status updates and estimated completion times</p>
          </div>
        </div>
      </section>

      <section>
        <h2 class="heading-2">Technical Details</h2>
        <div class="warning-box">
          <p class="info-box-text"><strong>Alignment Engine:</strong> razers3 for fast and accurate short-read alignment</p>
          <p class="info-box-text"><strong>Host Genomes:</strong> Pre-indexed reference genomes for hosts</p>
          <p class="info-box-text"><strong>Probe Length Range:</strong> 20-50 nucleotides (default: 30bp)</p>
          <p class="info-box-text"><strong>Mismatch Tolerance:</strong> 0-2 mismatches (default: 2)</p>
        </div>
      </section>

      <section>
        <h2 class="heading-2">Use Cases</h2>
        <ul class="unordered-list">
          <li>Design of highly specific probes for single-cell transcriptomics of host-associated microbiomes</li>
          <li>Spatial transcriptomics studies to localize microbial transcripts within host tissues</li>
          <li>Quality control and decontamination of metagenomic assemblies from clinical or environmental samples</li>
          <li>Decontamination of gut microbiome genome sequences</li>
          <li>Validation of microbial genome purity by excluding host-derived sequences</li>
          <li>Pre-processing and probe selection for comparative genomics and functional studies</li>
        </ul>
      </section>

      <section class="info-box-blue">
        <p class="blue-box-text"><strong>Citation:</strong> If you use this tool in your research, please cite</p>
      </section>

      <section>
        <h2 class="heading-2">Getting Started</h2>
        <p class="text-body" style="margin-bottom: 1rem;">
          Ready to design high-specificity probes for your microbiome study? 
          Go to the <strong>Design Probes</strong> tab to upload your 
          microbial gene sequences and set your analysis parameters. Probe filtering
          against the host genome will be completed within minutes.
        </p>
        <Button
          variant="gradient"
          size="md"
          onClick={() => handleNavigate({ detail: { view: 'upload' } })}
        >
          Start Analysis →
        </Button>
      </section>
    </div>
  {/if}
  
  <Bottom />
</main>