<script>
  import { onMount, onDestroy } from 'svelte';
  export let onGetStarted;

  let visible = false;

  onMount(() => {
    document.body.classList.add('landing-active');
    setTimeout(() => { visible = true; }, 100);
  });
  onDestroy(() => {
    document.body.classList.remove('landing-active');
  });
</script>

<svelte:head>
  <title>SHARP-FISH — Within-Organism & Microbial FISH Probe Design</title>
  <meta name="description" content="Design highly specific oligonucleotide FISH probes within an organism's transcriptome or for a microbial gene." />
</svelte:head>

<div class="sp-page" class:sp-visible={visible}>

  <!-- Hero -->
  <section class="sp-hero">
    <div class="sp-container sp-center">
      <h1 class="sp-hero-title">SHARP-FISH</h1>
      <p class="sp-hero-sub">
        Design highly specific oligonucleotide FISH probes within an organism's transcriptome or for a microbial gene
      </p>
    </div>
  </section>

  <!-- Tools -->
  <section id="tools" class="sp-section">
    <div class="sp-container">
      <div class="sp-section-head">
        <h2 class="sp-section-title">Probe Design Modes</h2>
        <p class="sp-section-sub">Choose where your target lives — within a single organism's transcriptome or in a microbial genome.</p>
      </div>

      <div class="sp-tool-grid">
        <button type="button" class="sp-tool-card sp-tool-host" on:click={() => onGetStarted('host')}>
          <span class="sp-tag sp-tag-host">Organism</span>
          <h3 class="sp-tool-title">Within-Organism Probe Design</h3>
          <p class="sp-tool-desc">
            Target a single gene's transcript and stay specific within the organism's transcriptome (human, mouse, zebrafish, fly, worm, frog, …).
          </p>
        </button>

        <button type="button" class="sp-tool-card sp-tool-microbe" on:click={() => onGetStarted('microbe')}>
          <span class="sp-tag sp-tag-microbe">Microbe</span>
          <h3 class="sp-tool-title">Microbial Probe Design</h3>
          <p class="sp-tool-desc">
            Target a microbial gene with optional cross-checks against host and co-residing microbiome
            catalogs (gut, oral, skin, vaginal, mouse-gut).
          </p>
        </button>
      </div>
    </div>
  </section>

</div>

<style>
  /* ─────────────────────────  Page shell  ─────────────────────────── */
  .sp-page {
    width: 100%;
    opacity: 0;
    transition: opacity 0.6s ease;
    background: transparent;
  }
  .sp-page.sp-visible { opacity: 1; }

  .sp-container {
    box-sizing: border-box;
    width: 100%;
    max-width: 1120px;
    margin: 0 auto;
    padding: 0 24px;
  }
  .sp-center { text-align: center; }

  /* ─────────────────────────  Hero  ────────────────────────────────── */
  .sp-hero {
    box-sizing: border-box;
    width: 100%;
    position: relative;
    padding: 80px 24px 96px;
    background: transparent;
  }
  .sp-hero-title {
    margin: 0 0 14px 0;
    font-size: 56px;
    font-weight: 800;
    line-height: 1.1;
    letter-spacing: -0.02em;
    background: linear-gradient(90deg, #2563eb 0%, #7c3aed 100%);
    -webkit-background-clip: text;
    background-clip: text;
    color: transparent;
  }
  .sp-hero-sub {
    margin: 0 auto;
    max-width: 720px;
    font-size: 19px;
    color: #475569;
    line-height: 1.55;
  }

  /* ─────────────────────────  Sections  ────────────────────────────── */
  .sp-section { padding: 72px 0 88px; }

  #tools {
    box-sizing: border-box;
    width: 100%;
    background: transparent;
  }

  .sp-section-head { text-align: center; margin-bottom: 40px; }
  .sp-section-title {
    margin: 0 0 10px 0;
    font-size: 30px;
    font-weight: 800;
    color: #0f172a;
    letter-spacing: -0.01em;
  }
  .sp-section-sub {
    margin: 0;
    color: #64748b;
    font-size: 16px;
    line-height: 1.55;
  }

  .sp-tool-grid {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(min(100%, 320px), 1fr));
    gap: 24px;
  }
  .sp-tool-card {
    box-sizing: border-box;
    position: relative;
    text-align: center;
    background: #ffffff;
    border: 1px solid #e2e8f0;
    border-radius: 16px;
    padding: 40px 28px 34px;
    cursor: pointer;
    transition: transform 0.2s ease, box-shadow 0.2s ease, border-color 0.2s ease;
    font: inherit;
    color: inherit;
    min-width: 0;
    overflow: hidden;
    word-wrap: break-word;
    overflow-wrap: anywhere;
    box-shadow: 0 1px 2px rgba(15, 23, 42, 0.04), 0 4px 14px rgba(15, 23, 42, 0.04);
  }
  /* Coloured accent bar across the top of each card */
  .sp-tool-card::before {
    content: '';
    position: absolute;
    top: 0;
    left: 0;
    right: 0;
    height: 5px;
  }
  .sp-tool-host::before {
    background: linear-gradient(90deg, #3b82f6 0%, #6366f1 50%, #7c3aed 100%);
  }
  .sp-tool-microbe::before {
    background: linear-gradient(90deg, #10b981 0%, #14b8a6 50%, #06b6d4 100%);
  }
  .sp-tool-card:hover, .sp-tool-card:focus-visible {
    transform: translateY(-4px);
    box-shadow: 0 14px 28px rgba(15, 23, 42, 0.1);
    outline: none;
  }

  .sp-tool-title {
    margin: 18px 0 12px 0;
    font-size: 22px;
    font-weight: 800;
    color: #0f172a;
  }
  .sp-tool-desc {
    margin: 0 auto;
    max-width: 42ch;
    font-size: 14.5px;
    color: #4b5563;
    line-height: 1.6;
  }
  .sp-tag {
    display: inline-block;
    padding: 5px 14px;
    border-radius: 999px;
    font-size: 11px;
    font-weight: 700;
    letter-spacing: 0.08em;
    text-transform: uppercase;
  }
  .sp-tag-host { background: #e0e7ff; color: #4338ca; }
  .sp-tag-microbe { background: #ccfbf1; color: #115e59; }

  @media (max-width: 640px) {
    .sp-hero { padding: 56px 0 36px; }
    .sp-hero-title { font-size: 40px; }
    .sp-hero-sub { font-size: 17px; }
    .sp-section { padding: 44px 0; }
    .sp-section-title { font-size: 26px; }
  }
</style>
