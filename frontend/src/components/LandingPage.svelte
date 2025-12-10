<script>
  import Button from './Button.svelte';
  export let onGetStarted;
</script>

<div class="landing-container">
  
  <div class="content-section">
    <h1 class="main-title">
       SoloMicrobe
    </h1>
    <h3 class="subtitle">Design highly specific oligonucleotide probes for microbial genomes</h3>
    
    <div class="image-container">
      <img src="/genome-mixing.png" alt="Genome mixing illustration" class="hero-image" />
      <div class="image-overlay">
        <p class="overlay-text">
          Filter out probes that cross-hybridize with host genomes. Perfect for single-cell and spatial transcriptomics studies of host-associated microbiomes.
        </p>
      </div>
    </div>
    
    <div class="cta-container">
      <Button 
        variant="gradient" 
        size="lg"
        onClick={onGetStarted}
      >
        Get Started →
      </Button>
    </div>
  </div>

  <div class="features-grid">
    <div class="feature-card">
      <div class="feature-icon"></div>
      <h3 class="feature-title">High Specificity</h3>
      <p class="feature-description">
        Filter probes based on alignment to host genomes with customizable mismatch tolerance
      </p>
    </div>

    <div class="feature-card">
      <div class="feature-icon"></div>
      <h3 class="feature-title">Fast Processing</h3>
      <p class="feature-description">
        Optimized algorithms handle large genomes efficiently with parallel processing
      </p>
    </div>

    <div class="feature-card">
      <div class="feature-icon"></div>
      <h3 class="feature-title">Interactive Results</h3>
      <p class="feature-description">
        Explore alignments with genome browser visualization and download filtered probe sets
      </p>
    </div>
  </div>
</div>