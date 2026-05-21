<script>
  export let showJobResults = false;
</script>

<div style="display:flex; flex-direction:column; gap:1rem;">
  <div>
    <h1 style="margin:0; font-size:2rem; font-weight:700; color:#111827;">
      SHARP-FISH
    </h1>
    <div style="margin-top:0.4rem; font-size:1.125rem; font-weight:400; color:#4b5563;">
      A webtool for designing host and microbial FISH probes with validated target specificity
    </div>
  </div>

  <p style="line-height:1.6; color:#4b5563; margin:0;">
    SHARP-FISH is a web-based tool for designing fluorescence
    in situ hybridization probes that target a specific gene within
    a chosen host (human or mouse) transcriptome or microbial genome. Each probe is
    validated against the host transcriptome and co-residing microbial genomes to
    prevent cross-hybridization. SHARP-FISH is freely available under GNU AGPL.
  </p>

  {#if showJobResults}
    <h1 style="margin:0; font-size:2rem; font-weight:700; color:#111827; margin-top:1rem;">
      Job Results
    </h1>
    <p style="line-height:1.6; color:#4b5563; margin:0;">
      Monitor your probe design job progress and download results when complete.
    </p>
  {/if}
</div>