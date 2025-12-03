<script>
  export let showJobResults = false;
</script>

<div style="display:flex; flex-direction:column; gap:1rem;">
  <h1 style="margin:0; font-size:2rem; font-weight:700; color:#1f2937;">
    Microbial Probe Designer
  </h1>
  
  <p style="line-height:1.6; color:#4b5563; margin:0;">
    Spatial transcriptomics using fluorescence in situ hybridization 
    enables visualization of microbial gene expression directly within host 
    tissues, providing critical insights into host-microbe interactions. 
    However, a major technical challenge arises from sequence similarity 
    between microbial genes and host transcripts. When designing probes to 
    target microbial transcripts, even partial sequence homology can cause 
    probes to cross-hybridize with host RNA, producing false-positive signals 
    that confound biological interpretation. This tool systematically validates 
    microbial probe specificity against host transcriptomes through mismatch-tolerant 
    alignment, identifies probes with potential host cross-hybridization, and ranks safe 
    probes by thermodynamic quality—enabling researchers to design validated probes that minimize 
    host-derived artifacts in FISH experiments.
  </p>

  {#if showJobResults}
    <h1 style="margin:0; font-size:2rem; font-weight:700; color:#1f2937; margin-top:1rem;">
      Job Results
    </h1>
    <p style="line-height:1.6; color:#4b5563; margin:0;">
      Monitor your probe design job progress and download results when complete. 
      The pipeline identifies host-specific probes to filter contaminating sequences 
      from microbiome assemblies.
    </p>
  {/if}
</div>