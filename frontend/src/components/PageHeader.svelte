<script>
  export let showJobResults = false;
</script>

<div style="display:flex; flex-direction:column; gap:1rem;">
  <h1 style="margin:0; font-size:2rem; font-weight:700; color:#1f2937;">
    Microbial-probe designer
  </h1>
  
  <p style="line-height:1.6; color:#4b5563; margin:0;">
    Microbiome genome assemblies are often biased by contamination 
    from the host genome. During sequencing, host and microbial DNA 
    are extracted together, and if host reads aren't properly filtered 
    before assembly, fragments of host DNA can be incorporated into 
    microbial contigs. This leads to chimeric assemblies where the supposed 
    microbial genome includes host sequences, distorting gene content. 
    When designing probes to detect microbial transcripts in host tissues, 
    one major problem is sequence similarity between microbial genes and host 
    transcripts. Even partial homology can lead a probe to bind to host RNA, 
    producing false positives. So this tool is designed to filter probes that have
    off-target with the host genome.
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