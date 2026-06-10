<script>
  import { createEventDispatcher, onMount } from 'svelte';
  import Button from './Button.svelte';
  const dispatch = createEventDispatcher();

  export let mode = 'microbe';
  let inputType = 'gene';

  // Organism registry loaded from /config/organisms (backed by
  // config/organisms.yml). Renders the host dropdown + microbiome
  // catalog list. Falls back to a minimal set if the fetch fails.
  let organisms = { hosts: [] };
  let organismsLoaded = false;
  onMount(async () => {
    try {
      const res = await fetch('/api/organisms');
      if (res.ok) {
        organisms = await res.json();
        const def = organisms.hosts.find(h => h.default) || organisms.hosts[0];
        if (def) species = def.id;
      }
    } catch (err) {
      console.warn('Failed to load organism registry:', err);
    } finally {
      organismsLoaded = true;
    }
  });
  $: currentHost = organisms.hosts.find(h => h.id === species);
  $: hostMicrobiomes = currentHost ? currentHost.microbiomes : [];

  // In microbial mode only hosts that have at least one microbiome catalog
  // are useful — drop the others. In host mode all organisms are eligible.
  $: visibleHosts = mode === 'microbe'
    ? organisms.hosts.filter(h => h.microbiomes && h.microbiomes.length > 0)
    : organisms.hosts;

  // If the user switches to microbial mode while currently on a host that
  // has no microbiomes (e.g. zebrafish), bounce them to the first valid one.
  $: if (mode === 'microbe' && organismsLoaded && visibleHosts.length > 0
        && !visibleHosts.some(h => h.id === species)) {
    species = visibleHosts[0].id;
  }
  const defaultGeneSequence = `>transcript:ENSB:3PgobK0mHtDbpdO CDS=1-37686
ATGAAGGGTTCCGACGGCACCTCGCCGCGCACCACGGACGCGCCGATCGCGGTCGTCGGA
CTGTCCTGCCGCCTTCCCGGAGCACCCGACCCCGCCACGTTCCGGCAACTGCTCCTCGAC
GGCGCCGACGCCATTACGGAGGCCCCCGAGGGACTGTGGGGAATGGGCACGGACGCGGAC
CTCGGCCGCCGGGGCGGATTCCTGGACCGGGACCGGATCGACCACTTCGACGCCGCCTTC
TTCGGCATATCGCCACGCGAGGCCGCGGCCATGGACCCCCAGCAGCGCCTGACCCTGGAA
CTGACCTGGGAGGCCCTCGAAGACGCCGGAATCATTCCGGACCGGCTCCGCGACAGCCGT
ACCGGCGTGTACATCGGAGTGATCGCGGACGACTACGCCACCCTGATCCGCCGGGGCGGC
CCGGCGGCCATCGACCGGCACAGCTTCACCGGACTCCACCGCGGCATCATCGCCAACCGC
GTCTCCTACCACCTCGGCCTGCGCGGCCCCAGCCTCACCCTCGACGCCGGCCAGGCATCA
TCGCTGGCGGCCATCCACCTGGCCTGCGAAAGCATCCGCCGCGGCGAGACGTCCCTCGCG
TCGCTGGCGGCCATCCACCTGGCCTGCGAAAGCATCCGCCGCGGCGAGACGTCCCTCGCG
ATCGCCGGCGGTGTCCATCTCAACCTCGCCGTCGAAAGCGGCGTCAGCGCAGAGCGGTTC
GGCGGTCTGTCGCCGGACGGCGTCACCTACACCTTCGACGCCCGCGCCAATGGCTTCGTA
CGCGGCGAGGGCGGCGGCGCCGTCGTCCTCAAGCCCCTCGCCGACGCCCTCGCCGACGGG
GACGCCGTGTACTGCGTCATCCGCGGCAGCGCGCTCAACAACGACGGTGGCGGGGACCAC
CTCACCACGCCCCACCAGGCCGCCCAGGAGGACCTCCTGCGGCGCGCCTACCGGCAGGCC
GGAGTCGACCCCGCCCGGGTCCAGTACGTGGAACTCCACGGCACCGGAACGAAGGTCGGC
GACCCGATCGAGGCCGCCGCCCTGGGCGAGGTGCTGGGCGCGGCCCGGCGGCCCACGGAT
GCCCCGCTGCTGGTGGGGTCGGCCAAGACCAACGTGGGCCATCTGGAGGGCGCGGCCGGT
GTCGTCGGCTTCATCAAGACGGCCCTGGGTCTCAAGCACGGCGAACTCTTCCCCAGCCTC`;

  const defaultProbeSequence = `>probe_0|start=1|end=20|transcript:ENSB:3PgobK0mHtDbpdO
ATGAAGGGTTCCGACGGCAC
>probe_1|start=97|end=116|transcript:ENSB:3PgobK0mHtDbpdO
ACGTTCCGGCAACTGCTCCT
>probe_2|start=100|end=119|transcript:ENSB:3PgobK0mHtDbpdO
TTCCGGCAACTGCTCCTCGA
>probe_3|start=211|end=230|transcript:ENSB:3PgobK0mHtDbpdO
GACCGGATCGACCACTTCGA
>probe_4|start=212|end=231|transcript:ENSB:3PgobK0mHtDbpdO
ACCGGATCGACCACTTCGAC
>probe_5|start=223|end=242|transcript:ENSB:3PgobK0mHtDbpdO
CACTTCGACGCCGCCTTCTT
>probe_6|start=224|end=243|transcript:ENSB:3PgobK0mHtDbpdO
ACTTCGACGCCGCCTTCTTC
>probe_7|start=230|end=249|transcript:ENSB:3PgobK0mHtDbpdO
ACGCCGCCTTCTTCGGCATA
>probe_8|start=231|end=250|transcript:ENSB:3PgobK0mHtDbpdO
CGCCGCCTTCTTCGGCATAT
>probe_9|start=232|end=251|transcript:ENSB:3PgobK0mHtDbpdO
GCCGCCTTCTTCGGCATATC`;

  let pastedText = '';
  $: placeholderText = inputType === 'gene' ? defaultGeneSequence : defaultProbeSequence;
  let species = 'human';
  let selectedMicrobiome = '';
  let kmerLength = '';
  let alignMicrobiome = false;
  let alignHost = true;
  $: if (mode === 'host') {
    alignHost = true;
    alignMicrobiome = false;
    selectedMicrobiome = '';
  }
  let probe_length = '';
  let max_mismatches = '';
  let tmRange = '';
  let gcRange = '';
  let uploading = false;
  let error = ''; 
  $: hostLabel = currentHost ? `${currentHost.display_name} Transcriptome` : 'Organism Transcriptome';
  let effectiveSpecies;
  $: {
    if (selectedMicrobiome && alignMicrobiome) {
      effectiveSpecies = selectedMicrobiome;
    } else {
      effectiveSpecies = species;
    }
  }
  async function submit(e) {
    e.preventDefault();
    error = '';

    console.log('=== FORM SUBMISSION DEBUG ===');
    console.log('Input type:', inputType);
    console.log('Pasted text length:', pastedText.length);
    console.log('Pasted text (first 100 chars):', pastedText.substring(0, 100));

    // Use default value if user didn't provide input
    const sequenceToSubmit = pastedText.trim().length === 0 
      ? (inputType === 'gene' ? defaultGeneSequence : defaultProbeSequence)
      : pastedText;

    console.log('Sequence to submit (first 100 chars):', sequenceToSubmit.substring(0, 100));

    if (probe_length !== '' && (probe_length < 20 || probe_length > 50)) {
      error = 'Probe length must be between 20 and 50 bp';
      return;
    }
    if (tmRange !== '') {
      const parts = tmRange.split('-').map(s => parseFloat(s.trim()));
      if (parts.some(isNaN) || parts.length < 1 || parts.length > 2) {
        error = 'Enter a single Tm (e.g. 45) or a range (e.g. 42-47)';
        return;
      }
      if (parts.some(v => v < 10 || v > 100)) {
        error = 'Please enter a valid Tm value';
        return;
      }
      if (parts.length === 2 && parts[0] > parts[1]) {
        error = 'Tm minimum cannot be greater than maximum';
        return;
      }
    }
    if (gcRange !== '') {
      const parts = gcRange.split('-').map(s => parseFloat(s.trim()));
      if (parts.some(isNaN) || parts.length < 1 || parts.length > 2) {
        error = 'Enter a single GC% (e.g. 50) or a range (e.g. 40-80)';
        return;
      }
      if (parts.some(v => v < 0 || v > 100)) {
        error = 'GC% must be between 0 and 100';
        return;
      }
      if (parts.length === 2 && parts[0] > parts[1]) {
        error = 'GC minimum cannot be greater than maximum';
        return;
      }
    }
    if (!alignMicrobiome && !alignHost) {
      error = 'Please select at least one target (Host Transcriptome or Additional Microbiome)';
      return;
    }
    if (alignMicrobiome && !selectedMicrobiome) {
      if (alignHost) {
        error = 'Please select a microbiome type or uncheck Additional Microbiome';
      } else {
        error = 'Please select a microbiome type, or uncheck Additional Microbiome and select a host transcriptome';
      }
      return;
    }

    if (inputType === 'gene') {
      const effectiveProbeLength = probe_length || 36;
      const lines = sequenceToSubmit.split('\n');
      let currentHeader = '';
      let currentSeq = '';
      const shortSequences = [];
      let headerCount = 0;

      for (const line of lines) {
        const trimmed = line.trim();
        if (trimmed.startsWith('>')) {
          headerCount++;
          if (currentHeader && currentSeq.length > 0 && currentSeq.length < effectiveProbeLength) {
            shortSequences.push({ header: currentHeader, length: currentSeq.length });
          }
          currentHeader = trimmed.substring(1).split(/\s/)[0] || 'unnamed';
          currentSeq = '';
        } else {
          currentSeq += trimmed.replace(/\s/g, '');
        }
      }
      if (currentHeader && currentSeq.length > 0 && currentSeq.length < effectiveProbeLength) {
        shortSequences.push({ header: currentHeader, length: currentSeq.length });
      }

      if (headerCount > 1) {
        error = `Only a single FASTA sequence is allowed per job. Found ${headerCount} sequences. Please submit one sequence at a time.`;
        return;
      }

      if (shortSequences.length > 0) {
        const details = shortSequences.map(s => `"${s.header}" (${s.length} bp)`).join(', ');
        error = `${shortSequences.length === 1 ? 'Sequence' : 'Sequences'} shorter than probe length (${effectiveProbeLength} bp): ${details}. Please remove short sequences or reduce the probe length.`;
        return;
      }
    }

    if (inputType === 'probe') {
      const lines = sequenceToSubmit.split('\n').map(l => l.trim()).filter(l => l.length > 0);
      const headersWithoutSeq = [];
      const seqsWithoutHeader = [];
      let lastWasHeader = false;
      let lastHeader = '';

      for (let i = 0; i < lines.length; i++) {
        const line = lines[i];
        if (line.startsWith('>')) {
          if (lastWasHeader) {
            headersWithoutSeq.push(lastHeader);
          }
          lastHeader = line.split(/\s/)[0].substring(1) || `line ${i + 1}`;
          lastWasHeader = true;
        } else if (/^[ATCGNatcgn]+$/.test(line)) {
          if (!lastWasHeader) {
            seqsWithoutHeader.push(`line ${i + 1}`);
          }
          lastWasHeader = false;
        }
      }
      if (lastWasHeader) {
        headersWithoutSeq.push(lastHeader);
      }

      if (headersWithoutSeq.length > 0) {
        error = `FASTA header without sequence: "${headersWithoutSeq[0]}". Every header must be followed by a nucleotide sequence.`;
        return;
      }
      if (seqsWithoutHeader.length > 0) {
        error = 'Sequence found without a FASTA header. Every sequence must be preceded by a header line (starting with ">").';
        return;
      }
      if (lines.length === 0 || !lines.some(l => l.startsWith('>'))) {
        error = 'Invalid FASTA format: no FASTA headers found. Each probe must have a header line (starting with ">") followed by its sequence.';
        return;
      }
    }

    // Host-token routing applies only when the user explicitly chose host probe
    // design. In microbial mode the input is treated as microbial regardless of
    // header content — any host alignment downstream is counted as cross-
    // reactivity, not as a routing decision here. The "Mode" chip at the top of
    // the form is the user-facing signal.
    if (mode === 'host') {
      try {
        const hostFd = new FormData();
        hostFd.append('fasta', sequenceToSubmit);
        hostFd.append('species', effectiveSpecies);
        const hostResp = await fetch('/validate-host-token', { method: 'POST', body: hostFd });
        if (hostResp.ok) {
          const hostJson = await hostResp.json();
          if (hostJson.is_host && hostJson.should_block) {
            error = hostJson.message;
            return;
          }
        }
      } catch (err) {
        console.warn('Host-token check failed (continuing — backend will re-check):', err);
      }
    }

    const form = new FormData();

    if (inputType === 'gene') {
      form.append('gene_sequence', sequenceToSubmit);
      console.log('Appending gene_sequence');
    } else {
      form.append('probe_sequence', sequenceToSubmit);
      console.log('Appending probe_sequence');
    }
    form.append('species', effectiveSpecies);
    form.append('mode', mode);
    form.append('align_microbiome', alignMicrobiome ? 'true' : 'false');
    form.append('align_host', alignHost ? 'true' : 'false');
    if (inputType === 'gene') {
      form.append('probe_length', String(probe_length || 36));
    }
    form.append('kmer_length', String(kmerLength || 18));
    form.append('max_mismatches', String(max_mismatches === '' || max_mismatches == null ? 2 : max_mismatches));
    form.append('tm_range', tmRange || '42-47');
    form.append('gc_range', gcRange || '40-80');

    console.log('Form data being sent:');
    for (let [key, value] of form.entries()) {
      console.log(`  ${key}:`, typeof value === 'string' ? value.substring(0, 50) + '...' : value);
    }

    uploading = true;
    try {
      console.log('Sending POST request to /jobs...');
      const res = await fetch('/jobs', { method: 'POST', body: form });
      console.log('Response status:', res.status);
      
      const contentType = res.headers.get('content-type') || '';
      let data;

      if (contentType.includes('application/json')) {
        data = await res.json();
      } else {
        const text = await res.text();
        console.log('Response text:', text);
        try {
          data = JSON.parse(text);
        } catch {
          data = { detail: text };
        }
      }

      console.log('Response data:', data);

      if (!res.ok) {
        throw new Error(data.detail || `HTTP ${res.status} ${res.statusText}`);
      }

      console.log('✅ Job submitted successfully:', data);
      dispatch('submitted', data);
    } catch (err) {
      console.error('❌ Submission error:', err);
      error = String(err);
    } finally {
      uploading = false;
    }
  }
</script>

<form on:submit|preventDefault={(e) => submit(e).catch(() => {})} class="form-container">
  <div class="form-section">
    <div class="label-text">Input Type</div>
    <div class="radio-group">
      <label class="radio-label">
        <input type="radio" bind:group={inputType} value="gene"> Gene Sequences (FASTA)
      </label>
      <label class="radio-label">
        <input type="radio" bind:group={inputType} value="probe"> Probe Sequences (FASTA)
      </label>
    </div>
  </div>
  {#if inputType === 'gene'}
    <div class="info-notice">
      <ul style="margin:0; padding-left:18px; list-style:disc;">
        {#if mode === 'host'}
          <li>For within-organism probe design, paste the <strong>transcript</strong> FASTA sequence (header should include a recognizable transcript ID or gene symbol — e.g. Ensembl <code>ENST…</code> / <code>ENSMUST…</code> / <code>ENSDART…</code> / <code>FBtr…</code>, NCBI <code>NM_…</code>, or a gene name like <code>EGFR</code>) and select the matching organism.</li>
        {:else}
          <li>For microbe probe design, paste the <strong>microbe gene</strong> FASTA sequence.</li>
        {/if}
        <li>The input should be a <strong>single</strong> FASTA sequence.</li>
      </ul>
    </div>
  {/if}
  <textarea
    bind:value={pastedText}
    placeholder={placeholderText}
    class="textarea-input"
  ></textarea>
    <div class="config-box">
      <h3 class="heading-3" style="margin-bottom: 16px;">Configuration</h3>
  
      <div class="form-grid">
        <div>
          <label for="species" class="label-text">{mode === 'host' ? 'Organism' : 'Host Organism'}</label>
          <select id="species" bind:value={species} on:change={() => { selectedMicrobiome = ''; alignMicrobiome = false; }} class="select-input">
            {#each visibleHosts as h}
              <option value={h.id}>{h.display_name}</option>
            {/each}
          </select>

          {#if mode !== 'host'}
          <div style="margin-top: 12px;">
            <div class="label-text">Select Target Transcriptome</div>
            <div style="border: 1px solid #d1d5db; border-radius: 6px; padding: 12px; margin-top: 4px; display: flex; flex-direction: column; gap: 6px;">

              <label class="radio-label">
                <input type="checkbox" bind:checked={alignHost}> {hostLabel}
                {#if currentHost?.reference_url}
                  <a href={currentHost.reference_url} target="_blank" rel="noopener" style="margin-left: 4px; font-size: 12px; color: #3b82f6;">(source)</a>
                {/if}
              </label>

              {#if hostMicrobiomes.length > 0}
                <label class="radio-label">
                  <input type="checkbox" bind:checked={alignMicrobiome} on:change={(e) => { if (!e.target.checked) selectedMicrobiome = ''; }}> Additional Microbiome
                </label>
                <div style="margin-left: 24px; display: flex; flex-direction: column; gap: 4px; opacity: {alignMicrobiome ? 1 : 0.5};">
                  {#each hostMicrobiomes as m}
                    <label class="radio-label">
                      <input type="radio" bind:group={selectedMicrobiome} value={m.id} disabled={!alignMicrobiome}> {m.display_name}
                      {#if m.source_url}
                        <a href={m.source_url} target="_blank" rel="noopener" style="margin-left: 4px; font-size: 12px; color: #3b82f6;">(source)</a>
                      {/if}
                    </label>
                  {/each}
                </div>
              {/if}
              {#if !alignHost && !alignMicrobiome}
                <div style="color: #dc2626; font-size: 13px; margin-top: 8px;">
                  ⚠️ Warning: At least one category must be selected.
                </div>
              {/if}
            </div>
          </div>
          {/if}
        </div>
        
        <div style="grid-column: 1 / -1;">
          <div class="label-text" style="color:#111827;">Probe Design</div>
          <div style="padding: 16px; border: 1px solid #e5e7eb; border-radius: 8px; background: white; display: grid; grid-template-columns: repeat({inputType === 'gene' ? 3 : 2}, 1fr); gap: 16px; margin-top: 4px;">
            {#if inputType === 'gene'}
              <div>
                <label for="probe_length" class="radio-label" style="display:block; margin-bottom:4px;">Probe Length (bp)</label>
                <input id="probe_length" type="number" bind:value={probe_length} placeholder="Range: 20-50 bp (default: 36)" min="20" max="50" class="number-input" />
              </div>
            {/if}
            <div>
              <label for="kmer_length" class="radio-label" style="display:block; margin-bottom:4px;">K-mer Length (bp)</label>
              <input id="kmer_length" type="text" bind:value={kmerLength} placeholder="Default: 18" class="number-input" />
            </div>
            <div>
              <label for="max_mismatches" class="radio-label" style="display:block; margin-bottom:4px;">Max Mismatches</label>
              <input id="max_mismatches" type="number" min="0" max="6" bind:value={max_mismatches} placeholder="Default: 2" class="number-input" />
            </div>
          </div>
        </div>
        <div style="grid-column: 1 / -1;">
          <div class="label-text" style="color:#111827;">Thermodynamic Constraints</div>
          <div style="padding: 16px; border: 1px solid #e5e7eb; border-radius: 8px; background: white; display: grid; grid-template-columns: repeat(2, 1fr); gap: 16px; margin-top: 4px;">
            <div>
              <label for="gc_range" class="radio-label" style="display:block; margin-bottom:4px;">GC Range (%)</label>
              <input id="gc_range" type="text" bind:value={gcRange} placeholder="e.g. 40-80 or 50 (default: 40-80)" class="number-input" />
            </div>
            <div>
              <label for="tm_range" class="radio-label" style="display:block; margin-bottom:4px;">Tm Range (°C)</label>
              <input id="tm_range" type="text" bind:value={tmRange} placeholder="e.g. 42-47 or 45 (default: 42-47)" class="number-input" />
            </div>
          </div>
        </div>
      </div>
    </div>

  <div class="form-actions">
    <Button 
      type="submit" 
      variant="gradient" 
      size="md"
      disabled={uploading}
      fullWidth={false}
    >
      {uploading ? '⏳ Processing...' : 'Submit'}
    </Button>

    {#if error}
      <div class="error-box">
        <strong>❌ Error:</strong> {error}
      </div>
    {/if}
    
    {#if uploading}
      <div class="loading-box">
        ⏳ Submitting job to pipeline...
      </div>
    {/if}
  </div>
</form>

<style>
  .info-notice {
    width: calc(100% - 20px);
    background: #f0f7ff;
    border-left: 4px solid #3b82f6;
    border-radius: 6px;
    padding: 12px 16px;
    font-size: 13px;
    color: #1e3a5f;
    line-height: 1.6;
    box-sizing: border-box;
  }

  input[type="number"] {
    appearance: textfield;
  }
  
  input[type="number"]::-webkit-inner-spin-button,
  input[type="number"]::-webkit-outer-spin-button {
    -webkit-appearance: none;
    margin: 0;
  }
</style>