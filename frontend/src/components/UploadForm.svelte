<script>
  import { createEventDispatcher } from 'svelte';
  import Button from './Button.svelte';
  const dispatch = createEventDispatcher();

  let inputType = 'gene';
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
  let probe_length = '';
  let max_mismatches = '';
  let tmRange = '';
  let uploading = false;
  let error = ''; 
  $: hostLabel = species === 'mouse' ? 'Mouse Transcriptome' : 'Human Transcriptome';
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

    const form = new FormData();

    if (inputType === 'gene') {
      form.append('gene_sequence', sequenceToSubmit);
      console.log('Appending gene_sequence');
    } else {
      form.append('probe_sequence', sequenceToSubmit);
      console.log('Appending probe_sequence');
    }
    form.append('species', effectiveSpecies);
    form.append('align_microbiome', alignMicrobiome ? 'true' : 'false');
    form.append('align_host', alignHost ? 'true' : 'false'); 
    if (inputType === 'gene') {
      form.append('probe_length', String(probe_length || 36));
    }
    form.append('kmer_length', String(kmerLength || 18));
    form.append('max_mismatches', String(max_mismatches ?? 2));
    form.append('tm_range', tmRange || '42-47');

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
        <input type="radio" bind:group={inputType} value="gene"> Transcript FASTA
      </label>
      <label class="radio-label">
        <input type="radio" bind:group={inputType} value="probe"> Probe FASTA
      </label>
    </div>
  </div>
  {#if inputType === 'gene'}
    <div class="info-notice">
      <div style="margin-bottom:6px;">
        <strong>⚠ Only DNA-form sequences accepted (A, T, C, G, N).</strong> RNA sequences containing U must be converted to T before submission.
      </div>
      <ul style="margin:0; padding-left:18px; list-style:disc;">
        <li>For host probe design, paste a <strong>transcript</strong> FASTA sequence.</li>
        <li>For microbe probe design, paste the <strong>microbe gene</strong> FASTA sequence.</li>
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
          <label for="species" class="label-text">Host Organism</label>
          <select id="species" bind:value={species} on:change={() => { selectedMicrobiome = ''; alignMicrobiome = false; }} class="select-input">
            <option value="human">Human</option>
            <option value="mouse">Mouse</option>
          </select>
          
          <div style="margin-top: 12px;">
            <div class="label-text">Select Target Transcriptome</div>
            <div style="border: 1px solid #d1d5db; border-radius: 6px; padding: 12px; margin-top: 4px; display: flex; flex-direction: column; gap: 6px;">
              
              <label class="radio-label">
                <input type="checkbox" bind:checked={alignHost}> {hostLabel}
                {#if species === 'human'}
                  <a href="https://www.gencodegenes.org/human/" target="_blank" rel="noopener" style="margin-left: 4px; font-size: 12px; color: #3b82f6;">(source)</a>
                {:else}
                  <a href="https://www.gencodegenes.org/mouse/" target="_blank" rel="noopener" style="margin-left: 4px; font-size: 12px; color: #3b82f6;">(source)</a>
                {/if}
              </label>

              <label class="radio-label">
                <input type="checkbox" bind:checked={alignMicrobiome} on:change={(e) => { if (!e.target.checked) selectedMicrobiome = ''; }}> Additional Microbiome
              </label>
              <div style="margin-left: 24px; display: flex; flex-direction: column; gap: 4px; opacity: {alignMicrobiome ? 1 : 0.5};">
                {#if species === 'human'}
                  <label class="radio-label">
                    <input type="radio" bind:group={selectedMicrobiome} value="gut-microbe" disabled={!alignMicrobiome}> Human Gut Microbiome
                    <a href="https://www.ebi.ac.uk/metagenomics/genome-catalogues/human-gut-v2-0-2" target="_blank" rel="noopener" style="margin-left: 4px; font-size: 12px; color: #3b82f6;">(source)</a>
                  </label>
                  <label class="radio-label">
                    <input type="radio" bind:group={selectedMicrobiome} value="human-oral-microbiome" disabled={!alignMicrobiome}> Human Oral Microbiome
                    <a href="https://www.ebi.ac.uk/metagenomics/genome-catalogues/human-oral-v1-0-1" target="_blank" rel="noopener" style="margin-left: 4px; font-size: 12px; color: #3b82f6;">(source)</a>
                  </label>
                  <label class="radio-label">
                    <input type="radio" bind:group={selectedMicrobiome} value="human-skin-microbiome" disabled={!alignMicrobiome}> Human Skin Microbiome
                    <a href="https://www.ebi.ac.uk/metagenomics/genome-catalogues/human-skin-v1-0" target="_blank" rel="noopener" style="margin-left: 4px; font-size: 12px; color: #3b82f6;">(source)</a>
                  </label>
                  <label class="radio-label">
                    <input type="radio" bind:group={selectedMicrobiome} value="human-vaginal-microbiome" disabled={!alignMicrobiome}> Human Vaginal Microbiome
                    <a href="https://www.ebi.ac.uk/metagenomics/genome-catalogues/human-vaginal-v1-0" target="_blank" rel="noopener" style="margin-left: 4px; font-size: 12px; color: #3b82f6;">(source)</a>
                  </label>
                {:else}
                  <label class="radio-label">
                    <input type="radio" bind:group={selectedMicrobiome} value="mouse-gut-microbiome" disabled={!alignMicrobiome}> Mouse Gut Microbiome
                    <a href="https://www.ebi.ac.uk/metagenomics/genome-catalogues/mouse-gut-v1-0" target="_blank" rel="noopener" style="margin-left: 4px; font-size: 12px; color: #3b82f6;">(source)</a>
                  </label>
                {/if}
              </div>
              {#if !alignHost && !alignMicrobiome}
                <div style="color: #dc2626; font-size: 13px; margin-top: 8px;">
                  ⚠️ Warning: At least one category must be selected.
                </div>
              {/if}
            </div>
          </div>
        </div>
        
       {#if inputType === 'gene'}
          <div>
            <label for="probe_length" class="label-text">Probe Length (bp)</label>
            <input id="probe_length" type="number" bind:value={probe_length} placeholder="Range: 20-50 bp (default: 36)" min="20" max="50" class="number-input" />        </div>
        {/if}
        <div>
          <label for="kmer_length" class="label-text">K-mer Length (bp)</label>
          <input id="kmer_length" type="text" bind:value={kmerLength} placeholder="Default: 18" class="number-input" />      
        </div>
        <div>
          <label for="max_mismatches" class="label-text">Max Mismatches</label>
          <input id="max_mismatches" type="number" min="0" max="2" bind:value={max_mismatches} placeholder="Default: 2" class="number-input" />
        </div>
        <div>
          <label for="tm_range" class="label-text">Tm Range (°C)</label>
          <input id="tm_range" type="text" bind:value={tmRange} placeholder="e.g. 42-47 or 45 (default: 42-47)" class="number-input" />
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