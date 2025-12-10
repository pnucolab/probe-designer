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
  let probe_length = 30;
  let max_mismatches = 2;
  let uploading = false;
  let error = '';

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

    if (probe_length <= 0) {
      error = 'Probe length must be a positive number';
      return;
    }

    const form = new FormData();

    if (inputType === 'gene') {
      form.append('gene_sequence', sequenceToSubmit);
      console.log('Appending gene_sequence');
    } else {
      form.append('probe_sequence', sequenceToSubmit);
      console.log('Appending probe_sequence');
    }

    form.append('species', species);
    if (inputType === 'gene') {
      form.append('probe_length', String(probe_length));
    }

    form.append('max_mismatches', String(max_mismatches));

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
        <input type="radio" bind:group={inputType} value="gene"> Gene Sequence
      </label>
      <label class="radio-label">
        <input type="radio" bind:group={inputType} value="probe"> Probe Sequence
      </label>
    </div>
  </div>
  <textarea
    bind:value={pastedText}
    placeholder={placeholderText}
    class="textarea-input"
  ></textarea>
  <div class="config-box">
    <h3 class="heading-3" style="margin-bottom: 16px;">Configuration</h3>

    <div class="form-grid">
      <div>
        <label for="species" class="label-text">Host Organism / Target Microbiome</label>
        <select id="species" bind:value={species} class="select-input">
          <option value="human">Human Transcriptome</option>
          <option value="gut-microbe">Human Gut Microbiome</option>
          <option value="human-oral-microbiome">Human Oral Microbiome</option>
          <option value="human-skin-microbiome">Human Skin Microbiome</option>
          <option value="human-vaginal-microbiome">Human Vaginal Microbiome</option>
          <option value="mouse-gut-microbiome">Mouse Gut Microbiome</option>
          <option value="mouse">Mouse Transcriptome</option>

        </select>
      </div>
     {#if inputType === 'gene'}
        <div>
          <label for="probe_length" class="label-text">Probe Length (bp)</label>
          <input id="probe_length" type="number" bind:value={probe_length} class="number-input" />
        </div>
      {/if}
      <div>
        <label for="max_mismatches" class="label-text">Max Mismatches</label>
        <input id="max_mismatches" type="number" min="0" bind:value={max_mismatches} class="number-input" />
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
  input[type="number"] {
    appearance: textfield;
  }
  
  input[type="number"]::-webkit-inner-spin-button,
  input[type="number"]::-webkit-outer-spin-button {
    -webkit-appearance: none;
    margin: 0;
  }
</style>