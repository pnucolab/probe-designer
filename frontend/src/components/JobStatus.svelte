<script>
  import { onMount } from 'svelte';
  import { fade, slide } from 'svelte/transition';
  import { quintOut } from 'svelte/easing';
  import JBrowseViewer from './JBrowseViewer.svelte';
  import Button from './Button.svelte';
  export let jobId;
  let status = null;
  let info = null;
  let files = [];
  let alignments = [];
  let totalAlignments = 0;
  let currentPage = 1;
  let pageSize = 25;
  let totalPages = 0;
  let error = null;
  let mismatchFilter = null; 
  let groupedAlignments = [];
  let sequenceLength = 0;
  let referenceId = 'reference';
  let expandedProbe = null;
  let expandedProbeAlignments = [];
  let loadingProbeAlignments = false;
  let loadingProgress = { current: 0, total: 0, percent: 0 };
  let searchProbeId = '';
  let searchedProbe = null;
  let searchedAlignments = [];
  let isSearching = false;
  let searchError = null;
  let inputType = null;
  
  function formatBytes(bytes) {
    if (!bytes && bytes !== 0) return '-';
    if (bytes === 0) return '0 B';
    const k = 1024;
    const sizes = ['B', 'KB', 'MB', 'GB', 'TB'];
    const i = Math.floor(Math.log(bytes) / Math.log(k));
    return parseFloat((bytes / Math.pow(k, i)).toFixed(2)) + ' ' + sizes[i];
  }

  function formatNumber(num) {
    if (num === null || num === undefined) return '-';
    return num.toLocaleString();
  }
  
  function calculateGC(sequence) {
    if (!sequence) return '-';
    const cleanSeq = sequence.replace(/-/g, '').toUpperCase();
    if (cleanSeq.length === 0) return '-';
    const gcCount = (cleanSeq.match(/[GC]/g) || []).length;
    return ((gcCount / cleanSeq.length) * 100).toFixed(1) + '%';
  }

  function groupAlignmentsByProbe(alignments) {
    const groups = [];
    let currentGroup = null;
    
    for (const aln of alignments) {
      if (!currentGroup || currentGroup.probe_id !== aln.probe_id) {
        currentGroup = {
          probe_id: aln.probe_id,
          sequence: aln.sequence,
          gc_content: calculateGC(aln.sequence),
          alignments: []
        };
        groups.push(currentGroup);
      }
      
      currentGroup.alignments.push({
        target_transcript: aln.target_transcript,
        mismatches: aln.mismatches,
        position: aln.position,
        strand: aln.strand
      });
    }
    
    return groups;
  }
  
  async function fetchAllProbeAlignments(probeId) {
    try {
      console.log('=== Fetching all alignments for probe:', probeId);
      
      const originalProbeId = probeId.trim();
      const baseProbeId = probeId.split('|')[0].trim();
      
      let allMatchingAlignments = [];
      let currentPage = 1;
      let totalPages = 1;
      
      const BATCH_SIZE = 5;
      const PAGE_SIZE = 100;
      
      while (currentPage <= totalPages) {
        const batchPromises = [];
        const batchStart = currentPage;
        const batchEnd = Math.min(currentPage + BATCH_SIZE - 1, totalPages);
        
        for (let batchPage = batchStart; batchPage <= batchEnd; batchPage++) {
          batchPromises.push(
            fetch(`/jobs/${jobId}/alignments?page=${batchPage}&page_size=${PAGE_SIZE}`)
              .then(res => res.ok ? res.json() : null)
              .then(data => ({ page: batchPage, data }))
          );
        }
        
        const batchResults = await Promise.all(batchPromises);
        
        for (const result of batchResults) {
          if (!result || !result.data) continue;
          
          const { data } = result;
          const pageAlignments = data.alignments || [];
          if (totalPages === 1) totalPages = data.total_pages || 1;
          
          const matching = pageAlignments.filter(aln => {
            const alnBase = aln.probe_id.split('|')[0].trim();
            return aln.probe_id.trim() === originalProbeId || alnBase === baseProbeId;
          });
          
          allMatchingAlignments.push(...matching);
        }
        
        currentPage = batchEnd + 1;
        
        loadingProgress = { 
          current: batchEnd, 
          total: totalPages, 
          percent: Math.round((batchEnd / totalPages) * 100)
        };
      }
      
      return allMatchingAlignments;
      
    } catch (e) {
      console.error('Failed to fetch all probe alignments:', e);
      return [];
    }
  }

    async function handleFeatureClick(featureData) {
    console.log('Feature clicked:', featureData);
    
    if (expandedProbe?.probe_id === featureData.probe_id) {
      expandedProbe = null;
      expandedProbeAlignments = [];
      return;
    }
    
    // Fetch the probe sequence from candidate probes file
    let probeSequence = null;
    try {
      const baseProbeId = featureData.probe_id.split('|')[0].trim();
      
      console.log('Available files:', files.map(f => f.filename));
      
      // Find all .fa or .fasta files that might contain probes
      const probeFiles = files.filter(f => {
        const name = f.filename.toLowerCase();
        return (name.endsWith('.fa') || name.endsWith('.fasta')) && 
               (name.includes('probe') || name.includes('candidate'));
      });
      
      console.log('Probe files to try:', probeFiles.map(f => f.filename));
      
      for (const file of probeFiles) {
        try {
          console.log(`Trying to fetch: ${file.filename}`);
          const res = await fetch(`/jobs/${jobId}/download/${file.filename}`);
          if (res.ok) {
            const fastaText = await res.text();
            const lines = fastaText.split('\n');
            
            for (let i = 0; i < lines.length; i++) {
              if (lines[i].startsWith('>')) {
                const header = lines[i].substring(1);
                const headerBase = header.split('|')[0].trim();
                
                if (header.trim() === featureData.probe_id || headerBase === baseProbeId) {
                  // Get the sequence on the next line
                  if (i + 1 < lines.length) {
                    probeSequence = lines[i + 1].trim();
                    console.log(`✓ Found probe sequence in ${file.filename}:`, probeSequence);
                    break;
                  }
                }
              }
            }
            
            if (probeSequence) break;
          }
        } catch (e) {
          console.log(`Could not fetch from ${file.filename}:`, e.message);
        }
      }
      
      if (!probeSequence) {
        console.warn('Probe sequence not found in candidate files');
        console.log('Files checked:', probeFiles.map(f => f.filename));
      }
    } catch (e) {
      console.error('Failed to fetch probe sequence:', e);
    }
    
    expandedProbe = { ...featureData, sequence: probeSequence };
    loadingProbeAlignments = true;
    expandedProbeAlignments = [];
    
    try {
      if (featureData.status === 'safe') {
        loadingProbeAlignments = false;
        return;
      }
      
      const allProbeAlignments = await fetchAllProbeAlignments(featureData.probe_id);
      
      if (allProbeAlignments.length > 0) {
        expandedProbeAlignments = allProbeAlignments.map((aln) => ({
          target_transcript: aln.target_transcript,
          mismatches: aln.mismatches,
          position: aln.position,
          strand: aln.strand,
          sequence: aln.sequence,
          gc_content: calculateGC(aln.sequence)
        }));
      }
    } catch (error) {
      console.error('Error in handleFeatureClick:', error);
    } finally {
      loadingProbeAlignments = false;
    }
  }
  async function searchProbeById() {
    if (!searchProbeId.trim()) {
      searchError = 'Please enter a probe ID';
      return;
    }
    
    console.log('\n========================================');
    console.log('=== SEARCH INITIATED ===');
    console.log('Searching for probe ID:', searchProbeId);
    console.log('========================================\n');
    
    isSearching = true;
    searchError = null;
    searchedAlignments = [];
    searchedProbe = null;
    
    try {
      searchedProbe = {
        probe_id: searchProbeId.trim(),
        status: 'unknown' 
      };
      
      console.log('🔄 Fetching alignments...');
      const allProbeAlignments = await fetchAllProbeAlignments(searchProbeId);
      
      console.log('📊 Fetch complete. Alignments found:', allProbeAlignments.length);
      
      if (allProbeAlignments.length > 0) {
        searchedAlignments = allProbeAlignments.map((aln, idx) => {
          const processed = {
            target_transcript: aln.target_transcript,
            mismatches: aln.mismatches,
            position: aln.position,
            strand: aln.strand,
            sequence: aln.sequence,
            gc_content: calculateGC(aln.sequence)
          };
          console.log(`  ${idx + 1}.`, processed.target_transcript, '- mismatches:', processed.mismatches);
          return processed;
        });
        
        const minMismatches = Math.min(...searchedAlignments.map(a => a.mismatches));
        if (minMismatches <= 1) {
          searchedProbe.status = 'high_risk';
          searchedProbe.mismatches = minMismatches;
        } else {
          searchedProbe.status = 'medium_risk';
          searchedProbe.mismatches = minMismatches;
        }
        
        console.log('✓ searchedAlignments set to:', searchedAlignments.length, 'items');
      } else {
        searchedProbe.status = 'safe';
        console.log('✓ Safe probe - no alignments found');
      }
    } catch (error) {
      console.error('❌ ERROR in searchProbeById:', error);
      searchError = 'Failed to fetch alignments: ' + error.message;
      searchedProbe = null;
    } finally {
      isSearching = false;
      console.log('\n✓ Search complete');
      console.log('Final state:');
      console.log('  - searchedProbe:', searchedProbe ? searchedProbe.probe_id : 'null');
      console.log('  - searchedAlignments.length:', searchedAlignments.length);
      console.log('========================================\n');
    }
  }

  function clearSearch() {
    searchProbeId = '';
    searchedProbe = null;
    searchedAlignments = [];
    searchError = null;
    loadingProgress = { current: 0, total: 0, percent: 0 };
  }
  
  function handleSearchKeypress(event) {
    if (event.key === 'Enter') {
      searchProbeById();
    }
  }

  function formatSequenceWithHighlights(sequence) {
    if (!sequence) return '';
    return sequence.split('').map((char, idx) => {
      const isLowercase = char === char.toLowerCase() && char !== char.toUpperCase();
      const isDash = char === '-';
      return { char, isLowercase, isDash, idx };
    });
  }

  async function fetchAlignments(page = 1) {
    try {
      let url = `/jobs/${jobId}/alignments?page=${page}&page_size=${pageSize}`;
      if (mismatchFilter !== null) url += `&mismatch=${mismatchFilter}`;
      const alignmentsRes = await fetch(url);
      if (alignmentsRes.ok) {
        const alignData = await alignmentsRes.json();
        alignments = alignData.alignments || [];
        groupedAlignments = groupAlignmentsByProbe(alignments);
        totalAlignments = alignData.total_alignments || 0;
        currentPage = alignData.page || 1;
        totalPages = alignData.total_pages || 0;
      }
    } catch (e) {
      console.error('Failed to fetch alignments:', e);
    }
  }

  function onMismatchFilterChange(value) {
    mismatchFilter = value === 'all' ? null : parseInt(value);
    fetchAlignments(1);
  }

  async function fetchStatus() {
    try {
      const res = await fetch(`/jobs/${jobId}`);
      if (!res.ok) throw new Error('Failed to get status');
      const data = await res.json();
      status = data.status;
      info = data.info;
      inputType = data.input_type || info?.input_type;

      if (status === 'SUCCESS' || status === 'FAILURE') {
        if (intervalId) {
          clearInterval(intervalId);
          intervalId = null;
        }
      }

      if (status === 'SUCCESS' || status === 'FAILURE') {
        const filesRes = await fetch(`/jobs/${jobId}/files`);
        if (filesRes.ok) files = await filesRes.json().then(r => r.files || []);
        
        await fetchAlignments(1);
        
        try {
          const faiRes = await fetch(`/jobs/${jobId}/download/reference.fasta.fai`);
          if (faiRes.ok) {
            const faiText = await faiRes.text();
            const faiLine = faiText.split('\n')[0];
            const parts = faiLine.split('\t');
            if (parts.length >= 2) {
              referenceId = parts[0];
              sequenceLength = parseInt(parts[1]);
            }
          }
        } catch (e) {
          console.error('Failed to fetch sequence length:', e);
        }
      }
    } catch (e) {
      error = String(e);
    }
  }

  function goToPage(page) {
    if (page >= 1 && page <= totalPages) {
      fetchAlignments(page);
    }
  }

  function nextPage() {
    if (currentPage < totalPages) {
      goToPage(currentPage + 1);
    }
  }

  function prevPage() {
    if (currentPage > 1) {
      goToPage(currentPage - 1);
    }
  }

  function getPageNumbers() {
    const pages = [];
    const maxVisible = 5;
    
    if (totalPages <= maxVisible) {
      for (let i = 1; i <= totalPages; i++) {
        pages.push(i);
      }
    } else {
      if (currentPage <= 3) {
        for (let i = 1; i <= 4; i++) pages.push(i);
        pages.push('...');
        pages.push(totalPages);
      } else if (currentPage >= totalPages - 2) {
        pages.push(1);
        pages.push('...');
        for (let i = totalPages - 3; i <= totalPages; i++) pages.push(i);
      } else {
        pages.push(1);
        pages.push('...');
        for (let i = currentPage - 1; i <= currentPage + 1; i++) pages.push(i);
        pages.push('...');
        pages.push(totalPages);
      }
    }
    
    return pages;
  }
  
  function getStatusLabel(probe) {
    if (probe.status === 'safe') {
      return 'Safe (No alignment)';
    } else if (probe.status === 'high_risk') {
      return `High Risk (${probe.mismatches} mismatch${probe.mismatches !== 1 ? 'es' : ''})`;
    } else if (probe.status === 'medium_risk') {
      return `Medium Risk (${probe.mismatches} mismatches)`;
    } else {
      return 'Unknown';
    }
  }
  
  function getStatusColor(status) {
    if (status === 'safe') return '#10b981';
    if (status === 'high_risk') return '#dc2626';
    if (status === 'medium_risk') return '#ea580c';
    return '#6b7280';
  }

  let intervalId;
  onMount(() => {
    fetchStatus();
    intervalId = setInterval(fetchStatus, 3000);
    return () => {
      clearInterval(intervalId);
    };
  });
</script>

<div style="margin-top:1.5rem; padding:12px; border:1px solid #e5e7eb; border-radius:6px; overflow:auto; box-sizing:border-box;">
  
  {#if error}
    <div style="color: red">{error}</div>
  {/if}

  <table style="width:100%; border-collapse: collapse; margin-top: 0.5rem;">
    <thead>
      <tr>
        <th style="text-align:left; border-bottom:1px solid #ddd; padding:12px 16px">Job ID</th>
        <th style="text-align:left; border-bottom:1px solid #ddd; padding:12px 16px">Status</th>
        <th style="text-align:right; border-bottom:1px solid #ddd; padding:12px 16px">Execution time (s)</th>
        <th style="text-align:left; border-bottom:1px solid #ddd; padding:12px 16px">Completed at (UTC)</th>
      </tr>
    </thead>
    <tbody>
      <tr>
        <td style="padding:12px 16px; vertical-align:top">{jobId}</td>
        <td style="padding:12px 16px; vertical-align:top">
          <div style="display:flex; align-items:center; gap:8px;">
            {#if status === 'PENDING' || status === 'PROGRESS'}
              <div class="spinner"></div>
              <span style="font-weight:500; color:#f59e0b;">In Progress</span>
            {:else if status === 'SUCCESS'}
              <svg width="16" height="16" viewBox="0 0 16 16" fill="none" style="flex-shrink:0;">
                <circle cx="8" cy="8" r="7" fill="#10b981"/>
                <path d="M5 8l2 2 4-4" stroke="white" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"/>
              </svg>
              <span style="font-weight:500; color:#10b981;">Completed</span>
            {:else if status === 'FAILURE'}
              <svg width="16" height="16" viewBox="0 0 16 16" fill="none" style="flex-shrink:0;">
                <circle cx="8" cy="8" r="7" fill="#dc2626"/>
                <path d="M6 6l4 4M10 6l-4 4" stroke="white" stroke-width="2" stroke-linecap="round"/>
              </svg>
              <span style="font-weight:500; color:#dc2626;">Failed</span>
            {:else}
              <span style="font-weight:500; color:#6b7280;">{status}</span>
            {/if}
          </div>
        </td>
        <td style="padding:12px 16px; vertical-align:top; text-align:right">{info && info.execution_time ? info.execution_time.toFixed(2) : '-'}</td>
        <td style="padding:12px 16px; vertical-align:top">{info && info.completed_at ? info.completed_at.split('T')[1].split('.')[0] : '-'}</td>
      </tr>
    </tbody>
  </table>
  {#if info && info.stats}
  <div
    style="
      margin-top:1.5rem;
      padding:16px;
      background:#f9fafb;
      border:1px solid #e5e7eb;
      border-radius:6px;
    "
  >
    <h3 style="margin:0 0 12px 0; font-weight:600; font-size:16px;">Result Summary</h3>

    <table
      style="
        width:100%;
        border-collapse:collapse;
        font-family:'Times New Roman', serif;
        font-size:14px;
        color:#111;
        border:1px solid #d1d5db;
      "
    >
      <thead>
        <tr style="background:#f3f4f6;">
          <th style="padding:10px; border:1px solid #d1d5db; text-align:center; font-weight:700;">Total Candidate Probes (40–60% GC content)</th>
          <th style="padding:10px; border:1px solid #d1d5db; text-align:center; font-weight:700;">Safe Probes (non-aligned)</th>
          <th style="padding:10px; border:1px solid #d1d5db; text-align:center; font-weight:700;">Filtered Alignments (&lt;=2 mismatches)</th>
        </tr>
      </thead>

      <tbody>
        <tr style="background:white;">
          <td style="padding:12px; border:1px solid #e5e7eb; text-align:center;">
            {formatNumber(info.stats.candidate_probes)}
          </td>
          <td style="padding:12px; border:1px solid #e5e7eb; text-align:center; line-height:1.6;">
            <div>{formatNumber(info.stats.non_aligned_probes)}</div>
            {#each files as f}
              {#if f.filename.toLowerCase() === 'non_aligned_probes.fa'}
                <a
                  href={`/jobs/${jobId}/download/${f.filename}`}
                  style="
                    display:inline-flex;
                    align-items:center;
                    gap:4px;
                    color:#1d4ed8;
                    font-size:13px;
                    text-decoration:none;
                    font-style:italic;
                    margin-top:3px;
                    padding:3px 8px;
                    border-radius:3px;
                    background:#dbeafe;
                    border:1px solid #bfdbfe;
                    transition:background 0.15s ease;
                  "
                  on:focus={() => (event.target.style.background = '#bfdbfe')}
                  on:blur={() => (event.target.style.background = '#dbeafe')}
                  on:mouseover={(e) => {
                    e.target.style.background = '#bfdbfe';
                  }}
                  on:mouseout={(e) => {
                    e.target.style.background = '#dbeafe';
                  }}
                >
                  <span>Download Probes</span>
                </a>
              {/if}
            {/each}
          </td>

          
          <td style="padding:12px; border:1px solid #e5e7eb; text-align:center; line-height:1.6;">
            <div>{formatNumber(info.stats.filtered_alignments)}</div>
            {#each files as f}
              {#if f.filename.toLowerCase() === 'probe_alignments.sam'}
                <a
                  href={`/jobs/${jobId}/download/${f.filename}`}
                  style="
                    display:inline-flex;
                    align-items:center;
                    gap:4px;
                    color:#1d4ed8;
                    font-size:13px;
                    text-decoration:none;
                    font-style:italic;
                    margin-top:3px;
                    padding:3px 8px;
                    border-radius:3px;
                    background:#dbeafe;
                    border:1px solid #bfdbfe;
                    transition:background 0.15s ease;
                  "
                  on:focus={() => (event.target.style.background = '#bfdbfe')}
                  on:blur={() => (event.target.style.background = '#dbeafe')}
                  on:mouseover={(e) => {
                    e.target.style.background = '#bfdbfe';
                  }}
                  on:mouseout={(e) => {
                    e.target.style.background = '#dbeafe';
                  }}
                >
                  Download Alignments
                </a>
              {/if}
            {/each}
          </td>
        </tr>
      </tbody>
    </table>
  </div>
{/if}
  
  {#if status === 'SUCCESS' && sequenceLength > 0 && inputType === 'gene_sequence'}
    <div style="margin-top:1.5rem;">
      <h3 style="margin:0 0 12px 0; font-weight:600; font-size:16px;">Probe Alignment Browser</h3>
      
      <JBrowseViewer 
        {jobId} 
        {referenceId}
        {sequenceLength}
        onFeatureClick={handleFeatureClick}
      />
      
      {#if expandedProbe}
        <div 
          transition:slide={{ duration: 300, easing: quintOut }}
          style="margin-top:16px; padding:20px; background:white; border:2px solid {getStatusColor(expandedProbe.status)}; border-radius:6px; box-shadow: 0 4px 12px rgba(0,0,0,0.1);"
        >
          <div style="display:flex; justify-content:flex-end; margin-bottom:8px;">
            <Button 
              variant="ghost" 
              size="sm"
              onClick={() => {expandedProbe = null; expandedProbeAlignments = [];}}
              style="transition: all 0.2s ease;"
              onMouseEnter={(e) => {
                e.currentTarget.style.transform = 'rotate(90deg) scale(1.2)';
                e.currentTarget.style.color = '#dc2626';
              }}
              onMouseLeave={(e) => {
                e.currentTarget.style.transform = 'rotate(0deg) scale(1)';
                e.currentTarget.style.color = '';
              }}
            >
              ×
            </Button>
          </div>

          <div style="margin-bottom:16px;">
            <h4 style="margin:0 0 8px 0; font-size:16px; font-weight:600; color:{getStatusColor(expandedProbe.status)};">
              {expandedProbe.probe_id}
            </h4>
            <div style="display:flex; gap:16px; flex-wrap:wrap; font-size:13px; color:#6b7280;">
              <span><strong>Position:</strong> {expandedProbe.start.toLocaleString()} - {expandedProbe.end.toLocaleString()} bp</span>
              <span><strong>Length:</strong> {(expandedProbe.end - expandedProbe.start).toLocaleString()} bp</span>
              <span><strong>Status:</strong> {getStatusLabel(expandedProbe)}</span>
            </div>
            {#if expandedProbe.sequence}
              <div style="margin-top:8px;">
                <div style="font-size:13px; color:#374151; font-weight:600; margin-bottom:6px;">Sequence</div>
                <div style="font-family:monospace; font-size:12px; line-height:1.5; background:#f8fafc; padding:8px; border-radius:6px; word-break:break-all; color:#111;">
                  {#each formatSequenceWithHighlights(expandedProbe.sequence) as {char, isLowercase, isDash}}
                    <span class:lowercase={isLowercase} class:dash={isDash}>{char}</span>
                  {/each}
                </div>
              </div>
            {:else}
              <div style="margin-top:8px; color:#9ca3af; font-size:13px;">
                <strong>Sequence:</strong> N/A
              </div>
            {/if}
          </div>

          {#if loadingProbeAlignments}
            <div style="text-align:center; padding:32px; color:#6b7280;">
              <div style="width:40px; height:40px; border:3px solid #e5e7eb; border-top-color:#3b82f6; border-radius:50%; margin:0 auto 12px; animation: spin 1s linear infinite;"></div>
              <p style="margin:0; font-size:14px;">Loading alignments...</p>
              {#if loadingProgress.total > 0}
                <div style="margin-top:8px; font-size:12px;">
                  {loadingProgress.current} / {loadingProgress.total} ({loadingProgress.percent}%)
                </div>
              {/if}
            </div>
          
          {:else if expandedProbe.status === 'safe'}
            <div style="padding:20px; background:#ecfdf5; border:1px solid #10b981; border-radius:6px; text-align:center; color:#065f46;">
              <div style="font-size:28px; margin-bottom:8px;">✓</div>
              <p style="margin:0; font-weight:600;">No off-target alignments</p>
            </div>
            
          {:else if expandedProbeAlignments.length === 0}
            <div style="padding:20px; background:#ecfdf5; border:1px solid #10b981; border-radius:6px; text-align:center; color:#065f46;">
              <div style="font-size:28px; margin-bottom:8px;">✓</div>
              <p style="margin:0; font-weight:600;">No off-target alignments</p>
            </div>
            
          {:else}
            <div style="overflow-x:auto;">
              <div style="margin-bottom:12px; font-size:14px; color:#374151; font-weight:600;">
                {expandedProbeAlignments.length} alignment{expandedProbeAlignments.length !== 1 ? 's' : ''}
              </div>
              <table style="width:100%; border-collapse:collapse; font-size:12px;">
                <thead>
                  <tr style="background:#f9fafb; border-bottom:2px solid #d1d5db;">
                    <th style="padding:8px; text-align:left;">#</th>
                    <th style="padding:8px; text-align:left;">Target</th>
                    <th style="padding:8px; text-align:center;">Mismatches</th>
                    <th style="padding:8px; text-align:center;">Position</th>
                    <th style="padding:8px; text-align:center;">Strand</th>
                    <th style="padding:8px; text-align:center;">GC%</th>
                    <th style="padding:8px; text-align:left;">Sequence</th>
                  </tr>
                </thead>
                <tbody>
                  {#each expandedProbeAlignments as aln, idx}
                    <tr style="border-bottom:1px solid #e5e7eb;">
                      <td style="padding:8px; color:#6b7280;">{idx + 1}</td>
                      <td style="padding:8px; font-family:monospace; color:#1f2937;">{aln.target_transcript}</td>
                      <td style="padding:8px; text-align:center;">
                        <span style="padding:3px 6px; border-radius:3px; font-weight:600; background:{aln.mismatches <= 1 ? '#fee2e2' : '#fff7ed'}; color:{aln.mismatches <= 1 ? '#991b1b' : '#9a3412'};">
                          {aln.mismatches}
                        </span>
                      </td>
                      <td style="padding:8px; text-align:center; font-family:monospace;">{aln.position?.toLocaleString() || '-'}</td>
                      <td style="padding:8px; text-align:center;">
                        <span style="padding:3px 6px; border-radius:3px; background:#f3f4f6;">{aln.strand}</span>
                      </td>
                      <td style="padding:8px; text-align:center; font-family:monospace;">{aln.gc_content}</td>
                      <td style="padding:8px;">
                        {#if aln.sequence}
                          <div style="font-family:monospace; font-size:10px; line-height:1.5; word-break:break-all;">
                            {#each formatSequenceWithHighlights(aln.sequence) as {char, isLowercase, isDash}}
                              <span class:lowercase={isLowercase} class:dash={isDash}>{char}</span>
                            {/each}
                          </div>
                        {:else}
                          <span style="color:#9ca3af;">N/A</span>
                        {/if}
                      </td>
                    </tr>
                  {/each}
                </tbody>
              </table>
            </div>
          {/if}
        </div>
      {/if}
      
    </div>
  {/if}
  {#if status === 'SUCCESS' && inputType === 'probe_sequence'}
    <div style="margin-top:1.5rem;">
      <h3 style="margin:0 0 12px 0; font-weight:600; font-size:16px;">Probe Alignment Browser</h3>
      
      <div style="margin-bottom:12px; display:flex; align-items:center; gap:12px;">
        <label for="mismatch-filter" style="font-size:14px; color:#374151;">Filter by mismatches:</label>
        <select 
          id="mismatch-filter"
          on:change={(e) => onMismatchFilterChange(e.target.value)}
          style="padding:6px 12px; border:1px solid #d1d5db; border-radius:4px; font-size:14px;"
        >
          <option value="all">All</option>
          <option value="0">0 mismatches</option>
          <option value="1">1 mismatch</option>
          <option value="2">2 mismatches</option>
          <option value="3">3 mismatches</option>
        </select>
        <span style="margin-left:auto; font-size:14px; color:#6b7280;">
          {totalAlignments.toLocaleString()} total alignments
        </span>
      </div>

      {#if groupedAlignments.length > 0}
        <div style="overflow-x:auto; border:1px solid #e5e7eb; border-radius:6px;">
          <table style="width:100%; border-collapse:collapse; font-size:13px;">
            <thead>
              <tr style="background:#f9fafb; border-bottom:2px solid #d1d5db;">
                <th style="padding:10px; text-align:left;">Probe ID</th>
                <th style="padding:10px; text-align:center;">GC%</th>
                <th style="padding:10px; text-align:left;">Sequence</th>
                <th style="padding:10px; text-align:left;">Target</th>
                <th style="padding:10px; text-align:center;">Mismatches</th>
                <th style="padding:10px; text-align:center;">Position</th>
                <th style="padding:10px; text-align:center;">Strand</th>
              </tr>
            </thead>
            <tbody>
              {#each groupedAlignments as group}
                {#each group.alignments as aln, idx}
                  <tr style="border-bottom:1px solid #e5e7eb;">
                    {#if idx === 0}
                      <td style="padding:10px; font-family:monospace; color:#1f2937; vertical-align:top;" rowspan={group.alignments.length}>
                        {group.probe_id}
                      </td>
                      <td style="padding:10px; text-align:center; font-family:monospace; vertical-align:top;" rowspan={group.alignments.length}>
                        {group.gc_content}
                      </td>
                      <td style="padding:10px; vertical-align:top;" rowspan={group.alignments.length}>
                        <div style="font-family:monospace; font-size:11px; line-height:1.5; word-break:break-all;">
                          {#each formatSequenceWithHighlights(group.sequence) as {char, isLowercase, isDash}}
                            <span class:lowercase={isLowercase} class:dash={isDash}>{char}</span>
                          {/each}
                        </div>
                      </td>
                    {/if}
                    <td style="padding:10px; font-family:monospace; color:#374151;">
                      {aln.target_transcript}
                    </td>
                    <td style="padding:10px; text-align:center;">
                      <span style="padding:3px 8px; border-radius:3px; font-weight:600; background:{aln.mismatches <= 1 ? '#fee2e2' : '#fff7ed'}; color:{aln.mismatches <= 1 ? '#991b1b' : '#9a3412'};">
                        {aln.mismatches}
                      </span>
                    </td>
                    <td style="padding:10px; text-align:center; font-family:monospace; color:#6b7280;">
                      {aln.position?.toLocaleString() || '-'}
                    </td>
                    <td style="padding:10px; text-align:center;">
                      <span style="padding:3px 8px; background:#e5e7eb; border-radius:3px; font-family:monospace;">
                        {aln.strand}
                      </span>
                    </td>
                  </tr>
                {/each}
              {/each}
            </tbody>
          </table>
        </div>
        <div style="margin-top:16px; display:flex; justify-content:center; align-items:center; gap:8px;">
          <Button 
            variant="primary" 
            size="sm"
            onClick={prevPage}
            disabled={currentPage === 1}
          >
            Previous
          </Button>
          
          {#each getPageNumbers() as pageNum}
            {#if pageNum === '...'}
              <span style="padding:6px 8px; color:#9ca3af;">...</span>
            {:else}
              <button
                on:click={() => goToPage(pageNum)}
                style="padding:6px 12px; border:1px solid #d1d5db; border-radius:4px; background:{pageNum === currentPage ? '#3b82f6' : 'white'}; color:{pageNum === currentPage ? 'white' : '#374151'}; cursor:pointer; font-weight:{pageNum === currentPage ? 600 : 400};"
              >
                {pageNum}
              </button>
            {/if}
          {/each}
          
          <Button 
            variant="primary" 
            size="sm"
            onClick={nextPage}
            disabled={currentPage === totalPages}
          >
            Next
          </Button>
        </div>
      {:else}
        <div style="padding:40px; text-align:center; color:#6b7280; background:#f9fafb; border-radius:6px;">
          No alignments found for the selected filter
        </div>
      {/if}
    </div>
  {/if}

</div>

<style>
  @keyframes spin {
    from { transform: rotate(0deg); }
    to { transform: rotate(360deg); }
  }

  .spinner {
    width: 16px;
    height: 16px;
    border: 2px solid #e5e7eb;
    border-top-color: #3b82f6;
    border-radius: 50%;
    animation: spin 0.8s linear infinite;
  }

  .lowercase {
    background-color: rgba(239, 68, 68, 0.3);
    color: #991b1b;
    font-weight: 700;
    padding: 1px 2px;
    border-radius: 2px;
  }

  .dash {
    background-color: rgba(251, 191, 36, 0.3);
    color: #92400e;
    padding: 1px 2px;
    border-radius: 2px;
  }
</style>