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
  let inputType = null;
  let safeProbes = [];
  let filteredSafeProbes = [];
  let showSafeProbes = false;
  let loadingSafeProbes = false;
  let scoreFilter = 'all';
  let tmFilter = 'all';
  let gcFilter = 'all';
  let homopolymerFilter = 'all';
  let statusFilter = 'all';
  let safeProbesPage = 1;
  let safeProbesPageSize = 25;
  let safeProbesPageSizeOptions = [10, 25, 50, 100];

  $: paginatedProbes = filteredSafeProbes.slice(
    (safeProbesPage - 1) * safeProbesPageSize,
    safeProbesPage * safeProbesPageSize
  );
  $: safeProbesTotalPages = Math.ceil(filteredSafeProbes.length / safeProbesPageSize) || 1;
  $: isPrevDisabled = safeProbesPage <= 1;
  $: isNextDisabled = safeProbesPage >= safeProbesTotalPages;
  $: alignmentsPrevDisabled = currentPage <= 1;
  $: alignmentsNextDisabled = currentPage >= totalPages;
  $: safeProbesPageNumbers = getSafeProbesPageNumbers(safeProbesPage, safeProbesTotalPages);
  $: pageNumbers = getPageNumbers(currentPage, totalPages);

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
        gene_id: aln.gene_id,
        mismatches: aln.mismatches,
        position: aln.position,
        strand: aln.strand,
        species: aln.species
      });
    }
    
    return groups;
  }
  
  async function fetchAllProbeAlignments(probeId) {
    try {
      console.log('Fetching alignments for probe:', probeId);
      const res = await fetch(`/jobs/${jobId}/alignments?probe_id=${encodeURIComponent(probeId)}`);
      if (res.ok) {
        const data = await res.json();
        console.log('Found', data.alignments?.length || 0, 'alignments');
        return data.alignments || [];
      }
      return [];
    } catch (e) {
      console.error('Failed to fetch probe alignments:', e);
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
  
    let probeSequence = null;
    try {
      const baseProbeId = featureData.probe_id.split('|')[0].trim();
      
      console.log('Available files:', files.map(f => f.filename));
      
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
                  if (i + 1 < lines.length) {
                    probeSequence = lines[i + 1].trim();
                    console.log(`Found probe sequence in ${file.filename}:`, probeSequence);
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
   
        console.log('Fetching scoring data for safe probe:', featureData.probe_id);
        
        try {
          const scoresRes = await fetch(`/jobs/${jobId}/download/safe_probes_scores.txt`);
          if (scoresRes.ok) {
            const scoresText = await scoresRes.text();
            const lines = scoresText.split('\n');
            const baseProbeId = featureData.probe_id.split('|')[0].trim();
            
            for (const line of lines) {
              // Extract probe ID from the line and check for exact match
              const parts = line.trim().split(/\s+/);
              if (parts.length >= 11) {
                const lineProbeId = parts[1]; // probe_id is in column 2
                const lineProbeBase = lineProbeId.split('|')[0].trim();
                
                // Extract probe numbers for exact numeric comparison
                const baseProbeNum = baseProbeId.match(/probe_(\d+)/)?.[1];
                const lineProbeNum = lineProbeBase.match(/probe_(\d+)/)?.[1];
                
                // Exact match: compare full probe_id OR exact probe number
                if (lineProbeId === featureData.probe_id || baseProbeNum === lineProbeNum) {
                  const scoringData = {
                    rank: parts[0],
                    probe_id: parts[1],
                    sequence: parts[2],
                    score: parseFloat(parts[3]),
                    tm: parseFloat(parts[4]),
                    gc_content: parseFloat(parts[5]),
                    length: parseInt(parts[6]),
                    complexity: parseFloat(parts[7]),
                    sec_struct: parseFloat(parts[8]),
                    homopoly: parts[9],
                    status: parts[10],
                    rejected: parts[11]
                  };
                  expandedProbeAlignments = [{
                    type: 'scoring_data',
                    data: scoringData
                  }];
                  
                  console.log('Found scoring data:', scoringData);
                  break;
                }
              }
            }
          }
        } catch (e) {
          console.error('Failed to fetch scoring data:', e);
        }
        
        loadingProbeAlignments = false;
        return;
      }
      
      const allProbeAlignments = await fetchAllProbeAlignments(featureData.probe_id);
      
      if (allProbeAlignments.length > 0) {
  // For microbiome: check if only self-alignment
  const isMicrobiome = info?.species === 'gut-microbe' || 
                       info?.species === 'human-oral-microbiome' || 
                       info?.species === 'human-skin-microbiome' || 
                       info?.species === 'human-vaginal-microbiome' || 
                       info?.species === 'mouse-gut-microbiome';
  
  if (isMicrobiome) {
    const perfectMatches = allProbeAlignments.filter(aln => aln.mismatches === 0);
    if (perfectMatches.length === 1 && allProbeAlignments.length === 1) {
      expandedProbeAlignments = []; // Safe - only self-alignment
    } else {
      expandedProbeAlignments = allProbeAlignments.map((aln) => ({
        target_transcript: aln.target_transcript,
        gene_id: aln.gene_id,
        mismatches: aln.mismatches,
        position: aln.position,
        strand: aln.strand,
        species: aln.species,
        sequence: aln.sequence,
        gc_content: calculateGC(aln.sequence)
      }));
    }
  } else {
    expandedProbeAlignments = allProbeAlignments.map((aln) => ({
      target_transcript: aln.target_transcript,
      gene_id: aln.gene_id,
      mismatches: aln.mismatches,
      position: aln.position,
      strand: aln.strand,
      species: aln.species,
      sequence: aln.sequence,
      gc_content: calculateGC(aln.sequence)
    }));
  }
}
    } catch (error) {
      console.error('Error in handleFeatureClick:', error);
    } finally {
      loadingProbeAlignments = false;
    }
  }

  async function fetchSafeProbes() {
    loadingSafeProbes = true;
    try {
      const res = await fetch(`/jobs/${jobId}/download/safe_probes_scores.txt`);
      if (res.ok) {
        const text = await res.text();
        const lines = text.split('\n');
        
        let dataStart = 0;
        for (let i = 0; i < lines.length; i++) {
          if (lines[i].trim().startsWith('-'.repeat(10))) {
            dataStart = i + 1;
            break;
          }
        }
       
        const probes = [];
        for (let i = dataStart; i < lines.length; i++) {
          const line = lines[i].trim();
          if (!line || line.startsWith('=')) continue;
          
          
          const parts = line.split(/\s+/);
          
          console.log('Parsing line:', line);
          console.log('Parts:', parts);
          
          if (parts.length >= 11) {
            const statusParts = parts.slice(10);
            const status = statusParts.join(' ');
            const probe = {
              rank: parseInt(parts[0]),
              probe_id: parts[1],
              sequence: parts[2],
              score: parseFloat(parts[3]),
              tm: parseFloat(parts[4]),
              gc_content: parseFloat(parts[5]),
              length: parseInt(parts[6]),
              complexity: parseFloat(parts[7]),
              sec_struct: parseFloat(parts[8]),
              homopolymer: parts[9],
              status: status
            };
            
            probes.push(probe);
          }
        }
        
        console.log('Total probes parsed:', probes.length);
        console.log('Sample probe:', probes[0]);
        
        safeProbes = probes;
        filteredSafeProbes = [...safeProbes];
      }
    } catch (e) {
      console.error('Failed to fetch safe probes:', e);
    } finally {
      loadingSafeProbes = false;
    }
  }

  function filterSafeProbes() {
    filteredSafeProbes = safeProbes.filter(probe => {
  
      if (scoreFilter !== 'all') {
        if (scoreFilter === '80+' && probe.score < 80) return false;
        if (scoreFilter === '60-79' && (probe.score < 60 || probe.score >= 80)) return false;
        if (scoreFilter === '40-59' && (probe.score < 40 || probe.score >= 60)) return false;
        if (scoreFilter === '<40' && probe.score >= 40) return false;
      }
      
      if (tmFilter !== 'all') {
        if (tmFilter === '60+' && probe.tm < 60) return false;
        if (tmFilter === '55-59' && (probe.tm < 55 || probe.tm >= 60)) return false;
        if (tmFilter === '50-54' && (probe.tm < 50 || probe.tm >= 55)) return false;
        if (tmFilter === '<50' && probe.tm >= 50) return false;
      }
      
      if (gcFilter !== 'all') {
        if (gcFilter === '40-50' && (probe.gc_content < 40 || probe.gc_content > 50)) return false;
        if (gcFilter === '50-60' && (probe.gc_content < 50 || probe.gc_content > 60)) return false;
        if (gcFilter === '60-70' && (probe.gc_content < 60 || probe.gc_content > 70)) return false;
        if (gcFilter === '70-80' && (probe.gc_content < 70 || probe.gc_content > 80)) return false;
      }
      
      if (homopolymerFilter !== 'all') {
        console.log('Homopolymer filter:', homopolymerFilter, 'Probe value:', probe.homopolymer);
        if (homopolymerFilter === 'no' && probe.homopolymer !== 'No') return false;
        if (homopolymerFilter === 'yes' && probe.homopolymer === 'No') return false;
      }
      
      if (statusFilter !== 'all' && probe.status.toLowerCase() !== statusFilter.toLowerCase()) {
        return false;
      }
      
      return true;
    });
    safeProbesPage = 1;
  }

  async function toggleSafeProbes() {
    showSafeProbes = !showSafeProbes;
    if (showSafeProbes && safeProbes.length === 0) {
      await fetchSafeProbes();
    }
  }
  function getSafeProbesPageData() {
    const startIdx = (safeProbesPage - 1) * safeProbesPageSize;
    const endIdx = startIdx + safeProbesPageSize;
    return filteredSafeProbes.slice(startIdx, endIdx);
  }

  function getSafeProbesTotalPages() {
    return Math.ceil(filteredSafeProbes.length / safeProbesPageSize) || 1;
  }

  function goToSafeProbePage(page) {
    const totalPages = getSafeProbesTotalPages();
    if (page >= 1 && page <= totalPages) {
      safeProbesPage = page;
    }
  }

  function nextSafeProbePage() {
    const totalPages = getSafeProbesTotalPages();
    if (safeProbesPage < totalPages) {
      safeProbesPage++;
    }
  }

  function prevSafeProbePage() {
    if (safeProbesPage > 1) {
      safeProbesPage--;
    }
  }

  function getSafeProbesPageNumbers(currentPage, totalPages) {
    const pages = [];
    const maxVisible = 5;
    
    if (totalPages <= maxVisible) {
      for (let i = 1; i <= totalPages; i++) {
        pages.push(i);
      }
      return pages;
    }
    
    pages.push(1);
    
    if (currentPage <= 3) {
      for (let i = 2; i <= 4; i++) {
        pages.push(i);
      }
      pages.push('...');
      pages.push(totalPages);
      
    } else if (currentPage >= totalPages - 2) {
      pages.push('...');
      for (let i = totalPages - 3; i <= totalPages; i++) {
        pages.push(i);
      }
      
    } else {
      pages.push('...');
      pages.push(currentPage - 1);
      pages.push(currentPage);
      pages.push(currentPage + 1);
      pages.push('...');
      pages.push(totalPages);
    }
    
    return pages;
  }

  function handlePageSizeChange(event) {
    safeProbesPageSize = parseInt(event.target.value);
    safeProbesPage = 1;
    filteredSafeProbes = [...filteredSafeProbes]; 
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
      if (res.status === 404) {
        if (intervalId) {
          clearInterval(intervalId);
          intervalId = null;
        }
        error = `Job ${jobId} not found. It may have expired or been deleted.`;
        return;
      }
      if (!res.ok) throw new Error('Failed to get status');
      const data = await res.json();
      status = data.status;
      info = data.info;
      console.log('info.species:', info?.species);
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
      if (intervalId) {
        clearInterval(intervalId);
        intervalId = null;
      }
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

  function getPageNumbers(currentPage, totalPages) {
    const pages = [];
    const maxVisible = 5;
    
    if (totalPages <= maxVisible) {
      for (let i = 1; i <= totalPages; i++) {
        pages.push(i);
      }
      return pages;
    } 
    
    pages.push(1);
    
    if (currentPage <= 3) {
      for (let i = 2; i <= 4; i++) {
        pages.push(i);
      }
      pages.push('...');
      pages.push(totalPages);
      
    } else if (currentPage >= totalPages - 2) {
      pages.push('...');
      for (let i = totalPages - 3; i <= totalPages; i++) {
        pages.push(i);
      }
      
    } else {
      pages.push('...');
      pages.push(currentPage - 1);
      pages.push(currentPage);
      pages.push(currentPage + 1);
      pages.push('...');
      pages.push(totalPages);
    }
    
    return pages;
  }
  
  function getStatusLabel(probe) {
    if (probe.status === 'safe') {
      return 'No off-target alignments';
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
    {#if info.stats.candidate_probes === 0}
      <div style="
        padding:16px;
        background:#fef2f2;
        border:1px solid #fca5a5;
        border-radius:6px;
        margin-bottom:16px;
      ">
        <div style="display:flex; align-items:start; gap:12px;">
          <svg width="24" height="24" viewBox="0 0 24 24" fill="none" style="flex-shrink:0; margin-top:2px;">
            <circle cx="12" cy="12" r="10" fill="#dc2626"/>
            <path d="M12 8v4M12 16h.01" stroke="white" stroke-width="2" stroke-linecap="round"/>
          </svg>
          <div>
            <div style="font-weight:600; color:#991b1b; font-size:15px; margin-bottom:4px;">
              None of the probes passed the GC content filter (40-60%). 
              This may occur if the input sequence has extreme GC content. 
              Consider adjusting the probe length or target a different region of the sequence.
            </div>
          </div>
        </div>
      </div>
    {/if}
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
          <th style="padding:10px; border:1px solid #d1d5db; text-align:center; font-weight:700;">Total Candidate Probes (40–80% GC content)</th>
          <th style="padding:10px; border:1px solid #d1d5db; text-align:center; font-weight:700;">Total Safe Probes Identified</th>
          <th style="padding:10px; border:1px solid #d1d5db; text-align:center; font-weight:700;">K-mer Safety Analysis Result (&gt;=14 consecutive matches)</th>
        </tr>
      </thead>

      <tbody>
        <tr style="background:white;">
          <td style="padding:12px; border:1px solid #e5e7eb; text-align:center;">
            {formatNumber(info.stats.candidate_probes)}
          </td>
          <td style="padding:12px; border:1px solid #e5e7eb; text-align:center; line-height:1.6;">
            <div>{formatNumber(info.stats.non_aligned_probes)}</div>
            {#if info.stats.non_aligned_probes > 0}
            {#each files as f}
              {#if f.filename.toLowerCase() === 'safe_probes_scores.txt'}
                <a
                  href={`/jobs/${jobId}/download/${f.filename}`}
                  class="download-link"
                >
                  <span>Download Probes</span>
                </a>
              {/if}
            {/each}
            {/if}
          </td>

          
          <td style="padding:12px; border:1px solid #e5e7eb; text-align:center; line-height:1.6;">
            
            {#each files as f}
              {#if f.filename.toLowerCase() === '14mer_matches_report.txt'}
                <a
                  href={`/jobs/${jobId}/download/${f.filename}`}
                  class="download-link"
                >
                  Download Report
                </a>
              {/if}
            {/each}
          </td>
        </tr>
      </tbody>
    </table>
  </div>
{/if}
  {#if info && info.stats && info.stats.non_aligned_probes > 0}
    <div style="margin-top:1.5rem;">
      <button
        on:click={toggleSafeProbes}
        class="collapsible-button"
      >
        <div style="display: flex; align-items: center; gap: 12px;">
          <div style="text-align: left;">
            <div style="font-weight: 600; font-size: 15px; font-family: inherit;">
              Safe Probe Scoring Details ({formatNumber(info.stats.non_aligned_probes)} probes)
            </div>
            <div style="font-size: 12px; color: #6b7280; font-weight: 400;">
              Click to view quality metrics and rankings
            </div>
          </div>
        </div>
        <span style="transform: rotate({showSafeProbes ? '180' : '0'}deg); transition: transform 0.3s;">▼</span>
      </button>

      {#if showSafeProbes}
        <div transition:slide={{ duration: 300, easing: quintOut }} style="margin-top:12px; border:1px solid #e5e7eb; border-radius:6px; overflow:hidden;">
        
          <div style="padding:16px; background:#f9fafb; border-bottom:1px solid #e5e7eb;">
            
            <div style="margin-top:12px; display:flex; justify-content:space-between; align-items:center;">
              <div style="font-size:13px; color:#6b7280;">
                Showing {((safeProbesPage - 1) * safeProbesPageSize) + 1}-{Math.min(safeProbesPage * safeProbesPageSize, filteredSafeProbes.length)} of {filteredSafeProbes.length} probes
              </div>
              <div style="display:flex; align-items:center; gap:8px;">
                <label for="page-size" style="font-size:13px; color:#374151;">Rows per page:</label>
                <select
                  id="page-size"
                  on:change={handlePageSizeChange}
                  value={safeProbesPageSize}
                  style="padding:4px 8px; border:1px solid #d1d5db; border-radius:4px; font-size:13px;"
                >
                  {#each safeProbesPageSizeOptions as option}
                    <option value={option}>{option}</option>
                  {/each}
                </select>
              </div>
            </div>
          </div>

          {#if loadingSafeProbes}
            <div style="padding:40px; text-align:center;">
              <div class="spinner" style="width:32px; height:32px; margin:0 auto 12px;"></div>
              <p style="color:#6b7280; margin:0;">Loading probe scores...</p>
            </div>
          {:else}
            <div style="overflow-x:auto;">
              <table style="width:100%; border-collapse:collapse; font-size:12px;">
                <thead>
                  <tr style="background:#f9fafb; border-bottom:2px solid #d1d5db;">
                    <th style="padding:10px; text-align:center;">Rank</th>
                    <th style="padding:10px; text-align:left;">Probe ID</th>
                    <th style="padding:10px; text-align:left;">Sequence</th>
                    <th style="padding:10px; text-align:center; position:relative;">
                      <div style="margin-bottom:4px;">Score</div>
                      <select
                        on:change={(e) => {scoreFilter = e.target.value; filterSafeProbes();}}
                        value={scoreFilter}
                        style="width:100%; padding:2px 4px; border:1px solid #d1d5db; border-radius:3px; font-size:11px; background:white;"
                      >
                        <option value="all">All</option>
                        <option value="80+">80+</option>
                        <option value="60-79">60-79</option>
                        <option value="40-59">40-59</option>
                        <option value="<40">&lt;40</option>
                      </select>
                    </th>
                    <th style="padding:10px; text-align:center; position:relative;">
                      <div style="margin-bottom:4px;">Tm (°C)</div>
                      <select
                        on:change={(e) => {tmFilter = e.target.value; filterSafeProbes();}}
                        value={tmFilter}
                        style="width:100%; padding:2px 4px; border:1px solid #d1d5db; border-radius:3px; font-size:11px; background:white;"
                      >
                        <option value="all">All</option>
                        <option value="60+">60+</option>
                        <option value="55-59">55-59</option>
                        <option value="50-54">50-54</option>
                        <option value="<50">&lt;50</option>
                      </select>
                    </th>
                    <th style="padding:10px; text-align:center; position:relative;">
                      <div style="margin-bottom:4px;">GC%</div>
                      <select
                        on:change={(e) => {gcFilter = e.target.value; filterSafeProbes();}}
                        value={gcFilter}
                        style="width:100%; padding:2px 4px; border:1px solid #d1d5db; border-radius:3px; font-size:11px; background:white;"
                      >
                        <option value="all">All</option>
                        <option value="40-50">40-50%</option>
                        <option value="50-60">50-60%</option>
                        <option value="60-70">60-70%</option>
                        <option value="70-80">70-80%</option>
                      </select>
                    </th>
                    <th style="padding:10px; text-align:center;">Length</th>
                    <th style="padding:10px; text-align:center;">Complexity</th>
                    <th style="padding:10px; text-align:center;">Sec. Struct</th>
                    <th style="padding:10px; text-align:center; position:relative;">
                      <div style="margin-bottom:4px;">Homopolymer</div>
                      <select
                        on:change={(e) => {homopolymerFilter = e.target.value; filterSafeProbes();}}
                        value={homopolymerFilter}
                        style="width:100%; padding:2px 4px; border:1px solid #d1d5db; border-radius:3px; font-size:11px; background:white;"
                      >
                        <option value="all">All</option>
                        <option value="yes">Yes</option>
                        <option value="no">No</option>
                      </select>
                    </th>
                    <th style="padding:10px; text-align:center; position:relative;">
                      <div style="margin-bottom:4px;">Status</div>
                      <select
                        on:change={(e) => {statusFilter = e.target.value; filterSafeProbes();}}
                        value={statusFilter}
                        style="width:100%; padding:2px 4px; border:1px solid #d1d5db; border-radius:3px; font-size:11px; background:white;"
                      >
                        <option value="all">All</option>
                        <option value="excellent">Excellent</option>
                        <option value="very good">Very Good</option>
                        <option value="good">Good</option>
                        <option value="acceptable">Acceptable</option>
                        <option value="marginal">Marginal</option>
                        <option value="rejected">Rejected</option>
                      </select>
                    </th>
                  </tr>
                </thead>
                <tbody>
                  {#each paginatedProbes as probe (probe.probe_id)}
                    <tr style="border-bottom:1px solid #e5e7eb;">
                      <td style="padding:10px; text-align:center; color:#6b7280;">{probe.rank}</td>
                      <td style="padding:10px; font-family:monospace; color:#1f2937;">{probe.probe_id}</td>
                      <td style="padding:10px; font-family:monospace; font-size:10px; word-break:break-all; max-width:200px;">{probe.sequence}</td>
                      <td style="padding:10px; text-align:center;">
                        <span style="
                          padding:4px 8px;
                          border-radius:4px;
                          font-weight:600;
                          background:{probe.score >= 80 ? '#d1fae5' : probe.score >= 60 ? '#fef3c7' : '#fee2e2'};
                          color:{probe.score >= 80 ? '#065f46' : probe.score >= 60 ? '#92400e' : '#991b1b'};
                        ">
                          {probe.score.toFixed(1)}
                        </span>
                      </td>
                      <td style="padding:10px; text-align:center; font-family:monospace;">{probe.tm.toFixed(1)}</td>
                      <td style="padding:10px; text-align:center; font-family:monospace;">{probe.gc_content.toFixed(1)}</td>
                      <td style="padding:10px; text-align:center;">{probe.length}</td>
                      <td style="padding:10px; text-align:center; font-family:monospace;">{probe.complexity.toFixed(2)}</td>
                      <td style="padding:10px; text-align:center; font-family:monospace;">{probe.sec_struct.toFixed(2)}</td>
                      <td style="padding:10px; text-align:center;">{probe.homopolymer}</td>
                      <td style="padding:10px; text-align:center;">
                        <span style="
                          padding:3px 8px;
                          border-radius:3px;
                          font-size:11px;
                          font-weight:600;
                          background:{probe.status.toLowerCase() === 'excellent' ? '#d1fae5' : 
                                    probe.status.toLowerCase() === 'very good' ? '#dbeafe' :
                                    probe.status.toLowerCase() === 'good' ? '#e0e7ff' :
                                    probe.status.toLowerCase() === 'acceptable' ? '#fef3c7' :
                                    probe.status.toLowerCase() === 'rejected' ? '#fee2e2' : '#f3f4f6'};
                          color:{probe.status.toLowerCase() === 'excellent' ? '#065f46' :
                                probe.status.toLowerCase() === 'very good' ? '#1e40af' :
                                probe.status.toLowerCase() === 'good' ? '#3730a3' :
                                probe.status.toLowerCase() === 'acceptable' ? '#92400e' :
                                probe.status.toLowerCase() === 'rejected' ? '#991b1b' : '#374151'};
                        ">
                          {probe.status}
                        </span>
                      </td>
                    </tr>
                  {/each}
                </tbody>
              </table>
            </div>
            
            <div style="padding:16px; display:flex; justify-content:center; align-items:center; gap:8px; border-top:1px solid #e5e7eb;">
              <button
                on:click={prevSafeProbePage}
                disabled={isPrevDisabled}
                style="
                  padding:6px 12px;
                  border:1px solid #d1d5db;
                  border-radius:4px;
                  background:{isPrevDisabled ? '#e5e7eb' : '#3b82f6'};
                  color:{isPrevDisabled ? '#9ca3af' : 'white'};
                  cursor:{isPrevDisabled ? 'not-allowed' : 'pointer'};
                  font-weight:500;
                  font-size:13px;
                  opacity:{isPrevDisabled ? 0.6 : 1};
                "
              >
                Previous
              </button>
              
              {#each safeProbesPageNumbers as pageNum}
                {#if pageNum === '...'}
                  <span style="padding:6px 8px; color:#9ca3af;">...</span>
                {:else}
                  <button
                    on:click={() => goToSafeProbePage(pageNum)}
                    style="
                      padding:6px 12px;
                      border:1px solid #d1d5db;
                      border-radius:4px;
                      background:{pageNum === safeProbesPage ? '#3b82f6' : 'white'};
                      color:{pageNum === safeProbesPage ? 'white' : '#374151'};
                      cursor:pointer;
                      font-weight:{pageNum === safeProbesPage ? 600 : 400};
                      font-size:13px;
                    "
                  >
                    {pageNum}
                  </button>
                {/if}
              {/each}
              
              <button
                on:click={nextSafeProbePage}
                disabled={isNextDisabled}
                style="
                  padding:6px 12px;
                  border:1px solid #d1d5db;
                  border-radius:4px;
                  background:{isNextDisabled ? '#e5e7eb' : '#3b82f6'};
                  color:{isNextDisabled ? '#9ca3af' : 'white'};
                  cursor:{isNextDisabled ? 'not-allowed' : 'pointer'};
                  font-weight:500;
                  font-size:13px;
                  opacity:{isNextDisabled ? 0.6 : 1};
                "
              >
                Next
              </button>
            </div>

          {/if}
        </div>
      {/if}
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
            {#if expandedProbeAlignments.length > 0 && expandedProbeAlignments[0].type === 'scoring_data'}
              <div style="overflow-x:auto;">
                <div style="margin-bottom:12px; font-size:14px; color:#374151; font-weight:600;">
                  Probe Quality Metrics
                </div>
                <table style="width:100%; border-collapse:collapse; font-size:12px;">
                  <thead>
                    <tr style="background:#f9fafb; border-bottom:2px solid #d1d5db;">
                      <th style="padding:8px; text-align:left;">Metric</th>
                      <th style="padding:8px; text-align:center;">Value</th>
                    </tr>
                  </thead>
                  <tbody>
                    <tr style="border-bottom:1px solid #e5e7eb;">
                      <td style="padding:8px; font-weight:600;">Score</td>
                      <td style="padding:8px; text-align:center;">{expandedProbeAlignments[0].data.score.toFixed(2)}</td>
                    </tr>
                    <tr style="border-bottom:1px solid #e5e7eb;">
                      <td style="padding:8px; font-weight:600;">Quality Status</td>
                      <td style="padding:8px; text-align:center;">
                        <span style="padding:3px 8px; border-radius:3px; background:#10b981; color:white; font-weight:600;">
                          {expandedProbeAlignments[0].data.status}
                        </span>
                      </td>
                    </tr>
                    <tr style="border-bottom:1px solid #e5e7eb;">
                      <td style="padding:8px; font-weight:600;">Tm (°C)</td>
                      <td style="padding:8px; text-align:center;">{expandedProbeAlignments[0].data.tm.toFixed(2)}</td>
                    </tr>
                    <tr style="border-bottom:1px solid #e5e7eb;">
                      <td style="padding:8px; font-weight:600;">GC Content (%)</td>
                      <td style="padding:8px; text-align:center;">{expandedProbeAlignments[0].data.gc_content.toFixed(1)}</td>
                    </tr>
                    <tr style="border-bottom:1px solid #e5e7eb;">
                      <td style="padding:8px; font-weight:600;">Length (bp)</td>
                      <td style="padding:8px; text-align:center;">{expandedProbeAlignments[0].data.length}</td>
                    </tr>
                    <tr style="border-bottom:1px solid #e5e7eb;">
                      <td style="padding:8px; font-weight:600;">Complexity</td>
                      <td style="padding:8px; text-align:center;">{expandedProbeAlignments[0].data.complexity.toFixed(3)}</td>
                    </tr>
                    <tr style="border-bottom:1px solid #e5e7eb;">
                      <td style="padding:8px; font-weight:600;">Secondary Structure</td>
                      <td style="padding:8px; text-align:center;">{expandedProbeAlignments[0].data.sec_struct.toFixed(3)}</td>
                    </tr>
                    <tr style="border-bottom:1px solid #e5e7eb;">
                      <td style="padding:8px; font-weight:600;">Homopolymer</td>
                      <td style="padding:8px; text-align:center;">{expandedProbeAlignments[0].data.homopoly}</td>
                    </tr>
                  </tbody>
                </table>
                
                <div style="margin-top:16px; padding:12px; background:#ecfdf5; border:1px solid #10b981; border-radius:6px; text-align:center; color:#065f46;">
                  <p style="margin:0; font-weight:600;">Safe Probe - No off-target alignments</p>
                </div>
              </div>
            {:else}
              
              <div style="padding:20px; background:#ecfdf5; border:1px solid #10b981; border-radius:6px; text-align:center; color:#065f46;">
                <p style="margin:0; font-weight:600;">No off-target alignments</p>
              </div>
            {/if}
            
          {:else if expandedProbeAlignments.length === 0}
            <div style="padding:20px; background:#ecfdf5; border:1px solid #10b981; border-radius:6px; text-align:center; color:#065f46;">
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
                    <th style="padding:8px; text-align:left;">Off-Target</th>
                    <th style="padding:8px; text-align:left;">Annotation</th>
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
                      {#if info?.species === 'human' || info?.species === 'mouse'}
                        <td style="padding:8px; color:#1f2937; font-weight:600;">
                          {aln.gene_id || '-'}
                        </td>
                      {:else}
                        <td style="padding:8px; font-family:monospace; color:#374151; font-size:12px;">
                          {aln.species || '-'}
                        </td>
                      {/if}
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
                {#if info?.species === 'human' || info?.species === 'mouse'}
                  <th style="padding:10px; text-align:left;">Gene ID</th>
                {:else}
                  <th style="padding:10px; text-align:left;">Species</th>
                {/if}
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
                    <td style="padding:8px; color:#1f2937; font-weight:600;">
                      {#if info?.species === 'human' || info?.species === 'mouse'}
                        {aln.gene_id || '-'}
                      {:else}
                        {aln.species || '-'}
                      {/if}
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
          <button
            on:click={prevPage}
            disabled={alignmentsPrevDisabled}
            style="
              padding:6px 12px;
              border:1px solid #d1d5db;
              border-radius:4px;
              background:{alignmentsPrevDisabled ? '#e5e7eb' : '#3b82f6'};
              color:{alignmentsPrevDisabled ? '#9ca3af' : 'white'};
              cursor:{alignmentsPrevDisabled ? 'not-allowed' : 'pointer'};
              font-weight:500;
              font-size:13px;
              opacity:{alignmentsPrevDisabled ? 0.6 : 1};
            "
          >
            Previous
          </button>
          
          {#each pageNumbers as pageNum}
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
          
          <button
            on:click={nextPage}
            disabled={alignmentsNextDisabled}
            style="
              padding:6px 12px;
              border:1px solid #d1d5db;
              border-radius:4px;
              background:{alignmentsNextDisabled ? '#e5e7eb' : '#3b82f6'};
              color:{alignmentsNextDisabled ? '#9ca3af' : 'white'};
              cursor:{alignmentsNextDisabled ? 'not-allowed' : 'pointer'};
              font-weight:500;
              font-size:13px;
              opacity:{alignmentsNextDisabled ? 0.6 : 1};
            "
          >
            Next
          </button>
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
  .download-link {
    display: inline-flex;
    align-items: center;
    gap: 4px;
    color: #1d4ed8;
    font-size: 13px;
    text-decoration: none;
    font-style: italic;
    margin-top: 3px;
    padding: 3px 8px;
    border-radius: 3px;
    background: #dbeafe;
    border: 1px solid #bfdbfe;
    transition: background 0.15s ease;
  }
  
  .download-link:hover,
  .download-link:focus {
    background: #bfdbfe;
  }
  .collapsible-button {
  width: 100%;
  padding: 12px 16px;
  background: #d1fae5;  
  border: 1px solid #10b981;  
  border-radius: 6px;
  cursor: pointer;
  display: flex;
  justify-content: space-between;
  align-items: center;
  font-weight: 600;
  font-size: 15px;
  color: #065f46; 
  transition: background 0.2s;
}

.collapsible-button:hover {
  background: #a7f3d0;  
}

  
</style>