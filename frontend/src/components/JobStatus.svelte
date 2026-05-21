<script>
  import { onMount } from 'svelte';
  import { fade, slide } from 'svelte/transition';
  import { quintOut } from 'svelte/easing';
  import JBrowseViewer from './JBrowseViewer.svelte';
  import Button from './Button.svelte';
  export let jobId;
  let status = null;
  let info = null;
  let submittedAt = null;
  let files = [];
  let alignments = [];
  let loadingAlignments = true;
  let totalAlignments = 0;
  const ALIGNMENT_PREVIEW_SIZE = 100;
  const ALIGNMENTS_PAGE_SIZE = 10;
  let alignmentsPage = 1;
  let error = null;
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
  let safeProbesFetched = false;
  let tmFilter = '';
  let gcFilter = '';
  let safeProbesPage = 1;
  let safeProbesPageSize = 10;
  let safeProbesPageSizeOptions = [10, 25, 50, 100];

  let kmerReportOpening = false;

  let probeRiskSummary = null;
  let probeRiskSummaryLoading = false;

  async function loadProbeRiskSummary() {
    if (probeRiskSummary || probeRiskSummaryLoading) return;
    probeRiskSummaryLoading = true;
    try {
      const res = await fetch(`/jobs/${jobId}/download/probes.gff3`);
      if (!res.ok) throw new Error(`HTTP ${res.status}`);
      const text = await res.text();
      const counts = {
        safe: 0,
        high_risk: 0,
        medium_risk_kmer: 0,
        medium_risk_mismatch: 0,
        no_alignment: 0,
        total: 0,
      };
      for (const line of text.split('\n')) {
        if (!line || line.startsWith('#')) continue;
        const cols = line.split('\t');
        if (cols.length < 9) continue;
        const attrs = {};
        for (const kv of cols[8].split(';')) {
          const i = kv.indexOf('=');
          if (i > 0) attrs[kv.slice(0, i)] = kv.slice(i + 1);
        }
        const risk = attrs.risk_level;
        const desc = (attrs.description || '').toLowerCase();
        counts.total++;
        if (risk === 'safe') counts.safe++;
        else if (risk === 'high_risk') counts.high_risk++;
        else if (risk === 'no_alignment') counts.no_alignment++;
        else if (risk === 'medium_risk') {
          if (desc.includes('k-mer')) counts.medium_risk_kmer++;
          else counts.medium_risk_mismatch++;
        }
      }
      probeRiskSummary = counts;
    } catch (err) {
      console.warn('Failed to load probe risk summary:', err);
    } finally {
      probeRiskSummaryLoading = false;
    }
  }

  $: if (status === 'SUCCESS' && info?.stats?.candidate_probes > 0 && !probeRiskSummary && !probeRiskSummaryLoading) {
    loadProbeRiskSummary();
  }

  async function openKmerReportInTab() {
    if (kmerReportOpening) return;
    kmerReportOpening = true;
    try {
      const res = await fetch(`/jobs/${jobId}/download/kmer_matches_report.txt`);
      if (!res.ok) throw new Error(`HTTP ${res.status}`);
      const blob = new Blob([await res.text()], { type: 'text/plain' });
      const url = URL.createObjectURL(blob);
      window.open(url, '_blank');
      // Revoke after the new tab has had a chance to load it.
      setTimeout(() => URL.revokeObjectURL(url), 60_000);
    } catch (err) {
      alert(`Failed to open report: ${err.message || err}`);
    } finally {
      kmerReportOpening = false;
    }
  }

  const OFFTARGET_PALETTE = [
    '#3b82f6', '#ef4444', '#10b981', '#f59e0b',
    '#8b5cf6', '#ec4899', '#14b8a6', '#f97316',
    '#6366f1', '#84cc16', '#06b6d4', '#a855f7'
  ];
  function offtargetColorFor(i) {
    if (i < OFFTARGET_PALETTE.length) return OFFTARGET_PALETTE[i];
    // Beyond the fixed palette, generate distinct hues via the golden angle.
    const hue = ((i - OFFTARGET_PALETTE.length) * 137.508) % 360;
    return `hsl(${hue.toFixed(1)}, 60%, 55%)`;
  }

  const OFFTARGET_TOP_N_MAX = 100;
  let offTargetTopN = 50;
  let offTargetCustomInput = 50;
  let offTargetSummary = { groups: [], totalGroups: 0, totalAlignments: 0, hiddenGroups: 0, excludedUnknown: 0, groupBy: null };
  let loadingOffTargetSummary = false;

  async function fetchOffTargetSummary() {
    loadingOffTargetSummary = true;
    try {
      const res = await fetch(`/jobs/${jobId}/offtarget-summary?top_n=${offTargetTopN}`);
      if (!res.ok) return;
      const data = await res.json();
      const groups = (data.groups || []).map((g, i) => ({
        label: g.label,
        count: g.count,
        color: offtargetColorFor(i)
      }));
      const totalAlignments = data.total_alignments || 0;
      const hiddenGroups = data.hidden_groups || 0;
      const topSum = groups.reduce((s, g) => s + g.count, 0);
      const othersCount = Math.max(0, totalAlignments - topSum);
      if (hiddenGroups > 0 && othersCount > 0) {
        groups.push({
          label: `Others (${hiddenGroups.toLocaleString()})`,
          count: othersCount,
          color: '#9ca3af'
        });
      }
      offTargetSummary = {
        groups,
        totalGroups: data.total_groups || 0,
        totalAlignments,
        hiddenGroups,
        excludedUnknown: data.excluded_unknown || 0,
        groupBy: data.group_by || null
      };
    } catch (e) {
      console.error('Failed to fetch off-target summary:', e);
    } finally {
      loadingOffTargetSummary = false;
    }
  }

  function applyOffTargetCustom() {
    let v = parseInt(offTargetCustomInput, 10);
    if (!Number.isFinite(v) || v <= 0) return;
    if (v > OFFTARGET_TOP_N_MAX) v = OFFTARGET_TOP_N_MAX;
    offTargetCustomInput = v;
    if (v === offTargetTopN) return;
    offTargetTopN = v;
    fetchOffTargetSummary();
  }

  function onOffTargetCustomKey(event) {
    if (event.key === 'Enter') applyOffTargetCustom();
  }

  function buildPieSlices(groups, radius, totalForPercent) {
    const sliceTotal = groups.reduce((s, g) => s + g.count, 0);
    if (sliceTotal === 0) return [];
    const denom = totalForPercent && totalForPercent > 0 ? totalForPercent : sliceTotal;
    if (groups.length === 1) {
      return [{ ...groups[0], path: null, isFullCircle: true, percent: (groups[0].count / denom) * 100 }];
    }
    let cumulative = 0;
    return groups.map((g) => {
      const startAngle = (cumulative / sliceTotal) * 2 * Math.PI;
      cumulative += g.count;
      const endAngle = (cumulative / sliceTotal) * 2 * Math.PI;
      const x1 = radius * Math.sin(startAngle);
      const y1 = -radius * Math.cos(startAngle);
      const x2 = radius * Math.sin(endAngle);
      const y2 = -radius * Math.cos(endAngle);
      const largeArc = endAngle - startAngle > Math.PI ? 1 : 0;
      const path = `M 0 0 L ${x1.toFixed(3)} ${y1.toFixed(3)} A ${radius} ${radius} 0 ${largeArc} 1 ${x2.toFixed(3)} ${y2.toFixed(3)} Z`;
      return { ...g, path, isFullCircle: false, percent: (g.count / denom) * 100 };
    });
  }

  $: offTargetSlices = buildPieSlices(offTargetSummary.groups, 100, offTargetSummary.totalAlignments);
  $: offTargetGroupBy = offTargetSummary.groupBy
    || ((info?.species === 'human' || info?.species === 'mouse') ? 'gene' : 'species');
  $: offTargetNoun = offTargetGroupBy === 'species' ? 'species' : offTargetGroupBy;
  $: offTargetNounPlural = offTargetGroupBy === 'species' ? 'species' : offTargetGroupBy + 's';

  $: paginatedProbes = filteredSafeProbes.slice(
    (safeProbesPage - 1) * safeProbesPageSize,
    safeProbesPage * safeProbesPageSize
  );
  $: safeProbesTotalPages = Math.ceil(filteredSafeProbes.length / safeProbesPageSize) || 1;
  $: isPrevDisabled = safeProbesPage <= 1;
  $: isNextDisabled = safeProbesPage >= safeProbesTotalPages;
  $: safeProbesPageNumbers = getSafeProbesPageNumbers(safeProbesPage, safeProbesTotalPages);

  $: groupedAlignments = groupAlignmentsByProbe(displayedAlignments);

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
  
  function alignmentLabel(aln) {
    const isHumanMouse = info?.species === 'human' || info?.species === 'mouse';
    const v = isHumanMouse ? aln.gene_id : aln.species;
    if (v && v !== 'Unknown') return v;
    return aln.target_transcript || '-';
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
        const rawSeq = aln.sequence ? aln.sequence.toUpperCase().replace(/-/g, '') : '';
        currentGroup = {
          probe_id: aln.probe_id,
          sequence: rawSeq,
          gc_content: calculateGC(rawSeq),
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
        species: aln.species,
        sequence: aln.sequence,
        target_sequence: aln.target_sequence
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
              if (parts.length >= 8) {
                const lineProbeId = parts[0]; // probe_id is in column 1
                const lineProbeBase = lineProbeId.split('|')[0].trim();

                // Extract probe numbers for exact numeric comparison
                const baseProbeNum = baseProbeId.match(/probe_(\d+)/)?.[1];
                const lineProbeNum = lineProbeBase.match(/probe_(\d+)/)?.[1];

                // Exact match: compare full probe_id OR exact probe number
                if (lineProbeId === featureData.probe_id || baseProbeNum === lineProbeNum) {
                  const scoringData = {
                    probe_id: parts[0],
                    sequence: parts[1],
                    tm: parseFloat(parts[2]),
                    gc_content: parseFloat(parts[3]),
                    length: parseInt(parts[4]),
                    complexity: parseFloat(parts[5]),
                    sec_struct: parseFloat(parts[6]),
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
      
      // If probe failed k-mer safety analysis, show the failure message instead of alignments
      if (featureData.status === 'medium_risk' && featureData.description && featureData.description.includes('k-mer')) {
        expandedProbeAlignments = [];
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
          
          if (parts.length >= 8) {
            const probe = {
              probe_id: parts[0],
              sequence: parts[1],
              tm: parseFloat(parts[2]),
              gc_content: parseFloat(parts[3]),
              length: parseInt(parts[4]),
              complexity: parseFloat(parts[5]),
              sec_struct: parseFloat(parts[6]),
            };

            probes.push(probe);
          }
        }
        
        safeProbes = probes;
        filteredSafeProbes = [...safeProbes];
      }
    } catch (e) {
      console.error('Failed to fetch safe probes:', e);
    } finally {
      loadingSafeProbes = false;
      safeProbesFetched = true;
    }
  }

  function filterSafeProbes() {
    filteredSafeProbes = safeProbes.filter(probe => {

      if (tmFilter.trim() !== '') {
        const val = tmFilter.trim();
        const rangeMatch = val.match(/^(\d+\.?\d*)\s*-\s*(\d+\.?\d*)$/);
        if (rangeMatch) {
          const lo = parseFloat(rangeMatch[1]);
          const hi = parseFloat(rangeMatch[2]);
          if (probe.tm < lo || probe.tm > hi) return false;
        } else {
          if (!probe.tm.toFixed(1).startsWith(val)) return false;
        }
      }

      if (gcFilter.trim() !== '') {
        const val = gcFilter.trim();
        const rangeMatch = val.match(/^(\d+\.?\d*)\s*-\s*(\d+\.?\d*)$/);
        if (rangeMatch) {
          const lo = parseFloat(rangeMatch[1]);
          const hi = parseFloat(rangeMatch[2]);
          if (probe.gc_content < lo || probe.gc_content > hi) return false;
        } else {
          if (!probe.gc_content.toFixed(1).startsWith(val)) return false;
        }
      }
      
      return true;
    });
    safeProbesPage = 1;
  }

  async function toggleSafeProbes() {
    if (safeProbes.length === 0) {
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

  async function fetchAlignments() {
    loadingAlignments = true;
    try {
      const url = `/jobs/${jobId}/alignments?page=1&page_size=${ALIGNMENT_PREVIEW_SIZE}`;
      const alignmentsRes = await fetch(url);
      if (alignmentsRes.ok) {
        const alignData = await alignmentsRes.json();
        alignments = alignData.alignments || [];
        totalAlignments = alignData.total_alignments || 0;
        alignmentsPage = 1;
      }
    } catch (e) {
      console.error('Failed to fetch alignments:', e);
    } finally {
      loadingAlignments = false;
    }
  }

  $: alignmentsTotalPages = Math.max(1, Math.ceil(alignments.length / ALIGNMENTS_PAGE_SIZE));
  $: displayedAlignments = alignments.slice(
    (alignmentsPage - 1) * ALIGNMENTS_PAGE_SIZE,
    alignmentsPage * ALIGNMENTS_PAGE_SIZE
  );

  function goToAlignmentsPage(p) {
    if (p < 1 || p > alignmentsTotalPages) return;
    alignmentsPage = p;
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
      submittedAt = data.submitted_at;
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

        if (status === 'SUCCESS' && inputType === 'probe_sequence') {
          fetchOffTargetSummary();
        }
        await fetchAlignments();
        
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

  function getStatusLabel(probe) {
    if (probe.status === 'safe') {
      return 'No off-target alignments';
    } else if (probe.status === 'high_risk') {
      return `High Risk (has exact matches with other gene)`;
    } else if (probe.status === 'medium_risk') {
      if (probe.mismatches === undefined || probe.mismatches === null) {
        return 'Medium Risk (Failed k-mer safety analysis)';
      }
      return `Medium Risk (off-target matches with 1-2 mismatches)`;
    } else if (probe.status === 'no_alignment') {
      return 'No alignment to source gene';
    } else {
      return 'Unknown';
    }
  }

  function getStatusColor(status) {
    if (status === 'safe') return '#10b981';
    if (status === 'high_risk') return '#dc2626';
    if (status === 'medium_risk') return '#ea580c';
    if (status === 'no_alignment') return '#6b7280';
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
  $: if (status === 'SUCCESS' && info?.stats?.non_aligned_probes > 0 && safeProbes.length === 0 && !safeProbesFetched) {
    fetchSafeProbes();
  }

  function downloadFilteredProbes() {
    if (filteredSafeProbes.length === 0) return;
    const isFiltered = tmFilter.trim() !== '' || gcFilter.trim() !== '';
    const headers = ['Probe ID', 'Sequence', 'Tm', 'GC%', 'Complexity', 'Sec. Struct'];
    const rows = filteredSafeProbes.map(p => [
      p.probe_id, p.sequence, p.tm.toFixed(1),
      p.gc_content.toFixed(1), p.complexity.toFixed(2), p.sec_struct.toFixed(2)
    ].join('\t'));
    const tsv = [headers.join('\t'), ...rows].join('\n');
    const blob = new Blob([tsv], { type: 'text/tab-separated-values' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `${isFiltered ? 'filtered' : 'safe'}_probes_${filteredSafeProbes.length}.tsv`;
    a.click();
    URL.revokeObjectURL(url);
  }
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
        <th style="text-align:left; border-bottom:1px solid #ddd; padding:12px 16px">Submitted at (UTC)</th>
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
        <td style="padding:12px 16px; vertical-align:top">{submittedAt ? submittedAt.split('T')[1].split('.')[0] : '-'}</td>
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
    {#if inputType === 'probe_sequence' && info.stats.input_probes !== null && info.stats.input_probes !== undefined}
      {@const rejections = [
        { label: 'GC', value: info.stats.gc_rejected ?? 0 },
        { label: 'Tm', value: info.stats.tm_rejected ?? 0 },
        { label: 'homopolymer', value: info.stats.homopolymer_rejected ?? 0 },
      ].filter(r => r.value > 0)}
      {#if rejections.length > 0}
        <p style="margin:0 0 12px 0; font-size:14px; color:#374151;">
          From <strong>{formatNumber(info.stats.input_probes)}</strong> input probes, {#each rejections as r, idx}<strong>{formatNumber(r.value)}</strong> {idx === 0 ? 'rejected ' : ''}by {r.label}{idx < rejections.length - 1 ? ', ' : ''}{/each}.
        </p>
      {/if}
    {/if}
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
              {#if info.stats.gc_passed === 0}
                None of the probes passed the GC content filter (40-80%).
                This may occur if the input sequence has extreme GC content.
                Consider adjusting the probe length or target a different region of the sequence.
              {:else if info.stats.tm_rejected > 0 && info.stats.tm_rejected === info.stats.gc_passed}
                None of the probes passed the Tm filter.
                All {info.stats.gc_passed} probes that passed the GC filter were rejected by the Tm range.
                Consider widening the Tm range or adjusting the probe length.
              {:else if info.stats.homopolymer_rejected > 0 && info.stats.homopolymer_rejected === info.stats.gc_passed}
                None of the probes passed the homopolymer filter.
                All {info.stats.gc_passed} probes that passed the GC filter contain homopolymer runs.
                Consider adjusting the probe length or target a different region of the sequence.
              {:else}
                No probes passed the filtering criteria.
                Consider adjusting the probe length, Tm range, or target a different region of the sequence.
              {/if}
            </div>
          </div>
        </div>
      </div>
    {:else}
    <table
      style="
        width:100%;
        border-collapse:collapse;
        font-size:14px;
        color:#111;
        border:1px solid #d1d5db;
      "
    >
      <thead>
        <tr style="background:#f3f4f6;">
          <th style="padding:10px; border:1px solid #d1d5db; text-align:center; font-weight:700;">Total Candidate Probes</th>
          <th style="padding:10px; border:1px solid #d1d5db; text-align:center; font-weight:700;">Total Potential On-target Probes</th>
          <th style="padding:10px; border:1px solid #d1d5db; text-align:center; font-weight:700;">Total Safe Probes </th>
          <th style="padding:10px; border:1px solid #d1d5db; text-align:center; font-weight:700;">K-mer Filter Result</th>
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
                {#if f.filename.toLowerCase() === 'safe_probes.fa'}
                  <a
                    href={`/jobs/${jobId}/download/${f.filename}`}
                    class="download-link"
                  >
                    <svg width="16" height="16" viewBox="0 0 24 24" fill="currentColor">
                      <path d="M19 9h-4V3H9v6H5l7 7 7-7zM5 18v2h14v-2H5z"/>
                    </svg>
                    <span>Download Potential On-target Probes</span>
                  </a>
                {/if}
              {/each}
            {/if}
          </td>
          <td style="padding:12px; border:1px solid #e5e7eb; text-align:center; line-height:1.6;">
            <div>{formatNumber(info.stats.safe_probes || 0)}</div>
            {#if info.stats.safe_probes > 0}
            {#each files as f}
              {#if f.filename.toLowerCase() === 'safe_probes_scores.txt'}
                <a
                  href={`/jobs/${jobId}/download/${f.filename}`}
                  class="download-link"
                >
                <svg width="16" height="16" viewBox="0 0 24 24" fill="currentColor">
                  <path d="M19 9h-4V3H9v6H5l7 7 7-7zM5 18v2h14v-2H5z"/>
                </svg>
                <span>Download Safe Probes</span>
                </a>
              {/if}
            {/each}
            {/if}
          </td>

          <td style="padding:12px; border:1px solid #e5e7eb; text-align:center; line-height:1.6;">
            {#if info.stats.non_aligned_probes > 0}
            {#each files as f}
              {#if f.filename.toLowerCase() === 'kmer_matches_report.txt'}
                <a
                  href={`/jobs/${jobId}/download/${f.filename}`}
                  class="download-link"
                >
                 <svg width="16" height="16" viewBox="0 0 24 24" fill="currentColor">
                  <path d="M19 9h-4V3H9v6H5l7 7 7-7zM5 18v2h14v-2H5z"/>
                </svg>
                <span>Download K-mer Analysis Report</span>
                </a>
                <button
                  type="button"
                  class="download-link"
                  style="border:none; cursor:pointer;"
                  disabled={kmerReportOpening}
                  on:click={openKmerReportInTab}
                >
                  <svg width="16" height="16" viewBox="0 0 24 24" fill="currentColor">
                    <path d="M14 3v2h3.59l-9.83 9.83 1.41 1.41L19 6.41V10h2V3h-7zM19 19H5V5h7V3H5c-1.11 0-2 .9-2 2v14c0 1.1.89 2 2 2h14c1.1 0 2-.9 2-2v-7h-2v7z"/>
                  </svg>
                  <span>{kmerReportOpening ? 'Opening…' : 'View K-mer Analysis Report'}</span>
                </button>
              {/if}
            {/each}
            {/if}
          </td>
        </tr>
      </tbody>
    </table>
    {/if}
  </div>
{/if}
  
  {#if info && info.stats && info.stats.non_aligned_probes > 0}
    
    <div style="margin-top:1.5rem;">
            <h3 style="margin:0 0 12px 0; font-weight:600; font-size:16px;">
        Safe Probe Scoring Details ({formatNumber(info.stats.safe_probes)} probes)
      </h3>

      {#if true}
        <div transition:slide={{ duration: 300, easing: quintOut }} style="margin-top:12px; border:1px solid #e5e7eb; border-radius:6px; overflow:hidden;">
        
          {#if loadingSafeProbes}
            <div style="padding:40px; text-align:center;">
              <div class="spinner" style="width:32px; height:32px; margin:0 auto 12px;"></div>
              <p style="color:#6b7280; margin:0;">Loading probe scores...</p>
            </div>
          {:else if safeProbesFetched && safeProbes.length === 0}
            <div style="padding:40px; text-align:center;">
              <p style="color:#6b7280; margin:0;">No safe probes found for scoring.</p>
              {#if info && info.stats && info.stats.non_aligned_probes > 0 && (info.stats.safe_probes === 0 || !info.stats.safe_probes)}
                <div style="
                  margin-top:16px;
                  padding:14px 20px;
                  background:#fffbeb;
                  border:1px solid #fbbf24;
                  border-radius:8px;
                  text-align:left;
                  font-size:13px;
                  color:#92400e;
                  line-height:1.6;
                ">
                  <strong>Note:</strong> {formatNumber(info.stats.non_aligned_probes)} potential on-target probes were identified from alignment,
                  but all were filtered out during the k-mer filter (k-mer matches found in other genes).
                  You can try re-running with a higher k-mer length to relax the filtering stringency.

                </div>
              {/if}
            </div>
          {:else}
            <div style="padding:16px; background:#f9fafb; border-bottom:1px solid #e5e7eb;">
              <div style="display:flex; justify-content:space-between; align-items:center;">
                <div style="font-size:13px; color:#6b7280;">
                  Showing {((safeProbesPage - 1) * safeProbesPageSize) + 1}-{Math.min(safeProbesPage * safeProbesPageSize, filteredSafeProbes.length)} of {filteredSafeProbes.length} probes
                </div>
                <div style="display:flex; align-items:center; gap:12px;">
                  <button
                    on:click={downloadFilteredProbes}
                    disabled={filteredSafeProbes.length === 0}
                    class="download-link"
                  >
                    <svg width="16" height="16" viewBox="0 0 24 24" fill="currentColor">
                      <path d="M19 9h-4V3H9v6H5l7 7 7-7zM5 18v2h14v-2H5z"/>
                    </svg>
                    <span>{tmFilter.trim() || gcFilter.trim() ? 'Download Filtered' : 'Download All'} ({filteredSafeProbes.length})</span>
                  </button>
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
            <div style="overflow-x:auto;">
              <table style="width:100%; border-collapse:collapse; font-size:12px;">
                <thead>
                  <tr style="background:#f9fafb; border-bottom:2px solid #d1d5db;">
                    <th style="padding:10px; text-align:left;">Probe ID</th>
                    <th style="padding:10px; text-align:left;">Sequence</th>
                    <th style="padding:10px; text-align:center; position:relative;">
                      <div style="margin-bottom:4px;">Tm (°C)</div>
                      <input
                        type="text"
                        placeholder="e.g. 46.9 or 45-48"
                        bind:value={tmFilter}
                        on:input={() => filterSafeProbes()}
                        style="width:100%; padding:2px 4px; border:1px solid #d1d5db; border-radius:3px; font-size:11px; background:white; box-sizing:border-box;"
                      />
                    </th>
                    <th style="padding:10px; text-align:center; position:relative;">
                      <div style="margin-bottom:4px;">GC%</div>
                      <input
                        type="text"
                        placeholder="e.g. 55.0 or 40-60"
                        bind:value={gcFilter}
                        on:input={() => filterSafeProbes()}
                        style="width:100%; padding:2px 4px; border:1px solid #d1d5db; border-radius:3px; font-size:11px; background:white; box-sizing:border-box;"
                      />
                    </th>
                    <th style="padding:10px; text-align:center;">Complexity</th>
                    <th style="padding:10px; text-align:center;">Sec. Struct</th>
                  </tr>
                </thead>
                <tbody>
                  {#each paginatedProbes as probe, idx (idx + safeProbesPage * safeProbesPageSize)}
                    <tr style="border-bottom:1px solid #e5e7eb;">
                      <td style="padding:10px; font-family:monospace; color:#1f2937;">{probe.probe_id}</td>
                      <td style="padding:10px; font-family:monospace; font-size:12px; white-space:nowrap;">{probe.sequence}</td>
                      <td style="padding:10px; text-align:center; font-family:monospace;">{probe.tm.toFixed(1)}</td>
                      <td style="padding:10px; text-align:center; font-family:monospace;">{probe.gc_content.toFixed(1)}</td>
                      <td style="padding:10px; text-align:center; font-family:monospace;">{probe.complexity.toFixed(2)}</td>
                      <td style="padding:10px; text-align:center; font-family:monospace;">{probe.sec_struct.toFixed(2)}</td>
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
  
  {#if status === 'SUCCESS' && sequenceLength > 0 && inputType === 'gene_sequence' && info.stats.candidate_probes > 0}
    <div style="margin-top:1.5rem;">
      <h3 style="margin:0 0 12px 0; font-weight:600; font-size:16px;">Probe Alignment Browser</h3>

      {#if probeRiskSummary}
        <div style="margin-bottom:12px;">
          <div style="margin:0 0 6px 0; font-size:13px; font-weight:500; color:#374151;">
            Probe Risk Summary
          </div>
          <div style="display:flex; flex-wrap:wrap; gap:8px; font-size:13px;">
            <span title="No off-target alignments and clean k-mer profile" style="padding:4px 10px; border-radius:12px; background:#d1fae5; color:#065f46; font-weight:600;">
              Safe: {probeRiskSummary.safe.toLocaleString()}
            </span>
            <span title="At least one k-mer matches off-target sequence" style="padding:4px 10px; border-radius:12px; background:#fed7aa; color:#9a3412; font-weight:600;">
              Medium risk — failed k-mer: {probeRiskSummary.medium_risk_kmer.toLocaleString()}
            </span>
            <span title="Off-target alignment with 1–2 mismatches" style="padding:4px 10px; border-radius:12px; background:#fed7aa; color:#9a3412; font-weight:600;">
              Medium risk — 1–2 mismatches: {probeRiskSummary.medium_risk_mismatch.toLocaleString()}
            </span>
            <span title="Exact off-target alignment to another gene" style="padding:4px 10px; border-radius:12px; background:#fecaca; color:#991b1b; font-weight:600;">
              High risk: {probeRiskSummary.high_risk.toLocaleString()}
            </span>
            <span title="Probe does not align to the source gene at all" style="padding:4px 10px; border-radius:12px; background:#e5e7eb; color:#374151; font-weight:600;">
              No alignment: {probeRiskSummary.no_alignment.toLocaleString()}
            </span>
            <span title="All probes shown in the browser" style="padding:4px 10px; border-radius:12px; background:#e5e7eb; color:#111827; font-weight:600;">
              Total: {probeRiskSummary.total.toLocaleString()}
            </span>
          </div>
          <p style="margin:6px 0 0 0; font-size:12px; color:#6b7280; line-height:1.5;">
            Probes are grouped by off-target risk. <strong>Safe</strong> = no off-target hits. <strong>Medium risk</strong> = either a failing k-mer match or an off-target alignment with 1–2 mismatches. <strong>High risk</strong> = exact off-target match to another gene. <strong>No alignment</strong> = probe doesn't bind the source gene at all.
          </p>
        </div>
      {/if}

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
              <span><strong>Length:</strong> {(expandedProbe.end - expandedProbe.start + 1).toLocaleString()} bp</span>
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
            {#if expandedProbe.status === 'medium_risk'}
              <div style="padding:20px; background:#fef3c7; border:1px solid #f59e0b; border-radius:6px; text-align:center; color:#92400e;">
                <p style="margin:0; font-weight:600;">Failed k-mer safety analysis</p>
                <p style="margin:8px 0 0 0; font-size:13px;">No full-length off-target matches, but contains short k-mer segments shared with other genes.</p>
              </div>
            {:else if expandedProbe.status === 'no_alignment'}
              <div style="padding:20px; background:#f3f4f6; border:1px solid #9ca3af; border-radius:6px; text-align:center; color:#374151;">
                <p style="margin:0; font-weight:600;">No alignment to source gene</p>
                <p style="margin:8px 0 0 0; font-size:13px;">This probe doesn't bind the target — not usable for hybridization.</p>
              </div>
            {:else}
              <div style="padding:20px; background:#ecfdf5; border:1px solid #10b981; border-radius:6px; text-align:center; color:#065f46;">
                <p style="margin:0; font-weight:600;">No off-target alignments</p>
              </div>
            {/if}
            
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
                          {alignmentLabel(aln)}
                        </td>
                      {:else}
                        <td style="padding:8px; font-family:monospace; color:#374151; font-size:12px;">
                          {alignmentLabel(aln)}
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
  {#if status === 'SUCCESS' && inputType === 'probe_sequence' && info.stats.candidate_probes > 0}
    {#if offTargetSummary.groups.length > 0}
      <div style="margin-top:1.5rem;">
        <div style="display:flex; align-items:center; justify-content:space-between; gap:12px; flex-wrap:wrap; margin-bottom:4px;">
          <h3 style="margin:0; font-weight:600; font-size:16px;">
            Off-target distribution
          </h3>
          <div style="display:flex; align-items:center; gap:8px; font-size:13px; color:#374151; flex-wrap:wrap;">
            <label for="offtarget-top-n">Show top:</label>
            <input
              id="offtarget-top-n"
              type="number"
              min="1"
              max={OFFTARGET_TOP_N_MAX}
              bind:value={offTargetCustomInput}
              on:keydown={onOffTargetCustomKey}
              disabled={loadingOffTargetSummary}
              style="width:72px; padding:4px 6px; border:1px solid #d1d5db; border-radius:4px; font-size:13px;"
            />
            <button
              on:click={applyOffTargetCustom}
              disabled={loadingOffTargetSummary || !Number.isFinite(parseInt(offTargetCustomInput, 10)) || parseInt(offTargetCustomInput, 10) <= 0 || parseInt(offTargetCustomInput, 10) === offTargetTopN}
              style="padding:4px 10px; border:1px solid #d1d5db; border-radius:4px; font-size:13px; background:white; cursor:pointer;"
            >
              Apply
            </button>
            {#if loadingOffTargetSummary}
              <span style="color:#6b7280;">loading…</span>
            {/if}
          </div>
        </div>
        <p style="margin:0 0 12px 0; font-size:13px; color:#6b7280;">
          {#if offTargetSummary.hiddenGroups > 0}
            {@const topShown = offTargetSummary.groups.length - 1}
            The pie shows the top {topShown} by alignment count; the remaining {offTargetSummary.hiddenGroups.toLocaleString()} {offTargetSummary.hiddenGroups === 1 ? offTargetNoun : offTargetNounPlural} are grouped as "Others".
          {/if}
        </p>
        <div style="display:flex; gap:24px; align-items:center; flex-wrap:wrap; padding:16px; background:#f9fafb; border:1px solid #e5e7eb; border-radius:6px;">
          <svg viewBox="-110 -110 220 220" width="220" height="220" style="flex-shrink:0;">
            {#each offTargetSlices as slice}
              {#if slice.isFullCircle}
                <circle cx="0" cy="0" r="100" fill={slice.color} stroke="white" stroke-width="1.5">
                  <title>{slice.label}: {slice.count.toLocaleString()} alignment{slice.count !== 1 ? 's' : ''} ({slice.percent.toFixed(2)}% of total)</title>
                </circle>
              {:else}
                <path d={slice.path} fill={slice.color} stroke="white" stroke-width="1.5">
                  <title>{slice.label}: {slice.count.toLocaleString()} alignment{slice.count !== 1 ? 's' : ''} ({slice.percent.toFixed(2)}% of total)</title>
                </path>
              {/if}
            {/each}
          </svg>
          <div style="flex:1; min-width:240px; display:grid; grid-template-columns:repeat(auto-fill, minmax(220px, 1fr)); gap:6px 16px;">
            {#each offTargetSlices as slice}
              <div style="display:flex; align-items:center; gap:8px; font-size:13px;">
                <span style="display:inline-block; width:12px; height:12px; border-radius:2px; background:{slice.color}; flex-shrink:0;"></span>
                <span style="flex:1; color:#1f2937; overflow:hidden; text-overflow:ellipsis; white-space:nowrap;" title={slice.label}>
                  {slice.label}
                </span>
                <span style="color:#6b7280; font-variant-numeric:tabular-nums;">
                  {slice.count.toLocaleString()} ({slice.percent.toFixed(2)}%)
                </span>
              </div>
            {/each}
          </div>
        </div>
      </div>
    {/if}
    <div style="margin-top:1.5rem;">
      <h3 style="margin:0 0 12px 0; font-weight:600; font-size:16px;">Probe Alignment Browser</h3>

      {#if alignments.length > 0}
        <div style="margin-bottom:12px; display:flex; align-items:center; gap:12px; flex-wrap:wrap;">
          <span style="font-size:14px; color:#374151;">
            Showing first {alignments.length.toLocaleString()} of {totalAlignments.toLocaleString()} alignments. Download the TSV for full data.
          </span>
          <a
            href={`/jobs/${jobId}/alignments/download`}
            class="download-link"
            style="margin-left:auto;"
          >
            <svg width="16" height="16" viewBox="0 0 24 24" fill="currentColor">
              <path d="M19 9h-4V3H9v6H5l7 7 7-7zM5 18v2h14v-2H5z"/>
            </svg>
            <span>Download all alignments ({totalAlignments.toLocaleString()})</span>
          </a>
        </div>
      {/if}

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
                    {/if}
                    <td style="padding:10px; vertical-align:top;">
                      <div style="font-family:monospace; font-size:11px; line-height:1.5; word-break:break-all;">
                        {#each formatSequenceWithHighlights(aln.target_sequence || aln.sequence) as {char, isLowercase, isDash}}
                          <span class:lowercase={isLowercase} class:dash={isDash}>{char}</span>
                        {/each}
                      </div>
                    </td>
                    <td style="padding:10px; font-family:monospace; color:#374151;">
                      {aln.target_transcript}
                    </td>
                    <td style="padding:8px; color:#1f2937; font-weight:600;">
                      {alignmentLabel(aln)}
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
        {#if alignmentsTotalPages > 1}
          {@const prevDisabled = alignmentsPage <= 1}
          {@const nextDisabled = alignmentsPage >= alignmentsTotalPages}
          <div style="margin-top:12px; display:flex; justify-content:center; align-items:center; gap:8px; flex-wrap:wrap; font-size:13px;">
            <button
              on:click={() => goToAlignmentsPage(alignmentsPage - 1)}
              disabled={prevDisabled}
              style="padding:6px 10px; border:1px solid #d1d5db; border-radius:4px; background:{prevDisabled ? '#e5e7eb' : '#3b82f6'}; color:{prevDisabled ? '#9ca3af' : 'white'}; cursor:{prevDisabled ? 'not-allowed' : 'pointer'}; font-weight:500;"
            >Previous</button>
            {#each Array(alignmentsTotalPages) as _, i}
              {@const num = i + 1}
              {@const active = num === alignmentsPage}
              <button
                on:click={() => goToAlignmentsPage(num)}
                style="padding:6px 12px; border:1px solid #d1d5db; border-radius:4px; background:{active ? '#3b82f6' : 'white'}; color:{active ? 'white' : '#374151'}; cursor:pointer; font-weight:{active ? 600 : 400};"
              >{num}</button>
            {/each}
            <button
              on:click={() => goToAlignmentsPage(alignmentsPage + 1)}
              disabled={nextDisabled}
              style="padding:6px 10px; border:1px solid #d1d5db; border-radius:4px; background:{nextDisabled ? '#e5e7eb' : '#3b82f6'}; color:{nextDisabled ? '#9ca3af' : 'white'}; cursor:{nextDisabled ? 'not-allowed' : 'pointer'}; font-weight:500;"
            >Next</button>
          </div>
        {/if}
      {:else if loadingAlignments}
        <div style="padding:40px; display:flex; justify-content:center; align-items:center; gap:10px; color:#6b7280; background:#f9fafb; border-radius:6px;">
          <div class="spinner"></div>
          <span>Loading alignments...</span>
        </div>
      {:else}
        <div style="padding:40px; text-align:center; color:#6b7280; background:#f9fafb; border-radius:6px;">
          No off-target alignment.
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
    gap: 6px;
    color: white;
    font-size: 13px;
    font-weight: 600;
    line-height: 1.2;
    text-decoration: none;
    margin-top: 8px;
    padding: 8px 16px;
    border-radius: 8px;
    background: linear-gradient(135deg, #3b82f6 0%, #2563eb 100%);
    border: none;
    cursor: pointer;
    transition: all 0.2s ease;
    box-shadow: 0 2px 8px rgba(59, 130, 246, 0.3);
  }
  .download-link:hover,
  .download-link:focus {
    transform: translateY(-2px);
    box-shadow: 0 4px 12px rgba(59, 130, 246, 0.4);
  }
  .download-link:disabled {
    background: #9ca3af;
    box-shadow: none;
    cursor: not-allowed;
  }
  .download-link:disabled:hover,
  .download-link:disabled:focus {
    transform: none;
    box-shadow: none;
  }
</style>