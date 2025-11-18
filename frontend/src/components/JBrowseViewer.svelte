<script>
  import { onMount, onDestroy } from 'svelte';
  import { createViewState, JBrowseLinearGenomeView } from '@jbrowse/react-linear-genome-view';
  import { reaction } from 'mobx';
  export let jobId;
  export let referenceId = 'reference';
  export let sequenceLength = 0;
  export let onFeatureClick = null;
  
  let container;
  let state;
  let reactRoot;
  $: assembly = {
    name: referenceId,
    sequence: {
      type: 'ReferenceSequenceTrack',
      trackId: `${referenceId}-ReferenceSequenceTrack`,
      adapter: {
        type: 'IndexedFastaAdapter',
        fastaLocation: {
          uri: `/jobs/${jobId}/download/reference.fasta`,
          locationType: 'UriLocation',
        },
        faiLocation: {
          uri: `/jobs/${jobId}/download/reference.fasta.fai`,
          locationType: 'UriLocation',
        },
      },
    },
  };
  
  $: tracks = [
    {
      type: 'FeatureTrack',
      trackId: 'probes-track',
      name: 'Probes (by Risk Level)',
      assemblyNames: [referenceId],
      adapter: {
        type: 'Gff3Adapter',
        gffLocation: {
          uri: `/jobs/${jobId}/download/probes.gff3`,
          locationType: 'UriLocation',
        },
      },
      displays: [
        {
          type: 'LinearBasicDisplay',
          displayId: 'probes-linear-display',
          renderer: {
            type: 'SvgFeatureRenderer',
            color1: 'jexl:get(feature, "color") || "#888"',
            labels: {
              name: 'jexl:get(feature, "name")',
              description: 'jexl:get(feature, "description")',
            },
          },
          mouseoverFeature: false,
        },
      ],
    }
  ];

  $: defaultSession = {
    name: 'Probe Visualization',
    view: {
      id: 'linearGenomeView',
      type: 'LinearGenomeView',
      tracks: [
        {
          id: 'probes-track-view',
          type: 'FeatureTrack',
          configuration: 'probes-track',
          displays: [
            {
              id: 'probes-linear-display-view',
              type: 'LinearBasicDisplay',
              configuration: 'probes-linear-display',
              height: 230,
            },
          ],
        }
      ],
    },
  };
  
  $: location = `${referenceId}:1-${sequenceLength}`;
  
  onMount(async () => {
    if (!container) return;
    
    try {
      const React = await import('react');
      const ReactDOM = await import('react-dom/client');
      
      console.log('Initializing JBrowse with config:', {
        assembly,
        tracks,
        location,
      });
      
      state = createViewState({
        assembly,
        tracks,
        location,
        defaultSession,
        configuration: {
          theme: {
            palette: {
              primary: {
                main: '#3b82f6',
              },
              secondary: {
                main: '#10b981',
              },
            },
          },
        },
        onChange: (patch) => {
          console.log('JBrowse state changed:', patch);
        },
      });
      const element = React.createElement(JBrowseLinearGenomeView, { viewState: state });
      
      reactRoot = ReactDOM.createRoot(container);
      reactRoot.render(element);
      
      console.log('JBrowse initialized successfully');
      
      if (onFeatureClick) {
        setTimeout(() => {
          try {
            const view = state.session.views[0];
            if (view) {
              reaction(
                () => {
                  const widgets = state.session.widgets;
                  const baseFeatureWidget = widgets.get('baseFeature');
                  return baseFeatureWidget?.featureData;
                },
                (featureData) => {
                  if (featureData) {
                    setTimeout(() => {
                      try {
                        const widget = state.session.widgets.get('baseFeature');
                        if (widget) {
                          state.session.hideWidget(widget);
                        }
                      } catch (e) {
                        console.warn('Could not hide widget:', e);
                      }
                    }, 0);
                    
                    console.log('Feature clicked:', featureData);
                    
                    const extractedData = {
                      probe_id: featureData.name || featureData.ID,
                      start: featureData.start,
                      end: featureData.end,
                      status: featureData.risk_level,
                      mismatches: featureData.mismatches,
                      description: featureData.description,
                    };
                    
                    console.log('Extracted feature data:', extractedData);
                    onFeatureClick(extractedData);
                  }
                }
              );
            }
          } catch (err) {
            console.warn('Could not set up feature click handler:', err);
          }
        }, 1000);
      }
      
    } catch (error) {
      console.error('Failed to initialize JBrowse:', error);
    }
  });
  
  onDestroy(() => {
    if (reactRoot) {
      reactRoot.unmount();
    }
  });
</script>

<div bind:this={container} class="jbrowse-container"></div>

<style>
  .jbrowse-container {
    width: 100%;
    height: 400px; 
    min-height: 400px;
    resize: vertical;
    border: 1px solid #d1d5db;
    border-radius: 6px;
    overflow: auto;
    background: #fff;
  }
  
  :global(.jbrowse-container .MuiPaper-root) {
    box-shadow: none !important;
  }
  :global(.MuiDialog-root),
  :global(.MuiPaper-root[role="dialog"]),
  :global([class*="BaseFeatureDetail"]) {
    display: none !important;
  }
</style>