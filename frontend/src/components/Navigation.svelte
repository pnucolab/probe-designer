<script>
  import { createEventDispatcher } from 'svelte';
  
  export let currentView = 'home';
  export let hasJob = false;
  export let activeMode = '';

  const dispatch = createEventDispatcher();

  let mobileMenuOpen = false;

  function navigate(view, mode) {
    dispatch('navigate', mode ? { view, mode } : { view });
    mobileMenuOpen = false;
  }
  
  function toggleMobileMenu() {
    mobileMenuOpen = !mobileMenuOpen;
  }
</script>

<nav class="nav-container">
  <div class="nav-desktop">
    <button
      on:click={() => navigate('home')}
      class="nav-brand"
      aria-label="SHARP-FISH home"
    >
      SHARP-FISH
    </button>

    <div class="nav-right">
      <button
        on:click={() => navigate('home')}
        class="nav-button"
        class:active={currentView === 'home'}
      >
       Home
      </button>

      <button
        on:click={() => navigate('upload', 'host')}
        class="nav-button"
        class:active={currentView === 'upload' && activeMode === 'host'}
      >
       Within-Organism Probe Design
      </button>

      <button
        on:click={() => navigate('upload', 'microbe')}
        class="nav-button"
        class:active={currentView === 'upload' && activeMode === 'microbe'}
      >
       Microbial Probe Design
      </button>

      {#if currentView !== 'home'}
        <button
          on:click={() => hasJob ? navigate('results') : null}
          disabled={!hasJob}
          class="nav-button"
          class:active={currentView === 'results'}
          class:disabled={!hasJob}
        >
         Results
        </button>
      {/if}

      <button
        on:click={() => navigate('about')}
        class="nav-button"
        class:active={currentView === 'about'}
      >
       About
      </button>
    </div>
  </div>

  <div class="nav-mobile">
    <button
      on:click={() => navigate('home')}
      class="nav-brand"
      aria-label="SHARP-FISH home"
    >
      SHARP-FISH
    </button>

    <button
      on:click={toggleMobileMenu}
      class="hamburger-button"
      aria-label="Toggle menu"
    >
      {#if mobileMenuOpen}
        ✕
      {:else}
        ☰
      {/if}
    </button>
  </div>

  {#if mobileMenuOpen}
    <div class="mobile-menu">
      <button
        on:click={() => navigate('home')}
        class="mobile-menu-item"
        class:active={currentView === 'home'}
      >
       Home
      </button>

      <button
        on:click={() => navigate('upload', 'host')}
        class="mobile-menu-item"
        class:active={currentView === 'upload' && activeMode === 'host'}
      >
       Within-Organism Probe Design
      </button>

      <button
        on:click={() => navigate('upload', 'microbe')}
        class="mobile-menu-item"
        class:active={currentView === 'upload' && activeMode === 'microbe'}
      >
       Microbial Probe Design
      </button>

      {#if currentView !== 'home'}
        <button
          on:click={() => hasJob ? navigate('results') : null}
          disabled={!hasJob}
          class="mobile-menu-item"
          class:active={currentView === 'results'}
          class:disabled={!hasJob}
        >
         Results
        </button>
      {/if}

      <button
        on:click={() => navigate('about')}
        class="mobile-menu-item"
        class:active={currentView === 'about'}
      >
       About
      </button>
    </div>
  {/if}
</nav>