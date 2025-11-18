<script>
  import { createEventDispatcher } from 'svelte';
  
  export let currentView = 'home';
  export let hasJob = false;
  
  const dispatch = createEventDispatcher();
  
  let mobileMenuOpen = false;
  
  function navigate(view) {
    dispatch('navigate', { view });
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
      class="nav-button"
      class:active={currentView === 'home'}
    >
     Home
    </button>
    
    <div class="nav-right">
      <button
        on:click={() => navigate('upload')}
        class="nav-button"
        class:active={currentView === 'upload'}
      >
       Design Probes
      </button>
      
      <button
        on:click={() => hasJob ? navigate('results') : null}
        disabled={!hasJob}
        class="nav-button"
        class:active={currentView === 'results'}
        class:disabled={!hasJob}
      >
       Results
      </button>
      
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
      class="nav-button-mobile"
      class:active={currentView === 'home'}
    >
     Home
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
        on:click={() => navigate('upload')}
        class="mobile-menu-item"
        class:active={currentView === 'upload'}
      >
       Design Probes
      </button>
      
      <button
        on:click={() => hasJob ? navigate('results') : null}
        disabled={!hasJob}
        class="mobile-menu-item"
        class:active={currentView === 'results'}
        class:disabled={!hasJob}
      >
       Results
      </button>
      
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