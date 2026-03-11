import { svelte } from '@sveltejs/vite-plugin-svelte';
import { defineConfig } from 'vite';

export default defineConfig({
  plugins: [svelte()],
  publicDir: 'static',
  resolve: {
    alias: {
      buffer: 'buffer/'
    }
  },
  optimizeDeps: {
    include: ['buffer']
  },
  server: {  
    port: 5173,
    proxy: {
      '/jobs': {
        target: 'http://localhost:8000',
        changeOrigin: true
      }
    }
  }
});