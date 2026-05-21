import { svelte } from '@sveltejs/vite-plugin-svelte';
import { defineConfig, loadEnv } from 'vite';

// Backend (FastAPI) URL for the Vite dev proxy. Override per-environment
// by setting one of these (in order of precedence):
//   1. env var BACKEND_URL  (e.g. BACKEND_URL=http://localhost:8002 npm run dev)
//   2. .env / .env.local file with BACKEND_URL=…  or  VITE_BACKEND_URL=…
//   3. defaults to http://localhost:8000
export default defineConfig(({ mode }) => {
  const env = loadEnv(mode, process.cwd(), '');
  const backend =
    process.env.BACKEND_URL ||
    env.BACKEND_URL ||
    env.VITE_BACKEND_URL ||
    'http://localhost:8002';

  const proxyTo = { target: backend, changeOrigin: true };

  return {
    plugins: [svelte()],
    publicDir: 'static',
    resolve: { alias: { buffer: 'buffer/' } },
    optimizeDeps: { include: ['buffer'] },
    server: {
      port: 5173,
      proxy: {
        '/jobs': proxyTo,
        '/api': proxyTo,
        '/validate-host-token': proxyTo
      }
    }
  };
});