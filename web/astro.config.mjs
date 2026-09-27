// @ts-check
import { defineConfig } from 'astro/config';

// Served from GitHub Pages at https://ejhong.github.io/sar/
export default defineConfig({
  site: 'https://ejhong.github.io',
  base: '/sar',
  trailingSlash: 'ignore',
  output: 'static',
  devToolbar: { enabled: false },
  build: { format: 'directory' },
  redirects: {
    '/underworld.html': '/sar/underworld/',
    '/field.html': '/sar/archive/field.html',
  },
  vite: {
    build: { chunkSizeWarningLimit: 1200 },
  },
});
