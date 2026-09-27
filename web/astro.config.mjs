// @ts-check
import { defineConfig } from 'astro/config';

// Served from GitHub Pages at https://ejhong.github.io/sar/, from the docs/ folder on main:
// the build is written there and committed (npm run build).
export default defineConfig({
  site: 'https://ejhong.github.io',
  base: '/sar',
  outDir: '../docs',
  trailingSlash: 'ignore',
  output: 'static',
  devToolbar: { enabled: false },
  build: { format: 'directory' },
  vite: {
    build: { chunkSizeWarningLimit: 1200 },
  },
});
