import type { GridBlock, SiteScene } from './types';

/** Bilinear sample of a grid block at (x, y), clamped at the edges. */
export function sampleGrid(g: GridBlock, x: number, y: number): number {
  const fx = Math.min(Math.max((x - g.x0) / g.dx, 0), g.nx - 1);
  const fy = Math.min(Math.max((y - g.y0) / g.dx, 0), g.ny - 1);
  const ix = Math.min(Math.floor(fx), g.nx - 2);
  const iy = Math.min(Math.floor(fy), g.ny - 2);
  const tx = fx - ix;
  const ty = fy - iy;
  const at = (i: number, j: number) => g.z[i + g.nx * j];
  return (
    at(ix, iy) * (1 - tx) * (1 - ty) + at(ix + 1, iy) * tx * (1 - ty) + at(ix, iy + 1) * (1 - tx) * ty + at(ix + 1, iy + 1) * tx * ty
  );
}

export function terrainHeight(s: SiteScene, x: number, y: number): number {
  return s.terrain.kind === 'flat' ? s.terrain.z : sampleGrid(s.terrain, x, y);
}

/** The material at (x, y, z) from the composition, ignoring features and structures. */
export function materialAt(s: SiteScene, x: number, y: number, z: number): string {
  const surf = terrainHeight(s, x, y);
  if (z > surf) return 'air';
  const depth = surf - z;
  let t = 0;
  for (const c of s.cover) {
    if (depth >= t && depth < t + c.thickness) return c.material;
    t += c.thickness;
  }
  let m = s.strata[0].material;
  for (const st of s.strata) {
    const top = st.top === 'surface' ? surf : Math.min(sampleGrid(st.top, x, y), surf);
    if (z <= top) m = st.material;
  }
  return m;
}
