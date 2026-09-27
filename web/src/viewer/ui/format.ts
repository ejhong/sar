import { terrainHeight } from '../data/sample';
import type { Feature, SiteScene } from '../data/types';

export const fmtM = (m: number) => (m >= 1000 ? `${(m / 1000).toFixed(m % 1000 ? 2 : 0)} km` : `${Math.round(m)} m`);

/** Depth of a feature's top and bottom below the ground above its centre. */
export function featureDepth(s: SiteScene, f: Feature): string {
  const sh = f.shape;
  let cx: number, cy: number, top: number, bottom: number;
  if (sh.type === 'prism') {
    cx = sh.polygon.reduce((a, p) => a + p[0], 0) / sh.polygon.length;
    cy = sh.polygon.reduce((a, p) => a + p[1], 0) / sh.polygon.length;
    top = sh.top;
    bottom = sh.bottom;
  } else {
    [cx, cy] = sh.centre;
    const half =
      sh.type === 'box'
        ? verticalHalf(sh.size, sh.pitch_deg ?? 0)
        : sh.type === 'cylinder'
          ? sh.height / 2
          : sh.type === 'sphere'
            ? sh.radius
            : sh.height / 2;
    top = sh.centre[2] + half;
    bottom = sh.centre[2] - half;
  }
  const g = terrainHeight(s, cx, cy);
  const dTop = g - top;
  const dBot = g - bottom;
  if (dBot <= 0) return `${fmtM(-dBot)} above ground`;
  if (dTop <= 0) return `to ${fmtM(dBot)} deep`;
  return `${fmtM(dTop)} to ${fmtM(dBot)} deep`;
}

function verticalHalf(size: [number, number, number], pitchDeg: number): number {
  const p = (pitchDeg * Math.PI) / 180;
  return (Math.abs(size[1] * Math.sin(p)) + Math.abs(size[2] * Math.cos(p))) / 2;
}
