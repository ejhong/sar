import type { SceneTheme } from '../engine/theme';
import { materialAt, terrainHeight } from '../data/sample';
import type { SiteScene } from '../data/types';

/**
 * A vertical section through the composition, drawn to a canvas the way a
 * geological block diagram is printed: each material in its colour with its
 * conventional pattern (limestone as bedded blocks, sand as stipple, marl as
 * dashes, igneous rock as crosses, masonry as large blocks), interfaces
 * inked, the water table dashed. Above the ground the canvas is clear.
 */
export interface SectionSpec {
  /** Start and end of the section line in site coordinates (x, y). */
  from: [number, number];
  to: [number, number];
  zMin: number;
  zMax: number;
  /** Pixels per metre along both axes (same scale, so patterns are square). */
  ppm: number;
}

const hex = (c: string) => {
  const v = parseInt(c.slice(1), 16);
  return [(v >> 16) & 255, (v >> 8) & 255, v & 255];
};

function tint(c: string, k: number, towards: number[]): string {
  const [r, g, b] = hex(c);
  const m = (a: number, t: number) => Math.round(a * k + t * (1 - k));
  return `rgb(${m(r, towards[0])},${m(g, towards[1])},${m(b, towards[2])})`;
}

export function drawSection(s: SiteScene, spec: SectionSpec, theme: SceneTheme, unit: number): HTMLCanvasElement {
  const len = Math.hypot(spec.to[0] - spec.from[0], spec.to[1] - spec.from[1]);
  const W = Math.max(8, Math.min(2048, Math.round(len * spec.ppm)));
  const H = Math.max(8, Math.min(2048, Math.round((spec.zMax - spec.zMin) * spec.ppm)));
  const cv = document.createElement('canvas');
  cv.width = W;
  cv.height = H;
  const g = cv.getContext('2d')!;
  const px2z = (py: number) => spec.zMax - (py + 0.5) / spec.ppm;
  const at = (px: number): [number, number] => {
    const t = (px + 0.5) / W;
    return [spec.from[0] + (spec.to[0] - spec.from[0]) * t, spec.from[1] + (spec.to[1] - spec.from[1]) * t];
  };
  const night = theme.name === 'night';
  const ground = night ? hex('#0f1214') : hex('#f4efe5');
  const colourOf = (m: string) => tint(s.materials[m]?.colour ?? '#888888', theme.wallDim, ground);
  const ink = night ? 'rgba(233,229,219,0.20)' : 'rgba(45,41,36,0.26)';

  // Column by column: the material runs, then fill each run and pattern it.
  const step = Math.max(1, Math.round(W / 512));
  const runsByCol: { px: number; runs: { m: string; y0: number; y1: number }[]; surfPy: number }[] = [];
  for (let px = 0; px < W; px += step) {
    const [x, y] = at(px);
    const surfPy = Math.max(0, (spec.zMax - terrainHeight(s, x, y)) * spec.ppm);
    const runs: { m: string; y0: number; y1: number }[] = [];
    let cur = '';
    const rowStep = Math.max(1, Math.round(H / 400));
    for (let py = Math.floor(surfPy); py < H; py += rowStep) {
      const m = materialAt(s, x, y, px2z(py));
      if (m !== cur) {
        if (runs.length) runs[runs.length - 1].y1 = py;
        runs.push({ m, y0: py, y1: H });
        cur = m;
      }
    }
    runsByCol.push({ px, runs, surfPy });
    for (const r of runs) {
      if (r.m === 'air') continue;
      g.fillStyle = colourOf(r.m);
      g.fillRect(px, r.y0, step, r.y1 - r.y0);
    }
  }

  // Patterns, clipped to each material's region.
  const u = Math.max(3, unit * spec.ppm);
  const families = new Set(Object.values(s.materials).map((m) => m.family));
  for (const fam of families) {
    const mats = Object.values(s.materials).filter((m) => m.family === fam).map((m) => m.id);
    g.save();
    g.beginPath();
    for (const col of runsByCol) for (const r of col.runs) if (mats.includes(r.m)) g.rect(col.px, r.y0, step, r.y1 - r.y0);
    g.clip();
    g.strokeStyle = ink;
    g.fillStyle = ink;
    g.lineWidth = Math.max(0.6, u / 14);
    pattern(g, fam, W, H, u);
    g.restore();
  }

  // Interfaces: where the material changes down a column.
  g.fillStyle = theme.interface;
  const lw = Math.max(1, Math.round(u / 10));
  for (const col of runsByCol) {
    for (let i = 1; i < col.runs.length; i++) g.fillRect(col.px, col.runs[i].y0 - lw / 2, step, lw);
    g.fillRect(col.px, col.surfPy - lw / 2, step, lw * 1.4);
  }

  // Depth ticks down the left edge, measured from the ground there.
  {
    const [x0, y0] = at(0);
    const ground = terrainHeight(s, x0, y0);
    const depth = ground - spec.zMin;
    const step = niceStep(depth / 5);
    g.fillStyle = night ? 'rgba(233,229,219,0.62)' : 'rgba(45,41,36,0.7)';
    g.strokeStyle = g.fillStyle;
    g.lineWidth = Math.max(1, lw * 0.8);
    const fs = Math.max(9, Math.min(26, u * 0.9));
    g.font = `500 ${fs}px ui-monospace, SFMono-Regular, Menlo, monospace`;
    g.textBaseline = 'middle';
    for (let d = step; d < depth - step * 0.3; d += step) {
      const py = (spec.zMax - (ground - d)) * spec.ppm;
      g.beginPath();
      g.moveTo(0, py);
      g.lineTo(u * 0.9, py);
      g.stroke();
      g.fillText(`${d} m`, u * 1.2, py);
    }
  }

  if (s.water_table !== null && s.water_table > spec.zMin && s.water_table < spec.zMax) {
    const py = (spec.zMax - s.water_table) * spec.ppm;
    g.save();
    g.beginPath();
    for (const col of runsByCol) if (col.surfPy < py) g.rect(col.px, py - 4, step, 8);
    g.clip();
    g.strokeStyle = night ? 'rgba(130,179,255,0.8)' : 'rgba(44,90,160,0.85)';
    g.lineWidth = lw;
    g.setLineDash([u * 0.9, u * 0.5]);
    g.beginPath();
    g.moveTo(0, py);
    g.lineTo(W, py);
    g.stroke();
    g.restore();
  }
  return cv;
}

function pattern(g: CanvasRenderingContext2D, family: string, W: number, H: number, u: number) {
  g.beginPath();
  if (family === 'limestone' || family === 'structure') {
    const ch = family === 'structure' ? u * 1.3 : u; // course height
    const bw = ch * 2.4;
    for (let r = 0, y = 0; y < H; r++, y += ch) {
      g.moveTo(0, y);
      g.lineTo(W, y);
      const off = (r % 2) * bw * 0.5;
      for (let x = -off; x < W; x += bw) {
        g.moveTo(x, y);
        g.lineTo(x, y + ch);
      }
    }
    g.stroke();
  } else if (family === 'unconsolidated') {
    const d = u * 0.55;
    const r = Math.max(0.6, u / 16);
    let k = 0;
    for (let y = d / 2; y < H; y += d) {
      for (let x = ((k++ % 2) * d) / 2; x < W; x += d) {
        const jx = Math.sin(x * 12.9898 + y * 78.233) * 0.3 * d;
        const jy = Math.cos(x * 4.1414 + y * 3.33) * 0.3 * d;
        g.moveTo(x + jx + r, y + jy);
        g.arc(x + jx, y + jy, r, 0, Math.PI * 2);
      }
    }
    g.fill();
  } else if (family === 'clastic') {
    const d = u * 0.7;
    let k = 0;
    for (let y = d / 2; y < H; y += d) {
      const off = (k++ % 2) * u;
      for (let x = -off; x < W; x += u * 2) {
        g.moveTo(x, y);
        g.lineTo(x + u * 1.2, y);
      }
    }
    g.stroke();
  } else if (family === 'igneous') {
    const d = u * 1.1;
    let k = 0;
    for (let y = d / 2; y < H; y += d) {
      for (let x = ((k++ % 2) * d) / 2; x < W; x += d) {
        const a = u * 0.22;
        g.moveTo(x - a, y);
        g.lineTo(x + a, y);
        g.moveTo(x, y - a);
        g.lineTo(x, y + a);
      }
    }
    g.stroke();
  } else if (family === 'fluid') {
    const d = u * 0.8;
    for (let y = d / 2; y < H; y += d) {
      g.moveTo(0, y);
      for (let x = 0; x < W; x += u / 4) g.lineTo(x, y + Math.sin((x / u) * Math.PI) * u * 0.12);
    }
    g.stroke();
  }
}

/** The block's four faces as section specs (south, east, north, west). */
export function faceSpecs(s: SiteScene, zMax: number, ppm: number) {
  const { x, y, z } = s.extent;
  const faces: { id: string; spec: SectionSpec }[] = [
    { id: 'south', spec: { from: [x[0], y[0]], to: [x[1], y[0]], zMin: z[0], zMax, ppm } },
    { id: 'east', spec: { from: [x[1], y[0]], to: [x[1], y[1]], zMin: z[0], zMax, ppm } },
    { id: 'north', spec: { from: [x[1], y[1]], to: [x[0], y[1]], zMin: z[0], zMax, ppm } },
    { id: 'west', spec: { from: [x[0], y[1]], to: [x[0], y[0]], zMin: z[0], zMax, ppm } },
  ];
  return faces;
}

function niceStep(raw: number): number {
  const p = Math.pow(10, Math.floor(Math.log10(raw)));
  const m = raw / p;
  return (m < 1.5 ? 1 : m < 3.5 ? 2 : m < 7.5 ? 5 : 10) * p;
}
