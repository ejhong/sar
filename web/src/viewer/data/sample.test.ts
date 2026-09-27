import { readFileSync, readdirSync } from 'node:fs';
import { join } from 'node:path';
import { describe, expect, it } from 'vitest';
import { materialAt, sampleGrid, terrainHeight } from './sample';
import type { SiteScene } from './types';

const DATA = join(__dirname, '../../../public/data/sites');
const scene = (id: string): SiteScene => JSON.parse(readFileSync(join(DATA, id, 'scene.json'), 'utf8'));

describe('exported sites', () => {
  const index = JSON.parse(readFileSync(join(DATA, 'index.json'), 'utf8')).sites as { id: string }[];

  it('lists every scene that was written', () => {
    const dirs = readdirSync(DATA, { withFileTypes: true }).filter((d) => d.isDirectory()).map((d) => d.name).sort();
    expect(index.map((s) => s.id).sort()).toEqual(dirs);
  });

  it('gives every feature a source and every material a colour', () => {
    for (const { id } of index) {
      const s = scene(id);
      for (const f of s.features) expect(f.source.length).toBeGreaterThan(0);
      for (const m of Object.values(s.materials)) expect(m.colour).toMatch(/^#[0-9a-f]{6}$/i);
    }
  });
});

describe('composition sampling', () => {
  it('reads the one-chamber bench as bare limestone under air', () => {
    const s = scene('bench-void');
    expect(terrainHeight(s, 0, 0)).toBe(0);
    expect(materialAt(s, 0, 0, 1)).toBe('air');
    expect(materialAt(s, 0, 0, -0.5)).toBe('limestone-mokattam');
    expect(materialAt(s, 0, 0, -20)).toBe('limestone-mokattam');
  });

  it('interpolates grids bilinearly and clamps at the edges', () => {
    const g = { x0: 0, y0: 0, dx: 1, nx: 2, ny: 2, z: [0, 1, 2, 3] };
    expect(sampleGrid(g, 0.5, 0.5)).toBeCloseTo(1.5);
    expect(sampleGrid(g, -5, -5)).toBe(0);
    expect(sampleGrid(g, 9, 9)).toBe(3);
  });

  it('puts the bedrock of the shaft field under its fractured layer', () => {
    const s = scene('bench-shafts');
    expect(materialAt(s, 0, 0, -30)).toBe('limestone-mokattam');
    expect(materialAt(s, 0, 0, -4)).toBe('limestone-fractured');
  });
});
