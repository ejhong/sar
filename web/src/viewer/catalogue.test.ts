import { readFileSync } from 'node:fs';
import { join } from 'node:path';
import { describe, expect, it } from 'vitest';
import { areaOf, defaultChoice, dimensions, findChoice, geophoneMethods, satelliteMethods, step, viewSet } from './catalogue';
import type { SiteScene } from './data/types';

const DATA = join(__dirname, '../../public/data/sites');
const scene = (id: string): SiteScene => JSON.parse(readFileSync(join(DATA, id, 'scene.json'), 'utf8'));

describe('the lab’s catalogue', () => {
  it('reads places out of titles', () => {
    expect(areaOf('What the gated reconstruction computes at Khufu')).toBe('Khufu');
    expect(areaOf('The gated reconstruction across Khafre, from the 2022 image')).toBe('Khafre');
    expect(areaOf('What the gated reconstruction computes over the bench')).toBeUndefined();
  });

  it('puts every volume of every site under exactly one instrument', () => {
    for (const id of ['bench-void', 'bench-khafre-claim', 'giza', 'sacsayhuaman']) {
      const s = scene(id);
      const ids = [...geophoneMethods(s), ...satelliteMethods(s)].flatMap((m) => m.choices.filter((c) => c.kind === 'volume').map((c) => c.id));
      expect(ids.sort()).toEqual(s.volumes.map((v) => v.id).sort());
    }
  });

  it('sets Sacsayhuamán’s walls beside its houses and fields, the fields as the control, every method under the same names', () => {
    const s = scene('sacsayhuaman');
    const sat = satelliteMethods(s);
    const paper = sat.find((m) => m.key === 'paper')!;
    expect(paper.choices.map((c) => c.area)).toEqual(['Zigzag walls', 'Rodadero outcrop', 'Fields north of the site', 'San Cristóbal houses', 'San Blas houses', 'City grid houses']);
    expect(paper.choices.filter((c) => c.control).map((c) => c.area)).toEqual(['Fields north of the site']);
    expect(paper.choices[0].stats).toMatch(/walls, the houses and the fields/);
    const gated = sat.find((m) => m.key === 'gated');
    for (const c of gated?.choices ?? []) expect(paper.choices.map((p) => p.area)).toContain(c.area);
  });

  it('gives the one-chamber bench its waves, five geophone methods and the satellite’s picture with its control', () => {
    const s = scene('bench-void');
    const g = geophoneMethods(s);
    expect(g[0].key).toBe('waves');
    expect(g.filter((m) => m.dot === 'recovered')).toHaveLength(5);
    const fwi = g.find((m) => m.name === 'Full-waveform inversion')!;
    expect(fwi.choices.map((c) => c.input).sort()).toEqual(['P speed', 'S speed']);
    const paper = satelliteMethods(s).find((m) => m.key === 'paper')!;
    expect(paper.choices.map((c) => [c.input, c.control])).toEqual([
      ['With the chamber', false],
      ['Without it', true],
    ]);
    expect(defaultChoice(g)?.id).toBe('tt-crosshole');
  });

  it('opens Giza’s satellite on the gated reconstruction’s real image, and steps between its pictures keeping the rest', () => {
    const s = scene('giza');
    const sat = satelliteMethods(s);
    const first = defaultChoice(sat)!;
    expect(first.input).toBe('Real image');
    expect(first.support).toBe(1);
    const gated = findChoice(sat, first.id)!.method;
    const dims = dimensions(gated).map((d) => d.key);
    expect(dims).toContain('input');
    expect(dims).toContain('support');
    const copy = step(gated, first, 'input', 'Motionless copy');
    expect([copy.input, copy.support, copy.control]).toEqual(['Motionless copy', 1, true]);
    const four = step(gated, copy, 'support', 4);
    expect([four.input, four.support]).toEqual(['Motionless copy', 4]);
    const paper = sat.find((m) => m.key === 'paper')!;
    expect(paper.choices.find((c) => c.control)?.area).toBe('Open plateau');
    // open plateau is a place of its own, a control fed the real image; every place is drawn with its like
    const plateau = gated.choices.find((c) => c.area === 'Open plateau')!;
    expect([plateau.input, plateau.control]).toEqual(['Real image', true]);
    expect(dimensions(gated).find((d) => d.key === 'input')?.values).not.toContain('Open plateau');
    const together = viewSet(gated, first).map((c) => c.area);
    expect(together).toContain('Khufu');
    expect(together).toContain('Open plateau');
  });

  it('shows the claimed underworld’s picture under the satellite, not the geophones', () => {
    const s = scene('bench-khafre-claim');
    expect(geophoneMethods(s)).toHaveLength(0);
    const sat = satelliteMethods(s);
    expect(sat).toHaveLength(1);
    expect(sat[0].choices).toHaveLength(1);
  });
});
