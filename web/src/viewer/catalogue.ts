import type { DepthScale, SiteScene, VolumeInfo } from './data/types';
import { GATED_CODE } from '../data/credits';

/**
 * What the lab offers on a site, arranged as it works: an instrument (geophones on the ground, or the satellite), the
 * methods that turn its records into a picture, and each picture's choices (where it was run, what went in, which of
 * them are controls). Built from the exported scene alone, so a new run appears without code here changing.
 */

export type Instrument = 'geophones' | 'satellite';

export interface Choice {
  /** The volume or wavefield shown. */
  id: string;
  kind: 'volume' | 'wave';
  /** Where the method was run, when it was run over several places. A method's pictures of every place with the same
   * pass, lines, input and support are shown together, one underworld. */
  area?: string;
  /** The radar pass it was made from, and how the reconstruction's lines were laid. */
  pass?: string;
  lines?: string;
  /** What went in: the real image, a motionless copy, the ground with the chamber, … */
  input: string;
  /** A control: a picture that should hold nothing, drawn beside the one that might. */
  control: boolean;
  /** The gated reconstruction's support: how many positions in a row must agree. */
  support?: number;
  /** One line about this picture. */
  sub: string;
  /** The run behind it, and what was run there. */
  note?: string;
  stats?: string;
  run?: string;
  /** Where to stand to see it, when it is small beside the site. */
  focus?: [number, number, number];
  radius_m?: number;
  /** A starting threshold for the volume's display (0..1 of its scale). */
  threshold: number;
  /** Order among a method's choices: the picture that might hold something first, its controls after. */
  rank?: number;
}

export interface Method {
  key: string;
  instrument: Instrument;
  name: string;
  /** The colour it is drawn in: the instrument's. */
  dot: 'recovered' | 'radar' | 'wave';
  /** What its voxels hold. */
  quantity: 'material property' | 'focused power' | 'fit score' | 'wave motion';
  /** What it is, in one line. */
  line: string;
  /** Why its picture looks the way it does. */
  why?: string;
  /** Code this lab runs but did not write, credited where its pictures are shown. */
  code?: { name: string; url: string; version: string };
  choices: Choice[];
}

const cap = (s: string) => s.charAt(0).toUpperCase() + s.slice(1);

/** Where a volume sits: the centre of its top, and how far to stand back to see it whole. */
function anchor(v: VolumeInfo): { focus: [number, number, number]; radius_m: number } {
  const w = (v.shape[0] - 1) * v.spacing;
  const h = (v.shape[1] - 1) * v.spacing;
  return { focus: [v.origin[0] + w / 2, v.origin[1] + h / 2, v.origin[2]], radius_m: 2.2 * Math.max(w, h) };
}
/** A caption's finding about the chamber first, its conditions after in brackets; the colour scale's note left out. */
const head = (caption?: string) => {
  const parts = (caption ?? '').split(' · ').filter(Boolean);
  const finding = parts.find((x) => /chamber/.test(x)) ?? parts[0] ?? '';
  const rest = parts.filter((x) => x !== finding && !/^brightest/.test(x));
  return rest.length ? `${finding} (${rest.join(', ')})` : finding;
};

/** The place a title names: "… at Khufu", "… across Khafre, from …". */
export function areaOf(title: string): string | undefined {
  const m = title.match(/\b(?:at|across|over|beneath|under)\s+(?:the\s+)?([A-Z][\w-]*(?:\s+[A-Z][\w-]*)*)/);
  return m?.[1];
}

const SUP = '⁰¹²³⁴⁵⁶⁷⁸⁹';
const power = (b: number) => `×10${String(Math.round(Math.log10(b))).replace(/\d/g, (d) => SUP[+d])}`;

/** Each input the gated reconstruction was fed: its name, whether it is a control, a line about it, and its order. */
function gatedInput(kind: string, boost?: number): [string, boolean, string, number] {
  switch (kind) {
    case 'real':
      return ['Real image', false, 'the real image, through the unchanged code', 0];
    case 'with':
      return ['With the chamber', false, 'a synthetic image of this ground, shaken as Giza shakes, the chamber’s imprint in it', 1];
    case 'without':
      return ['Without it', true, 'the same image with no chamber', 2];
    case 'null':
      return ['Random, same size', true, 'a random perturbation the size of the chamber’s imprint, at the same pixels', 3];
    case 'boosted':
      return [`Imprint ${power(boost ?? 1)}`, false, `the chamber’s imprint made ${(boost ?? 1).toLocaleString('en-US')} times stronger`, 4 + Math.log10(boost ?? 1) / 100];
    case 'twin':
      return ['Motionless copy', true, 'the same crop with nothing moving and nothing inside', 5];
    case 'plateau':
      return ['Open plateau', true, 'the same raster over open plateau, where no monument stands', 6];
    default:
      return [cap(kind), false, '', 7];
  }
}

/** Where the gated reconstruction's depths come from on a pass (P2-34): its turn and repeat, and the two levels they make. */
function scaleLine(ds?: DepthScale | null): string | undefined {
  if (!ds) return undefined;
  const n = (v: number) => v.toLocaleString('en-US');
  const head = `On this pass one turn of the fit spans ${ds.turn_m.toFixed(1)} m and the depth scale repeats every ${Math.round(ds.repeat_m)} m`;
  return ds.repeat_m < 300
    ? `${head}, inside the 300 m the method draws, so each fit is drawn twice: near the surface, and again near ${Math.round(ds.repeat_m)} m down where the scale folds back. ${n(ds.shallow)} of ${n(ds.positions)} positions have their best score in the shallow level, ${n(ds.mirror)} in the deep one: the deeper level is the shallow one again.`
    : `${head}, just past the 300 m the method draws, so its fits stand near the surface (${n(ds.shallow)} of ${n(ds.positions)} positions).`;
}

export const supportName = (p: number) => (p === 1 ? 'one position' : `${p} in a row`);

export function geophoneMethods(s: SiteScene): Method[] {
  const out: Method[] = [];
  const waves = s.wavefields ?? [];
  if (waves.length)
    out.push({
      key: 'waves',
      instrument: 'geophones',
      name: 'The waves they record',
      dot: 'wave',
      quantity: 'wave motion',
      line: waves.some((w) => w.id.startsWith('vibrator'))
        ? 'a hammer blow, a steady vibrator, and what the chamber sends back or changes'
        : 'a hammer blow spreading through the ground, and what the chamber sends back',
      choices: waves.map((w) => ({ id: w.id, kind: 'wave', input: w.label, control: false, sub: w.caption, threshold: 0 })),
    });
  const groups = new Map<string, VolumeInfo[]>();
  for (const v of s.volumes.filter((v) => v.status !== 'radar')) {
    const k = `${v.survey ?? ''}|${v.method}`;
    groups.set(k, [...(groups.get(k) ?? []), v]);
  }
  for (const [k, vs] of groups) {
    const several = vs.length > 1;
    out.push({
      key: `g:${k}`,
      instrument: 'geophones',
      name: several ? vs[0].label.split(', ')[0] : vs[0].label,
      dot: 'recovered',
      quantity: 'material property',
      line: several ? vs.map((v) => `${v.label.split(', ').slice(1).join(', ')}: ${head(v.caption)}`).join('; ') : head(vs[0].caption) || vs[0].method,
      choices: vs.map((v) => ({
        id: v.id,
        kind: 'volume',
        input: several ? cap(v.label.split(', ').slice(1).join(', ')) : v.label,
        control: false,
        sub: `${v.method}. ${cap(v.caption ?? '')}`,
        run: v.run,
        threshold: 0.18,
      })),
    });
  }
  return out;
}

export function satelliteMethods(s: SiteScene): Method[] {
  const out: Method[] = [];
  const r = s.radar;
  const vol = (id: string) => s.volumes.find((v) => v.id === id);
  // the pipeline as the 2022 paper describes it
  const paper: Choice[] = [];
  if (r && !Array.isArray(r.volumes)) {
    for (const [id, input, control] of [
      [r.volumes.with, 'With the chamber', false],
      [r.volumes.without, 'Without it', true],
    ] as const) {
      const v = vol(id);
      if (v) paper.push({ id, kind: 'volume', input, control, sub: cap(v.caption ?? ''), run: v.run, threshold: 0.6 });
    }
  } else {
    const ids = r && Array.isArray(r.volumes) ? r.volumes : s.volumes.filter((v) => v.status === 'radar' && v.tint !== 'gated').map((v) => v.id);
    for (const id of ids) {
      const v = vol(id);
      if (!v) continue;
      const where = v.label.replace(/^Satellite · the published method\s*/, '').replace(/^(at|over|on)\s+/, '');
      const control = /\(control\)/.test(where);
      const place = where.replace(/\s*\(control\)/, '').split(',')[0].replace(/\s+pyramid$/, '');
      const single = ids.length === 1;
      paper.push({
        ...(single ? {} : anchor(v)),
        id,
        kind: 'volume',
        area: single ? undefined : cap(place),
        input: single ? cap(place) : 'Real image',
        pass: !single && r?.acquisition.date ? `${r.acquisition.date.slice(0, 4)} pass` : undefined,
        control,
        sub: cap(v.caption ?? ''),
        run: v.run,
        // a site deep enough for the whole axis shows its thinner columns at a lower threshold
        threshold: s.extent.z[1] - s.extent.z[0] > 1000 ? 0.6 : 0.8,
      });
    }
  }
  const real = r?.kind !== 'bench';
  if (paper.length)
    out.push({
      key: 'paper',
      instrument: 'satellite',
      name: 'Paper-style pipeline',
      dot: 'radar',
      quantity: 'focused power',
      line: 'as the 2022 paper describes it: sub-aperture pairs registered and focused, no selection gates',
      why: !r
        ? 'Here the method’s depth axis is drawn whole, relabelled so that it repeats where the claim puts the bottoms of its shafts. Its pillars run the full depth, where the surface reading is noisiest, and a block sits at each repeat, where every steering phase coincides: faint as the paper states the method, bright with pairs taken close together in the band (P2-35). That is where the published pictures’ shafts and deep structure can come from. Over the Giza plateau the same kind of volume is cut at the block’s floor, well above the first repeat, and shares one brightness scale with open ground, so neither the repeat blocks nor one patch’s own stretch appear there.'
        : real
        ? 'Pillars: a pixel whose registration wanders is bright at every depth. Bands: along a pillar the power rises and falls once per step of the axis’s resolution. Blocks: at the surface and at each repeat depth every steering phase coincides. Open ground draws the same shapes.'
        : 'Pillars where a pixel’s registration wanders, bands at each step of the axis’s resolution, blocks where every steering phase coincides; none of it depends on what is below.',
      choices: paper.map((c) =>
        real && r?.stats
          ? {
              ...c,
              note: `The ${r.acquisition.date?.slice(0, 4) ?? ''} image through the pipeline as the 2022 paper describes it: 50 half-band pairs, no selection gates, focused power on a log scale, depth relabelled so it repeats at 648 m as the claim does, and smoothed for display.`,
              stats: r.stats.text,
            }
          : c,
      ),
    });
  // the stricter gated reconstruction: one study and any lab runs, each over its own area
  const gated: Choice[] = [];
  const fmt = (v: number | null | undefined) => (v === null || v === undefined ? '–' : v.toFixed(2));
  if (r?.gated) {
    const g = r.gated;
    const area = r.kind === 'bench' ? undefined : areaOf(g.title);
    const stats =
      [
        g.chambers
          ? `Inside the surveyed chambers and passages the real image scores ${fmt(g.chambers.real[0])} on average, against ${fmt(g.chambers.real[1])} at the same depths elsewhere; the motionless copy ${fmt(g.chambers.twin[0])} and ${fmt(g.chambers.twin[1])}.`
          : '',
        scaleLine(g.depth_scale) ?? '',
      ]
        .filter(Boolean)
        .join(' ') || undefined;
    for (const v of g.volumes) {
      const plateau = v.case === 'plateau';
      const [input, control, sub, rank] = plateau ? gatedInput('real') : gatedInput(v.case, v.boost);
      gated.push({
        rank,
        id: v.id,
        kind: 'volume',
        area: plateau ? 'Open plateau' : area,
        pass: r.kind === 'bench' ? (v.shaking === 'vibrator' ? 'A vibrator, 30 m away' : v.shaking ? 'Giza’s own shaking' : undefined) : `${g.date} pass`,
        lines: r.kind === 'bench' ? undefined : 'East–west lines',
        input,
        control: plateau || control,
        support: v.support,
        sub: plateau
          ? 'the same raster over open plateau, where no monument stands'
          : v.case === 'real'
            ? `the ${g.date} image, through the unchanged code`
            : v.shaking === 'vibrator'
              ? v.case === 'with'
                ? 'shaken by the vibrator at the force where the reflectors’ data hold the chamber; the chamber in it'
                : v.case === 'without'
                  ? 'the same shaking, no chamber'
                  : 'the same without the chamber, changed by random noise the size of the chamber’s own change'
              : sub,
        note: v.note ?? g.note,
        stats,
        run: vol(v.id)?.run,
        focus: r.kind === 'bench' ? undefined : (v.focus ?? g.focus),
        radius_m: r.kind === 'bench' ? undefined : g.radius_m,
        threshold: 0.08,
      });
    }
  }
  for (const lab of r?.lab ?? []) {
    for (const v of lab.volumes) {
      const [input, control, sub, rank] = gatedInput(v.case);
      gated.push({
        rank,
        id: v.id,
        kind: 'volume',
        area: lab.area ?? areaOf(lab.title) ?? cap(lab.name),
        pass: lab.pass === 'both' ? 'Both passes agree' : `${lab.pass ?? '2022'} pass`,
        lines: lab.lines === 'both' ? 'Both layouts agree' : lab.lines === 'ns' ? 'North–south lines' : 'East–west lines',
        input,
        control,
        support: v.support,
        sub,
        note: lab.note,
        stats: scaleLine(lab.depth_scale),
        run: lab.run,
        focus: lab.focus,
        radius_m: lab.radius_m,
        threshold: 0.08,
      });
    }
  }
  if (gated.length)
    out.push({
      key: 'gated',
      instrument: 'satellite',
      name: 'Gated reconstruction',
      dot: 'radar',
      quantity: 'fit score',
      line: 'a stricter version with selection gates and a depth fit, its public code run unchanged',
      code: { name: GATED_CODE.name, url: GATED_CODE.url, version: GATED_CODE.version },
      why:
        r?.kind === 'bench'
          ? 'A column wherever a position passed the gates, banded where the fit’s phase turns whole times across a window. The speckle decides where: the chamber’s imprint, even made a hundred million times stronger, moves the columns no nearer to it than random noise of its size does.'
          : 'A column hangs wherever a position passed the gates; along it the score peaks where the depth fit’s phase turns a whole number of times across a window, and again at the mirror of that depth, where the scale folds back. How many metres a turn spans, and where the scale folds, is set by the pass’s geometry, not by the ground: the depth is the frequency the gates chose, not a measured depth.',
      choices: ordered(gated),
    });
  return out;
}

/** Areas in the order they were run, then the picture before its controls, then support. */
function ordered(cs: Choice[]): Choice[] {
  // places in the order they were run, the controls' places last
  const isControl = (a: string) => cs.filter((c) => (c.area ?? '') === a).every((c) => c.control);
  const first = [...new Set(cs.map((c) => c.area ?? ''))];
  const areas = [...first.filter((a) => !isControl(a)), ...first.filter(isControl)];
  return [...cs].sort(
    (a, b) =>
      areas.indexOf(a.area ?? '') - areas.indexOf(b.area ?? '') || (a.rank ?? 9) - (b.rank ?? 9) || (a.support ?? 0) - (b.support ?? 0),
  );
}

export function methodsOf(s: SiteScene, instrument: Instrument): Method[] {
  return instrument === 'geophones' ? geophoneMethods(s) : satelliteMethods(s);
}

export function findChoice(methods: Method[], id: string | undefined): { method: Method; choice: Choice } | undefined {
  if (!id) return undefined;
  for (const method of methods) {
    const choice = method.choices.find((c) => c.id === id);
    if (choice) return { method, choice };
  }
  return undefined;
}

/** The instrument's first picture: a real image before its controls, a picture before the waves. */
export function defaultChoice(methods: Method[]): Choice | undefined {
  const pictures = methods.filter((m) => m.dot !== 'wave');
  const pool = pictures.length ? pictures : methods;
  const preferred = pool.find((m) => m.key === 'gated' && m.choices.some((c) => c.focus)) ?? pool[0];
  return preferred?.choices.find((c) => !c.control) ?? preferred?.choices[0];
}

/** The dimensions a method's choices vary over, each with its values in order of first appearance. Every value is
 * offered wherever it was run: choosing one the current place lacks moves to the places that have it. */
export type Dimension = 'area' | 'pass' | 'lines' | 'input' | 'support';

export function dimensions(m: Method): { key: Dimension; values: (string | number)[] }[] {
  const dims: { key: Dimension; values: (string | number)[] }[] = [];
  for (const key of ['area', 'pass', 'lines', 'input', 'support'] as const) {
    const values: (string | number)[] = [];
    for (const c of m.choices) {
      const v = c[key];
      if (v !== undefined && !values.includes(v)) values.push(v);
    }
    if (values.length > 1) dims.push({ key, values });
  }
  return dims;
}

/** Step from the current choice along one dimension, keeping the others where they are when such a picture exists. */
export function step(m: Method, current: Choice, key: Dimension, value: string | number): Choice {
  const pool = m.choices.filter((c) => c[key] === value);
  const score = (c: Choice) =>
    (c.area === current.area ? 16 : 0) +
    (c.pass === current.pass ? 8 : 0) +
    (c.lines === current.lines ? 4 : 0) +
    (c.input === current.input ? 2 : 0) +
    (c.support === current.support ? 1 : 0);
  return pool.reduce((best, c) => (score(c) > score(best) ? c : best), pool[0]);
}

/** Everything drawn with a choice: the method's pictures of every place made the same way. */
export function viewSet(m: Method, c: Choice): Choice[] {
  return m.choices.filter((x) => x.pass === c.pass && x.lines === c.lines && x.input === c.input && x.support === c.support);
}
