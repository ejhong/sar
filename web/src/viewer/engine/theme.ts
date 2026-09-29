/**
 * The underworld's two grounds. Night (deep water) is the default: a lamp-lit
 * slate-blue in which the only things that glow are the ones that carry
 * meaning. Day (mist) reads like a printed block diagram. Keep the roles in
 * step with src/styles/tokens.css.
 */
export type ThemeName = 'night' | 'day';

export interface SceneTheme {
  name: ThemeName;
  background: number;
  backgroundInner: string;
  fog: number;
  glow: number; // bloom strength; 0 disables bloom
  // roles
  void: number; // lapis: voids and chambers
  recovered: number; // verdigris: what an instrument recovered
  sensor: number; // ochre: sensors and sources
  claimed: number; // porphyry: claimed structures
  radar: number; // cinnabar: the radar
  gated: number; // gold: the gated reconstruction's fit scores
  // ground
  terrainLow: number;
  terrainHigh: number;
  terrainLight: number;
  contour: number;
  wallDim: number; // multiplier for material colours on the block walls (0..1)
  interface: string; // CSS colour for stratum interfaces on walls
  structure: number;
  structureEdge: number;
  frame: number; // block edges
  text: string;
}

export const THEMES: Record<ThemeName, SceneTheme> = {
  night: {
    name: 'night',
    background: 0x0a121a,
    backgroundInner: '#1d3043',
    fog: 0x0f1822,
    glow: 0.55,
    void: 0x82b3ff,
    recovered: 0x74d8bc,
    sensor: 0xe8b95c,
    claimed: 0xc0a6da,
    radar: 0xf0906f,
    gated: 0xf2c14e,
    terrainLow: 0x223140,
    terrainHigh: 0x56636f,
    terrainLight: 0xe8f0f7,
    contour: 0x7f93a6,
    wallDim: 0.3,
    interface: 'rgba(223,230,236,0.55)',
    structure: 0xdfe6ec,
    structureEdge: 0xe8eef3,
    frame: 0x5d7082,
    text: '#dfe6ec',
  },
  day: {
    name: 'day',
    background: 0xdbe3ea,
    backgroundInner: '#f8fafc',
    fog: 0xe8eef2,
    glow: 0,
    void: 0x2c5aa0,
    recovered: 0x16876e,
    sensor: 0xa8741f,
    claimed: 0x7a55a6,
    radar: 0xaa4b33,
    gated: 0x9a6b12,
    terrainLow: 0xc9c2b2,
    terrainHigh: 0xefe9dc,
    terrainLight: 0xffffff,
    contour: 0x8f8a7e,
    wallDim: 1,
    interface: 'rgba(28,40,51,0.55)',
    structure: 0xebe6dc,
    structureEdge: 0x55636f,
    frame: 0x6b7d8e,
    text: '#1c2833',
  },
};
