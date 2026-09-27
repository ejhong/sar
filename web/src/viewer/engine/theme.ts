/**
 * The underworld's two grounds. Night (basalt) is the default: the only
 * things that glow are the ones that carry meaning. Day (limestone) reads
 * like a printed block diagram. Keep the roles in step with
 * src/styles/tokens.css.
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
    background: 0x0b0d0f,
    backgroundInner: '#1b2024',
    fog: 0x0f1214,
    glow: 0.55,
    void: 0x82b3ff,
    recovered: 0x74d8bc,
    sensor: 0xe8b95c,
    claimed: 0xc0a6da,
    radar: 0xf0906f,
    terrainLow: 0x2a2c2b,
    terrainHigh: 0x5b5548,
    terrainLight: 0xf1e3c4,
    contour: 0x8d826c,
    wallDim: 0.26,
    interface: 'rgba(233,229,219,0.55)',
    structure: 0xd9ccb0,
    structureEdge: 0xe9dcc0,
    frame: 0x59625f,
    text: '#e9e5db',
  },
  day: {
    name: 'day',
    background: 0xece5d8,
    backgroundInner: '#faf7f0',
    fog: 0xf4efe5,
    glow: 0,
    void: 0x2c5aa0,
    recovered: 0x347f6f,
    sensor: 0xa8741f,
    claimed: 0x7f6396,
    radar: 0xaa4b33,
    terrainLow: 0xcbbd9f,
    terrainHigh: 0xefe6d2,
    terrainLight: 0xfffaf0,
    contour: 0x9d8f74,
    wallDim: 1,
    interface: 'rgba(45,41,36,0.55)',
    structure: 0xe6dcc6,
    structureEdge: 0x6b5f4e,
    frame: 0x8a7a64,
    text: '#2d2924',
  },
};
