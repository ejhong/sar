/** The shapes of web/public/data/sites/** as written by sim/katabasis/export. */

export interface GridBlock {
  x0: number;
  y0: number;
  dx: number;
  nx: number;
  ny: number;
  /** Heights, x fastest (index = ix + nx * iy). */
  z: number[];
}

export type Shape =
  | { type: 'box'; centre: [number, number, number]; size: [number, number, number]; yaw_deg?: number; pitch_deg?: number }
  | { type: 'cylinder'; centre: [number, number, number]; radius: number; height: number }
  | { type: 'sphere'; centre: [number, number, number]; radius: number }
  | { type: 'prism'; polygon: [number, number][]; bottom: number; top: number }
  | { type: 'pyramid'; centre: [number, number, number]; base: number; height: number; yaw_deg?: number };

export type FeatureStatus = 'truth' | 'surveyed' | 'claimed' | 'representative';

export interface Feature {
  id: string;
  name: string;
  kind: string;
  status: FeatureStatus;
  fill: string;
  shape: Shape;
  source: string;
  note?: string;
  group?: string;
  placement?: string;
}

export interface Structure {
  id: string;
  name: string;
  material: string;
  shape: Shape;
  source: string;
}

export interface MaterialSummary {
  id: string;
  name: string;
  family: string;
  colour: string;
  vp: { value: number; range: [number, number] | null; status: string; source: string | null };
  vs: { value: number; range: [number, number] | null; status: string; source: string | null };
  rho: { value: number; range: [number, number] | null; status: string; source: string | null };
}

export interface Stratum {
  name: string;
  material: string;
  top: 'surface' | GridBlock;
  heterogeneous: boolean;
}

export interface SurveyInfo {
  id: string;
  label: string;
  description: string;
  sources: [number, number, number][];
  receivers: [number, number, number][];
  boreholes?: [number, number][];
}

export interface VolumeInfo {
  id: string;
  label: string;
  method: string;
  survey?: string;
  status: 'truth' | 'tomogram' | 'radar';
  quantity: string;
  units: string;
  file: string;
  shape: [number, number, number];
  origin: [number, number, number];
  spacing: number;
  range: [number, number];
  run?: string;
  caption?: string;
}

export interface Extent {
  x: [number, number];
  y: [number, number];
  z: [number, number];
}

export interface SiteScene {
  id: string;
  name: string;
  kind: 'test' | 'real';
  summary: string;
  frame: Record<string, unknown>;
  extent: Extent;
  terrain: { kind: 'flat'; z: number } | ({ kind: 'grid'; source: string } & GridBlock);
  cover: { material: string; thickness: number }[];
  strata: Stratum[];
  water_table: number | null;
  materials: Record<string, MaterialSummary>;
  structures: Structure[];
  features: Feature[];
  landmarks: { id: string; name: string; position: [number, number, number]; source: string }[];
  notes: string[];
  sources: Record<string, string>;
  volumes: VolumeInfo[];
  surveys?: SurveyInfo[];
  wavefields?: import('../scene/Wavefield').WavefieldInfo[];
}

export interface SiteIndexEntry {
  id: string;
  name: string;
  kind: 'test' | 'real';
  summary: string;
  extent: Extent;
  features: Record<string, number>;
  volumes: number;
}
