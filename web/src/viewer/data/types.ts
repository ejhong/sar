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
  /** The depth step, where it is finer than the horizontal one (the paper-style volumes, so their bands show). */
  spacing_z?: number;
  range: [number, number];
  run?: string;
  caption?: string;
  /** 'gated': the stricter reconstruction's fit scores, drawn in gold rather than the radar's cinnabar. */
  tint?: 'gated';
}

/** The gated reconstruction's depth scale on the pass it read (P2-34): one turn of its fit, the depth at which the scale
 * repeats, and how many passing positions stand near the surface and how many at the mirror of their depth. */
export interface DepthScale {
  turn_m: number;
  repeat_m: number;
  positions: number;
  shallow: number;
  mirror: number;
}

export interface RadarSensors {
  east: number[];
  north: number[];
  looks_s: number[];
  truth_um_s: number[][];
  complex_um_s: number[][];
  magnitude_um_s: number[][];
  test_wave: { wavelength_m: number; f_hz: number; v_m_s: number; azimuth_deg: number };
  gains: { complex: number; magnitude: number };
  run: string;
}

export interface RadarInfo {
  /** 'bench': a simulated image on flat ground; 'real': a real product resampled onto the site's terrain. */
  kind: 'bench' | 'real';
  acquisition: {
    name: string;
    heading_deg: number;
    incidence_deg: number;
    los_enu: [number, number, number];
    aperture_s: number;
    track_km: number;
    slant_range_km: number;
    along_track_en: [number, number];
    ground_range_en: [number, number];
    satellite?: string;
    date?: string;
  };
  image:
    | { kind?: 'plane'; file: string; centre: [number, number]; width_m: number; height_m: number; rotation_deg: number }
    | { kind: 'ortho'; file: string; extent: { x: [number, number]; y: [number, number] } };
  sensors?: RadarSensors;
  /** A bench's volumes with and without its chamber, or a real site's list of the method's volumes. */
  volumes: { with: string; without: string } | string[];
  reach_m?: number;
  /** Start with the draped image hidden (a deep site, where it would cover what hangs beneath it). */
  image_hidden?: boolean;
  /** The numbers quoted beside a real pass's volumes, and the sentence the exporter quotes them in. */
  stats?: {
    monument_vs_control_profile_corr?: number;
    patch_profile_corr_range?: [number, number];
    profile_corr_range?: [number, number];
    pillar_power_vs_energy_min: number;
    text?: string;
  };
  run?: string;
  /** The gated reconstruction over one monument: its volumes (real image and motionless copy, by support), where to look,
   * and how its scores sit inside the surveyed chambers against the same depths elsewhere. */
  gated?: {
    title: string;
    date: string;
    volumes: {
      id: string;
      /** What went in: the real image or its motionless copy or open plateau; on a bench, the chamber, none, its imprint
       * boosted `boost` times, or a random perturbation of the imprint's size. */
      case: 'real' | 'twin' | 'plateau' | 'with' | 'without' | 'boosted' | 'null';
      support: number;
      boost?: number;
      /** On a bench: what shook it, Giza's ambient motion (P2-20) or a vibrator beside it (P2-24). */
      shaking?: 'ambient' | 'vibrator';
      /** This picture's own run, where it differs from the study's. */
      note?: string;
      focus?: [number, number, number];
    }[];
    focus: [number, number, number];
    radius_m: number;
    chambers?: { real: [number | null, number | null]; twin: [number | null, number | null] };
    depth_scale?: DepthScale | null;
    note: string;
  };
  /** Runs made with the lab's processing command (katabasis.lab), each over its own area. */
  lab?: {
    name: string;
    title: string;
    /** The place the run is named by, as the other methods name it. */
    area?: string;
    /** The pass's year ('both' where two passes are compared), and how its lines were laid ('ew', 'ns', or 'both'). */
    pass?: string;
    lines?: 'ew' | 'ns' | 'both';
    volumes: { id: string; case: 'real' | 'twin'; support: number }[];
    focus: [number, number, number];
    radius_m: number;
    note: string;
    run?: string;
    depth_scale?: DepthScale | null;
  }[];
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
  radar?: RadarInfo;
}

export interface SiteIndexEntry {
  id: string;
  name: string;
  kind: 'test' | 'real';
  summary: string;
  extent: Extent;
  features: Record<string, number>;
  volumes: number;
  /** How many pictures (and, for the geophones, recorded waves) each instrument has on the site. */
  instruments?: { geophones: number; satellite: number };
}
