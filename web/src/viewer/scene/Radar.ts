import {
  BoxGeometry,
  BufferAttribute,
  BufferGeometry,
  Color,
  ConeGeometry,
  DoubleSide,
  Group,
  LineBasicMaterial,
  LineDashedMaterial,
  LineSegments,
  Matrix4,
  Mesh,
  MeshBasicMaterial,
  MeshLambertMaterial,
  PlaneGeometry,
  Points,
  PointsMaterial,
  SRGBColorSpace,
  TextureLoader,
  Vector3,
  type Plane,
  type Texture,
} from 'three';
import type { SceneTheme } from '../engine/theme';
import type { RadarInfo } from '../data/types';
import { dataUrl } from '../data/load';

export type Reading = 'complex' | 'magnitude';

/** Seconds of wall time per look while the sensors play. */
const LOOK_SECONDS = 0.55;
/** Sensor stalk height (m) per µm/s, and the clip for wild readings. */
const STALK_M_PER_UM_S = 5 / 1400;
const CLIP_UM_S = 1500;

/**
 * The satellite over a site, in site coordinates: its spotlight beam, the
 * image it made draped on the ground, and a grid of the image's virtual
 * sensors. Each sensor carries two stalks: the true motion under it (ochre)
 * and what the image reports there (cinnabar). Nothing here is to scale in
 * the sky: the real satellite is hundreds of kilometres away.
 */
export class Radar {
  readonly group = new Group();
  private sat = new Group();
  private beam: LineSegments;
  private track: LineSegments;
  private truthLines: LineSegments;
  private readLines: LineSegments;
  private truthTips: Points;
  private readTips: Points;
  private mats: { mat: { color: Color }; role: 'sensor' | 'radar' | 'frame' }[] = [];
  private t = 0;
  playing = true;
  reading: Reading = 'complex';
  private readonly n: number;
  private readonly corners: Vector3[];
  private readonly satBase: Vector3;
  private readonly along: Vector3;
  private readonly span: number;

  constructor(
    readonly info: RadarInfo,
    siteId: string,
    size: number,
    ground: number,
    clip: Plane,
    theme: SceneTheme,
  ) {
    const a = info.acquisition;
    const im = info.image;
    const rot = (im.rotation_deg * Math.PI) / 180;
    const gr = new Vector3(Math.cos(rot), Math.sin(rot), 0);
    const at = new Vector3(-Math.sin(rot), Math.cos(rot), 0);
    const c = new Vector3(im.centre[0], im.centre[1], ground);
    this.corners = [
      [-1, -1],
      [1, -1],
      [1, 1],
      [-1, 1],
    ].map(([u, v]) =>
      c.clone().addScaledVector(gr, (u * im.width_m) / 2).addScaledVector(at, (v * im.height_m) / 2),
    );

    // the image, draped just above the ground
    const tex: Texture = new TextureLoader().load(dataUrl(`sites/${siteId}/${im.file}`));
    tex.colorSpace = SRGBColorSpace;
    const plane = new Mesh(
      new PlaneGeometry(im.width_m, im.height_m),
      new MeshBasicMaterial({ map: tex, side: DoubleSide, transparent: true, opacity: 0.92, clippingPlanes: [clip], depthWrite: false }),
    );
    plane.rotation.z = rot;
    plane.position.set(c.x, c.y, ground + 0.08);
    plane.renderOrder = 1;

    // the satellite: far up the line of sight (not to scale), moving along its track
    const los = new Vector3(...a.los_enu).normalize();
    const d = size * 0.8;
    this.satBase = c.clone().addScaledVector(los, d);
    this.along = new Vector3(a.along_track_en[0], a.along_track_en[1], 0).normalize();
    this.span = size * 0.55;
    const s = size / 60;
    const body = new Mesh(new BoxGeometry(1.4 * s, 1.4 * s, 2.2 * s), new MeshLambertMaterial({ color: 0xd9dee3 }));
    const panelMat = new MeshLambertMaterial({ color: 0x3a5f8a, side: DoubleSide });
    const wing = new BoxGeometry(4.2 * s, 1.5 * s, 0.08 * s);
    const w1 = new Mesh(wing, panelMat);
    const w2 = new Mesh(wing, panelMat);
    w1.position.x = 2.9 * s;
    w2.position.x = -2.9 * s;
    const dish = new Mesh(new ConeGeometry(0.9 * s, 0.9 * s, 20, 1, true), new MeshLambertMaterial({ color: 0xf2f2f2, side: DoubleSide }));
    dish.rotation.x = Math.PI / 2;
    dish.position.z = -1.5 * s;
    this.sat.add(body, w1, w2, dish);
    // face the ground: the dish (-z) along -los, the wings along the track
    const zAxis = los.clone();
    const xAxis = this.along.clone();
    const yAxis = zAxis.clone().cross(xAxis).normalize();
    xAxis.copy(yAxis.clone().cross(zAxis)).normalize();
    this.sat.quaternion.setFromRotationMatrix(new Matrix4().makeBasis(xAxis, yAxis, zAxis));

    const beamMat = new LineBasicMaterial({ color: theme.radar, transparent: true, opacity: 0.55 });
    this.mats.push({ mat: beamMat, role: 'radar' });
    this.beam = new LineSegments(new BufferGeometry(), beamMat);
    this.beam.geometry.setAttribute('position', new BufferAttribute(new Float32Array(8 * 3 + 6), 3));
    const trackMat = new LineDashedMaterial({ color: theme.frame, dashSize: size / 40, gapSize: size / 60, transparent: true, opacity: 0.8 });
    this.mats.push({ mat: trackMat, role: 'frame' });
    const t0 = this.satBase.clone().addScaledVector(this.along, -this.span * 1.3);
    const t1 = this.satBase.clone().addScaledVector(this.along, this.span * 1.3);
    this.track = new LineSegments(new BufferGeometry().setFromPoints([t0, t1]), trackMat);
    this.track.computeLineDistances();

    // the virtual sensors
    const sv = info.sensors;
    this.n = sv.east.length;
    const mkLines = (color: number, role: 'sensor' | 'radar', opacity: number) => {
      const mat = new LineBasicMaterial({ color, transparent: true, opacity, clippingPlanes: [clip] });
      this.mats.push({ mat, role });
      const g = new BufferGeometry();
      g.setAttribute('position', new BufferAttribute(new Float32Array(this.n * 6), 3));
      const l = new LineSegments(g, mat);
      l.renderOrder = 6;
      return l;
    };
    const mkTips = (color: number, role: 'sensor' | 'radar', px: number) => {
      const mat = new PointsMaterial({ color, size: px, sizeAttenuation: false, clippingPlanes: [clip] });
      this.mats.push({ mat, role });
      const g = new BufferGeometry();
      g.setAttribute('position', new BufferAttribute(new Float32Array(this.n * 3), 3));
      const p = new Points(g, mat);
      p.renderOrder = 7;
      return p;
    };
    this.truthLines = mkLines(theme.sensor, 'sensor', 0.9);
    this.readLines = mkLines(theme.radar, 'radar', 0.95);
    this.truthTips = mkTips(theme.sensor, 'sensor', 4);
    this.readTips = mkTips(theme.radar, 'radar', 4);
    this.ground = ground;

    this.group.add(plane, this.sat, this.beam, this.track, this.truthLines, this.readLines, this.truthTips, this.readTips);
    this.setTime(0);
  }

  private ground: number;

  get looks(): number {
    return this.info.sensors.looks_s.length;
  }

  /** The current look (fractional) and its time in the dwell. */
  get look(): { index: number; time_s: number } {
    const L = this.looks;
    const f = this.t / LOOK_SECONDS;
    const i = f % L;
    const ls = this.info.sensors.looks_s;
    const i0 = Math.floor(i);
    const w = i - i0;
    return { index: i, time_s: ls[i0] * (1 - w) + ls[(i0 + 1) % L] * w };
  }

  tick(dt: number) {
    if (!this.playing) return;
    this.setTime(this.t + dt);
  }

  setTime(t: number) {
    this.t = t;
    const L = this.looks;
    const f = (t / LOOK_SECONDS) % L;
    const i0 = Math.floor(f);
    const i1 = Math.min(i0 + 1, L - 1);
    const w = f - i0;
    const sv = this.info.sensors;
    const read = this.reading === 'complex' ? sv.complex_um_s : sv.magnitude_um_s;
    const off = Math.min(1.2, 60 / Math.sqrt(this.n)) * 0.35; // the two stalks side by side
    const tl = this.truthLines.geometry.getAttribute('position') as BufferAttribute;
    const rl = this.readLines.geometry.getAttribute('position') as BufferAttribute;
    const tt = this.truthTips.geometry.getAttribute('position') as BufferAttribute;
    const rt = this.readTips.geometry.getAttribute('position') as BufferAttribute;
    const z0 = this.ground + 0.15;
    const h = (v: number) => Math.max(-CLIP_UM_S, Math.min(CLIP_UM_S, v)) * STALK_M_PER_UM_S;
    for (let k = 0; k < this.n; k++) {
      const x = sv.east[k];
      const y = sv.north[k];
      const vt = sv.truth_um_s[k][i0] * (1 - w) + sv.truth_um_s[k][i1] * w;
      const vr = read[k][i0] * (1 - w) + read[k][i1] * w;
      tl.setXYZ(2 * k, x - off, y, z0);
      tl.setXYZ(2 * k + 1, x - off, y, z0 + h(vt));
      rl.setXYZ(2 * k, x + off, y, z0);
      rl.setXYZ(2 * k + 1, x + off, y, z0 + h(vr));
      tt.setXYZ(k, x - off, y, z0 + h(vt));
      rt.setXYZ(k, x + off, y, z0 + h(vr));
    }
    for (const a of [tl, rl, tt, rt]) a.needsUpdate = true;
    for (const o of [this.truthLines, this.readLines, this.truthTips, this.readTips]) o.geometry.computeBoundingSphere();
    // the satellite glides along its track across the loop; the beam stays on the scene (spotlight)
    const phase = (t / (LOOK_SECONDS * L)) % 1;
    const p = this.satBase.clone().addScaledVector(this.along, (phase * 2 - 1) * this.span);
    this.sat.position.copy(p);
    const b = this.beam.geometry.getAttribute('position') as BufferAttribute;
    this.corners.forEach((cn, j) => {
      b.setXYZ(2 * j, p.x, p.y, p.z);
      b.setXYZ(2 * j + 1, cn.x, cn.y, cn.z + 0.1);
    });
    const mid = this.corners.reduce((acc, v) => acc.add(v), new Vector3()).multiplyScalar(0.25);
    b.setXYZ(8, p.x, p.y, p.z);
    b.setXYZ(9, mid.x, mid.y, mid.z + 0.1);
    b.needsUpdate = true;
    this.beam.geometry.computeBoundingSphere();
  }

  setReading(r: Reading) {
    this.reading = r;
    this.setTime(this.t);
  }

  setTheme(t: SceneTheme) {
    for (const { mat, role } of this.mats) mat.color.set(role === 'sensor' ? t.sensor : role === 'radar' ? t.radar : t.frame);
  }

  dispose() {
    this.group.traverse((o) => {
      const m = o as Mesh;
      m.geometry?.dispose?.();
      const mat = m.material as { dispose?: () => void; map?: Texture } | undefined;
      mat?.map?.dispose();
      mat?.dispose?.();
    });
  }
}
