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
import type { RadarInfo, RadarSensors } from '../data/types';
import { dataUrl } from '../data/load';

export type Reading = 'complex' | 'magnitude';

/** Seconds of wall time per look while the sensors play. */
const LOOK_SECONDS = 0.55;
/** Sensor stalk height (m) per µm/s, and the clip for wild readings. */
const STALK_M_PER_UM_S = 5 / 1400;
const CLIP_UM_S = 1500;
/** Seconds for the satellite to glide once along its track when there are no sensors to pace it. */
const PASS_SECONDS = 9;

/**
 * The satellite over a site, in site coordinates: its spotlight beam and the image it made, laid on the ground. On a
 * bench the image is the simulated one on flat ground, with a grid of the image's virtual sensors, each carrying two
 * stalks: the true motion under it (ochre) and what the image reports there (cinnabar). On a real site the image is
 * the real product resampled onto the terrain, and the rim of the block marks how far down the wave can reach. The
 * satellite is not to scale in the sky: the real one is hundreds of kilometres away.
 */
export class Radar {
  readonly group = new Group();
  private sat = new Group();
  private beam: LineSegments;
  private track: LineSegments;
  private sensors?: {
    info: RadarSensors;
    n: number;
    truthLines: LineSegments;
    readLines: LineSegments;
    truthTips: Points;
    readTips: Points;
  };
  private image: Mesh;
  private mats: { mat: { color: Color }; role: 'sensor' | 'radar' | 'frame' }[] = [];
  private t = 0;
  playing = true;
  reading: Reading = 'complex';
  private readonly corners: Vector3[];
  private readonly satBase: Vector3;
  private readonly along: Vector3;
  private readonly span: number;

  constructor(
    readonly info: RadarInfo,
    siteId: string,
    size: number,
    private height: (x: number, y: number) => number,
    clip: Plane,
    theme: SceneTheme,
  ) {
    const a = info.acquisition;
    const im = info.image;
    const tex: Texture = new TextureLoader().load(dataUrl(`sites/${siteId}/${im.file}`));
    tex.colorSpace = SRGBColorSpace;
    let centre: Vector3;
    if (im.kind === 'ortho') {
      // the real image, resampled onto the ground: a mesh that follows the terrain
      const [x0, x1] = im.extent.x;
      const [y0, y1] = im.extent.y;
      centre = new Vector3((x0 + x1) / 2, (y0 + y1) / 2, height((x0 + x1) / 2, (y0 + y1) / 2));
      this.corners = [
        [x0, y0],
        [x1, y0],
        [x1, y1],
        [x0, y1],
      ].map(([x, y]) => new Vector3(x, y, height(x, y)));
      const nx = 161;
      const ny = Math.max(2, Math.round((nx - 1) * ((y1 - y0) / (x1 - x0))) + 1);
      const pos = new Float32Array(nx * ny * 3);
      const uv = new Float32Array(nx * ny * 2);
      for (let j = 0; j < ny; j++)
        for (let i = 0; i < nx; i++) {
          const k = i + nx * j;
          const x = x0 + ((x1 - x0) * i) / (nx - 1);
          const y = y0 + ((y1 - y0) * j) / (ny - 1);
          pos.set([x, y, height(x, y)], 3 * k);
          uv.set([i / (nx - 1), j / (ny - 1)], 2 * k);
        }
      const idx: number[] = [];
      for (let j = 0; j < ny - 1; j++)
        for (let i = 0; i < nx - 1; i++) {
          const p = i + nx * j;
          idx.push(p, p + 1, p + nx + 1, p, p + nx + 1, p + nx);
        }
      const geo = new BufferGeometry();
      geo.setAttribute('position', new BufferAttribute(pos, 3));
      geo.setAttribute('uv', new BufferAttribute(uv, 2));
      geo.setIndex(idx);
      this.image = new Mesh(
        geo,
        new MeshBasicMaterial({ map: tex, side: DoubleSide, clippingPlanes: [clip], polygonOffset: true, polygonOffsetFactor: -2, polygonOffsetUnits: -4 }),
      );
    } else {
      // the simulated image on flat ground, rotated to the pass
      const rot = (im.rotation_deg * Math.PI) / 180;
      const gr = new Vector3(Math.cos(rot), Math.sin(rot), 0);
      const at = new Vector3(-Math.sin(rot), Math.cos(rot), 0);
      const z = height(im.centre[0], im.centre[1]);
      centre = new Vector3(im.centre[0], im.centre[1], z);
      this.corners = [
        [-1, -1],
        [1, -1],
        [1, 1],
        [-1, 1],
      ].map(([u, v]) => centre.clone().addScaledVector(gr, (u * im.width_m) / 2).addScaledVector(at, (v * im.height_m) / 2));
      this.image = new Mesh(
        new PlaneGeometry(im.width_m, im.height_m),
        new MeshBasicMaterial({ map: tex, side: DoubleSide, transparent: true, opacity: 0.92, clippingPlanes: [clip], depthWrite: false }),
      );
      this.image.rotation.z = rot;
      this.image.position.set(centre.x, centre.y, z + 0.08);
    }
    this.image.renderOrder = 1;

    // the satellite: far up the line of sight (not to scale), moving along its track
    const los = new Vector3(...a.los_enu).normalize();
    this.satBase = centre.clone().addScaledVector(los, size * 0.8);
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
    this.beam.geometry.setAttribute('position', new BufferAttribute(new Float32Array(10 * 3), 3));
    const trackMat = new LineDashedMaterial({ color: theme.frame, dashSize: size / 40, gapSize: size / 60, transparent: true, opacity: 0.8 });
    this.mats.push({ mat: trackMat, role: 'frame' });
    const t0 = this.satBase.clone().addScaledVector(this.along, -this.span * 1.3);
    const t1 = this.satBase.clone().addScaledVector(this.along, this.span * 1.3);
    this.track = new LineSegments(new BufferGeometry().setFromPoints([t0, t1]), trackMat);
    this.track.computeLineDistances();
    this.group.add(this.image, this.sat, this.beam, this.track);

    // the virtual sensors, on a bench
    if (info.sensors) {
      const n = info.sensors.east.length;
      const mkLines = (color: number, role: 'sensor' | 'radar', opacity: number) => {
        const mat = new LineBasicMaterial({ color, transparent: true, opacity, clippingPlanes: [clip] });
        this.mats.push({ mat, role });
        const g = new BufferGeometry();
        g.setAttribute('position', new BufferAttribute(new Float32Array(n * 6), 3));
        const l = new LineSegments(g, mat);
        l.renderOrder = 6;
        return l;
      };
      const mkTips = (color: number, role: 'sensor' | 'radar', px: number) => {
        const mat = new PointsMaterial({ color, size: px, sizeAttenuation: false, clippingPlanes: [clip] });
        this.mats.push({ mat, role });
        const g = new BufferGeometry();
        g.setAttribute('position', new BufferAttribute(new Float32Array(n * 3), 3));
        const p = new Points(g, mat);
        p.renderOrder = 7;
        return p;
      };
      this.sensors = {
        info: info.sensors,
        n,
        truthLines: mkLines(theme.sensor, 'sensor', 0.9),
        readLines: mkLines(theme.radar, 'radar', 0.95),
        truthTips: mkTips(theme.sensor, 'sensor', 4),
        readTips: mkTips(theme.radar, 'radar', 4),
      };
      const sv = this.sensors;
      this.group.add(sv.truthLines, sv.readLines, sv.truthTips, sv.readTips);
    }
    this.setTime(0);
  }

  get hasSensors(): boolean {
    return !!this.sensors;
  }

  get looks(): number {
    return this.sensors?.info.looks_s.length ?? 0;
  }

  /** The current look (fractional) and its time in the dwell. */
  get look(): { index: number; time_s: number } {
    const sv = this.sensors;
    if (!sv) return { index: 0, time_s: 0 };
    const L = this.looks;
    const i = (this.t / LOOK_SECONDS) % L;
    const ls = sv.info.looks_s;
    const i0 = Math.floor(i);
    const w = i - i0;
    return { index: i, time_s: ls[i0] * (1 - w) + ls[(i0 + 1) % L] * w };
  }

  setImageVisible(on: boolean) {
    this.image.visible = on;
  }

  tick(dt: number) {
    if (!this.playing) return;
    this.setTime(this.t + dt);
  }

  setTime(t: number) {
    this.t = t;
    const sv = this.sensors;
    const loop = sv ? LOOK_SECONDS * this.looks : PASS_SECONDS;
    if (sv) {
      const L = this.looks;
      const f = (t / LOOK_SECONDS) % L;
      const i0 = Math.floor(f);
      const i1 = Math.min(i0 + 1, L - 1);
      const w = f - i0;
      const s = sv.info;
      const read = this.reading === 'complex' ? s.complex_um_s : s.magnitude_um_s;
      const off = Math.min(1.2, 60 / Math.sqrt(sv.n)) * 0.35; // the two stalks side by side
      const tl = sv.truthLines.geometry.getAttribute('position') as BufferAttribute;
      const rl = sv.readLines.geometry.getAttribute('position') as BufferAttribute;
      const tt = sv.truthTips.geometry.getAttribute('position') as BufferAttribute;
      const rt = sv.readTips.geometry.getAttribute('position') as BufferAttribute;
      const h = (v: number) => Math.max(-CLIP_UM_S, Math.min(CLIP_UM_S, v)) * STALK_M_PER_UM_S;
      for (let k = 0; k < sv.n; k++) {
        const x = s.east[k];
        const y = s.north[k];
        const z0 = this.height(x, y) + 0.15;
        const vt = s.truth_um_s[k][i0] * (1 - w) + s.truth_um_s[k][i1] * w;
        const vr = read[k][i0] * (1 - w) + read[k][i1] * w;
        tl.setXYZ(2 * k, x - off, y, z0);
        tl.setXYZ(2 * k + 1, x - off, y, z0 + h(vt));
        rl.setXYZ(2 * k, x + off, y, z0);
        rl.setXYZ(2 * k + 1, x + off, y, z0 + h(vr));
        tt.setXYZ(k, x - off, y, z0 + h(vt));
        rt.setXYZ(k, x + off, y, z0 + h(vr));
      }
      for (const a of [tl, rl, tt, rt]) a.needsUpdate = true;
      for (const o of [sv.truthLines, sv.readLines, sv.truthTips, sv.readTips]) o.geometry.computeBoundingSphere();
    }
    // the satellite glides along its track; the beam stays on the scene (spotlight)
    const phase = (t / loop) % 1;
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
