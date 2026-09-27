import {
  AdditiveBlending,
  BoxGeometry,
  BufferAttribute,
  BufferGeometry,
  CanvasTexture,
  Color,
  ConeGeometry,
  CylinderGeometry,
  DoubleSide,
  EdgesGeometry,
  FrontSide,
  ExtrudeGeometry,
  Group,
  LineBasicMaterial,
  LineDashedMaterial,
  LineSegments,
  Matrix4,
  Mesh,
  MeshBasicMaterial,
  MeshLambertMaterial,
  NormalBlending,
  Plane,
  Shape as ThreeShape,
  SphereGeometry,
  SRGBColorSpace,
  ShaderMaterial,
  Vector2,
  Vector3,
  type Material,
  type Object3D,
} from 'three';
import type { SceneTheme } from '../engine/theme';
import { terrainHeight } from '../data/sample';
import type { Feature, Shape, SiteScene, Structure } from '../data/types';
import { drawSection, faceSpecs, type SectionSpec } from './sections';

/**
 * One site as a block diagram. Everything inside `site` is in site
 * coordinates (x east, y north, z up, metres); `root` maps them to the
 * world: rotated so z is up, centred, scaled to a two-unit block, with the
 * vertical exaggerated on demand.
 */
export class Block {
  readonly root = new Group();
  readonly site = new Group();
  readonly layers = {
    terrain: new Group(),
    walls: new Group(),
    section: new Group(),
    structures: new Group(),
    features: new Group(),
    claimed: new Group(),
    frame: new Group(),
  };
  readonly pickables: Mesh[] = [];
  readonly clip = new Plane(new Vector3(0, 0, -1), 0);
  private terrainMat?: ShaderMaterial;
  private themed: { obj: Material; apply: (t: SceneTheme) => void }[] = [];
  private cut = 0; // 0..1 of the way north
  private scale = 1;
  private exaggeration = 1;
  readonly zTop: number;
  readonly unit: number; // pattern unit in metres
  readonly size: number; // largest horizontal side, metres

  constructor(
    readonly scene: SiteScene,
    private theme: SceneTheme,
  ) {
    const { x, y, z } = scene.extent;
    this.size = Math.max(x[1] - x[0], y[1] - y[0]);
    this.scale = 2 / this.size;
    this.zTop = this.maxGround();
    this.unit = Math.max(this.size / 48, 0.4);
    this.site.rotation.x = -Math.PI / 2;
    this.root.add(this.site);
    for (const g of Object.values(this.layers)) this.site.add(g);
    this.buildTerrain();
    this.buildWalls();
    this.buildStructures();
    this.buildFeatures();
    this.buildFrame();
    this.setExaggeration(this.defaultExaggeration());
    void z;
  }

  /** Deep, narrow sites are shown with less vertical scale so they fit. */
  defaultExaggeration(): number {
    const depth = this.zTop - this.scene.extent.z[0];
    return depth > this.size * 0.9 ? Math.max(0.5, (this.size * 0.9) / depth) : 1;
  }

  private maxGround(): number {
    const t = this.scene.terrain;
    if (t.kind === 'flat') return t.z;
    return Math.max(...t.z);
  }

  private get centre(): [number, number, number] {
    const { x, y, z } = this.scene.extent;
    return [(x[0] + x[1]) / 2, (y[0] + y[1]) / 2, (z[0] + this.zTop) / 2];
  }

  setExaggeration(k: number) {
    this.exaggeration = k;
    const s = this.scale;
    const [cx, cy, cz] = this.centre;
    this.root.scale.set(s, s * k, s);
    this.root.position.set(-cx * s, -cz * s * k, cy * s);
    this.updateClip();
  }

  get verticalExaggeration() {
    return this.exaggeration;
  }

  /** World height of the block (for camera framing). */
  get worldHeight() {
    return (this.zTop - this.scene.extent.z[0]) * this.scale * this.exaggeration;
  }

  /** World position of a site-coordinate point. */
  world(p: [number, number, number]): Vector3 {
    return new Vector3(p[0], p[1], p[2]).applyMatrix4(this.site.matrixWorld);
  }

  // ---------- section ----------

  setCut(f: number) {
    this.cut = Math.min(Math.max(f, 0), 0.98);
    this.updateClip();
    this.buildSection();
  }

  private cutY(): number {
    const { y } = this.scene.extent;
    return y[0] + (y[1] - y[0]) * this.cut;
  }

  private updateClip() {
    // keep y_site >= cut:  world z = root.z - s * y
    this.clip.constant = this.root.position.z - this.scale * this.cutY();
    if (this.cut <= 0) this.clip.constant += 1e3; // nothing clipped
  }

  private buildSection() {
    const g = this.layers.section;
    for (const c of [...g.children]) {
      g.remove(c);
      disposeTree(c);
    }
    if (this.cut <= 0) return;
    const { x, z } = this.scene.extent;
    const y = this.cutY();
    const spec: SectionSpec = { from: [x[0], y], to: [x[1], y], zMin: z[0], zMax: this.zTop, ppm: this.ppm() };
    g.add(this.wallMesh(spec, 0, 1));
    // an inked rim around the cut face
    const pts: Vector3[] = [];
    const n = 64;
    for (let i = 0; i <= n; i++) {
      const xx = x[0] + ((x[1] - x[0]) * i) / n;
      pts.push(new Vector3(xx, y, terrainHeight(this.scene, xx, y)));
    }
    pts.push(new Vector3(x[1], y, z[0]), new Vector3(x[0], y, z[0]), pts[0].clone());
    const segs: Vector3[] = [];
    for (let i = 0; i < pts.length - 1; i++) segs.push(pts[i], pts[i + 1]);
    const rim = new LineSegments(
      new BufferGeometry().setFromPoints(segs),
      new LineBasicMaterial({ color: this.theme.frame, transparent: true, opacity: 0.9 }),
    );
    g.add(rim);
  }

  private ppm(): number {
    const { x, y } = this.scene.extent;
    return Math.min(1600 / Math.max(x[1] - x[0], y[1] - y[0]), 1600 / (this.zTop - this.scene.extent.z[0]), 24);
  }

  // ---------- terrain ----------

  private buildTerrain() {
    const s = this.scene;
    const { x, y } = s.extent;
    let nx: number, ny: number, x0: number, y0: number, dx: number, dy: number;
    let h: (i: number, j: number) => number;
    if (s.terrain.kind === 'grid') {
      const t = s.terrain;
      nx = t.nx;
      ny = t.ny;
      x0 = t.x0;
      y0 = t.y0;
      dx = t.dx;
      dy = t.dx;
      h = (i, j) => t.z[i + t.nx * j];
    } else {
      nx = ny = 2;
      x0 = x[0];
      y0 = y[0];
      dx = x[1] - x[0];
      dy = y[1] - y[0];
      const z0 = s.terrain.z;
      h = () => z0;
    }
    const pos = new Float32Array(nx * ny * 3);
    for (let j = 0; j < ny; j++)
      for (let i = 0; i < nx; i++) {
        const k = (i + nx * j) * 3;
        pos[k] = Math.min(Math.max(x0 + i * dx, x[0]), x[1]);
        pos[k + 1] = Math.min(Math.max(y0 + j * dy, y[0]), y[1]);
        pos[k + 2] = h(i, j);
      }
    const idx: number[] = [];
    for (let j = 0; j < ny - 1; j++)
      for (let i = 0; i < nx - 1; i++) {
        const a = i + nx * j,
          b = a + 1,
          c = a + nx,
          d = c + 1;
        idx.push(a, b, d, a, d, c);
      }
    const geo = new BufferGeometry();
    geo.setAttribute('position', new BufferAttribute(pos, 3));
    geo.setIndex(idx);
    geo.computeVertexNormals();
    const zs = Array.from({ length: nx * ny }, (_, k) => pos[k * 3 + 2]);
    const zMin = Math.min(...zs);
    const zMax = Math.max(...zs);
    const relief = zMax - zMin;
    const interval = relief > 60 ? 5 : relief > 12 ? 1 : relief > 2 ? 0.25 : 0;
    this.terrainMat = new ShaderMaterial({
      uniforms: {
        uLow: { value: new Color() },
        uHigh: { value: new Color() },
        uLightCol: { value: new Color() },
        uContour: { value: new Color() },
        uZ: { value: new Vector2(zMin, Math.max(zMax, zMin + 1)) },
        uInterval: { value: interval },
        uGrid: { value: this.gridStep() },
        uOpacity: { value: 0.9 },
        uLight: { value: new Vector3(-0.55, 0.62, 0.56).normalize() },
      },
      vertexShader: TERRAIN_VS,
      fragmentShader: TERRAIN_FS,
      transparent: true,
      depthWrite: true,
      side: DoubleSide,
      clipping: true,
      clippingPlanes: [this.clip],
    });
    const mesh = new Mesh(geo, this.terrainMat);
    mesh.renderOrder = 1;
    this.layers.terrain.add(mesh);
    this.themed.push({
      obj: this.terrainMat,
      apply: (t) => {
        const u = this.terrainMat!.uniforms;
        u.uLow.value.set(t.terrainLow);
        u.uHigh.value.set(t.terrainHigh);
        u.uLightCol.value.set(t.terrainLight);
        u.uContour.value.set(t.contour);
      },
    });
  }

  private gridStep(): number {
    const s = this.size;
    return s <= 150 ? 10 : s <= 600 ? 50 : 100;
  }

  setGroundOpacity(o: number) {
    if (this.terrainMat) this.terrainMat.uniforms.uOpacity.value = o;
    this.layers.terrain.visible = o > 0.01;
  }

  // ---------- walls ----------

  private buildWalls() {
    for (const c of [...this.layers.walls.children]) {
      this.layers.walls.remove(c);
      disposeTree(c);
    }
    const ppm = this.ppm();
    const light: Record<string, number> = { south: 0.86, east: 1, north: 0.8, west: 0.72 };
    for (const f of faceSpecs(this.scene, this.zTop, ppm)) this.layers.walls.add(this.wallMesh(f.spec, 1, light[f.id]));
    // the floor of the block
    const { x, y, z } = this.scene.extent;
    const floor = new Mesh(
      new BufferGeometry().setFromPoints([
        new Vector3(x[0], y[0], z[0]),
        new Vector3(x[1], y[1], z[0]),
        new Vector3(x[1], y[0], z[0]),
        new Vector3(x[0], y[0], z[0]),
        new Vector3(x[0], y[1], z[0]),
        new Vector3(x[1], y[1], z[0]),
      ]),
      new MeshBasicMaterial({
        color: new Color(this.theme.background).lerp(new Color(this.theme.frame), 0.25),
        side: DoubleSide,
        clippingPlanes: [this.clip],
      }),
    );
    this.layers.walls.add(floor);
  }

  private wallMesh(spec: SectionSpec, clipped: number, light: number): Mesh {
    const cv = drawSection(this.scene, spec, this.theme, this.unit);
    const tex = new CanvasTexture(cv);
    tex.colorSpace = SRGBColorSpace;
    tex.anisotropy = 4;
    const n = 96;
    const pos: number[] = [];
    const uv: number[] = [];
    const col: number[] = [];
    const idx: number[] = [];
    const span = spec.zMax - spec.zMin;
    const night = this.theme.name === 'night';
    const deep = night ? 0.42 : 0.78; // shade at the floor relative to the top
    for (let i = 0; i <= n; i++) {
      const t = i / n;
      const xx = spec.from[0] + (spec.to[0] - spec.from[0]) * t;
      const yy = spec.from[1] + (spec.to[1] - spec.from[1]) * t;
      const top = Math.min(terrainHeight(this.scene, xx, yy), spec.zMax);
      pos.push(xx, yy, top, xx, yy, spec.zMin);
      uv.push(t, (top - spec.zMin) / span, t, 0);
      col.push(light, light, light, light * deep, light * deep, light * deep);
      if (i < n) {
        const a = i * 2;
        idx.push(a, a + 1, a + 3, a, a + 3, a + 2);
      }
    }
    const geo = new BufferGeometry();
    geo.setAttribute('position', new BufferAttribute(new Float32Array(pos), 3));
    geo.setAttribute('uv', new BufferAttribute(new Float32Array(uv), 2));
    geo.setAttribute('color', new BufferAttribute(new Float32Array(col), 3));
    geo.setIndex(idx);
    const mat = new MeshBasicMaterial({
      map: tex,
      vertexColors: true,
      side: FrontSide,
      transparent: false,
      alphaTest: 0.5,
      clippingPlanes: clipped ? [this.clip] : [],
    });
    const m = new Mesh(geo, mat);
    m.renderOrder = 0;
    return m;
  }

  // ---------- structures and features ----------

  private buildStructures() {
    for (const st of this.scene.structures) this.layers.structures.add(this.structureMesh(st));
  }

  private structureMesh(st: Structure): Group {
    const g = new Group();
    g.name = st.id;
    const geo = shapeGeometry(st.shape);
    const fill = new MeshLambertMaterial({ color: this.theme.structure, transparent: true, opacity: 0.12, clippingPlanes: [this.clip] });
    const edges = new LineBasicMaterial({ color: this.theme.structureEdge, transparent: true, opacity: 0.75, clippingPlanes: [this.clip] });
    g.add(new Mesh(geo, fill), new LineSegments(new EdgesGeometry(geo, 20), edges));
    this.themed.push({
      obj: fill,
      apply: (t) => {
        fill.color.set(t.structure);
        fill.opacity = t.name === 'night' ? 0.1 : 0.78;
        edges.color.set(t.structureEdge);
        edges.opacity = t.name === 'night' ? 0.7 : 0.55;
      },
    });
    return g;
  }

  private buildFeatures() {
    for (const f of this.scene.features) {
      const g = this.featureMesh(f);
      const plumb = this.plumbLine(f);
      if (plumb) g.add(plumb);
      (f.status === 'claimed' ? this.layers.claimed : this.layers.features).add(g);
    }
  }

  /** A dashed line from a buried feature up to the ground, and a small mark there. */
  private plumbLine(f: Feature): LineSegments | undefined {
    const box = shapeGeometry(f.shape);
    box.computeBoundingBox();
    const bb = box.boundingBox!;
    box.dispose();
    const cx = (bb.min.x + bb.max.x) / 2;
    const cy = (bb.min.y + bb.max.y) / 2;
    const ground = terrainHeight(this.scene, cx, cy);
    if (bb.max.z >= ground - this.unit * 0.8) return undefined;
    const r = Math.max((bb.max.x - bb.min.x) / 2, this.unit * 0.8);
    const pts = [new Vector3(cx, cy, bb.max.z), new Vector3(cx, cy, ground)];
    const n = 24;
    for (let i = 0; i < n; i++) {
      const a0 = (i / n) * Math.PI * 2;
      const a1 = ((i + 1) / n) * Math.PI * 2;
      pts.push(new Vector3(cx + r * Math.cos(a0), cy + r * Math.sin(a0), ground), new Vector3(cx + r * Math.cos(a1), cy + r * Math.sin(a1), ground));
    }
    const claimed = f.status === 'claimed';
    const mat = new LineDashedMaterial({
      color: claimed ? this.theme.claimed : this.theme.void,
      dashSize: this.unit * 0.7,
      gapSize: this.unit * 0.55,
      transparent: true,
      opacity: 0.55,
      depthTest: false,
      clippingPlanes: [this.clip],
    });
    const line = new LineSegments(new BufferGeometry().setFromPoints(pts), mat);
    line.computeLineDistances();
    line.renderOrder = 4;
    this.themed.push({ obj: mat, apply: (t) => mat.color.set(claimed ? t.claimed : t.void) });
    return line;
  }

  private featureMesh(f: Feature): Group {
    const g = new Group();
    g.name = f.id;
    const geo = shapeGeometry(f.shape);
    const claimed = f.status === 'claimed';
    const fill = new MeshBasicMaterial({
      color: claimed ? this.theme.claimed : this.theme.void,
      transparent: true,
      opacity: claimed ? 0.14 : 0.42,
      depthTest: false,
      depthWrite: false,
      side: DoubleSide,
      clippingPlanes: [this.clip],
    });
    const edgeGeo = new EdgesGeometry(geo, 25);
    const line = claimed
      ? new LineSegments(
          edgeGeo,
          new LineDashedMaterial({
            color: this.theme.claimed,
            dashSize: this.unit * 1.6,
            gapSize: this.unit * 1.1,
            transparent: true,
            opacity: 0.95,
            depthTest: false,
            clippingPlanes: [this.clip],
          }),
        )
      : new LineSegments(
          edgeGeo,
          new LineBasicMaterial({ color: this.theme.void, transparent: true, opacity: 0.95, depthTest: false, clippingPlanes: [this.clip] }),
        );
    if (claimed) line.computeLineDistances();
    const mesh = new Mesh(geo, fill);
    mesh.userData.feature = f;
    mesh.renderOrder = 5;
    line.renderOrder = 6;
    g.add(mesh, line);
    this.pickables.push(mesh);
    const lm = line.material as LineBasicMaterial;
    this.themed.push({
      obj: fill,
      apply: (t) => {
        const c = claimed ? t.claimed : t.void;
        fill.color.set(c);
        lm.color.set(c);
        const night = t.name === 'night';
        fill.blending = night ? AdditiveBlending : NormalBlending;
        fill.opacity = claimed ? (night ? 0.16 : 0.12) : night ? 0.5 : 0.55;
        fill.needsUpdate = true;
      },
    });
    return g;
  }

  // ---------- frame ----------

  private buildFrame() {
    const { x, y, z } = this.scene.extent;
    const b = z[0];
    const c = (xx: number, yy: number) => new Vector3(xx, yy, b);
    const v = (xx: number, yy: number) => [new Vector3(xx, yy, b), new Vector3(xx, yy, terrainHeight(this.scene, xx, yy))];
    const pts = [
      c(x[0], y[0]), c(x[1], y[0]), c(x[1], y[0]), c(x[1], y[1]), c(x[1], y[1]), c(x[0], y[1]), c(x[0], y[1]), c(x[0], y[0]),
      ...v(x[0], y[0]), ...v(x[1], y[0]), ...v(x[1], y[1]), ...v(x[0], y[1]),
    ];
    const mat = new LineBasicMaterial({ color: this.theme.frame, transparent: true, opacity: 0.8, clippingPlanes: [this.clip] });
    this.layers.frame.add(new LineSegments(new BufferGeometry().setFromPoints(pts), mat));
    // north arrow on the ground at the north-west corner
    const s = this.size * 0.045;
    const ax = x[0] + s * 1.6;
    const ay = y[1] - s * 2.6;
    const az = terrainHeight(this.scene, ax, ay) + this.size * 0.004;
    const arrow = new ThreeShape();
    arrow.moveTo(0, s);
    arrow.lineTo(s * 0.42, -s * 0.6);
    arrow.lineTo(0, -s * 0.28);
    arrow.lineTo(-s * 0.42, -s * 0.6);
    arrow.closePath();
    const am = new MeshBasicMaterial({ color: this.theme.frame, side: DoubleSide, transparent: true, opacity: 0.9 });
    const a = new Mesh(new ExtrudeGeometry(arrow, { depth: 0.001, bevelEnabled: false }), am);
    a.position.set(ax, ay, az);
    this.layers.frame.add(a);
    this.themed.push({
      obj: mat,
      apply: (t) => {
        mat.color.set(t.frame);
        am.color.set(t.frame);
      },
    });
  }

  // ---------- theme ----------

  setTheme(t: SceneTheme) {
    this.theme = t;
    for (const th of this.themed) th.apply(t);
    this.buildWalls();
    this.buildSection();
  }

  applyTheme() {
    for (const th of this.themed) th.apply(this.theme);
  }

  dispose() {
    disposeTree(this.root);
  }
}

// ---------- geometry ----------

export function shapeGeometry(shape: Shape): BufferGeometry {
  let g: BufferGeometry;
  switch (shape.type) {
    case 'box': {
      g = new BoxGeometry(shape.size[0], shape.size[1], shape.size[2]);
      const m = new Matrix4()
        .makeRotationZ(((shape.yaw_deg ?? 0) * Math.PI) / 180)
        .multiply(new Matrix4().makeRotationX(((shape.pitch_deg ?? 0) * Math.PI) / 180));
      g.applyMatrix4(m);
      g.translate(...shape.centre);
      break;
    }
    case 'cylinder':
      g = new CylinderGeometry(shape.radius, shape.radius, shape.height, 32, 1);
      g.rotateX(Math.PI / 2);
      g.translate(...shape.centre);
      break;
    case 'sphere':
      g = new SphereGeometry(shape.radius, 24, 16);
      g.translate(...shape.centre);
      break;
    case 'prism': {
      const s = new ThreeShape(shape.polygon.map(([x, y]) => new Vector2(x, y)));
      g = new ExtrudeGeometry(s, { depth: shape.top - shape.bottom, bevelEnabled: false });
      g.translate(0, 0, shape.bottom);
      break;
    }
    case 'pyramid': {
      g = new ConeGeometry(shape.base / Math.SQRT2, shape.height, 4, 1);
      g.rotateY(Math.PI / 4);
      g.translate(0, shape.height / 2, 0);
      g.rotateX(Math.PI / 2);
      g.rotateZ(((shape.yaw_deg ?? 0) * Math.PI) / 180);
      g.translate(...shape.centre);
      break;
    }
  }
  return g;
}

function disposeTree(o: Object3D) {
  o.traverse((c) => {
    const m = c as Mesh;
    m.geometry?.dispose?.();
    const mats = Array.isArray(m.material) ? m.material : m.material ? [m.material] : [];
    for (const mt of mats) {
      (mt as MeshBasicMaterial).map?.dispose();
      mt.dispose();
    }
  });
}

// ---------- terrain shader: raking light, elevation tint, contours ----------

const TERRAIN_VS = /* glsl */ `
#include <clipping_planes_pars_vertex>
varying vec3 vNormalW;
varying vec3 vSite;
void main() {
  vSite = position;
  vNormalW = normalize(transpose(inverse(mat3(modelMatrix))) * normal);
  vec4 mvPosition = modelViewMatrix * vec4(position, 1.0);
  #include <clipping_planes_vertex>
  gl_Position = projectionMatrix * mvPosition;
}`;

const TERRAIN_FS = /* glsl */ `
#include <clipping_planes_pars_fragment>
uniform vec3 uLow, uHigh, uLightCol, uContour, uLight;
uniform vec2 uZ;
uniform float uInterval, uGrid, uOpacity;
varying vec3 vNormalW;
varying vec3 vSite;
float lineAA(float v, float w) {
  float f = abs(fract(v - 0.5) - 0.5) / max(fwidth(v), 1e-5);
  return 1.0 - smoothstep(0.0, w, f);
}
void main() {
  #include <clipping_planes_fragment>
  vec3 n = normalize(vNormalW);
  if (!gl_FrontFacing) n = -n;
  float shade = clamp(dot(n, uLight), 0.0, 1.0);
  float t = clamp((vSite.z - uZ.x) / (uZ.y - uZ.x), 0.0, 1.0);
  vec3 base = mix(uLow, uHigh, smoothstep(0.0, 1.0, t));
  vec3 col = base * (0.55 + 0.6 * shade) + uLightCol * pow(shade, 6.0) * 0.12;
  if (uInterval > 0.0) {
    float c = lineAA(vSite.z / uInterval, 1.0);
    float major = lineAA(vSite.z / (uInterval * 5.0), 1.2);
    col = mix(col, uContour, max(c * 0.35, major * 0.6));
  }
  float g = max(lineAA(vSite.x / uGrid, 0.8), lineAA(vSite.y / uGrid, 0.8));
  col = mix(col, uContour, g * 0.18);
  gl_FragColor = vec4(col, uOpacity);
  #include <colorspace_fragment>
}`;
