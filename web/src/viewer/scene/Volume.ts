import {
  AdditiveBlending,
  BackSide,
  BoxGeometry,
  Color,
  Data3DTexture,
  DataTexture,
  RGBAFormat,
  LinearFilter,
  Mesh,
  NormalBlending,
  RedFormat,
  ShaderMaterial,
  UnsignedByteType,
  Vector3,
  type Matrix4,
} from 'three';
import type { SceneTheme } from '../engine/theme';
import type { VolumeInfo } from '../data/types';
import { MAGMA } from './magma';

/** The magma table as a texture, shared by every volume that asks for it. */
let magmaLut: DataTexture | null = null;
function magma(): DataTexture {
  if (magmaLut) return magmaLut;
  const d = new Uint8Array(MAGMA.length * 4);
  MAGMA.forEach(([r, g, b], i) => d.set([Math.round(r * 255), Math.round(g * 255), Math.round(b * 255), 255], i * 4));
  magmaLut = new DataTexture(d, MAGMA.length, 1, RGBAFormat);
  magmaLut.minFilter = magmaLut.magFilter = LinearFilter;
  magmaLut.needsUpdate = true;
  return magmaLut;
}

/**
 * A voxel volume (a tomogram) ray-marched inside its box, in site
 * coordinates. Values are 0..255: 0 is background, 255 the strongest
 * anomaly the exporter mapped. What an instrument recovered glows in
 * verdigris, what the satellite's processing computes in cinnabar, each with its
 * lightness following the value; the section plane cuts it like everything else.
 */
export class Volume {
  readonly mesh: Mesh;
  private mat: ShaderMaterial;

  constructor(
    readonly info: VolumeInfo,
    data: Uint8Array,
    theme: SceneTheme,
  ) {
    const [nx, ny, nz] = info.shape;
    const tex = new Data3DTexture(data, nx, ny, nz);
    tex.format = RedFormat;
    tex.type = UnsignedByteType;
    tex.minFilter = tex.magFilter = LinearFilter;
    tex.unpackAlignment = 1;
    tex.needsUpdate = true;
    const h = info.spacing;
    const hz = info.spacing_z ?? h; // a finer depth step where the volume has one
    const [x0, y0, ztop] = info.origin; // centre of the first cell; z of the top layer
    const lo = new Vector3(x0 - h / 2, y0 - h / 2, ztop - (nz - 0.5) * hz);
    const hi = new Vector3(x0 + (nx - 0.5) * h, y0 + (ny - 0.5) * h, ztop + hz / 2);
    const size = hi.clone().sub(lo);
    const geo = new BoxGeometry(size.x, size.y, size.z);
    geo.translate((lo.x + hi.x) / 2, (lo.y + hi.y) / 2, (lo.z + hi.z) / 2);
    this.mat = new ShaderMaterial({
      uniforms: {
        uTex: { value: tex },
        uLo: { value: lo },
        uHi: { value: hi },
        uColor: { value: new Color(tint(info, theme)) },
        uThreshold: { value: 0.18 },
        uDensity: { value: info.tint === 'gated' ? 6.0 : 1.6 }, // sparse fit-score columns need more opacity to read
        uCutY: { value: -1e9 },
        // enough steps along a ray to see every layer of a tall, finely layered volume (the shader stops at 512)
        uSteps: { value: Math.min(512, Math.max(160, Math.ceil(1.5 * Math.max(nx, ny, nz)))) },
        uCamLocal: { value: new Vector3() },
        uDay: { value: theme.name === 'day' ? 1 : 0 },
        uCmap: { value: info.cmap === 'magma' ? 1 : 0 },
        uLut: { value: magma() },
      },
      vertexShader: VS,
      fragmentShader: FS,
      side: BackSide,
      transparent: true,
      depthWrite: false,
      depthTest: false,
      blending: theme.name === 'day' ? NormalBlending : AdditiveBlending,
    });
    this.mesh = new Mesh(geo, this.mat);
    this.mesh.renderOrder = 3;
    this.mesh.onBeforeRender = (_r, _s, camera) => {
      const inv = (this.mesh.matrixWorld as Matrix4).clone().invert();
      this.mat.uniforms.uCamLocal.value.copy(camera.position).applyMatrix4(inv);
    };
  }

  setTheme(t: SceneTheme) {
    this.mat.uniforms.uColor.value.set(tint(this.info, t));
    this.mat.uniforms.uDay.value = t.name === 'day' ? 1 : 0;
    this.mat.blending = t.name === 'day' ? NormalBlending : AdditiveBlending;
    this.mat.needsUpdate = true;
  }

  setThreshold(v: number) {
    this.mat.uniforms.uThreshold.value = v;
  }

  setCut(y: number) {
    this.mat.uniforms.uCutY.value = y;
  }

  dispose() {
    this.mesh.geometry.dispose();
    (this.mat.uniforms.uTex.value as Data3DTexture).dispose();
    this.mat.dispose();
  }
}

/** A volume's colour is its instrument's: the satellite's in cinnabar whichever method drew it, the geophones' in
 * verdigris. What the voxels hold (a material property, focused power or a fit score) the lab says in words. */
function tint(info: VolumeInfo, t: SceneTheme): number {
  return info.status === 'radar' ? t.radar : t.recovered;
}

const VS = /* glsl */ `
varying vec3 vLocal;
void main() {
  vLocal = position;
  gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
}`;

const FS = /* glsl */ `
precision highp float;
precision highp sampler3D;
uniform sampler3D uTex;
uniform vec3 uLo, uHi, uColor, uCamLocal;
uniform float uThreshold, uDensity, uCutY, uDay, uCmap;
uniform sampler2D uLut;
uniform int uSteps;
varying vec3 vLocal;

vec2 hitBox(vec3 o, vec3 d) {
  vec3 inv = 1.0 / d;
  vec3 t0 = (uLo - o) * inv, t1 = (uHi - o) * inv;
  vec3 tmin = min(t0, t1), tmax = max(t0, t1);
  return vec2(max(max(tmin.x, tmin.y), tmin.z), min(min(tmax.x, tmax.y), tmax.z));
}

void main() {
  vec3 o = uCamLocal;
  vec3 d = normalize(vLocal - o);
  vec2 t = hitBox(o, d);
  t.x = max(t.x, 0.0);
  if (t.x >= t.y) discard;
  float len = t.y - t.x;
  float dt = len / float(uSteps);
  vec3 size = uHi - uLo;
  vec4 acc = vec4(0.0);
  float jitter = fract(sin(dot(gl_FragCoord.xy, vec2(12.9898, 78.233))) * 43758.5453);
  for (int i = 0; i < 512; i++) {
    if (i >= uSteps || acc.a > 0.97) break;
    vec3 p = o + d * (t.x + (float(i) + jitter) * dt);
    if (p.y < uCutY) continue;
    vec3 uvw = vec3((p.x - uLo.x) / size.x, (p.y - uLo.y) / size.y, (uHi.z - p.z) / size.z);
    float v = texture(uTex, uvw).r;
    float a = smoothstep(uThreshold, 1.0, v) * uDensity * dt / max(size.x, size.y) * 18.0;
    a = clamp(a, 0.0, 1.0);
    // one hue, its lightness the value: dim and deep for the lowest shown, near white (by day, near ink) for the highest
    float s = smoothstep(uThreshold, 1.0, v);
    vec3 c = uDay > 0.5 ? mix(mix(uColor, vec3(1.0), 0.55), uColor * 0.55, s) : mix(uColor * 0.45, mix(uColor, vec3(1.0), 0.6), s);
    // a volume with its own colour scale (the gated reconstruction's magma) takes its colour from the value itself
    if (uCmap > 0.5) c = texture2D(uLut, vec2((0.5 + clamp(v, 0.0, 1.0) * 32.0) / 33.0, 0.5)).rgb;
    acc.rgb += (1.0 - acc.a) * a * c;
    acc.a += (1.0 - acc.a) * a;
  }
  if (acc.a < 0.003) discard;
  gl_FragColor = uDay > 0.5 ? vec4(acc.rgb / max(acc.a, 1e-4), acc.a) : vec4(acc.rgb, acc.a);
  #include <colorspace_fragment>
}`;
