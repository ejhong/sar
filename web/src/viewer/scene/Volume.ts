import {
  AdditiveBlending,
  BackSide,
  BoxGeometry,
  Color,
  Data3DTexture,
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

/**
 * A voxel volume (a tomogram) ray-marched inside its box, in site
 * coordinates. Values are 0..255: 0 is background, 255 the strongest
 * anomaly the exporter mapped. What an instrument recovered glows in
 * verdigris; the section plane cuts it like everything else.
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
    const [x0, y0, ztop] = info.origin; // centre of the first cell; z of the top layer
    const lo = new Vector3(x0 - h / 2, y0 - h / 2, ztop - (nz - 0.5) * h);
    const hi = new Vector3(x0 + (nx - 0.5) * h, y0 + (ny - 0.5) * h, ztop + h / 2);
    const size = hi.clone().sub(lo);
    const geo = new BoxGeometry(size.x, size.y, size.z);
    geo.translate((lo.x + hi.x) / 2, (lo.y + hi.y) / 2, (lo.z + hi.z) / 2);
    this.mat = new ShaderMaterial({
      uniforms: {
        uTex: { value: tex },
        uLo: { value: lo },
        uHi: { value: hi },
        uColor: { value: new Color(theme.recovered) },
        uThreshold: { value: 0.18 },
        uDensity: { value: 1.6 },
        uCutY: { value: -1e9 },
        uSteps: { value: 160 },
        uCamLocal: { value: new Vector3() },
        uDay: { value: theme.name === 'day' ? 1 : 0 },
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
    this.mat.uniforms.uColor.value.set(t.recovered);
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
uniform float uThreshold, uDensity, uCutY, uDay;
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
    vec3 c = mix(uColor * 0.55, uColor * 1.35, smoothstep(uThreshold, 1.0, v));
    acc.rgb += (1.0 - acc.a) * a * c;
    acc.a += (1.0 - acc.a) * a;
  }
  if (acc.a < 0.003) discard;
  gl_FragColor = uDay > 0.5 ? vec4(acc.rgb / max(acc.a, 1e-4), acc.a) : vec4(acc.rgb, acc.a);
  #include <colorspace_fragment>
}`;
