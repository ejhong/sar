import {
  AdditiveBlending,
  Color,
  Data3DTexture,
  DoubleSide,
  LinearFilter,
  Mesh,
  NormalBlending,
  PlaneGeometry,
  RedFormat,
  ShaderMaterial,
  UnsignedByteType,
} from 'three';
import type { SceneTheme } from '../engine/theme';

export interface WavefieldInfo {
  id: string;
  label: string;
  caption: string;
  file: string;
  shape: [number, number, number]; // frames, nx, ny
  origin: [number, number]; // x, y of the first cell centre
  spacing: number;
  z: number; // height of the ground it is drawn on
  dt_ms: number;
  kind: 'motion' | 'echo';
  shot: [number, number, number];
}

/**
 * Vertical ground velocity played back on the surface as light: the waves a
 * hammer blow sends out, or the chamber's echo alone. Values are signed,
 * stored as 0..255 about 128 with a square-root gain so the faint echo is
 * visible beside the direct wave.
 */
export class Wavefield {
  readonly mesh: Mesh;
  private mat: ShaderMaterial;
  private t = 0;
  playing = true;
  speed = 1;

  constructor(
    readonly info: WavefieldInfo,
    data: Uint8Array,
    theme: SceneTheme,
    extent: { x: [number, number]; y: [number, number] },
  ) {
    const [nf, nx, ny] = info.shape;
    const tex = new Data3DTexture(data, nx, ny, nf);
    tex.format = RedFormat;
    tex.type = UnsignedByteType;
    tex.minFilter = tex.magFilter = LinearFilter;
    tex.unpackAlignment = 1;
    tex.needsUpdate = true;
    const w = extent.x[1] - extent.x[0];
    const hgt = extent.y[1] - extent.y[0];
    const geo = new PlaneGeometry(w, hgt, 1, 1);
    geo.translate((extent.x[0] + extent.x[1]) / 2, (extent.y[0] + extent.y[1]) / 2, info.z + 0.02);
    const h = info.spacing;
    this.mat = new ShaderMaterial({
      uniforms: {
        uTex: { value: tex },
        uFrame: { value: 0 },
        uFrames: { value: nf },
        uLo: { value: [info.origin[0] - h / 2, info.origin[1] - h / 2] },
        uSize: { value: [nx * h, ny * h] },
        uColor: { value: new Color(info.kind === 'echo' ? theme.void : 0xf3ead6) },
        uDay: { value: theme.name === 'day' ? 1 : 0 },
        uGain: { value: info.kind === 'echo' ? 1.4 : 1.0 },
      },
      vertexShader: VS,
      fragmentShader: FS,
      transparent: true,
      depthWrite: false,
      side: DoubleSide,
      blending: theme.name === 'day' ? NormalBlending : AdditiveBlending,
    });
    this.mesh = new Mesh(geo, this.mat);
    this.mesh.renderOrder = 2;
  }

  get frames() {
    return this.info.shape[0];
  }

  /** Advance by real seconds; the playback shows about 1.2 s per record. */
  tick(dt: number) {
    if (!this.playing) return;
    this.t = (this.t + (dt * this.speed * this.frames) / 2.4) % (this.frames + 12);
    this.mat.uniforms.uFrame.value = Math.min(this.t, this.frames - 1);
  }

  setFrame(f: number) {
    this.t = f;
    this.mat.uniforms.uFrame.value = f;
  }

  get time_ms() {
    return Math.min(this.t, this.frames - 1) * this.info.dt_ms;
  }

  setTheme(t: SceneTheme) {
    this.mat.uniforms.uColor.value.set(this.info.kind === 'echo' ? t.void : t.name === 'day' ? 0x3a3632 : 0xf3ead6);
    this.mat.uniforms.uDay.value = t.name === 'day' ? 1 : 0;
    this.mat.blending = t.name === 'day' ? NormalBlending : AdditiveBlending;
    this.mat.needsUpdate = true;
  }

  dispose() {
    this.mesh.geometry.dispose();
    (this.mat.uniforms.uTex.value as Data3DTexture).dispose();
    this.mat.dispose();
  }
}

const VS = /* glsl */ `
varying vec2 vXY;
void main() {
  vXY = position.xy;
  gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
}`;

const FS = /* glsl */ `
precision highp float;
precision highp sampler3D;
uniform sampler3D uTex;
uniform float uFrame, uFrames, uDay, uGain;
uniform vec2 uLo, uSize;
uniform vec3 uColor;
varying vec2 vXY;
void main() {
  vec2 uv = (vXY - uLo) / uSize;
  if (any(lessThan(uv, vec2(0.0))) || any(greaterThan(uv, vec2(1.0)))) discard;
  float w = (uFrame + 0.5) / uFrames;
  float v = (texture(uTex, vec3(uv, w)).r - 0.5) * 2.0 * uGain;
  float a = clamp(abs(v), 0.0, 1.0);
  // crests bright, troughs dimmer: a sheet of light that breathes with the wave
  float crest = v > 0.0 ? 1.0 : 0.55;
  vec3 c = uColor * crest;
  if (a < 0.02) discard;
  gl_FragColor = uDay > 0.5 ? vec4(c, a * 0.85) : vec4(c * a * 1.25, a);
  #include <colorspace_fragment>
}`;
