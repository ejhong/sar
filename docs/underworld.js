/* Underworld viewer: WebGL2 ray-marched tomogram volumes with the radar surface on top.
 *
 * No framework. One 3-D texture per patch, one 2-D texture for the surface, and an analytic
 * ray-box march in the fragment shader. Bundles are produced by experiments/export_voxels.py
 * and described by data/voxels/index.json, so a new site is new data, not new code.
 */
'use strict';

const VERT = `#version 300 es
in vec2 position;
out vec2 uv;
void main() { uv = position; gl_Position = vec4(position, 0.0, 1.0); }`;

const FRAG = `#version 300 es
precision highp float;
precision highp sampler3D;
in vec2 uv;
out vec4 frag;

uniform sampler3D volume;
uniform sampler2D surface;
uniform mat4 invViewProj;
uniform vec3 cameraPos;
uniform vec3 boxMin;
uniform vec3 boxMax;
uniform float threshold;
uniform float density;
uniform float depthClip;      // 0..1 of the box depth that is drawn
uniform float surfaceMix;
uniform float repeatFrac;     // where the volume repeats, as a fraction of box depth
uniform int   showRepeat;
uniform vec2  voidBand;       // known-void depth range as fractions of box depth; (0,0) = none
uniform int   showVoids;
uniform int   steps;

vec3 ramp(float t) {
  // dark blue -> teal -> amber -> white; monotonic in lightness so structure reads honestly
  const vec3 c0 = vec3(0.043, 0.075, 0.180);
  const vec3 c1 = vec3(0.090, 0.330, 0.450);
  const vec3 c2 = vec3(0.310, 0.620, 0.520);
  const vec3 c3 = vec3(0.880, 0.650, 0.270);
  const vec3 c4 = vec3(1.000, 0.960, 0.900);
  if (t < 0.25) return mix(c0, c1, t / 0.25);
  if (t < 0.50) return mix(c1, c2, (t - 0.25) / 0.25);
  if (t < 0.75) return mix(c2, c3, (t - 0.50) / 0.25);
  return mix(c3, c4, (t - 0.75) / 0.25);
}

bool hitBox(vec3 o, vec3 d, out float t0, out float t1) {
  vec3 inv = 1.0 / d;
  vec3 a = (boxMin - o) * inv;
  vec3 b = (boxMax - o) * inv;
  vec3 lo = min(a, b), hi = max(a, b);
  t0 = max(max(lo.x, lo.y), lo.z);
  t1 = min(min(hi.x, hi.y), hi.z);
  return t1 > max(t0, 0.0);
}

void main() {
  vec4 near = invViewProj * vec4(uv, -1.0, 1.0);
  vec4 far  = invViewProj * vec4(uv,  1.0, 1.0);
  vec3 origin = near.xyz / near.w;
  vec3 dir = normalize(far.xyz / far.w - origin);

  float t0, t1;
  if (!hitBox(origin, dir, t0, t1)) { frag = vec4(0.0); return; }
  t0 = max(t0, 0.0);

  vec3 span = boxMax - boxMin;
  // the surface sits at the top of the box (y = boxMax.y)
  float tSurface = (boxMax.y - origin.y) / dir.y;

  vec4 acc = vec4(0.0);
  float dt = (t1 - t0) / float(steps);

  // The surface is a single thin layer, so composite it once at its exact crossing rather
  // than hoping a march step lands on it; stepping past it leaves the ground stippled.
  float tStart = t0;
  if (surfaceMix > 0.0 && tSurface >= t0 && tSurface <= t1) {
    vec3 sg = (origin + dir * tSurface - boxMin) / span;
    if (all(greaterThanEqual(sg.xz, vec2(0.0))) && all(lessThanEqual(sg.xz, vec2(1.0)))) {
      float s = texture(surface, vec2(sg.x, 1.0 - sg.z)).r;
      vec3 sc = mix(vec3(0.06, 0.07, 0.08), vec3(0.99, 0.98, 0.95), s);
      acc.rgb += sc * surfaceMix;
      acc.a   += surfaceMix;
      if (dir.y < 0.0) tStart = max(t0, tSurface);   // looking down: start below the ground
    }
  }

  for (int i = 0; i < 512; i++) {
    if (i >= steps || acc.a > 0.985) break;
    float t = tStart + (float(i) + 0.5) * (t1 - tStart) / float(steps);
    vec3 p = origin + dir * t;
    vec3 g = (p - boxMin) / span;              // 0..1 in box space

    float depthFrac = 1.0 - g.y;               // 0 at the surface, 1 at the bottom
    if (depthFrac > depthClip) continue;

    float v = texture(volume, vec3(g.x, g.z, depthFrac)).r;
    float a = smoothstep(threshold, 1.0, v) * density;
    if (a <= 0.0) continue;
    vec3 c = ramp(v);
    if (showVoids == 1 && voidBand.y > voidBand.x &&
        depthFrac >= voidBand.x && depthFrac <= voidBand.y) {
      c = mix(c, vec3(0.30, 0.86, 0.62), 0.55);      // where the surveyed shafts actually are
      a = max(a, 0.05);
    }
    if (showRepeat == 1 && abs(depthFrac - repeatFrac) < 0.004) {
      c = vec3(1.0, 0.35, 0.25); a = max(a, 0.5);
    }
    a = 1.0 - pow(1.0 - a, ((t1 - tStart) / float(steps)) * 220.0);
    acc.rgb += (1.0 - acc.a) * c * a;
    acc.a   += (1.0 - acc.a) * a;
  }
  frag = acc;
}`;

function compile(gl, type, src) {
  const s = gl.createShader(type);
  gl.shaderSource(s, src);
  gl.compileShader(s);
  if (!gl.getShaderParameter(s, gl.COMPILE_STATUS)) throw new Error(gl.getShaderInfoLog(s));
  return s;
}

function mat4Multiply(a, b) {
  const out = new Float32Array(16);
  for (let i = 0; i < 4; i++)
    for (let j = 0; j < 4; j++) {
      let v = 0;
      for (let k = 0; k < 4; k++) v += a[k * 4 + j] * b[i * 4 + k];
      out[i * 4 + j] = v;
    }
  return out;
}

function perspective(fovy, aspect, near, far) {
  const f = 1 / Math.tan(fovy / 2);
  const out = new Float32Array(16);
  out[0] = f / aspect; out[5] = f; out[11] = -1;
  out[10] = (far + near) / (near - far);
  out[14] = (2 * far * near) / (near - far);
  return out;
}

function lookAt(eye, centre, up) {
  const z = normalize(sub(eye, centre));
  const x = normalize(cross(up, z));
  const y = cross(z, x);
  return new Float32Array([
    x[0], y[0], z[0], 0, x[1], y[1], z[1], 0, x[2], y[2], z[2], 0,
    -dot(x, eye), -dot(y, eye), -dot(z, eye), 1]);
}

const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
const normalize = (a) => { const l = Math.hypot(...a) || 1; return [a[0] / l, a[1] / l, a[2] / l]; };

function invert(m) {
  const inv = new Float32Array(16), a = m;
  inv[0] = a[5]*a[10]*a[15]-a[5]*a[11]*a[14]-a[9]*a[6]*a[15]+a[9]*a[7]*a[14]+a[13]*a[6]*a[11]-a[13]*a[7]*a[10];
  inv[4] = -a[4]*a[10]*a[15]+a[4]*a[11]*a[14]+a[8]*a[6]*a[15]-a[8]*a[7]*a[14]-a[12]*a[6]*a[11]+a[12]*a[7]*a[10];
  inv[8] = a[4]*a[9]*a[15]-a[4]*a[11]*a[13]-a[8]*a[5]*a[15]+a[8]*a[7]*a[13]+a[12]*a[5]*a[11]-a[12]*a[7]*a[9];
  inv[12] = -a[4]*a[9]*a[14]+a[4]*a[10]*a[13]+a[8]*a[5]*a[14]-a[8]*a[6]*a[13]-a[12]*a[5]*a[10]+a[12]*a[6]*a[9];
  inv[1] = -a[1]*a[10]*a[15]+a[1]*a[11]*a[14]+a[9]*a[2]*a[15]-a[9]*a[3]*a[14]-a[13]*a[2]*a[11]+a[13]*a[3]*a[10];
  inv[5] = a[0]*a[10]*a[15]-a[0]*a[11]*a[14]-a[8]*a[2]*a[15]+a[8]*a[3]*a[14]+a[12]*a[2]*a[11]-a[12]*a[3]*a[10];
  inv[9] = -a[0]*a[9]*a[15]+a[0]*a[11]*a[13]+a[8]*a[1]*a[15]-a[8]*a[3]*a[13]-a[12]*a[1]*a[11]+a[12]*a[3]*a[9];
  inv[13] = a[0]*a[9]*a[14]-a[0]*a[10]*a[13]-a[8]*a[1]*a[14]+a[8]*a[2]*a[13]+a[12]*a[1]*a[10]-a[12]*a[2]*a[9];
  inv[2] = a[1]*a[6]*a[15]-a[1]*a[7]*a[14]-a[5]*a[2]*a[15]+a[5]*a[3]*a[14]+a[13]*a[2]*a[7]-a[13]*a[3]*a[6];
  inv[6] = -a[0]*a[6]*a[15]+a[0]*a[7]*a[14]+a[4]*a[2]*a[15]-a[4]*a[3]*a[14]-a[12]*a[2]*a[7]+a[12]*a[3]*a[6];
  inv[10] = a[0]*a[5]*a[15]-a[0]*a[7]*a[13]-a[4]*a[1]*a[15]+a[4]*a[3]*a[13]+a[12]*a[1]*a[7]-a[12]*a[3]*a[5];
  inv[14] = -a[0]*a[5]*a[14]+a[0]*a[6]*a[13]+a[4]*a[1]*a[14]-a[4]*a[2]*a[13]-a[12]*a[1]*a[6]+a[12]*a[2]*a[5];
  inv[3] = -a[1]*a[6]*a[11]+a[1]*a[7]*a[10]+a[5]*a[2]*a[11]-a[5]*a[3]*a[10]-a[9]*a[2]*a[7]+a[9]*a[3]*a[6];
  inv[7] = a[0]*a[6]*a[11]-a[0]*a[7]*a[10]-a[4]*a[2]*a[11]+a[4]*a[3]*a[10]+a[8]*a[2]*a[7]-a[8]*a[3]*a[6];
  inv[11] = -a[0]*a[5]*a[11]+a[0]*a[7]*a[9]+a[4]*a[1]*a[11]-a[4]*a[3]*a[9]-a[8]*a[1]*a[7]+a[8]*a[3]*a[5];
  inv[15] = a[0]*a[5]*a[10]-a[0]*a[6]*a[9]-a[4]*a[1]*a[10]+a[4]*a[2]*a[9]+a[8]*a[1]*a[6]-a[8]*a[2]*a[5];
  let det = a[0]*inv[0]+a[1]*inv[4]+a[2]*inv[8]+a[3]*inv[12];
  if (!det) return inv;
  det = 1.0 / det;
  for (let i = 0; i < 16; i++) inv[i] *= det;
  return inv;
}

class Underworld {
  constructor(canvas, status) {
    this.canvas = canvas;
    this.status = status;
    const gl = canvas.getContext('webgl2', { antialias: false, alpha: false });
    if (!gl) throw new Error('WebGL2 is required for volume rendering.');
    this.gl = gl;
    const program = gl.createProgram();
    gl.attachShader(program, compile(gl, gl.VERTEX_SHADER, VERT));
    gl.attachShader(program, compile(gl, gl.FRAGMENT_SHADER, FRAG));
    gl.linkProgram(program);
    if (!gl.getProgramParameter(program, gl.LINK_STATUS)) throw new Error(gl.getProgramInfoLog(program));
    this.program = program;
    gl.useProgram(program);
    // Attribute state is global without vertex array objects, and the survey overlay binds its
    // own attributes, so the quad is rebound at the start of every frame rather than once.
    this.quad = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, this.quad);
    gl.bufferData(gl.ARRAY_BUFFER, new Float32Array([-1, -1, 3, -1, -1, 3]), gl.STATIC_DRAW);
    this.quadLoc = gl.getAttribLocation(program, 'position');
    this.u = {};
    for (const name of ['volume', 'surface', 'invViewProj', 'cameraPos', 'boxMin', 'boxMax',
                        'threshold', 'density', 'depthClip', 'surfaceMix', 'repeatFrac',
                        'showRepeat', 'steps', 'voidBand', 'showVoids']) {
      this.u[name] = gl.getUniformLocation(program, name);
    }
    this.volumeTex = gl.createTexture();
    this.surfaceTex = gl.createTexture();
    gl.uniform1i(this.u.volume, 0);
    gl.uniform1i(this.u.surface, 1);

    this.settings = { threshold: 0.42, density: 0.30, depthClip: 1.0, surfaceMix: 0.85,
                      exaggeration: 1.4, showRepeat: true, showVoids: true, showSurvey: true,
                      steps: 256 };
    this.camera = { azimuth: -0.65, elevation: 0.42, distance: 1.85 };
    this.header = null;
    this.attachControls();
    this.resize();
    window.addEventListener('resize', () => { this.resize(); this.draw(); });
  }

  attachControls() {
    const c = this.canvas;
    let dragging = false, lx = 0, ly = 0;
    const down = (e) => { dragging = true; lx = e.clientX; ly = e.clientY; c.setPointerCapture(e.pointerId); };
    const move = (e) => {
      if (!dragging) return;
      this.camera.azimuth -= (e.clientX - lx) * 0.008;
      this.camera.elevation = Math.max(-1.45, Math.min(1.45, this.camera.elevation + (e.clientY - ly) * 0.008));
      lx = e.clientX; ly = e.clientY;
      this.draw();
    };
    const up = (e) => { dragging = false; try { c.releasePointerCapture(e.pointerId); } catch (_) {} };
    c.addEventListener('pointerdown', down);
    c.addEventListener('pointermove', move);
    c.addEventListener('pointerup', up);
    c.addEventListener('pointercancel', up);
    c.addEventListener('wheel', (e) => {
      e.preventDefault();
      this.camera.distance = Math.max(1.15, Math.min(7, this.camera.distance * (1 + e.deltaY * 0.0012)));
      this.draw();
    }, { passive: false });
    c.tabIndex = 0;
    c.addEventListener('keydown', (e) => {
      const step = 0.09;
      if (e.key === 'ArrowLeft') this.camera.azimuth -= step;
      else if (e.key === 'ArrowRight') this.camera.azimuth += step;
      else if (e.key === 'ArrowUp') this.camera.elevation = Math.min(1.45, this.camera.elevation + step);
      else if (e.key === 'ArrowDown') this.camera.elevation = Math.max(-1.45, this.camera.elevation - step);
      else if (e.key === '+' || e.key === '=') this.camera.distance = Math.max(1.15, this.camera.distance - 0.2);
      else if (e.key === '-') this.camera.distance = Math.min(7, this.camera.distance + 0.2);
      else return;
      e.preventDefault();
      this.draw();
    });
  }

  resize() {
    const dpr = Math.min(window.devicePixelRatio || 1, 2);
    const w = Math.max(320, this.canvas.clientWidth), h = Math.max(280, this.canvas.clientHeight);
    this.canvas.width = Math.round(w * dpr);
    this.canvas.height = Math.round(h * dpr);
    this.gl.viewport(0, 0, this.canvas.width, this.canvas.height);
  }

  async load(bundle) {
    const gl = this.gl;
    this.status.textContent = `Loading ${bundle.label}…`;
    const base = `data/voxels/${bundle.site}/${bundle.patch}`;
    const [header, buffer, image] = await Promise.all([
      fetch(`${base}.json`).then((r) => r.json()),
      fetch(`${base}.bin`).then((r) => r.arrayBuffer()),
      new Promise((res, rej) => { const i = new Image(); i.onload = () => res(i); i.onerror = rej; i.src = `${base}.png`; }),
    ]);
    const [nx, ny, nz] = header.shape;
    const data = new Uint8Array(buffer);
    if (data.length !== nx * ny * nz) throw new Error(`bundle size ${data.length} does not match shape ${header.shape}`);
    gl.activeTexture(gl.TEXTURE0);
    gl.bindTexture(gl.TEXTURE_3D, this.volumeTex);
    gl.pixelStorei(gl.UNPACK_ALIGNMENT, 1);
    gl.texImage3D(gl.TEXTURE_3D, 0, gl.R8, nx, ny, nz, 0, gl.RED, gl.UNSIGNED_BYTE, data);
    for (const p of [gl.TEXTURE_WRAP_S, gl.TEXTURE_WRAP_T, gl.TEXTURE_WRAP_R]) gl.texParameteri(gl.TEXTURE_3D, p, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_3D, gl.TEXTURE_MIN_FILTER, gl.LINEAR);
    gl.texParameteri(gl.TEXTURE_3D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
    gl.activeTexture(gl.TEXTURE1);
    gl.bindTexture(gl.TEXTURE_2D, this.surfaceTex);
    gl.texImage2D(gl.TEXTURE_2D, 0, gl.R8, gl.RED, gl.UNSIGNED_BYTE, image);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
    this.header = header;
    this.setSurvey(this.surveyBoxes(header));
    this.status.textContent = '';
    this.draw();
    return header;
  }

  /* Map surveyed features from metres into the box frame used by the renderer. */
  surveyBoxes(header) {
    const survey = header.survey;
    if (!survey || !survey.features || !survey.features.length) return [];
    const ex = header.extent_m;
    const wide = Math.max(ex.along_azimuth, ex.across_slant_range);
    const fit = Math.min(1, wide / Math.max(ex.depth, 1e-6) * 1.6);
    const sy = (ex.depth / wide) * this.settings.exaggeration * fit;
    return survey.features.map((f) => ({
      kind: f.kind,
      name: f.name,
      centre: [f.centre_m[0] / wide,
               -(f.depth_below_apex_m / ex.depth) * sy,
               f.centre_m[1] / wide],
      size: [f.size_m[0] / wide,
             (f.size_m[2] / ex.depth) * sy,
             f.size_m[1] / wide],
    })).filter((b) => Math.abs(b.centre[1]) <= sy * 1.02);
  }

  bindQuad() {
    const gl = this.gl;
    gl.useProgram(this.program);
    gl.bindBuffer(gl.ARRAY_BUFFER, this.quad);
    gl.enableVertexAttribArray(this.quadLoc);
    gl.vertexAttribPointer(this.quadLoc, 2, gl.FLOAT, false, 0, 0);
  }

  draw() {
    const gl = this.gl, h = this.header;
    this.bindQuad();
    gl.clearColor(0.043, 0.047, 0.055, 1);
    gl.clear(gl.COLOR_BUFFER_BIT);
    if (!h) return;
    const ex = h.extent_m;
    const wide = Math.max(ex.along_azimuth, ex.across_slant_range);
    const sx = ex.along_azimuth / wide, sz = ex.across_slant_range / wide;
    // A deep box (the monument bundles reach past 200 m) would otherwise fill the view; damp
    // the exaggeration rather than making the user find the slider.
    const fit = Math.min(1, wide / Math.max(ex.depth, 1e-6) * 1.6);
    const sy = (ex.depth / wide) * this.settings.exaggeration * fit;
    const boxMin = [-sx / 2, -sy, -sz / 2], boxMax = [sx / 2, 0, sz / 2];
    const { azimuth, elevation, distance } = this.camera;
    const eye = [distance * Math.cos(elevation) * Math.sin(azimuth),
                 distance * Math.sin(elevation) + sy * 0.15,
                 distance * Math.cos(elevation) * Math.cos(azimuth)];
    const view = lookAt(eye, [0, -sy / 2, 0], [0, 1, 0]);
    const proj = perspective(0.85, this.canvas.width / this.canvas.height, 0.05, 40);
    gl.uniformMatrix4fv(this.u.invViewProj, false, invert(mat4Multiply(proj, view)));
    gl.uniform3fv(this.u.cameraPos, eye);
    gl.uniform3fv(this.u.boxMin, boxMin);
    gl.uniform3fv(this.u.boxMax, boxMax);
    gl.uniform1f(this.u.threshold, this.settings.threshold);
    gl.uniform1f(this.u.density, this.settings.density);
    gl.uniform1f(this.u.depthClip, this.settings.depthClip);
    gl.uniform1f(this.u.surfaceMix, this.settings.surfaceMix);
    gl.uniform1f(this.u.repeatFrac, h.depth_axis.repeat_period_m / ex.depth);
    gl.uniform1i(this.u.showRepeat, this.settings.showRepeat ? 1 : 0);
    const kv = h.known_voids;
    gl.uniform2f(this.u.voidBand, kv ? kv.depth_min_m / ex.depth : 0, kv ? kv.depth_max_m / ex.depth : 0);
    gl.uniform1i(this.u.showVoids, this.settings.showVoids && kv ? 1 : 0);
    gl.uniform1i(this.u.steps, this.settings.steps);
    gl.activeTexture(gl.TEXTURE0); gl.bindTexture(gl.TEXTURE_3D, this.volumeTex);
    gl.activeTexture(gl.TEXTURE1); gl.bindTexture(gl.TEXTURE_2D, this.surfaceTex);
    gl.drawArrays(gl.TRIANGLES, 0, 3);
    this.setSurvey(this.surveyBoxes(h));                      // exaggeration may have changed
    this.drawSurvey(mat4Multiply(proj, view), sy, boxMin, boxMax);
  }
}

window.Underworld = Underworld;

/* ---------------------------------------------------------------- survey overlay
 * Independently known subsurface geometry, drawn as wireframe boxes over the volume. These
 * come from excavation surveys and muon imaging, never from the radar, and the wireframe is
 * deliberate: it is a model, not a measurement.
 */
const LINE_VERT = `#version 300 es
in vec3 position;
in vec3 colour;
uniform mat4 viewProj;
out vec3 vColour;
void main() { vColour = colour; gl_Position = viewProj * vec4(position, 1.0); }`;

const LINE_FRAG = `#version 300 es
precision highp float;
in vec3 vColour;
out vec4 frag;
uniform float alpha;
void main() { frag = vec4(vColour, alpha); }`;

const KIND_COLOUR = {
  chamber: [0.30, 0.88, 0.62],
  passage: [0.42, 0.72, 0.98],
  void: [1.00, 0.70, 0.25],
};

function boxEdges(centre, size, colour, out) {
  const [cx, cy, cz] = centre, [sx, sy, sz] = size;
  const hx = sx / 2, hy = sy / 2, hz = sz / 2;
  const corner = (i, j, k) => [cx + (i ? hx : -hx), cy + (j ? hy : -hy), cz + (k ? hz : -hz)];
  const pairs = [
    [[0,0,0],[1,0,0]], [[1,0,0],[1,1,0]], [[1,1,0],[0,1,0]], [[0,1,0],[0,0,0]],
    [[0,0,1],[1,0,1]], [[1,0,1],[1,1,1]], [[1,1,1],[0,1,1]], [[0,1,1],[0,0,1]],
    [[0,0,0],[0,0,1]], [[1,0,0],[1,0,1]], [[1,1,0],[1,1,1]], [[0,1,0],[0,1,1]],
  ];
  for (const [a, b] of pairs) {
    for (const p of [corner(...a), corner(...b)]) {
      out.position.push(p[0], p[1], p[2]);
      out.colour.push(colour[0], colour[1], colour[2]);
    }
  }
}

Underworld.prototype.initSurvey = function () {
  const gl = this.gl;
  const program = gl.createProgram();
  gl.attachShader(program, compile(gl, gl.VERTEX_SHADER, LINE_VERT));
  gl.attachShader(program, compile(gl, gl.FRAGMENT_SHADER, LINE_FRAG));
  gl.linkProgram(program);
  if (!gl.getProgramParameter(program, gl.LINK_STATUS)) throw new Error(gl.getProgramInfoLog(program));
  this.lineProgram = program;
  this.lineBuffers = { position: gl.createBuffer(), colour: gl.createBuffer() };
  this.lineCount = 0;
  this.lineUniforms = { viewProj: gl.getUniformLocation(program, 'viewProj'),
                        alpha: gl.getUniformLocation(program, 'alpha') };
};

/* `features` are in metres relative to the patch centre: x along azimuth, y across range,
 * z positive downward from the surface, matching the volume box. */
Underworld.prototype.setSurvey = function (features) {
  const gl = this.gl;
  if (!this.lineProgram) this.initSurvey();
  const out = { position: [], colour: [] };
  for (const f of features || []) {
    boxEdges(f.centre, f.size, KIND_COLOUR[f.kind] || [1, 1, 1], out);
  }
  this.lineCount = out.position.length / 3;
  this.surveyFeatures = features || [];
  gl.bindBuffer(gl.ARRAY_BUFFER, this.lineBuffers.position);
  gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(out.position), gl.STATIC_DRAW);
  gl.bindBuffer(gl.ARRAY_BUFFER, this.lineBuffers.colour);
  gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(out.colour), gl.STATIC_DRAW);
};

Underworld.prototype.drawSurvey = function (viewProj, scale, boxMin, boxMax) {
  const gl = this.gl;
  if (!this.lineCount || !this.settings.showSurvey) return;
  gl.useProgram(this.lineProgram);
  const bind = (name, buffer) => {
    const loc = gl.getAttribLocation(this.lineProgram, name);
    gl.bindBuffer(gl.ARRAY_BUFFER, buffer);
    gl.enableVertexAttribArray(loc);
    gl.vertexAttribPointer(loc, 3, gl.FLOAT, false, 0, 0);
  };
  bind('position', this.lineBuffers.position);
  bind('colour', this.lineBuffers.colour);
  gl.lineWidth(1);
  gl.uniformMatrix4fv(this.lineUniforms.viewProj, false, viewProj);
  gl.uniform1f(this.lineUniforms.alpha, 0.95);
  gl.enable(gl.BLEND);
  gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA);
  gl.drawArrays(gl.LINES, 0, this.lineCount);
  gl.disable(gl.BLEND);
  for (const name of ['position', 'colour']) {
    gl.disableVertexAttribArray(gl.getAttribLocation(this.lineProgram, name));
  }
  this.bindQuad();
};
