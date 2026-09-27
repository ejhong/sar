import {
  CanvasTexture,
  HalfFloatType,
  NoToneMapping,
  PerspectiveCamera,
  Scene,
  SRGBColorSpace,
  Timer,
  Vector2,
  Vector3,
  WebGLRenderer,
  WebGLRenderTarget,
} from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { EffectComposer } from 'three/addons/postprocessing/EffectComposer.js';
import { RenderPass } from 'three/addons/postprocessing/RenderPass.js';
import { UnrealBloomPass } from 'three/addons/postprocessing/UnrealBloomPass.js';
import { OutputPass } from 'three/addons/postprocessing/OutputPass.js';
import { THEMES, type SceneTheme, type ThemeName } from './theme';

export interface FrameInfo {
  /** Seconds since the previous frame (clamped). */
  dt: number;
  /** Seconds since start. */
  time: number;
}

export interface CameraPose {
  position: [number, number, number];
  target: [number, number, number];
  fov?: number;
}

/** A soft radial ground behind the scene: lighter at the centre, as under a lamp. */
function vignette(inner: string, outer: number): CanvasTexture {
  const c = document.createElement('canvas');
  c.width = c.height = 256;
  const g = c.getContext('2d')!;
  const grad = g.createRadialGradient(128, 118, 10, 128, 128, 190);
  grad.addColorStop(0, inner);
  grad.addColorStop(1, '#' + outer.toString(16).padStart(6, '0'));
  g.fillStyle = grad;
  g.fillRect(0, 0, 256, 256);
  const t = new CanvasTexture(c);
  t.colorSpace = SRGBColorSpace;
  return t;
}

const easeInOut = (t: number) => (t < 0.5 ? 4 * t * t * t : 1 - Math.pow(-2 * t + 2, 3) / 2);

/**
 * Renderer, camera and loop (after the knots engine). Deliberately small:
 * every visual layer is a plain Three.js object added to `scene`, and every
 * system that animates registers an `onFrame` callback.
 */
export class Engine {
  readonly renderer: WebGLRenderer;
  readonly scene = new Scene();
  readonly camera: PerspectiveCamera;
  readonly controls: OrbitControls;
  theme: SceneTheme;

  readonly composer: EffectComposer;
  readonly bloom: UnrealBloomPass;
  private timer = new Timer();
  private frameCallbacks = new Set<(f: FrameInfo) => void>();
  private themeCallbacks = new Set<(t: SceneTheme) => void>();
  private running = false;
  private flight?: {
    from: CameraPose;
    to: CameraPose;
    t: number;
    duration: number;
    done?: () => void;
  };
  private resizeObserver: ResizeObserver;
  private idleTime = 0;
  /**
   * Shifts the rendered view by a fraction of the canvas (e.g. x = 0.18
   * places the block right of centre) without changing the perspective.
   */
  viewShift = { x: 0, y: 0 };
  /** Gentle turntable when the viewer is left alone. */
  autoRotate = false;
  autoRotateSpeed = 0.06;

  constructor(readonly canvas: HTMLCanvasElement, themeName: ThemeName = 'night') {
    this.renderer = new WebGLRenderer({
      canvas,
      antialias: true,
      alpha: false,
      powerPreference: 'high-performance',
      preserveDrawingBuffer: false,
    });
    this.renderer.outputColorSpace = SRGBColorSpace;
    this.renderer.toneMapping = NoToneMapping;
    this.renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));

    this.camera = new PerspectiveCamera(30, 1, 0.01, 60);
    this.camera.position.set(2.6, 1.9, 3.4);

    this.controls = new OrbitControls(this.camera, canvas);
    this.controls.enableDamping = true;
    this.controls.dampingFactor = 0.07;
    this.controls.rotateSpeed = 0.55;
    this.controls.zoomSpeed = 0.7;
    this.controls.panSpeed = 0.6;
    this.controls.minDistance = 0.4;
    this.controls.maxDistance = 12;
    this.controls.target.set(0, -0.2, 0);
    this.renderer.localClippingEnabled = true;
    this.controls.addEventListener('start', () => {
      this.idleTime = 0;
      this.flight = undefined;
    });

    const target = new WebGLRenderTarget(1, 1, { type: HalfFloatType, samples: 4 });
    this.composer = new EffectComposer(this.renderer, target);
    this.composer.addPass(new RenderPass(this.scene, this.camera));
    this.bloom = new UnrealBloomPass(new Vector2(256, 256), 0.55, 0.5, 0.62);
    this.composer.addPass(this.bloom);
    this.composer.addPass(new OutputPass());

    this.theme = THEMES[themeName];
    this.applyTheme();

    this.resizeObserver = new ResizeObserver(() => this.resize());
    this.resizeObserver.observe(canvas.parentElement ?? canvas);
    this.resize();
  }

  setTheme(name: ThemeName) {
    this.theme = THEMES[name];
    this.applyTheme();
  }

  private applyTheme() {
    (this.scene.background as { dispose?: () => void } | null)?.dispose?.();
    this.scene.background = vignette(this.theme.backgroundInner, this.theme.background);
    this.bloom.enabled = this.theme.glow > 0;
    this.bloom.strength = this.theme.glow;
    for (const cb of this.themeCallbacks) cb(this.theme);
  }

  onTheme(cb: (t: SceneTheme) => void) {
    this.themeCallbacks.add(cb);
    cb(this.theme);
    return () => this.themeCallbacks.delete(cb);
  }

  onFrame(cb: (f: FrameInfo) => void) {
    this.frameCallbacks.add(cb);
    return () => this.frameCallbacks.delete(cb);
  }

  resize() {
    const el = this.canvas.parentElement ?? this.canvas;
    const w = Math.max(1, el.clientWidth);
    const h = Math.max(1, el.clientHeight);
    this.renderer.setSize(w, h, false);
    this.composer.setPixelRatio(this.renderer.getPixelRatio());
    this.composer.setSize(w, h);
    this.camera.aspect = w / h;
    if (this.viewShift.x || this.viewShift.y)
      this.camera.setViewOffset(w, h, -this.viewShift.x * w, this.viewShift.y * h, w, h);
    else this.camera.clearViewOffset();
    this.camera.updateProjectionMatrix();
  }

  setViewShift(x: number, y = 0) {
    this.viewShift = { x, y };
    this.resize();
  }

  get pose(): CameraPose {
    const p = this.camera.position;
    const t = this.controls.target;
    return { position: [p.x, p.y, p.z], target: [t.x, t.y, t.z], fov: this.camera.fov };
  }

  /** Animated camera move; resolves when it lands. */
  flyTo(to: CameraPose, duration = 1.6): Promise<void> {
    return new Promise((resolve) => {
      if (duration <= 0) {
        this.setPose(to);
        resolve();
        return;
      }
      this.flight = { from: this.pose, to, t: 0, duration, done: resolve };
    });
  }

  setPose(p: CameraPose) {
    this.camera.position.set(...p.position);
    this.controls.target.set(...p.target);
    if (p.fov) {
      this.camera.fov = p.fov;
      this.camera.updateProjectionMatrix();
    }
    this.controls.update();
  }

  start() {
    if (this.running) return;
    this.running = true;
    this.timer.connect(document);
    this.renderer.setAnimationLoop((t) => this.frame(t));
  }

  stop() {
    this.running = false;
    this.renderer.setAnimationLoop(null);
  }

  private frame(timestamp?: number) {
    this.timer.update(timestamp);
    const dt = Math.min(this.timer.getDelta(), 1 / 20);
    const time = this.timer.getElapsed();

    if (this.flight) {
      const f = this.flight;
      f.t = Math.min(1, f.t + dt / f.duration);
      const k = easeInOut(f.t);
      const lerp = (a: number[], b: number[]) => a.map((x, i) => x + (b[i] - x) * k) as [number, number, number];
      // Arc the camera around the target rather than cutting through the ground.
      const target = lerp(f.from.target, f.to.target);
      const fromOff = new Vector3(...f.from.position).sub(new Vector3(...f.from.target));
      const toOff = new Vector3(...f.to.position).sub(new Vector3(...f.to.target));
      const r = fromOff.length() + (toOff.length() - fromOff.length()) * k;
      const dir = fromOff.clone().normalize().lerp(toOff.clone().normalize(), k);
      if (dir.lengthSq() < 1e-6) dir.copy(toOff).normalize();
      dir.normalize().multiplyScalar(r);
      this.camera.position.set(target[0] + dir.x, target[1] + dir.y, target[2] + dir.z);
      this.controls.target.set(...target);
      if (f.from.fov && f.to.fov) {
        this.camera.fov = f.from.fov + (f.to.fov - f.from.fov) * k;
        this.camera.updateProjectionMatrix();
      }
      if (f.t >= 1) {
        this.flight = undefined;
        f.done?.();
      }
    } else if (this.autoRotate) {
      this.idleTime += dt;
      if (this.idleTime > 4) {
        const t = this.controls.target;
        const off = this.camera.position.clone().sub(t);
        const a = this.autoRotateSpeed * dt * Math.min(1, (this.idleTime - 4) / 3);
        off.applyAxisAngle(new Vector3(0, 1, 0), a);
        this.camera.position.copy(t).add(off);
      }
    }

    this.controls.update();
    for (const cb of this.frameCallbacks) cb({ dt, time });
    this.composer.render(dt);
  }

  /** Marks user activity (pauses the idle turntable). */
  poke() {
    this.idleTime = 0;
  }

  dispose() {
    this.stop();
    this.resizeObserver.disconnect();
    this.controls.dispose();
    this.composer.dispose();
    this.renderer.dispose();
  }
}
