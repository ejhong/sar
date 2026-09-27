import { AmbientLight, DirectionalLight, Raycaster, TOUCH, Vector2, Vector3, type Mesh } from 'three';
import { Engine } from './engine/Engine';
import { THEMES, type ThemeName } from './engine/theme';
import { loadIndex, loadScene, loadVolume } from './data/load';
import type { Feature, SiteIndexEntry, SiteScene } from './data/types';
import { Block } from './scene/Block';
import type { Wavefield } from './scene/Wavefield';
import { featureDepth, fmtM } from './ui/format';

export type Mode = 'truth' | 'recovered' | 'waves';

export interface TourStop {
  site: string;
  mode?: Mode;
  item?: string;
  caption: string;
  seconds?: number;
}

export interface UnderworldOptions {
  root: HTMLElement;
  initial?: string;
  /** Embedded in a scrolling page: plain wheel and one-finger drags scroll the page. */
  embedded?: boolean;
  /** Keep the view in the URL hash (full page only). */
  hash?: boolean;
  /** A control-free tour through these stops (the front-page hero). */
  tour?: TourStop[];
}

const STATUS_LABEL: Record<string, string> = {
  truth: 'ground truth',
  surveyed: 'surveyed',
  representative: 'representative',
  claimed: 'claimed',
};

/**
 * The underworld viewer: one engine, many sites. The DOM it binds to is
 * rendered by UnderworldViewer.astro (or Hero.astro for the tour); everything
 * site-specific is filled in here from the exported scene.
 */
export class Underworld {
  readonly engine: Engine;
  private block?: Block;
  private scene?: SiteScene;
  private sites: SiteIndexEntry[] = [];
  private current = '';
  private mode: Mode = 'truth';
  private item: Partial<Record<Mode, string>> = {};
  private themeName: ThemeName = 'night';
  private ray = new Raycaster();
  private pointer = new Vector2(2, 2);
  private hovered?: Feature;
  private pinned?: Feature;
  private wave?: Wavefield;
  private down?: { x: number; y: number };
  private loading = 0;
  private scaleTick = 0;
  private $ = <T extends HTMLElement>(sel: string) => this.opts.root.querySelector<T>(sel);
  private $$ = <T extends HTMLElement>(sel: string) => [...this.opts.root.querySelectorAll<T>(sel)];

  constructor(private opts: UnderworldOptions) {
    const canvas = this.$<HTMLCanvasElement>('canvas')!;
    this.engine = new Engine(canvas, 'night');
    const still = matchMedia('(prefers-reduced-motion: reduce)').matches;
    this.engine.autoRotate = !still;
    this.engine.autoRotateSpeed = opts.tour ? 0.07 : 0.045;
    this.engine.setViewShift(opts.tour ? 0.12 : 0.06, 0.01);
    if (opts.embedded || opts.tour) this.yieldScrollToPage(canvas);
    const sun = new DirectionalLight(0xffffff, 1.6);
    sun.position.set(-3, 4, 2);
    this.engine.scene.add(sun, new AmbientLight(0xffffff, 0.9));
    canvas.addEventListener('pointermove', (e) => {
      const r = canvas.getBoundingClientRect();
      this.pointer.set(((e.clientX - r.left) / r.width) * 2 - 1, -((e.clientY - r.top) / r.height) * 2 + 1);
    });
    canvas.addEventListener('pointerleave', () => this.pointer.set(2, 2));
    canvas.addEventListener('pointerdown', (e) => (this.down = { x: e.clientX, y: e.clientY }));
    canvas.addEventListener('pointerup', (e) => {
      if (!this.down || opts.tour) return;
      const moved = Math.hypot(e.clientX - this.down.x, e.clientY - this.down.y);
      this.down = undefined;
      if (moved < 5) this.pin(this.hovered);
    });
    this.engine.onFrame((f) => {
      if (!opts.tour) this.pick();
      if (this.wave?.mesh.visible) {
        this.wave.tick(f.dt);
        this.updateTimeline();
      }
      if (++this.scaleTick % 6 === 0) this.updateScaleBar();
    });
    if (!opts.tour) this.bindControls();
  }

  /**
   * In a scrolling page the viewer must not steal the scroll: the wheel zooms
   * only with ctrl or ⌘ held (a trackpad pinch sends exactly that), and on
   * touch screens one finger scrolls the page while two turn and zoom.
   */
  private yieldScrollToPage(canvas: HTMLCanvasElement) {
    const c = this.engine.controls;
    c.enableZoom = false;
    canvas.style.touchAction = 'pan-y';
    c.touches = { ONE: -1 as unknown as TOUCH, TWO: TOUCH.DOLLY_ROTATE };
    const hint = this.$('[data-uw=hint]');
    let t = 0;
    canvas.addEventListener(
      'wheel',
      (e) => {
        const zoom = (e.ctrlKey || e.metaKey) && !this.opts.tour;
        c.enableZoom = zoom;
        if (!zoom && hint) {
          hint.hidden = false;
          clearTimeout(t);
          t = window.setTimeout(() => (hint.hidden = true), 1400);
        }
      },
      { capture: true, passive: true },
    );
  }

  async start() {
    this.engine.start();
    this.sites = await loadIndex();
    if (this.opts.tour) {
      (window as any).underworld = this;
      await this.runTour(this.opts.tour);
      return;
    }
    this.renderSiteList();
    const h = this.opts.hash ? parseHash() : null;
    const first = this.sites.find((s) => s.id === h?.site)?.id ?? this.opts.initial ?? this.sites[0].id;
    if (h?.mode) this.mode = h.mode;
    if (h?.item && h.mode) this.item[h.mode] = h.item;
    await this.show(first, false);
    if (this.opts.hash)
      addEventListener('hashchange', () => {
        const p = parseHash();
        if (!p) return;
        if (p.site !== this.current) void this.show(p.site);
        else if (p.mode && p.mode !== this.mode) this.setMode(p.mode, p.item);
      });
    (window as any).viewerReady = true;
    (window as any).underworld = this;
  }

  // ---------- sites ----------

  async show(id: string, animate = true) {
    if (!this.sites.find((s) => s.id === id) || id === this.current) return;
    this.current = id;
    this.setLoading(true);
    const scene = await loadScene(id);
    if (this.current !== id) return;
    const old = this.block;
    const block = new Block(scene, THEMES[this.themeName]);
    block.applyTheme();
    this.block = block;
    this.scene = scene;
    this.engine.scene.add(block.root);
    if (old) {
      this.engine.scene.remove(old.root);
      old.dispose();
    }
    await Promise.all([
      ...scene.volumes.map(async (v) => block.addVolume(v, await loadVolume(scene.id, v.file))),
      ...(scene.wavefields ?? []).map(async (w) => block.addWavefield(w, await loadVolume(scene.id, w.file))),
    ]);
    this.pin(undefined);
    this.applyLayerToggles();
    const cut = this.$<HTMLInputElement>('[data-uw=cut]');
    if (cut) cut.value = '0';
    const ex = this.$<HTMLInputElement>('[data-uw=exaggeration]');
    if (ex) ex.value = String(block.verticalExaggeration);
    this.renderSiteInfo(scene);
    this.markCurrent();
    const available = this.availableModes();
    this.setMode(available.includes(this.mode) ? this.mode : 'truth', this.item[this.mode]);
    this.frame(animate);
    this.setLoading(false);
  }

  private availableModes(): Mode[] {
    const s = this.scene;
    if (!s) return ['truth'];
    const m: Mode[] = ['truth'];
    if (s.volumes.length) m.push('recovered');
    if (s.wavefields?.length) m.push('waves');
    return m;
  }

  // ---------- modes ----------

  setMode(mode: Mode, item?: string) {
    const b = this.block;
    const s = this.scene;
    if (!b || !s) return;
    if (!this.availableModes().includes(mode)) mode = 'truth';
    this.mode = mode;
    const vols = s.volumes;
    const waves = s.wavefields ?? [];
    if (mode === 'recovered') {
      const id = vols.find((v) => v.id === item)?.id ?? this.item.recovered ?? vols[0].id;
      this.item.recovered = id;
      b.showVolume(id);
      const v = vols.find((x) => x.id === id);
      b.showSurvey(s.surveys?.find((sv) => sv.id === v?.survey) ?? null);
    } else {
      b.showVolume(null);
      b.showSurvey(null);
    }
    if (mode === 'waves') {
      const id = waves.find((w) => w.id === item)?.id ?? this.item.waves ?? waves[0].id;
      this.item.waves = id;
      this.wave = b.showWavefield(id);
      this.wave?.setFrame(0);
      if (this.wave) this.wave.playing = true;
    } else {
      b.showWavefield(null);
      this.wave = undefined;
    }
    b.setFeatureEmphasis(mode === 'truth' ? 1 : 0.35);
    this.renderModeUI();
    this.writeHash();
    this.engine.poke();
  }

  private renderModeUI() {
    const s = this.scene;
    if (!s || this.opts.tour) return;
    const available = this.availableModes();
    this.$$<HTMLButtonElement>('[data-mode]').forEach((b) => {
      const m = b.dataset.mode as Mode;
      b.setAttribute('aria-pressed', String(m === this.mode));
      b.disabled = !available.includes(m);
      b.title = b.disabled ? 'No runs on this site yet' : '';
    });
    const ctx = this.$('[data-uw=mode-context]');
    if (!ctx) return;
    if (this.mode === 'truth') {
      const counts = s.features.reduce<Record<string, number>>((a, f) => ((a[f.status] = (a[f.status] ?? 0) + 1), a), {});
      ctx.innerHTML = `<p class="uw-ctx-note">What is really in the ground: ${Object.entries(counts)
        .map(([k, v]) => `<span class="uw-inline"><span class="uw-dot ${k}"></span>${v} ${STATUS_LABEL[k] ?? k}</span>`)
        .join(' ')}. Click one to pin its details.</p>`;
      return;
    }
    if (this.mode === 'recovered') {
      ctx.innerHTML =
        s.volumes
          .map(
            (v) => `<label class="uw-radio"><input type="radio" name="vol" value="${v.id}" ${v.id === this.item.recovered ? 'checked' : ''}/>
          <span><span class="uw-opt-name">${esc(v.label)}</span><span class="uw-opt-sub">${esc(v.caption ?? v.method)}</span></span></label>`,
          )
          .join('') +
        `<label class="uw-range"><span>Show anomalies stronger than</span><input type="range" min="0.02" max="0.9" step="0.01" value="0.18" data-uw="threshold" /></label>`;
      ctx.querySelectorAll<HTMLInputElement>('input[name=vol]').forEach((el) =>
        el.addEventListener('change', () => this.setMode('recovered', el.value)),
      );
      const th = ctx.querySelector<HTMLInputElement>('[data-uw=threshold]')!;
      th.addEventListener('input', () => this.block?.setVolumeThreshold(Number(th.value)));
      this.block?.setVolumeThreshold(Number(th.value));
      return;
    }
    const waves = s.wavefields ?? [];
    ctx.innerHTML =
      waves
        .map(
          (w) => `<label class="uw-radio"><input type="radio" name="wave" value="${w.id}" ${w.id === this.item.waves ? 'checked' : ''}/>
        <span><span class="uw-opt-name">${esc(w.label)}</span><span class="uw-opt-sub">${esc(w.caption)}</span></span></label>`,
        )
        .join('') +
      `<div class="uw-timeline">
        <button type="button" class="uw-play" data-uw="play" aria-label="Pause">❚❚</button>
        <input type="range" min="0" max="1" step="0.001" value="0" data-uw="scrub" aria-label="Time" />
        <span class="uw-time" data-uw="time">0.0 ms</span>
      </div>`;
    ctx.querySelectorAll<HTMLInputElement>('input[name=wave]').forEach((el) =>
      el.addEventListener('change', () => this.setMode('waves', el.value)),
    );
    const play = ctx.querySelector<HTMLButtonElement>('[data-uw=play]')!;
    play.addEventListener('click', () => {
      if (!this.wave) return;
      this.wave.playing = !this.wave.playing;
      play.textContent = this.wave.playing ? '❚❚' : '▶';
      play.setAttribute('aria-label', this.wave.playing ? 'Pause' : 'Play');
    });
    const scrub = ctx.querySelector<HTMLInputElement>('[data-uw=scrub]')!;
    scrub.addEventListener('input', () => {
      if (!this.wave) return;
      this.wave.playing = false;
      play.textContent = '▶';
      this.wave.setFrame(Number(scrub.value) * (this.wave.frames - 1));
      this.updateTimeline();
    });
  }

  private updateTimeline() {
    const w = this.wave;
    if (!w) return;
    const t = this.$('[data-uw=time]');
    if (t) t.textContent = `${w.time_ms.toFixed(1)} ms`;
    const s = this.$<HTMLInputElement>('[data-uw=scrub]');
    if (s && w.playing) s.value = String(Math.min(1, w.time_ms / (w.info.dt_ms * (w.frames - 1))));
  }

  // ---------- the tour ----------

  private async runTour(stops: TourStop[]) {
    const cap = this.$('[data-uw=caption]');
    const still = matchMedia('(prefers-reduced-motion: reduce)').matches;
    let i = 0;
    for (;;) {
      const stop = stops[i % stops.length];
      if (cap) cap.classList.remove('on');
      await this.show(stop.site, i > 0);
      if (stop.mode) this.setMode(stop.mode, stop.item);
      if (cap) {
        cap.innerHTML = stop.caption;
        cap.classList.add('on');
      }
      (window as any).viewerReady = true;
      if (still) return;
      await new Promise((r) => setTimeout(r, (stop.seconds ?? 9) * 1000));
      i++;
      if (document.hidden) await new Promise((r) => document.addEventListener('visibilitychange', r, { once: true }));
    }
  }

  // ---------- camera ----------

  private frame(animate: boolean) {
    const b = this.block!;
    const h = b.worldHeight;
    const r = (this.opts.tour ? 5.0 : 5.5) + h * 1.25;
    const az = (38 * Math.PI) / 180; // from the south, toward the east
    const el = (21 * Math.PI) / 180;
    const pose = {
      position: [r * Math.cos(el) * Math.sin(az), r * Math.sin(el) - h * 0.1, r * Math.cos(el) * Math.cos(az)] as [number, number, number],
      target: [0, -h * 0.1, 0] as [number, number, number],
      fov: 28,
    };
    if (animate) void this.engine.flyTo(pose, 1.3);
    else this.engine.setPose(pose);
  }

  /** A scale bar measured on the ground plane at the block's centre, in nice metres. */
  private updateScaleBar() {
    const el = this.$('[data-uw=scalebar]');
    const b = this.block;
    if (!el || !b) return;
    const { x, y } = b.scene.extent;
    const cx = (x[0] + x[1]) / 2;
    const cy = (y[0] + y[1]) / 2;
    const z = b.zTop;
    const canvas = this.engine.canvas;
    const px = (p: Vector3) => {
      const v = p.clone().project(this.engine.camera);
      return [((v.x + 1) / 2) * canvas.clientWidth, ((1 - v.y) / 2) * canvas.clientHeight];
    };
    const span = (x[1] - x[0]) / 4;
    const a = px(b.world([cx, cy, z]));
    const c = px(b.world([cx + span, cy, z]));
    const ppm = Math.hypot(c[0] - a[0], c[1] - a[1]) / span;
    if (!isFinite(ppm) || ppm <= 0) return;
    const target = 90 / ppm;
    const p10 = Math.pow(10, Math.floor(Math.log10(target)));
    const len = [1, 2, 5, 10].map((k) => k * p10).reduce((best, v) => (Math.abs(v - target) < Math.abs(best - target) ? v : best));
    const bar = el.querySelector<HTMLElement>('i');
    const lab = el.querySelector<HTMLElement>('span');
    if (bar) bar.style.width = `${Math.round(len * ppm)}px`;
    if (lab) lab.textContent = fmtM(len);
  }

  // ---------- controls ----------

  private bindControls() {
    const on = (sel: string, ev: string, fn: (el: HTMLInputElement) => void) => {
      const el = this.$<HTMLInputElement>(sel);
      if (el) el.addEventListener(ev, () => fn(el));
    };
    this.$$<HTMLInputElement>('[data-layer]').forEach((el) => el.addEventListener('change', () => this.applyLayerToggles()));
    on('[data-uw=ground]', 'input', (el) => this.block?.setGroundOpacity(Number(el.value)));
    on('[data-uw=cut]', 'input', (el) => {
      this.block?.setCut(Number(el.value));
      this.engine.poke();
    });
    on('[data-uw=exaggeration]', 'input', (el) => {
      this.block?.setExaggeration(Number(el.value));
      this.updateScaleNote();
    });
    this.$$<HTMLButtonElement>('[data-theme]').forEach((b) =>
      b.addEventListener('click', () => this.setTheme(b.dataset.theme as ThemeName)),
    );
    this.$$<HTMLButtonElement>('[data-mode]').forEach((b) =>
      b.addEventListener('click', () => this.setMode(b.dataset.mode as Mode)),
    );
    this.$('[data-uw=reset]')?.addEventListener('click', () => this.frame(true));
    this.$('[data-uw=unpin]')?.addEventListener('click', () => this.pin(undefined));
    if (this.opts.hash) addEventListener('keydown', (e) => this.key(e));
  }

  private key(e: KeyboardEvent) {
    const t = e.target as HTMLElement;
    if (t && (t.tagName === 'INPUT' || t.tagName === 'TEXTAREA' || t.isContentEditable) && (t as HTMLInputElement).type !== 'range') return;
    if (e.metaKey || e.ctrlKey || e.altKey) return;
    const n = Number(e.key);
    if (n >= 1 && n <= this.sites.length) void this.show(this.sites[n - 1].id);
    else if (e.key === 't') this.setMode('truth');
    else if (e.key === 'v') this.setMode('recovered');
    else if (e.key === 'w') this.setMode('waves');
    else if (e.key === 'r') this.frame(true);
    else if (e.key === 'Escape') this.pin(undefined);
    else if (e.key === ' ' && this.wave) {
      e.preventDefault();
      this.$<HTMLButtonElement>('[data-uw=play]')?.click();
    } else return;
  }

  setTheme(name: ThemeName) {
    this.themeName = name;
    this.engine.setTheme(name);
    this.block?.setTheme(THEMES[name]);
    this.opts.root.dataset.theme = name;
    this.$$<HTMLButtonElement>('[data-theme]').forEach((b) => b.setAttribute('aria-pressed', String(b.dataset.theme === name)));
  }

  private applyLayerToggles() {
    const b = this.block;
    if (!b) return;
    this.$$<HTMLInputElement>('[data-layer]').forEach((el) => {
      const g = b.layers[el.dataset.layer as keyof Block['layers']];
      if (g) g.visible = el.checked;
    });
    const ground = this.$<HTMLInputElement>('[data-uw=ground]');
    if (ground) b.setGroundOpacity(Number(ground.value));
  }

  private writeHash() {
    if (!this.opts.hash || !this.current) return;
    const item = this.mode === 'truth' ? '' : this.item[this.mode] ?? '';
    const h = [this.current, this.mode === 'truth' ? '' : this.mode, item].filter(Boolean).join('/');
    if (location.hash.slice(1) !== h) history.replaceState(null, '', `#${h}`);
  }

  // ---------- picking ----------

  private pick() {
    const b = this.block;
    if (!b || this.pointer.x > 1) return this.setHover(undefined);
    this.ray.setFromCamera(this.pointer, this.engine.camera);
    const visible = b.pickables.filter((m) => isVisible(m));
    const hit = this.ray.intersectObjects(visible, false)[0];
    this.setHover(hit ? ((hit.object as Mesh).userData.feature as Feature) : undefined);
  }

  private setHover(f?: Feature) {
    if (f === this.hovered) return;
    this.hovered = f;
    this.engine.canvas.style.cursor = f ? 'pointer' : '';
    if (!this.pinned) this.renderCard(f, false);
  }

  private pin(f?: Feature) {
    this.pinned = f;
    this.renderCard(f ?? this.hovered, !!f);
  }

  private renderCard(f: Feature | undefined, pinned: boolean) {
    const box = this.$('[data-uw=card]');
    if (!box) return;
    if (!f || !this.scene) {
      box.hidden = true;
      return;
    }
    const s = this.scene;
    const sh = f.shape;
    const size =
      sh.type === 'box'
        ? sh.size.map((v) => fmtM(v)).join(' × ')
        : sh.type === 'cylinder'
          ? `${fmtM(2 * sh.radius)} across, ${fmtM(sh.height)} tall`
          : sh.type === 'sphere'
            ? `${fmtM(2 * sh.radius)} across`
            : '';
    box.hidden = false;
    box.classList.toggle('pinned', pinned);
    box.querySelector('[data-uw=card-body]')!.innerHTML = `
      <div class="uw-card-head"><span class="uw-dot ${f.status}"></span>${esc(f.name)}</div>
      <div class="uw-card-meta">${esc(STATUS_LABEL[f.status] ?? f.status)} · ${esc(f.kind)} · ${featureDepth(s, f)}${size ? ` · ${size}` : ''}${
        f.fill !== 'air' ? ` · filled with ${esc(s.materials[f.fill]?.name ?? f.fill)}` : ''
      }</div>
      ${f.placement === 'approximate' ? '<div class="uw-card-warn">placement approximate</div>' : ''}
      <div class="uw-card-src">${esc(f.source)}</div>
      ${pinned && f.note ? `<div class="uw-card-src">${esc(f.note)}</div>` : ''}`;
    this.engine.poke();
  }

  // ---------- DOM ----------

  private renderSiteList() {
    const list = this.$('[data-uw=sites]');
    if (!list) return;
    const groups: [string, SiteIndexEntry[]][] = [
      ['Test benches', this.sites.filter((s) => s.kind === 'test')],
      ['Real sites', this.sites.filter((s) => s.kind === 'real')],
    ];
    let n = 0;
    list.innerHTML = groups
      .filter(([, g]) => g.length)
      .map(
        ([label, g]) => `
        <div class="uw-group">${label}</div>
        ${g
          .map((s) => {
            n++;
            const runs = s.volumes ? `<span class="uw-site-runs" title="recovered volumes">${s.volumes}</span>` : '';
            return `<button class="uw-site" data-site="${s.id}" type="button">
              <span class="uw-site-key">${n}</span>
              <span class="uw-site-text"><span class="uw-site-name">${esc(s.name)}${runs}</span>
              <span class="uw-site-sub">${esc(firstSentence(s.summary))}</span></span>
            </button>`;
          })
          .join('')}`,
      )
      .join('');
    list.querySelectorAll<HTMLButtonElement>('[data-site]').forEach((b) =>
      b.addEventListener('click', () => void this.show(b.dataset.site!)),
    );
  }

  private markCurrent() {
    this.$$<HTMLButtonElement>('[data-site]').forEach((b) => b.setAttribute('aria-current', String(b.dataset.site === this.current)));
  }

  private renderSiteInfo(s: SiteScene) {
    const title = this.$('[data-uw=title]');
    if (title)
      title.innerHTML = `<span class="uw-kind">${s.kind === 'test' ? 'test bench · simulated ground' : 'real site'}</span>
        <h2>${esc(s.name)}</h2><p>${esc(s.summary)}</p>`;
    const legend = this.$('[data-uw=composition]');
    if (legend) {
      const used = new Set<string>([...s.cover.map((c) => c.material), ...s.strata.map((x) => x.material)]);
      const structural = new Set(s.structures.map((x) => x.material));
      const rows = [...used, ...[...structural].filter((m) => !used.has(m))]
        .map((id) => s.materials[id])
        .filter(Boolean)
        .map((m) => {
          const marks = (p: { status: string }) => (p.status === 'sourced' ? '' : p.status === 'derived' ? '′' : '*');
          return `<div class="uw-mat">
            <span class="uw-swatch" style="background:${m.colour}"></span>
            <span class="uw-mat-name">${esc(m.name)}</span>
            <span class="uw-mat-v">${fmtV(m.vp.value)}${marks(m.vp)} / ${fmtV(m.vs.value)}${marks(m.vs)}</span>
          </div>`;
        })
        .join('');
      legend.innerHTML = `${rows}<div class="uw-mat-foot">P / S wave speed, m/s · ′ derived from a source · * assumed</div>`;
    }
    const notes = this.$('[data-uw=notes]');
    if (notes) notes.innerHTML = s.notes.map((n) => `<p>${esc(n)}</p>`).join('');
    this.updateScaleNote();
  }

  private updateScaleNote() {
    const el = this.$('[data-uw=scale]');
    const b = this.block;
    if (!el || !b) return;
    const { x, y, z } = b.scene.extent;
    el.textContent = `block ${fmtM(x[1] - x[0])} × ${fmtM(y[1] - y[0])} × ${fmtM(b.zTop - z[0])} deep · vertical ×${b.verticalExaggeration.toFixed(2)}`;
  }

  private setLoading(on: boolean) {
    this.loading += on ? 1 : -1;
    const el = this.$('[data-uw=loading]');
    if (el) el.hidden = this.loading <= 0;
  }
}

function parseHash(): { site: string; mode?: Mode; item?: string } | null {
  const h = decodeURIComponent(location.hash.slice(1));
  if (!h) return null;
  const [site, mode, item] = h.split('/');
  return { site, mode: (['truth', 'recovered', 'waves'] as Mode[]).includes(mode as Mode) ? (mode as Mode) : undefined, item };
}

function isVisible(o: { visible: boolean; parent: any }): boolean {
  for (let p: any = o; p; p = p.parent) if (!p.visible) return false;
  return true;
}

const esc = (s: string) => s.replace(/[&<>"]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' })[c]!);
const firstSentence = (s: string) => (s.match(/^[^.:]+[.:]?/)?.[0] ?? s).replace(/[.:]$/, '');
const fmtV = (v: number) => (v === 0 ? '0' : Math.round(v).toLocaleString('en-US'));
