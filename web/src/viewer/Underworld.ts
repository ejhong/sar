import { AmbientLight, DirectionalLight, Raycaster, TOUCH, Vector2, Vector3, type Mesh } from 'three';
import { Engine } from './engine/Engine';
import { THEMES, type ThemeName } from './engine/theme';
import { loadIndex, loadScene, loadVolume } from './data/load';
import type { Feature, SiteIndexEntry, SiteScene } from './data/types';
import {
  defaultChoice,
  dimensions,
  findChoice,
  methodsOf,
  step,
  supportName,
  viewSet,
  type Choice,
  type Dimension,
  type Instrument,
  type Method,
} from './catalogue';
import { Block } from './scene/Block';
import type { Wavefield } from './scene/Wavefield';
import { featureDepth, fmtM } from './ui/format';

/** What to look with: the ground as it is, the geophones on it, or the satellite over it. */
export type Mode = 'truth' | 'geophones' | 'satellite';

const MODES: Mode[] = ['truth', 'geophones', 'satellite'];
/** Names the lab used before it was arranged by instrument, still honoured in links. */
const LEGACY: Record<string, Mode> = { recovered: 'geophones', waves: 'geophones', instruments: 'geophones', ground: 'truth' };
const INSTRUMENT_NAME: Record<Mode, string> = { truth: 'The ground', geophones: 'Geophones', satellite: 'Satellite' };
const DIM_LABEL: Record<Dimension, string> = { area: 'go to a place', pass: 'the pass', lines: 'the lines', input: 'what went in', support: 'support' };

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
  private methods: Record<Instrument, Method[]> = { geophones: [], satellite: [] };
  /** The display threshold, kept across one method's pictures so a picture and its control are drawn alike. */
  private threshold = { method: '', value: 0.18 };
  /** A mode named in an old link, resolved once the site's methods are known. */
  private legacy?: string;
  /** A link that names one picture opens on its place; otherwise the whole site. */
  private placeOnLoad = false;
  private linesShown: [number, number][] = [];
  private gridShown: import('./data/types').LineGrid[] = [];
  /** Labels on the stage for each place the picture covers, and where each stands. */
  private places: { el: HTMLButtonElement; at: [number, number, number] }[] = [];
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
      const radar = this.block?.radar;
      if (radar && this.mode === 'satellite') {
        radar.tick(f.dt);
        if (radar.hasSensors) this.updateLookLabel();
        this.engine.poke();
      }
      if (++this.scaleTick % 6 === 0) this.updateScaleBar();
      this.updatePlaces();
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
    this.legacy = h?.legacy;
    this.placeOnLoad = !!h?.item;
    await this.show(first, false);
    if (this.opts.hash)
      addEventListener('hashchange', () => {
        const p = parseHash();
        if (!p) return;
        if (p.site !== this.current) void this.show(p.site);
        else if (p.mode && p.mode !== this.mode) this.setMode(p.legacy ?? p.mode, p.item);
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
    block.addRadar();
    if (scene.radar?.image_hidden) block.radar?.setImageVisible(false);
    this.pin(undefined);
    this.applyLayerToggles();
    const cut = this.$<HTMLInputElement>('[data-uw=cut]');
    if (cut) cut.value = '0';
    const ex = this.$<HTMLInputElement>('[data-uw=exaggeration]');
    if (ex) ex.value = String(block.verticalExaggeration);
    this.methods = { geophones: methodsOf(scene, 'geophones'), satellite: methodsOf(scene, 'satellite') };
    this.threshold.method = '';
    this.renderSiteInfo(scene);
    this.markCurrent();
    this.setMode(this.legacy ?? this.mode, this.item[this.mode], false);
    this.legacy = undefined;
    // the whole site, with every place labelled; a link that names one picture opens on its place
    const sel = this.selected()?.choice;
    if (this.placeOnLoad && sel?.focus && sel.radius_m) this.frameOn(sel.focus, sel.radius_m, animate);
    else this.frame(animate);
    this.placeOnLoad = false;
    this.setLoading(false);
  }

  private availableModes(): Mode[] {
    const s = this.scene;
    if (!s) return ['truth'];
    const m: Mode[] = ['truth'];
    if (this.methods.geophones.length) m.push('geophones');
    if (this.methods.satellite.length || s.radar) m.push('satellite');
    return m;
  }

  /** The picture on show: its method and its choice, or nothing. */
  private selected(): { method: Method; choice: Choice } | undefined {
    if (this.mode === 'truth') return undefined;
    return findChoice(this.methods[this.mode], this.item[this.mode]);
  }

  // ---------- modes ----------

  setMode(mode: Mode | string, item?: string, reframe = true, place = false) {
    const b = this.block;
    const s = this.scene;
    if (!b || !s) return;
    let m: Mode = LEGACY[mode] ?? (MODES.includes(mode as Mode) ? (mode as Mode) : 'truth');
    // links from before the lab was arranged by instrument: a radar picture once sat under the instruments
    if (m === 'geophones' && item && findChoice(this.methods.satellite, item)) m = 'satellite';
    if (mode === 'recovered' && !this.methods.geophones.length && this.methods.satellite.length) m = 'satellite';
    if (mode === 'waves' && !findChoice(this.methods.geophones, item)) item = this.methods.geophones.find((x) => x.key === 'waves')?.choices[0]?.id;
    if (m === 'satellite' && s.radar && !Array.isArray(s.radar.volumes) && (item === 'with' || item === 'without')) item = s.radar.volumes[item];
    if (!this.availableModes().includes(m)) m = 'truth';
    const was = this.selected();
    const refit = (m === 'satellite') !== (this.mode === 'satellite') && !!b.radar;
    this.mode = m;
    let choice: Choice | undefined;
    if (m !== 'truth') {
      const methods = this.methods[m];
      if (m === 'satellite' && (item === 'none' || (item === undefined && this.item.satellite === 'none'))) this.item.satellite = 'none';
      else {
        choice = (findChoice(methods, item) ?? findChoice(methods, this.item[m]))?.choice ?? defaultChoice(methods);
        this.item[m] = choice?.id;
      }
    }
    const vol = choice?.kind === 'volume' ? s.volumes.find((v) => v.id === choice!.id) : undefined;
    const chosen = choice ? findChoice(this.methods[m as Instrument], choice.id) : undefined;
    const set = chosen && vol ? viewSet(chosen.method, chosen.choice).filter((c) => c.kind === 'volume') : [];
    b.showVolumes(set.map((c) => c.id));
    // the lines a lab run laid, on the ground, while one of its pictures is shown
    this.linesShown = b.showLines(m === 'satellite' ? chosen?.choice.grid : undefined);
    const gs = m === 'satellite' ? chosen?.choice.grid : undefined;
    this.gridShown = Array.isArray(gs) ? gs : gs ? [gs] : [];
    b.showSurvey(vol?.survey ? (s.surveys?.find((sv) => sv.id === vol.survey) ?? null) : null);
    if (choice?.kind === 'wave') {
      this.wave = b.showWavefield(choice.id);
      this.wave?.setFrame(0);
      if (this.wave) this.wave.playing = true;
    } else {
      b.showWavefield(null);
      this.wave = undefined;
    }
    b.showRadar(m === 'satellite');
    b.setFeatureEmphasis(m === 'geophones' ? 0.35 : 1);
    const sel = this.selected();
    if (sel && vol) {
      if (this.threshold.method !== sel.method.key) this.threshold = { method: sel.method.key, value: sel.choice.threshold };
      b.setVolumeThreshold(this.threshold.value);
    }
    this.renderModeUI();
    this.renderNow();
    this.renderPlaces(set.length > 1 ? set : []);
    this.writeHash();
    // the camera moves when asked to go to a place, or when the satellite comes or goes; otherwise it stays
    if (reframe && !this.opts.tour) {
      if (place && choice?.focus && choice.radius_m) this.frameOn(choice.focus, choice.radius_m);
      else if (refit) this.frame(true);
    }
    void was;
    this.engine.poke();
  }

  /** A label on the stage for each place the picture covers; a click goes there. */
  private renderPlaces(set: Choice[]) {
    const box = this.$('[data-uw=places]');
    this.places = [];
    if (!box) return;
    box.innerHTML = '';
    const inst = this.mode as Instrument;
    for (const c of set) {
      if (!c.focus || !c.area) continue;
      const el = document.createElement('button');
      el.type = 'button';
      el.className = `uw-place${c.control ? ' control' : ''}${c.id === this.item[this.mode] ? ' current' : ''}`;
      el.textContent = c.control ? `${c.area} · control` : c.area;
      el.title = `Go to ${c.area}`;
      el.addEventListener('click', () => this.setMode(inst, c.id, true, true));
      box.appendChild(el);
      this.places.push({ el, at: c.focus });
    }
    this.updatePlaces();
  }

  private updatePlaces() {
    const b = this.block;
    if (!b || !this.places.length) return;
    const canvas = this.engine.canvas;
    for (const p of this.places) {
      const v = b.world(p.at).project(this.engine.camera);
      const off = v.z > 1 || Math.abs(v.x) > 1.05 || Math.abs(v.y) > 1.05;
      p.el.style.display = off ? 'none' : '';
      p.el.style.transform = `translate(${((v.x + 1) / 2) * canvas.clientWidth}px, ${((1 - v.y) / 2) * canvas.clientHeight}px) translate(-50%, -130%)`;
    }
  }

  private renderModeUI() {
    const s = this.scene;
    if (!s || this.opts.tour) return;
    const available = this.availableModes();
    this.$$<HTMLButtonElement>('[data-mode]').forEach((b) => {
      const m = b.dataset.mode as Mode;
      b.setAttribute('aria-pressed', String(m === this.mode));
      b.disabled = !available.includes(m);
      b.title = b.disabled ? `No ${m === 'geophones' ? 'geophone' : 'satellite'} runs on this site yet` : '';
    });
    const ctx = this.$('[data-uw=mode-context]');
    if (!ctx) return;
    if (this.mode === 'truth') {
      const counts = s.features.reduce<Record<string, number>>((a, f) => ((a[f.status] = (a[f.status] ?? 0) + 1), a), {});
      ctx.innerHTML = `<p class="uw-ctx-note">What is really in the ground: ${Object.entries(counts)
        .map(([k, v]) => `<span class="uw-inline"><span class="uw-dot ${k}"></span>${v} ${STATUS_LABEL[k] ?? k}</span>`)
        .join(' ')}. Click one to pin its details. Every picture the instruments make is scored against this.</p>`;
      return;
    }
    this.renderInstrument(ctx, s, this.mode);
  }

  /** An instrument: what it is, then each method as a card; the open card picks where, what went in, and its control. */
  private renderInstrument(ctx: HTMLElement, s: SiteScene, inst: Instrument) {
    const methods = this.methods[inst];
    const sel = this.selected();
    const waves = methods.filter((m) => m.dot === 'wave');
    const pictures = methods.filter((m) => m.dot !== 'wave');
    const intro =
      inst === 'geophones'
        ? `<p class="uw-ctx-note">Geophones on the ground or down boreholes record waves from a hammer, or the ground’s own hum.
            Each method below turns those records into a picture of the rock, drawn in green.</p>`
        : this.passHTML(s);
    ctx.innerHTML = `${intro}
      ${waves.length ? `<div class="uw-subhead">What they record</div>${waves.map((m) => this.methodHTML(m, sel)).join('')}` : ''}
      <div class="uw-subhead">${inst === 'geophones' ? 'What each method recovers' : 'What each method draws from the image'}</div>
      ${pictures.map((m) => this.methodHTML(m, sel)).join('')}
      ${
        inst === 'satellite'
          ? `<p class="uw-more">${sel ? '<button type="button" class="uw-link" data-choice="none">Hide the picture</button> · ' : ''}more methods join as they are built</p>`
          : ''
      }
      ${
        sel?.choice.kind === 'volume'
          ? `<label class="uw-range"><span>Show values above · ${esc(sel.method.quantity)}</span><input type="range" min="0.02" max="0.95" step="0.01" value="${this.threshold.value}" data-uw="threshold" /></label>`
          : ''
      }`;
    ctx.querySelectorAll<HTMLButtonElement>('[data-open]').forEach((el) =>
      el.addEventListener('click', () => {
        const m = methods.find((x) => x.key === el.dataset.open);
        if (!m) return;
        const keep = sel?.method.key === m.key ? sel.choice : undefined;
        this.setMode(inst, (keep ?? m.choices.find((c) => !c.control) ?? m.choices[0]).id);
      }),
    );
    ctx.querySelectorAll<HTMLButtonElement>('[data-dim]').forEach((el) =>
      el.addEventListener('click', () => {
        if (!sel) return;
        const key = el.dataset.dim as Dimension;
        const value = key === 'support' ? Number(el.dataset.value) : el.dataset.value!;
        this.setMode(inst, step(sel.method, sel.choice, key, value).id, true, key === 'area');
      }),
    );
    ctx.querySelector('[data-uw=all-places]')?.addEventListener('click', () => this.frame(true));
    ctx.querySelector('[data-choice=none]')?.addEventListener('click', () => this.setMode(inst, 'none'));
    const th = ctx.querySelector<HTMLInputElement>('[data-uw=threshold]');
    th?.addEventListener('input', () => {
      this.threshold.value = Number(th.value);
      this.block?.setVolumeThreshold(this.threshold.value);
      this.engine.poke();
    });
    if (inst === 'satellite') this.bindPass(ctx);
    if (sel?.choice.kind === 'wave') this.bindTimeline(ctx);
  }

  private methodHTML(m: Method, sel?: { method: Method; choice: Choice }): string {
    const c = sel?.method.key === m.key ? sel.choice : undefined;
    const chips = c
      ? dimensions(m)
          .map(
            (d) =>
              `<div class="uw-chips" role="group" aria-label="${DIM_LABEL[d.key]}">${d.key === 'area' ? '<span class="uw-chips-label">go to</span><button type="button" class="uw-chip" data-uw="all-places">the whole site</button>' : ''}${d.values
                .map((v) => {
                  const control = m.choices.filter((x) => x[d.key] === v).every((x) => x.control);
                  const label = d.key === 'support' ? supportName(v as number) : String(v);
                  const pressed = d.key === 'area' ? false : c[d.key] === v;
                  return `<button type="button" class="uw-chip${control ? ' control' : ''}${d.key === 'area' ? ' place' : ''}" data-dim="${d.key}" data-value="${esc(String(v))}" aria-pressed="${pressed}"${
                    control ? ' title="a control: nothing here for the method to find"' : ''
                  }>${esc(label)}</button>`;
                })
                .join('')}</div>`,
          )
          .join('')
      : '';
    return `<div class="uw-method${c ? ' open' : ''}">
      <button type="button" class="uw-method-head" data-open="${esc(m.key)}" aria-expanded="${!!c}">
        <span class="uw-dot ${m.dot}"></span><span class="uw-method-name">${esc(m.name)}</span><span class="uw-badge">${m.quantity}</span>
      </button>
      <p class="uw-method-line">${esc(m.line)}</p>
      ${
        c
          ? `${chips}
        <p class="uw-ctx-note">${c.control ? '<span class="uw-tag">control</span> ' : ''}${esc(c.sub)}</p>
        ${c.note ? `<p class="uw-ctx-note">${esc(c.note)}</p>` : ''}
        ${c.stats ? `<p class="uw-ctx-note">${esc(c.stats)}</p>` : ''}
        ${m.why ? `<p class="uw-why"><b>Why it looks like this.</b> ${esc(m.why)}</p>` : ''}
        ${
          c.kind === 'wave'
            ? `<div class="uw-timeline">
          <button type="button" class="uw-play" data-uw="play" aria-label="Pause">❚❚</button>
          <input type="range" min="0" max="1" step="0.001" value="0" data-uw="scrub" aria-label="Time" />
          <span class="uw-time" data-uw="time">0.0 ms</span>
        </div>`
            : ''
        }
        ${m.code ? `<p class="uw-run">code <a href="${esc(m.code.url)}" target="_blank" rel="noopener">${esc(m.code.name)}</a> · ${esc(m.code.version)}, unchanged</p>` : ''}
        ${c.run ? `<p class="uw-run">run ${esc(c.run)}</p>` : ''}`
          : ''
      }
    </div>`;
  }

  /** The satellite's pass: what it is, and on a bench the image's virtual sensors beside the true motion. */
  private passHTML(s: SiteScene): string {
    const r = s.radar;
    if (!r)
      return `<p class="uw-ctx-note">No pass is drawn over this bench: its picture is the paper-style pipeline’s, from the real 2025
        image of Khafre, placed where the claim places what it shows and drawn to the claim’s depths.</p>`;
    const a = r.acquisition;
    if (r.sensors) {
      const sv = r.sensors;
      const pct = (g: number) => `${Math.round(g * 100)}%`;
      const reading = this.block?.radar?.reading ?? 'complex';
      return `<p class="uw-ctx-note">An ICEYE dwell on the real Giza geometry: ${a.aperture_s.toFixed(1)} s, ${a.track_km.toFixed(1)} km
          of track, ${Math.round(a.slant_range_km)} km away, looking ${Math.round(a.incidence_deg)}° from straight down (satellite and
          beam not to scale). On the ground lies the image it makes of this shaking desert.</p>
        <details class="uw-sub"><summary>The image’s virtual sensors: the true motion beside what it reports</summary>
          <p class="uw-ctx-note"><span class="uw-inline"><span class="uw-dot sensor"></span>the true motion</span> under each sensor
            (a test wave ${sv.test_wave.wavelength_m} m long, put into the image itself) beside
            <span class="uw-inline"><span class="uw-dot radar"></span>what the image reports</span> there, read by</p>
          <label class="uw-radio"><input type="radio" name="reading" value="complex" ${reading === 'complex' ? 'checked' : ''}/>
            <span><span class="uw-opt-name">Complex correlation</span><span class="uw-opt-sub">the published way: recovers ${pct(sv.gains.complex)} of the motion</span></span></label>
          <label class="uw-radio"><input type="radio" name="reading" value="magnitude" ${reading === 'magnitude' ? 'checked' : ''}/>
            <span><span class="uw-opt-name">Magnitudes</span><span class="uw-opt-sub">the best tracker tried: ${pct(sv.gains.magnitude)} on this ground</span></span></label>
          <div class="uw-timeline">
            <button type="button" class="uw-play" data-uw="look-play" aria-label="Pause">❚❚</button>
            <span class="uw-time" data-uw="look">look 1 of ${sv.looks_s.length}</span>
            <span></span>
          </div>
        </details>`;
    }
    return `<p class="uw-ctx-note">The real ${esc(a.satellite ?? 'ICEYE')} pass of ${esc(a.date ?? '')}: ${a.aperture_s.toFixed(1)} s,
        ${a.track_km.toFixed(1)} km of track, ${Math.round(a.slant_range_km)} km away, looking ${Math.round(a.incidence_deg)}° from straight
        down (satellite and beam not to scale). On the ground lies the image it made, resampled onto the terrain; it reaches about
        ${Math.round((r.reach_m ?? 0.3) * 100)} cm into the driest sand, less into rock.</p>
      <label class="uw-check"><input type="checkbox" data-uw="sat-image" ${r.image_hidden ? '' : 'checked'} /><span class="uw-dot radar"></span>Show the image</label>`;
  }

  private bindPass(ctx: HTMLElement) {
    ctx.querySelectorAll<HTMLInputElement>('input[name=reading]').forEach((el) =>
      el.addEventListener('change', () => {
        this.block?.radar?.setReading(el.value as 'complex' | 'magnitude');
        this.engine.poke();
      }),
    );
    const play = ctx.querySelector<HTMLButtonElement>('[data-uw=look-play]');
    play?.addEventListener('click', () => {
      const rd = this.block?.radar;
      if (!rd) return;
      rd.playing = !rd.playing;
      play.textContent = rd.playing ? '❚❚' : '▶';
      play.setAttribute('aria-label', rd.playing ? 'Pause' : 'Play');
    });
    const img = ctx.querySelector<HTMLInputElement>('[data-uw=sat-image]');
    img?.addEventListener('change', () => {
      this.block?.radar?.setImageVisible(img.checked);
      this.engine.poke();
    });
  }

  private bindTimeline(ctx: HTMLElement) {
    const play = ctx.querySelector<HTMLButtonElement>('[data-uw=play]');
    const scrub = ctx.querySelector<HTMLInputElement>('[data-uw=scrub]');
    if (!play || !scrub) return;
    play.addEventListener('click', () => {
      if (!this.wave) return;
      this.wave.playing = !this.wave.playing;
      play.textContent = this.wave.playing ? '❚❚' : '▶';
      play.setAttribute('aria-label', this.wave.playing ? 'Pause' : 'Play');
    });
    scrub.addEventListener('input', () => {
      if (!this.wave) return;
      this.wave.playing = false;
      play.textContent = '▶';
      this.wave.setFrame(Number(scrub.value) * (this.wave.frames - 1));
      this.updateTimeline();
    });
  }

  /** On the stage: what the picture is, in one line, and a legend of only what is drawn. */
  private renderNow() {
    const s = this.scene;
    if (!s) return;
    const sel = this.selected();
    const now = this.$('[data-uw=now]');
    if (now) {
      if (this.mode === 'truth') now.innerHTML = `<span class="uw-dot surveyed"></span><b>The ground</b><span>what is really there</span>`;
      else if (!sel) now.innerHTML = `<span class="uw-dot radar"></span><b>${INSTRUMENT_NAME[this.mode]}</b><span>the pass, no picture</span>`;
      else {
        const { method: m, choice: c } = sel;
        const places = new Set(m.choices.map((x) => x.area)).size;
        const lines = new Set(m.choices.map((x) => x.lines)).size;
        const parts = [
          m.name,
          places > 1 ? undefined : c.area,
          c.pass,
          lines > 1 ? c.lines : undefined,
          c.input !== m.name ? c.input : undefined,
          c.support ? supportName(c.support) : undefined,
        ].filter(Boolean) as string[];
        now.innerHTML = `<span class="uw-dot ${m.dot}"></span><b>${INSTRUMENT_NAME[this.mode]}</b><span>${parts.map(esc).join(' · ')}</span>${
          c.control ? '<span class="uw-tag">control</span>' : ''
        }${c.kind === 'volume' ? `<span class="uw-badge">${m.quantity}</span>` : ''}`;
      }
    }
    const legend = this.$('[data-uw=legend]');
    if (legend) {
      const items: [string, string][] = [];
      if (s.features.some((f) => f.status !== 'claimed')) items.push(['surveyed', 'chamber or void']);
      if (s.features.some((f) => f.status === 'claimed')) items.push(['claimed', 'claimed']);
      if (s.structures.length) items.push(['structure', 'monument']);
      const vol = sel?.choice.kind === 'volume' ? s.volumes.find((v) => v.id === sel.choice.id) : undefined;
      if (this.mode === 'geophones') {
        if (sel?.choice.kind === 'wave') items.push(['wave', 'the wave']);
        if (vol) items.push(['recovered', `what the method recovered: ${vol.quantity}`]);
        if (vol?.survey) items.push(['sensor', 'geophones and sources']);
      }
      if (this.mode === 'satellite') {
        if (s.radar?.sensors) items.push(['sensor', 'true motion']);
        items.push(['radar', vol ? `the satellite’s picture: ${sel!.method.quantity}` : s.radar?.sensors ? 'what the image reports' : 'the satellite']);
        // the lines the run laid: one set, or both for a map made of two layouts
        const sets = this.gridShown.map((gr, i) => ({ gr, of: this.linesShown[i]?.[1] ?? 0, drawn: this.linesShown[i]?.[0] ?? 0 }));
        if (sets.length && sets.every((x) => x.of)) {
          const dir = (gr: { direction: string }) => (gr.direction === 'ns' ? 'north–south' : 'east–west');
          const drawn = sets.reduce((t, x) => t + x.drawn, 0);
          const of = sets.reduce((t, x) => t + x.of, 0);
          items.push(['line', `its ${sets.map((x) => `${x.of} ${dir(x.gr)}`).join(' and ')} lines, ${sets[0].gr.step} m apart${drawn < of ? ` (${drawn} drawn)` : ''}`]);
        }
      }
      legend.innerHTML =
        items.map(([dot, label]) => `<span><span class="uw-dot ${dot}"></span>${esc(label)}</span>`).join('') +
        (vol
          ? `<span class="uw-ramp ${vol.cmap === 'magma' ? 'magma' : this.mode === 'satellite' ? 'radar' : 'recovered'}" title="${esc(vol.quantity)}"><i></i>${vol.cmap === 'magma' ? '0 · fit score · 1' : 'lower · higher'}</span>`
          : '');
    }
  }

  private updateLookLabel() {
    const rd = this.block?.radar;
    const el = this.$('[data-uw=look]');
    if (!rd || !el) return;
    const l = rd.look;
    el.textContent = `look ${Math.floor(l.index) + 1} of ${rd.looks} · ${l.time_s >= 0 ? '+' : '−'}${Math.abs(l.time_s).toFixed(2)} s`;
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
    // with the satellite shown, stand further back and look higher, so its beam and the block both fit
    // on a bench the satellite and its beam are the story; over a real site, the ground and what lies under it
    const sat = this.mode === 'satellite' && !!b.radar;
    const bench = this.scene?.radar?.kind === 'bench';
    const r = ((this.opts.tour ? 5.0 : 5.5) + h * 1.25) * (sat ? (bench ? 1.55 : 1.12) : 1);
    const az = (38 * Math.PI) / 180; // from the south, toward the east
    const el = ((sat ? (bench ? 16 : 24) : 21) * Math.PI) / 180;
    const lift = sat ? (bench ? 0.55 : 0.08) : -h * 0.1;
    const pose = {
      position: [r * Math.cos(el) * Math.sin(az), r * Math.sin(el) + lift, r * Math.cos(el) * Math.cos(az)] as [number, number, number],
      target: [0, lift, 0] as [number, number, number],
      fov: 28,
    };
    if (animate) void this.engine.flyTo(pose, 1.3);
    else this.engine.setPose(pose);
  }

  /** Stand off from one place in the site (site coordinates, metres), from the south-east and a little above. */
  private frameOn(focus: [number, number, number], radius_m: number, animate = true) {
    const b = this.block!;
    const c = b.world(focus);
    const r = b.world([focus[0] + radius_m, focus[1], focus[2]]).distanceTo(c);
    const az = (38 * Math.PI) / 180;
    const el = (22 * Math.PI) / 180;
    const pose = {
      position: [c.x + r * Math.cos(el) * Math.sin(az), c.y + r * Math.sin(el), c.z + r * Math.cos(el) * Math.cos(az)] as [number, number, number],
      target: [c.x, c.y, c.z] as [number, number, number],
      fov: 28,
    };
    if (animate) void this.engine.flyTo(pose, 1.3);
    else this.engine.setPose(pose);
  }

  /** A scale bar measured where the camera looks, along the ground's east-west, in nice metres. */
  private updateScaleBar() {
    const el = this.$('[data-uw=scalebar]');
    const b = this.block;
    if (!el || !b) return;
    const canvas = this.engine.canvas;
    const px = (p: Vector3) => {
      const v = p.clone().project(this.engine.camera);
      return [((v.x + 1) / 2) * canvas.clientWidth, ((1 - v.y) / 2) * canvas.clientHeight];
    };
    const at = new Vector3(...this.engine.pose.target);
    const [cx, cy, z] = b.local(at);
    const perM = b.world([cx + 1, cy, z]).distanceTo(at);
    const span = Math.max(1, this.engine.camera.position.distanceTo(at) / perM / 10);     // a tenth of the camera's distance
    const a = px(at);
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
    else if (e.key === 'g' || e.key === 'v') this.setMode('geophones');
    else if (e.key === 'w') this.setMode('waves');
    else if (e.key === 's') this.setMode('satellite');
    else if (e.key === 'r') this.frame(true);
    else if (e.key === 'Escape') this.pin(undefined);
    else if (e.key === ' ' && (this.wave || this.mode === 'satellite')) {
      e.preventDefault();
      this.$<HTMLButtonElement>(this.mode === 'satellite' ? '[data-uw=look-play]' : '[data-uw=play]')?.click();
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

  /** The sites as a drop-down, benches and real sites apart; under it, one line on the current site. */
  private renderSiteList() {
    const list = this.$('[data-uw=sites]');
    if (!list) return;
    const groups: [string, SiteIndexEntry[]][] = [
      ['Test benches', this.sites.filter((s) => s.kind === 'test')],
      ['Real sites', this.sites.filter((s) => s.kind === 'real')],
    ];
    list.innerHTML = `<div class="uw-select"><select data-uw="site-select" aria-label="Site">${groups
      .filter(([, g]) => g.length)
      .map(([label, g]) => `<optgroup label="${label}">${g.map((s) => `<option value="${s.id}">${esc(s.name)}</option>`).join('')}</optgroup>`)
      .join('')}</select></div>
      <p class="uw-site-meta" data-uw="site-meta"></p>`;
    list.querySelector<HTMLSelectElement>('select')!.addEventListener('change', (e) => void this.show((e.target as HTMLSelectElement).value));
  }

  private markCurrent() {
    const sel = this.$<HTMLSelectElement>('[data-uw=site-select]');
    if (sel) sel.value = this.current;
    const meta = this.$('[data-uw=site-meta]');
    const s = this.sites.find((x) => x.id === this.current);
    if (!meta || !s) return;
    const ins = s.instruments;
    const dots = ins
      ? [
          ins.geophones ? `<span class="uw-inline"><span class="uw-dot recovered"></span>geophones</span>` : '',
          ins.satellite ? `<span class="uw-inline"><span class="uw-dot radar"></span>satellite</span>` : '',
        ].join(' ')
      : '';
    meta.innerHTML = `${esc(firstSentence(s.summary))}.${dots ? ` <span class="uw-site-ins">${dots}</span>` : ''}`;
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

function parseHash(): { site: string; mode?: Mode; item?: string; legacy?: string } | null {
  const h = decodeURIComponent(location.hash.slice(1));
  if (!h) return null;
  const [site, mode, item] = h.split('/');
  const m = LEGACY[mode] ?? (MODES.includes(mode as Mode) ? (mode as Mode) : undefined);
  // an old link's picture resolves in setMode, which knows the site's methods
  return { site, mode: m, item: mode === 'waves' && !item ? undefined : item, legacy: mode in LEGACY ? mode : undefined };
}

function isVisible(o: { visible: boolean; parent: any }): boolean {
  for (let p: any = o; p; p = p.parent) if (!p.visible) return false;
  return true;
}

const esc = (s: string) => s.replace(/[&<>"]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' })[c]!);
const firstSentence = (s: string) => (s.match(/^[^.:]+[.:]?/)?.[0] ?? s).replace(/[.:]$/, '');
const fmtV = (v: number) => (v === 0 ? '0' : Math.round(v).toLocaleString('en-US'));
