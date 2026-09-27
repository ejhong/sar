import { AmbientLight, DirectionalLight, Raycaster, TOUCH, Vector2, type Mesh } from 'three';
import { Engine } from './engine/Engine';
import { THEMES, type ThemeName } from './engine/theme';
import { loadIndex, loadScene, loadVolume } from './data/load';
import type { Feature, SiteIndexEntry, SiteScene } from './data/types';
import { Block } from './scene/Block';
import { featureDepth, fmtM } from './ui/format';

export interface UnderworldOptions {
  root: HTMLElement;
  initial?: string;
  /** Embedded in a scrolling page: plain wheel and one-finger drags scroll the page. */
  embedded?: boolean;
  /** Update the URL hash when the site changes (full page only). */
  hash?: boolean;
}

const STATUS_LABEL: Record<string, string> = {
  truth: 'ground truth',
  surveyed: 'surveyed',
  representative: 'representative',
  claimed: 'claimed',
};

/**
 * The underworld viewer: one engine, many sites. The DOM it binds to is
 * rendered by UnderworldViewer.astro; everything site-specific is filled in
 * here from the exported scene.
 */
export class Underworld {
  readonly engine: Engine;
  private block?: Block;
  private sites: SiteIndexEntry[] = [];
  private current = '';
  private themeName: ThemeName = 'night';
  private ray = new Raycaster();
  private pointer = new Vector2(2, 2);
  private hovered?: Feature;
  private wave?: import('./scene/Wavefield').Wavefield;
  private $ = <T extends HTMLElement>(sel: string) => this.opts.root.querySelector<T>(sel);
  private loading = 0;

  constructor(private opts: UnderworldOptions) {
    const canvas = this.$<HTMLCanvasElement>('canvas')!;
    this.engine = new Engine(canvas, 'night');
    this.engine.autoRotate = !matchMedia('(prefers-reduced-motion: reduce)').matches;
    this.engine.autoRotateSpeed = 0.045;
    this.engine.setViewShift(0.06, 0.01);
    if (opts.embedded) this.yieldScrollToPage(canvas);
    const sun = new DirectionalLight(0xffffff, 1.6);
    sun.position.set(-3, 4, 2);
    this.engine.scene.add(sun, new AmbientLight(0xffffff, 0.9));
    canvas.addEventListener('pointermove', (e) => {
      const r = canvas.getBoundingClientRect();
      this.pointer.set(((e.clientX - r.left) / r.width) * 2 - 1, -((e.clientY - r.top) / r.height) * 2 + 1);
    });
    canvas.addEventListener('pointerleave', () => this.pointer.set(2, 2));
    this.engine.onFrame((f) => {
      this.pick();
      if (this.wave?.mesh.visible) {
        this.wave.tick(f.dt);
        const el = this.$('[data-uw=wave-t]');
        if (el) el.textContent = `${this.wave.time_ms.toFixed(1)} ms`;
      }
    });
    this.bindControls();
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
        const zoom = e.ctrlKey || e.metaKey;
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
    this.renderSiteList();
    const fromHash = this.opts.hash ? location.hash.slice(1) : '';
    const first = this.sites.find((s) => s.id === fromHash)?.id ?? this.opts.initial ?? this.sites[0].id;
    await this.show(first, false);
    if (this.opts.hash) addEventListener('hashchange', () => this.show(location.hash.slice(1)));
    (window as any).viewerReady = true;
    (window as any).underworld = this;
  }

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
    this.engine.scene.add(block.root);
    if (old) {
      this.engine.scene.remove(old.root);
      old.dispose();
    }
    this.applyLayerToggles();
    await this.loadVolumes(scene, block);
    await this.loadWaves(scene, block);
    const cut = this.$<HTMLInputElement>('[data-uw=cut]');
    if (cut) cut.value = '0';
    const ex = this.$<HTMLInputElement>('[data-uw=exaggeration]');
    if (ex) ex.value = String(block.verticalExaggeration);
    this.renderSiteInfo(scene);
    this.markCurrent();
    this.frame(animate);
    if (this.opts.hash && location.hash.slice(1) !== id) history.replaceState(null, '', `#${id}`);
    this.setLoading(false);
  }

  private async loadVolumes(s: SiteScene, block: Block) {
    const box = this.$('[data-uw=recovered]');
    const sec = box?.closest('section') as HTMLElement | null;
    if (!box || !sec) return;
    sec.hidden = !s.volumes.length;
    if (!s.volumes.length) {
      block.showSurvey(null);
      return;
    }
    await Promise.all(s.volumes.map(async (v) => block.addVolume(v, await loadVolume(s.id, v.file))));
    const first = s.volumes[0];
    box.innerHTML =
      s.volumes
        .map(
          (v, i) => `<label class="uw-radio"><input type="radio" name="vol-${s.id}" value="${v.id}" ${i === 0 ? 'checked' : ''}/>
          <span><span class="uw-site-name">${esc(v.label)}</span><span class="uw-site-sub">${esc(v.caption ?? v.method)}</span></span></label>`,
        )
        .join('') +
      `<label class="uw-radio"><input type="radio" name="vol-${s.id}" value="" /><span><span class="uw-site-name">None</span></span></label>`;
    const pick = (id: string) => {
      block.showVolume(id || null);
      const v = s.volumes.find((x) => x.id === id);
      block.showSurvey(s.surveys?.find((sv) => sv.id === v?.survey) ?? null);
      this.engine.poke();
    };
    box.querySelectorAll<HTMLInputElement>('input').forEach((el) => el.addEventListener('change', () => pick(el.value)));
    pick(first.id);
    const th = this.$<HTMLInputElement>('[data-uw=threshold]');
    if (th) block.setVolumeThreshold(Number(th.value));
  }

  private async loadWaves(s: SiteScene, block: Block) {
    const box = this.$('[data-uw=waves]');
    const sec = box?.closest('section') as HTMLElement | null;
    this.wave = undefined;
    if (!box || !sec) return;
    const list = s.wavefields ?? [];
    sec.hidden = !list.length;
    if (!list.length) return;
    await Promise.all(list.map(async (w) => block.addWavefield(w, await loadVolume(s.id, w.file))));
    box.innerHTML =
      list
        .map(
          (w) => `<label class="uw-radio"><input type="radio" name="wave-${s.id}" value="${w.id}"/>
          <span><span class="uw-site-name">${esc(w.label)}</span><span class="uw-site-sub">${esc(w.caption)}</span></span></label>`,
        )
        .join('') +
      `<label class="uw-radio"><input type="radio" name="wave-${s.id}" value="" checked/><span><span class="uw-site-name">Still</span></span></label>`;
    box.querySelectorAll<HTMLInputElement>('input').forEach((el) =>
      el.addEventListener('change', () => {
        this.wave = block.showWavefield(el.value || null);
        this.engine.poke();
      }),
    );
  }

  private frame(animate: boolean) {
    const b = this.block!;
    const h = b.worldHeight;
    const r = 5.5 + h * 1.25;
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

  // ---------- controls ----------

  private bindControls() {
    const on = (sel: string, ev: string, fn: (el: HTMLInputElement) => void) => {
      const el = this.$<HTMLInputElement>(sel);
      if (el) el.addEventListener(ev, () => fn(el));
    };
    this.opts.root.querySelectorAll<HTMLInputElement>('[data-layer]').forEach((el) =>
      el.addEventListener('change', () => this.applyLayerToggles()),
    );
    on('[data-uw=ground]', 'input', (el) => this.block?.setGroundOpacity(Number(el.value)));
    on('[data-uw=threshold]', 'input', (el) => this.block?.setVolumeThreshold(Number(el.value)));
    on('[data-uw=cut]', 'input', (el) => {
      this.block?.setCut(Number(el.value));
      this.engine.poke();
    });
    on('[data-uw=exaggeration]', 'input', (el) => {
      this.block?.setExaggeration(Number(el.value));
      this.updateScaleNote();
    });
    this.opts.root.querySelectorAll<HTMLButtonElement>('[data-theme]').forEach((b) =>
      b.addEventListener('click', () => this.setTheme(b.dataset.theme as ThemeName)),
    );
    this.$('[data-uw=reset]')?.addEventListener('click', () => this.frame(true));
  }

  setTheme(name: ThemeName) {
    this.themeName = name;
    this.engine.setTheme(name);
    this.block?.setTheme(THEMES[name]);
    this.opts.root.dataset.theme = name;
    this.opts.root.querySelectorAll<HTMLButtonElement>('[data-theme]').forEach((b) =>
      b.setAttribute('aria-pressed', String(b.dataset.theme === name)),
    );
  }

  private applyLayerToggles() {
    const b = this.block;
    if (!b) return;
    this.opts.root.querySelectorAll<HTMLInputElement>('[data-layer]').forEach((el) => {
      const g = b.layers[el.dataset.layer as keyof Block['layers']];
      if (g) g.visible = el.checked;
    });
    const ground = this.$<HTMLInputElement>('[data-uw=ground]');
    if (ground) b.setGroundOpacity(Number(ground.value));
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
    const box = this.$('[data-uw=hover]');
    if (!box) return;
    if (!f) {
      box.hidden = true;
      return;
    }
    const s = this.block!.scene;
    const d = featureDepth(s, f);
    box.hidden = false;
    box.innerHTML = `
      <div class="uw-hover-head"><span class="uw-dot ${f.status}"></span>${esc(f.name)}</div>
      <div class="uw-hover-meta">${esc(STATUS_LABEL[f.status] ?? f.status)} · ${esc(f.kind)} · ${d}${
        f.placement === 'approximate' ? ' · placement approximate' : ''
      }</div>
      <div class="uw-hover-src">${esc(f.source)}</div>`;
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
    list.innerHTML = groups
      .filter(([, g]) => g.length)
      .map(
        ([label, g]) => `
        <div class="uw-group">${label}</div>
        ${g
          .map(
            (s) => `<button class="uw-site" data-site="${s.id}" type="button">
              <span class="uw-site-name">${esc(s.name)}</span>
              <span class="uw-site-sub">${esc(firstSentence(s.summary))}</span>
            </button>`,
          )
          .join('')}`,
      )
      .join('');
    list.querySelectorAll<HTMLButtonElement>('[data-site]').forEach((b) =>
      b.addEventListener('click', () => void this.show(b.dataset.site!)),
    );
  }

  private markCurrent() {
    this.opts.root.querySelectorAll<HTMLButtonElement>('[data-site]').forEach((b) =>
      b.setAttribute('aria-current', String(b.dataset.site === this.current)),
    );
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
      const counts = s.features.reduce<Record<string, number>>((a, f) => ((a[f.status] = (a[f.status] ?? 0) + 1), a), {});
      legend.innerHTML = `${rows}
        <div class="uw-mat-foot">P / S wave speed, m/s · ′ derived · * assumed</div>
        <div class="uw-counts">${Object.entries(counts)
          .map(([k, v]) => `<span><span class="uw-dot ${k}"></span>${v} ${STATUS_LABEL[k] ?? k}</span>`)
          .join('')}</div>`;
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
    el.textContent = `block ${fmtM(x[1] - x[0])} × ${fmtM(y[1] - y[0])} × ${fmtM(b.zTop - z[0])} deep · vertical ×${b.verticalExaggeration.toFixed(
      2,
    )}`;
  }

  private setLoading(on: boolean) {
    this.loading += on ? 1 : -1;
    const el = this.$('[data-uw=loading]');
    if (el) el.hidden = this.loading <= 0;
  }
}

function isVisible(o: { visible: boolean; parent: any }): boolean {
  for (let p: any = o; p; p = p.parent) if (!p.visible) return false;
  return true;
}

const esc = (s: string) => s.replace(/[&<>"]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' })[c]!);
const firstSentence = (s: string) => (s.match(/^[^.:]+[.:]?/)?.[0] ?? s).replace(/[.:]$/, '');
const fmtV = (v: number) => (v === 0 ? '0' : Math.round(v).toLocaleString('en-US'));
