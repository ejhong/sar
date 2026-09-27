import type { SiteIndexEntry, SiteScene } from './types';

const base = import.meta.env.BASE_URL.replace(/\/$/, '');

export const dataUrl = (path: string) => `${base}/data/${path}`;

export async function loadIndex(): Promise<SiteIndexEntry[]> {
  const r = await fetch(dataUrl('sites/index.json'));
  if (!r.ok) throw new Error(`site index: ${r.status}`);
  return (await r.json()).sites;
}

const scenes = new Map<string, Promise<SiteScene>>();

export function loadScene(id: string): Promise<SiteScene> {
  let p = scenes.get(id);
  if (!p) {
    p = fetch(dataUrl(`sites/${id}/scene.json`)).then((r) => {
      if (!r.ok) throw new Error(`scene ${id}: ${r.status}`);
      return r.json();
    });
    scenes.set(id, p);
  }
  return p;
}

export async function loadVolume(siteId: string, file: string): Promise<Uint8Array> {
  const r = await fetch(dataUrl(`sites/${siteId}/${file}`));
  if (!r.ok) throw new Error(`volume ${file}: ${r.status}`);
  return new Uint8Array(await r.arrayBuffer());
}
