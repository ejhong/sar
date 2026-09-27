/**
 * Screenshot bench for visual work.
 *
 *   npx tsx scripts/shot.ts <path> <out.png> [--w 1400] [--h 900]
 *        [--eval "<js run after load>"] [--wait 1500] [--base http://127.0.0.1:4321]
 *        [--scroll <y>] [--full]
 *
 * Waits for `window.viewerReady` when present. Uses the GPU via ANGLE/Metal
 * where available so WebGL renders at real quality. SHOT_CHROMIUM=<path>
 * launches another Chromium than Playwright's own (cloud sessions set it).
 */
import { chromium } from 'playwright';

const args = process.argv.slice(2);
const flag = (name: string, dflt?: string) => {
  const i = args.indexOf(`--${name}`);
  return i >= 0 ? args[i + 1] : dflt;
};
const has = (name: string) => args.includes(`--${name}`);
const [path, out] = args.filter((a, i) => !a.startsWith('--') && !(i > 0 && args[i - 1].startsWith('--') && !['--full'].includes(args[i - 1])));

async function main() {
  if (!path || !out) throw new Error('usage: shot.ts <path> <out.png>');
  const base = flag('base', 'http://127.0.0.1:4321')!;
  const browser = await chromium.launch({
    headless: true,
    executablePath: process.env.SHOT_CHROMIUM || undefined,
    args: ['--use-angle=metal', '--enable-gpu', '--ignore-gpu-blocklist', '--enable-unsafe-swiftshader'],
  });
  const page = await browser.newPage({
    viewport: { width: Number(flag('w', '1400')), height: Number(flag('h', '900')) },
    deviceScaleFactor: Number(flag('dpr', '1')),
  });
  const logs: string[] = [];
  page.on('console', (m) => logs.push(`[${m.type()}] ${m.text()}`));
  page.on('pageerror', (e) => logs.push(`[pageerror] ${e.message}`));
  await page.goto(base + path, { waitUntil: 'load' });
  await page
    .waitForFunction(() => (window as any).viewerReady === true || !document.querySelector('canvas'), null, {
      timeout: 90_000,
    })
    .catch(() => logs.push('[shot] viewerReady timeout'));
  const code = flag('eval');
  if (code) {
    const r = await page.evaluate(`(async () => { ${code} })()`);
    if (r !== undefined) logs.push(`[eval] ${JSON.stringify(r)}`);
  }
  const scroll = flag('scroll');
  if (scroll) await page.evaluate((y) => window.scrollTo(0, Number(y)), scroll);
  await page.waitForTimeout(Number(flag('wait', '1200')));
  const series = Number(flag('series', '1'));
  if (series > 1) {
    const interval = Number(flag('interval', '700'));
    for (let i = 0; i < series; i++) {
      await page.screenshot({ path: out.replace(/\.png$/, `-${i}.png`), fullPage: has('full') });
      await page.waitForTimeout(interval);
    }
  } else await page.screenshot({ path: out, fullPage: has('full') });
  await browser.close();
  for (const l of logs) console.log(l);
}

main().catch((e) => {
  console.error(e);
  process.exit(1);
});
