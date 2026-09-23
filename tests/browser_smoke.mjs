// Optional live check of the built dashboard and field atlas. Requires Node 22+, a served
// docs/ directory, and headless Chrome on --remote-debugging-port=9231. No npm dependencies.
//   cd docs && python3 -m http.server 4175 &
//   chrome --headless=new --remote-debugging-port=9231 &
//   node tests/browser_smoke.mjs
import assert from 'node:assert/strict';
import fs from 'node:fs/promises';

const endpoint = process.env.SAR_BROWSER_ENDPOINT || 'http://127.0.0.1:9231';
const base = process.env.SAR_SITE_URL || 'http://127.0.0.1:4175/';
const screenshotDir = process.env.SAR_SCREENSHOTS || '/private/tmp/sar-final-screenshots';
await fs.mkdir(screenshotDir, {recursive: true});
const tab = await (await fetch(`${endpoint}/json/new?about:blank`, {method: 'PUT'})).json();
const ws = new WebSocket(tab.webSocketDebuggerUrl);
await new Promise((resolve, reject) => { ws.onopen = resolve; ws.onerror = reject; });
let id = 0;
const pending = new Map();
let errors = [];
ws.onmessage = event => {
  const message = JSON.parse(event.data);
  if (message.method === 'Runtime.exceptionThrown') errors.push(message.params.exceptionDetails.text);
  if (message.method === 'Network.responseReceived' && message.params.response.status >= 400) {
    errors.push(`${message.params.response.status} ${message.params.response.url}`);
  }
  if (message.id && pending.has(message.id)) { pending.get(message.id)(message); pending.delete(message.id); }
};
const send = (method, params = {}) => new Promise(resolve => {
  const message = {id: ++id, method, params};
  pending.set(message.id, resolve);
  ws.send(JSON.stringify(message));
});
const evaluate = async expression => {
  const {result: r} = await send('Runtime.evaluate', {expression, returnByValue: true, awaitPromise: true});
  if (r.exceptionDetails) throw new Error(r.exceptionDetails.text);
  return r.result.value;
};
const until = async (expression, timeout = 15000) => {
  const deadline = Date.now() + timeout;
  while (Date.now() < deadline) {
    if (await evaluate(expression)) return true;
    await new Promise(r => setTimeout(r, 120));
  }
  throw new Error(`Condition did not become true: ${expression}`);
};
const screenshot = async name => {
  await evaluate('new Promise(r => requestAnimationFrame(() => requestAnimationFrame(() => r(true))))');
  const shot = await send('Page.captureScreenshot', {format: 'png', captureBeyondViewport: false});
  await fs.writeFile(`${screenshotDir}/${name}.png`, Buffer.from(shot.result.data, 'base64'));
};

try {
  await send('Page.enable'); await send('Runtime.enable'); await send('Network.enable');
  await send('Network.setCacheDisabled', {cacheDisabled: true});
  await send('Emulation.setDeviceMetricsOverride', {width: 1440, height: 1000, deviceScaleFactor: 1, mobile: false});

  // ---------------- the dashboard ----------------
  await send('Page.navigate', {url: base});
  await until("document.readyState === 'complete'");
  await evaluate('document.fonts.ready.then(() => true)');
  assert.equal(await evaluate('document.title'), 'SAR Depth, Tested');

  const counts = await evaluate(`JSON.stringify({
    stats: document.querySelectorAll('.db-stat').length,
    steps: document.querySelectorAll('.db-step').length,
    cards: document.querySelectorAll('.db-card').length,
    figures: document.querySelectorAll('.db-card img').length})`);
  const n = JSON.parse(counts);
  assert.ok(n.stats >= 4, `expected stat tiles, saw ${n.stats}`);
  assert.ok(n.steps >= 6, `expected argument steps, saw ${n.steps}`);
  assert.ok(n.cards >= 24, `expected result cards, saw ${n.cards}`);
  assert.ok(n.figures >= 20, `expected figures, saw ${n.figures}`);

  // every result section the builder emitted is actually on the page
  for (const t of ['r01_giza_controls', 'r09_coherence', 'r11_array', 'r13_known_voids',
                   'r14_budget', 'r15_modal_gate', 't01_static_pyramid']) {
    assert.equal(await evaluate(`!!document.getElementById('${t}')`), true, `missing ${t}`);
  }

  // the viewer sits above the argument, and WebGL actually paints
  assert.equal(await evaluate(
    "document.getElementById('viewer').getBoundingClientRect().top < " +
    "document.getElementById('argument').getBoundingClientRect().top"), true, 'viewer must lead the page');
  // the volume renders wherever WebGL2 exists; skip those checks on a software-less browser
  const webgl2 = await evaluate("!!document.createElement('canvas').getContext('webgl2')");
  if (webgl2) {
    await until("(() => {const c = document.getElementById('db-canvas'); return c && c.width > 300;})()");
    await until("document.getElementById('db-readout').textContent.trim().length > 0");
    const readout = await evaluate("document.getElementById('db-readout').textContent");
    for (const want of ['ICEYE', 'Voxels', 'Repeats every']) {
      assert.ok(readout.includes(want), `readout should mention ${want}: ${readout}`);
    }
    // the volume controls are live
    await evaluate("const s = document.getElementById('db-threshold');" +
                   "s.value = String(Math.min(Number(s.max), Number(s.value) + Number(s.step || 1) * 3));" +
                   "s.dispatchEvent(new Event('input', {bubbles: true}))");
    assert.notEqual(await evaluate("document.getElementById('db-threshold-v').textContent.trim()"), '',
                    'threshold readout should show a value');
  } else {
    console.log('note: WebGL2 unavailable, viewer rendering not checked');
  }

  await evaluate("scrollTo({top: 0, behavior: 'instant'})");
  await screenshot('dashboard-top');

  // a result card carries its figures and its per-test disclosure opens
  assert.ok(await evaluate("document.querySelectorAll('#r15_modal_gate img').length") >= 1,
            'the modal-gate card should carry figures');
  await evaluate("document.getElementById('r15_modal_gate').scrollIntoView({behavior:'instant'})");
  await evaluate("const d = document.querySelector('#r15_modal_gate details'); if (d) d.open = true;");
  assert.ok((await evaluate("document.getElementById('r15_modal_gate').textContent")).includes('repeat'),
            'the modal-gate card should discuss the repeat');
  await screenshot('dashboard-modal-gate');

  // ---------------- the field atlas ----------------
  errors = errors.filter(e => !e.includes('favicon'));
  const dashboardErrors = errors.slice();
  errors = [];
  await send('Page.navigate', {url: new URL('field.html', base).href});
  await until("document.readyState === 'complete'");
  await until("!document.getElementById('site-content').hidden", 20000);
  assert.equal(await evaluate("!!document.querySelector('.db-mark')"), true,
               'the atlas should wear the dashboard header');
  assert.equal(await evaluate("!!document.getElementById('next-test')"), false,
               'the superseded roadmap section should be gone');
  await until("(() => {const c = document.getElementById('radar-canvas'); return c && c.width > 4;})()");
  assert.ok(await evaluate("document.querySelectorAll('.atlas-toc a').length") >= 3,
            'the atlas should keep its section index');
  await screenshot('field-atlas');

  const all = dashboardErrors.concat(errors.filter(e => !e.includes('favicon')));
  assert.deepEqual(all, [], `page errors: ${all.join(' | ')}`);
  console.log(`OK  stats ${n.stats}  steps ${n.steps}  cards ${n.cards}  figures ${n.figures}`);
  console.log(`screenshots in ${screenshotDir}`);
} finally {
  ws.close();
}
