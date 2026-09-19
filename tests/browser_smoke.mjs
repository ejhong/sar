// Optional live check. Requires Node 22+, a served docs/ directory, and headless
// Chrome started with --remote-debugging-port=9231. No npm dependencies.
import assert from 'node:assert/strict';
import fs from 'node:fs/promises';

const endpoint = process.env.SAR_BROWSER_ENDPOINT || 'http://127.0.0.1:9231';
const url = process.env.SAR_SITE_URL || 'http://127.0.0.1:4175/';
const screenshotDir = process.env.SAR_SCREENSHOTS || '/private/tmp/sar-final-screenshots';
await fs.mkdir(screenshotDir, {recursive:true});
const tab = await (await fetch(`${endpoint}/json/new?about:blank`, {method:'PUT'})).json();
const ws = new WebSocket(tab.webSocketDebuggerUrl);
await new Promise((resolve, reject) => { ws.onopen = resolve; ws.onerror = reject; });
let id = 0;
const pending = new Map();
const errors = [];
ws.onmessage = event => {
  const message = JSON.parse(event.data);
  if (message.method === 'Runtime.exceptionThrown') errors.push(message.params.exceptionDetails.text);
  if (message.method === 'Network.responseReceived' && message.params.response.status >= 400) {
    errors.push(`${message.params.response.status} ${message.params.response.url}`);
  }
  if (!message.id || !pending.has(message.id)) return;
  const task = pending.get(message.id); pending.delete(message.id); clearTimeout(task.timeout);
  if (message.error) task.reject(message.error); else task.resolve(message.result);
};
const send = (method, params={}) => new Promise((resolve, reject) => {
  const key = ++id;
  const timeout = setTimeout(() => { pending.delete(key); reject(new Error(`Timed out: ${method}`)); }, 15000);
  pending.set(key, {resolve, reject, timeout}); ws.send(JSON.stringify({id:key,method,params}));
});
const evaluate = async expression => {
  const result = await send('Runtime.evaluate', {expression,returnByValue:true,awaitPromise:true});
  if (result.exceptionDetails) throw new Error(result.exceptionDetails.text);
  return result.result.value;
};
const until = async expression => {
  for(let attempt=0; attempt<50; attempt++) {
    if(await evaluate(expression)) return;
    await new Promise(resolve => setTimeout(resolve,100));
  }
  throw new Error(`Condition did not become true: ${expression}`);
};
const screenshot = async name => {
  await evaluate('new Promise(resolve => requestAnimationFrame(() => requestAnimationFrame(() => resolve(true))))');
  const result = await send('Page.captureScreenshot',{format:'png',captureBeyondViewport:false});
  await fs.writeFile(`${screenshotDir}/${name}.png`,Buffer.from(result.data,'base64'));
};
try {
  await send('Page.enable'); await send('Runtime.enable'); await send('Network.enable');
  await send('Emulation.setDeviceMetricsOverride',{width:1440,height:1000,deviceScaleFactor:1,mobile:false});
  await send('Page.navigate',{url});
  await until("document.readyState === 'complete'");
  await evaluate('document.fonts.ready.then(() => true)');
  assert.equal(await evaluate('document.title'),'Doppler Tomography, Tested');
  assert.equal(await evaluate("parseFloat(document.querySelector('#predicted-depth').textContent)"),26.6);
  await evaluate("document.querySelector('#wavelength').value='0.96'; document.querySelector('#wavelength').dispatchEvent(new Event('input',{bubbles:true}))");
  assert.equal(await evaluate("parseFloat(document.querySelector('#predicted-depth').textContent)"),53.2);
  await evaluate("document.querySelector('#reset-scale').click()");
  assert.equal(await evaluate("document.querySelector('#wavelength').value"),'0.48');
  await evaluate("document.querySelector('#separation').focus()");
  await send('Input.dispatchKeyEvent',{type:'keyDown',key:'ArrowRight',code:'ArrowRight',windowsVirtualKeyCode:39});
  await send('Input.dispatchKeyEvent',{type:'keyUp',key:'ArrowRight',code:'ArrowRight',windowsVirtualKeyCode:39});
  assert.equal(await evaluate("document.querySelector('#separation').value"),'3.1');
  await evaluate("document.querySelector('#reset-scale').click(); document.activeElement.blur(); scrollTo({top:0,behavior:'instant'})");
  await screenshot('desktop-overview');
  await evaluate("location.hash='t03_mechanism'");
  await until("document.querySelector('#t03_mechanism').open");
  await evaluate("document.querySelector('#t03_mechanism [data-zoom]').click()");
  assert.equal(await evaluate("document.querySelector('#figure-dialog').open"),true);
  await send('Input.dispatchKeyEvent',{type:'keyDown',key:'Escape',code:'Escape',windowsVirtualKeyCode:27});
  await send('Input.dispatchKeyEvent',{type:'keyUp',key:'Escape',code:'Escape',windowsVirtualKeyCode:27});
  await until("!document.querySelector('#figure-dialog').open");
  await evaluate("document.querySelector('#expand-experiments').click()");
  assert.equal(await evaluate("[...document.querySelectorAll('.experiment')].every(d=>d.open)"),true);
  await evaluate("document.querySelector('#expand-experiments').click()");
  assert.equal(await evaluate("[...document.querySelectorAll('.experiment')].every(d=>!d.open)"),true);
  const dataStatus = await evaluate("Promise.all([...document.querySelectorAll('a[download]')].map(async a=>{const r=await fetch(a.href); if(a.pathname.endsWith('.json')) await r.json(); else if(!(await r.text()).includes(',')) return false; return r.ok}))");
  assert(dataStatus.length >= 9 && dataStatus.every(Boolean));
  for(const width of [1440,768,390,320]) {
    await send('Emulation.setDeviceMetricsOverride',{width,height:width<600?844:1000,deviceScaleFactor:1,mobile:width<600});
    await evaluate("document.querySelectorAll('details').forEach(d=>d.open=true); document.querySelectorAll('img[src]').forEach(i=>i.loading='eager')");
    await evaluate("Promise.all([...document.querySelectorAll('img[src]')].map(i=>i.decode())).then(()=>true)");
    const dimensions=await evaluate('({viewport:innerWidth,document:document.documentElement.scrollWidth})');
    assert(dimensions.document <= dimensions.viewport, `Overflow at ${width}: ${JSON.stringify(dimensions)}`);
    await evaluate("document.querySelectorAll('details').forEach(d=>d.open=false); scrollTo({top:0,behavior:'instant'})");
    await screenshot(`overview-${width}`);
    await evaluate("document.querySelector('#method').scrollIntoView({behavior:'instant',block:'start'})");
    await screenshot(`method-${width}`);
    await evaluate("document.querySelector('#phase-control').scrollIntoView({behavior:'instant',block:'start'})");
    await screenshot(`signal-${width}`);
  }
  assert.deepEqual(errors,[]);
  console.log(JSON.stringify({passed:true,widths:[1440,768,390,320],checks:['slider input and keyboard','reset','deep link','figure dialog and Escape','expand/collapse','JSON downloads','all images decode','no overflow with all details open','no browser or HTTP errors'],screenshots:screenshotDir},null,2));
} finally {
  await send('Page.close'); ws.close();
}
