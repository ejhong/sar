/* Progressive enhancement: every experiment, figure and download works without JS. */
const one = selector => document.querySelector(selector);
const all = selector => [...document.querySelectorAll(selector)];

const separation = one('#separation');
const wavelength = one('#wavelength');
function updateScale() {
  const x = Number(separation.value);
  const wave = Number(wavelength.value);
  const depth = x * wave * Math.sin(35 * Math.PI / 180) / (299792458 / 9.65e9);
  one('#separation-value').textContent = `${x.toFixed(1)} m`;
  one('#wavelength-value').textContent = `${wave.toFixed(2)} m`;
  one('#predicted-depth').innerHTML = `${depth.toFixed(1)} <small>m</small>`;
  separation.setAttribute('aria-valuetext', `${x.toFixed(1)} metres`);
  wavelength.setAttribute('aria-valuetext', `${wave.toFixed(2)} metres`);
}
separation.addEventListener('input', updateScale);
wavelength.addEventListener('input', updateScale);
one('#reset-scale').addEventListener('click', () => {
  separation.value = '3'; wavelength.value = '0.48'; updateScale();
});
updateScale();

const scenarioData = JSON.parse(one('#scenario-data').textContent);
function updateScenario() {
  const snr = Number(one('#snr').value);
  const result = scenarioData.find(row => row.snr_db === snr && row.family === one('#scenario').value && row.estimator === one('#estimator').value);
  if (!result) return;
  const percent = value => `${(100 * value).toFixed(0)}%`;
  one('#snr-value').textContent = `${snr} dB`;
  one('#snr').setAttribute('aria-valuetext', `${snr} decibels`);
  one('#scenario-localized').textContent = percent(result.detected_within_20m_rate);
  one('#scenario-false-alarm').textContent = percent(result.false_alarm_rate);
  one('#localized-meter').value = result.detected_within_20m_rate;
  one('#false-alarm-meter').value = result.false_alarm_rate;
  one('#scenario-interval').textContent = `False-alarm interval: ${(100 * result.false_alarm_ci[0]).toFixed(1)}–${(100 * result.false_alarm_ci[1]).toFixed(1)}% (95% Wilson).`;
  one('#scenario-interpretation').textContent = result.false_alarm_rate >= .1
    ? 'Frequent false alarms: a high cavity hit rate cannot establish a useful detector in this condition.'
    : result.detected_within_20m_rate >= .8
      ? 'Strong recovery under these supplied assumptions. This is conditional synthetic performance.'
      : result.detected_within_20m_rate >= .3
        ? 'Partial recovery. Many cavity trials remain missed or incorrectly localized.'
        : 'Low noise-calibrated recovery in this tested condition.';
}
['#snr', '#scenario', '#estimator'].forEach(selector => one(selector).addEventListener('input', updateScenario));
updateScenario();

function revealHash(scroll = false) {
  let hash;
  try { hash = decodeURIComponent(location.hash.slice(1)); } catch { return; }
  const target = document.getElementById(hash);
  if (!target) return;
  if (target.matches('details')) target.open = true;
  for (let parent = target.parentElement; parent; parent = parent.parentElement) {
    if (parent.matches('details')) parent.open = true;
  }
  if (scroll) requestAnimationFrame(() => target.scrollIntoView({block: 'start'}));
}
window.addEventListener('hashchange', () => revealHash(true));
revealHash(Boolean(location.hash));

const experiments = all('.experiment');
const expand = one('#expand-experiments');
function updateExpand() {
  const opened = experiments.every(item => item.open);
  expand.textContent = opened ? 'Collapse all experiments −' : 'Expand all experiments +';
  expand.setAttribute('aria-expanded', String(opened));
}
expand.addEventListener('click', () => {
  const shouldOpen = !experiments.every(item => item.open);
  experiments.forEach(item => { item.open = shouldOpen; });
  updateExpand();
});
experiments.forEach(item => item.addEventListener('toggle', updateExpand));
updateExpand();

const dialog = one('#figure-dialog');
all('[data-zoom]').forEach(link => link.addEventListener('click', event => {
  if (event.ctrlKey || event.metaKey || event.shiftKey || event.altKey || typeof dialog.showModal !== 'function') return;
  event.preventDefault();
  const img = link.querySelector('img');
  one('#enlarged-figure').src = link.href;
  one('#enlarged-figure').alt = img.alt;
  one('#enlarged-caption').textContent = link.closest('figure').querySelector('figcaption').textContent;
  one('#original-figure').href = link.href;
  dialog.showModal();
}));
one('#close-figure').addEventListener('click', () => dialog.close());
dialog.addEventListener('click', event => { if (event.target === dialog) dialog.close(); });

const navLinks = all('.toc a');
const sections = navLinks.map(link => one(link.getAttribute('href'))).filter(Boolean);
let scheduled = false;
function updateNavigation() {
  let active = sections[0];
  for (const section of sections) if (section.getBoundingClientRect().top <= 150) active = section;
  navLinks.forEach(link => {
    if (link.hash === `#${active.id}`) link.setAttribute('aria-current', 'location');
    else link.removeAttribute('aria-current');
  });
  scheduled = false;
}
window.addEventListener('scroll', () => {
  if (!scheduled) { scheduled = true; requestAnimationFrame(updateNavigation); }
}, {passive: true});
updateNavigation();

let printState = [];
window.addEventListener('beforeprint', () => {
  printState = all('details').map(item => [item, item.open]);
  printState.forEach(([item]) => { item.open = true; });
});
window.addEventListener('afterprint', () => printState.forEach(([item, opened]) => { item.open = opened; }));
