import {SurveyViewer} from './model-viewer.js?v=84a96b52f022';

const $=id=>document.getElementById(id);
const esc=value=>String(value).replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const fmt=(v,n=3)=>Number(v).toFixed(n);
const comma=v=>Number(v).toLocaleString('en-US');
const pair=(a,b,n=4)=>`(${fmt(a,n)}, ${fmt(b,n)})`;
const json=async url=>{const r=await fetch(url);if(!r.ok)throw new Error(`Could not load ${url} (${r.status})`);return r.json();};
const setLink=(id,url)=>$(id).href=url;
let sites=[],current=null,currentReport=null,currentSurveys=null,generation=0,acquisitionRequest=0;
const cache=new Map();
const load=path=>{if(!cache.has(path))cache.set(path,json(path));return cache.get(path);};
const surveyViewer=new SurveyViewer($('model-canvas'));
$('model-canvas').addEventListener('viewchange',event=>document.querySelectorAll('[data-model-view]').forEach(b=>b.setAttribute('aria-pressed',String(b.dataset.modelView===event.detail))));

class RadarViewer {
  constructor(canvas){
    this.canvas=canvas;this.ctx=canvas.getContext('2d');this.zoom=1;this.pan=[0,0];this.image=null;this.points=[];this.showPoints=false;this.heightShift=0;this.loadNumber=0;
    new ResizeObserver(()=>this.draw()).observe(canvas);
    let drag=null;
    canvas.addEventListener('pointerdown',e=>{drag=[e.clientX,e.clientY];canvas.setPointerCapture(e.pointerId);canvas.focus({preventScroll:true});});
    canvas.addEventListener('pointermove',e=>{if(!drag)return;this.pan[0]+=e.clientX-drag[0];this.pan[1]+=e.clientY-drag[1];drag=[e.clientX,e.clientY];this.draw();});
    for(const name of ['pointerup','pointercancel','lostpointercapture'])canvas.addEventListener(name,()=>drag=null);
    canvas.addEventListener('wheel',e=>{e.preventDefault();this.zoomBy(Math.exp(-Math.sign(e.deltaY)*.13));},{passive:false});
    canvas.addEventListener('keydown',e=>{
      if(!['+','=','-','Home','ArrowLeft','ArrowRight','ArrowUp','ArrowDown'].includes(e.key))return;
      e.preventDefault();
      if(e.key==='+'||e.key==='=')this.zoomBy(1.3);if(e.key==='-')this.zoomBy(1/1.3);if(e.key==='Home')this.reset();
      if(e.key==='ArrowLeft')this.pan[0]-=25;if(e.key==='ArrowRight')this.pan[0]+=25;
      if(e.key==='ArrowUp')this.pan[1]-=25;if(e.key==='ArrowDown')this.pan[1]+=25;this.draw();
    });
  }
  async setView(view,base,points){
    const number=++this.loadNumber;
    const im=new Image();im.src=base+view.image;await im.decode();
    if(number!==this.loadNumber)return;
    this.view=view;this.image=im;this.points=points||[];this.reset();
    $('image-caption').textContent=view.description;
    setLink('image-download',im.src);
    this.canvas.setAttribute('aria-label',`${view.label}. Real satellite SAR intensity, not a depth image. Native rows ${view.rows.join(' to ')}, columns ${view.cols.join(' to ')}. Drag or arrow keys pan; plus and minus zoom; Home resets.`);
  }
  zoomBy(factor){this.zoom=Math.max(1,Math.min(10,this.zoom*factor));this.draw();}
  reset(){this.zoom=1;this.pan=[0,0];this.draw();}
  draw(){
    if(!this.image)return;
    const c=this.canvas,ctx=this.ctx,w=c.clientWidth,h=c.clientHeight;
    if(!w||!h)return;
    const dpr=Math.min(devicePixelRatio||1,2);
    if(c.width!==Math.round(w*dpr)||c.height!==Math.round(h*dpr)){c.width=Math.round(w*dpr);c.height=Math.round(h*dpr);}
    ctx.setTransform(dpr,0,0,dpr,0,0);ctx.fillStyle='#111c1a';ctx.fillRect(0,0,w,h);
    const scale=Math.min(w/this.image.width,h/this.image.height)*this.zoom;
    const iw=this.image.width*scale,ih=this.image.height*scale;
    this.pan[0]=Math.max(-iw/2,Math.min(iw/2,this.pan[0]));this.pan[1]=Math.max(-ih/2,Math.min(ih/2,this.pan[1]));
    const x=(w-iw)/2+this.pan[0],y=(h-ih)/2+this.pan[1];ctx.drawImage(this.image,x,y,iw,ih);
    const v=this.view;
    if(this.showPoints){
      for(const point of this.points){
        const dh=this.heightShift,delta=point.height_minus_plus_5m_shifts_px[dh<0?0:1];
        const row=point.projected_native_row_col[0]+delta[0]*Math.abs(dh)/5;
        const col=point.projected_native_row_col[1]+delta[1]*Math.abs(dh)/5;
        if(row<v.rows[0]||row>=v.rows[1]||col<v.cols[0]||col>=v.cols[1])continue;
        let py=(row-v.rows[0]+.5)/(v.rows[1]-v.rows[0]);if(v.flip_rows)py=1-py;
        const px=x+(col-v.cols[0]+.5)/(v.cols[1]-v.cols[0])*iw;py=y+py*ih;
        ctx.strokeStyle='#ffc976';ctx.lineWidth=1.5;ctx.beginPath();ctx.moveTo(px-5,py);ctx.lineTo(px+5,py);ctx.moveTo(px,py-5);ctx.lineTo(px,py+5);ctx.stroke();
        if(this.zoom>1.5){ctx.font='10px "IBM Plex Mono",monospace';ctx.lineWidth=3;ctx.strokeStyle='#17251be0';ctx.strokeText(point.id,px+8,py-7);ctx.fillStyle='#ffd99b';ctx.fillText(point.id,px+8,py-7);}
      }
    }
    ctx.fillStyle='#14261bd4';ctx.fillRect(12,12,130,32);ctx.fillStyle='#d9e4d5';ctx.font='9px "IBM Plex Mono",monospace';ctx.textAlign='left';ctx.fillText(`Azimuth ${v.flip_rows?'↑':'↓'}  Range →`,22,32);
    const metresPerPixel=v.approximate_ground_width_m/iw;
    const available=80*metresPerPixel,power=10**Math.floor(Math.log10(available));
    const length=[1,2,5,10].map(n=>n*power).filter(n=>n<=available).at(-1)||power;
    $('radar-scale').style.width=`${length/metresPerPixel}px`;
    $('radar-scale').style.minWidth='0';$('radar-scale').textContent=`≈ ${length<1?fmt(length,1):comma(length)} m`;
    c.dataset.view=v.id;c.dataset.zoom=this.zoom.toFixed(3);c.dataset.landmarks=String(this.showPoints);
  }
}
const radar=new RadarViewer($('radar-canvas'));


function renderChecks(report,base){
  const timing=report.timing,groups=report.translation.groups;
  const pass=groups.filter(g=>g.meets_0_02px_p95_tolerance).length;
  $('check-summary').innerHTML=`<div><span class="eyebrow">Complex samples checked</span><strong>${fmt(report.native_shape[0]*report.native_shape[1]/1e9,2)} billion</strong><small>Every I/Q value read; all finite</small></div><div><span class="eyebrow">Shift-control groups</span><strong>${pass} / ${groups.length}</strong><small>Meet the preset 0.02 px tolerance</small></div><div><span class="eyebrow">Survey matches scored</span><strong>0</strong><small>Registration and depth test still open</small></div>`;
  const projection=report.survey_projection;
  $('registration-title').textContent=projection?'Height datum changes the target position.':'Surface registration needs a survey reference.';
  if(projection){
    $('registration-body').innerHTML=`<p>The independent Glen Dash Foundation survey supplies surface coordinates and elevations above mean sea level. The native radar projection needs height above the WGS84 ellipsoid. A local EGM96 grid gives a provisional conversion of <strong>about +${fmt(projection.points[0].geoid_undulation_m,1)} m</strong>. The survey’s precise vertical tie still needs verification.</p><div class="registration-metrics"><div><span>Omitting the datum conversion</span><strong>${fmt(Math.abs(projection.mean_datum_range_shift_px),0)} pixels</strong><small>Mean range displacement of these ${projection.points.length} projected points</small></div><div><span>Approximate flat-ground equivalent</span><strong>${fmt(Math.abs(projection.approximate_mean_datum_ground_shift_m),0)} m</strong><small>Far larger than the reference chamber widths</small></div></div><p>${esc(projection.interpretation)}</p><p class="small"><a href="${esc(projection.reference)}">Glen Dash Foundation survey ↗</a> · <a href="https://raw.githubusercontent.com/OSGeo/PROJ-data/master/us_nga/us_nga_README.txt">NGA / PROJ geoid grid ↗</a> · <a href="https://geographiclib.sourceforge.io/C++/doc/geoid.html">Height conversion convention ↗</a></p><details><summary>Inspect all projected landmarks <span aria-hidden="true">+</span></summary><div class="table-scroll" tabindex="0" role="region" aria-label="Survey surface points and unfitted native projections"><table class="field-table"><thead><tr><th scope="col">Point</th><th scope="col">Published description</th><th scope="col">MSL height m</th><th scope="col">Projected row / col</th></tr></thead><tbody>${projection.points.map(p=>`<tr><th scope="row">${esc(p.id)}</th><td>${esc(p.description)}</td><td>${fmt(p.orthometric_height_m,3)}</td><td>${pair(...p.projected_native_row_col,1)}</td></tr>`).join('')}</tbody></table></div><p>No manually accepted image correspondence or held-out residual is assigned. These points are not cavity labels.</p></details>`;
  }else{
    $('registration-body').innerHTML='<p>The supplied native geometry is audited, including terrain-height sensitivity and axis conventions. An independently surveyed underground reference with usable coordinates and elevations is not yet registered here.</p><p>Terrain height and radar layover can shift an apparent feature. A clear-looking image does not establish placement accuracy. No target mask, known-negative label or depth error is assigned.</p>';
  }
  $('timing-comparison').innerHTML=`<div><strong>${fmt(timing.collection_duration_s,2)} s</strong><span>First to last transmitted pulse</span></div><div><strong>${fmt(timing.zero_doppler_image_span_s,3)} s</strong><span>Zero-Doppler sweep across image rows</span></div>`;
  $('timing-explanation').textContent=`These describe different aspects of the focused product. Dividing the processed bandwidth by the magnitude of the midpoint Doppler rate gives a nominal ${fmt(timing.nominal_bandwidth_over_rate_s[1],2)} s—not the short row-time span. This consistency check does not yet validate the time axis for extracting motion from this dwell acquisition.`;
  $('spectrum-image').src=base+report.figures.spectrum;
  $('spectrum-mobile').srcset=base+report.figures.spectrum_mobile;
  setLink('spectrum-link',base+report.figures.spectrum);
  setLink('spectrum-download',base+report.figures.spectrum);
  $('spectrum-details').innerHTML=`<p>Acquisition PRF: ${fmt(timing.acquisition_prf_hz,1)} Hz. Processing PRF: ${fmt(timing.processing_prf_hz,1)} Hz. Their difference reflects the focused output sampling; extra image samples are not independent transmitted pulses.</p><p>The Doppler-rate polynomial is evaluated relative to mid-range time, following the <a href="https://github.com/ngageoint/sarpy/blob/master/sarpy/io/complex/iceye.py">NGA SarPy ICEYE reader</a>. Near, centre and far rates are ${timing.doppler_rate_near_center_far_hz_s.map(v=>fmt(v,1)).join(', ')} Hz/s. Evaluating that polynomial at the absolute range delay would be a different calculation.</p><p>Fractions of measured strip power outside a zero-centred metadata band: ${report.spectra.map(s=>fmt(100*s.outside_zero_centered_metadata_band_fraction,2)+'%').join(', ')}. No per-strip recentering, deramp, taper or fitted bandwidth is applied. The scene contains distributed surface returns. These are diagnostic Fourier spectra, not recovered vibration histories.</p><p><a href="https://sar.iceye.com/6.0.0/productFormats/metadata/">ICEYE metadata definitions ↗</a>. The companion SICD dimensions and image plane must also be reconciled with the native HDF5 product before a motion adapter is accepted.</p>`;
  $('translation-state').textContent=pass===groups.length?'Tolerance met':'Limits found';
  $('translation-state').className=`check-state ${pass===groups.length?'':'fail'}`;
  const maxError=Math.max(...groups.map(g=>g.p95_error_px));
  $('translation-finding').textContent=`${pass} of ${groups.length} input groups meet the declared tolerance. The largest group’s 95th-percentile error is ${fmt(maxError,4)} px. All groups remain in the table; finite patches and scene texture affect the result.`;
  $('translation-rows').innerHTML=groups.map(g=>`<tr><th scope="row">${pair(g.injected_row_px,g.injected_col_px,3)}</th><td>${pair(g.median_measured_row_px,g.median_measured_col_px)}</td><td>${fmt(g.p95_error_px,5)}</td><td class="${g.meets_0_02px_p95_tolerance?'check-pass':'check-fail'}">${g.meets_0_02px_p95_tolerance?'Yes':'No'}</td></tr>`).join('');
  setLink('translation-download',base+report.translation.trials_csv);
}

function selectModel(id){
  const model=currentSurveys.models.find(m=>m.id===id)||currentSurveys.models[0];
  if(!model)return;
  $('model-select').value=model.id;surveyViewer.setModel(model);
  $('model-designation').textContent=model.designation;$('model-name').textContent=model.name;
  $('model-metrics').innerHTML=model.metrics.map(m=>`<div><dt>${esc(m.label)}</dt><dd>${esc(m.value)}</dd></div>`).join('');
  $('model-datum').textContent=model.datum;$('model-state').textContent=model.acquisition_state;
  $('model-assumptions').innerHTML=model.assumptions.map(a=>`<li>${esc(a)}</li>`).join('');
  $('model-source').textContent=model.source_title+' ↗';setLink('model-source',model.source_url);
  setLink('model-obj',`data/field/${current.id}/${model.id}.obj`);
  $('depth-guide').max=model.depth_m;$('depth-guide').value=0;$('guide-value').textContent='0.00 m';
  document.querySelectorAll('[data-model-view]').forEach(b=>b.setAttribute('aria-pressed',String(b.dataset.modelView==='orbit')));
  const url=new URL(location);url.searchParams.set('model',model.id);history.replaceState(null,'',url);
}

async function selectAcquisition(acquisition,token){
  const request=++acquisitionRequest;
  const reportPath='data/field/'+acquisition.report,base=reportPath.slice(0,reportPath.lastIndexOf('/')+1);
  const report=await load(reportPath);if(token!==generation||request!==acquisitionRequest)return;
  currentReport=report;
  $('acquisition-facts').innerHTML=[['Platform',acquisition.platform],['Collection',acquisition.date+' · '+fmt(report.timing.collection_duration_s,1)+' s'],['Native complex array',comma(report.native_shape[0])+' × '+comma(report.native_shape[1])],['Field depth results','None validated']].map(([k,v])=>`<div><dt>${esc(k)}</dt><dd>${esc(v)}</dd></div>`).join('');
  $('view-tabs').innerHTML=report.views.map((v,i)=>`<button type="button" data-image-view="${esc(v.id)}" aria-pressed="${i===0}">${esc(v.label)}</button>`).join('');
  const setView=async id=>{
    const view=report.views.find(v=>v.id===id)||report.views[0];
    document.querySelectorAll('[data-image-view]').forEach(b=>b.setAttribute('aria-pressed',String(b.dataset.imageView===view.id)));
    await radar.setView(view,base,report.survey_projection?.points);
  };
  document.querySelectorAll('[data-image-view]').forEach(b=>b.addEventListener('click',()=>setView(b.dataset.imageView).catch(showError)));
  $('landmark-controls').hidden=!report.survey_projection;
  $('show-landmarks').checked=false;radar.showPoints=false;radar.heightShift=0;$('height-shift').value=0;$('height-shift-value').textContent='0 m';
  renderChecks(report,base);
  setLink('report-download',reportPath);setLink('design-download','data/field/'+acquisition.design);setLink('audit-download','data/'+acquisition.geometry_audit);
  await setView(report.views[0].id);
}

async function selectSite(id){
  const token=++generation;current=sites.find(s=>s.id===id)||sites[0];
  $('load-status').textContent=`Loading ${current.name}…`;$('site-content').hidden=true;
  document.querySelectorAll('[data-site]').forEach(b=>b.setAttribute('aria-pressed',String(b.dataset.site===current.id)));
  const url=new URL(location);url.searchParams.set('site',current.id);if(!url.searchParams.has('model'))url.searchParams.set('model','idu');history.replaceState(null,'',url);
  const surveys=await load('data/field/'+current.surveys);if(token!==generation)return;
  currentSurveys=surveys;
  $('site-country').textContent=current.country+' / Research site';$('site-name').textContent=current.name;$('site-description').textContent=current.description;
  $('acquisition-select').innerHTML=current.acquisitions.map(a=>`<option value="${esc(a.id)}">${esc(a.date+' · '+a.platform+' · Complex SLC')}</option>`).join('');
  $('survey-workspace').hidden=!surveys.models.length;$('survey-empty').hidden=Boolean(surveys.models.length);
  if(surveys.models.length){
    $('model-select').innerHTML=surveys.models.map(m=>`<option value="${esc(m.id)}">${esc(m.name+' · '+m.designation.split(' · ')[0])}</option>`).join('');
    selectModel(new URL(location).searchParams.get('model'));
  }else{
    $('survey-empty').innerHTML=`<h3>A surveyed reference is still needed.</h3><p>${esc(surveys.status)}</p>`;
    const next=new URL(location);next.searchParams.delete('model');history.replaceState(null,'',next);
  }
  setLink('model-json','data/field/'+current.surveys);
  $('additional-surveys').innerHTML=surveys.additional_references.map(r=>`<article><h3>${esc(r.name)}</h3><p>${esc(r.description)} <a href="${esc(r.source_url)}">${esc(r.source_title)} ↗</a></p><p>${esc(r.qualification)}</p></article>`).join('');
  await selectAcquisition(current.acquisitions[0],token);if(token!==generation)return;
  $('load-status').textContent='';$('site-content').hidden=false;
  radar.draw();surveyViewer.draw();document.body.dataset.site=current.id;
}

function showError(error){$('load-status').textContent=`The saved atlas could not load. ${error.message}. You can still read the main research page and download the data.`;console.error(error);}
$('image-zoom-in').addEventListener('click',()=>radar.zoomBy(1.4));$('image-zoom-out').addEventListener('click',()=>radar.zoomBy(1/1.4));$('image-reset').addEventListener('click',()=>radar.reset());
$('show-landmarks').addEventListener('change',()=>{radar.showPoints=$('show-landmarks').checked;radar.draw();});
$('height-shift').addEventListener('input',()=>{radar.heightShift=Number($('height-shift').value);$('height-shift-value').textContent=`${radar.heightShift>0?'+':''}${radar.heightShift} m`;radar.draw();});
$('model-select').addEventListener('change',()=>selectModel($('model-select').value));
document.querySelectorAll('[data-model-view]').forEach(b=>b.addEventListener('click',()=>{document.querySelectorAll('[data-model-view]').forEach(s=>s.setAttribute('aria-pressed',String(s===b)));surveyViewer.view(b.dataset.modelView);}));
$('model-zoom-in').addEventListener('click',()=>surveyViewer.zoomBy(1.3));$('model-zoom-out').addEventListener('click',()=>surveyViewer.zoomBy(1/1.3));
$('depth-guide').addEventListener('input',()=>{surveyViewer.guide=Number($('depth-guide').value);$('guide-value').textContent=`${fmt(surveyViewer.guide,2)} m`;surveyViewer.draw();});
$('acquisition-select').addEventListener('change',()=>{const a=current.acquisitions.find(a=>a.id===$('acquisition-select').value);selectAcquisition(a,generation).catch(showError);});

try{
  ({sites}=await load('data/field/catalog.json'));
  $('site-tabs').innerHTML=sites.map(s=>`<button type="button" class="site-tab" data-site="${esc(s.id)}" aria-pressed="false"><span><strong>${esc(s.name)}</strong><small>${esc(s.country)} · ${s.acquisitions.length} acquisition${s.acquisitions.length===1?'':'s'}</small></span><span aria-hidden="true">↗</span></button>`).join('');
  document.querySelectorAll('[data-site]').forEach(b=>b.addEventListener('click',()=>selectSite(b.dataset.site).catch(showError)));
  await selectSite(new URL(location).searchParams.get('site'));
}catch(error){showError(error);}
