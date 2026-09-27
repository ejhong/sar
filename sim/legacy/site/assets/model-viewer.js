// Orthographic metric survey viewer. Geometry comes only from the survey JSON.
// Canvas keeps it portable and usable without WebGL, CDNs or a graphics driver.
export class SurveyViewer {
  constructor(canvas) {
    this.canvas = canvas; this.ctx = canvas.getContext('2d');
    this.model = null; this.guide = 0; this.mode = 'orbit'; this.zoom = 1;
    this.yaw = -.63; this.pitch = .43;
    this.observer = new ResizeObserver(() => this.draw()); this.observer.observe(canvas);
    let drag = null;
    canvas.addEventListener('pointerdown', e => {
      drag = [e.clientX,e.clientY]; canvas.setPointerCapture(e.pointerId); canvas.focus({preventScroll:true});
    });
    canvas.addEventListener('pointermove', e => {
      if(!drag) return;
      this.orbitInteraction();
      this.yaw += (e.clientX-drag[0])*.009;
      this.pitch = Math.max(-.35,Math.min(1.56,this.pitch+(e.clientY-drag[1])*.007));
      drag = [e.clientX,e.clientY]; this.draw();
    });
    for(const event of ['pointerup','pointercancel','lostpointercapture']) canvas.addEventListener(event,()=>drag=null);
    canvas.addEventListener('wheel',e=>{e.preventDefault();this.zoomBy(Math.exp(-Math.sign(e.deltaY)*.11));},{passive:false});
    canvas.addEventListener('keydown',e=>{
      if(!['ArrowLeft','ArrowRight','ArrowUp','ArrowDown','+','=','-','Home'].includes(e.key))return;
      e.preventDefault();
      if(e.key.startsWith('Arrow'))this.orbitInteraction();
      if(e.key==='ArrowLeft')this.yaw-=.13;
      if(e.key==='ArrowRight')this.yaw+=.13;
      if(e.key==='ArrowUp')this.pitch=Math.min(1.56,this.pitch+.1);
      if(e.key==='ArrowDown')this.pitch=Math.max(-.35,this.pitch-.1);
      if(e.key==='+'||e.key==='=')this.zoomBy(1.2);
      if(e.key==='-')this.zoomBy(1/1.2);
      if(e.key==='Home')this.view('orbit');
      this.draw();
    });
  }
  setModel(model) { this.model=model; this.guide=0; this.view('orbit'); }
  orbitInteraction() {
    if(this.mode==='plan'||this.mode==='section'){
      this.mode='orbit';this.canvas.dispatchEvent(new CustomEvent('viewchange',{detail:this.mode}));
    }
  }
  zoomBy(factor) { this.zoom=Math.max(.5,Math.min(6,this.zoom*factor));this.draw(); }
  view(mode) {
    this.mode=mode;this.zoom=1;
    this.yaw=mode==='section'?Math.PI/2:mode==='plan'?0:-.63;
    this.pitch=mode==='plan'?Math.PI/2:mode==='section'?0:.43;
    this.canvas.dispatchEvent(new CustomEvent('viewchange',{detail:mode}));
    this.draw();
  }
  draw() {
    if(!this.model)return;
    const canvas=this.canvas,ctx=this.ctx,w=canvas.clientWidth,h=canvas.clientHeight;
    if(!w||!h)return;
    const ratio=Math.min(devicePixelRatio||1,2);
    if(canvas.width!==Math.round(w*ratio)||canvas.height!==Math.round(h*ratio)){
      canvas.width=Math.round(w*ratio);canvas.height=Math.round(h*ratio);
    }
    ctx.setTransform(ratio,0,0,ratio,0,0);ctx.clearRect(0,0,w,h);
    const model=this.model;
    const center=this.mode==='chamber'?model.chamber_focus:[0,-1,-model.depth_m/2];
    const transform=p=>{
      const x=p[0]-center[0],y=p[1]-center[1],z=p[2]-center[2];
      const a=x*Math.cos(this.yaw)-y*Math.sin(this.yaw),b=x*Math.sin(this.yaw)+y*Math.cos(this.yaw);
      return [a,b*Math.sin(this.pitch)-z*Math.cos(this.pitch),b*Math.cos(this.pitch)+z*Math.sin(this.pitch)];
    };
    const all=model.parts.flatMap(p=>p.polygon.flatMap(([x,y])=>[[x,y,p.top_m],[x,y,p.bottom_m]]));
    const fit=this.mode==='chamber'?model.parts.filter(p=>p.kind==='chamber').flatMap(p=>p.polygon.flatMap(([x,y])=>[[x,y,p.top_m],[x,y,p.bottom_m]])):all;
    const mapped=fit.map(transform),xs=mapped.map(p=>p[0]),ys=mapped.map(p=>p[1]);
    const scale=Math.min((w-110)/(Math.max(...xs)-Math.min(...xs)+3),(h-100)/(Math.max(...ys)-Math.min(...ys)+1.2))*this.zoom;
    const bx=(Math.max(...xs)+Math.min(...xs))/2,by=(Math.max(...ys)+Math.min(...ys))/2;
    const project=p=>{const t=transform(p);return [w/2+(t[0]-bx)*scale,h/2+(t[1]-by)*scale,t[2]];};
    const path=vertices=>{ctx.beginPath();vertices.forEach((v,i)=>{const p=project(v);i?ctx.lineTo(p[0],p[1]):ctx.moveTo(p[0],p[1]);});};
    const line=(a,b,color,width=1)=>{path([a,b]);ctx.strokeStyle=color;ctx.lineWidth=width;ctx.stroke();};
    const text=(word,p,color='#665842',align='left')=>{
      const [x,y]=project(p);ctx.font='10px "IBM Plex Mono", monospace';ctx.textAlign=align;ctx.fillStyle=color;
      ctx.fillText(word,x,y);
    };
    const extent=Math.max(4,...all.map(p=>Math.max(Math.abs(p[0]),Math.abs(p[1]))));
    // Ground is a reference plane, not a terrain or masonry reconstruction.
    for(let k=-Math.ceil(extent);k<=Math.ceil(extent);k++){
      line([k,-extent,0],[k,extent,0],'#7b958125',.7);
      line([-extent,k,0],[extent,k,0],'#7b958125',.7);
    }
    const faces=[];
    for(const part of model.parts){
      const top=part.polygon.map(([x,y])=>[x,y,part.top_m]);
      const bottom=part.polygon.map(([x,y])=>[x,y,part.bottom_m]);
      const polygons=[top,bottom,...top.map((p,i)=>[p,top[(i+1)%top.length],bottom[(i+1)%top.length],bottom[i]])];
      polygons.forEach((polygon,i)=>faces.push({polygon,depth:polygon.reduce((s,p)=>s+transform(p)[2],0)/polygon.length,kind:part.kind,side:i}));
    }
    faces.sort((a,b)=>a.depth-b.depth);
    for(const f of faces){
      path(f.polygon);ctx.closePath();
      ctx.fillStyle=f.kind==='chamber'?(f.side===1?'#b78e5850':'#c9a27038'):'#67947e23';ctx.fill();
      ctx.strokeStyle=f.kind==='chamber'?'#916632c0':'#537966ad';ctx.lineWidth=f.kind==='chamber'?1.25:.95;ctx.stroke();
    }
    // Depth guide is a ruler plane, never a slice of measured radar values.
    if(this.guide>0){
      const z=-this.guide,e=Math.min(extent,4);
      path([[-e,-e,z],[e,-e,z],[e,e,z],[-e,e,z]]);ctx.closePath();
      ctx.fillStyle='#1e756014';ctx.fill();ctx.strokeStyle='#28715c77';ctx.setLineDash([3,4]);ctx.stroke();ctx.setLineDash([]);
      const gx=Math.max(64,Math.min(...[[-e,-e,z],[e,-e,z],[e,e,z],[-e,e,z]].map(p=>project(p)[0]))-9);
      ctx.font='10px "IBM Plex Mono",monospace';ctx.textAlign='right';ctx.fillStyle='#24665b';
      ctx.fillText(`${this.guide.toFixed(2)} m`,gx,project([0,0,z])[1]);
    }
    if(this.mode!=='plan'&&this.mode!=='chamber'){
      const rulerX=Math.min(w-48,Math.max(...all.map(p=>project(p)[0]))+32);
      const shaft=model.parts.find(p=>p.kind==='shaft'),n=shaft.polygon.length;
      const cx=shaft.polygon.reduce((s,p)=>s+p[0],0)/n,cy=shaft.polygon.reduce((s,p)=>s+p[1],0)/n;
      const py=depth=>project([cx,cy,-depth])[1];
      ctx.beginPath();ctx.moveTo(rulerX,py(0));ctx.lineTo(rulerX,py(model.depth_m));ctx.strokeStyle='#8b9a8970';ctx.lineWidth=1;ctx.stroke();
      const step=model.depth_m>15?5:2;
      for(let depth=0;depth<=model.depth_m;depth+=step){
        ctx.beginPath();ctx.moveTo(rulerX-3,py(depth));ctx.lineTo(rulerX+3,py(depth));ctx.strokeStyle='#7b8f7f';ctx.stroke();
        ctx.font='10px "IBM Plex Mono",monospace';ctx.textAlign='left';ctx.fillStyle='#5d6f61';ctx.fillText(`${depth} m`,rulerX+7,py(depth)+3);
      }
    }
    if(Math.abs(this.pitch)>.1){
      line([-extent,-extent,0],[-extent+1.4,-extent,0],'#587866',1.5);
      line([-extent,-extent,0],[-extent,-extent+1.4,0],'#587866',1.5);
      text('E',[-extent+1.65,-extent,0],'#416650');
      text('N',[-extent,-extent+1.65,0],'#416650');
    }else if(this.mode==='section'){
      ctx.font='10px "IBM Plex Mono", monospace';ctx.textAlign='left';ctx.fillStyle='#587866';
      ctx.fillText('N ←  Section  → S',18,49);
    }
    ctx.font='10px "IBM Plex Mono", monospace';ctx.textAlign='right';ctx.fillStyle='#5e7164';
    const barMetres=this.mode==='chamber'?1:model.depth_m>15?5:2;
    const barWidth=barMetres*scale;
    if(barWidth<w*.35){ctx.beginPath();ctx.moveTo(w-22-barWidth,h-29);ctx.lineTo(w-22,h-29);ctx.strokeStyle='#708474';ctx.lineWidth=2;ctx.stroke();ctx.fillText(`${barMetres} m`,w-22,h-37);}
    canvas.dataset.model=model.id;canvas.dataset.view=this.mode;canvas.dataset.zoom=this.zoom.toFixed(3);
  }
}
