(function(root){
'use strict';const {math:M,quality}=Topologia;
class Renderer{
 constructor(canvas){this.canvas=canvas;this.ctx=canvas.getContext('2d');this.yaw=.65;this.pitch=.55;this.zoom=1;this.polygons=[];this.dpr=1;this.options={source:true,candidates:false,color:'growth'};}
 reset(){this.yaw=.65;this.pitch=.55;this.zoom=1;}
 resize(){const r=this.canvas.getBoundingClientRect(),d=Math.min(devicePixelRatio||1,2);this.w=r.width;this.h=r.height;this.dpr=d;if(this.canvas.width!==Math.round(r.width*d)||this.canvas.height!==Math.round(r.height*d)){this.canvas.width=Math.round(r.width*d);this.canvas.height=Math.round(r.height*d);}this.ctx.setTransform(d,0,0,d,0,0);}
 project(v){const c=Math.cos(this.yaw),s=Math.sin(this.yaw),cp=Math.cos(this.pitch),sp=Math.sin(this.pitch),x=c*v[0]+s*v[2],z=-s*v[0]+c*v[2],y=cp*v[1]-sp*z,depth=sp*v[1]+cp*z,factor=Math.min(this.w,this.h)*.36*this.zoom*4/(4-depth);return [this.w/2+x*factor,this.h/2-y*factor,depth];}
 path(p){const ctx=this.ctx;ctx.beginPath();ctx.moveTo(p[0][0],p[0][1]);for(let i=1;i<p.length;i++)ctx.lineTo(p[i][0],p[i][1]);ctx.closePath();}
 render(sim){this.resize();const stamp=[this.w,this.h,this.canvas.width,this.canvas.height,this.yaw,this.pitch,this.zoom].join(':');if(sim?.isPreview&&this.previewMesh===sim.source&&this.previewStamp===stamp)return;this.previewMesh=sim?.isPreview?sim.source:null;this.previewStamp=stamp;const c=this.ctx;c.clearRect(0,0,this.w,this.h);if(!sim)return;const pv=sim.vertices.map(v=>this.project(v)),sv=sim.source.vertices.map(v=>this.project(v));
  // Depth-buffer equivalent for opaque faces: painter order on compact mesh patches.
  const list=[];if(this.options.source||sim.isPreview)sim.source.faces.forEach((f,i)=>list.push({points:f.map(j=>sv[j]),source:true,depth:f.reduce((s,j)=>s+sv[j][2],0)/f.length,id:i}));
  sim.faces.forEach((f,i)=>{if(sim.solid[i]||this.options.candidates)list.push({points:f.map(j=>pv[j]),source:false,depth:f.reduce((s,j)=>s+pv[j][2],0)/4+.001,id:i});});list.sort((a,b)=>a.depth-b.depth);this.polygons=[];
  for(const item of list){const {points:p,id}=item;this.path(p);if(item.source){const shade=sim.isPreview?.32+.68*Math.abs(M.dot(M.normal(sim.source.vertices,sim.source.faces[id]),M.unit([-.4,.8,.6]))):0;c.fillStyle=sim.isPreview?`hsl(165 13% ${22+shade*34}%)`:'#10212beb';c.strokeStyle=sim.isPreview?'#92aea72b':'#47607270';c.lineWidth=.6;c.fill();if(sim.isPreview||!this.options.candidates)c.stroke();continue;}
   if(!sim.solid[id]){c.fillStyle='#173b43';c.strokeStyle='#81d8df';c.lineWidth=.9;c.fill();c.stroke();continue;}const f=sim.faces[id],n=M.normal(sim.vertices,f),light=.5+.5*Math.abs(M.dot(n,M.unit([-.5,1,.7])));let hue=158,sat=55,lum=28+light*24;
   if(this.options.color==='phase')hue=(sim.phase[id]/(2*Math.PI)*360+360)%360;else if(this.options.color==='quality')hue=quality(sim.vertices,f)*155;else if(this.options.color==='nuclei')hue=(sim.owner[id]*137.508)%360;else if(this.options.color==='error')hue=Math.max(0,155-200*sim.errorField[id]);else if(this.options.color==='vertex')hue=f.reduce((s,i)=>s+sim.vertexWeight[i],0)/4*155;else if(this.options.color==='area')hue=Math.max(0,155-200*Math.abs(Math.log(Topologia.faceArea(sim.vertices,f)/sim.targetArea[id])));else if(this.options.color==='mass')hue=210-Math.min(1,sim.mass[id])*160;else if(sim.tick-sim.birth[id]<45){hue=38;sat=75;lum=45+light*20;}
   c.fillStyle=`hsl(${hue} ${sat}% ${lum}%)`;c.strokeStyle=`hsl(${hue} ${sat*.7}% ${Math.min(85,lum+18)}%)`;c.lineWidth=.65;c.fill();c.stroke();if(this.options.color==='vertex'){for(let k=0;k<4;k++){const a=f[k],b=f[(k+1)%4],e=sim.relationLookup.get(a<b?`${a}:${b}`:`${b}:${a}`);c.beginPath();c.moveTo(p[k][0],p[k][1]);c.lineTo(p[(k+1)%4][0],p[(k+1)%4][1]);c.strokeStyle=`hsl(${e.weight*155} 65% 65%)`;c.lineWidth=.5+1.8*e.weight;c.stroke();}}this.polygons.push(item);
  }
  // Draw candidate seeds even when hidden by the source surface? Keep depth ordering;
  // only visible active polygons are selectable, and source/candidate faces can seed.
  c.fillStyle='#7b98aa';c.font='10px Consolas';c.fillText('Y',22,this.h-58);c.strokeStyle='#486476';c.lineWidth=1;c.beginPath();c.moveTo(25,this.h-65);c.lineTo(25,this.h-88);c.moveTo(25,this.h-65);c.lineTo(47,this.h-65);c.stroke();c.fillText('X',51,this.h-62);
  this.all=sim.faces.map((f,i)=>({id:i,points:f.map(j=>pv[j]),depth:f.reduce((s,j)=>s+pv[j][2],0)/4})).sort((a,b)=>b.depth-a.depth);
 }
 pick(x,y){for(const face of this.all||[]){const p=face.points;let inside=false;for(let i=0,j=p.length-1;i<p.length;j=i++){const a=p[i],b=p[j];if((a[1]>y)!==(b[1]>y)&&x<(b[0]-a[0])*(y-a[1])/(b[1]-a[1])+a[0])inside=!inside;}if(inside)return face.id;}return -1;}
}
root.MeshRenderer=Renderer;
})(globalThis);
