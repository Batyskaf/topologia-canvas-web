(function(root){
'use strict';
const T=typeof module!=='undefined'&&module.exports?require('./engine.js'):root.Topologia,{mean}=T.math;
const key=(a,b)=>a<b?`${a}:${b}`:`${b}:${a}`;
function convert(mesh){
 if(mesh.faces.some(f=>f.length!==3&&f.length!==4))throw Error('Konwersja oczekuje trójkątów i quadów.');
 const topo=T.topology(mesh),split=new Set(),odd=new Set(),vertices=mesh.vertices.map(v=>v.slice()),vertexSources=(mesh.vertexSources||mesh.vertices.map(()=>[])).map(a=>a.slice()),faces=[],faceSources=[];
 for(const f of mesh.faces)if(f.length===3)for(let k=0;k<3;k++)split.add(key(f[k],f[(k+1)%3]));
 const edgeIds=mesh.faces.map(f=>f.map((v,k)=>key(v,f[(k+1)%f.length])));
 mesh.faces.forEach((f,i)=>{if(f.length===4&&edgeIds[i].filter(e=>split.has(e)).length%2)odd.add(i);});
 // Pair odd quad boundaries along unsplit shared edges. Each path adds one
 // split at its endpoints and two at intermediate faces: no T-junctions.
 while(odd.size){
  const start=odd.values().next().value,queue=[start],parent=new Map([[start,null]]);let end=null,boundary=null;
  for(let head=0;head<queue.length&&end===null;head++){const id=queue[head];for(const e of edgeIds[id]){
   if(split.has(e))continue;const row=topo.edges.get(e),other=row.faces.find(i=>i!==id);
   if(other===undefined){end=id;boundary=e;break;}
   if(parent.has(other))continue;parent.set(other,{face:id,edge:e});if(odd.has(other)){end=other;break;}queue.push(other);
  }}
  if(end===null)throw Error('Nie można uzgodnić parzystości granic quadów.');
  if(boundary)split.add(boundary);else odd.delete(end);
  let cursor=end;while(cursor!==start){const p=parent.get(cursor);split.add(p.edge);cursor=p.face;}odd.delete(start);
 }
 const mids=new Map();for(const e of split){const row=topo.edges.get(e),id=vertices.length;vertices.push(mean([vertices[row.a],vertices[row.b]]));vertexSources.push([...new Set(row.faces.flatMap(i=>mesh.faceSources?.[i]||[]))]);mids.set(e,id);}
 const emit=(f,sources)=>{faces.push(f);faceSources.push(sources.slice());};
 mesh.faces.forEach((f,id)=>{
  const sources=mesh.faceSources?.[id]||[],count=edgeIds[id].filter(e=>split.has(e)).length;
  if(f.length===4&&count===0){emit(f.slice(),sources);return;}
  const boundary=[];f.forEach((v,k)=>{boundary.push(v);const mid=mids.get(edgeIds[id][k]);if(mid!==undefined)boundary.push(mid);});
  if(f.length===4&&count===2){
   // Opposite splits allow two quads; adjacent splits need a three-quad fan.
   let best=null,score=0;for(let i=0;i<3;i++){const a=[0,1,2,3].map(k=>boundary[(i+k)%6]),b=[3,4,5,0].map(k=>boundary[(i+k)%6]),q=Math.min(T.quality(vertices,a),T.quality(vertices,b));if(q>score){score=q;best=[a,b];}}
   if(best&&score>1e-8){for(const quad of best)emit(quad,sources);return;}
  }
  const center=vertices.length;vertices.push(mean(f.map(i=>vertices[i])));vertexSources.push(sources.slice());let best=null,bestQuality=-1;
  for(const offset of [0,1]){const fan=[];for(let k=0;k<boundary.length;k+=2)fan.push([boundary[(offset+k)%boundary.length],boundary[(offset+k+1)%boundary.length],boundary[(offset+k+2)%boundary.length],center]);const q=Math.min(...fan.map(face=>T.quality(vertices,face)));if(q>bestQuality){bestQuality=q;best=fan;}}
  if(bestQuality<=1e-10)throw Error('Nie można utworzyć poprawnych quadów w lokalnym fragmencie.');for(const quad of best)emit(quad,sources);
 });
 return {vertices,faces,vertexSources,faceSources};
}
const api={convert};if(typeof module!=='undefined'&&module.exports)module.exports=api;else root.TopologiaQuads=api;
})(globalThis);
