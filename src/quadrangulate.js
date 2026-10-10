(function(root){
'use strict';
const T=typeof module!=='undefined'&&module.exports?require('./engine.js'):root.Topologia,{mean}=T.math;
const key=(a,b)=>a<b?`${a}:${b}`:`${b}:${a}`;
function convert(mesh,minPairQuality=.2){
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
   if(best&&score>minPairQuality){for(const quad of best)emit(quad,sources);return;}
  }
  const center=vertices.length;vertices.push(mean(f.map(i=>vertices[i])));vertexSources.push(sources.slice());let best=null,bestQuality=-1;
  for(const offset of [0,1]){const fan=[];for(let k=0;k<boundary.length;k+=2)fan.push([boundary[(offset+k)%boundary.length],boundary[(offset+k+1)%boundary.length],boundary[(offset+k+2)%boundary.length],center]);const q=Math.min(...fan.map(face=>T.quality(vertices,face)));if(q>bestQuality){bestQuality=q;best=fan;}}
  if(bestQuality<=1e-10)throw Error('Nie można utworzyć poprawnych quadów w lokalnym fragmencie.');for(const quad of best)emit(quad,sources);
 });
 return {vertices,faces,vertexSources,faceSources};
}
function regularize(mesh,source,surface,sourceIds){
 const M=T.math,A=typeof module!=='undefined'&&module.exports?require('./audit.js'):root.TopologiaAudit;
 const {add,sub,mul,dot,cross,len,unit,mean,normal}=M;
 if(mesh.faces.some(f=>f.length!==4||new Set(f).size!==4))throw Error('Wynik musi zawierać wyłącznie quady z czterema różnymi wierzchołkami.');
 const v=mesh.vertices,srcTopo=T.topology(source),srcNormals=source.faces.map(f=>normal(source.vertices,f)),protectedSource=new Set();
 for(const e of srcTopo.edges.values())if(e.faces.length===1||dot(srcNormals[e.faces[0]],srcNormals[e.faces[1]])<Math.cos(Math.PI*50/180)){protectedSource.add(e.a);protectedSource.add(e.b);}
 let rotated=0,moved=0;const qualityBefore=mesh.faces.reduce((s,f)=>s+T.quality(v,f),0)/mesh.faces.length;
 // Rotating the shared edge of two planar quads keeps the same six-sided patch.
 // Convexity and coplanarity make the replacement cover exactly the old area.
 for(let pass=0;pass<2;pass++){
  const topo=T.topology(mesh),used=new Set();
  for(const e of topo.edges.values()){
   if(e.faces.length!==2||e.faces.some(i=>used.has(i)))continue;
   const [a,b]=e.faces,fa=mesh.faces[a],fb=mesh.faces[b],n=normal(v,fa);
   if(dot(n,normal(v,fb))<.999999)continue;
   const edges=[];for(const f of [fa,fb])for(let k=0;k<4;k++)if(key(f[k],f[(k+1)%4])!==key(e.a,e.b))edges.push([f[k],f[(k+1)%4]]);
   const next=new Map(edges),ring=[edges[0][0]];for(let k=1;k<6;k++)ring.push(next.get(ring.at(-1)));
   if(new Set(ring).size!==6||ring.some(i=>i===undefined)||next.get(ring[5])!==ring[0])continue;
   const scale=Math.max(...ring.map(i=>len(sub(v[i],v[ring[0]]))));
   if(ring.some(i=>Math.abs(dot(sub(v[i],v[ring[0]]),n))>scale*1e-7)||ring.some((i,k)=>dot(cross(sub(v[ring[(k+1)%6]],v[i]),sub(v[ring[(k+2)%6]],v[ring[(k+1)%6]])),n)<-scale*scale*1e-10))continue;
   let best=null,score=Math.min(T.quality(v,fa),T.quality(v,fb));
   for(let k=0;k<3;k++){
    const f=[0,1,2,3].map(j=>ring[(k+j)%6]),g=[3,4,5,0].map(j=>ring[(k+j)%6]),diagonal=key(f[0],f[3]);
    if(topo.edges.has(diagonal))continue;
    const q=Math.min(T.quality(v,f),T.quality(v,g));if(q>score+1e-5){score=q;best=[f,g];}
   }
   if(best){mesh.faces[a]=best[0];mesh.faces[b]=best[1];const refs=[...new Set([...(mesh.faceSources?.[a]||[]),...(mesh.faceSources?.[b]||[])])];if(mesh.faceSources){mesh.faceSources[a]=refs.slice();mesh.faceSources[b]=refs.slice();}for(const i of ring)mesh.vertexSources[i]=[...new Set([...(mesh.vertexSources[i]||[]),...refs])];used.add(a);used.add(b);rotated++;}
  }
 }
 const topo=T.topology(mesh),pinned=new Set();for(const e of topo.edges.values())if(e.faces.length===1){pinned.add(e.a);pinned.add(e.b);}
 const neighbors=v.map(()=>[]);for(const e of topo.edges.values()){neighbors[e.a].push(e.b);neighbors[e.b].push(e.a);}
 const charts=v.map((_,i)=>{const ids=new Set(mesh.vertexSources[i]||[]);for(const id of [...ids])for(const j of srcTopo.adj[id]||[])ids.add(j);return ids;});
 for(let i=0;i<v.length;i++){if(protectedSource.has(sourceIds[i]))pinned.add(i);const ns=[...charts[i]].map(j=>srcNormals[j]);if(ns.some(n=>ns.some(m=>dot(n,m)<Math.cos(Math.PI*50/180))))pinned.add(i);}
 const collider=new A.CollisionIndex(mesh);
 for(let pass=0;pass<3;pass++)for(let i=0;i<v.length;i++){
  if(pinned.has(i)||!neighbors[i].length)continue;
  const fs=topo.incident[i],before=v[i],qs=fs.map(j=>T.quality(v,mesh.faces[j])),floor=Math.min(...qs),oldCost=qs.reduce((s,q)=>s+(1-q)**2,0),ns=fs.map(j=>normal(v,mesh.faces[j])),direction=sub(mean(neighbors[i].map(j=>v[j])),before);
  for(const rate of [.5,.25,.125]){
   const hit=surface.sample(add(before,mul(direction,rate)),j=>charts[i].has(j));if(hit.face<0)continue;
   v[i]=hit.point;const next=fs.map(j=>T.quality(v,mesh.faces[j])),cost=next.reduce((s,q)=>s+(1-q)**2,0);
   if(Math.min(...next)>=floor-1e-9&&cost<oldCost-1e-9&&fs.every((j,k)=>dot(normal(v,mesh.faces[j]),ns[k])>.5)&&collider.canMove(fs)){moved++;break;}
   v[i]=before;collider.refreshFaces(fs);
  }
 }
 return {rotatedEdges:rotated,movedVertices:moved,qualityBefore,qualityAfter:mesh.faces.reduce((s,f)=>s+T.quality(v,f),0)/mesh.faces.length};
}
const api={convert,regularize};if(typeof module!=='undefined'&&module.exports)module.exports=api;else root.TopologiaQuads=api;
})(globalThis);
