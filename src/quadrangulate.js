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
// Replace three-quad star fans by two quads across their six-vertex boundary.
function mergeFans(mesh,source,surface,tolerance){
 const {normal,dot,mean}=T.math,A=typeof module!=='undefined'&&module.exports?require('./audit.js'):root.TopologiaAudit,backup=JSON.parse(JSON.stringify(mesh)),topo=T.topology(mesh),used=new Set(),removed=new Set();
 for(let center=0;center<mesh.vertices.length;center++){
  const region=topo.incident[center];if(region.length!==3||region.some(i=>used.has(i)))continue;
  const edges=[];for(const i of region){const f=mesh.faces[i];for(let k=0;k<4;k++){const e=topo.edges.get(key(f[k],f[(k+1)%4]));if(e.faces.filter(j=>region.includes(j)).length===1)edges.push([f[k],f[(k+1)%4]]);}}
  if(edges.length!==6)continue;const next=new Map(edges),ring=[edges[0][0]];for(let k=1;k<6;k++)ring.push(next.get(ring.at(-1)));if(new Set(ring).size!==6||ring.some(i=>i===undefined)||next.get(ring[5])!==ring[0])continue;
  const v=mesh.vertices,n=normal(v,mesh.faces[region[0]]),refs=[...new Set(region.flatMap(i=>mesh.faceSources[i]||[]))],old=region.map(i=>T.quality(v,mesh.faces[i])),limit=Math.min(.025*Math.sqrt(region.reduce((s,i)=>s+T.faceArea(v,mesh.faces[i]),0)),2*tolerance);
  if(region.some(i=>dot(n,normal(v,mesh.faces[i]))<.995))continue;
  let best=null,score=Math.max(.35,Math.min(...old));for(let k=0;k<3;k++){const pair=[[0,1,2,3],[3,4,5,0]].map(indices=>indices.map(j=>ring[(j+k)%6])),q=Math.min(...pair.map(f=>T.quality(v,f)));if(q<=score||pair.some(f=>dot(normal(v,f),n)<.995))continue;let hit=false;for(const f of pair)for(const p of [mean(f.map(i=>v[i])),...T.triangulate(v,f).map(t=>mean(t.map(i=>v[i])))])if(surface.sample(p,i=>refs.includes(i)).distance>limit)hit=true;if(!hit){score=q;best=pair;}}
  if(!best)continue;mesh.faces[region[0]]=best[0];mesh.faces[region[1]]=best[1];mesh.faceSources[region[0]]=refs.slice();mesh.faceSources[region[1]]=refs.slice();for(const i of ring)mesh.vertexSources[i]=[...new Set([...(mesh.vertexSources[i]||[]),...refs])];region.forEach(i=>used.add(i));removed.add(region[2]);
 }
 if(!removed.size)return 0;
 mesh.faces=mesh.faces.filter((_,i)=>!removed.has(i));mesh.faceSources=mesh.faceSources.filter((_,i)=>!removed.has(i));
 // Keep unused vertices until the shared compaction stage; audit a compact copy.
 const ids=[...new Set(mesh.faces.flat())],map=new Map(ids.map((i,j)=>[i,j])),check={vertices:ids.map(i=>mesh.vertices[i]),faces:mesh.faces.map(f=>f.map(i=>map.get(i)))},report=A.inspect(check);
 if(!report.valid||report.chi!==topo.chi||report.components!==topo.components){Object.assign(mesh,backup);return 0;}return removed.size;
}
// Merge four quads into a curvature-controlled patch. Retain shared edge
// midpoints next to finer patches; reconcile parity before quadrangulating.
function coarsen(mesh,source,surface,sourceIds,tolerance=.04){
 const A=typeof module!=='undefined'&&module.exports?require('./audit.js'):root.TopologiaAudit,{dot,normal,mean}=T.math;
 const original=JSON.parse(JSON.stringify(mesh)),oldIds=sourceIds.slice();
 const fans=mergeFans(mesh,source,surface,tolerance);if(fans){const ids=[...new Set(mesh.faces.flat())].sort((a,b)=>a-b),map=new Map(ids.map((v,i)=>[v,i]));mesh.faces=mesh.faces.map(f=>f.map(i=>map.get(i)));mesh.vertices=ids.map(i=>mesh.vertices[i]);mesh.vertexSources=ids.map(i=>mesh.vertexSources[i]);const mapped=ids.map(i=>sourceIds[i]);sourceIds.length=mapped.length;mapped.forEach((v,i)=>sourceIds[i]=v);}
 const topo=T.topology(mesh),v=mesh.vertices,occupied=new Set(),patches=[],owner=new Map();
 for(let center=0;center<v.length;center++){
  const region=topo.incident[center];if(region.length!==4||region.some(i=>occupied.has(i)))continue;
  const edges=[];for(const i of region)for(let k=0;k<4;k++){const f=mesh.faces[i],e=topo.edges.get(key(f[k],f[(k+1)%4]));if(e.faces.filter(j=>region.includes(j)).length===1)edges.push([f[k],f[(k+1)%4]]);}
  if(edges.length!==8)continue;const next=new Map(edges),ring=[edges[0][0]];for(let k=1;k<8;k++)ring.push(next.get(ring.at(-1)));if(new Set(ring).size!==8||ring.some(i=>i===undefined)||next.get(ring[7])!==ring[0])continue;
  const radial=new Set(region.flatMap(i=>{const f=mesh.faces[i],k=f.indexOf(center);return [f[(k+1)%4],f[(k+3)%4]];}));if(radial.size!==4)continue;
  const corners=ring.filter(i=>!radial.has(i));if(corners.length!==4||T.quality(v,corners)<.5)continue;
  const weights=region.map(i=>T.faceArea(v,mesh.faces[i])),average=T.math.unit(region.reduce((sum,i,k)=>T.math.add(sum,T.math.mul(normal(v,mesh.faces[i]),weights[k])),[0,0,0]));
  const n=normal(v,corners),refs=[...new Set(region.flatMap(i=>mesh.faceSources[i]||[]))],area=region.reduce((s,i)=>s+T.faceArea(v,mesh.faces[i]),0);
  // Compare the merged boundary with the area-weighted mean orientation;
  // reject patches whose constituent normals vary by more than 12 degrees.
  if(dot(n,average)<.995||region.some(i=>dot(normal(v,mesh.faces[i]),average)<Math.cos(Math.PI/15))||refs.some(i=>dot(normal(source.vertices,source.faces[i]),n)<Math.cos(Math.PI/15)))continue;
  const limit=Math.min(2*tolerance,Math.sqrt(area)*.025),flat=new T.Surface({vertices:corners.map(i=>v[i]),faces:[[0,1,2,3]]});
  if([center,...ring].some(i=>flat.sample(v[i]).distance>limit)||[mean(corners.map(i=>v[i])),...T.triangulate(v,corners).map(f=>mean(f.map(i=>v[i])))].some(p=>surface.sample(p,i=>refs.includes(i)).distance>limit))continue;
  const mids=corners.map((a,k)=>{const pos=ring.indexOf(a);return ring[(pos+1)%8];});
  const id=patches.length;patches.push({corners,mids,refs,region});for(const i of region){occupied.add(i);owner.set(i,id);}
 }
 if(!patches.length)return {collapsedQuads:fans,before:original.faces.length,after:mesh.faces.length,mergedFans:fans};
 const splits=patches.map(p=>p.mids.map((mid,k)=>{const a=p.corners[k],b=p.corners[(k+1)%4],outside=[...new Set([...(topo.edges.get(key(a,mid))?.faces||[]),...(topo.edges.get(key(mid,b))?.faces||[])])].filter(i=>!p.region.includes(i));return outside.length===0||outside.some(i=>!owner.has(i))||new Set(outside.map(i=>owner.get(i))).size!==1;}));
 const adj=new Map();patches.forEach((p,i)=>p.corners.forEach((a,k)=>{const e=key(a,p.corners[(k+1)%4]);if(!adj.has(e))adj.set(e,[]);adj.get(e).push([i,k]);}));
 // Odd boundaries are joined through coarse neighbors, adding a shared split.
 const odd=new Set(splits.map((s,i)=>s.filter(Boolean).length%2?i:-1).filter(i=>i>=0));
 while(odd.size){const start=odd.values().next().value,queue=[start],parent=new Map([[start,null]]);let end=null;
  for(let h=0;h<queue.length&&end===null;h++){const i=queue[h];for(let k=0;k<4;k++){if(splits[i][k])continue;const rows=adj.get(key(patches[i].corners[k],patches[i].corners[(k+1)%4])),other=rows.find(r=>r[0]!==i);if(!other||parent.has(other[0]))continue;parent.set(other[0],{i,k,otherK:other[1]});if(odd.has(other[0])){end=other[0];break;}queue.push(other[0]);}}
  if(end===null)return {collapsedQuads:fans,before:original.faces.length,after:mesh.faces.length,mergedFans:fans,rejected:true};
  odd.delete(start);odd.delete(end);let i=end;while(i!==start){const p=parent.get(i);splits[i][p.otherK]=true;splits[p.i][p.k]=true;i=p.i;}
 }
 const faces=[],faceSources=[];const emit=(f,refs)=>{faces.push(f);faceSources.push(refs.slice());};mesh.faces.forEach((f,i)=>{if(!occupied.has(i))emit(f.slice(),mesh.faceSources[i]||[]);});
 patches.forEach((p,i)=>{const ring=[];p.corners.forEach((a,k)=>{ring.push(a);if(splits[i][k])ring.push(p.mids[k]);});
  if(ring.length===4){emit(ring,p.refs);return;}
  if(ring.length===6){let best=null,q=.2;for(let k=0;k<3;k++){const a=[0,1,2,3].map(j=>ring[(k+j)%6]),b=[3,4,5,0].map(j=>ring[(k+j)%6]),score=Math.min(T.quality(v,a),T.quality(v,b));if(score>q){q=score;best=[a,b];}}if(best){best.forEach(f=>emit(f,p.refs));return;}}
  const center=v.length;v.push(mean(p.corners.map(j=>v[j])));mesh.vertexSources.push(p.refs.slice());sourceIds.push(-1);let best=null,score=-1;for(const offset of [0,1]){const fan=[];for(let k=0;k<ring.length;k+=2)fan.push([ring[(k+offset)%ring.length],ring[(k+offset+1)%ring.length],ring[(k+offset+2)%ring.length],center]);const q=Math.min(...fan.map(f=>T.quality(v,f)));if(q>score){score=q;best=fan;}}best.forEach(f=>emit(f,p.refs));
 });
 faces.forEach((f,i)=>f.forEach(id=>{mesh.vertexSources[id]=[...new Set([...(mesh.vertexSources[id]||[]),...faceSources[i]])];}));
 const ids=[...new Set(faces.flat())].sort((a,b)=>a-b),map=new Map(ids.map((v,i)=>[v,i]));mesh.faces=faces.map(f=>f.map(i=>map.get(i)));mesh.faceSources=faceSources;mesh.vertices=ids.map(i=>v[i]);mesh.vertexSources=ids.map(i=>mesh.vertexSources[i]);const mapped=ids.map(i=>sourceIds[i]);sourceIds.length=mapped.length;mapped.forEach((id,i)=>sourceIds[i]=id);
 const qualityFloor=Math.min(.3,...original.faces.map(f=>T.quality(original.vertices,f)));
 const report=A.inspect(mesh);if(mesh.faces.some(f=>T.quality(mesh.vertices,f)<qualityFloor-1e-8)||!report.valid||report.chi!==topo.chi||report.components!==topo.components){Object.assign(mesh,original);sourceIds.length=oldIds.length;oldIds.forEach((id,i)=>sourceIds[i]=id);return {collapsedQuads:0,before:original.faces.length,after:original.faces.length,rejected:true,reason:report.errors.join(' ')||'Scalenie pogorszyłoby jakość quadów.'};}
 return {collapsedQuads:original.faces.length-mesh.faces.length,before:original.faces.length,after:mesh.faces.length,mergedPatches:patches.length,mergedFans:fans};
}
// Exchange a shared edge during growth, including gently curved patches.
// Face count and face IDs stay stable so diffusion and nucleus fields survive.
function evolve(mesh,surface,active,charts,cursor=0,tolerance=.04){
 const A=typeof module!=='undefined'&&module.exports?require('./audit.js'):root.TopologiaAudit,{normal,dot,mean}=T.math,topo=T.topology(mesh),edges=[...topo.edges.values()].filter(e=>e.faces.length===2),v=mesh.vertices;
 for(let index=0;index<Math.min(48,edges.length);index++){
  const e=edges[(cursor+index)%edges.length],[a,b]=e.faces;if(!active[a]||!active[b])continue;
  const fa=mesh.faces[a],fb=mesh.faces[b],n=normal(v,fa);if(dot(n,normal(v,fb))<.866)continue;
  const border=[];for(const f of [fa,fb])for(let k=0;k<4;k++)if(key(f[k],f[(k+1)%4])!==key(e.a,e.b))border.push([f[k],f[(k+1)%4]]);
  const next=new Map(border),ring=[border[0][0]];for(let k=1;k<6;k++)ring.push(next.get(ring.at(-1)));if(new Set(ring).size!==6||ring.some(i=>i===undefined)||next.get(ring[5])!==ring[0])continue;
  const refs=new Set(ring.flatMap(i=>[...charts[i]])),floor=Math.min(T.quality(v,fa),T.quality(v,fb)),oldArea=T.faceArea(v,fa)+T.faceArea(v,fb);
  for(let k=0;k<3;k++){
   const pair=[[0,1,2,3],[3,4,5,0]].map(ids=>ids.map(j=>ring[(k+j)%6]));if(topo.edges.has(key(pair[0][0],pair[0][3])))continue;
   if(Math.min(...pair.map(f=>T.quality(v,f)))<floor+.015||pair.some(f=>dot(normal(v,f),n)<.8)||Math.abs(pair.reduce((s,f)=>s+T.faceArea(v,f),0)/oldArea-1)>.04)continue;
   const samples=pair.flatMap(f=>[mean(f.map(i=>v[i])),...T.triangulate(v,f).map(t=>mean(t.map(i=>v[i])))]);if(samples.some(p=>surface.sample(p,i=>refs.has(i)).distance>Math.min(2*tolerance,.04*Math.sqrt(oldArea))))continue;
   mesh.faces[a]=pair[0];mesh.faces[b]=pair[1];const audit=A.inspect(mesh);
   if(audit.valid&&audit.chi===topo.chi&&audit.components===topo.components)return {changedEdges:1,faces:[a,b],vertices:ring,cursor:cursor+index+1};
   mesh.faces[a]=fa;mesh.faces[b]=fb;return {changedEdges:0,cursor:cursor+index+1};
  }
 }
 return {changedEdges:0,cursor:cursor+48};
}
const api={convert,regularize,coarsen,evolve};if(typeof module!=='undefined'&&module.exports)module.exports=api;else root.TopologiaQuads=api;
})(globalThis);
