(function(root){
'use strict';
// Cut ambiguous non-manifold fans on an in-memory copy. Positions and polygons
// are preserved; new boundaries are explicitly reported to the user.
function manifoldCopy(input,depth=0){
 const vertices=input.vertices.map(v=>v.slice()),faces=input.faces.map(f=>f.slice()),key=(a,b)=>a<b?`${a}:${b}`:`${b}:${a}`;
 const edges=new Map();faces.forEach((f,i)=>f.forEach((a,k)=>{const id=key(a,f[(k+1)%f.length]);if(!edges.has(id))edges.set(id,[]);edges.get(id).push(i);}));
 const detach=new Set();let badEdges=0;for(const ids of edges.values())if(ids.length>2){badEdges++;for(const i of ids.slice(2))detach.add(i);}
 for(const id of detach)faces[id]=faces[id].map(v=>{const n=vertices.length;vertices.push(vertices[v].slice());return n;});
 const e2=new Map(),incident=vertices.map(()=>[]);faces.forEach((f,i)=>{for(const v of f)incident[v].push(i);f.forEach((a,k)=>{const id=key(a,f[(k+1)%f.length]);if(!e2.has(id))e2.set(id,[]);e2.get(id).push(i);});});
 const adjacent=faces.map(()=>[]);for(const ids of e2.values())if(ids.length===2){adjacent[ids[0]].push(ids[1]);adjacent[ids[1]].push(ids[0]);}
 let splitVertices=0;for(let v=0;v<incident.length;v++){const fs=incident[v],allowed=new Set(fs),seen=new Set();let group=0;for(const first of fs)if(!seen.has(first)){const fan=[],todo=[first];seen.add(first);while(todo.length){const i=todo.pop();fan.push(i);for(const j of adjacent[i])if(allowed.has(j)&&!seen.has(j)&&faces[j].includes(v)){seen.add(j);todo.push(j);}}if(group++){const next=vertices.length;vertices.push(vertices[v].slice());for(const i of fan)faces[i]=faces[i].map(j=>j===v?next:j);splitVertices++;}}}
 const orientedEdges=new Map(),constraints=faces.map(()=>[]);faces.forEach((f,i)=>f.forEach((a,k)=>{const b=f[(k+1)%f.length],id=key(a,b);if(!orientedEdges.has(id))orientedEdges.set(id,[]);orientedEdges.get(id).push({face:i,direction:a<b?1:-1});}));for(const rows of orientedEdges.values())if(rows.length===2){const[a,b]=rows,r=-a.direction*b.direction;constraints[a.face].push([b.face,r]);constraints[b.face].push([a.face,r]);}
 const signs=new Int8Array(faces.length),conflicts=new Set();for(let i=0;i<faces.length;i++)if(!signs[i]){signs[i]=1;const todo=[i];while(todo.length){const a=todo.pop();for(const[b,r]of constraints[a]){if(signs[b]&&signs[b]!==signs[a]*r)conflicts.add(Math.max(a,b));else if(!signs[b]){signs[b]=signs[a]*r;todo.push(b);}}}}
 if(conflicts.size){if(depth>=6)throw Error('Nie można rozciąć konfliktów orientacji źródła.');for(const id of conflicts)faces[id]=faces[id].map(v=>{const n=vertices.length;vertices.push(vertices[v].slice());return n;});const next=manifoldCopy({vertices,faces},depth+1);next.repair.badEdges+=badEdges;next.repair.detachedFaces+=detach.size+conflicts.size;next.repair.splitVertices+=splitVertices;next.repair.orientationCuts=(next.repair.orientationCuts||0)+conflicts.size;return next;}
 const used=[...new Set(faces.flat())].sort((a,b)=>a-b),map=new Map(used.map((id,i)=>[id,i]));return {vertices:used.map(id=>vertices[id]),faces:faces.map(f=>f.map(id=>map.get(id))),repair:{badEdges,detachedFaces:detach.size,splitVertices,orientationCuts:0}};
}
function previewCopy(input){const lo=[Infinity,Infinity,Infinity],hi=[-Infinity,-Infinity,-Infinity];for(const v of input.vertices){if(v.length!==3||!v.every(Number.isFinite))throw Error('Nieprawidłowe współrzędne wejścia.');for(let k=0;k<3;k++){lo[k]=Math.min(lo[k],v[k]);hi[k]=Math.max(hi[k],v[k]);}}const center=lo.map((x,k)=>(x+hi[k])/2),scale=Math.max(...hi.map((x,k)=>x-lo[k]))/2;if(!(scale>0))throw Error('Zerowy rozmiar mesha.');return {vertices:input.vertices.map(v=>v.map((x,k)=>(x-center[k])/scale)),faces:input.faces.map(f=>f.slice()),transform:{center,scale}};}
function repairLocalIntersections(source,tolerance,initialReport=null){
 const T=typeof module!=='undefined'&&module.exports?require('./engine.js'):root.Topologia,A=typeof module!=='undefined'&&module.exports?require('./audit.js'):root.TopologiaAudit;
 const topo=T.topology(source),vertices=source.vertices.map(v=>v.slice()),neighbors=vertices.map(()=>new Set()),normals=source.faces.map(f=>T.math.normal(source.vertices,f)),cap=Math.min(.02,2*tolerance*.25);
 for(const e of topo.edges.values()){neighbors[e.a].add(e.b);neighbors[e.b].add(e.a);}
 const pinned=new Set();for(const e of topo.edges.values())if(e.faces.length===1){pinned.add(e.a);pinned.add(e.b);}
 let report=initialReport||A.inspect(source,1500000);const before=report.intersectionPairs.length;
 for(let round=0;round<6&&!report.valid&&report.intersectionPairs.length;round++){
  const selected=new Set(report.intersectionPairs.flatMap(pair=>pair.flatMap(id=>source.faces[id])));for(const id of [...selected])for(const j of neighbors[id])selected.add(j);
  const changes=[];for(const id of selected){if(pinned.has(id)||!neighbors[id].size)continue;const avg=T.math.mean([...neighbors[id]].map(i=>vertices[i])),next=T.math.add(vertices[id],T.math.mul(T.math.sub(avg,vertices[id]),.15));if(T.math.len(T.math.sub(next,source.vertices[id]))>cap)continue;
   if(topo.incident[id].every(face=>{const f=source.faces[face],ps=f.map(i=>i===id?next:vertices[i]);return T.math.dot(T.math.normal(ps,ps.map((_,i)=>i)),normals[face])>.8;}))changes.push([id,next]);
  }
  if(!changes.length)break;for(const[id,p]of changes)vertices[id]=p;report=A.inspect({vertices,faces:source.faces},1500000);
 }
 const mesh={vertices,faces:source.faces.map(f=>f.slice()),transform:source.transform};
 if(!report.valid||mesh.faces.some((f,i)=>T.math.dot(T.math.normal(vertices,f),normals[i])<.5))throw Object.assign(Error('Nie udało się naprawić przecięć lokalnie w granicach tolerancji. Model pozostaje w podglądzie.'),{code:'INPUT_INTERSECTIONS',details:{pairs:before,remaining:report.intersectionPairs.length}});
 const fit=A.compare(source,mesh,1);if(Math.max(fit.sourceToResult.max,fit.resultToSource.max)>2*tolerance)throw Error('Naprawa przekroczyłaby tolerancję kształtu.');
 const moved=vertices.map((v,i)=>T.math.len(T.math.sub(v,source.vertices[i])));return {mesh,repair:{intersectionPairsBefore:before,movedVertices:moved.filter(d=>d>1e-12).length,maxDisplacement:moved.reduce((maximum,d)=>Math.max(maximum,d),0),fit}};
}
// Excise one side of each unresolved crossing. This is an open-surface repair,
// not hole filling: no invented cap may intersect the remaining surface.
function repairIntersections(source,tolerance,initialReport=null,onProgress=()=>{}){
 const T=typeof module!=='undefined'&&module.exports?require('./engine.js'):root.Topologia,A=typeof module!=='undefined'&&module.exports?require('./audit.js'):root.TopologiaAudit;
 let report=initialReport||A.inspect(source,1500000);
 if(report.intersectionPairs.length<20){try{return repairLocalIntersections(source,tolerance,report);}catch(e){if(e.code!=='INPUT_INTERSECTIONS')throw e;}}
 let mesh=source,removedFaces=0;const initialPairs=report.intersectionPairs.length,initialComplete=report.complete,maximum=Math.max(8,Math.floor(source.faces.length*.01));
 for(let round=0;round<128&&report.intersectionPairs.length;round++){
  onProgress(`Naprawa przecięć: rozcinanie powierzchni · ${removedFaces} usuniętych ścian…`);
  const remove=new Set();for(const[a,b]of report.intersectionPairs){if(remove.has(a)||remove.has(b))continue;remove.add(T.faceArea(mesh.vertices,mesh.faces[a])<=T.faceArea(mesh.vertices,mesh.faces[b])?a:b);}
  if(removedFaces+remove.size>maximum)throw Error('Naprawa przecięć wymaga usunięcia zbyt dużej części powierzchni (limit 1%). Model pozostaje w podglądzie.');
  removedFaces+=remove.size;const copy=manifoldCopy({vertices:mesh.vertices,faces:mesh.faces.filter((_,i)=>!remove.has(i))});
  if(!copy.faces.length)throw Error('Naprawa usunęłaby całą powierzchnię.');
  const oriented=T.prepare(copy);mesh={vertices:oriented.vertices.map(v=>T.math.add(T.math.mul(v,oriented.transform.scale),oriented.transform.center)),faces:oriented.faces,transform:source.transform};report=A.inspect(mesh,1500000);
 }
 if(!report.valid)throw Error('Naprawa przecięć nie przeszła pełnej kontroli: '+report.errors.join(' '));
 // T.prepare may normalize again after cuts. Compare in the original frame.
 const aligned={vertices:mesh.vertices,faces:mesh.faces};
 const fit=A.compare(source,aligned,1);if(Math.max(fit.sourceToResult.max,fit.resultToSource.max)>2*tolerance)throw Error('Rozcięcie przecięć przekroczyłoby tolerancję kształtu. Zwiększ tolerancję lub napraw mesh w edytorze.');
 const repaired={vertices:aligned.vertices,faces:aligned.faces,transform:source.transform};
 return {mesh:repaired,repair:{method:'excision',intersectionPairsBefore:initialPairs,initialPairsLowerBound:!initialComplete,removedFaces,movedVertices:0,maxDisplacement:0,boundaryLoops:report.boundaryLoops,boundaryEdges:report.boundaryEdges,fit}};
}
function build(input,levels,params,options={}){
 const T=typeof module!=='undefined'&&module.exports?require('./engine.js'):root.Topologia,A=typeof module!=='undefined'&&module.exports?require('./audit.js'):root.TopologiaAudit,progress=options.onProgress||(()=>{});
 const raw=input.transform?{vertices:input.vertices.map(v=>T.math.add(T.math.mul(v,input.transform.scale),input.transform.center)),faces:input.faces}:input,original=previewCopy(raw);let prepared,repair=null,geometryRepair=null;
 progress('Sprawdzanie topologii wejścia…');try{prepared=T.prepare(raw);}catch(e){if(!options.repairTopology||!/non-manifold|nieorientowalna/.test(e.message))throw e;const copy=manifoldCopy(raw);prepared=T.prepare(copy);repair=copy.repair;}
 progress('Sprawdzanie przecięć wejścia…');const diagnostic=A.inspect(prepared,1500000);
 if(!diagnostic.valid){
  if(!options.repairIntersections||!diagnostic.intersectionPairs.length)throw Object.assign(Error(diagnostic.intersectionPairs.length?`Wejście zawiera ${diagnostic.complete?'':'co najmniej '}${diagnostic.intersectionPairs.length} par przecinających się ścian. Włącz naprawę drobnych przecięć na kopii.`:diagnostic.errors.join(' ')),{code:'INPUT_GEOMETRY',details:diagnostic});
  progress('Naprawa przecięć na kopii…');const fixed=repairIntersections(prepared,input.denseInput?(params.maxSurfaceError-input.denseInput.errorBound)*input.denseInput.originalExtent/(2*prepared.transform.scale):params.maxSurfaceError,diagnostic,progress);prepared=fixed.mesh;geometryRepair=fixed.repair;
 }
 const dense=input.denseInput;let effectiveParams=params;if(dense){const remaining=params.maxSurfaceError-dense.errorBound;if(!(remaining>0))throw Error('Przygotowanie gęstego wejścia wykorzystało całą tolerancję. Zwiększ tolerancję lub użyj lżejszego modelu.');effectiveParams={...params,maxSurfaceError:remaining*dense.originalExtent/(2*prepared.transform.scale)};}
 progress('Redukcja z ochroną cech i kontrolą kolizji…');const sim=new T.Simulation(prepared,levels,effectiveParams);if(dense){sim.denseInput=dense;sim.requestedMaxSurfaceError=params.maxSurfaceError;}
 if(geometryRepair){sim.inputOriginal=original;sim.importRepair=geometryRepair;progress('Sprawdzanie wyniku względem oryginalnego pliku…');const fit=sim.compareOriginal({vertices:sim.vertices,faces:sim.faces},1);if(Math.max(fit.sourceToResult.max,fit.resultToSource.max)>2*(dense?(params.maxSurfaceError-dense.errorBound)*dense.originalExtent/(2*original.transform.scale):params.maxSurfaceError))throw Error('Wynik przekracza tolerancję względem oryginału sprzed naprawy. Zwiększ procent pozostawionych poligonów.');sim.originalFit=fit;}
 return {sim,repair,geometryRepair,diagnostic};
}
const api={manifoldCopy,previewCopy,repairLocalIntersections,repairIntersections,build};if(typeof module!=='undefined'&&module.exports)module.exports=api;else root.TopologiaImport=api;
})(globalThis);
