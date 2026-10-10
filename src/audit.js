(function(root){
'use strict';const T=typeof module!=='undefined'&&module.exports?require('./engine.js'):root.Topologia,{sub,dot,cross,len,unit,mean,normal}=T.math;
const EPS=1e-9;
function boxOverlap(a,b){return a.lo.every((x,i)=>x<=b.hi[i]+EPS&&a.hi[i]+EPS>=b.lo[i]);}
function triangleBox(p){return {lo:[0,1,2].map(k=>Math.min(...p.map(v=>v[k]))),hi:[0,1,2].map(k=>Math.max(...p.map(v=>v[k])))};}
function segmentHits(a,b,t){const d=sub(b,a),u=sub(t[1],t[0]),v=sub(t[2],t[0]),p=cross(d,v),det=dot(u,p);if(Math.abs(det)<EPS*len(d)*len(u)*len(v))return false;const inv=1/det,s=sub(a,t[0]),x=dot(s,p)*inv,q=cross(s,u),y=dot(d,q)*inv,z=dot(v,q)*inv;return x>EPS&&y>EPS&&x+y<1-EPS&&z>EPS&&z<1-EPS;}
function trianglesCross(a,b){
 const n=unit(cross(sub(a[1],a[0]),sub(a[2],a[0]))),m=unit(cross(sub(b[1],b[0]),sub(b[2],b[0])));
 if(Math.abs(dot(n,m))>1-1e-8&&b.every(p=>Math.abs(dot(sub(p,a[0]),n))<EPS)){
  const axis=n.map(Math.abs).indexOf(Math.max(...n.map(Math.abs))),axes=[0,1,2].filter(k=>k!==axis),aa=a.map(v=>axes.map(k=>v[k])),bb=b.map(v=>axes.map(k=>v[k]));
  const orient=(p,q,r)=>(q[0]-p[0])*(r[1]-p[1])-(q[1]-p[1])*(r[0]-p[0]),inside=(p,t)=>{const s=Math.sign(orient(...t));return t.every((v,i)=>s*orient(v,t[(i+1)%3],p)>EPS);};
  for(let i=0;i<3;i++)for(let j=0;j<3;j++){const p=aa[i],q=aa[(i+1)%3],r=bb[j],s=bb[(j+1)%3];if(orient(p,q,r)*orient(p,q,s)<-(EPS**2)&&orient(r,s,p)*orient(r,s,q)<-(EPS**2))return true;}
  return aa.some(p=>inside(p,bb))||bb.some(p=>inside(p,aa))||inside(mean(a).filter((_,k)=>k!==axis),bb)||inside(mean(b).filter((_,k)=>k!==axis),aa);
 }
 for(let i=0;i<3;i++)if(segmentHits(a[i],a[(i+1)%3],b)||segmentHits(b[i],b[(i+1)%3],a))return true;return false;
}
function inspect(mesh,limit=1500000){
 const report={valid:false,complete:false,faces:mesh.faces.length,vertices:mesh.vertices.length,triangles:mesh.faces.reduce((s,f)=>s+f.length-2,0),errors:[],intersectionPairs:[],testedPairs:0};let topo;
 try{topo=T.topology(mesh);}catch(e){report.errors.push(e.message);return report;}
 const used=new Set(mesh.faces.flat());if(used.size!==mesh.vertices.length)report.errors.push('Nieużywane wierzchołki.');if(mesh.faces.some(f=>T.faceArea(mesh.vertices,f)<1e-14||(f.length===4&&T.quality(mesh.vertices,f)<=0)))report.errors.push('Zerowe pole lub odwrócone quady.');
 for(const e of topo.edges.values())if(e.faces.length===2){const dir=f=>{const k=f.indexOf(e.a);return f[(k+1)%f.length]===e.b?1:-1;};if(dir(mesh.faces[e.faces[0]])===dir(mesh.faces[e.faces[1]])){report.errors.push('Niezgodna orientacja przy wspólnej krawędzi.');break;}}
 const boundaryGraph=new Map(),degree=new Int32Array(mesh.vertices.length);for(const e of topo.edges.values()){degree[e.a]++;degree[e.b]++;if(e.faces.length===1){for(const[a,b]of [[e.a,e.b],[e.b,e.a]]){if(!boundaryGraph.has(a))boundaryGraph.set(a,[]);boundaryGraph.get(a).push(b);}}}let loops=0;const seen=new Set();for(const id of boundaryGraph.keys())if(!seen.has(id)){loops++;const stack=[id];seen.add(id);while(stack.length)for(const v of boundaryGraph.get(stack.pop()))if(!seen.has(v)){seen.add(v);stack.push(v);}}if([...boundaryGraph.values()].some(ns=>ns.length!==2))report.errors.push('Granica nie tworzy rozłącznych pętli.');
 Object.assign(report,{chi:topo.chi,components:topo.components,boundaryLoops:loops,boundaryEdges:topo.boundary,poles:Array.from(degree).reduce((s,d,i)=>s+(!boundaryGraph.has(i)&&d!==4?1:0),0)});
 let surface;try{surface=new T.Surface(mesh);}catch(e){report.errors.push('Nie można spójnie triangulować ścian: '+e.message);report.complete=true;return report;}const tri=surface.tri.map(f=>f.map(i=>mesh.vertices[i])),boxes=tri.map(triangleBox);
 const overlapping=(node,box,out)=>{if(!boxOverlap(node,box))return;if(node.ids)out.push(...node.ids);else{overlapping(node.left,box,out);overlapping(node.right,box,out);}};
 let stopped=false;for(let i=0;i<tri.length&&!stopped;i++){const candidates=[];overlapping(surface.tree,boxes[i],candidates);for(const j of candidates){if(j<=i||surface.triangleFaces[i]===surface.triangleFaces[j]||!boxOverlap(boxes[i],boxes[j]))continue;if(++report.testedPairs>limit){stopped=true;report.errors.push('Limit kontroli przecięć: wynik niezweryfikowany.');break;}if(trianglesCross(tri[i],tri[j])){const pair=[surface.triangleFaces[i],surface.triangleFaces[j]];if(!report.intersectionPairs.some(p=>p[0]===pair[0]&&p[1]===pair[1]))report.intersectionPairs.push(pair);if(report.intersectionPairs.length>=20){stopped=true;break;}}}}
 if(report.intersectionPairs.length)report.errors.push('Wykryto przecinające się ściany.');report.complete=!stopped;report.valid=report.complete&&!report.errors.length;return report;
}
function compare(source,result,scale=1){
 const inputSurface=new T.Surface(source),outputSurface=new T.Surface(result);
 const sample=(a,b,surface,ownSurface)=>{let total=0,weighted=0,normalSum=0,max=0;for(const f of ownSurface.tri){const area=T.faceArea(a.vertices,f),p=mean(f.map(i=>a.vertices[i])),hit=surface.sample(p),n=normal(a.vertices,f),m=normal(b.vertices,surface.tri[hit.triangle]);weighted+=area*hit.distance**2;normalSum+=area*(1-Math.max(-1,Math.min(1,dot(n,m))));total+=area;max=Math.max(max,hit.distance);}for(const p of a.vertices)max=Math.max(max,surface.sample(p).distance);return {rms:Math.sqrt(weighted/total)*scale,max:max*scale,normalError:normalSum/total,samples:ownSurface.tri.length+a.vertices.length};};
 return {sourceToResult:sample(source,result,outputSurface,inputSurface),resultToSource:sample(result,source,inputSurface,outputSurface),note:'Próbkowany błąd obustronny (środki ścian ważone polem i maksima na wierzchołkach); nie certyfikowana odległość Hausdorffa.'};
}
class CollisionIndex{
 constructor(mesh){this.mesh=mesh;this.surface=new T.Surface(mesh);this.faceTriangles=mesh.faces.map(()=>[]);this.surface.triangleFaces.forEach((face,id)=>this.faceTriangles[face].push(id));this.leaf=[];const visit=(node,parent=null)=>{node.parent=parent;if(node.ids)for(const id of node.ids)this.leaf[id]=node;else{visit(node.left,node);visit(node.right,node);}};visit(this.surface.tree);}
 refreshFaces(faces){
  const replacements=[];try{for(const face of faces){const parts=T.triangulate(this.mesh.vertices,this.mesh.faces[face]);if(parts.length!==this.faceTriangles[face].length)return false;parts.forEach((part,k)=>replacements.push([this.faceTriangles[face][k],part]));}}catch{return false;}
  const leaves=new Set();for(const[id,part]of replacements){this.surface.tri[id]=part;leaves.add(this.leaf[id]);}
  for(const leaf of leaves){leaf.lo=[Infinity,Infinity,Infinity];leaf.hi=[-Infinity,-Infinity,-Infinity];for(const id of leaf.ids)for(const v of this.surface.tri[id])for(let k=0;k<3;k++){leaf.lo[k]=Math.min(leaf.lo[k],this.mesh.vertices[v][k]);leaf.hi[k]=Math.max(leaf.hi[k],this.mesh.vertices[v][k]);}}
  for(const leaf of leaves){let node=leaf.parent;while(node){node.lo=node.left.lo.map((x,k)=>Math.min(x,node.right.lo[k]));node.hi=node.left.hi.map((x,k)=>Math.max(x,node.right.hi[k]));node=node.parent;}}return true;
 }
 canMove(faces){
  if(!this.refreshFaces(faces))return false;const s=this.surface,points=id=>s.tri[id].map(i=>this.mesh.vertices[i]),overlap=(node,box,out)=>{if(!boxOverlap(node,box))return;if(node.ids)out.push(...node.ids);else{overlap(node.left,box,out);overlap(node.right,box,out);}};
  for(const face of faces)for(const id of this.faceTriangles[face]){const a=points(id),box=triangleBox(a),ids=[];overlap(s.tree,box,ids);for(const other of ids)if(s.triangleFaces[other]!==face&&(!this.mesh.active||this.mesh.active[s.triangleFaces[other]])){const shared=s.tri[id].filter(i=>s.tri[other].includes(i));if(shared.length===2){const direction=tri=>tri[(tri.indexOf(shared[0])+1)%3]===shared[1]?1:-1;if(direction(s.tri[id])!==direction(s.tri[other])){const b=points(other),na=cross(sub(a[1],a[0]),sub(a[2],a[0])),nb=cross(sub(b[1],b[0]),sub(b[2],b[0]));if(dot(na,nb)>0)continue;}}const b=points(other);if(boxOverlap(box,triangleBox(b))&&trianglesCross(a,b))return false;}}return true;
 }
}
const api={inspect,compare,trianglesCross,CollisionIndex};if(typeof module!=='undefined'&&module.exports)module.exports=api;else root.TopologiaAudit=api;
})(globalThis);
