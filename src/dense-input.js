(function(root){
'use strict';
function prepare(mesh,options={}){
 if(mesh.faces.length<=150000)return mesh;
 const progress=options.onProgress||(()=>{});progress('Przygotowanie gęstego wejścia: sąsiedztwo i chronione krawędzie…');const inputVertices=mesh.vertices.length,inputTriangles=mesh.faces.length,v=mesh.vertices,f=mesh.faces,valid=new Uint8Array(f.length).fill(1),parent=Uint32Array.from({length:v.length},(_,i)=>i),incident=v.map(()=>new Set()),pinned=new Uint8Array(v.length),lo=[Infinity,Infinity,Infinity],hi=[-Infinity,-Infinity,-Infinity];
 for(const p of v)for(let k=0;k<3;k++){lo[k]=Math.min(lo[k],p[k]);hi[k]=Math.max(hi[k],p[k]);}
 const extent=Math.max(...hi.map((x,k)=>x-lo[k])),tolerance=options.tolerance??.04,cell=Math.min(.012,tolerance/(2*Math.sqrt(3)))*extent;
 if(!(cell>0))throw Error('Dla gęstego wejścia ustaw dodatnią tolerancję geometrii.');
 const grid=Math.ceil(extent/cell)+2,codes=Float64Array.from(v,p=>{const a=p.map((x,k)=>Math.floor((x-lo[k])/cell));return (a[0]*grid+a[1])*grid+a[2];});
 f.forEach((face,id)=>{if(face.length!==3)throw Error('Przygotowanie gęstego wejścia oczekuje trójkątów.');for(const i of face)incident[i].add(id);});
 const normal=face=>{const a=v[face[0]],b=v[face[1]],c=v[face[2]],x=b[0]-a[0],y=b[1]-a[1],z=b[2]-a[2],u=c[0]-a[0],w=c[1]-a[1],t=c[2]-a[2],n=[y*t-z*w,z*u-x*t,x*w-y*u],d=Math.hypot(...n);return d>1e-7*Math.max(x*x+y*y+z*z,u*u+w*w+t*t,extent*extent*1e-12)?n.map(x=>x/d):[0,0,0];},dot=(a,b)=>a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
 let inputEdges=0,badEdges=0;for(let id=0;id<f.length;id++)for(let k=0;k<3;k++){const a=f[id][k],b=f[id][(k+1)%3],shared=[...incident[a]].filter(j=>f[j].includes(b));if(shared[0]!==id)continue;inputEdges++;if(shared.length!==2){pinned[a]=pinned[b]=1;if(shared.length>2)badEdges++;}else if(dot(normal(f[shared[0]]),normal(f[shared[1]]))<Math.cos(50*Math.PI/180))pinned[a]=pinned[b]=1;}
 const neighbors=i=>{const out=new Set();for(const id of incident[i])if(valid[id])for(const j of f[id])if(j!==i)out.add(j);return out;};
 let count=f.length,collapses=0;
 for(let pass=0;pass<3;pass++){progress('Przygotowanie gęstego wejścia: bezpieczne kontrakcje · przebieg '+(pass+1)+'…');const before=count;for(let id=0;id<f.length&&count>80000;id++){
  if(!valid[id])continue;
  for(let k=0;k<3&&valid[id]&&count>80000;k++){
   let a=f[id][k],b=f[id][(k+1)%3];if(codes[a]!==codes[b]||pinned[a]&&pinned[b])continue;if(pinned[b])[a,b]=[b,a];
   const shared=[...incident[a]].filter(j=>valid[j]&&f[j].includes(b));if(shared.length!==2)continue;
   const na=neighbors(a),nb=neighbors(b),opposite=new Set(shared.flatMap(j=>f[j].filter(i=>i!==a&&i!==b))),common=[...na].filter(i=>nb.has(i));if(common.length!==opposite.size||common.some(i=>!opposite.has(i)))continue;
   const region=new Set([...incident[a],...incident[b]]),changes=[];let safe=true;
   for(const j of region){if(!valid[j]||shared.includes(j))continue;const next=f[j].map(i=>i===b?a:i);if(dot(normal(next),normal(f[j]))<.8){safe=false;break;}changes.push([j,next]);}
   if(!safe)continue;
   for(const j of region)if(valid[j])for(const i of f[j])incident[i].delete(j);
   for(const j of shared){valid[j]=0;count--;}
   for(const[j,face]of changes){f[j]=face;for(const i of face)incident[i].add(j);}parent[b]=a;incident[b].clear();collapses++;
  }
 }if(before===count)break;}
 const rootOf=i=>{let p=i;while(parent[p]!==p)p=parent[p];while(parent[i]!==i){const next=parent[i];parent[i]=p;i=next;}return p;};let displacement=0;
 for(let i=0;i<v.length;i++){const p=v[rootOf(i)];displacement=Math.max(displacement,Math.hypot(v[i][0]-p[0],v[i][1]-p[1],v[i][2]-p[2]));}
 const used=new Set();for(let id=0;id<f.length;id++)if(valid[id])for(const i of f[id])used.add(i);const ids=[...used],map=new Map(ids.map((id,i)=>[id,i])),faces=[];
 for(let id=0;id<f.length;id++)if(valid[id])faces.push(f[id].map(i=>map.get(i)));
 return {vertices:ids.map(i=>v[i]),faces,importInfo:mesh.importInfo,transform:mesh.transform,denseInput:{inputVertices,inputTriangles,proxyVertices:ids.length,proxyTriangles:faces.length,collapses,errorBound:displacement/extent,originalExtent:extent,inputChi:inputVertices-inputEdges+inputTriangles,inputNonManifoldEdges:badEdges,method:'feature-pinned edge contractions with link condition'}};
}
const api={prepare};if(typeof module!=='undefined'&&module.exports)module.exports=api;else root.TopologiaDense=api;
})(globalThis);
