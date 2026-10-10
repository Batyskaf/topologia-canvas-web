importScripts('engine.js','audit.js','reduction.js','quadrangulate.js','import-mesh.js','dense-input.js','../vendor/fbx-loader.js');
self.onmessage=({data})=>{
 try{
  if(data.action==='parse'){
   let mesh=data.format==='fbx'?TopologiaFBX.parseFBX(data.content,data.options):data.format==='obj'?Topologia.parseOBJ(data.content):data.format==='stl'?Topologia.parseSTL(data.content):JSON.parse(data.content);
   mesh=TopologiaDense.prepare(mesh,{tolerance:data.options?.denseTolerance??.04,onProgress:progress=>self.postMessage({progress})});self.postMessage({mesh});return;
  }
  const result=TopologiaImport.build(data.mesh,data.levels,data.params,{repairTopology:data.repair,repairIntersections:data.repairIntersections,onProgress:progress=>self.postMessage({progress})});self.postMessage(result);
 }catch(e){self.postMessage({error:e.message,code:e.code,details:e.details});}
};
