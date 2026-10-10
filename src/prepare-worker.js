importScripts('engine.js?v=1.3.7','audit.js?v=1.3.7','reduction.js?v=1.3.7','quadrangulate.js?v=1.3.7','import-mesh.js?v=1.3.7','dense-input.js?v=1.3.7','../vendor/fbx-loader.js?v=1.3.7');
self.onmessage=({data})=>{
 try{
  if(data.action==='parse'){
   let mesh=data.format==='fbx'?TopologiaFBX.parseFBX(data.content,data.options):data.format==='obj'?Topologia.parseOBJ(data.content):data.format==='stl'?Topologia.parseSTL(data.content):JSON.parse(data.content);
   mesh=TopologiaDense.prepare(mesh,{tolerance:data.options?.denseTolerance??.04,onProgress:progress=>self.postMessage({progress})});self.postMessage({mesh});return;
  }
  const result=TopologiaImport.build(data.mesh,data.levels,data.params,{repairTopology:data.repair,repairIntersections:data.repairIntersections,onProgress:progress=>self.postMessage({progress})});self.postMessage(result);
 }catch(e){self.postMessage({error:e.message,code:e.code,details:e.details});}
};
