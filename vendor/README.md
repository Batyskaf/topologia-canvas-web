# FBX importer (offline)

`fbx-loader.js` bundles Three.js 0.180.0 FBXLoader, NURBSCurve and fflate 0.8.2.
Sources: https://github.com/mrdoob/three.js/tree/r180 and https://github.com/101arrowz/fflate/tree/v0.8.2.
Licenses: THREE-LICENSE.txt and FFLATE-LICENSE.txt (MIT).

Entry source: `../src/fbx-adapter.mjs`. Runtime uses this existing IIFE script, including under file://; no CDN or module loader is needed. Textures are skipped, without external requests. FBX imports static/base geometry with object world transforms. Animation, skin deformation and material/texture transfer are excluded.

To rebuild only (unnecessary for running the app): install dev dependencies with `npm install`, then `npm run build:fbx`. The output is bundled to this folder. Tests exercise both ASCII FBX and compressed binary arrays in FBX 7400 / 7500.

The build postprocessor skips FBX image, texture and material parsing entirely. Those optional metadata blocks must not prevent importing geometry; regression tests cover non-string FBX image filenames.

Binary leaf records containing a trailing null record are recognized as single-property nodes. The official parser previously left arrays wrapped without the expected .a value for this encoding. Regression fixtures cover FBX 7400 and 7500.
