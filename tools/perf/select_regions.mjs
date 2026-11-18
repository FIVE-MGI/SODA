// Choose circles on a Visium HD 8 um layer that hold ~N bins the UMAP will use (Draw a region timings).
//   node select_regions.mjs <export dir> <layer> <out json> 1000 10000 50000
// Reads (locally, read only) expr/<layer>/positions.bin (float32 x,y per bin row) and the hvg
// chunks (per bin: row, total UMI), then, around one fixed centre (the median of the bins with
// >= 10 UMI), finds for each N the largest radius whose 72-gon (as viewer.js circlePts / inPoly)
// holds <= N bins with >= minUmi UMI: the count lasso_worker.js keeps (no QC filter).
import { readFileSync, readdirSync, writeFileSync } from "node:fs";
import { join } from "node:path";

const [dir, layer, out, ...ns] = process.argv.slice(2);
const N = ns.map(Number);
const MIN_UMI = 10;
const pb = readFileSync(join(dir, "expr", layer, "positions.bin"));
const pos = new Float32Array(pb.buffer, pb.byteOffset, pb.length / 4);
const nBins = pos.length / 2;
const lib = new Float32Array(nBins).fill(-1);
const hvgDir = join(dir, "expr", layer, "hvg");
const index = JSON.parse(readFileSync(join(hvgDir, "index.json"), "utf8"));
for (const f of readdirSync(hvgDir).filter((f) => f.endsWith(".bin"))) {
  const b = readFileSync(join(hvgDir, f));
  const buf = b.buffer.slice(b.byteOffset, b.byteOffset + b.length);
  const n = new DataView(buf).getUint32(0, true);
  const row = new Uint32Array(buf.slice(4, 4 + 4 * n)), libs = new Float32Array(buf.slice(4 + 4 * n, 4 + 8 * n));
  for (let i = 0; i < n; i++) lib[row[i]] = libs[i];
}
const ok = [];
for (let i = 0; i < nBins; i++) if (lib[i] >= MIN_UMI) ok.push(i);
const med = (a) => { const s = Float64Array.from(a).sort(); return s[Math.floor(s.length / 2)]; };
const cx = med(ok.map((i) => pos[2 * i])), cy = med(ok.map((i) => pos[2 * i + 1]));

function circlePts(c, r) {   // = viewer.js circlePts
  const pts = [];
  for (let k = 0; k < 72; k++) { const a = (k / 72) * 2 * Math.PI; pts.push([c[0] + r * Math.cos(a), c[1] + r * Math.sin(a)]); }
  return pts;
}
function inPoly(x, y, pts) {  // = viewer.js inPoly
  let inside = false;
  for (let i = 0, j = pts.length - 1; i < pts.length; j = i++) {
    const [xi, yi] = pts[i], [xj, yj] = pts[j];
    if ((yi > y) !== (yj > y) && x < xi + ((xj - xi) * (y - yi)) / (yj - yi)) inside = !inside;
  }
  return inside;
}
function count(r, given) {
  const pts = given || circlePts([cx, cy], r);
  let all = 0, kept = 0, chunks = new Set();
  for (let i = 0; i < nBins; i++) {
    const x = pos[2 * i], y = pos[2 * i + 1];
    if (Math.abs(x - cx) > r || Math.abs(y - cy) > r || !inPoly(x, y, pts)) continue;
    all++; if (lib[i] >= MIN_UMI) kept++;
    chunks.add(`${Math.floor(x / index.chunkPx)}_${Math.floor(y / index.chunkPx)}`);
  }
  return { all, kept, chunks: chunks.size };
}
const regions = [];
for (const target of N) {
  let lo = 1, hi = 20000;
  for (let it = 0; it < 40; it++) { const m = (lo + hi) / 2; if (count(m).kept <= target) lo = m; else hi = m; }
  const r = Math.floor(lo * 100) / 100;     // the viewer saves points to 0.01 px
  const pts = circlePts([cx, cy], r).map((p) => [Math.round(p[0] * 100) / 100, Math.round(p[1] * 100) / 100]);
  const c = count(r + 1, pts);              // counted on the rounded points actually handed to the viewer
  regions.push({ target, centre: [cx, cy], radiusPx: r, binsInShape: c.all, expectedInUmap: c.kept, chunks: c.chunks, points: pts });
  console.log(`N=${target}: radius ${r} px, ${c.all} bins in shape, ${c.kept} with >= ${MIN_UMI} UMI, ${c.chunks} chunks`);
}
writeFileSync(out, JSON.stringify({ dir, layer, nBins, nOk: ok.length, chunkPx: index.chunkPx, regions }, null, 1));
