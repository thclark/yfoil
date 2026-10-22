// The native-versus-WebAssembly bit-identity check (run by scripts/wasm-identity.sh).
//
// Every comparison is of bytes: the JSON records the native `yfoil` CLI writes against the ones
// the wasm module returns for the same inputs (serde_json writes each f64 as its shortest
// round-trip text, so equal text is equal bits), and the raw f64 bits of the maths library.
//
//   node identity.mjs --pkg DIR --yfoil BIN --native BIN --work DIR [--quick] [--samples N]
//
// --quick runs only the CI reference case and fewer samples. Exits 1 on any difference.

import { createRequire } from "node:module";
import { execFileSync } from "node:child_process";
import fs from "node:fs";
import path from "node:path";

const args = Object.fromEntries(
  process.argv.slice(2).reduce((acc, a, i, all) => {
    if (a.startsWith("--")) acc.push([a.slice(2), all[i + 1] && !all[i + 1].startsWith("--") ? all[i + 1] : true]);
    return acc;
  }, []),
);
const quick = args.quick === true;
const samples = Number(args.samples ?? (quick ? 1 << 14 : 1 << 18));
const work = path.resolve(args.work);
fs.mkdirSync(work, { recursive: true });
const yfoil = createRequire(import.meta.url)(path.resolve(args.pkg, "yfoil_wasm.js"));

const REF_CASE = "tests/fixtures/xfoil/naca0012_n60_a2_re1e6/panels.json";
let failures = 0;
const report = (label, ok, detail = "") => {
  if (!ok) failures++;
  console.log(`${ok ? "identical" : "DIFFERENT"}  ${label}${detail ? "  " + detail : ""}`);
};
const native = (cmd) => execFileSync(path.resolve(args.yfoil), cmd, { stdio: ["ignore", "ignore", "inherit"] });
const firstDifference = (a, b) => {
  let i = 0;
  while (i < a.length && a[i] === b[i]) i++;
  return `first difference at byte ${i}: native …${a.slice(Math.max(0, i - 40), i + 40)}… wasm …${b.slice(Math.max(0, i - 40), i + 40)}…`;
};
const compareText = (label, nativeFile, wasmText) => {
  const nativeText = fs.readFileSync(nativeFile, "utf8");
  const ok = nativeText === wasmText;
  report(label, ok, ok ? `(${nativeText.length} bytes)` : firstDifference(nativeText, wasmText));
};
const timed = (f) => {
  const t = process.hrtime.bigint();
  const r = f();
  return [r, Number(process.hrtime.bigint() - t) / 1e6];
};

// 1. geometry generators: `yfoil geometry naca` against naca()
const designations = quick ? ["0012"] : ["0012", "4412", "23012", "23112", "0012-34", "16-212", "63-415", "64A010"];
for (const d of designations) {
  const file = path.join(work, `naca_${d}.json`);
  native(["geometry", "naca", d, "-o", file]);
  compareText(`geometry naca ${d}`, file, yfoil.naca(d, "", NaN));
}

// 2. analyses: the CLI's records against analyse() and polar(), on every tracked panel set
const panelSets = quick
  ? [REF_CASE]
  : fs
      .readdirSync("tests/fixtures/xfoil")
      .map((c) => `tests/fixtures/xfoil/${c}/panels.json`)
      .filter((p) => fs.existsSync(p))
      .sort();
// [label, CLI arguments, binding call, binding options]
const runs = [
  ["analyse a2 (lagged closures)", ["analyse", "-a", "2", "--include-lagged-closures", "--allow-invalid"], "analyse",
    { alpha: 2, include_lagged_closures: true, allow_invalid: true }],
  ["analyse cl0.5 re3e5 m0.3", ["analyse", "--cl", "0.5", "-r", "3e5", "-m", "0.3", "--allow-invalid"], "analyse",
    { cl: 0.5, reynolds: 3e5, mach: 0.3, allow_invalid: true }],
  ["analyse a6 re1e4", ["analyse", "-a", "6", "-r", "1e4", "--allow-invalid"], "analyse",
    { alpha: 6, reynolds: 1e4, allow_invalid: true }],
  ["analyse inviscid a4 m0.6", ["analyse", "-a", "4", "-m", "0.6", "--inviscid", "--allow-invalid"], "analyse",
    { alpha: 4, mach: 0.6, inviscid: true, allow_invalid: true }],
  ["polar -10..20 iter50 (distributions)",
    ["polar", "--alpha-min", "-10", "--alpha-max", "20", "--alpha-step", "1", "--max-iterations", "50",
      "--distributions", "--allow-invalid", "--json"],
    "polar", { alpha_min: -10, alpha_max: 20, alpha_step: 1, max_iterations: 50, distributions: true, allow_invalid: true }],
];
const timing = { native: 0, wasm: 0 };
for (const panels of panelSets) {
  const geometry = fs.readFileSync(panels, "utf8");
  const caseName = path.basename(path.dirname(panels));
  for (const [label, cli, call, options] of quick ? runs.slice(0, 1) : runs) {
    const out = path.join(work, `${caseName}__${label.replace(/[^a-z0-9.]+/gi, "_")}.json`);
    const [, tn] = timed(() => native([cli[0], panels, ...cli.slice(1), "-o", out]));
    // the CLI names the foil after the file stem
    const [text, tw] = timed(() => yfoil[call](geometry, JSON.stringify({ ...options, name: "panels" })));
    timing.native += tn;
    timing.wasm += tw;
    compareText(`${caseName}: ${label}`, out, text);
  }
}

// 3. the maths library, bit for bit, on random and edge-case arguments
let state = 0x9e3779b97f4a7c15n;
const nextBits = () => {
  state ^= (state << 13n) & 0xffffffffffffffffn;
  state ^= state >> 7n;
  state ^= (state << 17n) & 0xffffffffffffffffn;
  return state;
};
const scratch = new DataView(new ArrayBuffer(8));
const fromBits = (b) => (scratch.setBigUint64(0, b), scratch.getFloat64(0));
const uniform = () => Number(nextBits() >> 11n) / 2 ** 53;
const randomBits = () => fromBits(nextBits());
const edges = [0, -0, Infinity, -Infinity, NaN, 1, -1, 0.5, 2, Number.MIN_VALUE, -Number.MIN_VALUE, 2.2250738585072014e-308,
  Number.MAX_VALUE, -Number.MAX_VALUE, Math.PI, Math.PI / 2, Math.PI / 4, 1e-300, 1e300, 709.782712893384, -745.1332191019411,
  1 + 2 ** -52, 1 - 2 ** -53, 1e22, 2 ** 1023, 3 * 2 ** 1022];
// each generator returns one argument tuple; a quarter are arbitrary bit patterns
const domain = {
  exp: () => [uniform() * 1460 - 750],
  ln: () => [uniform() < 0.5 ? Math.abs(randomBits()) : 0.5 + uniform()],
  log10: () => [uniform() < 0.5 ? Math.abs(randomBits()) : 0.5 + uniform()],
  tanh: () => [(uniform() - 0.5) * (uniform() < 0.5 ? 40 : 1e-6)],
  sin: () => [(uniform() - 0.5) * (uniform() < 0.8 ? 20 : 1e20)],
  cos: () => [(uniform() - 0.5) * (uniform() < 0.8 ? 20 : 1e20)],
  "sin_cos.sin": () => [(uniform() - 0.5) * (uniform() < 0.8 ? 20 : 1e20)],
  "sin_cos.cos": () => [(uniform() - 0.5) * (uniform() < 0.8 ? 20 : 1e20)],
  atan: () => [(uniform() - 0.5) * 2 ** (uniform() * 200 - 100)],
  asin: () => [uniform() * 2 - 1],
  pow: () => (uniform() < 0.5 ? [uniform() * 10, (uniform() - 0.5) * 40] : [0.9 + uniform() * 0.2, (uniform() - 0.5) * 2000]),
  atan2: () => [(uniform() - 0.5) * 10, (uniform() - 0.5) * 10],
  hypot: () => [(uniform() - 0.5) * 2 ** (uniform() * 200 - 100), (uniform() - 0.5) * 2 ** (uniform() * 200 - 100)],
  // products that cancel against the addend, so the result depends on the unrounded product
  fma: () => {
    const x = (uniform() - 0.5) * 2 ** (uniform() * 40 - 20), y = (uniform() - 0.5) * 2 ** (uniform() * 40 - 20);
    return [x, y, uniform() < 0.5 ? -(x * y) : (uniform() - 0.5) * 2 ** (uniform() * 80 - 40)];
  },
};
const arity = { pow: 2, atan2: 2, hypot: 2, fma: 3 };
for (const [name, gen] of Object.entries(domain)) {
  const n = samples + edges.length ** ((arity[name] ?? 1) > 1 ? 2 : 1);
  const xyz = [new Float64Array(n), new Float64Array(n), new Float64Array(n)];
  let i = 0;
  for (const a of edges) for (const b of (arity[name] ?? 1) > 1 ? edges : [1]) xyz[0][i] = a, xyz[1][i] = b, xyz[2][i++] = a;
  for (; i < n; i++) {
    const v = uniform() < 0.25 ? [randomBits(), randomBits(), randomBits()] : gen();
    for (let k = 0; k < 3; k++) xyz[k][i] = v[k] ?? 1;
  }
  const input = path.join(work, `${name}.in.bin`), output = path.join(work, `${name}.out.bin`);
  const header = Buffer.alloc(4);
  header.writeUInt32LE(n);
  fs.writeFileSync(input, Buffer.concat([header, ...xyz.map((a) => Buffer.from(a.buffer))]));
  execFileSync(path.resolve(args.native), ["evaluate", name, input, output], { stdio: "inherit" });
  const nativeBits = new BigUint64Array(new Uint8Array(fs.readFileSync(output)).buffer);
  const wasmBits = new BigUint64Array(yfoil.evaluate(name, xyz[0], xyz[1], xyz[2]).buffer);
  let differ = 0, nanOnly = 0;
  const examples = [];
  const hex = (b) => b.toString(16).padStart(16, "0");
  const inputBits = xyz.map((a) => new BigUint64Array(a.buffer));
  for (let j = 0; j < n; j++) {
    if (nativeBits[j] !== wasmBits[j]) {
      differ++;
      if (Number.isNaN(fromBits(nativeBits[j])) && Number.isNaN(fromBits(wasmBits[j]))) nanOnly++;
      if (examples.length < 3) {
        const args = inputBits.slice(0, arity[name] ?? 1).map((a) => hex(a[j])).join(",");
        examples.push(`(${args}) native ${hex(nativeBits[j])} wasm ${hex(wasmBits[j])}`);
      }
    }
  }
  report(`${name} on ${n} arguments`, differ === 0,
    differ ? `${differ} differ (${nanOnly} of them NaN against NaN), e.g. ${examples.join("; ")}` : "");
}

// 4. speed: the CI reference case's analysis, repeated (native timed by the example binary)
const options = path.join(work, "bench_options.json");
fs.writeFileSync(options, JSON.stringify({ alpha: 2 }));
const repeats = quick ? 5 : 50;
const nativeMs = Number(execFileSync(path.resolve(args.native), ["bench", REF_CASE, options, String(repeats)]).toString());
const geometry = fs.readFileSync(REF_CASE, "utf8");
yfoil.analyse(geometry, '{"alpha":2}');
const [, total] = timed(() => { for (let i = 0; i < repeats; i++) yfoil.analyse(geometry, '{"alpha":2}'); });
console.log(`\nspeed  reference case analysis: native ${nativeMs.toFixed(2)} ms, wasm ${(total / repeats).toFixed(2)} ms ` +
  `(x${(total / repeats / nativeMs).toFixed(2)}); all CLI runs: native ${(timing.native / 1e3).toFixed(1)} s ` +
  `(process start included), wasm ${(timing.wasm / 1e3).toFixed(1)} s`);
console.log(failures ? `\n${failures} comparison(s) DIFFERENT` : "\nevery comparison identical");
process.exit(failures ? 1 : 0);
