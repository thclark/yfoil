// Smoke test of packaged release assets (scripts/wasm.sh package): the checksums hold, and each
// package, unpacked, loads the way its users load it and analyses the CI reference case.
//
//   node smoke.mjs DIST_DIR [--expect RECORD.json]
//
// DIST_DIR holds yfoil-wasm-<version>-{web,nodejs}.tar.gz and SHA256SUMS. --expect is the record
// `yfoil analyse <REF_CASE panels.json> -a 2 -o RECORD.json` writes, which both packages must
// reproduce byte for byte. Exits 1 on any failure.

import { createRequire } from "node:module";
import { execFileSync } from "node:child_process";
import crypto from "node:crypto";
import fs from "node:fs";
import os from "node:os";
import path from "node:path";
import { pathToFileURL } from "node:url";

const dist = path.resolve(process.argv[2] ?? "");
const expectAt = process.argv.indexOf("--expect");
const expected = expectAt > 0 ? fs.readFileSync(process.argv[expectAt + 1], "utf8") : null;
const repoRoot = new URL("../../../", import.meta.url);
const geometry = fs.readFileSync(new URL("tests/fixtures/xfoil/naca0012_n60_a2_re1e6/panels.json", repoRoot), "utf8");
let failures = 0;
const check = (label, ok, detail = "") => {
  if (!ok) failures++;
  console.log(`${ok ? "ok  " : "FAIL"}  ${label}${detail ? "  " + detail : ""}`);
};

// the checksums
const sums = fs.readFileSync(path.join(dist, "SHA256SUMS"), "utf8").trim().split("\n").map((l) => l.split(/\s+\*?/));
check("SHA256SUMS names two archives", sums.length === 2, sums.map((s) => s[1]).join(", "));
for (const [sum, file] of sums) {
  const actual = crypto.createHash("sha256").update(fs.readFileSync(path.join(dist, file))).digest("hex");
  check(`sha256 ${file}`, actual === sum);
}

const work = fs.mkdtempSync(path.join(os.tmpdir(), "yfoil-wasm-smoke-"));
const unpack = (target) => {
  const archive = sums.map((s) => s[1]).find((f) => f.endsWith(`-${target}.tar.gz`));
  execFileSync("tar", ["-xzf", path.join(dist, archive), "-C", work]);
  const dir = path.join(work, archive.replace(/\.tar\.gz$/, ""));
  for (const f of ["yfoil_wasm.js", "yfoil_wasm_bg.wasm", "yfoil_wasm.d.ts", "LICENSE", "LICENSE-core-math"]) {
    check(`${target}: ${f} present`, fs.existsSync(path.join(dir, f)));
  }
  return dir;
};
const exercise = (target, yfoil) => {
  const record = yfoil.analyse(geometry, JSON.stringify({ alpha: 2, name: "panels" }));
  const parsed = JSON.parse(record);
  check(`${target}: analyse reference case`, parsed.results.status === "Valid",
    `CL ${parsed.results.values?.cl}, version ${yfoil.version()}`);
  if (expected !== null) check(`${target}: record identical to native`, record === expected);
  const section = JSON.parse(yfoil.naca("4412", "", NaN));
  check(`${target}: naca 4412`, section.x.length === 160, `${section.x.length} nodes`);
};

// nodejs: require(), as Node users load it
const nodeDir = unpack("nodejs");
exercise("nodejs", createRequire(import.meta.url)(path.join(nodeDir, "yfoil_wasm.js")));

// web: an ES module initialised with the module's bytes, as a page does with fetch()
const webDir = unpack("web");
const web = await import(pathToFileURL(path.join(webDir, "yfoil_wasm.js")).href);
web.initSync({ module: fs.readFileSync(path.join(webDir, "yfoil_wasm_bg.wasm")) });
exercise("web", web);

fs.rmSync(work, { recursive: true, force: true });
console.log(failures ? `\n${failures} check(s) FAILED` : "\nevery check passed");
process.exit(failures ? 1 : 0);
