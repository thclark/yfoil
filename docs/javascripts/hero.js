/*
 * yFoil landing-page hero: streamlines advecting past an aerofoil.
 *
 * The flow is not decorative noise. Particles are advected through a real
 * incompressible potential-flow field with the Kutta condition applied at the
 * trailing edge, obtained from the Joukowski conformal map of flow past a
 * circular cylinder. It is the same flow XFOIL's inviscid panel solution
 * converges to, computed exactly rather than by panels.
 *
 * The velocity field sits behind a deliberately narrow interface:
 *
 *     field.outline()          -> [[x, y], ...]   the section, in world units
 *     field.velocityAt(x, y)   -> [u, v]          physical velocity, or null inside
 *     field.setAlpha(radians)
 *
 * Nothing below the field boundary knows the flow is analytic. When yFoil is
 * available as WASM, `joukowskiField` is replaced by a field backed by a real
 * yFoil solution — panel geometry, the converged inviscid velocities, and in
 * time the viscous edge velocity — and the renderer is untouched.
 */

(() => {
  "use strict";

  // -------------------------------------------------------------- constants --

  const CHORD_FRACTION = 0.46;  // of the canvas width the section spans
  const CHORD_FRACTION_TALL = 0.82;  // ... on a portrait canvas, where width is scarce
  const CHORD_CEILING = 0.9;    // ... but never more than this fraction of the height
  const FOIL_X = 0.60;          // section centre, as a fraction of canvas width
  const FOIL_Y = 0.42;          // section centre, as a fraction of canvas height
  const PARTICLE_AREA = 1250;   // one particle per this many CSS pixels
  const PARTICLE_RANGE = [420, 1900];
  const TRAIL = 18;             // positions retained per particle
  const SPEED = 1.35;           // world units per second at freestream
  const MAX_SPEED_DRAWN = 2.2;  // colour/alpha saturate here (TE is singular)
  const ALPHA = 6.0 * Math.PI / 180;   // fixed incidence

  // --------------------------------------------------------- complex helpers --
  //
  // Complex numbers are [re, im] pairs rather than objects: this runs a few
  // hundred thousand times a second and allocation is the whole cost.

  const cMul = (a, b) => [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]];

  const cDiv = (a, b) => {
    const d = b[0] * b[0] + b[1] * b[1];
    return [(a[0] * b[0] + a[1] * b[1]) / d, (a[1] * b[0] - a[0] * b[1]) / d];
  };

  const cSqrt = (a) => {
    const m = Math.hypot(a[0], a[1]);
    const re = Math.sqrt(Math.max(0, (m + a[0]) / 2));
    const im = Math.sqrt(Math.max(0, (m - a[0]) / 2));
    return [re, a[1] < 0 ? -im : im];
  };

  // ------------------------------------------------------------- flow field --

  /**
   * Flow past a Joukowski aerofoil.
   *
   * The map z = ζ + a²/ζ carries the exterior of a circle of radius R centred
   * at ζ0 = (−thickness, camber) onto the exterior of an aerofoil. The circle
   * is required to pass through ζ = a, which places a sharp trailing edge at
   * z = 2a. Circulation is set by the Kutta condition, so the rear stagnation
   * point sits exactly on that trailing edge — the same condition XFOIL
   * enforces on its panel solution.
   */
  function joukowskiField({ a = 1, thickness = 0.11, camber = 0.055 } = {}) {
    const z0 = [-thickness, camber];
    const R = Math.hypot(a + thickness, camber);
    const beta = Math.asin(Math.min(1, camber / R));

    let alpha = 0;
    let cosA = 1;
    let sinA = 0;
    let gamma = 0;

    function setAlpha(value) {
      alpha = value;
      cosA = Math.cos(alpha);
      sinA = Math.sin(alpha);
      // Kutta condition: Γ = 4πUR sin(α + β)
      gamma = 4 * Math.PI * R * Math.sin(alpha + beta);
    }

    /** Inverse map: the exterior root of ζ² − zζ + a² = 0. */
    function toCircle(z) {
      const zz = cMul(z, z);
      const root = cSqrt([zz[0] - 4 * a * a, zz[1]]);
      const p = [(z[0] + root[0]) / 2, (z[1] + root[1]) / 2];
      const q = [(z[0] - root[0]) / 2, (z[1] - root[1]) / 2];
      const dp = Math.hypot(p[0] - z0[0], p[1] - z0[1]);
      const dq = Math.hypot(q[0] - z0[0], q[1] - z0[1]);
      return dp >= dq ? p : q;
    }

    function velocityAt(x, y) {
      const zeta = toCircle([x, y]);
      const d = [zeta[0] - z0[0], zeta[1] - z0[1]];
      if (Math.hypot(d[0], d[1]) < R * 1.0005) return null;  // inside the section

      // dw/dζ = U e^{-iα} − U R² e^{iα} / (ζ−ζ0)² + iΓ / (2π(ζ−ζ0))
      const dd = cMul(d, d);
      const term2 = cDiv([R * R * cosA, R * R * sinA], dd);
      const term3 = cDiv([0, gamma / (2 * Math.PI)], d);
      const dwdz = [
        cosA - term2[0] + term3[0],
        -sinA - term2[1] + term3[1],
      ];

      // dz/dζ = 1 − a²/ζ², vanishing at the trailing edge where dw/dζ also
      // vanishes by the Kutta condition. The ratio is finite but numerically
      // delicate, so it is clamped rather than special-cased.
      const zz = cMul(zeta, zeta);
      const dzdz = [1 - cDiv([a * a, 0], zz)[0], -cDiv([a * a, 0], zz)[1]];
      if (Math.hypot(dzdz[0], dzdz[1]) < 1e-4) return [1, 0];

      const w = cDiv(dwdz, dzdz);
      const u = w[0];
      const v = -w[1];  // velocity is the conjugate of the complex potential
      const speed = Math.hypot(u, v);
      if (!Number.isFinite(speed)) return [1, 0];
      if (speed > 6) return [(u / speed) * 6, (v / speed) * 6];
      return [u, v];
    }

    function outline(samples = 260) {
      const points = [];
      for (let i = 0; i <= samples; i++) {
        const t = (i / samples) * 2 * Math.PI;
        const zeta = [z0[0] + R * Math.cos(t), z0[1] + R * Math.sin(t)];
        const inv = cDiv([a * a, 0], zeta);
        points.push([zeta[0] + inv[0], zeta[1] + inv[1]]);
      }
      return points;
    }

    setAlpha(0);
    return { setAlpha, velocityAt, outline, get alpha() { return alpha; } };
  }

  // ---------------------------------------------------------------- colours --

  function readPalette(element) {
    const style = getComputedStyle(element);
    const pick = (name, fallback) => (style.getPropertyValue(name).trim() || fallback);
    return {
      slow: pick("--yf-hero-glow-a", "#4f46e5"),
      fast: pick("--yf-hero-glow-b", "#06b6d4"),
      ink: pick("--md-default-fg-color", "#1a1a1a"),
      bg: pick("--yf-hero-ground", pick("--md-default-bg-color", "#ffffff")),
      streakBase: parseFloat(pick("--yf-hero-streak-base", "0.13")),
      streakGain: parseFloat(pick("--yf-hero-streak-gain", "0.50")),
    };
  }

  /** Pre-mix the speed ramp once per theme, so the draw loop only indexes it. */
  function buildRamp(from, to, steps = 24) {
    const probe = document.createElement("canvas").getContext("2d");
    const parse = (colour) => {
      probe.fillStyle = colour;
      const hex = probe.fillStyle;
      if (hex.startsWith("#")) {
        const n = parseInt(hex.slice(1), 16);
        return [(n >> 16) & 255, (n >> 8) & 255, n & 255];
      }
      return (hex.match(/\d+/g) || [80, 80, 200]).slice(0, 3).map(Number);
    };
    const [r0, g0, b0] = parse(from);
    const [r1, g1, b1] = parse(to);
    const ramp = [];
    for (let i = 0; i < steps; i++) {
      const t = i / (steps - 1);
      ramp.push([
        Math.round(r0 + (r1 - r0) * t),
        Math.round(g0 + (g1 - g0) * t),
        Math.round(b0 + (b1 - b0) * t),
      ]);
    }
    return ramp;
  }

  // ----------------------------------------------------------------- runner --

  function start(canvas) {
    const context = canvas.getContext("2d", { alpha: true });
    const field = joukowskiField();
    const outline = field.outline();
    const reduced = window.matchMedia("(prefers-reduced-motion: reduce)").matches;

    let palette = readPalette(canvas);
    let ramp = buildRamp(palette.slow, palette.fast);
    let width = 0;
    let height = 0;
    let scale = 1;
    let originX = 0;
    let originY = 0;
    let view = { xMin: -3, xMax: 3, yMin: -1, yMax: 1 };

    // The camera is rotated by −α, so the freestream runs horizontally across the
    // canvas and the section is seen to pitch. That is a change of viewpoint and
    // nothing else: the flow field is untouched, and a rotation of a solution of
    // the potential-flow equations is a solution of the same equations.
    let camCos = 1;
    let camSin = 0;

    const viewX = (x, y) => x * camCos + y * camSin;
    const viewY = (x, y) => -x * camSin + y * camCos;
    const worldX = (vx, vy) => vx * camCos - vy * camSin;
    const worldY = (vx, vy) => vx * camSin + vy * camCos;

    const toScreenX = (x, y) => originX + viewX(x, y) * scale;
    const toScreenY = (x, y) => originY - viewY(x, y) * scale;

    function respawn(index, anywhere) {
      // Seeded in view coordinates, then carried back to world coordinates, so
      // the inflow edge stays flush with the canvas at every incidence.
      const vx = anywhere
        ? view.xMin + Math.random() * (view.xMax - view.xMin)
        : view.xMin - Math.random() * 0.25;
      const vy = view.yMin + Math.random() * (view.yMax - view.yMin);
      px[index] = worldX(vx, vy);
      py[index] = worldY(vx, vy);
      age[index] = 0;
      life[index] = 3.5 + Math.random() * 4.5;
      filled[index] = 0;
    }

    // Particle state, held in flat arrays: position, age, and a ring buffer of
    // the last TRAIL world positions used to draw the streak. The population is
    // sized to the canvas area, so a full-viewport hero on a large display is
    // seeded as densely as a small one, and reallocated when the canvas resizes.
    let count = 0;
    let px = new Float32Array(0);
    let py = new Float32Array(0);
    let age = new Float32Array(0);
    let life = new Float32Array(0);
    let trail = new Float32Array(0);
    let filled = new Uint8Array(0);

    function allocate(n) {
      count = n;
      px = new Float32Array(n);
      py = new Float32Array(n);
      age = new Float32Array(n);
      life = new Float32Array(n);
      trail = new Float32Array(n * TRAIL * 2);
      filled = new Uint8Array(n);
      for (let i = 0; i < n; i++) respawn(i, true);
    }

    function resize() {
      const rect = canvas.getBoundingClientRect();
      if (rect.width < 1 || rect.height < 1) return false;
      const dpr = Math.min(window.devicePixelRatio || 1, 2);
      width = rect.width;
      height = rect.height;
      canvas.width = Math.round(width * dpr);
      canvas.height = Math.round(height * dpr);
      context.setTransform(dpr, 0, 0, dpr, 0, 0);

      // The section spans 4a in world units. It is fitted to a fraction of the
      // canvas width, capped against the height so that a short, wide viewport
      // does not push it off the top and bottom of the panel.
      const chordFraction = width < height ? CHORD_FRACTION_TALL : CHORD_FRACTION;
      scale = Math.min(width * chordFraction, height * CHORD_CEILING) / 4;
      originX = width * FOIL_X;
      originY = height * FOIL_Y;
      view = {
        xMin: -originX / scale,
        xMax: (width - originX) / scale,
        yMin: -(height - originY) / scale,
        yMax: originY / scale,
      };

      const target = Math.round(
        Math.min(
          PARTICLE_RANGE[1],
          Math.max(PARTICLE_RANGE[0], (width * height) / PARTICLE_AREA),
        ),
      );
      if (target !== count) allocate(target);
      return true;
    }

    function step(dt) {
      for (let i = 0; i < count; i++) {
        age[i] += dt;

        // Midpoint (RK2) integration of dz/dt = V(z).
        const v1 = field.velocityAt(px[i], py[i]);
        if (!v1) { respawn(i, false); continue; }
        const hx = px[i] + v1[0] * SPEED * dt * 0.5;
        const hy = py[i] + v1[1] * SPEED * dt * 0.5;
        const v2 = field.velocityAt(hx, hy) || v1;
        const nx = px[i] + v2[0] * SPEED * dt;
        const ny = py[i] + v2[1] * SPEED * dt;

        const vx = viewX(nx, ny);
        const vy = viewY(nx, ny);
        const out =
          vx > view.xMax + 0.3 || vx < view.xMin - 0.6 ||
          vy > view.yMax + 0.4 || vy < view.yMin - 0.4;
        if (out || age[i] > life[i]) { respawn(i, false); continue; }

        px[i] = nx;
        py[i] = ny;

        // Shift the trail ring buffer and record the new head.
        const base = i * TRAIL * 2;
        for (let k = TRAIL - 1; k > 0; k--) {
          trail[base + k * 2] = trail[base + (k - 1) * 2];
          trail[base + k * 2 + 1] = trail[base + (k - 1) * 2 + 1];
        }
        trail[base] = nx;
        trail[base + 1] = ny;
        if (filled[i] < TRAIL) filled[i]++;
      }
    }

    function draw() {
      context.clearRect(0, 0, width, height);
      context.lineCap = "round";
      context.lineJoin = "round";

      for (let i = 0; i < count; i++) {
        const n = filled[i];
        if (n < 3) continue;
        const base = i * TRAIL * 2;

        const v = field.velocityAt(px[i], py[i]);
        if (!v) continue;
        const speed = Math.hypot(v[0], v[1]);
        const t = Math.min(1, speed / MAX_SPEED_DRAWN);
        const [r, g, b] = ramp[Math.min(ramp.length - 1, (t * ramp.length) | 0)];

        // Fade in on birth and out on death, so respawns are not a flicker.
        const fade = Math.min(1, age[i] / 0.6) * Math.min(1, (life[i] - age[i]) / 0.9);
        const alpha = palette.streakBase + palette.streakGain * t * t;

        context.strokeStyle = `rgba(${r}, ${g}, ${b}, ${(alpha * fade).toFixed(3)})`;
        context.lineWidth = 0.7 + 1.5 * t;
        context.beginPath();
        context.moveTo(toScreenX(trail[base], trail[base + 1]),
                       toScreenY(trail[base], trail[base + 1]));
        for (let k = 1; k < n; k++) {
          const tx = trail[base + k * 2];
          const ty = trail[base + k * 2 + 1];
          context.lineTo(toScreenX(tx, ty), toScreenY(tx, ty));
        }
        context.stroke();
      }

      // The section itself, drawn last so streaks never cross it.
      context.beginPath();
      context.moveTo(toScreenX(outline[0][0], outline[0][1]),
                     toScreenY(outline[0][0], outline[0][1]));
      for (let i = 1; i < outline.length; i++) {
        context.lineTo(toScreenX(outline[i][0], outline[i][1]),
                       toScreenY(outline[i][0], outline[i][1]));
      }
      context.closePath();
      context.fillStyle = palette.bg;
      context.fill();

      const gradient = context.createLinearGradient(
        toScreenX(-2, 0), toScreenY(-2, 0), toScreenX(2, 0), toScreenY(2, 0),
      );
      gradient.addColorStop(0, palette.fast);
      gradient.addColorStop(1, palette.slow);
      context.strokeStyle = gradient;
      context.lineWidth = 1.6;
      context.stroke();
    }

    // ------------------------------------------------------------ animation --

    let raf = 0;
    let last = 0;
    let visible = true;

    field.setAlpha(ALPHA);
    camCos = Math.cos(ALPHA);
    camSin = Math.sin(ALPHA);

    function frame(now) {
      raf = requestAnimationFrame(frame);
      if (!visible) { last = now; return; }

      const dt = Math.min(0.05, last ? (now - last) / 1000 : 0.016);
      last = now;

      step(dt);
      draw();
    }

    if (!resize()) return () => {};

    if (reduced) {
      // No animation loop: integrate a still frame of streaklines and stop.
      for (let i = 0; i < 260; i++) step(0.016);
      draw();
    } else {
      raf = requestAnimationFrame(frame);
    }

    const resizeObserver = new ResizeObserver(() => {
      if (resize() && reduced) draw();
    });
    resizeObserver.observe(canvas);

    const intersectionObserver = new IntersectionObserver(
      ([entry]) => { visible = entry.isIntersecting; },
      { threshold: 0 },
    );
    intersectionObserver.observe(canvas);

    const onVisibility = () => { if (document.hidden) visible = false; };
    const themeObserver = new MutationObserver(() => {
      palette = readPalette(canvas);
      ramp = buildRamp(palette.slow, palette.fast);
      if (reduced) draw();
    });
    themeObserver.observe(document.body, {
      attributes: true,
      attributeFilter: ["data-md-color-scheme"],
    });

    document.addEventListener("visibilitychange", onVisibility);

    return () => {
      cancelAnimationFrame(raf);
      resizeObserver.disconnect();
      intersectionObserver.disconnect();
      themeObserver.disconnect();
      document.removeEventListener("visibilitychange", onVisibility);
    };
  }

  // ------------------------------------------------------------------ setup --
  //
  // Instant navigation swaps the document body without a page load, so the
  // previous canvas has to be torn down explicitly on every navigation.

  let teardown = null;

  function mount() {
    if (teardown) { teardown(); teardown = null; }
    const host = document.querySelector(".mdx-hero");
    if (!host) return;
    let canvas = host.querySelector("canvas");
    if (!canvas) {
      canvas = document.createElement("canvas");
      canvas.setAttribute("aria-hidden", "true");
      host.appendChild(canvas);
    }
    teardown = start(canvas);
  }

  if (typeof document$ !== "undefined") {
    document$.subscribe(mount);
  } else if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", mount);
  } else {
    mount();
  }
})();
