"use strict";

/* ====================================================================
   Welford's recurrences and sampling mean / variance of hackers' scores
   ==================================================================== */

// ---------- random numbers (seeded, reproducible) ----------
function mulberry32(seed) {
  let a = seed >>> 0;
  return function () {
    a = (a + 0x6D2B79F5) >>> 0;
    let t = a;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

// ====================== Exercise 1: Welford ======================

// Welford on an array. The core logic is the 3 lines inside the loop.
function welford(xs) {
  let mean = 0, ss = 0, k = 0;
  const rows = [];
  for (const x of xs) {
    k++;
    const delta = x - mean;        // δ_k = x_k - M_{k-1}
    mean += delta / k;             // M_k = M_{k-1} + δ_k / k
    ss += delta * (x - mean);      // SS_k = SS_{k-1} + δ_k (x_k - M_k)
    rows.push({ k, x, delta, mean, ss, s2: k > 1 ? ss / (k - 1) : NaN });
  }
  return { mean, ss, n: k, variance: k > 1 ? ss / (k - 1) : NaN, popVariance: ss / k, rows };
}

// two-pass reference: first the mean, then the sum of squared deviations
function twoPass(xs) {
  const mean = xs.reduce((a, b) => a + b, 0) / xs.length;
  const ss = xs.reduce((a, b) => a + (b - mean) * (b - mean), 0);
  return { mean, ss };
}

// textbook one-pass formula: SS = sum(x^2) - (sum x)^2 / n  (numerically unstable)
function naiveOnePass(xs) {
  let s = 0, s2 = 0;
  for (const x of xs) { s += x; s2 += x * x; }
  return { mean: s / xs.length, ss: s2 - s * s / xs.length };
}

// ====================== Exercise 2: experiment ======================

// pmf over scores 0..10 and its theoretical moments
function makePmf(kind, p, weights) {
  let w = new Array(11).fill(0);
  if (kind === "uniform") w = w.map(() => 1);
  else if (kind === "custom") w = weights.slice(0, 11).map((v) => Math.max(0, v || 0));
  else { // Binomial(10, p)
    let c = 1;
    for (let j = 0; j <= 10; j++) {
      w[j] = c * Math.pow(p, j) * Math.pow(1 - p, 10 - j);
      c = c * (10 - j) / (j + 1);
    }
  }
  const total = w.reduce((a, b) => a + b, 0) || 1;
  const q = w.map((v) => v / total);
  let mu = 0;
  q.forEach((v, j) => { mu += j * v; });
  let s2 = 0, m4 = 0;
  q.forEach((v, j) => { s2 += (j - mu) ** 2 * v; m4 += (j - mu) ** 4 * v; });
  const cum = [];
  q.reduce((acc, v, j) => { cum[j] = acc + v; return cum[j]; }, 0);
  cum[10] = 1;
  return { q, cum, mu, sigma2: s2, mu4: m4 };
}

function sampler(pmf, rng) {
  const cum = pmf.cum;
  return function () {
    const u = rng();
    let j = 0;
    while (u >= cum[j]) j++;
    return j;
  };
}

// one run: m hackers, each with n scores; returns the sampling mean and variances of each hacker
function runExperiment(pmf, m, n, rng) {
  const draw = sampler(pmf, rng);
  const means = new Float64Array(m), s2 = new Float64Array(m), biased = new Float64Array(m);
  for (let h = 0; h < m; h++) {
    let mean = 0, ss = 0;
    for (let k = 1; k <= n; k++) {
      const x = draw();
      const delta = x - mean;
      mean += delta / k;
      ss += delta * (x - mean);
    }
    means[h] = mean;
    s2[h] = ss / (n - 1);
    biased[h] = ss / n;
  }
  return { means, s2, biased };
}

function avg(a) {
  let s = 0;
  for (const v of a) s += v;
  return s / a.length;
}

function varOf(a) {
  const mu = avg(a);
  let s = 0;
  for (const v of a) s += (v - mu) * (v - mu);
  return s / (a.length - 1);
}

// theoretical variance of the sample variance s^2 (divisor n-1)
function theoryVarS2(pmf, n) {
  return pmf.mu4 / n - pmf.sigma2 * pmf.sigma2 * (n - 3) / (n * (n - 1));
}

function summarise(pmf, n, res) {
  return {
    n,
    meanOfMeans: [avg(res.means), pmf.mu],
    varOfMeans: [varOf(res.means), pmf.sigma2 / n],
    meanOfVars: [avg(res.s2), pmf.sigma2],
    meanOfBiased: [avg(res.biased), pmf.sigma2 * (n - 1) / n],
    varOfVars: [varOf(res.s2), theoryVarS2(pmf, n)]
  };
}

// running Welford estimates for one hacker, stored for every k
function runningPath(pmf, nmax, rng) {
  const draw = sampler(pmf, rng);
  const mean = new Float64Array(nmax + 1), vari = new Float64Array(nmax + 1);
  let mu = 0, ss = 0;
  for (let k = 1; k <= nmax; k++) {
    const x = draw();
    const delta = x - mu;
    mu += delta / k;
    ss += delta * (x - mu);
    mean[k] = mu;
    vari[k] = k > 1 ? ss / (k - 1) : 0;
  }
  return { mean, vari };
}

// ====================== plotting helpers (SVG strings) ======================
function fmt(x) {
  if (x === 0) return "0";
  if (!isFinite(x)) return "–";
  const a = Math.abs(x);
  if (a >= 1e5 || a < 1e-3) return x.toExponential(2);
  return String(Number(x.toPrecision(4)));
}

function tableHtml(head, rows) {
  return "<tr>" + head.map((h) => `<th>${h}</th>`).join("") + "</tr>" +
    rows.map((r) => "<tr>" + r.map((c) => `<td>${c}</td>`).join("") + "</tr>").join("");
}

function histogramSvg(values, lo, hi, bins, centre, variance, caption) {
  const W = 260, H = 140, L = 8, R = 8, T = 14, B = 22;
  const bw = (hi - lo) / bins, counts = new Array(bins).fill(0);
  for (const v of values) {
    const b = Math.floor((v - lo) / bw);
    if (b >= 0 && b < bins) counts[b]++;
  }
  const dens = counts.map((c) => c / (values.length * bw));
  const pdf = (x) => Math.exp(-((x - centre) ** 2) / (2 * variance)) / Math.sqrt(2 * Math.PI * variance);
  const yMax = 1.1 * Math.max(...dens, pdf(centre));
  const X = (x) => L + (x - lo) / (hi - lo) * (W - L - R);
  const Y = (y) => H - B - y / yMax * (H - B - T);
  let svg = `<svg viewBox="0 0 ${W} ${H}" role="img" aria-label="${caption}"><text x="${L}" y="10">${caption}</text>`;
  dens.forEach((d, i) => {
    const x0 = X(lo + i * bw), x1 = X(lo + (i + 1) * bw);
    svg += `<rect class="bar" x="${x0.toFixed(1)}" y="${Y(d).toFixed(1)}" width="${Math.max(0.6, x1 - x0 - 0.4).toFixed(1)}" height="${(H - B - Y(d)).toFixed(1)}"/>`;
  });
  let path = "";
  for (let i = 0; i <= 120; i++) {
    const x = lo + (hi - lo) * i / 120;
    path += (i ? "L" : "M") + X(x).toFixed(1) + " " + Y(pdf(x)).toFixed(1);
  }
  svg += `<path class="curve" d="${path}"/><line class="axis" x1="${L}" y1="${H - B}" x2="${W - R}" y2="${H - B}"/>`;
  [lo, (lo + hi) / 2, hi].forEach((v, i) => {
    svg += `<text x="${X(v)}" y="${H - 8}" text-anchor="${i === 0 ? "start" : i === 2 ? "end" : "middle"}">${fmt(v)}</text>`;
  });
  return svg + "</svg>";
}

function pathsSvg(paths, ks, spec) {
  const W = 300, H = 190, L = 38, R = 8, T = 16, B = 24;
  let ymin = Infinity, ymax = -Infinity;
  const upd = (y) => { ymin = Math.min(ymin, y); ymax = Math.max(ymax, y); };
  for (const path of paths) for (const k of ks) upd(path[k]);
  for (const k of ks) { upd(spec.lo(k)); upd(spec.hi(k)); }
  const pad = 0.04 * (ymax - ymin || 1);
  ymin -= pad; ymax += pad;
  const kmax = ks[ks.length - 1], lg = Math.max(1e-9, Math.log10(kmax));
  const X = (k) => L + Math.log10(k) / lg * (W - L - R);
  const Y = (y) => H - B - (y - ymin) / (ymax - ymin) * (H - B - T);
  let svg = `<svg viewBox="0 0 ${W} ${H}" role="img" aria-label="${spec.title}"><text x="${L}" y="11">${spec.title}</text>`;
  let poly = ks.map((k) => `${X(k).toFixed(1)},${Y(spec.lo(k)).toFixed(1)}`).join(" ");
  poly += " " + ks.slice().reverse().map((k) => `${X(k).toFixed(1)},${Y(spec.hi(k)).toFixed(1)}`).join(" ");
  svg += `<polygon class="band" points="${poly}"/>`;
  svg += `<polyline class="meanline" points="${ks.map((k) => `${X(k).toFixed(1)},${Y(spec.mean).toFixed(1)}`).join(" ")}"/>`;
  paths.forEach((path, i) => {
    svg += `<polyline fill="none" stroke="hsl(${i * 60},70%,50%)" stroke-width="1.1" points="${ks.map((k) => `${X(k).toFixed(1)},${Y(path[k]).toFixed(1)}`).join(" ")}"/>`;
  });
  svg += `<line class="axis" x1="${L}" y1="${H - B}" x2="${W - R}" y2="${H - B}"/><line class="axis" x1="${L}" y1="${T}" x2="${L}" y2="${H - B}"/>`;
  for (let j = 0; j <= Math.floor(lg + 1e-9); j++) {
    svg += `<text x="${X(Math.pow(10, j))}" y="${H - 10}" text-anchor="middle">${Math.pow(10, j)}</text>`;
  }
  svg += `<text x="${W - R}" y="${H - 1}" text-anchor="end">k (log scale)</text>`;
  [ymin + pad, ymax - pad].forEach((v) => { svg += `<text x="${L - 3}" y="${Y(v) + 3}" text-anchor="end">${fmt(v)}</text>`; });
  return svg + "</svg>";
}

function logSpaced(kmin, kmax, count) {
  const set = new Set([kmin, kmax]);
  const a = Math.log10(kmin), b = Math.log10(kmax);
  for (let i = 0; i <= count; i++) set.add(Math.max(kmin, Math.round(Math.pow(10, a + (b - a) * i / count))));
  return [...set].sort((x, y) => x - y);
}

// ====================== user interface ======================
function initUI() {
  const $ = (id) => document.getElementById(id);
  const clamp = (v, a, b) => Math.min(b, Math.max(a, v));
  const pair = (a) => `${fmt(a[0])} <span class="mut">(${fmt(a[1])})</span>`;

  // ----- Exercise 1: step by step -----
  function showSteps() {
    const xs = $("data").value.split(",").map((s) => parseFloat(s.trim())).filter((v) => !isNaN(v));
    if (xs.length < 2) { $("tSteps").innerHTML = ""; $("nSteps").textContent = "Enter at least two numbers."; return; }
    const w = welford(xs), ref = twoPass(xs);
    $("tSteps").innerHTML = tableHtml(["k", "x<sub>k</sub>", "δ<sub>k</sub> = x<sub>k</sub> − M<sub>k−1</sub>", "M<sub>k</sub>", "SS<sub>k</sub>", "s²<sub>k</sub> = SS<sub>k</sub>/(k−1)"],
      w.rows.map((r) => [r.k, fmt(r.x), fmt(r.delta), fmt(r.mean), fmt(r.ss), isNaN(r.s2) ? "–" : fmt(r.s2)]));
    $("nSteps").textContent = `Two-pass check: mean = ${fmt(ref.mean)}, SS = ${fmt(ref.ss)}, variance = ${fmt(ref.ss / (xs.length - 1))}. Welford gives mean = ${fmt(w.mean)}, SS = ${fmt(w.ss)}, variance = ${fmt(w.variance)}.`;
  }
  function showStability() {
    const base = $("data").value.split(",").map((s) => parseFloat(s.trim())).filter((v) => !isNaN(v));
    const off = parseFloat($("offset").value) || 0;
    if (base.length < 2) return;
    const xs = base.map((v) => v + off), truth = twoPass(base);
    const rows = [["True value (data without offset)", fmt(truth.ss), "–"],
      ["Textbook formula Σx² − (Σx)²/n", fmt(naiveOnePass(xs).ss), fmt(naiveOnePass(xs).ss - truth.ss)],
      ["Two-pass", fmt(twoPass(xs).ss), fmt(twoPass(xs).ss - truth.ss)],
      ["Welford", fmt(welford(xs).ss), fmt(welford(xs).ss - truth.ss)]];
    $("tStab").innerHTML = tableHtml(["method (data + offset)", "SS", "error"], rows);
  }
  $("steps").onclick = () => { showSteps(); showStability(); };
  $("stab").onclick = showStability;

  // ----- Exercise 2 -----
  function render() {
    const kind = $("dist").value;
    const p = clamp(parseFloat($("p").value) || 0.3, 0.01, 0.99);
    const weights = $("weights").value.split(",").map((s) => parseFloat(s.trim()));
    let m = clamp(parseInt($("m").value, 10) || 2000, 100, 5000);
    const seed = parseInt($("seed").value, 10) || 1;
    let ns = [...new Set($("ns").value.split(",").map((s) => parseInt(s.trim(), 10)).filter((v) => v >= 2 && v <= 10000))].sort((a, b) => a - b).slice(0, 6);
    if (!ns.length) ns = [5, 10, 50, 100, 1000];
    const total = ns.reduce((a, b) => a + b, 0);
    let note = "";
    if (m * total > 1e8) { m = Math.floor(1e8 / total); note = ` (m reduced to ${m} to keep the simulation fast)`; }
    const pmf = makePmf(kind, p, weights);
    const rng = mulberry32(seed);

    $("theory").innerHTML = `Theoretical distribution of one score: μ = <b>${fmt(pmf.mu)}</b>, σ² = <b>${fmt(pmf.sigma2)}</b> (σ = ${fmt(Math.sqrt(pmf.sigma2))}), fourth central moment μ<sub>4</sub> = ${fmt(pmf.mu4)}. Probabilities of scores 0…10: ${pmf.q.map((v) => fmt(v)).join(", ")}.`;

    // running Welford estimates for 6 hackers
    const nmax = ns[ns.length - 1];
    const paths = Array.from({ length: 6 }, () => runningPath(pmf, nmax, rng));
    const ks = logSpaced(1, nmax, 160), ks2 = ks.filter((k) => k >= 2);
    const sdVar = (k) => Math.sqrt(Math.max(0, theoryVarS2(pmf, k)));
    $("runs").innerHTML =
      pathsSvg(paths.map((r) => r.mean), ks, { title: "running mean M_k", mean: pmf.mu, lo: (k) => pmf.mu - 2 * Math.sqrt(pmf.sigma2 / k), hi: (k) => pmf.mu + 2 * Math.sqrt(pmf.sigma2 / k) }) +
      pathsSvg(paths.map((r) => r.vari), ks2, { title: "running variance s²_k", mean: pmf.sigma2, lo: (k) => pmf.sigma2 - 2 * sdVar(k), hi: (k) => pmf.sigma2 + 2 * sdVar(k) });

    // experiment for each n (fresh scores for every n)
    const results = ns.map((n) => ({ n, res: runExperiment(pmf, m, n, rng) }));
    const sums = results.map((r) => summarise(pmf, r.n, r.res));
    const shared = $("shared").checked;
    const sdMean = (n) => Math.sqrt(pmf.sigma2 / n);
    const rangeMean = (list) => [Math.min(...list.map((n) => pmf.mu - 4 * sdMean(n))), Math.max(...list.map((n) => pmf.mu + 4 * sdMean(n)))];
    const rangeVar = (list) => [Math.max(0, Math.min(...list.map((n) => pmf.sigma2 - 4 * Math.sqrt(theoryVarS2(pmf, n))))),
      Math.max(...list.map((n) => pmf.sigma2 + 4 * Math.sqrt(theoryVarS2(pmf, n))))];
    $("hMean").innerHTML = results.map((r) => {
      const [lo, hi] = shared ? rangeMean(ns) : rangeMean([r.n]);
      return histogramSvg(r.res.means, lo, hi, 41, pmf.mu, pmf.sigma2 / r.n, `n = ${r.n}`);
    }).join("");
    $("hVar").innerHTML = results.map((r) => {
      const [lo, hi] = shared ? rangeVar(ns) : rangeVar([r.n]);
      return histogramSvg(r.res.s2, lo, hi, 41, pmf.sigma2, theoryVarS2(pmf, r.n), `n = ${r.n}`);
    }).join("");

    $("tRes").innerHTML = tableHtml(
      ["n", "average of sampling means (μ)", "Var of sampling means (σ²/n)", "average of s², divisor n−1 (σ²)", "average of SS/n (σ²(n−1)/n)", "Var of s² (theory)"],
      sums.map((s) => [s.n, pair(s.meanOfMeans), pair(s.varOfMeans), pair(s.meanOfVars), pair(s.meanOfBiased), pair(s.varOfVars)]));
    $("nRes").textContent = `Empirical value (theory in grey) from m = ${m} hackers per n. With m hackers the average of the means has itself a standard error σ/√(nm), so small differences from μ are noise.`;
    $("status").textContent = `Done${note}. Seed ${seed}.`;
  }

  $("run").onclick = () => { $("status").textContent = "Simulating…"; setTimeout(render, 20); };
  $("shared").onchange = render;
  showSteps();
  showStability();
  render();
}

if (typeof document !== "undefined") initUI();