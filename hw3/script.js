"use strict";

/* ====================================================================
   Penetration process: LLN and CLT simulation
   Bernoulli and Rademacher increments, S_n, S_n/n and S_n/sqrt(n)
   ==================================================================== */

// ---------- random numbers (seeded, so results are reproducible) ----------
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

// ---------- models ----------
// mu = E[X], sigma2 = Var[X], step = spacing of the values of S_n
function modelInfo(model, p) {
  if (model === "bernoulli") {
    return { model, p, mu: p, sigma2: p * (1 - p), step: 1, low: 0, high: 1 };
  }
  return { model, p: 0.5, mu: 0, sigma2: 1, step: 2, low: -1, high: 1 };
}

// ---------- simulation ----------
// Simulates m attackers over nmax systems. Returns S_n for every checkpoint n in ns
// and the full trajectories of the first nPaths attackers.
function simulate(info, m, ns, seed, nPaths) {
  const rng = mulberry32(seed);
  const nmax = ns[ns.length - 1];
  const bernoulli = info.model === "bernoulli";
  const values = ns.map(() => new Float64Array(m));
  const paths = [];
  for (let j = 0; j < m; j++) {
    const keep = j < nPaths ? new Float64Array(nmax + 1) : null;
    let s = 0;
    let c = 0;
    for (let i = 1; i <= nmax; i++) {
      if (bernoulli) s += rng() < info.p ? 1 : 0;
      else s += rng() < 0.5 ? 1 : -1;
      if (keep) keep[i] = s;
      if (i === ns[c]) { values[c][j] = s; c++; }
    }
    if (keep) paths.push(keep);
  }
  return { values, paths, nmax };
}

// ---------- statistics ----------
function mean(a) {
  let s = 0;
  for (const v of a) s += v;
  return s / a.length;
}

function variance(a) { // unbiased sample variance
  const mu = mean(a);
  let s = 0;
  for (const v of a) s += (v - mu) * (v - mu);
  return s / (a.length - 1);
}

function normalCdf(x) { // Abramowitz-Stegun 7.1.26, error < 1.5e-7
  const z = Math.abs(x) / Math.SQRT2;
  const t = 1 / (1 + 0.3275911 * z);
  const poly = ((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t + 0.254829592;
  const erf = 1 - poly * t * Math.exp(-z * z);
  return 0.5 * (1 + (x >= 0 ? erf : -erf));
}

function ksDistance(values, cdf) { // Kolmogorov-Smirnov distance to a continuous cdf
  const sorted = Float64Array.from(values).sort();
  const m = sorted.length;
  let d = 0;
  for (let i = 0; i < m; i++) {
    const f = cdf(sorted[i]);
    d = Math.max(d, Math.abs((i + 1) / m - f), Math.abs(i / m - f));
  }
  return d;
}

function skewKurtosis(z) { // sample skewness and excess kurtosis
  const mu = mean(z);
  let m2 = 0, m3 = 0, m4 = 0;
  for (const v of z) {
    const d = v - mu;
    m2 += d * d; m3 += d * d * d; m4 += d * d * d * d;
  }
  const n = z.length;
  m2 /= n; m3 /= n; m4 /= n;
  return { skew: m3 / Math.pow(m2, 1.5), kurt: m4 / (m2 * m2) - 3 };
}

// theory for the standardized sum Z_n = (S_n - n mu)/(sigma sqrt n)
function theorySkewKurt(info, n) {
  const pq = info.sigma2;
  if (info.model === "bernoulli") {
    return { skew: (1 - 2 * info.p) / Math.sqrt(n * pq), kurt: (1 - 6 * pq) / (n * pq) };
  }
  return { skew: 0, kurt: -2 / n };
}

// ---------- plot helpers (SVG strings) ----------
function fmt(x) {
  if (x === 0) return "0";
  const a = Math.abs(x);
  if (a >= 1e5 || a < 1e-3) return x.toExponential(1);
  return String(Number(x.toPrecision(3)));
}

// Quantities: S, A (average), N (normalized). scale(n) divides S_n.
const QUANTITIES = {
  S: { scale: () => 1, name: "S_n" },
  A: { scale: (n) => n, name: "S_n/n" },
  N: { scale: (n) => Math.sqrt(n), name: "S_n/sqrt(n)" }
};

// theoretical mean and variance of each quantity
function quantityMoments(info, key, n) {
  const sc = QUANTITIES[key].scale(n);
  return { mean: n * info.mu / sc, variance: n * info.sigma2 / (sc * sc) };
}

// horizontal range for a histogram (4 standard deviations around the mean)
function quantityRange(info, key, nList) {
  let lo = Infinity, hi = -Infinity;
  for (const n of nList) {
    const q = quantityMoments(info, key, n);
    const sd = Math.sqrt(q.variance);
    lo = Math.min(lo, q.mean - 4 * sd);
    hi = Math.max(hi, q.mean + 4 * sd);
  }
  if (info.model === "bernoulli" && key === "A") { lo = Math.max(lo, 0); hi = Math.min(hi, 1); }
  if (info.model === "bernoulli" && key === "S") lo = Math.max(lo, 0);
  if (info.model === "rademacher" && key === "A") { lo = Math.max(lo, -1); hi = Math.min(hi, 1); }
  return [lo, hi];
}

// Histogram of S-values seen through quantity q = S/scale. Bins are centred on the lattice of S.
function histogramSvg(sValues, n, key, info, range, bins, caption) {
  const W = 260, H = 140, L = 8, R = 8, T = 14, B = 22;
  const scale = QUANTITIES[key].scale(n);
  const stepS = info.step;
  const anchor = info.model === "rademacher" ? n % 2 : 0;
  const loS = range[0] * scale, hiS = range[1] * scale;
  let c = Math.ceil((hiS - loS) / bins / stepS);
  if (c % 2 === 0) c++;
  c = Math.max(1, c);
  const wS = stepS * c;
  const base = anchor - stepS / 2;
  const e0 = base + wS * Math.floor((loS - base) / wS);
  const nb = Math.max(1, Math.ceil((hiS - e0) / wS));
  const counts = new Array(nb).fill(0);
  for (const s of sValues) {
    const b = Math.floor((s - e0) / wS);
    if (b >= 0 && b < nb) counts[b]++;
  }
  const m = sValues.length;
  const widthQ = wS / scale;
  const dens = counts.map((k) => k / (m * widthQ));
  const q = quantityMoments(info, key, n);
  const pdf = (x) => Math.exp(-((x - q.mean) ** 2) / (2 * q.variance)) / Math.sqrt(2 * Math.PI * q.variance);
  const xmin = e0 / scale, xmax = (e0 + nb * wS) / scale;
  const yMax = 1.1 * Math.max(...dens, pdf(q.mean));
  const X = (x) => L + (x - xmin) / (xmax - xmin) * (W - L - R);
  const Y = (y) => H - B - y / yMax * (H - B - T);
  let svg = `<svg viewBox="0 0 ${W} ${H}" role="img" aria-label="${caption}">`;
  svg += `<text x="${L}" y="10">${caption}</text>`;
  dens.forEach((d, i) => {
    const x0 = X((e0 + i * wS) / scale), x1 = X((e0 + (i + 1) * wS) / scale);
    svg += `<rect class="bar" x="${x0.toFixed(1)}" y="${Y(d).toFixed(1)}" width="${Math.max(0.6, x1 - x0 - 0.4).toFixed(1)}" height="${(H - B - Y(d)).toFixed(1)}"/>`;
  });
  let path = "";
  for (let i = 0; i <= 120; i++) {
    const x = xmin + (xmax - xmin) * i / 120;
    path += (i ? "L" : "M") + X(x).toFixed(1) + " " + Y(pdf(x)).toFixed(1);
  }
  svg += `<path class="curve" d="${path}"/>`;
  svg += `<line class="axis" x1="${L}" y1="${H - B}" x2="${W - R}" y2="${H - B}"/>`;
  [xmin, (xmin + xmax) / 2, xmax].forEach((v, i) => {
    svg += `<text x="${X(v)}" y="${H - 8}" text-anchor="${i === 0 ? "start" : i === 2 ? "end" : "middle"}">${fmt(v)}</text>`;
  });
  return svg + "</svg>";
}

// Trajectory plot with logarithmic horizontal axis.
function pathsSvg(paths, ks, spec) {
  const W = 300, H = 190, L = 38, R = 8, T = 16, B = 24;
  let ymin = Infinity, ymax = -Infinity;
  const upd = (y) => { ymin = Math.min(ymin, y); ymax = Math.max(ymax, y); };
  for (const path of paths) for (const k of ks) upd(spec.f(k, path[k]));
  for (const k of ks) { upd(spec.lo(k)); upd(spec.hi(k)); }
  const pad = 0.04 * (ymax - ymin || 1);
  ymin -= pad; ymax += pad;
  const kmax = ks[ks.length - 1];
  const lg = Math.max(1e-9, Math.log10(kmax));
  const X = (k) => L + Math.log10(k) / lg * (W - L - R);
  const Y = (y) => H - B - (y - ymin) / (ymax - ymin) * (H - B - T);
  let svg = `<svg viewBox="0 0 ${W} ${H}" role="img" aria-label="${spec.title}"><text x="${L}" y="11">${spec.title}</text>`;
  let poly = ks.map((k) => `${X(k).toFixed(1)},${Y(spec.lo(k)).toFixed(1)}`).join(" ");
  poly += " " + ks.slice().reverse().map((k) => `${X(k).toFixed(1)},${Y(spec.hi(k)).toFixed(1)}`).join(" ");
  svg += `<polygon class="band" points="${poly}"/>`;
  svg += `<polyline class="meanline" points="${ks.map((k) => `${X(k).toFixed(1)},${Y(spec.mean(k)).toFixed(1)}`).join(" ")}"/>`;
  paths.forEach((path, i) => {
    const pts = ks.map((k) => `${X(k).toFixed(1)},${Y(spec.f(k, path[k])).toFixed(1)}`).join(" ");
    svg += `<polyline fill="none" stroke="hsl(${i * 45},70%,50%)" stroke-width="1.1" points="${pts}"/>`;
  });
  svg += `<line class="axis" x1="${L}" y1="${H - B}" x2="${W - R}" y2="${H - B}"/><line class="axis" x1="${L}" y1="${T}" x2="${L}" y2="${H - B}"/>`;
  for (let j = 0; j <= Math.floor(lg + 1e-9); j++) {
    svg += `<text x="${X(Math.pow(10, j))}" y="${H - 10}" text-anchor="middle">${Math.pow(10, j)}</text>`;
  }
  svg += `<text x="${W - R}" y="${H - 1}" text-anchor="end">k (log scale)</text>`;
  [ymin + pad, ymax - pad].forEach((v) => { svg += `<text x="${L - 3}" y="${Y(v) + 3}" text-anchor="end">${fmt(v)}</text>`; });
  return svg + "</svg>";
}

function logSpaced(kmax, count) {
  const set = new Set([1, kmax]);
  for (let i = 0; i <= count; i++) set.add(Math.max(1, Math.round(Math.pow(10, Math.log10(kmax) * i / count))));
  return [...set].sort((a, b) => a - b);
}

function tableHtml(head, rows) {
  return "<tr>" + head.map((h) => `<th>${h}</th>`).join("") + "</tr>" +
    rows.map((r) => "<tr>" + r.map((c) => `<td>${c}</td>`).join("") + "</tr>").join("");
}

// ---------- analysis of one simulation ----------
function analyse(info, ns, values, eps) {
  const m = values[0].length;
  const out = { mom: [], lln: [], clt: [] };
  ns.forEach((n, idx) => {
    const S = values[idx];
    const avg = S.map((s) => s / n);
    const nor = S.map((s) => s / Math.sqrt(n));
    const sd = Math.sqrt(info.sigma2);
    out.mom.push({
      n,
      ES: [mean(S), n * info.mu], VS: [variance(S), n * info.sigma2],
      EA: [mean(avg), info.mu], VA: [variance(avg), info.sigma2 / n],
      EN: [mean(nor), Math.sqrt(n) * info.mu], VN: [variance(nor), info.sigma2]
    });
    let exceed = 0, maxDev = 0;
    for (const a of avg) {
      const dev = Math.abs(a - info.mu);
      if (dev >= eps) exceed++;
      maxDev = Math.max(maxDev, dev);
    }
    const range = info.high - info.low;
    out.lln.push({
      n, p: exceed / m, maxDev,
      cheb: Math.min(1, info.sigma2 / (n * eps * eps)),
      hoef: Math.min(1, 2 * Math.exp(-2 * n * eps * eps / (range * range)))
    });
    const z = S.map((s) => (s - n * info.mu) / (sd * Math.sqrt(n)));
    const sk = skewKurtosis(z), th = theorySkewKurt(info, n);
    out.clt.push({ n, ks: ksDistance(z, normalCdf), skew: [sk.skew, th.skew], kurt: [sk.kurt, th.kurt] });
  });
  return out;
}

// ---------- user interface ----------
function initUI() {
  const $ = (id) => document.getElementById(id);
  const clamp = (v, a, b) => Math.min(b, Math.max(a, v));
  const pair = (a) => `${fmt(a[0])} <span class="mut">(${fmt(a[1])})</span>`;

  function readNs() {
    const list = $("ns").value.split(",").map((s) => parseInt(s.trim(), 10)).filter((v) => v >= 1 && v <= 100000);
    return [...new Set(list)].sort((a, b) => a - b).slice(0, 6);
  }

  function render() {
    const model = $("model").value;
    const p = clamp(parseFloat($("p").value) || 0.3, 0.01, 0.99);
    const m = clamp(parseInt($("m").value, 10) || 2000, 100, 5000);
    const seed = parseInt($("seed").value, 10) || 1;
    const eps = clamp(parseFloat($("eps").value) || 0.1, 0.001, 1);
    let ns = readNs();
    if (!ns.length) ns = [10, 100, 1000, 10000];
    const info = modelInfo(model, p);
    const sim = simulate(info, m, ns, seed, 8);
    const res = analyse(info, ns, sim.values, eps);

    // trajectories
    const ks = logSpaced(sim.nmax, 160);
    const sd = Math.sqrt(info.sigma2), mu = info.mu;
    const specs = [
      { title: "Total score S_k", f: (k, s) => s, mean: (k) => mu * k, lo: (k) => mu * k - 2 * sd * Math.sqrt(k), hi: (k) => mu * k + 2 * sd * Math.sqrt(k) },
      { title: "Average S_k/k", f: (k, s) => s / k, mean: () => mu, lo: (k) => mu - 2 * sd / Math.sqrt(k), hi: (k) => mu + 2 * sd / Math.sqrt(k) },
      { title: "Normalized S_k/sqrt(k)", f: (k, s) => s / Math.sqrt(k), mean: (k) => mu * Math.sqrt(k), lo: (k) => mu * Math.sqrt(k) - 2 * sd, hi: (k) => mu * Math.sqrt(k) + 2 * sd }
    ];
    $("paths").innerHTML = specs.map((sp) => pathsSvg(sim.paths, ks, sp)).join("");

    // histograms
    const shared = $("shared").checked;
    for (const [id, key] of [["hS", "S"], ["hA", "A"], ["hN", "N"]]) {
      const common = quantityRange(info, key, ns);
      $(id).innerHTML = ns.map((n, i) => {
        const range = shared ? common : quantityRange(info, key, [n]);
        return histogramSvg(sim.values[i], n, key, info, range, 41, `n = ${n}`);
      }).join("");
    }

    // tables
    $("tMom").innerHTML = tableHtml(
      ["n", "E[S<sub>n</sub>]", "Var[S<sub>n</sub>]", "E[S<sub>n</sub>/n]", "Var[S<sub>n</sub>/n]", "E[S<sub>n</sub>/√n]", "Var[S<sub>n</sub>/√n]"],
      res.mom.map((r) => [r.n, pair(r.ES), pair(r.VS), pair(r.EA), pair(r.VA), pair(r.EN), pair(r.VN)]));
    $("nMom").innerHTML = `Empirical value (theoretical value in grey), from m = ${m} attackers. Var[S<sub>n</sub>] = nσ² grows, Var[S<sub>n</sub>/n] = σ²/n → 0, Var[S<sub>n</sub>/√n] = σ² = ${fmt(info.sigma2)} stays constant.`;

    $("tLLN").innerHTML = tableHtml(
      ["n", `P(|S<sub>n</sub>/n − μ| ≥ ε), ε = ${eps}`, "Chebyshev bound σ²/(nε²)", "Hoeffding bound", "max deviation over the m attackers"],
      res.lln.map((r) => [r.n, fmt(r.p), fmt(r.cheb), fmt(r.hoef), fmt(r.maxDev)]));
    $("nLLN").textContent = `μ = ${fmt(mu)}. The empirical probability goes to 0 as n grows, and stays below both bounds (bounds are capped at 1). Chebyshev is loose; Hoeffding decays exponentially.`;

    $("tCLT").innerHTML = tableHtml(
      ["n", "KS distance of Z<sub>n</sub> to N(0,1)", "noise level 1.36/√m", "skewness (theory)", "excess kurtosis (theory)"],
      res.clt.map((r) => [r.n, fmt(r.ks), fmt(1.358 / Math.sqrt(m)), pair(r.skew), pair(r.kurt)]));
    $("nCLT").innerHTML = `Z<sub>n</sub> = (S<sub>n</sub> − nμ)/(σ√n). For the Rademacher model Z<sub>n</sub> = S<sub>n</sub>/√n exactly. ` +
      (model === "bernoulli"
        ? "For Bernoulli, S<sub>n</sub>/√n itself has mean √n·p, which diverges, so only the centred version Z<sub>n</sub> converges to N(0,1)."
        : "The KS distance of a lattice variable is at least about 0.4/√n, so it decreases until it reaches the sampling noise level.");
    $("status").textContent = `Done: m = ${m}, n = ${ns.join(", ")}, seed = ${seed}.`;
  }

  $("model").onchange = () => { $("p").disabled = $("model").value === "rademacher"; };
  $("run").onclick = () => {
    $("status").textContent = "Simulating…";
    setTimeout(render, 20);
  };
  $("shared").onchange = render;
  render();
}

if (typeof document !== "undefined") initUI();