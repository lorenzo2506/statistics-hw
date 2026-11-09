"use strict";

/* ====================================================================
   Exercise 1: Riemann vs Lebesgue numerical integration
   Exercise 2: numerical computation of the standard normal CDF
   ==================================================================== */

const PI = Math.PI;
const now = () => (typeof performance !== "undefined" ? performance.now() : Date.now());

// ====================== Exercise 1 ======================

// Si(pi) = integral of sin(x)/x on [0, pi], by its (rapidly convergent) power series
function siPi() {
  let term = PI, sum = 0; // term_k = pi^(2k+1)/(2k+1)!
  for (let k = 0; k < 40; k++) {
    sum += (k % 2 ? -1 : 1) * term / (2 * k + 1);
    term *= PI * PI / ((2 * k + 2) * (2 * k + 3));
  }
  return sum;
}
const sinc = (x) => (x === 0 ? 1 : Math.sin(x) / x);

// muGE(t) = Lebesgue measure of the preimage {x : f(x) >= t}, known analytically for each f
const FUNS = [
  { name: "sin x on [0, π]", a: 0, b: PI, f: Math.sin, exact: 2, lo: 0, hi: 1,
    muGE: (t) => (t <= 0 ? PI : t > 1 ? 0 : PI - 2 * Math.asin(t)) },
  { name: "x² on [0, 1]", a: 0, b: 1, f: (x) => x * x, exact: 1 / 3, lo: 0, hi: 1,
    muGE: (t) => (t <= 0 ? 1 : t > 1 ? 0 : 1 - Math.sqrt(t)) },
  { name: "1 if x<0.5, 3 if x≥0.5", a: 0, b: 1, f: (x) => (x < 0.5 ? 1 : 3), exact: 2, lo: 1, hi: 3,
    muGE: (t) => (t <= 1 ? 1 : t <= 3 ? 0.5 : 0) },
  { name: "|x − 0.5| on [0, 1]", a: 0, b: 1, f: (x) => Math.abs(x - 0.5), exact: 0.25, lo: 0, hi: 0.5,
    muGE: (t) => (t <= 0 ? 1 : t > 0.5 ? 0 : 1 - 2 * t) },
  { name: "1/√x on (0, 1] (unbounded)", a: 0, b: 1, f: (x) => 1 / Math.sqrt(x), exact: 2, lo: 1, hi: 1e6, geom: true,
    muGE: (t) => (t <= 1 ? 1 : 1 / (t * t)) },
  { name: "1 if x<0.5, 0 if x≥0.5", a: 0, b: 1, f: (x) => (x < 0.5 ? 1 : 0), exact: 0.5, lo: 0, hi: 1,
    muGE: (t) => (t <= 0 ? 1 : t <= 1 ? 0.5 : 0) },
  { name: "sin x / x on [0, π]", a: 0, b: PI, f: sinc, exact: siPi(), lo: 0, hi: 1,
    muGE: (t) => { // sinc is decreasing: {sinc >= t} = [0, x_t], found by bisection
      if (t <= 0) return PI; if (t > 1) return 0;
      let l = 0, r = PI;
      for (let i = 0; i < 60; i++) { const m = (l + r) / 2; if (sinc(m) >= t) l = m; else r = m; }
      return l; } }
];

// Riemann sum with N equal subintervals: left or midpoint tags
function riemann(fn, N, mode) {
  const h = (fn.b - fn.a) / N;
  let s = 0;
  for (let i = 0; i < N; i++) s += fn.f(fn.a + (mode === "mid" ? i + 0.5 : i) * h);
  return s * h;
}

// partition of the RANGE [lo, hi] into K levels (uniform, or geometric for unbounded functions)
function levels(lo, hi, K, geom) {
  const y = [];
  for (let j = 0; j <= K; j++) y.push(geom ? lo * Math.pow(hi / lo, j / K) : lo + (hi - lo) * j / K);
  y[K] = hi;
  return y;
}

// Lebesgue sum of the simple function  s = sum_j y_j * 1_{E_j},  E_j = {y_j <= f < y_(j+1)}  (s <= f)
function lebesgueSum(muGE, lo, hi, K, geom) {
  const y = levels(lo, hi, K, geom);
  let s = 0;
  for (let j = 0; j < K; j++) s += y[j] * (muGE(y[j]) - muGE(y[j + 1])); // y_j * mu(E_j)
  return s + hi * muGE(hi);                                                // the set {f >= hi}
}

// Lebesgue sum with the measure of the preimages estimated from M samples of f (sort + binary search)
function lebesgueSampled(fn, M, K) {
  const h = (fn.b - fn.a) / M, v = new Float64Array(M);
  for (let i = 0; i < M; i++) v[i] = fn.f(fn.a + (i + 0.5) * h);
  v.sort();
  const muGE = (t) => {
    let l = 0, r = M;
    while (l < r) { const m = (l + r) >> 1; if (v[m] < t) l = m + 1; else r = m; }
    return (M - l) * h;
  };
  return lebesgueSum(muGE, v[0], v[M - 1], K, fn.geom);
}

function timed(f) { const t = now(); const v = f(); return { v, ms: now() - t }; }

// all four estimates for one function with N subintervals / K levels / M samples
function estimates(fn, N, withSampled) {
  const r = {};
  r.mid = timed(() => riemann(fn, N, "mid"));
  r.left = timed(() => riemann(fn, N, "left"));
  r.leb = timed(() => lebesgueSum(fn.muGE, fn.lo, fn.hi, N, fn.geom));
  r.smp = withSampled ? timed(() => lebesgueSampled(fn, N, N)) : null;
  return r;
}

// Thomae: smallest denominator q of a fraction p/q in [i/N, (i+1)/N]; sup of f there is 1/q
function thomaeUpper(N) {
  let U = 0;
  for (let i = 0; i < N; i++) {
    let q = 1;
    while (Math.floor(((i + 1) * q) / N) < Math.floor((i * q + N - 1) / N)) q++;
    U += 1 / q;
  }
  return U / N;
}
// number of fractions p/q in [0,1] with q <= k (Farey sequence length = 1 + sum of Euler's phi)
function fareyCount(k) {
  const phi = Array.from({ length: k + 1 }, (_, i) => i);
  for (let p = 2; p <= k; p++) if (phi[p] === p) for (let m = p; m <= k; m += p) phi[m] -= phi[m] / p;
  let c = 1;
  for (let q = 1; q <= k; q++) c += phi[q];
  return c;
}

// ====================== Exercise 2 ======================
const SQRT2 = Math.SQRT2, SQRTPI = Math.sqrt(PI), INV_SQRT_2PI = 1 / Math.sqrt(2 * PI);
const pdf = (t) => INV_SQRT_2PI * Math.exp(-0.5 * t * t);

// 1) numerical integration of the density on [0, |x|], then Phi = 1/2 +- integral
function phiTrapezoid(x, N) {
  const h = Math.abs(x) / N;
  let s = 0.5 * (pdf(0) + pdf(Math.abs(x)));
  for (let i = 1; i < N; i++) s += pdf(i * h);
  return x >= 0 ? 0.5 + s * h : 0.5 - s * h;
}
function phiSimpson(x, N) { // N even
  const h = Math.abs(x) / N;
  let s = pdf(0) + pdf(Math.abs(x));
  for (let i = 1; i < N; i++) s += (i % 2 ? 4 : 2) * pdf(i * h);
  s *= h / 3;
  return x >= 0 ? 0.5 + s : 0.5 - s;
}

// 2) power series (Hogg-McKean-Craig): alternating terms, cancellation for large |x|
function phiSeries(x) {
  let term = x, sum = 0, k = 0, big = 0; // term = (-1)^k x^(2k+1) / (2^k k!)
  for (; k < 600; k++) {
    const a = term / (2 * k + 1);
    sum += a; big = Math.max(big, Math.abs(a));
    if (Math.abs(a) < 1e-18 && k > Math.abs(x)) break;
    term *= -x * x / (2 * (k + 1));
  }
  return { value: 0.5 + INV_SQRT_2PI * sum, terms: k + 1, big: INV_SQRT_2PI * big };
}
// same integral with all-positive terms:  Phi(x) = 1/2 + pdf(x) * sum x^(2k+1) / (1*3*5*...*(2k+1))
function phiStableSeries(x) {
  let term = x, sum = 0, k = 0;
  for (; k < 2000; k++) {
    sum += term;
    if (Math.abs(term) < 1e-17 * Math.abs(sum)) break;
    term *= x * x / (2 * k + 3);
  }
  return { value: 0.5 + pdf(x) * sum, terms: k + 1 };
}

// 3) through the error function:  Phi(x) = (1 + erf(x/sqrt 2)) / 2
// 3a) Abramowitz-Stegun 7.1.26: fast rational approximation, |error| < 1.5e-7
function erfAS(z) {
  const t = 1 / (1 + 0.3275911 * Math.abs(z));
  const y = 1 - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t * Math.exp(-z * z);
  return z >= 0 ? y : -y;
}
const phiAS = (x) => 0.5 * (1 + erfAS(x / SQRT2));
// 3b) accurate: series for small arguments, continued fraction for the tail (erfc)
function erfSeries(z) { // z >= 0, positive terms
  let term = z, sum = z;
  for (let n = 1; n < 300; n++) { term *= 2 * z * z / (2 * n + 1); sum += term; if (term < 1e-17 * sum) break; }
  return 2 / SQRTPI * Math.exp(-z * z) * sum;
}
function erfcFraction(z) { // z > 0:  erfc z = e^(-z²)/sqrt(pi) / (z + (1/2)/(z + (2/2)/(z + (3/2)/(z + ...))))
  let t = z;
  for (let k = 80; k >= 1; k--) t = z + (k / 2) / t;
  return Math.exp(-z * z) / (SQRTPI * t);
}
function phiErf(x) {
  const z = Math.abs(x) / SQRT2;
  const tail = z < 3 ? 0.5 * (1 - erfSeries(z)) : 0.5 * erfcFraction(z); // P(Z > |x|)
  return x >= 0 ? 1 - tail : tail;
}

// trusted references found on the web (see the page for sources)
const REFS = [
  { x: 0, v: 0.5, src: "symmetry" },
  { x: 0.5, v: 0.691462, src: "Univ. Ulm table (6 decimals)" },
  { x: 1, v: (1 + 0.682689492137) / 2, src: "Wikipedia, 68–95–99.7 rule (erf, 12 digits)" },
  { x: 1.5, v: 0.933193, src: "Univ. Ulm table (6 decimals)" },
  { x: 1.96, v: 0.9750021048517796, src: "cdflib documentation (docs.rs)" },
  { x: 2, v: (1 + 0.954499736104) / 2, src: "Wikipedia, 68–95–99.7 rule" },
  { x: 3, v: (1 + 0.997300203937) / 2, src: "Wikipedia, 68–95–99.7 rule" },
  { x: -3, v: (1 - 0.997300203937) / 2, src: "Wikipedia, 68–95–99.7 rule (lower tail)" },
  { x: 4, v: (1 + 0.999936657516) / 2, src: "Wikipedia, 68–95–99.7 rule" },
  { x: -5, v: 2.866516e-7, src: "cdflib documentation: upper tail at 5 (7 digits)" }
];

// ====================== formatting ======================
function sci(x) { return x === 0 ? "0" : !isFinite(x) ? "∞" : x.toExponential(2); }
function fix(x, d) { return isFinite(x) ? x.toFixed(d) : "∞"; }
function ms(x) { return x < 0.01 ? "<0.01" : x.toFixed(x < 10 ? 2 : 0); }
function tableHtml(head, rows) {
  return "<tr>" + head.map((h) => `<th>${h}</th>`).join("") + "</tr>" +
    rows.map((r) => "<tr>" + r.map((c) => `<td>${c}</td>`).join("") + "</tr>").join("");
}

// ====================== user interface ======================
function initUI() {
  const $ = (id) => document.getElementById(id);
  const err = (v, ex) => (isFinite(v) ? sci(Math.abs(v - ex)) : "∞ / undefined");

  // ----- Exercise 1: overview for all functions -----
  function overview() {
    const N = Math.max(10, Math.min(200000, parseInt($("n").value, 10) || 1000));
    $("tOver").innerHTML = tableHtml(
      ["function", "exact value", `Riemann midpoint (N=${N})`, `Riemann left (N=${N})`, `Lebesgue, exact preimages (K=${N})`, `Lebesgue, sampled measure (M=K=${N})`],
      FUNS.map((fn) => { const e = estimates(fn, N, true);
        return [fn.name, fix(fn.exact, 8), err(e.mid.v, fn.exact), err(e.left.v, fn.exact), err(e.leb.v, fn.exact), err(e.smp.v, fn.exact)]; }));
    convergence();
  }

  // ----- convergence of one function -----
  function convergence() {
    const fn = FUNS[parseInt($("fsel").value, 10)];
    const Ns = [10, 100, 1000, 1e4, 1e5, 1e6], rows = [], E = [];
    for (const N of Ns) {
      const e = estimates(fn, N, N <= 1e5);
      const errs = [e.mid, e.left, e.leb, e.smp].map((m) => (m ? Math.abs(m.v - fn.exact) : NaN));
      E.push(errs);
      rows.push([N.toExponential(0).replace("e+", "e"), sci(errs[0]) + ` <span class="mut">${ms(e.mid.ms)} ms</span>`, isFinite(errs[1]) ? sci(errs[1]) : "∞",
        sci(errs[2]) + ` <span class="mut">${ms(e.leb.ms)} ms</span>`, e.smp ? sci(errs[3]) + ` <span class="mut">${ms(e.smp.ms)} ms</span>` : "–"]);
    }
    const order = (c, a, b) => { const p = E[a][c], q = E[b][c]; return p > 1e-14 && q > 1e-14 ? fix(Math.log10(p / q), 2) : "–"; };
    rows.push(["observed order p (error ∝ N<sup>−p</sup>, N=10³→10⁴)", order(0, 2, 3), order(1, 2, 3), order(2, 2, 3), order(3, 2, 3)]);
    $("tConv").innerHTML = tableHtml(["N", "Riemann midpoint: error, time", "Riemann left: error", "Lebesgue (exact preimages): error, time", "Lebesgue (sampled): error, time"], rows);
  }

  // ----- pathological functions -----
  function pathological() {
    const Ns = [10, 50, 100, 500, 1000, 3000];
    $("tThomae").innerHTML = tableHtml(["N subintervals", "lower Darboux sum", "upper Darboux sum"],
      Ns.map((N) => [N, "0", fix(thomaeUpper(N), 5)]));
    $("tFarey").innerHTML = tableHtml(["k", "number of points in {f ≥ 1/k}", "Lebesgue measure of {f ≥ 1/k}"],
      [1, 5, 10, 50, 100, 1000].map((k) => [k, fareyCount(k), "0"]));
    $("tDir").innerHTML = tableHtml(["ε", "total length of the cover of the first 40 rationals", "bound ε"],
      [1, 0.1, 0.001, 1e-6].map((e) => { let s = 0; for (let n = 1; n <= 40; n++) s += e / Math.pow(2, n); return [e, sci(s), e]; }));
  }

  // ----- Exercise 2 -----
  function normal() {
    const methods = [
      ["Trapezoid N=1000", (x) => phiTrapezoid(x, 1000)],
      ["Simpson N=1000", (x) => phiSimpson(x, 1000)],
      ["Alternating series", (x) => phiSeries(x).value],
      ["Positive-term series", (x) => phiStableSeries(x).value],
      ["erf, A&S 7.1.26", phiAS],
      ["erf, series + fraction", phiErf]
    ];
    const cell = (v, r) => { const e = Math.abs(v - r.v); return sci(e) + (r.v < 1e-3 ? ` <span class="mut">(rel ${sci(e / r.v)})</span>` : ""); };
    $("tRef").innerHTML = tableHtml(["x", "reference Φ(x)", "source"].concat(methods.map((m) => m[0] + ": abs. error")),
      REFS.map((r) => [r.x, r.v.toPrecision(12), r.src].concat(methods.map((m) => cell(m[1](r.x), r))))) ;

    const x0 = 1.96, ref = 0.9750021048517796, rows = []; let pt = NaN, ps = NaN;
    for (const N of [2, 4, 8, 16, 32, 64, 128, 256, 512]) {
      const et = Math.abs(phiTrapezoid(x0, N) - ref), es = Math.abs(phiSimpson(x0, N) - ref);
      rows.push([N, sci(et), fix(pt / et, 1), sci(es), fix(ps / es, 1)]); pt = et; ps = es;
    }
    $("tInt").innerHTML = tableHtml(["N", "trapezoid error", "ratio to previous (≈4 if error ∝ h²)", "Simpson error", "ratio to previous (≈16 if error ∝ h⁴)"], rows);

    $("tSer").innerHTML = tableHtml(["x", "alternating series: terms, largest term", "alternating error", "positive series: terms", "positive error"],
      [1, 2, 3, 4, 5, 6, 8].map((x) => { const a = phiSeries(x), b = phiStableSeries(x), r = phiErf(x);
        return [x, `${a.terms}, ${a.big.toPrecision(2)}`, sci(Math.abs(a.value - r)), b.terms, sci(Math.abs(b.value - r))]; }));

    const pts = Array.from({ length: 4001 }, (_, i) => -6 + 12 * i / 4000), tr = [];
    for (const [name, f] of methods) {
      const t = timed(() => { let s = 0; for (let rep = 0; rep < 5; rep++) for (const x of pts) s += f(x); return s; });
      let mx = 0; for (const x of pts) mx = Math.max(mx, Math.abs(f(x) - phiErf(x)));
      tr.push([name, sci(mx), (t.ms * 1e6 / (5 * pts.length)).toFixed(0) + " ns"]);
    }
    $("tTime").innerHTML = tableHtml(["method", "max abs. error on [−6, 6] (vs erf series + fraction)", "time per evaluation"], tr);
  }
  $("xgo").onclick = () => {
    const x = parseFloat($("xin").value); if (isNaN(x)) return;
    $("tCalc").innerHTML = tableHtml(["method", "Φ(x)"], [["Trapezoid N=1000", phiTrapezoid(x, 1000)], ["Simpson N=1000", phiSimpson(x, 1000)],
      ["Alternating series", phiSeries(x).value], ["Positive-term series", phiStableSeries(x).value], ["erf, A&S 7.1.26", phiAS(x)], ["erf, series + fraction", phiErf(x)]].map((r) => [r[0], r[1].toPrecision(15)]));
  };

  $("fsel").innerHTML = FUNS.map((f, i) => `<option value="${i}">${f.name}</option>`).join("");
  $("fsel").onchange = convergence;
  $("run").onclick = () => { $("status").textContent = "Running…"; setTimeout(() => { overview(); $("status").textContent = "Done."; }, 20); };
  overview(); pathological(); normal(); $("xgo").onclick();
}

if (typeof document !== "undefined") initUI();