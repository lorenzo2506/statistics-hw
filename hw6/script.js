"use strict";

/* ====================================================================
   Statistical independence: contingency table, expected frequencies,
   chi-square statistic, p-value (and a permutation p-value)
   ==================================================================== */

// ---------- random numbers (seeded) ----------
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

// ---------- example datasets (simulated) ----------
// groups: [label of X, size, probabilities of the categories of F]
function simulate(groups, fLabels, seed) {
  const rng = mulberry32(seed), lines = [];
  for (const [x, n, probs] of groups) {
    for (let i = 0; i < n; i++) {
      let u = rng(), j = 0;
      while (j < probs.length - 1 && u >= probs[j]) { u -= probs[j]; j++; }
      lines.push(`${x},${fLabels[j]}`);
    }
  }
  return lines.join("\n");
}
const EXAMPLES = {
  dependent: () => simulate([["No training", 120, [0.40, 0.60]], ["Basic training", 100, [0.25, 0.75]], ["Advanced training", 80, [0.12, 0.88]]],
    ["Clicked", "Did not click"], 11),
  independent: () => simulate([["Windows", 150, [0.5, 0.3, 0.2]], ["Linux", 90, [0.5, 0.3, 0.2]], ["macOS", 60, [0.5, 0.3, 0.2]]],
    ["Low", "Medium", "High"], 7),
  small: () => simulate([["Team A", 10, [0.6, 0.3, 0.1]], ["Team B", 10, [0.3, 0.4, 0.3]], ["Team C", 10, [0.2, 0.3, 0.5]]],
    ["Low", "Medium", "High"], 3)
};

// ---------- data handling ----------
function parsePairs(text, hasHeader) {
  const lines = text.split(/\r?\n/).map((s) => s.trim()).filter((s) => s.length);
  const pairs = [];
  lines.forEach((line, i) => {
    if (hasHeader && i === 0) return;
    const parts = line.split(/[,;\t]/).map((s) => s.trim());
    if (parts.length >= 2 && parts[0] !== "" && parts[1] !== "") pairs.push([parts[0], parts[1]]);
  });
  return pairs;
}

// contingency table: categories in order of first appearance
function crosstab(pairs) {
  const xs = [], fs = [], xi = [], fi = [];
  for (const [x, f] of pairs) {
    let a = xs.indexOf(x); if (a < 0) { xs.push(x); a = xs.length - 1; }
    let b = fs.indexOf(f); if (b < 0) { fs.push(f); b = fs.length - 1; }
    xi.push(a); fi.push(b);
  }
  const counts = xs.map(() => new Array(fs.length).fill(0));
  for (let k = 0; k < xi.length; k++) counts[xi[k]][fi[k]]++;
  return { xs, fs, xi, fi, counts };
}

// marginals, conditional frequencies, expected counts, chi-square
function analyse(tab) {
  const { counts } = tab, r = counts.length, c = counts[0].length;
  const rowT = counts.map((row) => row.reduce((a, b) => a + b, 0));
  const colT = counts[0].map((_, j) => counts.reduce((a, row) => a + row[j], 0));
  const N = rowT.reduce((a, b) => a + b, 0);
  const expected = counts.map((row, i) => row.map((_, j) => rowT[i] * colT[j] / N)); // e_ij = n_i. n_.j / N
  let chi2 = 0, minExp = Infinity, low = 0;
  const contrib = counts.map((row, i) => row.map((o, j) => {
    const e = expected[i][j], v = (o - e) * (o - e) / e;
    chi2 += v; minExp = Math.min(minExp, e); if (e < 5) low++;
    return v;
  }));
  const resid = counts.map((row, i) => row.map((o, j) => (o - expected[i][j]) / Math.sqrt(expected[i][j])));
  const df = (r - 1) * (c - 1);
  const cramerV = Math.sqrt(chi2 / (N * Math.min(r - 1, c - 1)));
  return { r, c, N, rowT, colT, expected, contrib, resid, chi2, df, cramerV, minExp, lowShare: low / (r * c) };
}

// ---------- chi-square distribution: p-value through the regularized incomplete gamma ----------
function lgamma(x) { // Lanczos approximation
  const g = [76.18009172947146, -86.50532032941677, 24.01409824083091, -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5];
  let y = x, t = x + 5.5;
  t -= (x + 0.5) * Math.log(t);
  let s = 1.000000000190015;
  for (const v of g) s += v / ++y;
  return -t + Math.log(2.5066282746310005 * s / x);
}
function gammaQ(a, x) { // upper regularized incomplete gamma Q(a, x)
  if (x <= 0) return 1;
  if (x < a + 1) { // series for P, then Q = 1 - P
    let ap = a, sum = 1 / a, del = sum;
    for (let n = 0; n < 1000; n++) { ap++; del *= x / ap; sum += del; if (Math.abs(del) < Math.abs(sum) * 1e-15) break; }
    return 1 - sum * Math.exp(-x + a * Math.log(x) - lgamma(a));
  }
  let b = x + 1 - a, c = 1e300, d = 1 / b, h = d; // continued fraction (Lentz)
  for (let i = 1; i < 1000; i++) {
    const an = -i * (i - a);
    b += 2; d = an * d + b; if (Math.abs(d) < 1e-300) d = 1e-300;
    c = b + an / c; if (Math.abs(c) < 1e-300) c = 1e-300;
    d = 1 / d; const del = d * c; h *= del;
    if (Math.abs(del - 1) < 1e-15) break;
  }
  return Math.exp(-x + a * Math.log(x) - lgamma(a)) * h;
}
const chi2PValue = (x, df) => gammaQ(df / 2, x / 2);
function chi2Critical(df, alpha) { // value with P(chi2 > value) = alpha, by bisection
  let lo = 0, hi = 1000;
  for (let i = 0; i < 200; i++) { const m = (lo + hi) / 2; if (chi2PValue(m, df) > alpha) lo = m; else hi = m; }
  return (lo + hi) / 2;
}

// ---------- permutation p-value: shuffle F, margins stay fixed ----------
function permutationP(tab, an, B, seed) {
  const rng = mulberry32(seed), f = tab.fi.slice(), x = tab.xi, n = f.length, c = an.c, r = an.r;
  let hits = 0;
  for (let b = 0; b < B; b++) {
    for (let i = n - 1; i > 0; i--) { const j = Math.floor(rng() * (i + 1)); const t = f[i]; f[i] = f[j]; f[j] = t; }
    const cnt = new Array(r * c).fill(0);
    for (let k = 0; k < n; k++) cnt[x[k] * c + f[k]]++;
    let s = 0;
    for (let i = 0; i < r; i++) for (let j = 0; j < c; j++) { const e = an.expected[i][j]; s += (cnt[i * c + j] - e) ** 2 / e; }
    if (s >= an.chi2 - 1e-9) hits++;
  }
  return (hits + 1) / (B + 1);
}

// ---------- formatting ----------
const pct = (v) => (100 * v).toFixed(1) + "%";
const f3 = (v) => v.toFixed(3);
function pText(p) { return p < 1e-4 ? p.toExponential(2) : p.toFixed(4); }
function table(head, rows) {
  return "<tr>" + head.map((h) => `<th>${h}</th>`).join("") + "</tr>" +
    rows.map((r) => "<tr>" + r.map((c, i) => `<${i ? "td" : "th"}>${c}</${i ? "td" : "th"}>`).join("") + "</tr>").join("");
}
const esc = (s) => String(s).replace(/&/g, "&amp;").replace(/</g, "&lt;");

// stacked horizontal bars: conditional distribution of F for each X, plus the marginal one
function barsSvg(tab, an) {
  const W = 520, rowH = 26, L = 130, rows = an.r + 1, H = rows * rowH + 40;
  const colors = tab.fs.map((_, j) => `hsl(${(j * 67 + 20) % 360},65%,55%)`);
  let svg = `<svg viewBox="0 0 ${W} ${H}" role="img" aria-label="Conditional distributions">`;
  const draw = (label, probs, y, bold) => {
    svg += `<text x="${L - 6}" y="${y + 15}" text-anchor="end"${bold ? ' style="font-weight:bold"' : ""}>${esc(label.slice(0, 22))}</text>`;
    let x = L;
    probs.forEach((p, j) => {
      const w = p * (W - L - 10);
      svg += `<rect x="${x.toFixed(1)}" y="${y + 3}" width="${w.toFixed(1)}" height="${rowH - 8}" fill="${colors[j]}"/>`;
      if (w > 28) svg += `<text x="${(x + w / 2).toFixed(1)}" y="${y + 16}" text-anchor="middle" style="fill:#000">${Math.round(p * 100)}%</text>`;
      x += w;
    });
  };
  tab.xs.forEach((x, i) => draw(`${x}`, tab.counts[i].map((o) => o / an.rowT[i]), i * rowH + 4, false));
  draw("all (marginal)", an.colT.map((v) => v / an.N), an.r * rowH + 14, true);
  let lx = L;
  tab.fs.forEach((f, j) => { svg += `<rect x="${lx}" y="${H - 18}" width="10" height="10" fill="${colors[j]}"/><text x="${lx + 14}" y="${H - 9}">${esc(f.slice(0, 18))}</text>`; lx += 24 + 6 * Math.min(18, f.length); });
  return svg + "</svg>";
}

// ====================== user interface ======================
function initUI() {
  const $ = (id) => document.getElementById(id);

  function run() {
    const pairs = parsePairs($("data").value, $("header").checked);
    const alpha = Math.min(0.5, Math.max(0.0001, parseFloat($("alpha").value) || 0.05));
    if (pairs.length < 4) { $("out").hidden = true; $("status").textContent = "Enter at least 4 rows 'X,F'."; return; }
    const tab = crosstab(pairs);
    if (tab.xs.length < 2 || tab.fs.length < 2) { $("out").hidden = true; $("status").textContent = "Both variables need at least 2 different categories."; return; }
    const an = analyse(tab), { r, c, N } = an;
    $("out").hidden = false;
    $("status").textContent = `${N} observations, ${r} × ${c} table.`;

    // observed counts with marginals
    $("tObs").innerHTML = table(["X \\ F"].concat(tab.fs.map(esc), ["total n<sub>i·</sub> (f<sub>i·</sub>)"]),
      tab.xs.map((x, i) => [esc(x)].concat(tab.counts[i], [`${an.rowT[i]} (${pct(an.rowT[i] / N)})`]))
        .concat([["total n<sub>·j</sub> (f<sub>·j</sub>)"].concat(an.colT.map((v) => `${v} (${pct(v / N)})`), [`${N} (100%)`])]));
    // joint relative frequencies
    $("tJoint").innerHTML = table(["X \\ F"].concat(tab.fs.map(esc)),
      tab.xs.map((x, i) => [esc(x)].concat(tab.counts[i].map((o) => f3(o / N)))));
    // conditional frequencies
    $("tRow").innerHTML = table(["X \\ F"].concat(tab.fs.map(esc), ["sum"]),
      tab.xs.map((x, i) => [esc(x)].concat(tab.counts[i].map((o) => pct(o / an.rowT[i])), ["100%"]))
        .concat([["marginal f<sub>·j</sub>"].concat(an.colT.map((v) => pct(v / N)), ["100%"])]));
    $("tCol").innerHTML = table(["X \\ F"].concat(tab.fs.map(esc), ["marginal f<sub>i·</sub>"]),
      tab.xs.map((x, i) => [esc(x)].concat(tab.counts[i].map((o, j) => pct(o / an.colT[j])), [pct(an.rowT[i] / N)]))
        .concat([["sum"].concat(tab.fs.map(() => "100%"), ["–"])]));
    $("bars").innerHTML = barsSvg(tab, an);
    // expected counts, contributions and residuals
    $("tExp").innerHTML = table(["X \\ F"].concat(tab.fs.map(esc)),
      tab.xs.map((x, i) => [esc(x)].concat(an.expected[i].map((e, j) => `${e.toFixed(2)} <span class="mut">(obs ${tab.counts[i][j]})</span>`))));
    $("tChi").innerHTML = table(["X \\ F"].concat(tab.fs.map(esc), ["row sum"]),
      tab.xs.map((x, i) => [esc(x)].concat(an.contrib[i].map((v, j) => `${f3(v)} <span class="mut">(res ${an.resid[i][j].toFixed(2)})</span>`), [f3(an.contrib[i].reduce((a, b) => a + b, 0))]))
        .concat([["χ²"].concat(tab.fs.map(() => ""), [`<b>${f3(an.chi2)}</b>`])]));

    // test
    const p = chi2PValue(an.chi2, an.df), crit = chi2Critical(an.df, alpha);
    const B = Math.max(200, Math.min(3000, Math.floor(2e7 / N)));
    const pPerm = permutationP(tab, an, B, 12345);
    const strength = an.cramerV < 0.1 ? "negligible" : an.cramerV < 0.3 ? "weak to moderate" : an.cramerV < 0.5 ? "moderate to strong" : "strong";
    let worst = [0, 0];
    for (let i = 0; i < r; i++) for (let j = 0; j < c; j++) if (Math.abs(an.resid[i][j]) > Math.abs(an.resid[worst[0]][worst[1]])) worst = [i, j];
    $("stats").innerHTML =
      `<p>χ² = Σ (O − E)²/E = <b>${f3(an.chi2)}</b>, degrees of freedom (r−1)(c−1) = <b>${an.df}</b>, critical value at α = ${alpha}: <b>${crit.toFixed(3)}</b>.<br>
       p-value (χ² distribution) = <b>${pText(p)}</b> &nbsp;|&nbsp; permutation p-value (${B} reshuffles) = ${pText(pPerm)}<br>
       Cramér’s V = ${f3(an.cramerV)} (${strength} association)</p>`;
    const reject = p < alpha;
    let txt = reject
      ? `<p><b>Decision:</b> p = ${pText(p)} &lt; α = ${alpha} (χ² = ${f3(an.chi2)} &gt; ${crit.toFixed(3)}): we <b>reject independence</b>. The data are not compatible with X and F being independent: the distribution of F changes with X (compare the conditional rows with the marginal row above). The strength of the association measured by Cramér’s V is ${strength} (${f3(an.cramerV)}).</p>`
      : `<p><b>Decision:</b> p = ${pText(p)} ≥ α = ${alpha} (χ² = ${f3(an.chi2)} ≤ ${crit.toFixed(3)}): we <b>do not reject independence</b>. The observed differences between conditional and marginal frequencies are of the size expected from random sampling. This does not prove independence: it only says that the data show no evidence of association (with a small sample the test has little power).</p>`;
    txt += `<p>Largest standardized residual: X = “${esc(tab.xs[worst[0]])}”, F = “${esc(tab.fs[worst[1]])}”: observed ${tab.counts[worst[0]][worst[1]]}, expected ${an.expected[worst[0]][worst[1]].toFixed(2)} (residual ${an.resid[worst[0]][worst[1]].toFixed(2)}); residuals beyond ±2 point at the cells that depart most from independence.</p>`;
    if (an.minExp < 5 && an.lowShare > 0.2 || an.minExp < 1) txt += `<p><b>Warning:</b> ${(100 * an.lowShare).toFixed(0)}% of the expected counts are below 5 (minimum ${an.minExp.toFixed(2)}), so the χ² approximation is unreliable. Trust the permutation p-value more, or merge categories / collect more data.</p>`;
    $("verdict").innerHTML = txt;
  }

  function load(name) { $("data").value = EXAMPLES[name](); $("header").checked = false; run(); }
  $("ex").onchange = () => load($("ex").value);
  $("go").onclick = run;
  load("dependent");
}

if (typeof document !== "undefined") initUI();