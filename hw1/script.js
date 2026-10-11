/* ==========================================================
   Toy elliptic curve: y^2 = x^3 + A*x + B (mod p), here A = 0, B = 7
   ========================================================== */

const A = 0;
const B = 7;

/* Small DOM helper */
const $ = (id) => document.getElementById(id);

/* ----------------------------------------------------------
   Modular arithmetic
   ---------------------------------------------------------- */

/** Always-positive modulo. */
const mod = (a, p) => ((a % p) + p) % p;

/** Modular inverse via extended Euclid. */
function inv(a, p) {
  let [r0, r1] = [mod(a, p), p];
  let [s0, s1] = [1, 0];

  while (r1) {
    const q = Math.floor(r0 / r1);
    [r0, r1] = [r1, r0 - q * r1];
    [s0, s1] = [s1, s0 - q * s1];
  }
  return mod(s0, p);
}

/* ----------------------------------------------------------
   Elliptic-curve group law
   null represents the point at infinity O.
   `info` (optional) is filled with the kind of operation and λ.
   ---------------------------------------------------------- */

function add(P, Q, p, info) {
  if (!P) return Q;
  if (!Q) return P;

  const [x1, y1] = P;
  const [x2, y2] = Q;

  // P + (-P) = O
  if (x1 === x2 && mod(y1 + y2, p) === 0) {
    if (info) info.kind = 'O';
    return null;
  }

  let lambda;
  if (x1 === x2 && y1 === y2) {
    // Point doubling
    lambda = mod(3 * x1 * x1 + A, p) * inv(2 * y1, p);
    if (info) info.kind = 'double';
  } else {
    // General addition
    lambda = mod(y2 - y1, p) * inv(x2 - x1, p);
    if (info) info.kind = 'add';
  }
  lambda = mod(lambda, p);

  const x3 = mod(lambda * lambda - x1 - x2, p);
  const y3 = mod(lambda * (x1 - x3) - y1, p);

  if (info) info.lam = lambda;
  return [x3, y3];
}

/** Double-and-add scalar multiplication; also returns a trace of the steps. */
function scalarMul(k, P, p) {
  let R = null;
  const trace = [];
  const bits = k.toString(2);

  for (const bit of bits) {
    R = add(R, R, p);
    const afterDouble = R;
    if (bit === '1') R = add(R, P, p);
    trace.push({ bit, afterDouble, result: R });
  }
  return { R, trace, bits };
}

/* ----------------------------------------------------------
   Helpers
   ---------------------------------------------------------- */

const fmt = (P) => (P ? `(${P[0]}, ${P[1]})` : 'O (infinity)');

const onCurve = (P, p) =>
  !P || mod(P[1] * P[1], p) === mod(P[0] ** 3 + A * P[0] + B, p);

const samePoint = (P, Q) => P && Q && P[0] === Q[0] && P[1] === Q[1];

/** All affine points of the curve over Z_p. */
function curvePoints(p) {
  const pts = [];
  for (let x = 0; x < p; x++) {
    for (let y = 0; y < p; y++) {
      if (mod(y * y, p) === mod(x ** 3 + B, p)) pts.push([x, y]);
    }
  }
  return pts;
}

/* ==========================================================
   Part 1: table of multiples G … kG
   ========================================================== */

const P17 = 17;
const G17 = [1, 5];

let steps = [];
let cur = 0;
let timer = null;

function compute() {
  const k = Math.max(1, Math.min(60, parseInt($('kIn').value) || 6));

  $('setup').innerHTML =
    `Curve: y² = x³ + 7 (mod 17) &nbsp;|&nbsp; G = (1, 5) &nbsp;|&nbsp; private key k = ${k}<br>` +
    `Check G on curve: 5² mod 17 = ${mod(25, 17)}, 1³+7 = 8 → ${onCurve(G17, P17)}`;

  // Build the list of multiples nG, n = 1 … k
  steps = [];
  let R = null;
  for (let n = 1; n <= k; n++) {
    const info = {};
    const prev = R;
    R = n === 1 ? G17 : add(G17, prev, P17, info);
    steps.push({ n, prev, P: R, info, onC: onCurve(R, P17) });
  }

  renderTable(k);
  renderFastCheck(k);

  cur = 0;
  draw();
}

function describeOperation(s) {
  if (s.n === 1) return 'G';
  if (s.info.kind === 'double') return 'double';
  if (s.info.kind === 'O') return 'P + (−P) = O';
  return `G + ${s.n - 1}G`;
}

function renderTable(k) {
  let html = '<tr><th>n</th><th>nG</th><th>operation</th><th>λ</th><th>on curve?</th></tr>';

  steps.forEach((s) => {
    html +=
      `<tr class="${s.n === k ? 'pub' : ''}">` +
      `<td>${s.n}</td>` +
      `<td>${fmt(s.P)}</td>` +
      `<td>${describeOperation(s)}</td>` +
      `<td>${s.info.lam ?? '–'}</td>` +
      `<td>${s.onC ? '✓' : '✗'}</td>` +
      `</tr>`;
  });

  $('tbl').innerHTML = html;
}

function renderFastCheck(k) {
  const last = steps[steps.length - 1];
  const dd = scalarMul(k, G17, P17);

  const same = dd.R && last.P ? samePoint(dd.R, last.P) : dd.R === last.P;
  const ones = [...dd.bits].filter((b) => b === '1').length;

  $('fast').innerHTML =
    `<b>TOY PUBLIC KEY</b> Q = ${k}·G = <b>${fmt(last.P)}</b><br>` +
    `Double-and-add check: k = ${k} = ${dd.bits}₂ → ${fmt(dd.R)} ${same ? '✓ same result' : '✗'} ` +
    `(${dd.bits.length} doublings, ${ones - 1} extra additions instead of ${k - 1} additions).`;
}

/* ==========================================================
   Part 2: step-by-step visualisation
   ========================================================== */

function drawPlot() {
  const svg = $('plot');
  const p = P17;
  const S = 24; // cell size
  const M = 24; // margin
  const px = (x) => M + x * S;
  const py = (y) => M + (16 - y) * S;

  let h = '';

  // Grid + axis labels
  for (let i = 0; i < p; i++) {
    h += `<line x1="${px(i)}" y1="${py(0)}" x2="${px(i)}" y2="${py(16)}" stroke="var(--line)" stroke-width=".5"/>`;
    h += `<line x1="${px(0)}" y1="${py(i)}" x2="${px(16)}" y2="${py(i)}" stroke="var(--line)" stroke-width=".5"/>`;
    h += `<text x="${px(i)}" y="${py(0) + 16}" font-size="9" text-anchor="middle" fill="var(--mut)">${i}</text>`;
    h += `<text x="${px(0) - 8}" y="${py(i) + 3}" font-size="9" text-anchor="end" fill="var(--mut)">${i}</text>`;
  }

  // All curve points (grey)
  curvePoints(p).forEach(([x, y]) => {
    h += `<circle cx="${px(x)}" cy="${py(y)}" r="4" fill="var(--pt)"/>`;
  });

  // Multiples computed so far
  steps.slice(0, cur).forEach((s, i) => {
    if (!s.P) return;
    const isCurrent = i === cur - 1;

    if (isCurrent && s.prev) {
      h += `<line x1="${px(s.prev[0])}" y1="${py(s.prev[1])}" x2="${px(s.P[0])}" y2="${py(s.P[1])}" ` +
           `stroke="var(--acc)" stroke-dasharray="4 3"/>`;
    }
    h += `<circle cx="${px(s.P[0])}" cy="${py(s.P[1])}" r="${isCurrent ? 8 : 6}" ` +
         `fill="${isCurrent ? 'var(--acc)' : 'none'}" stroke="var(--acc)" stroke-width="2"/>`;
    h += `<text x="${px(s.P[0]) + 9}" y="${py(s.P[1]) - 7}" font-size="10" fill="var(--fg)">${s.n}G</text>`;
  });

  svg.innerHTML = h;
}

function explanationText() {
  if (cur === 0) {
    return 'Press "Next" to start. Step n computes nG = G + (n−1)G.';
  }

  const s = steps[cur - 1];

  if (s.n === 1) {
    return '1G = G = (1, 5). Starting point.';
  }

  if (s.info.kind === 'O') {
    return `${s.n}G = G + ${s.n - 1}G = ${fmt(s.prev)} + (1,5): same x, opposite y → vertical line → ` +
           `result is the point at infinity O. The cycle closes: order of G = ${s.n}.`;
  }

  const [x1, y1] = G17;
  const [x2, y2] = s.prev;

  const lambdaLine =
    s.info.kind === 'double'
      ? `λ = 3x₁² / (2y₁) = 3·${x1}² / ${2 * y1}  → ${mod(3 * x1 * x1, P17)} · inv(${mod(2 * y1, P17)}) = ${s.info.lam}\n`
      : `λ = (y₂−y₁)/(x₂−x₁) = (${y2}−${y1})/(${x2}−${x1}) = ${mod(y2 - y1, P17)} · inv(${mod(x2 - x1, P17)}) = ${s.info.lam}   (mod 17)\n`;

  return (
    `${s.n}G = G + ${s.n - 1}G = (1,5) + ${fmt(s.prev)}\n` +
    lambdaLine +
    `x₃ = λ² − x₁ − x₂ = ${s.info.lam}² − ${x1} − ${x2} ≡ ${s.P[0]}\n` +
    `y₃ = λ(x₁ − x₃) − y₁ = ${s.info.lam}(${x1} − ${s.P[0]}) − ${y1} ≡ ${s.P[1]}\n` +
    `→ ${s.n}G = ${fmt(s.P)}   on curve: ${s.onC ? 'yes' : 'NO'}` +
    (s.n === steps.length ? '\n★ This is the TOY PUBLIC KEY' : '')
  );
}

function draw() {
  drawPlot();
  $('explain').textContent = explanationText();
}

/* ----------------------------------------------------------
   Navigation controls
   ---------------------------------------------------------- */

function stepTo(delta) {
  cur = Math.max(0, Math.min(steps.length, cur + delta));
  draw();
}

function stopTimer() {
  clearInterval(timer);
  timer = null;
}

function togglePlay() {
  if (timer) {
    stopTimer();
    return;
  }
  if (cur >= steps.length) cur = 0;

  timer = setInterval(() => {
    if (cur >= steps.length) stopTimer();
    else stepTo(1);
  }, 1500);
}

function resetPlot() {
  stopTimer();
  cur = 0;
  draw();
}

/* ==========================================================
   Part 3: statistical exploration
   ========================================================== */

/** Number of 1-bits in n. */
function popcnt(n) {
  let count = 0;
  while (n) {
    count += n & 1;
    n >>= 1;
  }
  return count;
}

/** Draw a 16-bin histogram of `vals` (range 0…p-1) into the SVG with the given id. */
function hist(id, vals, p, label) {
  const bins = 16;
  const counts = new Array(bins).fill(0);
  vals.forEach((v) => counts[Math.min(bins - 1, Math.floor((v * bins) / p))]++);

  const max = Math.max(...counts, 1);
  let h = `<text x="4" y="12" font-size="11" fill="var(--fg)">${label}: histogram (16 bins over 0…p−1)</text>`;

  counts.forEach((n, i) => {
    const barHeight = (100 * n) / max;
    h += `<rect class="bar" x="${10 + i * 25}" y="${130 - barHeight}" width="21" height="${barHeight}"/>`;
    h += `<text x="${20 + i * 25}" y="144" font-size="8" text-anchor="middle" fill="var(--mut)">${n}</text>`;
  });

  $(id).innerHTML = h;
}

/** Pearson correlation coefficient of two equal-length arrays. */
function pearson(a, b) {
  const n = a.length;
  const meanA = a.reduce((s, v) => s + v, 0) / n;
  const meanB = b.reduce((s, v) => s + v, 0) / n;

  let sab = 0;
  let saa = 0;
  let sbb = 0;
  for (let i = 0; i < n; i++) {
    sab += (a[i] - meanA) * (b[i] - meanB);
    saa += (a[i] - meanA) ** 2;
    sbb += (b[i] - meanB) ** 2;
  }
  return sab / Math.sqrt(saa * sbb || 1);
}

/** Find a generator: fixed for p = 17, otherwise the first valid point found. */
function findGenerator(p) {
  if (p === 17) return [1, 5];

  for (let x = 1; x < p; x++) {
    const rhs = mod(x ** 3 + B, p);
    for (let y = 1; y < p; y++) {
      if (mod(y * y, p) === rhs) return [x, y];
    }
  }
}

function runStats() {
  const p = parseInt($('curveSel').value);
  const G = findGenerator(p);

  // pts[k] = kG, with pts[0] = O, until the cycle closes
  const pts = [null, G];
  let R = G;
  while (true) {
    R = add(R, G, p);
    if (!R) break;
    pts.push(R);
  }
  const n = pts.length; // order of G

  // Collect samples for k = 1 … n-1
  const ks = [];
  const xs = [];
  const ys = [];
  for (let k = 1; k < n; k++) {
    ks.push(k);
    xs.push(pts[k][0]);
    ys.push(pts[k][1]);
  }

  // Bit statistics
  const bitsLen = Math.ceil(Math.log2(p));
  const onesX = xs.reduce((sum, v) => sum + popcnt(v), 0);
  const onesY = ys.reduce((sum, v) => sum + popcnt(v), 0);

  // x-value frequency
  const freq = new Map();
  xs.forEach((v) => freq.set(v, (freq.get(v) || 0) + 1));
  const distinctX = freq.size;
  const maxFreq = Math.max(...freq.values());

  // Avalanche tests: (k vs k xor 1) and (k vs k+1)
  let xorX = 0, xorY = 0, xorCount = 0;
  let incX = 0, incY = 0, incCount = 0;
  for (let k = 1; k < n; k++) {
    const k2 = k ^ 1;
    if (k2 >= 1 && k2 < n) {
      xorX += popcnt(pts[k][0] ^ pts[k2][0]);
      xorY += popcnt(pts[k][1] ^ pts[k2][1]);
      xorCount++;
    }
    if (k + 1 < n) {
      incX += popcnt(pts[k][0] ^ pts[k + 1][0]);
      incY += popcnt(pts[k][1] ^ pts[k + 1][1]);
      incCount++;
    }
  }

  // Symmetry: kG and (n-k)G share x
  const symmetric = pts.filter(
    (q, k) => k && k < n && pts[n - k] && pts[n - k][0] === q[0]
  ).length;

  const frac = (value, count) => (value / (count * bitsLen)).toFixed(3);

  $('statsOut').innerHTML =
    `p = ${p}, G = ${fmt(G)}, <b>order of G = ${n}</b> (so nG = O; k and k+${n} give the same point)<br>` +
    `Points sampled: k = 1…${n - 1}. Distinct x values: ${distinctX}; max frequency of a single x: ${maxFreq}.<br>` +
    `Bit statistics (${bitsLen} bits per coordinate): fraction of 1-bits in x = ${frac(onesX, xs.length)}, ` +
    `in y = ${frac(onesY, ys.length)} (uniform ideal ≈ 0.5, slightly lower because values < p).<br>` +
    `Pearson correlation: corr(k, x) = ${pearson(ks, xs).toFixed(3)}, corr(k, y) = ${pearson(ks, ys).toFixed(3)} ` +
    `(≈ 0 = no linear correlation).<br>` +
    `Avalanche test, flipping the last bit of k (k vs k⊕1): fraction of changed bits ` +
    `x = ${frac(xorX, xorCount)}, y = ${frac(xorY, xorCount)}; ` +
    `(k vs k+1): x = ${frac(incX, incCount)}, y = ${frac(incY, incCount)}. ` +
    `A hash would give ≈ 0.5 for every flip with no structure.<br>` +
    `Symmetry: kG and (n−k)G always share x and have opposite y (${symmetric} of ${n - 1} points) – ` +
    `a structural property no hash has.`;

  hist('histX', xs, p, 'x');
  hist('histY', ys, p, 'y');
}

/* ==========================================================
   Event wiring + initial render
   ========================================================== */

$('calc').onclick = compute;
$('next').onclick = () => stepTo(1);
$('prev').onclick = () => stepTo(-1);
$('reset').onclick = resetPlot;
$('play').onclick = togglePlay;
$('runStats').onclick = runStats;

compute();
runStats();
