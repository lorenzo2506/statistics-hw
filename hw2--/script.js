// ===================== core: secp256k1 over BigInt =====================
const P = 2n ** 256n - 2n ** 32n - 977n;
const GX = 0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798n;
const GY = 0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8n;
const G = [GX, GY];
const mod = (a) => { const r = a % P; return r >= 0n ? r : r + P; };
function inv(a) {
  let [r0, r1] = [mod(a), P], [s0, s1] = [1n, 0n];
  while (r1) { const q = r0 / r1; [r0, r1] = [r1, r0 - q * r1]; [s0, s1] = [s1, s0 - q * s1]; }
  return mod(s0);
}
function add(A, B) { // null = point at infinity
  if (!A) return B; if (!B) return A;
  if (A[0] === B[0]) {
    if (mod(A[1] + B[1]) === 0n) return null;
    const l = mod(3n * A[0] * A[0] * inv(2n * A[1]));
    const x = mod(l * l - 2n * A[0]);
    return [x, mod(l * (A[0] - x) - A[1])];
  }
  const l = mod((B[1] - A[1]) * inv(B[0] - A[0]));
  const x = mod(l * l - A[0] - B[0]);
  return [x, mod(l * (A[0] - x) - A[1])];
}
const neg = (A) => A && [A[0], mod(-A[1])];
function mul(k, A) { // double-and-add
  let R = null, Q = A;
  while (k > 0n) { if (k & 1n) R = add(R, Q); Q = add(Q, Q); k >>= 1n; }
  return R;
}
const same = (A, B) => (!A && !B) || (A && B && A[0] === B[0] && A[1] === B[1]);
const key = (A) => A ? A[0].toString(16) + (A[1] & 1n) : 'O';

// table m*G for m = 0..255 (message encoding)
const TAB = [null]; for (let m = 1; m < 256; m++) TAB.push(add(TAB[m - 1], G));
const LOOK = new Map(TAB.map((pt, m) => [key(pt), m]));
const rand128 = () => { const a = new Uint32Array(4); (globalThis.crypto || require('crypto').webcrypto).getRandomValues(a); return [...a].reduce((s, v) => (s << 32n) | BigInt(v), 0n) || 1n; };
const pub = (d) => mul(BigInt(d), G);
function encChar(m, Q) { const r = rand128(); return [mul(r, G), add(TAB[m], mul(r, Q))]; }
function decChar(ct, d) { const S = mul(BigInt(d), ct[0]); const M = add(ct[1], neg(S)); return LOOK.get(key(M)); }

// attacks on the restricted key space (1..65535)
const KMAX = 65535;
function bruteForce(Q, done, prog) {
  const t0 = Date.now(); let R = G, i = 1;
  (function loop() {
    const end = Math.min(i + 3000, KMAX + 1);
    for (; i < end; i++) { if (same(R, Q)) { done({ d: i, steps: i, ms: Date.now() - t0 }); return; } R = add(R, G); }
    if (i > KMAX) { done(null); return; }
    if (prog) prog(i); setTimeout(loop, 0);
  })();
}
function bsgs(Q) {
  const t0 = Date.now(), m = 256, baby = new Map(); let B = null, ops = 0;
  for (let j = 0; j < m; j++) { baby.set(key(B), j); B = add(B, G); ops++; }
  const step = neg(mul(BigInt(m), G)); let Y = Q;
  for (let i = 0; i <= 256; i++) {
    const j = baby.get(key(Y)); ops++;
    if (j !== undefined) { const d = i * m + j; if (d >= 1 && same(pub(d), Q)) return { d, steps: ops, ms: Date.now() - t0 }; }
    Y = add(Y, step);
  }
  return null;
}

// ===================== Caesar and statistics =====================
const A0 = 65;
function caesar(s, k) { return s.replace(/[A-Za-z]/g, (c) => { const b = c <= 'Z' ? 65 : 97; return String.fromCharCode(((c.charCodeAt(0) - b + k) % 26 + 26) % 26 + b); }); }
const letters = (s) => s.toUpperCase().replace(/[^A-Z]/g, '');
const norm = (a) => { const f = a.map(v => Math.max(v, 0.05) / 100), t = f.reduce((x, y) => x + y, 0); return f.map(v => v / t); };
const LANGS = {
  en: { name: 'English', ref: norm([8.167,1.492,2.782,4.253,12.702,2.228,2.015,6.094,6.966,0.153,0.772,4.025,2.406,6.749,7.507,1.929,0.095,5.987,6.327,9.056,2.758,0.978,2.360,0.150,1.974,0.074]),
    words: ['THE', 'AND', 'THAT', 'WITH', 'ARE', 'FOR', 'NOT', 'THIS', 'FROM', 'HAVE', 'ION', 'ING'],
    A: "Statistics is the science of learning from data. When we collect observations, we want to understand what they tell us about the world, and how much we can trust the conclusions that we draw. In cryptography, the same ideas appear in a different form. An attacker who sees an encrypted message tries to find patterns in it, because a good cipher should hide every pattern, while a weak cipher leaves traces of the original language. Letters in English are not used equally often: the letter e appears much more than the letter z, and some pairs of letters, such as th and he, are very common. These regularities are the starting point of statistical analysis, and they have been used for centuries to break simple ciphers. In the modern world we rely on mathematics instead of secrecy of the method, and the security of a system depends only on the key that is kept private.",
    B: "The history of writing secret messages is very old. Soldiers and kings used to send letters that nobody else could read, and they trusted simple rules to protect their plans. A famous example is the shift of every letter by a fixed number of places in the alphabet. This method looks clever, but there are only a few possible keys, and a patient person can try them all in a short time. Later, people noticed that the frequency of the letters gives away the secret, because the shift moves every letter in the same way. Today computers can test millions of keys every second, so a modern system must have a huge number of possible keys. It must also hide the statistics of the message, so that the encrypted text looks like random noise to anyone who does not know the key. These two requirements explain why the design of a good cipher is a difficult task.",
    sample: 'Statistical analysis of letter frequencies recovers the key of a Caesar cipher when the ciphertext is long enough to reveal the pattern of the language' },
  it: { name: 'Italian', ref: norm([11.74,0.92,4.50,3.73,11.79,0.95,1.64,1.54,11.28,0,0,6.51,2.51,6.88,9.83,3.05,0.51,6.37,4.98,5.62,3.01,2.10,0,0,0,0.49]),
    words: ['CHE', 'DEL', 'DELLA', 'NON', 'UNA', 'PER', 'CON', 'SONO', 'QUESTO', 'ZIONE', 'MENTE', 'ANTE'],
    A: "La statistica studia i dati per capire come funziona il mondo e ci dice quanto possiamo fidarci delle conclusioni che ricaviamo. Anche nella crittografia compaiono le stesse idee. Chi intercetta un messaggio cifrato cerca degli schemi ricorrenti, dato che un buon cifrario dovrebbe nascondere ogni traccia, mentre un cifrario debole lascia intravedere la lingua originale. In italiano le lettere non vengono usate con la stessa frequenza: la vocale a, la vocale e e la vocale i compaiono molto spesso, mentre la lettera z compare raramente, e alcune coppie di lettere, come di, che e re, sono molto comuni. Queste ripetizioni sono il punto di partenza della statistica applicata ai cifrari, e da secoli servono per rompere i cifrari semplici. Oggi la sicurezza di un sistema non dipende dalla segretezza del metodo ma soltanto dalla chiave che resta privata.",
    B: "Scrivere messaggi segreti ha una storia molto antica. Soldati e sovrani mandavano lettere che nessun altro poteva leggere e si affidavano a regole semplici per proteggere i loro piani. Un esempio famoso consiste nello spostare ogni lettera di un numero fisso di posti nell ordine alfabetico. Questo metodo sembra ingegnoso, ma le chiavi possibili sono poche e una persona paziente riesce a provarle tutte in poco tempo. Con il tempo gli studiosi notarono che la frequenza delle lettere tradisce il segreto, dato che lo spostamento muove ogni lettera nello stesso modo. Oggi i computer provano milioni di chiavi al secondo, quindi un sistema moderno deve avere moltissime chiavi possibili. Deve anche nascondere le statistiche del messaggio, in modo che il testo cifrato sembri rumore casuale a chiunque non conosca la chiave. Queste due condizioni mostrano quanto sia difficile progettare un buon cifrario.",
    sample: 'L analisi statistica delle frequenze delle lettere permette di ricavare la chiave di un cifrario di Cesare quando il testo cifrato ha abbastanza lettere da mostrare lo schema della lingua' }
};
let CUR = LANGS.en;
function counts(s) { const c = new Array(26).fill(0); for (const ch of s) c[ch.charCodeAt(0) - A0]++; return c; }
function chi2(s) { const c = counts(s), n = s.length || 1; let x = 0; for (let i = 0; i < 26; i++) { const e = n * CUR.ref[i]; x += (c[i] - e) ** 2 / e; } return x; }
function buildBigram(text) {
  const t = letters(text), c = Array.from({ length: 26 }, () => new Array(26).fill(0.5));
  for (let i = 0; i + 1 < t.length; i++) c[t.charCodeAt(i) - A0][t.charCodeAt(i + 1) - A0]++;
  return c.map(row => { const s = row.reduce((a, b) => a + b, 0); return row.map(v => Math.log(v / s)); });
}
function bigramLL(s, M) { let x = 0; for (let i = 0; i + 1 < s.length; i++) x += M[s.charCodeAt(i) - A0][s.charCodeAt(i + 1) - A0]; return x; }
function wordHits(s) { let n = 0; for (const w of CUR.words) { let p = -1; while ((p = s.indexOf(w, p + 1)) >= 0) n++; } return n; }
function scoreAll(cipher, M) { // for each key k: candidate = decrypt with shift k
  const c = letters(cipher); const out = [];
  for (let k = 0; k < 26; k++) { const cand = caesar(c, -k); out.push({ k, cand, chi: chi2(cand), bi: bigramLL(cand, M), w: wordHits(cand) }); }
  return out;
}
const argBest = (rows, f, sign) => { const best = Math.max(...rows.map(r => sign * f(r))); return rows.filter(r => Math.abs(sign * f(r) - best) < 1e-9).map(r => r.k); };
function successTest(M, trials) {
  const B = letters(CUR.B), lens = [5, 10, 20, 40, 80, 160], res = [];
  for (const L of lens) {
    let a = 0, b = 0, c = 0;
    for (let t = 0; t < trials; t++) {
      const st = Math.floor(Math.random() * (B.length - L)), k = Math.floor(Math.random() * 26);
      const rows = scoreAll(caesar(B.substr(st, L), k), M);
      const pick = (f, s) => { const ks = argBest(rows, f, s); return ks.includes(k) ? 1 / ks.length : 0; };
      a += pick(r => r.chi, -1); b += pick(r => r.bi, 1); c += pick(r => r.w, 1);
    }
    res.push({ L, chi: a / trials, bi: b / trials, w: c / trials });
  }
  return res;
}
// statistics used to compare ciphertexts
function ic(s) { const c = counts(s), n = s.length; let x = 0; for (const v of c) x += v * (v - 1); return x / (n * (n - 1)); } // index of coincidence
function chiUniform(s) { const c = counts(s), e = s.length / 26; return c.reduce((a, v) => a + (v - e) ** 2 / e, 0); }
// toy ECC ciphertext reduced to letters: (x of C2) mod 26 for each plaintext letter
function eccLetters(text, Q, done, prog) {
  const t = letters(text), out = []; let i = 0;
  (function loop() {
    const end = Math.min(i + 10, t.length);
    for (; i < end; i++) { const ct = encChar(t.charCodeAt(i), Q); out.push(ct[1] ? String.fromCharCode(65 + Number(ct[1][0] % 26n)) : 'A'); }
    if (i < t.length) { if (prog) prog(i, t.length); setTimeout(loop, 0); } else done(out.join(''));
  })();
}
function attackSummary(c, M) {
  const rows = scoreAll(c, M), chis = rows.map(r => r.chi).sort((a, b) => a - b);
  const k = argBest(rows, r => r.chi, -1)[0], best = rows.find(r => r.k === k);
  return { k, cand: best.cand, min: chis[0], med: (chis[12] + chis[13]) / 2, words: best.w };
}

// ===================== UI =====================
function initUI() {
  const $ = (id) => document.getElementById(id);
  const hx = (v) => v === undefined ? '–' : '0x' + v.toString(16).padStart(64, '0').slice(0, 12) + '…';
  const px = (A) => A ? hx(A[0]) : 'O';
  let cts = [];
  const getD = () => Math.max(1, Math.min(KMAX, parseInt($('d').value) || 1));
  function info() { const d = getD(), Q = pub(d); $('eccInfo').innerHTML = `Curve secp256k1. Private key d = ${d} (restricted range). Public key Q = d·G:<br>x = 0x${Q[0].toString(16)}<br>y = 0x${Q[1].toString(16)}`; return Q; }
  function table(rows, head) { return '<tr>' + head.map(h => `<th>${h}</th>`).join('') + '</tr>' + rows.map(r => '<tr>' + r.map(c => `<td>${c}</td>`).join('') + '</tr>').join(''); }
  $('enc').onclick = () => {
    const Q = info(), msg = $('msg').value; cts = []; const rows = [];
    for (const ch of msg) { const m = ch.charCodeAt(0) < 256 ? ch.charCodeAt(0) : 63; const ct = encChar(m, Q); cts.push(ct); rows.push([ch === ' ' ? '␣' : ch, m, px(ct[0]), px(ct[1])]); }
    $('eccTbl').innerHTML = table(rows, ['char', 'm', 'C₁.x = (r·G).x', 'C₂.x = (m·G + r·Q).x']);
    $('eccOut').innerHTML = `<p>Encrypted ${cts.length} characters. Each is a pair of points; press Decrypt.</p>`;
  };
  $('dec').onclick = () => {
    if (!cts.length) { $('eccOut').textContent = 'Encrypt something first.'; return; }
    const d = getD(); const out = cts.map(ct => { const m = decChar(ct, d); return m === undefined ? '?' : String.fromCharCode(m); }).join('');
    $('eccOut').innerHTML = `<p>Decrypted with d = ${d}: <b>${out.replace(/</g, '&lt;')}</b></p>`;
  };
  $('demoA').onclick = () => {
    const Q = info(); const rows = [], c = [];
    for (let i = 0; i < 8; i++) { const ct = encChar(65, Q); rows.push(['A', 65, px(ct[0]), px(ct[1])]); c.push(caesar('A', parseInt($('ck').value) || 3)); }
    $('eccTbl').innerHTML = table(rows, ['char', 'm', 'C₁.x', 'C₂.x']);
    $('eccOut').innerHTML = `<p>Eight encryptions of the same letter A: every ciphertext is different. With Caesar (key ${parseInt($('ck').value) || 3}) the same message would give <b>${c.join('')}</b>, the same letter every time, which is what frequency analysis exploits.</p>`;
  };
  $('cgo').onclick = () => {
    const m = $('cmsg').value, k = ((parseInt($('ck').value) || 0) % 26 + 26) % 26, c = caesar(m, k);
    $('cOut').innerHTML = `<p>Plaintext: <b>${m}</b><br>Ciphertext: <b>${c}</b><br>Decrypted: <b>${caesar(c, -k)}</b></p>`;
  };
  function bars(id, counts, n) {
    let h = ''; const w = 18, mx = Math.max(...counts.map(c => c / (n || 1)), ...CUR.ref);
    counts.forEach((c, i) => { const f = c / (n || 1), bh = 100 * f / mx, rh = 100 * CUR.ref[i] / mx;
      h += `<rect class="bar" x="${10 + i * 19.5}" y="${110 - bh}" width="${w}" height="${bh}"/><rect class="ref" x="${10 + i * 19.5}" y="${110 - rh}" width="${w}" height="${rh}"/><text x="${19 + i * 19.5}" y="124" font-size="9" text-anchor="middle" fill="var(--mut)">${String.fromCharCode(65 + i)}</text>`; });
    $(id).innerHTML = h;
  }
  $('an').onclick = () => {
    const M = buildBigram($('corpus').value), c = $('cipher').value, rows = scoreAll(c, M);
    const bc = argBest(rows, r => r.chi, -1), bb = argBest(rows, r => r.bi, 1), bw = argBest(rows, r => r.w, 1);
    let h = '<tr><th>key k</th><th>candidate plaintext (first 28 letters)</th><th>χ² (min)</th><th>bigram log-lik (max)</th><th>word hits (max)</th></tr>';
    rows.forEach(r => { h += `<tr class="${bb.includes(r.k) ? 'best' : ''}"><td>${r.k}</td><td>${r.cand.slice(0, 28)}</td><td>${r.chi.toFixed(1)}${bc.includes(r.k) ? ' ★' : ''}</td><td>${r.bi.toFixed(1)}${bb.includes(r.k) ? ' ★' : ''}</td><td>${r.w}${bw.includes(r.k) ? ' ★' : ''}</td></tr>`; });
    $('anOut').innerHTML = `<p>Estimated key — χ²: <b>${bc.join(',')}</b>; bigram: <b>${bb.join(',')}</b>; word search: <b>${bw.join(',')}</b> (★ = best for that method; highlighted row = bigram choice).</p><div class="scroll"><table>${h}</table></div>`;
    const L = letters(c); bars('freq', counts(L), L.length);
  };
  $('test').onclick = () => {
    $('testNote').textContent = 'Running…';
    setTimeout(() => {
      const res = successTest(buildBigram($('corpus').value), 400);
      $('testTbl').innerHTML = table(res.map(r => [r.L, (100 * r.chi).toFixed(0) + '%', (100 * r.bi).toFixed(0) + '%', (100 * r.w).toFixed(0) + '%']), ['letters', 'χ² letter frequency', 'bigram score', 'word search']);
      $('testNote').textContent = 'Success rate = fraction of trials where the true key is recovered (ties split randomly). Limitations: short texts do not show the language statistics; the corpus must match the language and style; non-English texts need a corpus in that language; all methods assume a plain shift (monoalphabetic) cipher.';
    }, 30);
  };
  $('bf').onclick = () => {
    const Q = info(); $('atkOut').textContent = 'Running brute force…';
    bruteForce(Q, (r) => { $('atkOut').innerHTML = r ? `Brute force found d = <b>${r.d}</b> after ${r.steps} point additions (${r.ms} ms).${recover(r.d)}` : 'Key not found.'; }, (i) => { $('atkOut').textContent = `Running brute force… tried ${i} keys`; });
  };
  $('bsgs').onclick = () => {
    const Q = info(); $('atkOut').textContent = 'Running…';
    setTimeout(() => { const r = bsgs(Q); $('atkOut').innerHTML = r ? `Baby-step giant-step found d = <b>${r.d}</b> with about ${r.steps} group operations (${r.ms} ms) instead of up to ${KMAX}.${recover(r.d)}` : 'Key not found.'; }, 20);
  };
  function recover(d) {
    if (!cts.length) return '<br>Encrypt a message in section 1 first to see it decrypted with the recovered key.';
    return '<br>Message decrypted with the recovered key: <b>' + cts.map(ct => { const m = decChar(ct, d); return m === undefined ? '?' : String.fromCharCode(m); }).join('').replace(/</g, '&lt;') + '</b>';
  }
  function setLang(v) { CUR = LANGS[v]; $('corpus').value = CUR.A; $('cipher').value = caesar(CUR.sample, 11); $('testTbl').innerHTML = ''; $('testNote').textContent = ''; $('an').onclick(); }
  $('lang').onchange = () => setLang($('lang').value);
  $('cmp').onclick = () => {
    const t = letters(CUR.B).slice(0, 150), k = ((parseInt($('ck').value) || 0) % 26 + 26) % 26, Q = info(), M = buildBigram($('corpus').value);
    $('cmpOut').textContent = 'Encrypting with ECC…';
    eccLetters(t, Q, (e) => {
      const c = caesar(t, k), f = (v, n) => v.toFixed(n);
      const sets = [['Plaintext', t], ['Caesar ciphertext', c], ['ECC ciphertext (C₂.x mod 26)', e]];
      let h1 = table(sets.map(([n, x]) => [n, f(ic(x), 4), f(chiUniform(x), 1) + (chiUniform(x) < 37.65 ? ' (uniform-like)' : ' (not uniform)'), f(chi2(x), 1)]), ['text (150 letters)', 'index of coincidence', 'χ² vs uniform (df 25, 5% limit 37.65)', 'χ² vs language frequencies']);
      const aC = attackSummary(c, M), aE = attackSummary(e, M);
      let h2 = table([['Caesar', aC.k + (aC.k === k ? ' (correct)' : ' (wrong)'), f(aC.min, 1), f(aC.med, 1), f(aC.med / aC.min, 1), aC.cand.slice(0, 30), aC.words],
        ['ECC', aE.k + ' (meaningless)', f(aE.min, 1), f(aE.med, 1), f(aE.med / aE.min, 1), aE.cand.slice(0, 30), aE.words]],
        ['ciphertext', 'key chosen by χ²', 'min χ² over 26 shifts', 'median χ²', 'median / min', 'best candidate plaintext', 'word hits']);
      $('cmpOut').innerHTML = `<p>True Caesar key: ${k}. Plain text length: ${t.length} letters.</p><div class="scroll"><table>${h1}</table></div><p><b>Frequency attack on both ciphertexts</b> (a clear key shows as min χ² much smaller than the median):</p><div class="scroll"><table>${h2}</table></div>`;
      bars('freqC', counts(c), c.length); bars('freqE', counts(e), e.length);
    }, (i, n) => { $('cmpOut').textContent = `Encrypting with ECC… ${i}/${n}`; });
  };
  info(); $('cgo').onclick(); setLang('en');
}
if (typeof document !== 'undefined') initUI();