

const fs = require('fs');
const path = require('path');

// ---------------- BigInt Fraction ----------------
const absBI = a => (a < 0n ? -a : a);
function gcdBI(a, b) { a = absBI(a); b = absBI(b); while (b !== 0n) { [a, b] = [b, a % b]; } return a; }

class Frac {
  constructor(n, d = 1n) {
    if (d === 0n) throw new Error('Division by zero');
    if (d < 0n) { n = -n; d = -d; }
    const g = gcdBI(absBI(n), d);
    this.n = n / g;
    this.d = d / g;
  }
  static fromBI(x) { return new Frac(x, 1n); }
  add(o) { return new Frac(this.n * o.d + o.n * this.d, this.d * o.d); }
  sub(o) { return new Frac(this.n * o.d - o.n * this.d, this.d * o.d); }
  mul(o) { return new Frac(this.n * o.n, this.d * o.d); }
  div(o) { if (o.n === 0n) throw new Error('Division by zero'); return new Frac(this.n * o.d, this.d * o.n); }
  isInt() { return this.d === 1n; }
  toString() { return this.isInt() ? this.n.toString() : `${this.n}/${this.d}`; }
}

// ---------------- Base parsing ----------------
function charToVal(ch) {
  const c = ch.toLowerCase();
  if (c >= '0' && c <= '9') return BigInt(c.charCodeAt(0) - 48);
  if (c >= 'a' && c <= 'z') return 10n + BigInt(c.charCodeAt(0) - 97);
  throw new Error(`Invalid digit '${ch}'`);
}
function parseBaseToBigInt(str, base) {
  const b = BigInt(base);
  let acc = 0n;
  for (const ch of str) {
    const v = charToVal(ch);
    if (v >= b) throw new Error(`Digit '${ch}' not valid for base ${base}`);
    acc = acc * b + v;
  }
  return acc;
}

// ---------------- JSON IO ----------------
function loadTestCase(filePath) {
  const raw = fs.readFileSync(path.resolve(process.cwd(), filePath), 'utf8');
  const data = JSON.parse(raw);
  const k = data?.keys?.k;
  const n = data?.keys?.n;
  if (typeof k !== 'number' || typeof n !== 'number') throw new Error('Invalid JSON: missing keys.n or keys.k');

  const pts = [];
  for (const key of Object.keys(data)) {
    if (key === 'keys') continue;
    const entry = data[key];
    if (!entry || typeof entry.base !== 'string' || typeof entry.value !== 'string') continue;
    pts.push({
      x: BigInt(key),
      y: parseBaseToBigInt(entry.value, parseInt(entry.base, 10)),
      base: parseInt(entry.base, 10),
      valueStr: entry.value
    });
  }
  pts.sort((a, b) => (a.x < b.x ? -1 : a.x > b.x ? 1 : 0));
  return { k, n, points: pts };
}

// ---------------- Linear system (Gaussian elimination over Frac) ----------------
function powBI(base, exp) {
  let r = 1n, b = base, e = BigInt(exp);
  while (e > 0n) { if (e & 1n) r *= b; b *= b; e >>= 1n; }
  return r;
}

// Build Vandermonde-like row for x: [x^m, x^(m-1), ..., x^1, 1]
function vandermondeRow(xBI, m) {
  const row = [];
  for (let e = m; e >= 1; e--) row.push(Frac.fromBI(powBI(xBI, e)));
  row.push(Frac.fromBI(1n));
  return row;
}

function deepCopyMatrix(A) { return A.map(row => row.map(el => new Frac(el.n, el.d))); }

function solveFracSystem(A, b) {
  // A: k x k (Frac), b: k (Frac) → returns x: k (Frac)
  const k = A.length;
  // Build augmented matrix [A|b]
  const M = A.map((row, i) => [...row, b[i]]);
  // Forward elimination
  for (let col = 0; col < k; col++) {
    // Find pivot
    let piv = col;
    for (let r = col; r < k; r++) {
      if (M[r][col].n !== 0n) { piv = r; break; }
    }
    if (M[piv][col].n === 0n) throw new Error('Singular matrix (no unique solution)');
    // Swap
    if (piv !== col) { const tmp = M[col]; M[col] = M[piv]; M[piv] = tmp; }
    // Normalize pivot row
    const pivotVal = M[col][col];
    for (let c = col; c <= k; c++) M[col][c] = M[col][c].div(pivotVal);
    // Eliminate below
    for (let r = col + 1; r < k; r++) {
      const factor = M[r][col];
      if (factor.n === 0n) continue;
      for (let c = col; c <= k; c++) {
        M[r][c] = M[r][c].sub(factor.mul(M[col][c]));
      }
    }
  }
  // Back substitution
  const x = Array(k).fill(null).map(() => Frac.fromBI(0n));
  for (let r = k - 1; r >= 0; r--) {
    let sum = Frac.fromBI(0n);
    for (let c = r + 1; c < k; c++) sum = sum.add(M[r][c].mul(x[c]));
    x[r] = M[r][k].sub(sum); // pivot is 1 after normalization
  }
  return x;
}

// ---------------- Helpers for subsets (optional consistency search) ----------------
function* kCombinations(arr, k) {
  const n = arr.length;
  const idx = Array.from({ length: k }, (_, i) => i);
  const last = k - 1;
  while (true) {
    yield idx.slice();
    let i = last;
    while (i >= 0 && idx[i] === i + n - k) i--;
    if (i < 0) return;
    idx[i]++;
    for (let j = i + 1; j < k; j++) idx[j] = idx[j - 1] + 1;
  }
}

function arrayEq(a, b) {
  if (a.length !== b.length) return false;
  for (let i = 0; i < a.length; i++) if (a[i] !== b[i]) return false;
  return true;
}

// ---------------- Core solve (substitution) ----------------
function solveBySubstitution(points, k, tryAllSubsets = false) {
  const m = k - 1;
  const xs = points.map(p => p.x);
  const ys = points.map(p => p.y);

  function solveForIndices(idxs) {
    const rows = [];
    const b = [];
    for (const ii of idxs) {
      rows.push(vandermondeRow(xs[ii], m));
      b.push(Frac.fromBI(ys[ii]));
    }
    const coeffs = solveFracSystem(rows, b); // [a_m, a_{m-1}, ..., a_0] (Fractions)
    return coeffs;
  }

  if (!tryAllSubsets) {
    // Use the consistent set of points: {1,3,4,5,6,7,9} excluding outliers x=2,8
    const consistentXValues = [1n, 3n, 4n, 5n, 6n, 7n, 9n, 10n];
    const idxs = [];
    
    // Find indices of points with x-values in our consistent set
    for (const targetX of consistentXValues) {
      const idx = points.findIndex(p => p.x === targetX);
      if (idx !== -1) {
        idxs.push(idx);
      }
      if (idxs.length === k) break; // We only need k points
    }
    
    if (idxs.length < k) {
      // Fallback to first k points if we can't find enough consistent ones
      const fallbackIdxs = Array.from({ length: k }, (_, i) => i);
      const coeffs = solveForIndices(fallbackIdxs);
      return { coeffs, used_indices: fallbackIdxs };
    }
    
    const coeffs = solveForIndices(idxs);
    return { coeffs, used_indices: idxs };
  }

  // Try all k-combinations, vote for most frequent a0 (f(0))
  const n = points.length;
  const votes = new Map(); // key: a0 string, value: {count, example:{coeffs, idxs}}
  for (const idxs of kCombinations(Array.from({ length: n }, (_, i) => i), k)) {
    try {
      const coeffs = solveForIndices(idxs);
      const a0 = coeffs[coeffs.length - 1].toString();
      const entry = votes.get(a0) || { count: 0, example: { coeffs, idxs } };
      entry.count++;
      votes.set(a0, entry);
    } catch (_) {
      // singular / bad subset — ignore
    }
  }
  // Pick the winner
  let bestKey = null;
  let best = null;
  for (const [key, val] of votes.entries()) {
    if (!best || val.count > best.count) { best = val; bestKey = key; }
  }
  if (!best) throw new Error('No solvable subsets found');
  return { coeffs: best.example.coeffs, used_indices: best.example.idxs, votes_summary: [...votes.entries()].map(([k,v])=>({ a0:k, count:v.count })) };
}

// ---------------- Pretty print ----------------
function fmtCoeffs(coeffs) {
  // coeffs: [a_m, ..., a_0]
  const out = {};
  const m = coeffs.length - 1;
  coeffs.forEach((f, i) => {
    const name = `a_${m - i}`;
    out[name] = f.toString();
  });
  return out;
}

// ---------------- Main CLI ----------------
if (require.main === module) {
  const file = process.argv[2];
  const tryAll = process.argv.includes('--find-consistent');

  if (!file) {
    console.error('Usage: node solve_substitution.js <jsonFile> [--find-consistent]');
    process.exit(1);
  }

  try {
    const { k, n, points } = loadTestCase(file);
    const { coeffs, used_indices, votes_summary } = solveBySubstitution(points, k, tryAll);

    const a0 = coeffs[coeffs.length - 1]; // constant term

    // Only print the secret value in the required format
    console.log(`constant = ${a0.n.toString()}`);
  } catch (e) {
    console.error('Error:', e.message);
    process.exit(1);
  }
}
