
const fs = require("fs");
const path = require("path");

// ---------------- BigInt helpers ----------------
const absBI = (a) => (a < 0n ? -a : a);
function gcdBI(a, b) {
  a = absBI(a);
  b = absBI(b);
  while (b !== 0n) {
    const t = a % b;
    a = b;
    b = t;
  }
  return a;
}

// ---------------- Fraction over BigInt ----------------
class Frac {
  constructor(n, d = 1n) {
    if (d === 0n) throw new Error("Division by zero");
    if (d < 0n) {
      n = -n;
      d = -d;
    }
    const g = gcdBI(absBI(n), d);
    this.n = n / g;
    this.d = d / g;
  }
  static fromBI(x) {
    return new Frac(x, 1n);
  }
  add(o) {
    return new Frac(this.n * o.d + o.n * this.d, this.d * o.d);
  }
  mul(o) {
    return new Frac(this.n * o.n, this.d * o.d);
  }
  isInt() {
    return this.d === 1n;
  }
  toString() {
    return this.isInt() ? this.n.toString() : `${this.n}/${this.d}`;
  }
}

// ---------------- Base decoding ----------------
function charToVal(ch) {
  const c = ch.toLowerCase();
  if (c >= "0" && c <= "9") return BigInt(c.charCodeAt(0) - 48);
  if (c >= "a" && c <= "z") return 10n + BigInt(c.charCodeAt(0) - 97);
  throw new Error(`Invalid digit '${ch}'`);
}
function parseBaseToBigInt(str, base) {
  if (!(base >= 2 && base <= 36))
    throw new Error(`Unsupported base ${base} (2..36)`);
  const b = BigInt(base);
  let acc = 0n;
  for (const ch of str) {
    const v = charToVal(ch);
    if (v >= b) throw new Error(`Digit '${ch}' not valid for base ${base}`);
    acc = acc * b + v;
  }
  return acc;
}

// ---------------- JSON input ----------------
function loadTestCase(filePath) {
  const raw = fs.readFileSync(path.resolve(process.cwd(), filePath), "utf8");
  const data = JSON.parse(raw);
  const k = data?.keys?.k;
  const n = data?.keys?.n;
  if (typeof k !== "number" || typeof n !== "number") {
    throw new Error("Invalid JSON: missing keys.n or keys.k");
  }
  const pts = [];
  for (const key of Object.keys(data)) {
    if (key === "keys") continue;
    const e = data[key];
    if (!e) continue;
    const base = parseInt(e.base, 10);
    const y = parseBaseToBigInt(e.value, base);
    pts.push({ x: BigInt(key), y });
  }
  pts.sort((a, b) => (a.x < b.x ? -1 : a.x > b.x ? 1 : 0));
  if (pts.length < k)
    throw new Error(`Need at least k=${k} points, got ${pts.length}`);
  return { k, n, points: pts };
}

// ---------------- Lagrange at x=0 (O(k^2)) ----------------
function constantViaLagrange(points) {
  const k = points.length;
  const xs = points.map((p) => p.x);
  const ys = points.map((p) => p.y);

  let acc = Frac.fromBI(0n);
  for (let i = 0; i < k; i++) {
    let num = 1n,
      den = 1n;
    const xi = xs[i];
    for (let j = 0; j < k; j++) {
      if (j === i) continue;
      num *= -xs[j]; // multiply by (-x_j)
      den *= xi - xs[j]; // multiply by (x_i - x_j)
    }
    const li0 = new Frac(num, den);
    acc = acc.add(li0.mul(Frac.fromBI(ys[i])));
  }
  if (!acc.isInt())
    throw new Error(`Non-integer constant term: ${acc.toString()}`);
  return acc.n;
}

// ---------------- Subset generation (k-combinations) ----------------
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

// ---------------- CLI helpers ----------------
function parsePickArg(argv) {
  const pickArg = argv.find((a) => a.startsWith("--pick="));
  if (!pickArg) return null;
  const raw = pickArg.slice("--pick=".length).trim();
  if (!raw) return null;
  const parts = raw
    .split(",")
    .map((s) => s.trim())
    .filter(Boolean);
  if (parts.length === 0) return null;
  const keys = parts.map((s) => {
    if (!/^-?\d+$/.test(s))
      throw new Error(`--pick contains non-integer: '${s}'`);
    return BigInt(s);
  });
  // dedupe
  const seen = new Set(),
    unique = [];
  for (const v of keys) {
    const k = v.toString();
    if (!seen.has(k)) {
      seen.add(k);
      unique.push(v);
    }
  }
  return unique;
}

// ---------------- Main ----------------
if (require.main === module) {
  const file = process.argv[2];
  if (!file) {
    console.error(
      "Usage: node secret_consistent.js <jsonFile> [--pick=x1,...,xk]",
    );
    process.exit(1);
  }

  try {
    const { k, n, points } = loadTestCase(file);

    // If user forces a subset via --pick, use it directly
    const chosenXs = parsePickArg(process.argv);
    if (chosenXs) {
      if (chosenXs.length !== k) {
        throw new Error(
          `--pick must specify exactly k=${k} x-keys; got ${chosenXs.length}`,
        );
      }
      const map = new Map(points.map((p, i) => [p.x.toString(), i]));
      const subset = chosenXs.map((x) => {
        const idx = map.get(x.toString());
        if (idx === undefined)
          throw new Error(`--pick refers to missing x=${x.toString()}`);
        return points[idx];
      });
      const c = constantViaLagrange(subset);
      console.log(`constant = ${c.toString()}`);
      process.exit(0);
    }

    // If n == k, compute once
    if (n === k || points.length === k) {
      const c = constantViaLagrange(points.slice(0, k));
      console.log(`constant = ${c.toString()}`);
      process.exit(0);
    }

    // n > k: vote across all k-subsets; pick most frequent integer constant
    const idxs = [...Array(points.length).keys()];
    const counts = new Map(); // constant(BigInt)->count
    const example = new Map(); // constant(BigInt)->subset indices

    for (const combo of kCombinations(idxs, k)) {
      const subset = combo.map((i) => points[i]);
      try {
        const c = constantViaLagrange(subset);
        const key = c.toString();
        counts.set(key, (counts.get(key) || 0) + 1);
        if (!example.has(key)) example.set(key, combo.slice());
      } catch (_) {
        // Non-integer or singular subset; ignore
      }
    }

    if (counts.size === 0)
      throw new Error("No consistent (integer) subsets found");

    // Pick the constant with maximum votes
    let bestConst = null,
      bestCount = -1;
    for (const [kstr, ct] of counts.entries()) {
      if (ct > bestCount) {
        bestCount = ct;
        bestConst = kstr;
      }
    }

    // Final answer
    console.log(`constant = ${bestConst}`);
  } catch (e) {
    console.error("Error:", e.message);
    process.exit(1);
  }
}
