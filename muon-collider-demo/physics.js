/* JavaScript port of the mucsim formulas (display only; the Python package is the reference).
   Conventions: opc = 1 + cos(theta), omc = 1 - cos(theta); GeV units; tree level, massless final fermions. */
(function (root) {
  const PI = Math.PI, GEV2PB = 0.3893793721e9, MMU = 0.1056583755, ALPHA = 0.00729735;
  const FERM = { e: [-0.5, -1, 1], mu: [-0.5, -1, 1], tau: [-0.5, -1, 1], u: [0.5, 2 / 3, 3], d: [-0.5, -1 / 3, 3],
                 s: [-0.5, -1 / 3, 3], c: [0.5, 2 / 3, 3], b: [-0.5, -1 / 3, 3] };
  const DEFAULT = { alpha: ALPHA, mz: 91.1876, gz: 2.4952, sw2: 0.2312, sw2eff: 0.2315 };

  function coup(f, p) {
    const [T3, Q, Nc] = FERM[f];
    const sc = Math.sqrt(p.sw2 * (1 - p.sw2));
    const se = p.sw2eff == null ? p.sw2 : p.sw2eff;
    return [(T3 - Q * se) / sc, (-Q * se) / sc, Q, Nc];
  }
  function chi(s, p) {                                  // s / (s - mz^2 + i mz gz)  as [re, im]
    const d = s - p.mz * p.mz, w = p.mz * p.gz, den = d * d + w * w;
    return [s * d / den, -s * w / den];
  }
  function ffSums(s, init, fin, p) {
    const [gLi, gRi, Qi] = coup(init, p), [gLf, gRf, Qf, Nc] = coup(fin, p), [cr, ci] = chi(s, p);
    const A2 = (gi, gf) => { const re = Qi * Qf + gi * gf * cr, im = gi * gf * ci; return re * re + im * im; };
    return { same: A2(gLi, gLf) + A2(gRi, gRf), opp: A2(gLi, gRf) + A2(gRi, gLf), Nc };
  }
  const sigmaFF = (s, fin, p, init = 'mu') => { const q = ffSums(s, init, fin, p);
    return GEV2PB * q.Nc * PI * p.alpha ** 2 / (3 * s) * (q.same + q.opp); };
  const afbFF = (s, fin, p, init = 'mu') => { const q = ffSums(s, init, fin, p); return 0.75 * (q.same - q.opp) / (q.same + q.opp); };
  const dsFF = (s, opc, omc, fin, p, init = 'mu') => { const q = ffSums(s, init, fin, p);
    return GEV2PB * q.Nc * PI * p.alpha ** 2 / (8 * s) * (q.same * opc * opc + q.opp * omc * omc); };

  function dsBhabha(s, opc, omc, p, mode = 'full') {
    const [gL, gR, Q] = coup('mu', p);
    const t = -0.5 * s * omc, um = -0.5 * s * opc, [cr, ci] = chi(s, p), mz2 = p.mz * p.mz;
    const As = (g1, g2) => [Q * Q + g1 * g2 * cr, g1 * g2 * ci];
    const At = (g1, g2) => Q * Q + g1 * g2 * t / (t - mz2);
    const comb = (g) => {
      const [ar, ai] = As(g, g), sr = ar / s, si = ai / s, at = At(g, g) / t;
      return mode === 'full' ? (sr + at) ** 2 + si ** 2 : sr * sr + si * si + at * at;
    };
    let S = um * um * (comb(gL) + comb(gR));
    const a = As(gL, gR); S += 2 * t * t * (a[0] ** 2 + a[1] ** 2) / (s * s);
    const b = At(gL, gR); S += 2 * s * s * b * b / (t * t);
    return GEV2PB * PI * p.alpha ** 2 * S / (2 * s);
  }
  function sigmaBhabha(s, thmin, p, mode = 'full', n = 4000) {     // sigma(theta in [thmin, pi - thmin]) [pb]
    const wmin = 2 * Math.sin(0.5 * thmin) ** 2, lo = Math.log(wmin);
    let acc = 0;                                         // midpoint rule in ln(1-cos), forward + backward halves
    for (let i = 0; i < n; i++) {
      const w = Math.exp(lo + (i + 0.5) * (-lo) / n);
      acc += (dsBhabha(s, 2 - w, w, p, mode) + dsBhabha(s, w, 2 - w, p, mode)) * w;
    }
    return acc * (-lo) / n;
  }

  /* ---- muon decay in flight ---- */
  const G = (x, sp) => sp === 'numu' ? 3 * x * x - (4 / 3) * x ** 3 : 6 * x * x - 4 * x ** 3;
  const mcdf = (x, sp) => sp === 'numu' ? 2 * x ** 3 - x ** 4 : 4 * x ** 3 - 3 * x ** 4;
  const eEta = (g) => { const E = g * MMU, pp = Math.sqrt((E - MMU) * (E + MMU)); return (E + pp) / MMU; };
  function labPdf(y, g, sp) {
    const beta = Math.sqrt(1 - 1 / (g * g)), omb = 1 / (g * g * (1 + beta));
    const xlo = 2 * y / (1 + beta), xhi = Math.min(1, 2 * y / omb);
    if (!(y > 0) || xlo >= 1) return 0;
    return (G(xhi, sp) - G(xlo, sp)) / beta;
  }
  function coneFrac(theta, g) { const k = (eEta(g) * Math.tan(0.5 * theta)) ** 2; return k / (1 + k); }
  function sampleX(sp) { const r = Math.random(); let lo = 0, hi = 1;
    for (let i = 0; i < 40; i++) { const m = 0.5 * (lo + hi); if (mcdf(m, sp) < r) lo = m; else hi = m; } return 0.5 * (lo + hi); }
  function sampleNu(g, n, sp) {
    const ee = eEta(g), y = new Float64Array(n), tg = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      const x = sampleX(sp), Es = 0.5 * MMU * x, a = Math.random(), flip = Math.random() < 0.5;
      const opc = flip ? a : 2 - a, omc = flip ? 2 - a : a;
      const pp = Es * ee * opc, pm = Es * omc / ee;
      y[i] = 0.5 * (pp + pm) / (g * MMU);
      tg[i] = 2 * Math.atan2(Math.sqrt(pm), Math.sqrt(pp)) * g;
    }
    return { y, tg };
  }
  function sampleFF(s, n, fin, p) {                      // exact inversion sampling of cos(theta)
    const q = ffSums(s, 'mu', fin, p), ps = q.same / (q.same + q.opp), c = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      const r = Math.random(), cb = Math.cbrt(r), big = 2 * cb, small = 2 * (1 - r) / (1 + cb + cb * cb);
      const same = Math.random() < ps, opc = same ? big : small, omc = same ? small : big;
      c[i] = opc < omc ? opc - 1 : 1 - omc;
    }
    return c;
  }

  /* ---- precision lab: naive f32 / naive f64 / stable f64 / 80-digit reference ---- */
  const f32 = Math.fround;
  function kernels(sqrtS, eps, side, D) {
    D.set({ precision: 80 });
    const gamma = sqrtS / 2 / MMU, small = Math.pow(10, eps);
    const out = [];
    // K1: 1 - beta
    {
      const naive = (r) => { const s = r(r(sqrtS) * r(sqrtS)), m = r(MMU); return r(1 - r(Math.sqrt(r(1 - r(r(4 * r(m * m)) / s))))); };
      const rr = 4 * MMU * MMU / (sqrtS * sqrtS), stable = rr / (1 + Math.sqrt(1 - rr));
      const dm = new D(MMU), ds = new D(sqrtS), dr = dm.mul(dm).mul(4).div(ds.mul(ds));
      out.push({ name: '1 − β of a beam muon', f32: naive(f32), f64: naive((x) => x), stable, ref: new D(1).minus(new D(1).minus(dr).sqrt()) });
    }
    // K2: m^2 reconstructed from boosted (E, pz)
    {
      const naive = (r) => { const g = r(gamma), m = r(MMU), beta = r(Math.sqrt(r(1 - r(1 / r(g * g))))), E = r(g * m), pz = r(r(g * beta) * m);
        return r(r(E * E) - r(pz * pz)); };
      const E = gamma * MMU, ee = (E + Math.sqrt((E - MMU) * (E + MMU))) / MMU, stable = (MMU * ee) * (MMU / ee);
      out.push({ name: 'm² of a boosted muon from (E, p_z)', f32: naive(f32), f64: naive((x) => x), stable, ref: new D(MMU).pow(2) });
    }
    const opc = side === 'bwd' ? small : 2 - small, omc = side === 'bwd' ? 2 - small : small;
    const cExact = side === 'bwd' ? new D(small).minus(1) : new D(1).minus(small);
    const cNaive = side === 'bwd' ? opc - 1 : 1 - omc;
    const dg = new D(gamma), dbeta = new D(1).minus(new D(1).div(dg.mul(dg))).sqrt();
    // K3: lab E_nu / E_mu for x = 1/2
    {
      const x = 0.5;
      const naive = (r) => { const g = r(gamma), beta = r(Math.sqrt(r(1 - r(1 / r(g * g))))); return r(r(0.5 * x) * r(1 + r(beta * r(cNaive)))); };
      const beta = Math.sqrt(1 - 1 / (gamma * gamma)), omb = 1 / (gamma * gamma * (1 + beta));
      out.push({ name: 'lab E_ν/E_μ  (x = 1/2)', f32: naive(f32), f64: naive((v) => v), stable: 0.5 * x * (omb + beta * opc),
                 ref: new D(x).div(2).mul(new D(1).plus(dbeta.mul(cExact))) });
    }
    // K4: lab polar angle
    {
      const naive = (r) => { const g = r(gamma), beta = r(Math.sqrt(r(1 - r(1 / r(g * g))))), c = r(cNaive);
        const cl = r(r(c + beta) / r(1 + r(beta * c))); return Math.acos(Math.max(-1, Math.min(1, cl))); };
      const E = gamma * MMU, ee = (E + Math.sqrt((E - MMU) * (E + MMU))) / MMU, stable = 2 * Math.atan2(Math.sqrt(omc), ee * Math.sqrt(opc));
      out.push({ name: 'lab polar angle θ of a decay ν [rad]', f32: naive(f32), f64: naive((v) => v), stable,
                 ref: new D(cExact.plus(dbeta)).div(new D(1).plus(dbeta.mul(cExact))).acos() });
    }
    return out;
  }
  function relErr(v, ref, D) {
    if (!isFinite(v)) return 1;
    const e = new D(v).minus(ref).abs().div(ref.abs()).toNumber();
    return Math.min(e, 1);
  }

  root.MUC = { DEFAULT, MMU, FERM, sigmaFF, afbFF, dsFF, dsBhabha, sigmaBhabha, labPdf, coneFrac, eEta, sampleNu, sampleFF, kernels, relErr };
})(typeof window !== 'undefined' ? window : globalThis);
