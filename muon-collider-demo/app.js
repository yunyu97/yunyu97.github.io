/* UI for the muon-collider demo. Physics lives in physics.js. */
(function () {
  const M = window.MUC;
  const $ = (id) => document.getElementById(id);
  const css = (v) => getComputedStyle(document.documentElement).getPropertyValue(v).trim();
  let DATA = null, mc = { c: null, nu: null }, charts = {};

  const FLAB = { e: ['e⁺e⁻', 'e⁻'], tau: ['τ⁺τ⁻', 'τ⁻'], u: ['u ū', 'u'], d: ['d d̄', 'd'], s: ['s s̄', 's'], c: ['c c̄', 'c'], b: ['b b̄', 'b'] };
  const fmtE = (x, d = 3) => (x === 0 ? '0' : x.toExponential(d));
  const fmtG = (x, d = 4) => Number(x.toPrecision(d)).toString();
  const params = () => ({ alpha: M.DEFAULT.alpha, mz: +$('mz').value, gz: +$('gz').value, sw2: M.DEFAULT.sw2, sw2eff: +$('sw').value });
  const sqrtS = () => 1000 * Math.pow(10, +$('tev').value);               // GeV

  function mkChart(id, opts) {
    const ctx = $(id).getContext('2d');
    if (charts[id]) charts[id].destroy();
    Chart.defaults.color = css('--text');
    Chart.defaults.borderColor = css('--line');
    charts[id] = new Chart(ctx, { type: 'scatter', data: { datasets: opts.datasets }, options: {
      animation: false, responsive: true, maintainAspectRatio: false,
      plugins: { legend: { labels: { boxWidth: 14, font: { size: 11 } } } },
      scales: { x: { type: opts.xlog ? 'logarithmic' : 'linear', title: { display: true, text: opts.xt }, ...(opts.xr || {}) },
                y: { type: opts.ylog ? 'logarithmic' : 'linear', title: { display: true, text: opts.yt },
                     ...(opts.ylog ? { ticks: { callback: (v) => { const e = Math.log10(v); return Math.abs(e - Math.round(e)) < 1e-9 ? '1e' + Math.round(e) : ''; } } } : {}),
                     ...(opts.yr || {}) } } } });
  }
  const line = (label, pts, color, dash) => ({ label, data: pts, showLine: true, pointRadius: 0, borderColor: color, borderWidth: 2, borderDash: dash || [], tension: 0 });
  const dots = (label, pts, color, style) => ({ label, data: pts, showLine: false, pointRadius: 3.5, pointStyle: style || 'circle',
      backgroundColor: style === 'rect' ? 'transparent' : color, borderColor: color, borderWidth: style === 'rect' ? 1.5 : 0 });
  const hist = (vals, lo, hi, nb, norm) => { const h = new Array(nb).fill(0), w = (hi - lo) / nb;
    for (const v of vals) { const k = Math.floor((v - lo) / w); if (k >= 0 && k < nb) h[k]++; }
    return h.map((c, i) => ({ x: lo + (i + 0.5) * w, y: c / (vals.length * w) })); };
  const histCounts = (counts, edges) => { const n = counts.reduce((a, b) => a + b, 0);
    return counts.map((c, i) => ({ x: 0.5 * (edges[i] + edges[i + 1]), y: c / (n * (edges[i + 1] - edges[i])) })); };

  function update() {
    const p = params(), sq = sqrtS(), s = sq * sq, fin = $('fin').value, thmin = +$('thmin').value;
    const g = sq / 2 / M.MMU;
    $('o_tev').textContent = sq < 1000 ? sq.toPrecision(3) + ' GeV' : (sq / 1000).toPrecision(3) + ' TeV';
    $('o_mz').textContent = (+$('mz').value).toFixed(4); $('o_gz').textContent = (+$('gz').value).toFixed(4);
    $('o_sw').textContent = (+$('sw').value).toFixed(4); $('o_th').textContent = thmin.toFixed(3);
    const bh = M.sigmaBhabha(s, thmin, p, 'full'), bhi = M.sigmaBhabha(s, thmin, p, 'incoherent');
    const r = 4 * M.MMU ** 2 / s;
    $('s_sig').textContent = fmtE(M.sigmaFF(s, fin, p)); $('s_afb').textContent = fmtG(M.afbFF(s, fin, p));
    $('s_bh').textContent = fmtE(bh); $('s_int').textContent = (100 * (bh / bhi - 1)).toFixed(2) + ' %';
    $('s_g').textContent = fmtG(g, 5); $('s_omb').textContent = fmtE(r / (1 + Math.sqrt(1 - r)));

    // chart A: cos theta
    const [fl, fpart] = FLAB[fin], bf = M.betaF(s, fin), mf = M.FMASS[fin] || 0;
    const massTxt = mf ? `m = ${mf} GeV, β_f = ${bf.toFixed(4)}` : 'massless quark';
    const sig = M.sigmaFF(s, fin, p), ana = [], ana0 = [];
    const sig0 = M.sigmaFF(s, fin, p, 'mu', false);
    for (let i = 0; i <= 200; i++) { const c = -1 + i / 100; ana.push({ x: c, y: M.dsFF(s, 1 + c, 1 - c, fin, p) / sig });
      ana0.push({ x: c, y: M.dsFF(s, 1 + c, 1 - c, fin, p, 'mu', false) / sig0 }); }
    const dsA = [line(`analytic, μ⁺μ⁻ → ${fl} (${massTxt})`, ana, css('--c4'))];
    if (mf && Math.abs(bf - 1) > 2e-4) dsA.push(line('analytic, massless limit', ana0, css('--c3'), [6, 4]));
    if (mc.c && mc.cKey === [sq, p.mz, p.gz, p.sw2eff, fin].join()) dsA.push(dots(`browser MC (${fl}), N = ${mc.c.length.toLocaleString()}`, hist(mc.c, -1, 1, 20), css('--c1')));
    const pk = DATA && DATA.pythia.ff_hist, sel = pk && (Math.abs(sq - 1000) < 1 ? '1000' : Math.abs(sq - 10000) < 1 ? '10000' : null);
    const pythOK = sel && fin === 'e' && Math.abs(p.mz - 91.1876) < 1e-6 && Math.abs(p.gz - 2.4952) < 1e-6 && Math.abs(p.sw2eff - 0.2315) < 1e-6;
    if ($('pyth').checked && pythOK) dsA.push(dots('Pythia 8, μ⁺μ⁻ → e⁺e⁻ (200k events)', histCounts(pk[sel].counts, pk[sel].edges), css('--c2'), 'rect'));
    $('tA').innerHTML = `cosθ of ${fpart} in μ⁺μ⁻ → ${fl} <span class="badge">${massTxt}</span>`;
    $('l_sig').textContent = `σ(μ⁺μ⁻ → ${fl}) [pb]`;
    $('s_beta').textContent = mf ? bf.toFixed(6) : '1 (massless)';
    const dsig = 100 * (sig0 / sig - 1);
    $('fnote').innerHTML = (fin === 'e' || fin === 'tau')
      ? `e, μ and τ have identical electroweak couplings, so μ⁺μ⁻ → e⁺e⁻ and μ⁺μ⁻ → τ⁺τ⁻ differ <b>only through the final-state mass</b>. ` +
        `For ${fl} here, neglecting the mass would change σ by ${dsig.toPrecision(2)} % and A_FB by ${(M.afbFF(s, fin, p, 'mu', false) - M.afbFF(s, fin, p)).toExponential(1)}; ` +
        `it matters near threshold (see the chart below), not at TeV energies.`
      : `Quark final states are leading-order partonic estimates with <b>massless quarks and no QCD corrections</b>` +
        (sq < 200 ? ` <b>— not reliable at this low energy</b> (thresholds, resonances, hadronisation).` : `; they are not compared with Pythia8.`);
    mkChart('chA', { datasets: dsA, xt: `cos θ of ${fpart} w.r.t. the μ⁻ beam`, yt: '(1/σ) dσ/dcosθ', yr: { min: 0 } });

    // chart C: mu mu -> mu mu vs angle
    const full = [], inc = [], thLo = Math.max(thmin, 0.01);
    for (let i = 0; i <= 120; i++) { const th = thLo * Math.pow((Math.PI / 2) / thLo, i / 120), w = 2 * Math.sin(th / 2) ** 2, o = 2 * Math.cos(th / 2) ** 2;
      full.push({ x: th, y: M.dsBhabha(s, o, w, p, 'full') }); inc.push({ x: th, y: M.dsBhabha(s, o, w, p, 'incoherent') }); }
    mkChart('chC', { datasets: [line('s+t with interference', full, css('--c1')), line('incoherent |M_s|²+|M_t|²', inc, css('--c2'), [6, 4])],
                     xlog: true, ylog: true, xt: 'θ [rad]', yt: 'dσ/dcosθ [pb]' });

    // charts B: neutrino from decay
    const ys = [], ts = [];
    for (let i = 1; i < 400; i++) { const y = i / 400; ys.push({ x: y, y: M.labPdf(y, g, 'numu') }); }
    for (let i = 1; i <= 250; i++) { const t = 0.02 * i; ts.push({ x: t, y: (M.coneFrac((t + 0.01) / g, g) - M.coneFrac((t - 0.01) / g, g)) / 0.02 }); }
    const dB1 = [line('exact analytic', ys, css('--c4'))], dB2 = [line('analytic d/dt [κ/(1+κ)]', ts, css('--c4'))];
    if (mc.nu && mc.nuKey === g) { dB1.push(dots('browser MC', hist(mc.nu.y, 0, 1, 40), css('--c1'))); dB2.push(dots('browser MC', hist(mc.nu.tg, 0, 5, 40), css('--c1'))); }
    if ($('pyth').checked && DATA) { const d = DATA.pythia.decay.numu;
      dB1.push(dots('Pythia 8 (5 TeV μ)', histCounts(d.y_hist, d.y_edges), css('--c2'), 'rect'));
      dB2.push(dots('Pythia 8 (5 TeV μ)', histCounts(d.theta_gamma_hist, d.theta_gamma_edges).map((q) => q), css('--c2'), 'rect')); }
    mkChart('chB1', { datasets: dB1, xt: 'E_ν / E_μ', yt: 'pdf', yr: { min: 0 } });
    mkChart('chB2', { datasets: dB2, xt: 'θ × γ', yt: 'pdf', yr: { min: 0 } });
    $('pynote').textContent = pythOK ? 'Pythia 8 overlay for μ⁺μ⁻ → e⁺e⁻ is shown because the parameters match the reference run (1 or 10 TeV, default electroweak inputs).'
      : 'Pythia 8 cross-section histograms exist only for e⁺e⁻ at 1 and 10 TeV with the default electroweak inputs. The neutrino overlay (5 TeV muons) is valid at any energy because the shape depends on θ·γ and E/E_μ only for γ ≫ 1.';
    lab();
  }

  function generate() {
    const p = params(), sq = sqrtS(), s = sq * sq, fin = $('fin').value, n = +$('nev').value, g = sq / 2 / M.MMU;
    $('genmsg').textContent = 'generating…';
    setTimeout(() => {
      const t0 = performance.now();
      mc.c = M.sampleFF(s, n, fin, p); mc.cKey = [sq, p.mz, p.gz, p.sw2eff, fin].join();
      mc.nu = M.sampleNu(g, n, 'numu'); mc.nuKey = g;
      $('genmsg').textContent = `${n.toLocaleString()} events in ${((performance.now() - t0) / 1000).toFixed(2)} s`;
      update();
    }, 20);
  }

  // ---------------- precision lab ----------------
  function errCell(e) { const c = e < 1e-12 ? 'good' : e < 1e-6 ? 'warn' : 'bad'; return `<td class="num ${c}">${e >= 1 ? '1 (none)' : fmtE(e, 2)}</td>`; }
  function lab() {
    const eps = +$('eps').value; $('o_eps').textContent = eps.toFixed(1);
    const rows = M.kernels(sqrtS(), eps, $('side').value, window.Decimal);
    let h = '<tr><th>quantity</th><th>reference (80 digits)</th><th>naive f32</th><th>naive f64</th><th>stable f64</th>' +
            '<th>err naive f32</th><th>err naive f64</th><th>err stable f64</th></tr>';
    for (const r of rows) {
      h += `<tr><td>${r.name}</td><td class="num">${r.ref.toSignificantDigits(15).toExponential()}</td>` +
           `<td class="num">${fmtE(r.f32, 6)}</td><td class="num">${fmtE(r.f64, 12)}</td><td class="num">${fmtE(r.stable, 12)}</td>` +
           errCell(M.relErr(r.f32, r.ref, Decimal)) + errCell(M.relErr(r.f64, r.ref, Decimal)) + errCell(M.relErr(r.stable, r.ref, Decimal)) + '</tr>';
    }
    $('plab').innerHTML = h;
  }

  function tauChart() {
    const pp = { alpha: 0.00729735, mz: 91.1876, gz: DATA.pythia.settings.GammaZ, sw2: 0.2312, sw2eff: 0.2315 }, pts = [], n = 80;
    for (let i = 0; i <= n; i++) { const sq = 3.6 * Math.pow(1000 / 3.6, i / n), s = sq * sq;
      pts.push({ x: sq, y: Math.max(Math.abs(M.sigmaFF(s, 'tau', pp, 'mu', false) / M.sigmaFF(s, 'tau', pp) - 1), 1e-12) }); }
    const py = DATA.pythia.ff_lowE.map((e) => ({ x: e.ecm, y: Math.max(Math.abs(M.sigmaFF(e.ecm * e.ecm, 'tau', pp, 'mu', false) / e.sigma_pb - 1), 1e-12) }));
    const pyMass = DATA.pythia.ff_lowE.map((e) => ({ x: e.ecm, y: Math.max(Math.abs(M.sigmaFF(e.ecm * e.ecm, 'tau', pp) / e.sigma_pb - 1), 1e-12) }));
    mkChart('chD', { datasets: [line('massless vs finite mass (analytic)', pts, css('--c1')), dots('massless analytic vs Pythia 8', py, css('--c2'), 'rect'),
      dots('finite-mass analytic vs Pythia 8', pyMass, css('--c3'))],
      xlog: true, ylog: true, xt: '√s [GeV]', yt: '|relative error of σ(μ⁺μ⁻→τ⁺τ⁻)|', yr: { min: 1e-6, max: 2 } });
  }

  function scanChart() {
    const k = $('kernel').value, K = DATA.precision.kernels[k], keys = Object.keys(K.naive_f64).sort((a, b) => Number(a) - Number(b));
    const xs = keys.map(Number);
    const mk = (impl, label, color, dash) => line(label, keys.map((key) => ({ x: Number(key), y: Math.max(K[impl][key].max, 1e-18) })), color, dash);
    mkChart('chP', { datasets: [mk('naive_f32', 'naive float32', css('--c2'), [6, 4]), mk('naive_f64', 'naive float64', '#E69F00'), mk('stable_f64', 'stable float64', css('--c1')),
      line('float64 epsilon', [{ x: xs[0], y: DATA.precision.eps64 }, { x: xs[xs.length - 1], y: DATA.precision.eps64 }], css('--text'), [2, 3])],
      xlog: true, ylog: true, xt: '√s [TeV]', yt: 'max relative error', yr: { min: 1e-18, max: 5 } });
  }

  // ---------------- verification tables ----------------
  function tables() {
    const d = DATA, s = d.validation_summary;
    $('vstats').innerHTML = [[`${s.n_passed} / ${s.n_checks}`, 'verification checks passed'],
      [`${d.mutation.n_detected} / ${d.mutation.n_mutants}`, 'injected defects detected'],
      [fmtE(Math.max(...Object.values(d.precision.summary).map((r) => r.naive_f64_max))), 'worst naive float64 error at 10 TeV'],
      [fmtE(Math.max(...Object.values(d.precision.summary).map((r) => r.stable_f64_max))), 'worst stable float64 error at 10 TeV']]
      .map(([a, b]) => `<div class="stat"><b>${a}</b><span>${b}</span></div>`).join('');
    $('vchecks').innerHTML = '<tr><th>group</th><th>check</th><th>metric</th><th>value</th><th>criterion</th><th></th></tr>' +
      d.checks.map((c) => `<tr><td>${c.group}</td><td>${c.name}</td><td>${c.metric}</td><td class="num">${fmtG(c.value, 3)}</td>` +
        `<td class="num">${c.op} ${fmtG(c.threshold, 3)}</td><td class="${c.passed ? 'good' : c.informational ? 'warn' : 'bad'}">${c.passed ? 'PASS' : c.informational ? 'baseline' : 'FAIL'}</td></tr>`).join('');
    $('vmut').innerHTML = '<tr><th>injected defect</th><th>caught by</th><th>first failing check</th></tr>' +
      d.mutation.mutants.map((m) => `<tr><td>${m.description}</td><td>${m.layers.join(', ') || 'NOT DETECTED'}</td><td>${m.first_failure || '–'}</td></tr>`).join('');
    $('vsur').innerHTML = '<tr><th>target</th><th>interp. median</th><th>interp. p95</th><th>extrap. median</th><th>extrap. p95</th><th>σ error @10 TeV</th><th>σ error @100 TeV</th><th>ns/point (NN vs closed form)</th></tr>' +
      Object.entries(d.surrogate).map(([k, m]) => `<tr><td>${k}</td><td class="num">${fmtG(m.interp.median_rel_err, 3)}</td><td class="num">${fmtG(m.interp.p95_rel_err, 3)}</td>` +
        `<td class="num">${fmtG(m.extrap.median_rel_err, 3)}</td><td class="num">${fmtG(m.extrap.p95_rel_err, 3)}</td><td class="num">${fmtG(m.sigma_forward['10TeV'].rel_err, 3)}</td>` +
        `<td class="num">${fmtG(m.sigma_forward['100TeV'].rel_err, 3)}</td><td class="num">${fmtG(m.ns_per_eval_surrogate, 3)} vs ${fmtG(m.ns_per_eval_analytic, 3)}</td></tr>`).join('');
    const g = d.codegen_gate;
    $('vgate').innerHTML = `<tr><th>kernel</th><th>candidates until pass</th><th>rejected classes</th><th>final max rel. error</th></tr>` +
      g.runs.map((r) => `<tr><td>${r.task}</td><td>${r.iterations}</td><td>${r.defects_caught.join(', ') || 'none'}</td><td class="num">${fmtE(r.final_max_rel_err, 2)}</td></tr>`).join('') +
      `<tr><td colspan="4">${g.n_passed} of ${g.n_tasks} kernels passed; ${g.total_rejected} of ${g.total_candidates} candidates rejected. Backend: ${g.backend}.</td></tr>`;
    $('vdev').innerHTML = '<tr><th>defect</th><th>found by</th><th>evidence</th></tr>' +
      d.dev_log.map((e) => `<tr><td>${e.defect}</td><td>${e.found_by}</td><td>${e.evidence}</td></tr>`).join('');
    $('kernel').innerHTML = Object.entries(d.precision.kernels).map(([k, v]) => `<option value="${k}">${v.label}</option>`).join('');
    $('stress').textContent = 'stress domain: ' + d.precision.stress_domain;
    $('kernel').addEventListener('change', scanChart);
    scanChart();
    tauChart();
  }

  ['tev', 'mz', 'gz', 'sw', 'thmin', 'fin', 'pyth'].forEach((id) => $(id).addEventListener('input', update));
  ['eps', 'side'].forEach((id) => $(id).addEventListener('input', lab));
  $('gen').addEventListener('click', generate);
  fetch('data/demo_data.json').then((r) => r.json()).then((d) => { DATA = d; tables(); update(); })
    .catch((e) => { console.error(e); update(); });
  window.addEventListener('load', () => { if (!DATA) update(); });
  window.__muc = { update, lab };
})();
