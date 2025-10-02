/* =========================
   Constants (pure data)
   ========================= */
const MATERIAL_PROPERTIES = {
  "Custom": { "rho": 0, "E": 0, "nu": 0 },
  "Concrete": { "rho": 2.3e3, "E": 2.1e10, "nu": 0.005 },
  "Gypsum Board": { "rho": 0.8e3, "E": 0.18e10, "nu": 0.005 },
  "Plywood":      { "rho": 0.6e3, "E": 0.5e10,  "nu": 0.30  },
  "Glass":        { "rho": 2500,  "E": 70.0e9,  "nu": 0.23  }
};

const FREQS = [16,20,25,31.5,40,50,63,80,100,125,160,200,250,315,400,500,630,800,1000,1250,1600,2000,2500,3150,4000,5000,6300,8000];
const THIRD_OCTAVE_TICKS = [16, 31.5, 63, 125, 250, 500, 1000, 2000, 4000, 8000];

/* =========================
   State (set after DOM ready)
   ========================= */
let charts = {};
let ui = {}; // will be filled in DOMContentLoaded

/* =========================
   Utilities
   ========================= */
class Complex {
  constructor(re=0, im=0) { this.re = re; this.im = im; }
  add(z) { return new Complex(this.re + z.re, this.im + z.im); }
  sub(z) { return new Complex(this.re - z.re, this.im - z.im); }
  mul(z) { return new Complex(this.re*z.re - this.im*z.im, this.re*z.im + this.im*z.re); }
  div(z) { const d = z.re*z.re + z.im*z.im; if (d===0) return new Complex(Infinity, Infinity); return new Complex((this.re*z.re + this.im*z.im)/d, (this.im*z.re - this.re*z.im)/d); }
  get magnitude() { return Math.sqrt(this.re*this.re + this.im*this.im); }
  get phase() { return Math.atan2(this.im, this.re); }
  static i() { return new Complex(0, 1); }
}

function createTable(tableElement, headers, data) {
  if (!tableElement) return;
  tableElement.innerHTML = `<thead><tr>${headers.map(h => `<th>${h}</th>`).join('')}</tr></thead>
    <tbody>${data.map(row => `<tr>${row.map(cell => `<td>${cell}</td>`).join('')}</tr>`).join('')}</tbody>`;
}

function openTab(evt, name) {
  document.querySelectorAll('.tab-content').forEach(t => t.classList.remove('active'));
  document.querySelectorAll('.tab-button').forEach(b => b.classList.remove('active'));
  const target = document.getElementById(name);
  if (target) target.classList.add('active');
  if (evt && evt.currentTarget) evt.currentTarget.classList.add('active');
}
window.openTab = openTab; // expose to inline onclick

/* =========================
   DOM Ready
   ========================= */
document.addEventListener('DOMContentLoaded', () => {
  // Late binding of UI refs
  ui = {
    materialSelect: document.getElementById('material-select'),
    thickInput:     document.getElementById('thick-input'),
    lxInput:        document.getElementById('lx-input'),
    lyInput:        document.getElementById('ly-input'),
    eInput:         document.getElementById('e-input'),
    rhoInput:       document.getElementById('rho-input'),
    nuInput:        document.getElementById('nu-input'),
    baffleCond:     document.getElementById('baffle-cond'),
    posX:           document.getElementById('pos-x'),
    posY:           document.getElementById('pos-y'),
    basicTable:     document.getElementById('basic-results-table'),
    infImpedanceTable: document.getElementById('inf-impedance-table'),
    impedanceTable: document.getElementById('impedance-table'),
    stlTable:       document.getElementById('stl-table'),
    radTable:       document.getElementById('rad-table'),
    basicFormula:   document.getElementById('basic-formula'),
    impFormula:     document.getElementById('impedance-formula'),
    stlFormula:     document.getElementById('stl-formula'),
    radFormula:     document.getElementById('rad-formula'),
  };

  // Populate materials
  if (ui.materialSelect) {
    ui.materialSelect.innerHTML = '';
    Object.keys(MATERIAL_PROPERTIES).forEach(name => {
      const opt = document.createElement('option');
      opt.value = name;
      opt.textContent = name;
      ui.materialSelect.appendChild(opt);
    });
    ui.materialSelect.value = "Gypsum Board"; // default
    ui.materialSelect.addEventListener('change', onMaterialSelect);
    onMaterialSelect();
  }

  // First run
  calculateAllTabs();
});

/* =========================
   Input handling
   ========================= */
function onMaterialSelect() {
  if (!ui || !ui.materialSelect) return;
  const selected = ui.materialSelect.value;
  const props = MATERIAL_PROPERTIES[selected];
  if (!props) return;

  if (ui.eInput)   ui.eInput.value   = (props.E   ?? '').toExponential ? props.E.toExponential(2) : (props.E ?? '');
  if (ui.rhoInput) ui.rhoInput.value = props.rho ?? '';
  if (ui.nuInput)  ui.nuInput.value  = props.nu  ?? '';
  [ui.eInput, ui.rhoInput, ui.nuInput].forEach(el => { if (el) el.readOnly = false; });
}

function getCommonInputs() {
  const h   = parseFloat(ui.thickInput?.value) / 1000.0;
  const E   = parseFloat(ui.eInput?.value);
  const rho = parseFloat(ui.rhoInput?.value);
  const nu  = parseFloat(ui.nuInput?.value);

  if (!(h > 0) || !(E > 0) || !(rho > 0)) throw new Error("Thickness, Young's Modulus, and Density must be positive.");

  const Lx = parseFloat(ui.lxInput?.value) / 1000.0;
  const Ly = parseFloat(ui.lyInput?.value) / 1000.0;
  if (!(Lx > 0) || !(Ly > 0)) throw new Error("Plate length and width must be positive.");

  const c0 = 343.0;
  const D  = (E * Math.pow(h, 3)) / (12 * (1 - nu*nu));
  const fc = (c0*c0/(2*Math.PI*h)) * Math.sqrt(12 * rho * (1 - nu*nu) / E);

  return { h, E, rho, nu, Lx, Ly, c0, D, fc };
}

/* =========================
   Top-level driver
   ========================= */
function calculateAllTabs() {
  try {
    calculateBasic();
    calculateImpedance();
    calculateSTL();
    calculateRadiation();
  } catch (e) {
    alert(e.message);
  }
}
window.calculateBasic     = calculateBasic;
window.calculateImpedance = calculateImpedance;
window.calculateSTL       = calculateSTL;
window.calculateRadiation = calculateRadiation;

/* =========================
   Basic
   ========================= */
function calculateBasic() {
  try {
    const { h, E, rho, nu, D, fc, Lx, Ly } = getCommonInputs();
    const C_L_p = Math.sqrt(E / (rho * (1 - nu*nu)));
    const S = Lx * Ly;
    const mode_density = (S * Math.sqrt(3)) / (h * C_L_p); // SS plate asymptotic

    createTable(ui.basicTable, ["Property", "Value", "Unit"], [
      ["Bending Stiffness (D)", D.toExponential(2), "N·m"],
      ["Critical Frequency (Fc)", fc.toFixed(2), "Hz"],
      ["Longitudinal Plate Wave Speed (C_L,p)", C_L_p.toFixed(2), "m/s"],
      ["Mode Density (n(f))", mode_density.toFixed(3), "modes/Hz"]
    ]);

    if (ui.basicFormula && window.MathJax) {
      ui.basicFormula.innerHTML = `<h5>Equations for Basic Properties</h5>
        <p><b>Bending Stiffness (D):</b> $$ D = \\frac{E h^3}{12(1-\\nu^2)} $$</p>
        <p><b>Longitudinal Plate Wave Speed (c_{L,p}):</b> $$ c_{L,p} = \\sqrt{\\frac{E}{\\rho(1-\\nu^2)}} $$</p>
        <p><b>Critical Frequency (f_c):</b> $$ f_c = \\frac{c_0^2 \\sqrt{3}}{\\pi h c_{L,p}} $$</p>
        <p><b>Mode Density (n(f)):</b> $$ n(f) = \\frac{S \\sqrt{3}}{h c_{L,p}} $$</p>`;
      MathJax.typeset();
    }
  } catch (e) { alert(e.message); }
}

/* =========================
   Helper: SS f_pq for p,q=1..3
   ========================= */
function ssModalFreqs_1to3(E,rho,nu,h,Lx,Ly){
  const rho_s = rho*h;
  const D = (E*h**3)/(12*(1-nu**2));
  const out = [];
  for (let p=1;p<=3;p++){
    for (let q=1;q<=3;q++){
      const f = (Math.PI/2) * Math.sqrt(D/rho_s) * ((p/Lx)**2 + (q/Ly)**2);
      out.push({p,q,f});
    }
  }
  out.sort((a,b)=>a.f-b.f);
  return out;
}

/* =========================
   Impedance (SS) — with requests applied
   ========================= */
function calculateImpedance() {
  try {
    const { E, rho, h, D, Lx, Ly, nu, fc } = getCommonInputs();
    const posX = parseFloat(ui.posX?.value) / 1000.0;
    const posY = parseFloat(ui.posY?.value) / 1000.0;
    const S = Lx * Ly;
    const rho_s = rho * h;

    // ref Z_inf and dB
    const C_L = Math.sqrt(E/rho);
    const Z_inf_ref = 2.3 * rho * C_L * h * h;
    const Z_INF_DB = (Z_inf_ref>0 && isFinite(Z_inf_ref)) ? 20*Math.log10(Z_inf_ref) : null;

    // ==== Property table: dB only ====
    createTable(ui.infImpedanceTable, ["Property", "Value"], [
      ["Infinite Plate Impedance (Z_inf, ref) [dB]", (Z_INF_DB==null?'NaN':Z_INF_DB.toFixed(2) + " dB re 1 Ns/m³")]
    ]);

    // ==== Impedance calculation up to 1000 Hz only ====
    const F_IMP = FREQS.filter(f => f <= 1000);

    const results = { mag_dB: [], real_dB: [], phase: [] };
    const tableData = [];
    for (const f of F_IMP) {
      const omega = 2 * Math.PI * f;
      const eta_f = 0.005 + 0.3 / Math.sqrt(Math.max(1e-6, f));

      let Y_sum = new Complex(0,0);
      const Pmax = 8, Qmax = 8;
      for (let p = 1; p <= Pmax; p++) {
        for (let q = 1; q <= Qmax; q++) {
          const f_pq = (Math.PI / 2) * Math.sqrt(D / rho_s) * ((p / Lx)**2 + (q / Ly)**2);
          const omega_pq = 2 * Math.PI * f_pq;
          const psi = Math.sin(p * Math.PI * posX / Lx) * Math.sin(q * Math.PI * posY / Ly);
          const denom = (new Complex(omega_pq**2, omega_pq**2 * eta_f)).sub(new Complex(omega**2, 0));
          Y_sum = Y_sum.add((new Complex(psi * psi, 0)).div(denom));
        }
      }
      const Y_dp = Complex.i().mul(new Complex(4 * omega / (rho_s * S), 0)).mul(Y_sum);

      // Z = 1/Y
      const d = Y_dp.re*Y_dp.re + Y_dp.im*Y_dp.im;
      let mag_dB = null, reZ_dB = null, ph = null;
      if (d > 0 && isFinite(d)) {
        const Zr =  Y_dp.re / d;
        const Zi = -Y_dp.im / d;
        const mag = Math.hypot(Zr, Zi);
        mag_dB = 20 * Math.log10(Math.max(mag, 1e-12));
        // Re(Z) を dB 表示（絶対値）
        reZ_dB = 20 * Math.log10(Math.max(Math.abs(Zr), 1e-12));
        ph = Math.atan2(Zi, Zr) * 180 / Math.PI;
      }

      results.mag_dB.push(mag_dB);
      results.real_dB.push(reZ_dB);
      results.phase.push(ph);

      tableData.push([
        f,
        mag_dB !== null ? mag_dB.toFixed(2) : 'NaN',
        reZ_dB !== null ? reZ_dB.toFixed(2) : 'NaN',
        ph !== null ? ph.toFixed(1) : 'NaN'
      ]);
    }

    // datasets
    const datasets = [
      { label: '|Z| [dB re 1 Ns/m³]', data: results.mag_dB, yAxisID: 'y',  borderWidth:2, pointRadius:2, borderColor: 'rgba(255,99,132,1)' },
      { label: 'Re(Z) [dB]',          data: results.real_dB, yAxisID: 'y1', borderWidth:2, pointRadius:2, borderColor: 'rgba(54,162,235,1)' },
      { label: 'Phase(Z) [deg]',      data: results.phase,   yAxisID: 'y2', borderWidth:2, pointRadius:2, borderColor: 'rgba(255,206,86,1)' }
    ];

    if (Z_INF_DB !== null) {
      datasets.push({
        label: 'Z_inf_ref [dB]',
        data: F_IMP.map(() => Z_INF_DB),
        yAxisID: 'y',
        borderColor: 'rgba(0,160,0,0.9)',
        borderDash: [6,3],
        fill: false,
        borderWidth: 2,
        pointRadius: 0
      });

      const modalSmall = ssModalFreqs_1to3(E,rho,nu,h,Lx,Ly)
        .filter(m => m.f <= 1000); // 1000Hzまでに制限
      const modalDots = modalSmall.map(m => ({ x: m.f, y: Z_INF_DB, _mn: `f(${m.p},${m.q})=${m.f.toFixed(1)} Hz` }));
      datasets.push({
        type: 'scatter',
        label: 'f11..f33 on Z_inf',
        data: modalDots,
        yAxisID: 'y',
        pointRadius: 4,
        pointHoverRadius: 6,
        showLine: false
      });
    }

    // draw chart
    const canvas = document.getElementById('impedance-chart');
    if (charts['impedance-chart']) charts['impedance-chart'].destroy();
    charts['impedance-chart'] = new Chart(canvas.getContext('2d'), {
      type: 'line',
      data: { labels: F_IMP, datasets },
      options: {
        responsive: true, maintainAspectRatio: false,
        interaction: { mode: 'index', intersect: false },
        plugins: {
          title: { display: true, text: 'Driving-Point Impedance (SS, up to 1 kHz)' },
          tooltip: {
            callbacks: {
              label: function(ctx){
                const ds = ctx.dataset;
                if (ds.type === 'scatter' && ds.label === 'f11..f33 on Z_inf' && ds.data && ds.data[ctx.dataIndex]) {
                  const pt = ds.data[ctx.dataIndex];
                  return `${ds.label}: ${pt._mn}`;
                }
                const def = Chart.defaults.plugins.tooltip.callbacks.label;
                return def ? def(ctx) : `${ds.label}: ${ctx.formattedValue}`;
              }
            }
          }
        },
        scales: {
          x: { type: 'logarithmic', title: { display: true, text: 'Frequency (Hz)' },
               ticks: { callback: v => [16,20,25,31.5,40,50,63,80,100,125,160,200,250,315,400,500,630,800,1000].includes(v) ? v : null },
               min: 16, max: 1000 },
          y:  { position: 'left',  title: { display: true, text: '|Z| [dB]' } },
          y1: { position: 'right', title: { display: true, text: 'Re(Z) [dB]' }, grid: { drawOnChartArea: false } },
          y2: { position: 'right', title: { display: true, text: 'Phase [deg]' }, grid: { drawOnChartArea: false }, offset: true }
        }
      }
    });

    // table under the figure (Re(Z) in dB)
    createTable(ui.impedanceTable, ['Frequency (Hz)', '|Z| dB', 'Re(Z) dB', 'Phase (deg)'], tableData);

    if (ui.impFormula && window.MathJax) {
      ui.impFormula.innerHTML = `<h5>Driving-Point Impedance (SS, modal sum)</h5>
        <p><b>Mobility:</b> $$ Y_{dp} = i\\omega \\frac{4}{\\rho_s S} \\sum_{p=1}^P \\sum_{q=1}^Q
        \\frac{\\psi_{p,q}^2(x,y)}{\\omega_{p,q}^2(1+i\\eta) - \\omega^2} $$</p>
        <p><b>Impedance:</b> $$ Z_{dp} = \\frac{1}{Y_{dp}} $$</p>
        <p><b>表示:</b> \( |Z|, \\ \\Re\\{Z\\} \) はともに \(20\\log_{10}(\\cdot)\\) の dB で描画。</p>`;
      MathJax.typeset();
    }

  } catch (e) { alert(e.message); }
}

/* =========================
   STL（そのまま）
   ========================= */
function calculateSTL() {
  try {
    const { rho, h, fc } = getCommonInputs();
    const eta = 0.01;
    const m_dp = rho * h;
    const R_mass = FREQS.map(f => 20 * Math.log10(m_dp * f) - 42.5);
    const R_coincidence = R_mass.map((r, i) => {
      const f = FREQS[i], fr = f / fc;
      if (Math.abs(fr - 1) < 1e-6) return r;
      const c = (2 * eta) / (Math.PI * fr) * (1 / Math.pow(Math.abs(1 - fr * fr), 2));
      return r - 10 * Math.log10(Math.abs(c) + 1);
    });

    const canvas = document.getElementById('stl-chart');
    if (charts['stl-chart']) charts['stl-chart'].destroy();
    charts['stl-chart'] = new Chart(canvas.getContext('2d'), {
      type: 'line',
      data: { labels: FREQS, datasets: [
        { label: 'Mass Law (Normal Incidence)', data: R_mass, borderWidth:2, pointRadius:2 },
        { label: 'With Coincidence Dip', data: R_coincidence, borderWidth:2, pointRadius:2 }
      ]},
      options: {
        responsive: true, maintainAspectRatio: false,
        plugins: { title: { display: true, text: 'Sound Transmission Loss (STL)' } },
        scales: {
          x: { type: 'logarithmic', title:{display:true, text:'Frequency (Hz)'},
               ticks: { callback: v => THIRD_OCTAVE_TICKS.includes(v) ? v : null },
               min: Math.min(...FREQS), max: Math.max(...FREQS) },
          y: { position: 'left', title: { display: true, text: 'STL (dB)' } }
        }
      }
    });

    createTable(ui.stlTable, ['Frequency (Hz)', 'Mass Law (dB)', 'With Coincidence (dB)'],
      FREQS.map((f, i) => [f, R_mass[i].toFixed(1), R_coincidence[i].toFixed(1)]));

    if (ui.stlFormula && window.MathJax) {
      ui.stlFormula.innerHTML =
        `<h5>Sound Transmission Loss (Single Plate)</h5>
         <p><b>Surface Density:</b> $$ m'' = \\rho h $$</p>
         <p><b>Mass Law (Normal Incidence):</b> $$ TL_0 = 20\\log_{10}(m'' f) - 42.5 \\; (dB) $$</p>
         <p>Coincidence correction is applied around \( f_c \).</p>`;
      MathJax.typeset();
    }

  } catch (e) { alert(e.message); }
}

/* =========================
   Radiation（そのまま）
   ========================= */
function calculateRadiation() {
  try {
    const { h, Lx, Ly, c0, fc } = getCommonInputs();
    const C_BC = 1; // simply supported
    const C_OB = parseFloat(ui.baffleCond?.value);
    const l = 2 * (Lx + Ly), S = Lx * Ly, lambda_c = c0 / fc, L1 = Math.min(Lx, Ly), L2 = Math.max(Lx, Ly), k_fc = 2 * Math.PI * fc / c0;

    const sigma_simple = FREQS.map(f => {
      let s;
      if (f > fc) s = 1.0;
      else if (Math.abs(f - fc) < 1e-6) s = 0.45 * Math.sqrt(l / lambda_c);
      else s = (l * lambda_c / (Math.PI ** 2 * S)) * Math.sqrt(f / fc);
      return Math.max(1e-4, Math.min(s, 1.0));
    });

    const sigma_leppington = FREQS.map(f => {
      let s;
      const mu = Math.sqrt(fc / f);
      if (f > fc) s = 1 / Math.sqrt(1 - fc / f);
      else if (Math.abs(f - fc) < 1e-6) s = (0.5 - 0.15 * L1 / L2) * Math.sqrt(k_fc * Math.sqrt(L1));
      else {
        const k = 2 * Math.PI * f / c0;
        const mu2_1 = mu ** 2 - 1;
        if (mu2_1 <= 0 || k <= 0) return 1e-4;
        const term1 = l / (2 * Math.PI * k * S * Math.sqrt(mu2_1));
        const term2 = 2 * Math.atanh(1 / mu);
        const term3 = 2 * mu / mu2_1;
        const term4 = C_BC * C_OB - Math.pow(mu, -8) * (C_BC * C_OB - 1);
        s = term1 * (term2 + term3) * term4;
      }
      return Math.max(1e-4, Math.min(s, 1.0));
    });

    const sigma_wallace = FREQS.map(f => {
      const k = 2 * Math.PI * f / c0;
      let sigma_sum = 0;
      for (let p = 1; p <= 5; p++) {
        for (let q = 1; q <= 5; q++) {
          let s_pq;
          const p_odd = p % 2 !== 0, q_odd = q % 2 !== 0;
          if (p_odd && q_odd) {
            s_pq = (32 * k**2 * Lx * Ly / (Math.PI**5 * p**2 * q**2)) * (1 - (k**2 * Lx * Ly / 12) * ((1 - 8 / (p * Math.PI)**2) * Lx / Ly + (1 - 8 / (q * Math.PI)**2) * Ly / Lx));
          } else if (p_odd !== q_odd) {
            if (p_odd && !q_odd) {
              s_pq = (8 * k**4 * Lx**3 * Ly / (3 * Math.PI**5 * p**2 * q**2)) * (1 - (k**2 * Lx * Ly / 20) * ((1 - 8 / (p * Math.PI)**2) * Lx / Ly + (1 - 24 / (q * Math.PI)**2) * Ly / Lx));
            } else {
              s_pq = (8 * k**4 * Ly**3 * Lx / (3 * Math.PI**5 * q**2 * p**2)) * (1 - (k**2 * Ly * Lx / 20) * ((1 - 8 / (q * Math.PI)**2) * Ly / Lx + (1 - 24 / (p * Math.PI)**2) * Lx / Ly));
            }
          } else {
            s_pq = (2 * k**6 * Lx**3 * Ly**3 / (15 * Math.PI**5 * p**2 * q**2)) * (1 - (5 * k**2 * Lx**2 / 64) * ((1 - 24 / (p * Math.PI)**2) * Lx / Ly + (1 - 24 / (q * Math.PI)**2) * Ly / Lx));
          }
          sigma_sum += s_pq > 0 ? s_pq : 0;
        }
      }
      return Math.max(1e-4, Math.min(sigma_sum, 1.0));
    });

    const canvas = document.getElementById('rad-chart');
    if (charts['rad-chart']) charts['rad-chart'].destroy();
    charts['rad-chart'] = new Chart(canvas.getContext('2d'), {
      type: 'line',
      data: { labels: FREQS, datasets: [
        { label: 'Simple Model', data: sigma_simple, borderWidth:2, pointRadius:2 },
        { label: 'Leppington', data: sigma_leppington, borderWidth:2, pointRadius:2 },
        { label: 'Wallace (Modes Sum)', data: sigma_wallace, borderWidth:2, pointRadius:2 }
      ]},
      options: {
        responsive: true, maintainAspectRatio: false,
        plugins: { title: { display: true, text: 'Radiation Efficiency (σ)' } },
        scales: {
          x: { type: 'logarithmic', title:{display:true, text:'Frequency (Hz)'},
               ticks: { callback: v => THIRD_OCTAVE_TICKS.includes(v) ? v : null },
               min: Math.min(...FREQS), max: Math.max(...FREQS) },
          y: { position: 'left', title: { display: true, text: 'Coefficient (σ)' }, min:0, max:1.1 }
        }
      }
    });

    createTable(ui.radTable, ['Frequency (Hz)', 'σ (Simple)', 'σ (Leppington)', 'σ (Wallace)'],
      FREQS.map((f, i) => [f, sigma_simple[i].toFixed(3), sigma_leppington[i].toFixed(3), sigma_wallace[i].toExponential(2)]));

    if (ui.radFormula && window.MathJax) {
      ui.radFormula.innerHTML =
        `<h5>Radiation Efficiency (σ) — models</h5>
         <p>Three approximate models are shown for comparison.</p>`;
      MathJax.typeset();
    }

  } catch (e) { alert(e.message); }
}
