/* ---------- constants / UI ---------- */
const MATERIAL_PROPERTIES = {
  "Custom": { "rho": 0, "E": 0, "nu": 0 },
  "Concrete": { "rho": 2.3e3, "E": 2.1e10, "nu": 0.005 },
  "Gypsum Board": { "rho": 0.8e3, "E": 0.18e10, "nu": 0.005 },
  "Plywood": { "rho": 0.6e3, "E": 0.5e10, "nu": 0.30 },
  "Glass": { "rho": 2500, "E": 70.0e9, "nu": 0.23 }
};
const FREQS = [16,20,25,31.5,40,50,63,80,100,125,160,200,250,315,400,500,630,800,1000,1250,1600,2000,2500,3150,4000,5000,6300,8000];
const THIRD_OCTAVE_TICKS = [16, 31.5, 63, 125, 250, 500, 1000, 2000, 4000, 8000];
let charts = {};
const ui = {
  materialSelect: document.getElementById('material-select'),
  thickInput: document.getElementById('thick-input'),
  lxInput: document.getElementById('lx-input'),
  lyInput: document.getElementById('ly-input'),
  eInput: document.getElementById('e-input'),
  rhoInput: document.getElementById('rho-input'),
  nuInput: document.getElementById('nu-input'),
  baffleCond: document.getElementById('baffle-cond'),
  posX: document.getElementById('pos-x'),
  posY: document.getElementById('pos-y'),
  basicTable: document.getElementById('basic-results-table'),
  infImpedanceTable: document.getElementById('inf-impedance-table'),
  impedanceTable: document.getElementById('impedance-table'),
  stlTable: document.getElementById('stl-table'),
  radTable: document.getElementById('rad-table'),
};

/* ---------- small Complex class ---------- */
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

/* ---------- init ---------- */
document.addEventListener('DOMContentLoaded', () => {
  for (const name in MATERIAL_PROPERTIES) ui.materialSelect.add(new Option(name, name));
  ui.materialSelect.value = "Gypsum Board";
  ui.materialSelect.addEventListener('change', onMaterialSelect);
  onMaterialSelect();
  calculateAllTabs();
});

/* ---------- overall calculation driver ---------- */
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

/* ---------- input handling ---------- */
function getCommonInputs() {
  const h = parseFloat(ui.thickInput.value) / 1000.0;
  const E = parseFloat(ui.eInput.value);
  const rho = parseFloat(ui.rhoInput.value);
  const nu = parseFloat(ui.nuInput.value);
  if (!(h > 0) || !(E > 0) || !(rho > 0)) throw new Error("Thickness, Young's Modulus, and Density must be positive.");
  const Lx = parseFloat(ui.lxInput.value) / 1000.0;
  const Ly = parseFloat(ui.lyInput.value) / 1000.0;
  if (!(Lx > 0) || !(Ly > 0)) throw new Error("Plate length and width must be positive.");
  const c0 = 343.0;
  const D = (E * Math.pow(h, 3)) / (12 * (1 - nu*nu));
  const fc = (c0*c0/(2*Math.PI*h)) * Math.sqrt(12 * rho * (1 - nu*nu) / E);
  return { h, E, rho, nu, Lx, Ly, c0, D, fc };
}

/* ---------- Basic ---------- */
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
    const formula = document.getElementById('basic-formula');
    if (formula) {
      formula.innerHTML = `<h5>Equations for Basic Properties</h5>
        <p><b>Bending Stiffness (D):</b> $$ D = \\frac{E h^3}{12(1-\\nu^2)} $$</p>
        <p><b>Longitudinal Plate Wave Speed (c_{L,p}):</b> $$ c_{L,p} = \\sqrt{\\frac{E}{\\rho(1-\\nu^2)}} $$</p>
        <p><b>Critical Frequency (f_c):</b> $$ f_c = \\frac{c_0^2 \\sqrt{3}}{\\pi h c_{L,p}} $$</p>
        <p><b>Mode Density (n(f)):</b> $$ n(f) = \\frac{S \\sqrt{3}}{h c_{L,p}} $$</p>`;
      if (window.MathJax) MathJax.typeset();
    }
  } catch (e) { alert(e.message); }
}

/* ---------- SS modal frequencies f_pq for p,q=1..3 ---------- */
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

/* ---------- Impedance (SS, modal sum) ---------- */
function calculateImpedance() {
  try {
    const { E, rho, h, D, Lx, Ly, nu, fc } = getCommonInputs();
    const posX = parseFloat(ui.posX.value) / 1000.0;
    const posY = parseFloat(ui.posY.value) / 1000.0;
    const S = Lx * Ly;
    const rho_s = rho * h;

    // infinite-plate reference impedance
    const C_L = Math.sqrt(E/rho);
    const Z_inf_ref = 2.3 * rho * C_L * h * h;
    const Z_inf_ref_dB = (Z_inf_ref>0 && isFinite(Z_inf_ref)) ? 20*Math.log10(Z_inf_ref) : null;

    // show both linear and dB in the table
    createTable(ui.infImpedanceTable, ["Property", "Value"], [
      ["Infinite Plate Impedance (Z_inf, ref)", `${Z_inf_ref.toExponential(2)} Ns/m³`],
      ["Infinite Plate Impedance (Z_inf, ref) [dB]", (Z_inf_ref_dB==null?'NaN':Z_inf_ref_dB.toFixed(2) + " dB re 1 Ns/m³")]
    ]);

    // modal sum of mobility (SS)
    const results = { mag_dB: [], real: [], phase: [] };
    const tableData = [];
    for (const f of FREQS) {
      const omega = 2 * Math.PI * f;
      const eta_f = 0.005 + 0.3 / Math.sqrt(Math.max(1e-6, f)); // freq-dependent loss

      let Y_sum = new Complex(0,0);
      const Pmax = 8, Qmax = 8; // truncation
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
      const d = Y_dp.re*Y_dp.re + Y_dp.im*Y_dp.im;

      let mag_dB = null, reZ = null, ph = null;
      if (d > 0 && isFinite(d)) {
        const Zr =  Y_dp.re / d;
        const Zi = -Y_dp.im / d;
        const mag = Math.hypot(Zr, Zi);
        mag_dB = 20 * Math.log10(Math.max(mag, 1e-12));
        reZ = Zr;
        ph = Math.atan2(Zi, Zr) * 180 / Math.PI;
      }

      results.mag_dB.push(mag_dB);
      results.real.push(reZ);
      results.phase.push(ph);

      tableData.push([
        f,
        mag_dB !== null ? mag_dB.toFixed(2) : 'NaN',
        reZ !== null ? reZ.toExponential(2) : 'NaN',
        ph !== null ? ph.toFixed(1) : 'NaN'
      ]);
    }

    // datasets: |Z| dB, Re(Z), Phase, + Z_inf_ref horizontal + modal dots
    const datasets = [
      { label: '|Z| [dB re 1 Ns/m³]', data: results.mag_dB, yAxisID: 'y', borderColor: 'rgba(255,99,132,1)' },
      { label: 'Re(Z) [Ns/m³]',        data: results.real,   yAxisID: 'y1'
