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
  // Late binding of UI refs (確実に要素が存在してから取得)
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

  // Populate materials safely
  if (ui.materialSelect) {
    ui.materialSelect.innerHTML = ''; // clear just in case
    Object.keys(MATERIAL_PROPERTIES).forEach(name => {
      const opt = document.createElement('option');
      opt.value = name;
      opt.textContent = name;
      ui.materialSelect.appendChild(opt);
    });
    ui.materialSelect.value = "Gypsum Board"; // default
    ui.materialSelect.addEventListener('change', onMaterialSelect);
    onMaterialSelect(); // set fields
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
  if (!props) return; // safety

  // 値をセット（指数表記は parseFloat 可）
  if (ui.eInput)   ui.eInput.value   = (props.E   ?? '').toExponential ? props.E.toExponential(2) : (props.E ?? '');
  if (ui.rhoInput) ui.rhoInput.value = props.rho ?? '';
  if (ui.nuInput)  ui.nuInput.value  = props.nu  ?? '';

  // 入力可（Custom でも固定でもすべて手編集可能に）
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
    //console.error(e);
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
function ssModalFreq
