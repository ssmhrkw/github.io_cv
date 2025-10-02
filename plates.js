/* =========================
   Constants / UI handles
========================= */
const MATERIAL_PROPERTIES = {
  "Custom":       { rho: 0,     E: 0,       nu: 0    },
  "Concrete":     { rho: 2.3e3, E: 2.1e10,  nu: 0.005},
  "Gypsum Board": { rho: 0.8e3, E: 0.18e10, nu: 0.005},
  "Plywood":      { rho: 0.6e3, E: 0.5e10,  nu: 0.30 },
  "Glass":        { rho: 2500,  E: 70.0e9,  nu: 0.23 }
};

const FREQS = [16,20,25,31.5,40,50,63,80,100,125,160,200,250,315,400,500,630,800,1000,1250,1600,2000,2500,3150,4000,5000,6300,8000];

let charts = {};
const ui = {
  materialSelect:    document.getElementById('material-select'),
  thickInput:        document.getElementById('thick-input'),
  lxInput:           document.getElementById('lx-input'),
  lyInput:           document.getElementById('ly-input'),
  eInput:            document.getElementById('e-input'),
  rhoInput:          document.getElementById('rho-input'),
  nuInput:           document.getElementById('nu-input'),
  baffleCond:        document.getElementById('baffle-cond'),
  posX:              document.getElementById('pos-x'),
  posY:              document.getElementById('pos-y'),
  basicTable:        document.getElementById('basic-results-table'),
  naturalFreqTable:  document.getElementById('natural-freq-table'),
  infImpedanceTable: document.getElementById('inf-impedance-table'),
  impedanceTable:    document.getElementById('impedance-table'),
  stlTable:          document.getElementById('stl-table'),
  radTable:          document.getElementById('rad-table'),
};

/* =========================
   Complex (minimal)
========================= */
class Complex {
  constructor(re=0, im=0){ this.re = re; this.im = im; }
  static inv(z){
    const d = z.re*z.re + z.im*z.im;
    if (d === 0) return { re: Infinity, im: Infinity };
    return { re:  z.re/d, im: -z.im/d };
  }
  static mag(z){ return Math.hypot(z.re, z.im); }
  static phase(z){ return Math.atan2(z.im, z.re); }
}

/* =========================
   Utilities
========================= */
const clamp = (x, lo, hi) => Math.min(Math.max(x, lo), hi);
function lossFactor(f){ // η(f) with clamp for stability
  const eta = 0.005 + 0.3/Math.sqrt(Math.max(1, f));
  return clamp(eta, 1e-4, 0.2);
}
function mag2dB(x){
  const v = Math.max(x, 1e-12);
  return 20*Math.log10(v);
}
function toMeters(mm){ return mm/1000; }

/* =========================
   Input / derived quantities
========================= */
function getCommonInputs(){
  const h   = toMeters(parseFloat(ui.thickInput.value));
  const E   = parseFloat(ui.eInput.value);
  const rho = parseFloat(ui.rhoInput.value);
  const nu  = parseFloat(ui.nuInput.value);
  if (!(h>0) || !(E>0) || !(rho>0)) throw new Error("Thickness, Young's Modulus, and Density must be positive.");

  const Lx  = toMeters(parseFloat(ui.lxInput.value));
  const Ly  = toMeters(parseFloat(ui.lyInput.value));
  if (!(Lx>0) || !(Ly>0)) throw new Error("Plate length and width must be positive.");

  const c0 = 343.0;
  const D  = (E*h**3)/(12*(1-nu**2));
  const fc = (c0*c0/(2*Math.PI*h)) * Math.sqrt(12*rho*(1-nu**2)/E);

  return { h, E, rho, nu, Lx, Ly, c0, D, fc };
}

/* =========================
   UI helpers
========================= */
document.addEventListener('DOMContentLoaded', () => {
  for (const name in MATERIAL_PROPERTIES) ui.materialSelect.add(new Option(name, name));
  ui.materialSelect.value = "Gypsum Board";
  ui.materialSelect.addEventListener('change', onMaterialSelect);
  onMaterialSelect();
  calculateAllTabs();
});

function onMaterialSelect(){
  const selected = ui.materialSelect.value;
  const props = MATERIAL_PROPERTIES[selected];
  const isCustom = selected === "Custom";
  ui.eInput.value   = isCustom ? "" : props.E.toExponential(2);
  ui.rhoInput.value = isCustom ? "" : props.rho;
  ui.nuInput.value  = isCustom ? "" : props.nu;
  [ui.eInput, ui.rhoInput, ui.nuInput].forEach(i => i.readOnly = false);
}

function openTab(evt, tabName){
  document.querySelectorAll('.tab-content').forEach(tab => tab.classList.remove('active'));
  document.querySelectorAll('.tab-button').forEach(btn => btn.classList.remove('active'));
  document.getElementById(tabName).classList.add('active');
  evt.currentTarget.classList.add('active');
}

function createTable(container, headers, data){
  if (!container) return;
  container.innerHTML =
    `<table><thead><tr>${headers.map(h=>`<th>${h}</th>`).join('')}</tr></thead>
      <tbody>${data.map(row=>`<tr>${row.map(cell=>`<td>${cell}</td>`).join('')}</tr>`).join('')}</tbody></table>`;
}

/* =========================
   Plot (with 1/3-oct bands, fc line, zinf)
========================= */
function ThirdOctavePlugin(opts={}){
  const { showBands=true, fc=null, zInfDb=null, yAxisId='y' } = opts;
  const pow6 = Math.pow(2, 1/6);

  return {
    id: 'thirdOctavePlugin',
    afterDraw(chart){
      const {ctx, chartArea:area, scales} = chart;
      const xScale = scales.x;
      const yScale = scales[yAxisId];
      if (!xScale || !area) return;

      const labels = chart.data.labels.map(Number);

      // shaded 1/3-octave bands
      if (showBands){
        for (let i=0;i<labels.length;i++){
          const fc0 = labels[i];
          const fl = fc0 / pow6, fu = fc0 * pow6;
          const x0 = xScale.getPixelForValue(fl);
          const x1 = xScale.getPixelForValue(fu);
          if (x1 <= area.left || x0 >= area.right) continue;
          ctx.save();
          ctx.fillStyle = (i%2===0) ? 'rgba(200,200,200,0.05)' : 'rgba(180,180,180,0.03)';
          const left  = Math.max(x0, area.left);
          const right = Math.min(x1, area.right);
          ctx.fillRect(left, area.top, right-left, area.bottom-area.top);
          ctx.restore();
        }
      }

      // custom x labels at every center frequency
      ctx.save();
      ctx.font = '11px Arial';
      ctx.fillStyle = '#222';
      ctx.textAlign = 'center';
      const yText = area.bottom + 14;
      for (const f of labels){
        const x = xScale.getPixelForValue(f);
        if (x < area.left-10 || x > area.right+10) continue;
        ctx.fillText(String(f), x, yText);
      }
      ctx.restore();

      // vertical line at fc
      if (fc){
        const xfc = xScale.getPixelForValue(fc);
        if (xfc >= area.left && xfc <= area.right){
          ctx.save();
          ctx.strokeStyle = 'rgba(60,60,60,0.9)';
          ctx.setLineDash([4,4]);
          ctx.lineWidth = 1;
          ctx.beginPath();
          ctx.moveTo(xfc, area.top);
          ctx.lineTo(xfc, area.bottom);
          ctx.stroke();
          ctx.restore();
        }
      }

      // horizontal Z_inf_ref (dB)
      if (typeof zInfDb === 'number' && yScale){
        const y = yScale.getPixelForValue(zInfDb);
        if (y >= area.top && y <= area.bottom){
          ctx.save();
          ctx.strokeStyle = 'rgba(0,140,0,0.9)';
          ctx.setLineDash([6,3]);
          ctx.lineWidth = 1;
          ctx.beginPath();
          ctx.moveTo(area.left, y);
          ctx.lineTo(area.right, y);
          ctx.stroke();
          ctx.restore();
        }
      }
    }
  };
}

function createPlot(chartId, labels, datasets, title, yAxisLabel='dB', yAxisOptions={}, pluginOpts={}){
  const canvas = document.getElementById(chartId);
  if (!canvas) return;
  const ctx = canvas.getContext('2d');
  if (charts[chartId]) charts[chartId].destroy();

  const numericLabels = labels.map(Number);
  const scales = {
    x: {
      type: 'logarithmic',
      title: { display: true, text: 'Frequency (Hz)' },
      ticks: { display: false },
      min: Math.min(...numericLabels),
      max: Math.max(...numericLabels)
    },
    y:   { position: 'left',  title: { display: true, text: yAxisLabel }, ...yAxisOptions },
    y1:  { position: 'right', title: { display: true, text: 'Re(Z) [Ns/m³]' }, grid: { drawOnChartArea: false } },
    y2:  { position: 'right', title: { display: true, text: 'Phase [deg]' }, grid: { drawOnChartArea: false }, offset: true },
    yMd: { position: 'left',  title: { display: true, text: 'Modal Density [modes/Hz]' }, grid: { drawOnChartArea: false }, offset: true }
  };

  charts[chartId] = new Chart(ctx, {
    type: 'line',
    data: {
      labels: numericLabels,
      datasets: (datasets||[]).map(ds => ({
        ...ds,
        fill: false,
        spanGaps: true,
        borderWidth: 2,
        pointRadius: 2,
        data: (ds.data||[]).map(v => (v==null||!Number.isFinite(+v)) ? null : +v)
      }))
    },
    options: {
      responsive: true,
      maintainAspectRatio: false,
      interaction: { mode: 'index', intersect: false },
      plugins: { title: { display: true, text: title }, legend: { position: 'top' } },
      scales
    },
    plugins: [ ThirdOctavePlugin({ ...pluginOpts, yAxisId: 'y' }) ]
  });
}

/* =========================
   PlateModel with cache
========================= */
function buildModalCache({Lx, Ly, h, E, rho, nu, px, py, P=12, Q=12}){
  const rho_s = rho*h;
  const D = (E*h**3)/(12*(1-nu**2));
  const S = Lx*Ly;

  const N = P*Q;
  const psi2 = new Float64Array(N);
  const w2   = new Float64Array(N);

  let k=0;
  for (let p=1; p<=P; p++){
    const sp = Math.sin(p*Math.PI*px/Lx);
    const a  = (p/Lx)**2;
    for (let q=1; q<=Q; q++){
      const sq = Math.sin(q*Math.PI*py/Ly);
      psi2[k] = (sp*sq)*(sp*sq);
      const b = (q/Ly)**2;
      const f_pq = (Math.PI/2) * Math.sqrt(D/rho_s) * (a + b);
      w2[k] = (2*Math.PI*f_pq)**2;
      k++;
    }
  }
  // 昇順に並べ替え（任意：トランケーションに使える）
  const idx = [...Array(N).keys()].sort((i,j)=>w2[i]-w2[j]);
  const psi2_sorted = new Float64Array(N);
  const w2_sorted   = new Float64Array(N);
  for (let i=0;i<N;i++){ psi2_sorted[i] = psi2[idx[i]]; w2_sorted[i] = w2[idx[i]]; }

  return { rho_s, S, psi2: psi2_sorted, w2: w2_sorted, P, Q };
}

function mobilityFromCache(omega, eta, cache, eps=1e-10){
  const { psi2, w2, rho_s, S } = cache;
  let Yr=0, Yi=0; // sum (psi^2 / denom) * (Ar - i Ai)
  // トランケーション：寄与が微小になったらbreak（緩い停止条件）
  let accum=0;
  for (let k=0;k<psi2.length;k++){
    const w2k = w2[k];
    const Ar = (w2k - omega*omega);
    const Ai = (w2k * eta);
    const denom = Ar*Ar + Ai*Ai;
    const scale = psi2[k] / denom;
    const termAbs = scale * Math.hypot(Ar, Ai);
    accum += termAbs;
    Yr += scale * Ar;
    Yi -= scale * Ai;
    if (termAbs < eps && accum > 100*eps) break;
  }
  const C = 4/(rho_s*S);
  const Ydp_r = - omega*C*Yi; // Re{i*ω*C*(Yr+iYi)} = -ω*C*Yi
  const Ydp_i =   omega*C*Yr; // Im{...}           =  ω*C*Yr
  return { re: Ydp_r, im: Ydp_i };
}

class PlateModel {
  constructor(params){ this.setParams(params); }
  setParams({Lx,Ly,h,E,rho,nu,px,py,P=12,Q=12}){
    Object.assign(this, {Lx,Ly,h,E,rho,nu,px,py,P,Q});
    // posが端で sin=0 になり過ぎるのを避ける微小オフセット
    const eps = 1e-9;
    this.px = clamp(px, eps, Lx-eps);
    this.py = clamp(py, eps, Ly-eps);
    this.cache = buildModalCache(this);
  }
  impedanceAt(f){
    const omega = 2*Math.PI*f;
    const eta = lossFactor(f);
    const Y = mobilityFromCache(omega, eta, this.cache);
    const d = Y.re*Y.re + Y.im*Y.im;
    if (d === 0) return { Zr: null, Zi: null, ZdB: null, phaseDeg: null };
    const Zr =  Y.re/d, Zi = -Y.im/d;
    return {
      Zr, Zi,
      ZdB: mag2dB(Math.hypot(Zr, Zi)),
      phaseDeg: Math.atan2(Zi, Zr)*180/Math.PI
    };
  }
}

/* =========================
   Calculations
========================= */
function calculateAllTabs(){
  try{
    calculateBasic();
    calculateNaturalFreq();
    calculateImpedance();
    calculateSTL();
    calculateRadiation();
  }catch(e){ alert(e.message); }
}

/* --- Basic --- */
function calculateBasic(){
  try{
    const { h, E, rho, nu, D, fc, Lx, Ly } = getCommonInputs();
    const C_L_p = Math.sqrt(E/(rho*(1-nu**2)));
    const S = Lx*Ly;
    const mode_density = (S*Math.sqrt(3))/(h*C_L_p);
    createTable(ui.basicTable,
      ["Property","Value","Unit"],
      [
        ["Bending Stiffness (D)", D.toExponential(2), "N·m"],
        ["Critical Frequency (Fc)", fc.toFixed(2), "Hz"],
        ["Longitudinal Plate Wave Speed (C_L,p)", C_L_p.toFixed(2), "m/s"],
        ["Mode Density (n(f))", mode_density.toFixed(3), "modes/Hz"]
      ]
    );
    const elem = document.getElementById('basic-formula');
    if (elem){
      elem.innerHTML =
        `<h5>Equations for Basic Properties</h5>
         <p><b>Bending Stiffness (D):</b> $$ D = \\frac{E h^3}{12(1-\\nu^2)} $$</p>
         <p><b>Longitudinal Plate Wave Speed (c_{L,p}):</b> $$ c_{L,p} = \\sqrt{\\frac{E}{\\rho(1-\\nu^2)}} $$</p>
         <p><b>Critical Frequency (f_c):</b> $$ f_c = \\frac{c_0^2 \\sqrt{3}}{\\pi h c_{L,p}} $$</p>
         <p><b>Mode Density (n(f)):</b> $$ n(f) = \\frac{S \\sqrt{3}}{h c_{L,p}} $$</p>`;
      if (window.MathJax) MathJax.typeset();
    }
  }catch(e){ alert(e.message); }
}

/* --- Natural Frequencies --- */
function calculateNaturalFreq(){
  try{
    const { rho, h, D, Lx, Ly } = getCommonInputs();
    const rows = [];
    const nf = document.getElementById('natural-freq-formula');
    if (nf){
      nf.innerHTML =
        `<h5>Natural Frequencies (Simply Supported)</h5>
         <p>$$ f_{m,n} = \\frac{\\pi}{2} \\sqrt{\\frac{D}{\\rho h}}
         \\left( \\left(\\frac{m}{L_x}\\right)^2 + \\left(\\frac{n}{L_y}\\right)^2 \\right) $$</p>`;
      if (window.MathJax) MathJax.typeset();
    }
    for (let m=1;m<=5;m++){
      for (let n=1;n<=5;n++){
        const fmn = (Math.PI/2) * Math.sqrt(D/(rho*h)) * ((m/Lx)**2 + (n/Ly)**2);
        rows.push([`f(${m},${n})`, fmn.toFixed(2)]);
      }
    }
    createTable(ui.naturalFreqTable, ["Mode (m,n)","Frequency (Hz)"], rows);
  }catch(e){ alert(e.message); }
}

/* --- Impedance (fast & stable) --- */
function calculateImpedance(){
  try{
    const { E, rho, h, D, Lx, Ly, fc } = getCommonInputs();
    const px = toMeters(parseFloat(ui.posX.value));
    const py = toMeters(parseFloat(ui.posY.value));

    // Infinite-plate ref (as-is)
    const C_L = Math.sqrt(E/rho);
    const Z_inf_ref = 2.3 * rho * C_L * h*h;
    const Z_inf_ref_dB = Number.isFinite(Z_inf_ref) && Z_inf_ref>0 ? 20*Math.log10(Z_inf_ref) : null;
    createTable(ui.infImpedanceTable, ["Property","Value"],
      [["Infinite Plate Impedance (Z_inf, ref)", `${Z_inf_ref.toExponential(2)} Ns/m³`]]
    );

    // Build model/cache once (12x12 modes; adjust if needed)
    const model = new PlateModel({ Lx, Ly, h, E, rho, nu: parseFloat(ui.nuInput.value||0.3), px, py, P:12, Q:12 });

    const mag_dB = [], realZ = [], phase = [], modal_density = [];
    const S = Lx*Ly, C_L_p = Math.sqrt(E/(rho*(1-(parseFloat(ui.nuInput.value||0.3))**2)));
    const md = (S*Math.sqrt(3))/(h*C_L_p);

    const tableRows = [];
    for (const f of FREQS){
      const { Zr, Zi, ZdB, phaseDeg } = model.impedanceAt(f);
      mag_dB.push(ZdB);
      realZ.push(Zr);
      phase.push(phaseDeg);
      modal_density.push(md);
      tableRows.push([
        f,
        ZdB!=null ? ZdB.toFixed(2) : 'NaN',
        Zr!=null  ? Zr.toExponential(2) : 'NaN',
        phaseDeg!=null ? phaseDeg.toFixed(1) : 'NaN'
      ]);
    }

    const datasets = [
      { label: '|Z| [dB re 1 Ns/m³]', data: mag_dB, yAxisID:'y',   borderColor:'rgba(255,99,132,1)' },
      { label: 'Re(Z) [Ns/m³]',       data: realZ,  yAxisID:'y1',  borderColor:'rgba(54,162,235,1)' },
      { label: 'Phase(Z) [deg]',      data: phase,  yAxisID:'y2',  borderColor:'rgba(255,206,86,1)' },
      { label: 'Modal density n(f)',  data: modal_density, yAxisID:'yMd', borderColor:'rgba(0,120,120,0.85)', borderDash:[5,3] }
    ];
    if (Z_inf_ref_dB != null){
      datasets.push({ label:'Z_inf_ref [dB]', data: FREQS.map(()=>Z_inf_ref_dB), yAxisID:'y', borderColor:'rgba(0,160,0,0.9)', borderDash:[6,3] });
    }

    createPlot('impedance-chart', FREQS, datasets,
               'Driving-Point Impedance (point)', '|Z| [dB]',
               {}, { showBands:true, fc, zInfDb: Z_inf_ref_dB });

    createTable(ui.impedanceTable, ['Frequency (Hz)','|Z| dB','Re(Z)','Phase (deg)'], tableRows);

    const impF = document.getElementById('impedance-formula');
    if (impF){
      impF.innerHTML =
        `<h5>Driving-Point Impedance (modal sum)</h5>
         <p><b>Mobility:</b> $$ Y_{dp}(x,y) = i\\omega \\frac{4}{\\rho_s S}
         \\sum_{p=1}^{\\infty}\\sum_{q=1}^{\\infty}
         \\frac{\\psi_{pq}^2(x,y)}{\\omega_{pq}^2(1+i\\eta)-\\omega^2} $$</p>
         <p><b>Impedance:</b> $$ Z_{dp} = \\frac{1}{Y_{dp}} $$</p>`;
      if (window.MathJax) MathJax.typeset();
    }
  }catch(e){ alert(e.message); }
}

/* --- STL --- */
function calculateSTL(){
  try{
    const { rho, h, fc } = getCommonInputs();
    const eta = 0.01;
    const m_dp = rho*h;

    const R_mass = FREQS.map(f => 20*Math.log10(m_dp*f) - 42.5);
    const R_coin = R_mass.map((r,i)=>{
      const f = FREQS[i], fr = f/fc;
      if (Math.abs(fr-1) < 1e-6) return r;
      const c = (2*eta)/(Math.PI*fr) * (1/Math.pow(Math.abs(1-fr*fr),2));
      return r - 10*Math.log10(Math.abs(c)+1);
    });

    createPlot('stl-chart', FREQS,
      [
        { label:'Mass Law (Normal Incidence)', data:R_mass, yAxisID:'y', borderColor:'rgba(100,100,255,1)' },
        { label:'With Coincidence Dip',        data:R_coin, yAxisID:'y', borderColor:'rgba(255,100,100,1)' }
      ],
      'Sound Transmission Loss (STL)', 'STL (dB)',
      {}, { showBands: true, fc }
    );

    createTable(ui.stlTable, ['Frequency (Hz)','Mass Law (dB)','With Coincidence (dB)'],
      FREQS.map((f,i)=>[f, R_mass[i].toFixed(1), R_coin[i].toFixed(1)])
    );

    const stlF = document.getElementById('stl-formula');
    if (stlF){
      stlF.innerHTML =
        `<h5>Sound Transmission Loss</h5>
         <p>Surface density: $$ m'' = \\rho h $$</p>
         <p>Mass law (normal incidence): $$ TL_0 = 20\\log_{10}(m'' f) - 42.5 \\;\\text{(dB)} $$</p>
         <p>Coincidence correction: simple dip term around \(f_c\).</p>`;
      if (window.MathJax) MathJax.typeset();
    }
  }catch(e){ alert(e.message); }
}

/* --- Radiation Efficiency --- */
function calculateRadiation(){
  try{
    const { h, Lx, Ly, c0, fc } = getCommonInputs();
    const C_BC = 1;
    const C_OB = parseFloat(ui.baffleCond.value);
    const l = 2*(Lx+Ly), S=Lx*Ly, lambda_c=c0/fc, L1=Math.min(Lx,Ly), L2=Math.max(Lx,Ly), k_fc=2*Math.PI*fc/c0;

    const sigma_simple = FREQS.map(f=>{
      let s;
      if (f>fc) s = 1.0;
      else if (Math.abs(f-fc)<1e-6) s = 0.45*Math.sqrt(l/lambda_c);
      else s = (l*lambda_c/(Math.PI**2*S)) * Math.sqrt(f/fc);
      return clamp(s, 1e-4, 1.0);
    });

    const sigma_leppington = FREQS.map(f=>{
      let s;
      const mu = Math.sqrt(fc/f);
      if (f>fc) s = 1/Math.sqrt(1 - fc/f);
      else if (Math.abs(f-fc)<1e-6) s = (0.5 - 0.15*L1/L2)*Math.sqrt(k_fc*Math.sqrt(L1));
      else {
        const k = 2*Math.PI*f/c0;
        const mu2_1 = mu*mu - 1;
        if (mu2_1<=0 || k<=0) return 1e-4;
        const term1 = l/(2*Math.PI*k*S*Math.sqrt(mu2_1));
        const term2 = 2*Math.atanh(1/mu);
        const term3 = 2*mu/mu2_1;
        const term4 = C_BC*C_OB - Math.pow(mu,-8)*(C_BC*C_OB - 1);
        s = term1*(term2 + term3)*term4;
      }
      return clamp(s, 1e-4, 1.0);
    });

    const sigma_wallace = FREQS.map(f=>{
      const k = 2*Math.PI*f/c0;
      let sigma_sum = 0;
      for (let p=1;p<=5;p++){
        for (let q=1;q<=5;q++){
          let s_pq;
          const p_odd = (p%2)!==0, q_odd = (q%2)!==0;
          if (p_odd && q_odd){
            s_pq = (32*k**2*Lx*Ly/(Math.PI**5*p**2*q**2))*(1-(k**2*Lx*Ly/12)*((1-8/(p*Math.PI)**2)*Lx/Ly+(1-8/(q*Math.PI)**2)*Ly/Lx));
          } else if (p_odd !== q_odd){
            if (p_odd && !q_odd){
              s_pq = (8*k**4*Lx**3*Ly/(3*Math.PI**5*p**2*q**2))*(1-(k**2*Lx*Ly/20)*((1-8/(p*Math.PI)**2)*Lx/Ly+(1-24/(q*Math.PI)**2)*Ly/Lx));
            } else {
              s_pq = (8*k**4*Ly**3*Lx/(3*Math.PI**5*q**2*p**2))*(1-(k**2*Ly*Lx/20)*((1-8/(q*Math.PI)**2)*Ly/Lx+(1-24/(p*Math.PI)**2)*Lx/Ly));
            }
          } else {
            s_pq = (2*k**6*Lx**3*Ly**3/(15*Math.PI**5*p**2*q**2))*(1-(5*k**2*Lx**2/64)*((1-24/(p*Math.PI)**2)*Lx/Ly+(1-24/(q*Math.PI)**2)*Ly/Lx));
          }
          sigma_sum += Math.max(s_pq, 0);
        }
      }
      return clamp(sigma_sum, 1e-4, 1.0);
    });

    createPlot('rad-chart', FREQS,
      [
        { label:'Simple Model',           data:sigma_simple,    borderColor:'rgba(60,160,60,1)' },
        { label:'Leppington',             data:sigma_leppington,borderColor:'rgba(160,60,160,1)' },
        { label:'Wallace (Modes Sum)',    data:sigma_wallace,   borderColor:'rgba(60,60,160,1)' }
      ],
      'Radiation Efficiency (σ)', 'Coefficient (σ)', { min:0, max:1.1 }, { showBands:true, fc }
    );

    createTable(ui.radTable, ['Frequency (Hz)','σ (Simple)','σ (Leppington)','σ (Wallace)'],
      FREQS.map((f,i)=>[f, sigma_simple[i].toFixed(3), sigma_leppington[i].toFixed(3), sigma_wallace[i].toExponential(2)])
    );

    const radF = document.getElementById('rad-formula');
    if (radF){
      radF.innerHTML = `<h5>Radiation Efficiency (σ)</h5>
        <p>Three models shown: simple subcritical/supercritical approx., Leppington, and Wallace modal sum.</p>`;
      if (window.MathJax) MathJax.typeset();
    }
  }catch(e){ alert(e.message); }
}
