/* ---------------- constants ---------------- */
const MATERIAL_PROPERTIES = {
  "Custom": { rho: 0, E: 0, nu: 0.30 },
  "Concrete": { rho: 2300, E: 2.1e10, nu: 0.2 },
  "Gypsum Board": { rho: 800, E: 0.18e10, nu: 0.25 },
  "Plywood": { rho: 600, E: 0.5e10, nu: 0.30 },
  "Glass": { rho: 2500, E: 70.0e9, nu: 0.23 }
};
const FREQS = [16,20,25,31.5,40,50,63,80,100,125,160,200,250,315,400,500,630,800,1000];
let charts = {};
let ui = {};

/* ---------------- small Complex ---------------- */
class Complex {
  constructor(re=0, im=0){ this.re = re; this.im = im; }
  add(z){ return new Complex(this.re+z.re, this.im+z.im); }
  sub(z){ return new Complex(this.re-z.re, this.im-z.im); }
  mul(z){ return new Complex(this.re*z.re - this.im*z.im, this.re*z.im + this.im*z.re); }
  div(z){
    const d = z.re*z.re + z.im*z.im;
    if(d===0) return new Complex(Infinity, Infinity);
    return new Complex((this.re*z.re + this.im*z.im)/d, (this.im*z.re - this.re*z.im)/d);
  }
  get magnitude(){ return Math.hypot(this.re, this.im); }
  get phase(){ return Math.atan2(this.im, this.re); }
  static i(){ return new Complex(0,1); }
}

/* ---------------- helpers ---------------- */
function createTable(el, headers, rows){
  if (!el) return;
  el.innerHTML = `<thead><tr>${headers.map(h=>`<th>${h}</th>`).join('')}</tr></thead>
  <tbody>${rows.map(r=>`<tr>${r.map(c=>`<td>${c}</td>`).join('')}</tr>`).join('')}</tbody>`;
}
function createPlot(chartId, datasets, title, yTitle, extra={}){
  const canvas = document.getElementById(chartId);
  if(!canvas) return;
  const ctx = canvas.getContext('2d');
  if(charts[chartId]) charts[chartId].destroy();

  // background bands plugin + labels + vertical modal dots
  const bandsPlugin = {
    id: 'bands_'+chartId,
    afterDraw(chart){
      const {ctx, chartArea, scales} = chart;
      const x = scales.x, y = scales.y;
      if (!x || !chartArea) return;

      // 1/3-octave shaded
      const pow6 = Math.pow(2,1/6);
      for (let i=0;i<FREQS.length;i++){
        const fc = Number(FREQS[i]);
        const x0 = x.getPixelForValue(fc/pow6);
        const x1 = x.getPixelForValue(fc*pow6);
        const left = Math.max(x0, chartArea.left), right = Math.min(x1, chartArea.right);
        if (right<=chartArea.left || left>=chartArea.right) continue;
        ctx.save();
        ctx.fillStyle = (i%2===0) ? 'rgba(200,200,200,0.05)' : 'rgba(180,180,180,0.03)';
        ctx.fillRect(left, chartArea.top, right-left, chartArea.bottom-chartArea.top);
        ctx.restore();
      }
      // x labels for each 1/3-oct center
      ctx.save();
      ctx.font = '11px Arial'; ctx.fillStyle='#222'; ctx.textAlign='center';
      const yText = chartArea.bottom + 14;
      FREQS.forEach(fc=>{
        const xp = x.getPixelForValue(fc);
        if (xp>=chartArea.left && xp<=chartArea.right) ctx.fillText(String(fc), xp, yText);
      });
      ctx.restore();

      // modal dots provided in extra.modalDots
      if (extra.modalDots && Array.isArray(extra.modalDots) && typeof extra.zInfDb==='number'){
        ctx.save();
        ctx.fillStyle = '#222'; ctx.textAlign='center';
        extra.modalDots.forEach(md=>{
          const xp = x.getPixelForValue(md.f);
          const yp = y.getPixelForValue(extra.zInfDb);
          if (xp>=chartArea.left && xp<=chartArea.right && yp>=chartArea.top && yp<=chartArea.bottom){
            ctx.beginPath(); ctx.arc(xp, yp, 3, 0, Math.PI*2); ctx.fill();
            ctx.font = '11px Arial';
            ctx.fillText(`f${md.n}${md.m}`, xp, yp-8);
          }
        });
        ctx.restore();
      }
    }
  };

  const chart = new Chart(ctx, {
    type:'line',
    data:{
      datasets: datasets.map(ds=>({
        ...ds,
        borderWidth: 2,
        pointRadius: 0,
        spanGaps: true,
        parsing: false,              // we provide {x,y}
        data: ds.data
      }))
    },
    options:{
      responsive:true,
      maintainAspectRatio:false,
      plugins:{
        title:{display:true, text:title, font:{size:16}},
        legend:{position:'top'}
      },
      scales:{
        x:{
          type:'logarithmic',
          min:16, max:1000,
          grid:{display:false},
          ticks:{display:false}
        },
        y:{
          type:'linear', position:'left',
          title:{display:true, text:yTitle}
        },
        yR:{
          type:'linear', position:'right',
          grid:{drawOnChartArea:false},
          title:{display:true, text:'aux'}
        }
      }
    },
    plugins:[bandsPlugin]
  });
  charts[chartId] = chart;
  return chart;
}

/* ---------------- physics core ---------------- */
function getInputs(){
  const h  = parseFloat(ui.thickInput.value)/1000;
  const E  = parseFloat(ui.eInput.value);
  const rho= parseFloat(ui.rhoInput.value);
  const nu = parseFloat(ui.nuInput.value);
  const Lx = parseFloat(ui.lxInput.value)/1000;
  const Ly = parseFloat(ui.lyInput.value)/1000;
  if(!(h>0&&E>0&&rho>0&&Lx>0&&Ly>0)) throw new Error('Invalid inputs');
  const c0 = 343;
  const D  = (E*Math.pow(h,3))/(12*(1-nu*nu));
  const C_Lp = Math.sqrt(E/(rho*(1-nu*nu)));
  const fc = (c0*c0/(2*Math.PI*h))*Math.sqrt(12*rho*(1-nu*nu)/E);
  const S = Lx*Ly;
  const n_const = (S*Math.sqrt(3))/(h*C_Lp); // approx modal density
  return {h,E,rho,nu,Lx,Ly,S,c0,D,fc,C_Lp,n_const};
}

// natural freq SS plate
function f_mn_SS(m,n, D, rho, h, Lx, Ly){
  return (Math.PI/2)*Math.sqrt(D/(rho*h))*((m/Lx)**2 + (n/Ly)**2);
}

function calcBasicInline(){
  const {h,E,rho,nu,D,fc,C_Lp,S,n_const} = getInputs();
  createTable(ui.basicInline,
    ['Property','Value','Unit'],
    [
      ['Bending Stiffness (D)', D.toExponential(3), 'N·m'],
      ['Critical Frequency (fc)', fc.toFixed(1), 'Hz'],
      ['Longitudinal Plate Wave (c_L,p)', C_Lp.toFixed(1), 'm/s'],
      ['Modal density n(f) ~ const', n_const.toFixed(3), 'modes/Hz']
    ]
  );
}

// mobility via modal sum at a point (SS), or grid-average
function mobilityY_at(f, inputs, x, y){
  const {D,rho,h,Lx,Ly} = inputs;
  const rho_s = rho*h;
  const omega = 2*Math.PI*f;
  const eta_f = 0.005 + 0.3/Math.sqrt(Math.max(1e-6,f));

  let sum = new Complex(0,0);
  const Pmax = 8, Qmax = 8; // decent default
  for(let p=1;p<=Pmax;p++){
    for(let q=1;q<=Qmax;q++){
      const fpq = f_mn_SS(p,q, D, rho, h, Lx, Ly);
      const w0  = 2*Math.PI*fpq;
      const psi = Math.sin(p*Math.PI*x/Lx)*Math.sin(q*Math.PI*y/Ly);
      const denom = (new Complex(w0*w0, w0*w0*eta_f)).sub(new Complex(omega*omega,0));
      sum = sum.add( new Complex(psi*psi,0).div(denom) );
    }
  }
  // Y = i*omega * 4/(rho_s*S) * sum
  return Complex.i().mul(new Complex(4*omega/(rho_s*inputs.S),0)).mul(sum);
}

function toXY(freqs, arr){ return freqs.map((f,i)=>({x:f, y: arr[i]})); }

function calcImpedance(){
  const ins = getInputs();

  // Z_inf (reference) & its dB (shown in table; horizontal line in figure)
  const C_L = Math.sqrt(ins.E/ins.rho);
  const Z_inf_ref = 2.3 * ins.rho * C_L * Math.pow(ins.h,2);
  const Z_inf_ref_dB = 20*Math.log10(Z_inf_ref);

  createTable(ui.infImpTable, ['Property','Value'],
    [['Infinite Plate Impedance (Z_inf, ref) [dB]', Z_inf_ref_dB.toFixed(2)]]);

  // datasets
  const mag_dB = [], re_dB = [], ph_deg = [];
  const tableRows = [];

  // grid-average?
  const useGrid = ui.useMesh.checked;
  const Nx = Math.max(1, parseInt(ui.nxInput.value||'1',10));
  const Ny = Math.max(1, parseInt(ui.nyInput.value||'1',10));

  for(const f of FREQS){
    let Z;
    if (useGrid){
      // arithmetic mean of complex Z over grid
      let acc = new Complex(0,0);
      for (let ix=0; ix<Nx; ix++){
        for (let iy=0; iy<Ny; iy++){
          const x = (ix+0.5)*ins.Lx/Nx;
          const y = (iy+0.5)*ins.Ly/Ny;
          const Y = mobilityY_at(f, ins, x, y);
          const Zi = (new Complex(1,0)).div(Y);
          acc = acc.add(Zi);
        }
      }
      Z = new Complex(acc.re/(Nx*Ny), acc.im/(Nx*Ny));
    }else{
      const x = parseFloat(ui.posX.value)/1000;
      const y = parseFloat(ui.posY.value)/1000;
      const Y = mobilityY_at(f, ins, x, y);
      Z = (new Complex(1,0)).div(Y);
    }

    const mag = Z.magnitude;
    const re  = Math.abs(Z.re); // dB化するので絶対値に
    const ph  = Z.phase*180/Math.PI;
    mag_dB.push( (mag>0)? 20*Math.log10(mag) : null );
    re_dB .push( (re>0)?  20*Math.log10(re)  : null );
    ph_deg.push( isFinite(ph)? ph : null );

    tableRows.push([
      f, 
      mag>0 ? (20*Math.log10(mag)).toFixed(2) : 'NaN',
      re>0  ? (20*Math.log10(re)).toFixed(2)  : 'NaN',
      isFinite(ph) ? ph.toFixed(1) : 'NaN'
    ]);
  }

  createTable(ui.impTable, ['Frequency (Hz)', '|Z| [dB]', 'Re(Z) [dB]', 'Phase [deg]'], tableRows);

  // modal dots (f11,f12,f21,f13,f31,f23,f32,f33)
  const dots = [];
  const pairs = [[1,1],[1,2],[2,1],[1,3],[3,1],[2,3],[3,2],[3,3]];
  pairs.forEach(([m,n])=>{
    const fmn = f_mn_SS(m,n, ins.D, ins.rho, ins.h, ins.Lx, ins.Ly);
    if (fmn>=FREQS[0] && fmn<=FREQS[FREQS.length-1]) dots.push({f:fmn, m, n});
  });

  // datasets build (Infinite only by default; checkboxes add more)
  const sets = [
    {label:'Z_inf_ref [dB]', data: toXY(FREQS, FREQS.map(()=>Z_inf_ref_dB)), borderColor:'rgba(0,160,0,0.9)', borderDash:[6,3]}
  ];
  if (ui.showMag.checked) sets.push({label:'|Z| [dB]', data: toXY(FREQS, mag_dB), borderColor:'rgba(220,40,80,1)'});
  if (ui.showRe .checked) sets.push({label:'Re(Z) [dB]', data: toXY(FREQS, re_dB),  borderColor:'rgba(40,120,220,1)', yAxisID:'y'});
  if (ui.showPh .checked) sets.push({label:'Phase [deg]', data: toXY(FREQS, ph_deg), borderColor:'rgba(240,180,40,1)', yAxisID:'yR'});

  createPlot('impedance-chart', sets, 'Driving-Point Impedance (SS, 16–1000 Hz)', 'Level / dB', {
    zInfDb: Z_inf_ref_dB,
    modalDots: dots
  });

  // formulas (brief)
  ui.impFormula.innerHTML = `<p><b>Y<sub>dp</sub></b> = iω·4/(ρh·S)·Σ<sub>p,q</sub> { ψ² / [ω<sub>pq</sub>²(1+iη)-ω²] }, &nbsp; Z = 1/Y</p>`;
}

function calcSTL(){
  const {rho,h,fc} = getInputs();
  const eta = 0.01;
  const mpp = rho*h;
  const R_mass = FREQS.map(f=>20*Math.log10(mpp*f)-42.5);
  const R_coin = R_mass.map((r, i)=>{
    const f = FREQS[i], fr=f/fc;
    if (Math.abs(fr-1)<1e-6) return r;
    const c = (2*eta)/(Math.PI*fr)*(1/Math.pow(1-fr*fr,2));
    return r - 10*Math.log10(Math.abs(c)+1);
  });

  createPlot('stl-chart', [
    {label:'Mass law', data: toXY(FREQS, R_mass), borderColor:'rgba(60,140,220,1)'},
    {label:'With coincidence', data: toXY(FREQS, R_coin), borderColor:'rgba(220,90,120,1)'}
  ], 'STL (16–1000 Hz)', 'STL [dB]');
  createTable(ui.stlTable, ['f [Hz]','Mass law','With coincidence'],
    FREQS.map((f,i)=>[f, R_mass[i].toFixed(1), R_coin[i].toFixed(1)]));
  ui.stlFormula.innerHTML = `<p>TL ≈ 20log₁₀(m''f)−42.5, coincidence correction near f<sub>c</sub></p>`;
}

function calcRadiation(){
  const {h,Lx,Ly,c0,fc} = getInputs();
  const l = 2*(Lx+Ly), S=Lx*Ly, lambda_c = c0/fc;
  const sigma_simple = FREQS.map(f=>{
    if (f>fc) return 1;
    if (Math.abs(f-fc)<1e-6) return Math.min(1, 0.45*Math.sqrt(l/lambda_c));
    return Math.min(1, (l*lambda_c/(Math.PI**2*S))*Math.sqrt(f/fc));
  });

  createPlot('rad-chart', [
    {label:'σ (simple)', data: toXY(FREQS, sigma_simple), borderColor:'rgba(20,160,120,1)'}
  ], 'Radiation Efficiency (16–1000 Hz)', 'σ [-]');
  createTable(ui.radTable, ['f [Hz]','σ (simple)'],
    FREQS.map((f,i)=>[f, sigma_simple[i].toFixed(3)]));
  ui.radFormula.innerHTML = `<p>Simple σ model with sub/coincidence behavior.</p>`;
}

/* ---------------- UI init & wiring ---------------- */
function onMaterialSelect(){
  const name = ui.materialSelect.value;
  const props = MATERIAL_PROPERTIES[name];
  if (!props) return;
  ui.eInput.value   = props.E.toExponential(3);
  ui.rhoInput.value = props.rho;
  ui.nuInput.value  = props.nu;
}

function setAveragingUIState(){
  const use = ui.useMesh.checked;
  ui.posX.disabled = use;
  ui.posY.disabled = use;
}

function calculateAll(){
  calcBasicInline();
  calcImpedance();
  calcSTL();
  calcRadiation();
}

function toggleImpedanceDataset(kind){
  calcImpedance();
}

/* tabs */
function wireTabs(){
  const btns = document.querySelectorAll('.tab-btn');
  btns.forEach(b=>{
    b.addEventListener('click', ()=>{
      btns.forEach(x=>x.classList.remove('active'));
      b.classList.add('active');
      const key = b.getAttribute('data-tab');
      document.querySelectorAll('.tab').forEach(t=>t.classList.remove('active'));
      document.getElementById('tab-'+key).classList.add('active');
    });
  });
}

/* boot */
document.addEventListener('DOMContentLoaded', () => {
  // refs
  ui = {
    materialSelect: document.getElementById('material-select'),
    thickInput: document.getElementById('thick-input'),
    lxInput: document.getElementById('lx-input'),
    lyInput: document.getElementById('ly-input'),
    eInput: document.getElementById('e-input'),
    rhoInput: document.getElementById('rho-input'),
    nuInput: document.getElementById('nu-input'),
    posX: document.getElementById('pos-x'),
    posY: document.getElementById('pos-y'),
    useMesh: document.getElementById('use-mesh'),
    nxInput: document.getElementById('nx-input'),
    nyInput: document.getElementById('ny-input'),
    baffleCond: document.getElementById('baffle-cond'),
    basicInline: document.getElementById('basic-inline'),
    impChart: document.getElementById('impedance-chart'),
    impTable: document.getElementById('impedance-table'),
    infImpTable: document.getElementById('inf-impedance-table'),
    stlChart: document.getElementById('stl-chart'),
    stlTable: document.getElementById('stl-table'),
    radChart: document.getElementById('rad-chart'),
    radTable: document.getElementById('rad-table'),
    impFormula: document.getElementById('impedance-formula'),
    stlFormula: document.getElementById('stl-formula'),
    radFormula: document.getElementById('rad-formula'),
    showMag: document.getElementById('show-mag'),
    showRe: document.getElementById('show-re'),
    showPh: document.getElementById('show-ph')
  };

  function buildMaterialOptions() {
    if (!ui.materialSelect) return;
    const opts = Object.keys(MATERIAL_PROPERTIES)
      .map(name => `<option value="${name}">${name}</option>`).join('');
    ui.materialSelect.innerHTML = opts;
  }

  function setDefaults() {
    ui.materialSelect.value = 'Concrete';
    ui.thickInput.value = 180;
    ui.lxInput.value = 3000;
    ui.lyInput.value = 4000;
    ui.posX.value = 1500;
    ui.posY.value = 2000;
    if (ui.useMesh) ui.useMesh.checked = true;
    if (ui.nxInput) ui.nxInput.value = 5;
    if (ui.nyInput) ui.nyInput.value = 6;
  }

  function applyMaterialToFields() {
    const props = MATERIAL_PROPERTIES[ui.materialSelect.value];
    if (!props) return;
    ui.eInput.value   = Number(props.E).toExponential(3);
    ui.rhoInput.value = props.rho;
    ui.nuInput.value  = props.nu;
  }

  function setAveragingUIStateLocal() {
    const use = ui.useMesh?.checked;
    if (ui.posX) ui.posX.disabled = !!use;
    if (ui.posY) ui.posY.disabled = !!use;
  }

  buildMaterialOptions();
  setDefaults();
  applyMaterialToFields();
  setAveragingUIStateLocal();
  wireTabs();                 // tabs wired

  // 初回計算
  calculateAll();

  // ====== イベント束ね ======
  ui.materialSelect.addEventListener('change', () => {
    applyMaterialToFields();
    calculateAll();
  });

  ['thickInput','lxInput','lyInput','eInput','rhoInput','nuInput','baffleCond']
    .forEach(k => {
      ui[k]?.addEventListener('input',  calculateAll);
      ui[k]?.addEventListener('change', calculateAll);
    });

  ui.useMesh?.addEventListener('change', () => { setAveragingUIStateLocal(); calculateImpedance(); });
  ui.nxInput?.addEventListener('input',  calculateImpedance);
  ui.nxInput?.addEventListener('change', calculateImpedance);
  ui.nyInput?.addEventListener('input',  calculateImpedance);
  ui.nyInput?.addEventListener('change', calculateImpedance);

  // Grid OFF のときだけ位置で即時更新
  const posInstant = () => { if (!ui.useMesh?.checked) calculateImpedance(); };
  ui.posX?.addEventListener('input',  posInstant);
  ui.posX?.addEventListener('change', posInstant);
  ui.posY?.addEventListener('input',  posInstant);
  ui.posY?.addEventListener('change', posInstant);

  // 表示トグル
  ui.showMag?.addEventListener('change', () => calculateImpedance());
  ui.showRe ?.addEventListener('change', () => calculateImpedance());
  ui.showPh ?.addEventListener('change', () => calculateImpedance());
});
