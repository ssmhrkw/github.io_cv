/* ===== Constants ===== */
const MATERIAL_PROPERTIES = {
  "Custom":       { rho: 0,     E: 0,       nu: 0     },
  "Concrete":     { rho: 2.3e3, E: 2.1e10,  nu: 0.005 },
  "Gypsum Board": { rho: 0.8e3, E: 0.18e10, nu: 0.005 },
  "Plywood":      { rho: 0.6e3, E: 0.5e10,  nu: 0.30  },
  "Glass":        { rho: 2500,  E: 70e9,    nu: 0.23  }
};
const FREQS = [16,20,25,31.5,40,50,63,80,100,125,160,200,250,315,400,500,630,800,1000,1250,1600,2000,2500,3150,4000,5000,6300,8000];
const THIRD_OCT = [16,31.5,63,125,250,500,1000,2000,4000,8000];
const F_IMP = FREQS.filter(f => f <= 1000); // Impedance up to 1 kHz

/* ===== State/UI ===== */
let charts = {};
let ui = {};

/* ===== Complex ===== */
class Complex {
  constructor(re=0, im=0){ this.re=re; this.im=im; }
  add(z){ return new Complex(this.re+z.re, this.im+z.im); }
  sub(z){ return new Complex(this.re-z.re, this.im-z.im); }
  mul(z){ return new Complex(this.re*z.re - this.im*z.im, this.re*z.im + this.im*z.re); }
  div(z){
    const d = z.re*z.re + z.im*z.im;
    if (d===0) return new Complex(Infinity,Infinity);
    return new Complex((this.re*z.re + this.im*z.im)/d, (this.im*z.re - this.re*z.im)/d);
  }
  get magnitude(){ return Math.hypot(this.re, this.im); }
  get phase(){ return Math.atan2(this.im, this.re); }
  static i(){ return new Complex(0,1); }
}

/* ===== Helpers ===== */
function openTab(evt, id){
  document.querySelectorAll('.tab-content').forEach(t=>t.classList.remove('active'));
  document.querySelectorAll('.tab-button').forEach(b=>b.classList.remove('active'));
  document.getElementById(id).classList.add('active');
  if (evt && evt.currentTarget) evt.currentTarget.classList.add('active');
}
window.openTab = openTab;

function createTable(el, headers, rows){
  if (!el) return;
  el.innerHTML = `<thead><tr>${headers.map(h=>`<th>${h}</th>`).join('')}</tr></thead>` +
                 `<tbody>${rows.map(r=>`<tr>${r.map(c=>`<td>${c}</td>`).join('')}</tr>`).join('')}</tbody>`;
}

function ssModalFreq(p,q,{E,rho,nu,h,Lx,Ly}){
  const rho_s = rho*h;
  const D = (E*h**3)/(12*(1-nu**2));
  return (Math.PI/2)*Math.sqrt(D/rho_s) * ((p/Lx)**2 + (q/Ly)**2);
}

function ssModesOrder8(inputs){
  const pairs = [[1,1],[1,2],[2,1],[1,3],[3,1],[2,3],[3,2],[3,3]];
  return pairs.map(([p,q]) => ({p,q,f:ssModalFreq(p,q,inputs)}));
}

// Unicode subscript
const subDigit = d => ({'0':'₀','1':'₁','2':'₂','3':'₃','4':'₄','5':'₅','6':'₆','7':'₇','8':'₈','9':'₉'}[d]||d);
const sub = s => String(s).split('').map(subDigit).join('');

/* 1/3-oct band background plugin factory */
function bandPlugin(centers){
  const colors = [
    'rgba(230,159,0,0.06)','rgba(86,180,233,0.06)','rgba(0,158,115,0.06)',
    'rgba(240,228,66,0.06)','rgba(0,114,178,0.06)','rgba(213,94,0,0.06)',
    'rgba(204,121,167,0.06)','rgba(160,160,160,0.06)','rgba(120,190,190,0.06)',
    'rgba(190,120,190,0.06)','rgba(190,190,120,0.06)','rgba(120,120,190,0.06)'
  ];
  const pow6 = Math.pow(2,1/6);
  return {
    id: 'band_coloring_' + Math.random().toString(36).slice(2),
    afterDraw(chart){
      const {ctx, chartArea:area, scales:{x}} = chart;
      if (!x) return;
      ctx.save();
      centers.forEach((fc,i)=>{
        const fl = fc / pow6, fu = fc * pow6;
        const x0 = x.getPixelForValue(fl), x1 = x.getPixelForValue(fu);
        const left = Math.max(x0, area.left), right = Math.min(x1, area.right);
        if (right>left){
          ctx.fillStyle = colors[i % colors.length];
          ctx.fillRect(left, area.top, right-left, area.bottom-area.top);
        }
      });
      ctx.restore();
    }
  };
}

/* ===== DOM Ready ===== */
document.addEventListener('DOMContentLoaded', () => {
  // ---- UI refs ----
  ui = {
    materialSelect: document.getElementById('material-select'),
    thickInput:     document.getElementById('thick-input'),
    lxInput:        document.getElementById('lx-input'),
    lyInput:        document.getElementById('ly-input'),
    eInput:         document.getElementById('e-input'),
    rhoInput:       document.getElementById('rho-input'),
    nuInput:        document.getElementById('nu-input'),
    posX:           document.getElementById('pos-x'),
    posY:           document.getElementById('pos-y'),
    posBlock:       document.getElementById('pos-block'),
    useMesh:        document.getElementById('use-mesh'),
    nxInput:        document.getElementById('nx-input'),
    nyInput:        document.getElementById('ny-input'),
    nxBadge:        document.getElementById('nx-badge'),
    nyBadge:        document.getElementById('ny-badge'),
    baffleCond:     document.getElementById('baffle-cond'),
    basicInline:    document.getElementById('basic-inline'),
    impChart:       document.getElementById('impedance-chart'),
    impTable:       document.getElementById('impedance-table'),
    infImpTable:    document.getElementById('inf-impedance-table'),
    stlChart:       document.getElementById('stl-chart'),
    stlTable:       document.getElementById('stl-table'),
    radChart:       document.getElementById('rad-chart'),
    radTable:       document.getElementById('rad-table'),
    impFormula:     document.getElementById('impedance-formula'),
    stlFormula:     document.getElementById('stl-formula'),
    radFormula:     document.getElementById('rad-formula'),
    showMag:        document.getElementById('show-mag'),
    showRe:         document.getElementById('show-re'),
    showPh:         document.getElementById('show-ph')
  };

  // フォールバック付の「全再計算」呼び出し
  const runAll = () => {
    if (typeof calculateAll === 'function') { calculateAll(); }
    else if (typeof calculateAllTabs === 'function') { calculateAllTabs(); }
    else {
      // 最低限インピーダンスだけでも
      if (typeof calculateImpedance === 'function') calculateImpedance();
      if (typeof calculateSTL === 'function')        calculateSTL();
      if (typeof calculateRadiation === 'function')  calculateRadiation();
    }
  };

  // --- 材料選択肢を投入（未投入なら） ---
  if (ui.materialSelect && ui.materialSelect.options.length === 0 && typeof MATERIAL_PROPERTIES === 'object') {
    Object.keys(MATERIAL_PROPERTIES).forEach(name => {
      const opt = document.createElement('option');
      opt.value = name; opt.textContent = name;
      ui.materialSelect.appendChild(opt);
    });
  }

  // --- 既定値（ご要望どおり）を先にセット ---
  if (ui.materialSelect) ui.materialSelect.value = 'Concrete';
  if (ui.thickInput)     ui.thickInput.value     = 180;   // mm
  if (ui.lxInput)        ui.lxInput.value        = 3000;  // mm
  if (ui.lyInput)        ui.lyInput.value        = 4000;  // mm
  if (ui.posX)           ui.posX.value           = 1500;  // 任意
  if (ui.posY)           ui.posY.value           = 2000;
  if (ui.useMesh)        ui.useMesh.checked      = true;  // Grid平均を既定ONにするなら
  if (ui.nxInput)        ui.nxInput.value        = 5;
  if (ui.nyInput)        ui.nyInput.value        = 6;

  // --- 材料→E,ρ,ν を反映（これを先にやらないと計算が空振りする） ---
  if (typeof onMaterialSelect === 'function') onMaterialSelect();

  // --- Averaging UI（位置入力の有効/無効など） ---
  if (typeof setAveragingUIState === 'function') setAveragingUIState();

  // --- 初回描画（開いた瞬間に図・表を表示） ---
  runAll();

  // ---- イベント：共通ハンドラ ----
  const onAllChange = () => runAll();
  const onImpOnly   = () => { if (typeof calculateImpedance === 'function') calculateImpedance(); };

  // --- Material 変更：反映 → 全再計算 ---
  ui.materialSelect?.addEventListener('change', () => {
    onMaterialSelect?.();
    runAll();
  });

  // --- 主要入力の変更で全再計算 ---
  ['thickInput','lxInput','lyInput','eInput','rhoInput','nuInput','baffleCond'].forEach(k => {
    if (!ui[k]) return;
    ui[k].addEventListener('input',  onAllChange);
    ui[k].addEventListener('change', onAllChange);
  });

  // --- Grid average ON/OFF：UI切替 → インピーダンスのみ即再描画 ---
  ui.useMesh?.addEventListener('change', () => {
    setAveragingUIState?.();
    onImpOnly();
  });

  // --- Nx/Ny：上下△で変えたら即インピーダンス再描画（バッジ同期も） ---
  const syncBadges = () => {
    if (ui.nxBadge && ui.nxInput) ui.nxBadge.textContent = ui.nxInput.value;
    if (ui.nyBadge && ui.nyInput) ui.nyBadge.textContent = ui.nyInput.value;
  };
  ui.nxInput?.addEventListener('input',  () => { syncBadges(); onImpOnly(); });
  ui.nxInput?.addEventListener('change', () => { syncBadges(); onImpOnly(); });
  ui.nyInput?.addEventListener('input',  () => { syncBadges(); onImpOnly(); });
  ui.nyInput?.addEventListener('change', () => { syncBadges(); onImpOnly(); });

  // --- Position X/Y：Grid平均OFF のときだけ即インピーダンス再描画 ---
  const posInstantUpdate = () => { if (!ui.useMesh?.checked) onImpOnly(); };
  ui.posX?.addEventListener('input',  posInstantUpdate);
  ui.posX?.addEventListener('change', posInstantUpdate);
  ui.posY?.addEventListener('input',  posInstantUpdate);
  ui.posY?.addEventListener('change', posInstantUpdate);

  // --- 表示トグル（Mag/Re/Phase） ---
  ui.showMag?.addEventListener('change', () => toggleImpedanceDataset?.('mag'));
  ui.showRe ?.addEventListener('change', () => toggleImpedanceDataset?.('re'));
  ui.showPh ?.addEventListener('change', () => toggleImpedanceDataset?.('ph'));
});

function setAveragingUIState(){
  const use = ui.useMesh.checked;
  [...ui.posBlock.querySelectorAll('input')].forEach(inp=>{
    inp.disabled = use;
    inp.classList.toggle('disabled', use);
  });
  ui.nxInput.disabled = !use;
  ui.nyInput.disabled = !use;
  ui.nxBadge.classList.toggle('disabled', !use);
  ui.nyBadge.classList.toggle('disabled', !use);
}

function onMaterialSelect(){
  const props = MATERIAL_PROPERTIES[ui.materialSelect.value] || {};
  ui.eInput.value   = props.E?.toExponential ? props.E.toExponential(2) : (props.E ?? '');
  ui.rhoInput.value = props.rho ?? '';
  ui.nuInput.value  = props.nu ?? '';
  [ui.eInput,ui.rhoInput,ui.nuInput].forEach(el=>el.readOnly=false);
}

/* ===== Core calcs ===== */
function getCommonInputs(){
  const h   = parseFloat(ui.thickInput.value)/1000;
  const E   = parseFloat(ui.eInput.value);
  const rho = parseFloat(ui.rhoInput.value);
  const nu  = parseFloat(ui.nuInput.value);
  if (!(h>0 && E>0 && rho>0)) throw new Error('Thickness, E, ρ must be positive.');
  const Lx = parseFloat(ui.lxInput.value)/1000;
  const Ly = parseFloat(ui.lyInput.value)/1000;
  if (!(Lx>0 && Ly>0)) throw new Error('Lx and Ly must be positive.');
  const c0 = 343;
  const D  = (E*h**3)/(12*(1-nu**2));
  const C_Lp = Math.sqrt(E/(rho*(1-nu**2)));
  const fc = (c0*c0/(2*Math.PI*h)) * Math.sqrt(12*rho*(1-nu**2)/E);
  return {h,E,rho,nu,Lx,Ly,c0,D,fc,CLp:C_Lp};
}

function calculateAll(){
  try{
    renderBasicInline();
    renderImpedance();
    renderSTL();
    renderRadiation();
  }catch(e){ alert(e.message); }
}

/* ===== Basic inline ===== */
function renderBasicInline(){
  const {h,E,rho,nu,Lx,Ly,D,fc,CLp} = getCommonInputs();
  const S = Lx*Ly;
  const n_asym = (S*Math.sqrt(3))/(h*CLp);
  const tbl = `<table>
    <thead><tr><th>Property</th><th>Value</th><th>Unit</th></tr></thead>
    <tbody>
      <tr><td>Bending Stiffness (D)</td><td>${D.toExponential(2)}</td><td>N·m</td></tr>
      <tr><td>Critical Frequency (f<sub>c</sub>)</td><td>${fc.toFixed(2)}</td><td>Hz</td></tr>
      <tr><td>Longitudinal Plate Wave Speed (c<sub>L,p</sub>)</td><td>${CLp.toFixed(2)}</td><td>m/s</td></tr>
      <tr><td>Modal Density n(f) (asym.)</td><td>${n_asym.toFixed(3)}</td><td>modes/Hz</td></tr>
    </tbody></table>`;
  ui.basicInline.innerHTML = tbl;
  if (window.MathJax) MathJax.typeset();
}

/* ===== Impedance ===== */
let impDatasets = null;

function renderImpedance(){
  const {E,rho,nu,h,Lx,Ly,fc} = getCommonInputs();
  const S = Lx*Ly, rho_s = rho*h;

  const C_L = Math.sqrt(E/rho);
  const Z_inf_ref = 2.3 * rho * C_L * h*h;
  const Z_INF_DB = (Z_inf_ref>0 && isFinite(Z_inf_ref)) ? 20*Math.log10(Z_inf_ref) : null;

  createTable(ui.infImpTable, ['Property','Value'], [
    ['Infinite Plate Impedance (Z_inf, ref) [dB]', (Z_INF_DB==null?'NaN':Z_INF_DB.toFixed(2)+' dB re 1 Ns/m³')]
  ]);

  let points = [];
  if (ui.useMesh.checked){
    const Nx = Math.max(1, parseInt(ui.nxInput.value||'1',10));
    const Ny = Math.max(1, parseInt(ui.nyInput.value||'1',10));
    for (let ix=1; ix<=Nx; ix++){
      const x = (ix-0.5)/Nx * Lx;
      for (let iy=1; iy<=Ny; iy++){
        const y = (iy-0.5)/Ny * Ly;
        points.push({x,y});
      }
    }
  }else{
    const x = Math.max(0, Math.min(Lx, parseFloat(ui.posX.value)/1000));
    const y = Math.max(0, Math.min(Ly, parseFloat(ui.posY.value)/1000));
    points.push({x,y});
  }

  const results = { mag_dB: [], re_dB: [], ph: [], table: [] };
  const Pmax=8, Qmax=8;

  for (const f of F_IMP){
    const omega = 2*Math.PI*f;
    const eta_f = 0.005 + 0.3/Math.sqrt(Math.max(1e-6,f));

    let Zsum_re = 0, Zsum_im = 0, Np = points.length;

    for (const pt of points){
      let Ysum = new Complex(0,0);
      for (let p=1;p<=Pmax;p++){
        for (let q=1;q<=Qmax;q++){
          const fpq = ssModalFreq(p,q,{E,rho,nu,h,Lx,Ly});
          const opq = 2*Math.PI*fpq;
          const psi = Math.sin(p*Math.PI*pt.x/Lx)*Math.sin(q*Math.PI*pt.y/Ly);
          const denom = (new Complex(opq*opq, opq*opq*eta_f)).sub(new Complex(omega*omega,0));
          Ysum = Ysum.add( (new Complex(psi*psi,0)).div(denom) );
        }
      }
      const Ydp = Complex.i().mul(new Complex(4*omega/(rho_s*S),0)).mul(Ysum);
      const d = Ydp.re*Ydp.re + Ydp.im*Ydp.im;
      let Zr=Infinity, Zi=Infinity;
      if (d>0 && isFinite(d)){
        Zr =  Ydp.re/d;
        Zi = -Ydp.im/d;
      }
      Zsum_re += Zr;
      Zsum_im += Zi;
    }

    let Zr_avg = Zsum_re / Np;
    let Zi_avg = Zsum_im / Np;

    let mag_dB=null, re_dB=null, ph=null;
    if (isFinite(Zr_avg) && isFinite(Zi_avg)){
      const mag = Math.hypot(Zr_avg, Zi_avg);
      mag_dB = 20*Math.log10(Math.max(mag,1e-12));
      re_dB  = 20*Math.log10(Math.max(Math.abs(Zr_avg),1e-12));
      ph = Math.atan2(Zi_avg, Zr_avg)*180/Math.PI;
    }
    results.mag_dB.push(mag_dB);
    results.re_dB.push(re_dB);
    results.ph.push(ph);
    results.table.push([f,
      mag_dB!=null?mag_dB.toFixed(2):'NaN',
      re_dB!=null?re_dB.toFixed(2):'NaN',
      ph!=null?ph.toFixed(1):'NaN'
    ]);
  }
  createTable(ui.impTable, ['Frequency (Hz)','|Z| dB','Re(Z) dB','Phase (deg)'], results.table);

  const modes8raw = ssModesOrder8({E,rho,nu,h,Lx,Ly});
  const modes8 = modes8raw
    .filter(m => m.f>=F_IMP[0] && m.f<=F_IMP[F_IMP.length-1])
    .map(m => ({x:m.f, y:Z_INF_DB, label:`f${sub(m.p)}${sub(m.q)}`}));

  impDatasets = [
    { key:'mag', label:'|Z|', data:results.mag_dB, yAxisID:'y',
      borderColor:'rgba(255,99,132,1)', borderWidth:2, pointRadius:2, hidden: true },
    { key:'re',  label:'Re(Z)', data:results.re_dB, yAxisID:'y1',
      borderColor:'rgba(54,162,235,1)', borderWidth:2, pointRadius:2, hidden: true },
    { key:'ph',  label:'Phase', data:results.ph, yAxisID:'y2',
      borderColor:'rgba(255,206,86,1)', borderWidth:2, pointRadius:2, hidden: true },
    { key:'zinf',label:'Z_inf_ref', data:F_IMP.map(()=>Z_INF_DB), yAxisID:'y',
      borderColor:'rgba(0,160,0,0.95)', borderDash:[6,3], borderWidth:2, pointRadius:0, hidden:false },
    { key:'dots', type:'scatter', label:'f_nm', data:modes8, yAxisID:'y',
      pointRadius:4, pointHoverRadius:6, showLine:false, hidden:false }
  ];

  const labelPlugin = {
    id:'imp_labels_' + Math.random().toString(36).slice(2),
    afterDraw(chart){
      const {ctx, chartArea:area} = chart;
      const ds = chart.data.datasets.find(d=>d.key==='dots');
      if (!ds || ds.hidden) return;
      const meta = chart.getDatasetMeta(chart.data.datasets.indexOf(ds));
      ctx.save();
      ctx.font = '11px Arial';
      ctx.fillStyle = '#333';
      ds.data.forEach((pt,i)=>{
        const el = meta.data[i];
        if (!el) return;
        const {x:xp,y:yp} = el.getProps(['x','y'], true);
        if (xp>area.left && xp<area.right && yp>area.top && yp<area.bottom){
          ctx.fillText(pt.label, xp+6, yp-6);
        }
      });
      ctx.restore();
    }
  };

  if (charts.imp) charts.imp.destroy();
  charts.imp = new Chart(ui.impChart.getContext('2d'), {
    type:'line',
    data:{ labels:F_IMP, datasets: impDatasets },
    options:{
      responsive:true, maintainAspectRatio:false,
      plugins:{
        title:{display:true,text:'Driving-Point Impedance (SS, up to 1 kHz)'},
        legend:{position:'top'}
      },
      interaction:{mode:'index',intersect:false},
      scales:{
        x:{type:'logarithmic', title:{display:true,text:'Frequency (Hz)'},
           ticks:{callback:v=>[16,20,25,31.5,40,50,63,80,100,125,160,200,250,315,400,500,630,800,1000].includes(v)?v:null},
           min:16, max:1000},
        y:{position:'left', title:{display:true,text:'Level [dB]'}},
        y1:{position:'right', title:{display:true,text:'Re(Z) [dB]'}, grid:{drawOnChartArea:false}},
        y2:{position:'right', title:{display:true,text:'Phase [deg]'}, grid:{drawOnChartArea:false}, offset:true}
      }
    },
    plugins:[ bandPlugin(F_IMP), labelPlugin ]
  });

  // checkbox handlers (ensure sync even after re-render)
  document.getElementById('show-mag').onchange = ()=>toggleImpedanceDataset('mag');
  document.getElementById('show-re').onchange  = ()=>toggleImpedanceDataset('re');
  document.getElementById('show-ph').onchange  = ()=>toggleImpedanceDataset('ph');

  ui.impFormula.innerHTML = `<span>SS modal sum。平均は「線形の複素インピーダンス」を算術平均してから dB/位相へ変換。初期表示は Z_inf_ref と f<sub>nm</sub>。</span>`;
}

function toggleImpedanceDataset(kind){
  if (!charts.imp || !impDatasets) return;
  const key = ({mag:'mag',re:'re',ph:'ph'})[kind];
  const idx = charts.imp.data.datasets.findIndex(d=>d.key===key);
  if (idx<0) return;
  const chk = kind==='mag'?document.getElementById('show-mag')
            : kind==='re' ?document.getElementById('show-re')
            : document.getElementById('show-ph');
  charts.imp.data.datasets[idx].hidden = !chk.checked;
  charts.imp.update();
}

/* ===== STL ===== */
function renderSTL(){
  const {rho,h,fc} = getCommonInputs();
  const m = rho*h, eta=0.01;
  const Rm = FREQS.map(f=>20*Math.log10(m*f)-42.5);
  const Rc = Rm.map((r,i)=>{
    const f=FREQS[i], fr=f/fc;
    if (Math.abs(fr-1)<1e-6) return r;
    const c=(2*eta)/(Math.PI*fr) * (1/Math.pow(Math.abs(1-fr*fr),2));
    return r - 10*Math.log10(Math.abs(c)+1);
  });

  if (charts.stl) charts.stl.destroy();
  charts.stl = new Chart(ui.stlChart.getContext('2d'),{
    type:'line',
    data:{labels:FREQS,datasets:[
      {label:'Mass Law (Normal Incidence)',data:Rm,borderWidth:2,pointRadius:2},
      {label:'With Coincidence Dip',data:Rc,borderWidth:2,pointRadius:2}
    ]},
    options:{
      responsive:true,maintainAspectRatio:false,
      plugins:{title:{display:true,text:'Sound Transmission Loss (STL)'}},
      scales:{
        x:{type:'logarithmic',title:{display:true,text:'Frequency (Hz)'},
           ticks:{callback:v=>THIRD_OCT.includes(v)?v:null},
           min:Math.min(...FREQS),max:Math.max(...FREQS)},
        y:{position:'left',title:{display:true,text:'STL (dB)'}}
      }
    },
    plugins:[ bandPlugin(FREQS) ]
  });

  createTable(ui.stlTable, ['Frequency (Hz)','Mass Law (dB)','With Coincidence (dB)'],
    FREQS.map((f,i)=>[f,Rm[i].toFixed(1),Rc[i].toFixed(1)]));

  ui.stlFormula.textContent = 'm″=ρh, TL₀=20log₁₀(m″f)−42.5。Coincidence 簡易補正。';
}

/* ===== Radiation ===== */
function renderRadiation(){
  const {h,Lx,Ly,c0,fc} = getCommonInputs();
  const C_BC=1, C_OB=parseFloat(document.getElementById('baffle-cond').value);
  const l=2*(Lx+Ly), S=Lx*Ly, lambda_c=c0/fc, L1=Math.min(Lx,Ly), L2=Math.max(Lx,Ly), kfc=2*Math.PI*fc/c0;

  const sigma_simple = FREQS.map(f=>{
    let s; if (f>fc) s=1; else if (Math.abs(f-fc)<1e-6) s=0.45*Math.sqrt(l/lambda_c);
    else s=(l*lambda_c/(Math.PI**2*S))*Math.sqrt(f/fc);
    return Math.max(1e-4, Math.min(s,1));
  });

  const sigma_lepp = FREQS.map(f=>{
    let s; const mu=Math.sqrt(fc/f);
    if (f>fc) s=1/Math.sqrt(1-fc/f);
    else if (Math.abs(f-fc)<1e-6) s=(0.5-0.15*L1/L2)*Math.sqrt(kfc*Math.sqrt(L1));
    else {
      const k=2*Math.PI*f/c0; const mu2_1=mu*mu-1; if (mu2_1<=0||k<=0) return 1e-4;
      const term1 = l/(2*Math.PI*k*S*Math.sqrt(mu2_1));
      const term2 = 2*Math.atanh(1/mu);
      const term3 = 2*mu/mu2_1;
      const term4 = C_BC*C_OB - Math.pow(mu,-8)*(C_BC*C_OB-1);
      s = term1*(term2+term3)*term4;
    }
    return Math.max(1e-4, Math.min(s,1));
  });

  const sigma_wallace = FREQS.map(f=>{
    const k=2*Math.PI*f/c0; let sum=0;
    for (let p=1;p<=5;p++){
      for (let q=1;q<=5;q++){
        let spq; const p_odd=p%2!==0, q_odd=q%2!==0;
        if (p_odd && q_odd){
          spq=(32*k**2*Lx*Ly/(Math.PI**5*p**2*q**2))*(1-(k**2*Lx*Ly/12)*
             ((1-8/(p*Math.PI)**2)*Lx/Ly+(1-8/(q*Math.PI)**2)*Ly/Lx));
        }else if (p_odd!==q_odd){
          if (p_odd){
            spq=(8*k**4*Lx**3*Ly/(3*Math.PI**5*p**2*q**2))*(1-(k**2*Lx*Ly/20)*
               ((1-8/(p*Math.PI)**2)*Lx/Ly+(1-24/(q*Math.PI)**2)*Ly/Lx));
          }else{
            spq=(8*k**4*Ly**3*Lx/(3*Math.PI**5*q**2*p**2))*(1-(k**2*Ly*Lx/20)*
               ((1-8/(q*Math.PI)**2)*Ly/Lx+(1-24/(p*Math.PI)**2)*Lx/Ly));
          }
        }else{
          spq=(2*k**6*Lx**3*Ly**3/(15*Math.PI**5*p**2*q**2))*(1-(5*k**2*Lx**2/64)*
             ((1-24/(p*Math.PI)**2)*Lx/Ly+(1-24/(q*Math.PI)**2)*Ly/Lx));
        }
        sum += Math.max(0, spq);
      }
    }
    return Math.max(1e-4, Math.min(sum,1));
  });

  if (charts.rad) charts.rad.destroy();
  charts.rad = new Chart(ui.radChart.getContext('2d'),{
    type:'line',
    data:{labels:FREQS,datasets:[
      {label:'Simple Model',data:sigma_simple,borderWidth:2,pointRadius:2},
      {label:'Leppington',data:sigma_lepp,borderWidth:2,pointRadius:2},
      {label:'Wallace (Modes Sum)',data:sigma_wallace,borderWidth:2,pointRadius:2}
    ]},
    options:{
      responsive:true,maintainAspectRatio:false,
      plugins:{title:{display:true,text:'Radiation Efficiency (σ)'}},
      scales:{
        x:{type:'logarithmic',title:{display:true,text:'Frequency (Hz)'},
           ticks:{callback:v=>THIRD_OCT.includes(v)?v:null},
           min:Math.min(...FREQS),max:Math.max(...FREQS)},
        y:{position:'left',title:{display:true,text:'Coefficient (σ)'},min:0,max:1.1}
      }
    },
    plugins:[ bandPlugin(FREQS) ]
  });

  createTable(ui.radTable, ['Frequency (Hz)','σ (Simple)','σ (Leppington)','σ (Wallace)'],
    FREQS.map((f,i)=>[f,sigma_simple[i].toFixed(3),sigma_lepp[i].toFixed(3),sigma_wallace[i].toExponential(2)]));

  ui.radFormula.textContent = 'σの近似式（Simple, Leppington, Wallace）を比較表示。';
}
