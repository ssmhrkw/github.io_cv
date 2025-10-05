# Write a standalone JS file (plates.js) extracted from the delivered HTML's <script> with slight modularization.
js_code = r"""/* plates.js — Rectangular plate driving-point impedance (live) */
/* Requires: Plotly 2.x loaded before this script. */
/* Expects the following element IDs to exist in HTML:
   material,E,rho,nu,thick,Lx,Ly,x0,y0,eta,fmin,fmax,mx,ny,
   showZ,showY,logx,logy,reim,avgOn,avgN,seed,
   file,exportPng,reset,plot,kpiD,kpimm,kpif11
*/
(function(){
  "use strict";

  // --- Small utilities ------------------------------------------------------
  const $ = (id) => document.getElementById(id);
  const clamp = (v, lo, hi) => Math.min(hi, Math.max(lo, v));
  const EPS = 1e-30;

  function seededRandom(seed){
    let s = (seed >>> 0) || 1;
    return function(){
      // xorshift32
      s ^= s << 13; s ^= s >>> 17; s ^= s << 5;
      return (s >>> 0) / 4294967296;
    }
  }

  // debounce to avoid excessive Plotly.react calls
  function debounce(fn, wait=50){
    let t;
    return (...args)=>{
      clearTimeout(t);
      t = setTimeout(()=>fn(...args), wait);
    };
  }

  // Minimal complex arithmetic
  function c(re, im){ return {re, im} }
  function cadd(a,b){ return c(a.re+b.re, a.im+b.im) }
  function cdiv(a,b){
    const d = b.re*b.re + b.im*b.im || EPS;
    return c( (a.re*b.re + a.im*b.im)/d, (a.im*b.re - a.re*b.im)/d );
  }
  function cinv(a){
    const d = a.re*a.re + a.im*a.im || EPS;
    return c(a.re/d, -a.im/d);
  }
  function cabs(a){ return Math.hypot(a.re, a.im) }
  const creal = (a)=>a.re, cimag = (a)=>a.im;

  // --- State ----------------------------------------------------------------
  const state = {
    material: "concrete",
    E: 30e9, rho: 2400, nu: 0.20, h: 0.150,
    Lx: 3.6, Ly: 2.7,
    x0: 1.8, y0: 1.35,
    eta: 0.01,
    fmin: 20, fmax: 5000,
    mx: 30, ny: 30,
    showZ: true, showY: false, logx: true, logy: true, reim: false,
    avgOn: false, avgN: 9, seed: 42,
    overlay: null
  };

  const presets = {
    concrete:{E:30e9, rho:2400, nu:0.20, h_mm:150},
    gypsum:{E:3.0e9, rho:800, nu:0.30, h_mm:12.5},
    wood:{E:10e9, rho:600, nu:0.30, h_mm:18}
  };

  function applyPreset(name){
    if(!presets[name]) return;
    const p = presets[name];
    $("E").value = (p.E/1e9).toFixed(1);
    $("rho").value = p.rho;
    $("nu").value = p.nu.toFixed(2);
    $("thick").value = p.h_mm;
    pullFromUI(); // update state & replot
  }

  // --- Physics core ---------------------------------------------------------
  function bendingStiffness(E, nu, h){ return E*h*h*h/(12*(1-nu*nu)); } // D [N*m]
  function massPerArea(rho, h){ return rho*h } // m' [kg/m2]

  function fmn(m,n,Lx,Ly,E,nu,rho,h){
    const D = bendingStiffness(E,nu,h);
    const mp = massPerArea(rho,h);
    const lam2 = Math.PI*Math.PI*((m*m)/(Lx*Lx) + (n*n)/(Ly*Ly));
    const omega = Math.sqrt(D/mp) * (lam2*lam2); // ω_mn = sqrt(D/mp) * λ^2^2
    return omega/(2*Math.PI);
  }

  function mobilityAtPoint(freqs, params){
    const {E,nu,rho,h,Lx,Ly,x0,y0,eta,mx,ny} = params;
    const C = 4/(rho*h*Lx*Ly);
    const Y = new Array(freqs.length).fill(0).map(()=>c(0,0));
    for(let m=1; m<=mx; m++){
      const sx2 = Math.sin(m*Math.PI*x0/Lx); // squared later
      const mm = m*m/(Lx*Lx);
      for(let n=1; n<=ny; n++){
        const sy2 = Math.sin(n*Math.PI*y0/Ly);
        const nn = n*n/(Ly*Ly);
        const lam2 = Math.PI*Math.PI*(mm + nn);
        const omega_mn = Math.sqrt(bendingStiffness(E,nu,h)/massPerArea(rho,h)) * (lam2*lam2);
        const phi2 = (sx2*sx2) * (sy2*sy2); // (sin)^2 * (sin)^2
        for(let i=0; i<freqs.length; i++){
          const w = 2*Math.PI*freqs[i];
          // jω / (ω_mn^2 - ω^2 + j η ω_mn ω)
          const den = c( (omega_mn*omega_mn - w*w), (eta*omega_mn*w) );
          const num = c(0, w);
          const frac = cdiv(num, den);
          Y[i] = c( Y[i].re + C*phi2*frac.re, Y[i].im + C*phi2*frac.im );
        }
      }
    }
    return Y;
  }

  function spatialAverageMobility(freqs, params){
    const {avgN, seed, Lx, Ly} = params;
    const rng = seededRandom(seed);
    let acc = new Array(freqs.length).fill(0).map(()=>c(0,0));
    for(let k=0;k<avgN;k++){
      const p = {...params};
      p.x0 = rng()*Lx;
      p.y0 = rng()*Ly;
      const Yk = mobilityAtPoint(freqs,p);
      for(let i=0;i<freqs.length;i++){
        acc[i] = c( acc[i].re + Yk[i].re, acc[i].im + Yk[i].im );
      }
    }
    return acc.map(z=>c(z.re/avgN, z.im/avgN));
  }

  // --- UI sync --------------------------------------------------------------
  const plotDebounced = debounce(plot, 10);

  function pullFromUI(){
    state.material = $("material").value;
    state.E = parseFloat($("E").value)*1e9;
    state.rho = parseFloat($("rho").value);
    state.nu = clamp(parseFloat($("nu").value), 0.0, 0.49);
    state.h = parseFloat($("thick").value)/1000;
    state.Lx = parseFloat($("Lx").value);
    state.Ly = parseFloat($("Ly").value);
    state.x0 = clamp(parseFloat($("x0").value), 0, state.Lx);
    state.y0 = clamp(parseFloat($("y0").value), 0, state.Ly);
    state.eta = Math.max(0.0001, parseFloat($("eta").value)/100);
    state.fmin = Math.max(1, parseFloat($("fmin").value));
    state.fmax = Math.max(state.fmin+10, parseFloat($("fmax").value));
    state.mx = Math.max(1, parseInt($("mx").value));
    state.ny = Math.max(1, parseInt($("ny").value));
    state.showZ = $("showZ").checked;
    state.showY = $("showY").checked;
    state.logx = $("logx").checked;
    state.logy = $("logy").checked;
    state.reim = $("reim").checked;
    state.avgOn = $("avgOn").checked;
    state.avgN = Math.max(2, parseInt($("avgN").value));
    state.seed = Math.max(0, parseInt($("seed").value));

    // KPIs
    const D = bendingStiffness(state.E, state.nu, state.h);
    const mp = massPerArea(state.rho, state.h);
    $("kpiD").textContent = D.toExponential(3);
    $("kpimm").textContent = mp.toFixed(1);
    const f11 = fmn(1,1,state.Lx,state.Ly,state.E,state.nu,state.rho,state.h);
    $("kpif11").textContent = f11.toFixed(1);

    plotDebounced();
  }

  function attachEvents(){
    $("material").addEventListener("change",(e)=>{
      const v = e.target.value;
      if(v!=="custom"){ applyPreset(v); }
    });
    // All inputs trigger pullFromUI WITHOUT resetting material selection
    ["E","rho","nu","thick","Lx","Ly","x0","y0","eta","fmin","fmax","mx","ny",
     "showZ","showY","logx","logy","reim","avgOn","avgN","seed"].forEach(id=>{
      $(id).addEventListener("input", pullFromUI);
      $(id).addEventListener("change", pullFromUI);
    });

    $("reset").addEventListener("click", ()=>{
      $("material").value = "concrete";
      applyPreset("concrete");
      $("Lx").value = 3.6; $("Ly").value = 2.7;
      $("x0").value = 1.8; $("y0").value = 1.35;
      $("eta").value = 1.0;
      $("fmin").value = 20; $("fmax").value = 5000;
      $("mx").value = 30; $("ny").value = 30;
      $("showZ").checked = true; $("showY").checked = false;
      $("logx").checked = true; $("logy").checked = true; $("reim").checked=false;
      $("avgOn").checked = false; $("avgN").value = 9; $("seed").value = 42;
      state.overlay = null;
      pullFromUI();
    });

    $("exportPng").addEventListener("click", ()=>{
      Plotly.downloadImage("plot", {format:"png", filename:"plate_impedance"});
    });

    $("file").addEventListener("change", handleFileOverlay, false);
  }

  async function handleFileOverlay(evt){
    const file = evt.target.files?.[0];
    if(!file) return;
    try{
      const text = await file.text();
      const lines = text.split(/\r?\n/).filter(x=>x.trim().length>0);
      let f=[], zr=[], zi=[];
      for(let i=0;i<lines.length;i++){
        const row = lines[i].trim();
        if(i===0 && /f.*Zreal.*Zimag/i.test(row)) continue; // header
        const cols = row.split(/[,;\t]/).map(s=>s.trim());
        if(cols.length < 3) continue;
        const ff = parseFloat(cols[0]), r = parseFloat(cols[1]), im = parseFloat(cols[2]);
        if(Number.isFinite(ff) && Number.isFinite(r) && Number.isFinite(im)){
          f.push(ff); zr.push(r); zi.push(im);
        }
      }
      if(f.length>0){
        state.overlay = {f, zr, zi};
        plotDebounced();
      }else{
        alert("CSVを解釈できませんでした。ヘッダ: f,Zreal,Zimag で数値を含めてください。");
        state.overlay = null;
      }
    }catch(err){
      console.error(err);
      alert("CSV読み込みでエラーが発生しました。");
    }
  }

  // --- Plot -----------------------------------------------------------------
  function plot(){
    const N = 800;
    const fmin = state.fmin, fmax = state.fmax;
    const freqs = Array.from({length:N}, (_,i)=> fmin*Math.pow(fmax/fmin, i/(N-1)) );

    const params = {...state};
    const Y = (state.avgOn ? spatialAverageMobility(freqs, params) : mobilityAtPoint(freqs, params));
    const Z = Y.map(y => cinv(y));

    const traces = [];

    if(state.showZ){
      traces.push({
        x: freqs, y: Z.map(cabs),
        type: "scatter", mode:"lines",
        name: "|Z| [N·s/m]"
      });
    }
    if(state.showY){
      traces.push({
        x: freqs, y: Y.map(cabs),
        type: "scatter", mode:"lines",
        name: "|Y| [m/N·s]"
      });
    }
    if(state.reim){
      traces.push({ x: freqs, y: Z.map(creal), type:"scatter", mode:"lines", name:"Re{Z}" });
      traces.push({ x: freqs, y: Z.map(cimag), type:"scatter", mode:"lines", name:"Im{Z}" });
    }

    if(state.overlay){
      const {f, zr, zi} = state.overlay;
      const Zm = zr.map((r, i)=> Math.hypot(r, zi[i] ?? 0));
      traces.push({
        x: f, y: Zm, type:"scatter", mode:"lines",
        name: "Overlay |Z| (CSV)", line:{dash:"dot"}
      });
    }

    const layout = {
      title: {text: "駆動点インピーダンス / モビリティ（SSSS板, モード和）", font:{size:16}},
      xaxis: {
        title: "Frequency [Hz]",
        type: state.logx ? "log" : "linear",
        gridcolor: "#e5e7eb",
        zerolinecolor: "#e5e7eb"
      },
      yaxis: {
        title: state.showY && !state.showZ ? "Mobility |Y| [m/N·s]" : "Impedance / Mobility (magnitude)",
        type: state.logy ? "log" : "linear",
        gridcolor: "#e5e7eb",
        zerolinecolor: "#e5e7eb",
        rangemode: "tozero"
      },
      legend:{orientation:"h", x:0, y:1.06},
      plot_bgcolor:"#ffffff",
      paper_bgcolor:"#ffffff",
      margin:{l:70,r:20,t:60,b:60}
    };

    Plotly.react("plot", traces, layout, {responsive:true, displaylogo:false});
  }

  // --- Init -----------------------------------------------------------------
  function init(){
    // Guard for required nodes
    const ids = ["material","E","rho","nu","thick","Lx","Ly","x0","y0","eta","fmin","fmax","mx","ny",
      "showZ","showY","logx","logy","reim","avgOn","avgN","seed",
      "file","exportPng","reset","plot","kpiD","kpimm","kpif11"];
    for(const id of ids){
      if(!$(id)){
        console.warn("[plates.js] Missing element:", id);
      }
    }
    attachEvents();
    $("material").value = "concrete";
    applyPreset("concrete"); // calls pullFromUI() -> plotDebounced()
  }

  if(document.readyState === "loading"){
    document.addEventListener("DOMContentLoaded", init);
  }else{
    init();
  }

})();"""
path = "/mnt/data/plates.js"
with open(path, "w", encoding="utf-8") as f:
    f.write(js_code)

path
