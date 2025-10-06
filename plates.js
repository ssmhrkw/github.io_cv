/* plates.js — spatial averaging + structural damping + 1/1-octave shading */
(function(){
  "use strict";
  const $ = (id) => document.getElementById(id);
  const clamp = (v, lo, hi) => Math.min(hi, Math.max(lo, v));
  const EPS = 1e-30;

  function c(re, im){ return {re, im}; }
  function cdiv(a,b){ const d=b.re*b.re+b.im*b.im||EPS; return {re:(a.re*b.re+a.im*b.im)/d, im:(a.im*b.re-a.re*b.im)/d}; }
  function cinv(a){ const d=a.re*a.re+a.im*a.im||EPS; return {re:a.re/d, im:-a.im/d}; }
  const cabs = (z)=> Math.hypot(z.re, z.im);

  const state = {
    material:"concrete",
    E:30e9, rho:2400, nu:0.20, h:0.200,
    Lx:3.6, Ly:2.7,
    x0mm:1800, y0mm:1350,
    ptmode:"single", Nx:3, Ny:3, avgMode:"linear",
    eta:0.01, fmin:20, fmax:5000,
    mx:30, ny:30,
    showZ:true, showY:false, logx:true, db:true,
    showOct:true, showInf:true, showModes:true,
    bandBg:true
  };

  const presets = {
    concrete:{E:30e9,rho:2400,nu:0.20,h_mm:200},
    gypsum:{E:3.0e9,rho:800,nu:0.30,h_mm:12.5},
    wood:{E:10e9,rho:600,nu:0.30,h_mm:18}
  };

  const OCT_CENTERS = [16, 31.5, 63, 125, 250, 500, 1000, 2000];
  const BAND_FACTOR = Math.SQRT2;

  function makeBandShapes(fmin, fmax){
    const shapes = [];
    let alt = true;
    for(const fc of OCT_CENTERS){
      const lo = fc/BAND_FACTOR, hi = fc*BAND_FACTOR;
      if(hi < fmin || lo > fmax) continue;
      shapes.push({
        type:"rect", xref:"x", yref:"paper",
        x0: Math.max(lo, fmin), x1: Math.min(hi, fmax),
        y0: 0, y1: 1,
        layer: "below",
        fillcolor: alt ? "rgba(99,102,241,0.10)" : "rgba(16,185,129,0.10)",
        line:{width:0}
      });
      alt = !alt;
    }
    return shapes;
  }

  const D  = (E,nu,h)=> E*h*h*h/(12*(1-nu*nu));
  const mp = (rho,h)=> rho*h;

  // --- 修正：ωmn = √(D/ρh) * k^2  （k^2 を“さらに二乗”しない）
  function fmn(m,n,Lx,Ly,E,nu,rho,h){
    const lam2 = Math.PI*Math.PI*((m*m)/(Lx*Lx) + (n*n)/(Ly*Ly)); // k^2
    const omg = Math.sqrt(D(E,nu,h)/mp(rho,h)) * lam2;            // ωmn
    return omg/(2*Math.PI);
  }

  // --- 修正：構造減衰 (1+iη) → 分母 = (ωmn^2 - ω^2) + i(η ωmn^2)
  function mobilityAtPoint(freqs, params){
    const {E,nu,rho,h,Lx,Ly,x0,y0,eta,mx,ny} = params;
    const C = 4/(rho*h*Lx*Ly); // = 4/(ρ_s S)
    const Y = new Array(freqs.length).fill(0).map(()=>({re:0,im:0}));
    for(let m=1; m<=mx; m++){
      const sx  = Math.sin(m*Math.PI*x0/Lx);
      const sx2 = sx*sx;
      for(let n=1; n<=ny; n++){
        const sy  = Math.sin(n*Math.PI*y0/Ly);
        const sy2 = sy*sy;
        const lam2   = Math.PI*Math.PI*( (m*m)/(Lx*Lx) + (n*n)/(Ly*Ly) ); // k^2
        const omg_mn = Math.sqrt(D(E,nu,h)/mp(rho,h)) * lam2;            // ωmn
        const omg2   = omg_mn*omg_mn;
        const phi2   = sx2 * sy2;
        for(let i=0;i<freqs.length;i++){
          const w = 2*Math.PI*freqs[i];
          const den = {re:(omg2 - w*w), im:(eta*omg2)};   // (1+iη)ωmn^2 - ω^2
          const num = {re:0, im:w};                       // iω
          const frac = cdiv(num, den);
          Y[i].re += C*phi2*frac.re;
          Y[i].im += C*phi2*frac.im;
        }
      }
    }
    return Y;
  }

  const Zinf_const = (E,nu,rho,h)=> 8*Math.sqrt( D(E,nu,h) * (rho*h) );

  function toDB20(arr, ref){ return arr.map(v => 20*Math.log10( Math.max(v, EPS) / ref )); }

  function interpY(xarr, yarr, x){
    if(x<=xarr[0]) return yarr[0];
    if(x>=xarr[xarr.length-1]) return yarr[yarr.length-1];
    let lo=0, hi=xarr.length-1;
    while(hi-lo>1){ const mid=(lo+hi)>>1; if(xarr[mid] > x) hi=mid; else lo=mid; }
    const x0=xarr[lo], x1=xarr[hi], y0=yarr[lo], y1=yarr[hi];
    const t = (Math.log(x) - Math.log(x0))/(Math.log(x1)-Math.log(x0));
    return y0*(1-t)+y1*t;
  }

  function enableGridInputs(){
    const grid = state.ptmode==="grid";
    $("Nx").disabled = !grid; $("Ny").disabled = !grid;
    $("x0mm").disabled = grid; $("y0mm").disabled = grid;
  }

  function getPositions(){
    if(state.ptmode==="single"){
      return [[state.x0mm/1000, state.y0mm/1000]];
    }else{
      const xs = Array.from({length:state.Nx}, (_,i)=> (i+0.5)/state.Nx * state.Lx);
      const ys = Array.from({length:state.Ny}, (_,j)=> (j+0.5)/state.Ny * state.Ly);
      const pts = [];
      for(const x of xs) for(const y of ys) pts.push([x,y]);
      return pts;
    }
  }

  function computeAveraged(freqs, params){
    const pts = getPositions();
    const n = pts.length;
    const refZ = 1.0, refY = 1.0;
    let accZ, accY;
    if(state.avgMode==="db"){
      accZ = new Array(freqs.length).fill(0);
      accY = new Array(freqs.length).fill(0);
      for(const [x,y] of pts){
        const Y = mobilityAtPoint(freqs, {...params, x0:x, y0:y});
        const Z = Y.map(cinv);
        const mz = Z.map(z=>20*Math.log10(Math.max(cabs(z),EPS)/refZ));
        const my = Y.map(yc=>20*Math.log10(Math.max(cabs(yc),EPS)/refY));
        for(let i=0;i<freqs.length;i++){ accZ[i]+=mz[i]; accY[i]+=my[i]; }
      }
      for(let i=0;i<freqs.length;i++){ accZ[i]/=n; accY[i]/=n; }
      return { yZ: accZ, yY: accY, alreadyDB: true };
    }else{
      accZ = new Array(freqs.length).fill(0);
      accY = new Array(freqs.length).fill(0);
      for(const [x,y] of pts){
        const Y = mobilityAtPoint(freqs, {...params, x0:x, y0:y});
        const Z = Y.map(cinv);
        for(let i=0;i<freqs.length;i++){ accZ[i]+=cabs(Z[i]); accY[i]+=cabs(Y[i]); }
      }
      for(let i=0;i<freqs.length;i++){ accZ[i]/=n; accY[i]/=n; }
      return { yZ: accZ, yY: accY, alreadyDB: false };
    }
  }

  function pullFromUI(){
    state.material = $("material").value;
    state.E   = parseFloat($("E").value)*1e9;
    state.rho = parseFloat($("rho").value);
    state.nu  = clamp(parseFloat($("nu").value), 0.0, 0.49);
    state.h   = parseFloat($("thick").value)/1000;
    state.Lx  = parseFloat($("Lx").value);
    state.Ly  = parseFloat($("Ly").value);
    state.x0mm = Math.max(0, parseFloat($("x0mm").value));
    state.y0mm = Math.max(0, parseFloat($("y0mm").value));
    state.Nx = Math.max(1, parseInt($("Nx").value));
    state.Ny = Math.max(1, parseInt($("Ny").value));
    state.eta = Math.max(0.0001, parseFloat($("eta").value)/100);
    state.fmin= Math.max(1, parseFloat($("fmin").value));
    state.fmax= Math.max(state.fmin+10, parseFloat($("fmax").value));
    state.mx  = Math.max(1, parseInt($("mx").value));
    state.ny  = Math.max(1, parseInt($("ny").value));
    state.showZ = $("showZ").checked;
    state.showY = $("showY").checked;
    state.logx  = $("logx").checked;
    state.db    = $("db").checked;
    state.bandBg= $("bandbg") ? $("bandbg").checked : true;
    state.ptmode = $("pt_grid").checked ? "grid" : "single";
    state.avgMode = $("avg_db").checked ? "db" : "linear";
    enableGridInputs();

    const d = D(state.E,state.nu,state.h);
    const mps = mp(state.rho,state.h);
    $("kpiD").textContent = d.toExponential(3);
    $("kpimm").textContent = mps.toFixed(1);
    const f11 = fmn(1,1,state.Lx,state.Ly,state.E,state.nu,state.rho,state.h);
    $("kpif11").textContent = f11.toFixed(1);

    plot();
  }

  function attachEvents(){
    $("material").addEventListener("change",(e)=>{
      const v = e.target.value;
      if(v!=="custom"){ applyPreset(v); }
    });
    ["E","rho","nu","thick","Lx","Ly","x0mm","y0mm","Nx","Ny","eta","fmin","fmax","mx","ny",
     "showZ","showY","logx","db","bandbg","pt_single","pt_grid","avg_linear","avg_db"].forEach(id=>{
      const el = $(id);
      if(el){ el.addEventListener("input", pullFromUI); el.addEventListener("change", pullFromUI); }
    });
    $("toggleOct").addEventListener("click", ()=>{ state.showOct = !state.showOct; plot(); });
    $("toggleInf").addEventListener("click", ()=>{ state.showInf = !state.showInf; plot(); });
    $("toggleModes").addEventListener("click", ()=>{ state.showModes = !state.showModes; plot(); });
    $("exportPng").addEventListener("click", ()=>{
      Plotly.downloadImage("plot",{format:"png", filename:"plate_impedance"});
    });
    $("reset").addEventListener("click", ()=>{
      $("material").value = "concrete"; applyPreset("concrete");
      $("Lx").value = 3.6; $("Ly").value = 2.7;
      $("x0mm").value = 1800; $("y0mm").value = 1350;
      $("Nx").value = 3; $("Ny").value = 3;
      $("pt_single").checked = true; $("pt_grid").checked = false;
      $("avg_linear").checked = true; $("avg_db").checked = false;
      $("eta").value = 1.0;
      $("fmin").value = 20; $("fmax").value = 5000;
      $("mx").value = 30; $("ny").value = 30;
      $("showZ").checked=true; $("showY").checked=false;
      $("logx").checked=true; $("db").checked=true;
      if($("bandbg")) $("bandbg").checked=true;
      state.showOct = true; state.showInf = true; state.showModes = true; state.bandBg = true;
      pullFromUI();
    });
  }

  function plot(){
    const N = 900;
    const fmin = state.fmin, fmax = state.fmax;
    const freqs = Array.from({length:N}, (_,i)=> fmin*Math.pow(fmax/fmin, i/(N-1)) );

    const params = {
      E:state.E, nu:state.nu, rho:state.rho, h:state.h,
      Lx:state.Lx, Ly:state.Ly,
      eta: state.eta, mx: state.mx, ny: state.ny
    };

    const avg = computeAveraged(freqs, params);
    let yZ = avg.yZ, yY = avg.yY;
    const alreadyDB = avg.alreadyDB;

    if(state.db && !alreadyDB){
      yZ = toDB20(yZ, 1.0);
      yY = toDB20(yY, 1.0);
    }
    // if !state.db && alreadyDB: 逆変換できないのでそのまま

    const traces = [];
    if(state.showZ){
      traces.push({x:freqs, y:yZ, type:"scatter", mode:"lines",
        name: state.db? "|Z| [dB re 1 Ns/m]" : "|Z| [N·s/m]"});
    }
    if(state.showY){
      traces.push({x:freqs, y:yY, type:"scatter", mode:"lines",
        name: state.db? "|Y| [dB re 1 m/N·s]" : "|Y| [m/N·s]"});
    }

    const primaryY = (state.showZ || !state.showY) ? yZ : yY;

    // ★ Octave means
    if(state.showOct){
      const xs=[], ys=[];
      for(const fc of OCT_CENTERS){
        const lo = fc/BAND_FACTOR, hi = fc*BAND_FACTOR;
        if(hi < fmin || lo > fmax) continue;
        let acc=0, cnt=0;
        for(let i=0;i<freqs.length;i++){
          const f = freqs[i];
          if(f>=lo && f<=hi){ acc += primaryY[i]; cnt++; }
        }
        if(cnt>0){ xs.push(fc); ys.push(acc/cnt); }
      }
      if(xs.length){
        traces.push({x:xs, y:ys, type:"scatter", mode:"markers", name:"Octave mean (★)",
          marker:{symbol:"star", size:10}});
      }
    }

    // ● Infinite-plate markers
    if(state.showInf){
      const zinf = Zinf_const(state.E,state.nu,state.rho,state.h);
      const val = state.db ? 20*Math.log10(zinf/1.0) : zinf;
      const xs=[], ys=[];
      for(const fc of OCT_CENTERS){
        if(fc < fmin || fc > fmax) continue;
        xs.push(fc); ys.push(val);
      }
      if(xs.length){
        traces.push({x:xs, y:ys, type:"scatter", mode:"markers", name:"Infinite plate (●)",
          marker:{symbol:"circle", size:8}});
      }
    }

    // ▲ Mode markers up to f32
    if(state.showModes){
      const modes = [[1,1],[1,2],[2,1],[2,2],[3,1],[3,2]];
      const xm=[], ym=[], txt=[];
      modes.forEach(([m,n])=>{
        const f = fmn(m,n, state.Lx,state.Ly, state.E,state.nu,state.rho,state.h);
        if(f>=fmin && f<=fmax){
          xm.push(f);
          ym.push(interpY(freqs, primaryY, f));
          txt.push(`f_${m}${n} = ${f.toFixed(1)} Hz`);
        }
      });
      if(xm.length){
        traces.push({x:xm, y:ym, type:"scatter", mode:"markers+text", name:"Modes (▲)",
          marker:{symbol:"triangle-up", size:9}, text:txt, textposition:"top center", hoverinfo:"text+x+y"});
      }
    }

    const layout = {
      title: {text: "駆動点インピーダンス / モビリティ（SSSS板, モード和）", font:{size:16}},
      xaxis: { title: "Frequency [Hz]", type: state.logx ? "log" : "linear",
               gridcolor: "#e5e7eb", zerolinecolor: "#e5e7eb" },
      yaxis: {
        title: state.db
          ? (state.showZ||!state.showY ? "Level [dB re 1 Ns/m]" : "Level [dB re 1 m/N·s]")
          : (state.showZ||!state.showY ? "Impedance |Z|" : "Mobility |Y|"),
        type: state.db ? "linear" : "log",
        gridcolor: "#e5e7eb", zerolinecolor: "#e5e7eb", rangemode: "tozero"
      },
      legend:{orientation:"h", x:0, y:1.06},
      plot_bgcolor:"#ffffff", paper_bgcolor:"#ffffff",
      margin:{l:80,r:30,t:60,b:60},
      shapes: state.bandBg ? makeBandShapes(fmin, fmax) : []
    };
    Plotly.react("plot", traces, layout, {responsive:true, displaylogo:false});
  }

  function init(){
    $("material").value = "concrete";
    applyPreset("concrete");
    enableGridInputs();
  }
  if(document.readyState === "loading"){
    document.addEventListener("DOMContentLoaded", ()=>{ attachEvents(); init(); });
  }else{
    attachEvents(); init();
  }
})();
