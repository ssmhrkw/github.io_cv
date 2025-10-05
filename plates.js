/* plates_v2.js — Rectangular plate DP impedance (live) with octave stars & mode markers */
(function(){
  "use strict";
  const $ = (id) => document.getElementById(id);
  const clamp = (v, lo, hi) => Math.min(hi, Math.max(lo, v));
  const EPS = 1e-30;

  // ---------- Complex helpers ----------
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
  const cabs = (z)=> Math.hypot(z.re, z.im);
  const creal = (z)=> z.re, cimag = (z)=> z.im;

  // ---------- State ----------
  const state = {
    material: "concrete",
    E: 30e9, rho: 2400, nu: 0.20, h: 0.150,
    Lx: 3.6, Ly: 2.7,
    x0mm: 1800, y0mm: 1350,
    eta: 0.01,
    fmin: 20, fmax: 5000,
    mx: 30, ny: 30,
    showZ: true, showY: false, logx: true, db: false,
    octType: 1, showOct: true, markN: 6, showModes: true,
  };

  const presets = {
    concrete:{E:30e9, rho:2400, nu:0.20, h_mm:150},
    gypsum:{E:3.0e9, rho:800, nu:0.30, h_mm:12.5},
    wood:{E:10e9, rho:600, nu:0.30, h_mm:18}
  };

  function applyPreset(name){
    const p = presets[name]; if(!p) return;
    $("E").value = (p.E/1e9).toFixed(1);
    $("rho").value = p.rho;
    $("nu").value = p.nu.toFixed(2);
    $("thick").value = p.h_mm;
    pullFromUI();
  }

  // ---------- Physics core ----------
  const D = (E,nu,h)=> E*h*h*h/(12*(1-nu*nu));       // [N*m]
  const mp = (rho,h)=> rho*h;                         // [kg/m2]

  function fmn(m,n,Lx,Ly,E,nu,rho,h){
    const lam2 = Math.PI*Math.PI*((m*m)/(Lx*Lx) + (n*n)/(Ly*Ly));
    const omg = Math.sqrt(D(E,nu,h)/mp(rho,h)) * (lam2*lam2);
    return omg/(2*Math.PI);
  }

  function mobilityAtPoint(freqs, params){
    const {E,nu,rho,h,Lx,Ly,x0,y0,eta,mx,ny} = params;
    const C = 4/(rho*h*Lx*Ly);
    const Y = new Array(freqs.length).fill(0).map(()=>c(0,0));
    for(let m=1; m<=mx; m++){
      const sx2 = Math.sin(m*Math.PI*x0/Lx); // will square later
      const mm = m*m/(Lx*Lx);
      for(let n=1; n<=ny; n++){
        const sy2 = Math.sin(n*Math.PI*y0/Ly);
        const nn = n*n/(Ly*Ly);
        const lam2 = Math.PI*Math.PI*(mm + nn);
        const omg_mn = Math.sqrt(D(E,nu,h)/mp(rho,h)) * (lam2*lam2);
        const phi2 = (sx2*sx2) * (sy2*sy2);
        for(let i=0;i<freqs.length;i++){
          const w = 2*Math.PI*freqs[i];
          const den = c( (omg_mn*omg_mn - w*w), (eta*omg_mn*w) );
          const num = c(0, w);
          const frac = cdiv(num, den);
          Y[i] = c( Y[i].re + C*phi2*frac.re, Y[i].im + C*phi2*frac.im );
        }
      }
    }
    return Y;
  }

  // ---------- Octave utilities ----------
  function octaveCenters(type, fmin, fmax){
    const list = [];
    if(type===1){
      // Start around 31.5 Hz and double
      let fc = 31.5;
      while(fc < fmin) fc *= 2;
      for(; fc <= fmax*1.01; fc*=2){ list.push(fc); }
    }else{
      // 1/3 octave based on 31.5 Hz nominal
      const k = Math.pow(2, 1/3);
      // start so that we land near within range
      let fc = 31.5;
      while(fc/k > fmin) fc /= k; // move down to just below fmin
      while(fc < fmin) fc *= k;
      for(; fc <= fmax*1.01; fc*=k){ list.push(fc); }
    }
    return list;
  }
  function bandEdges(fc, type){
    const factor = (type===1)? Math.SQRT2 : Math.pow(2, 1/6);
    return [fc/factor, fc*factor];
  }

  // ---------- Plot helpers ----------
  function toDB20(arr, ref){
    return arr.map(v => 20*Math.log10( Math.max(v, EPS) / ref ));
  }
  function interpY(xarr, yarr, x){
    // log-spaced data; find nearest two
    if(x<=xarr[0]) return yarr[0];
    if(x>=xarr[xarr.length-1]) return yarr[yarr.length-1];
    let lo=0, hi=xarr.length-1;
    while(hi-lo>1){
      const mid=(lo+hi)>>1;
      if(xarr[mid] > x) hi = mid; else lo = mid;
    }
    const x0=xarr[lo], x1=xarr[hi], y0=yarr[lo], y1=yarr[hi];
    const t = (Math.log(x) - Math.log(x0))/(Math.log(x1)-Math.log(x0));
    return y0*(1-t)+y1*t;
  }

  // ---------- UI sync ----------
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
    state.eta = Math.max(0.0001, parseFloat($("eta").value)/100);
    state.fmin= Math.max(1, parseFloat($("fmin").value));
    state.fmax= Math.max(state.fmin+10, parseFloat($("fmax").value));
    state.mx  = Math.max(1, parseInt($("mx").value));
    state.ny  = Math.max(1, parseInt($("ny").value));
    state.showZ = $("showZ").checked;
    state.showY = $("showY").checked;
    state.logx = $("logx").checked;
    state.db   = $("db").checked;
    state.octType = parseInt($("octType").value);
    state.markN = Math.max(2, Math.min(20, parseInt($("markN").value)||6));

    // KPIs
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
      if(v!=="custom"){ applyPreset(v); } else { /* keep manual */ }
    });
    ["E","rho","nu","thick","Lx","Ly","x0mm","y0mm","eta","fmin","fmax","mx","ny",
     "showZ","showY","logx","db","octType","markN"].forEach(id=>{
      $(id).addEventListener("input", pullFromUI);
      $(id).addEventListener("change", pullFromUI);
    });
    $("toggleOct").addEventListener("click", ()=>{ state.showOct = !state.showOct; plot(); });
    $("toggleModes").addEventListener("click", ()=>{ state.showModes = !state.showModes; plot(); });
    $("exportPng").addEventListener("click", ()=>{
      Plotly.downloadImage("plot",{format:"png", filename:"plate_impedance_v2"});
    });
    $("reset").addEventListener("click", ()=>{
      $("material").value = "concrete"; applyPreset("concrete");
      $("Lx").value = 3.6; $("Ly").value = 2.7;
      $("x0mm").value = 1800; $("y0mm").value = 1350;
      $("eta").value = 1.0;
      $("fmin").value = 20; $("fmax").value = 5000;
      $("mx").value = 30; $("ny").value = 30;
      $("showZ").checked=true; $("showY").checked=false;
      $("logx").checked=true; $("db").checked=false;
      $("octType").value = "1"; $("markN").value="6";
      state.showOct = true; state.showModes = true;
      pullFromUI();
    });
  }

  // ---------- Plot ----------
  function plot(){
    const N = 900;
    const fmin = state.fmin, fmax = state.fmax;
    const freqs = Array.from({length:N}, (_,i)=> fmin*Math.pow(fmax/fmin, i/(N-1)) );

    const params = {
      E:state.E, nu:state.nu, rho:state.rho, h:state.h,
      Lx:state.Lx, Ly:state.Ly,
      x0: state.x0mm/1000, y0: state.y0mm/1000,
      eta: state.eta, mx: state.mx, ny: state.ny
    };
    const Y = mobilityAtPoint(freqs, params);
    const Z = Y.map(cinv);

    // choose curve(s)
    let traces = [];
    const yRefZ = 1.0; // 1 Ns/m
    const yRefY = 1.0; // 1 m/Ns
    const magZ = Z.map(cabs);
    const magY = Y.map(cabs);
    let yZ = state.db ? toDB20(magZ, yRefZ) : magZ;
    let yY = state.db ? toDB20(magY, yRefY) : magY;

    if(state.showZ){
      traces.push({x:freqs, y:yZ, type:"scatter", mode:"lines", name: state.db? "|Z| [dB re 1 Ns/m]" : "|Z| [N·s/m]"});
    }
    if(state.showY){
      traces.push({x:freqs, y:yY, type:"scatter", mode:"lines", name: state.db? "|Y| [dB re 1 m/N·s]" : "|Y| [m/N·s]"});
    }

    // --- Octave stars (current primary curve: prefer Z if shown else Y) ---
    const primaryY = state.showZ ? yZ : yY;
    if(state.showOct){
      const centers = octaveCenters(state.octType, fmin, fmax);
      const xs=[], ys=[];
      centers.forEach(fc=>{
        const [lo,hi] = bandEdges(fc, state.octType);
        // sample indices inside band
        let acc=0, cnt=0;
        for(let i=0;i<freqs.length;i++){
          const f = freqs[i];
          if(f>=lo && f<=hi){ acc += primaryY[i]; cnt++; }
        }
        if(cnt>0){
          xs.push(fc);
          ys.push(acc/cnt);
        }
      });
      if(xs.length){
        traces.push({x:xs, y:ys, type:"scatter", mode:"markers", name:"Octave mean (★)",
          marker:{symbol:"star", size:10}});
      }
    }

    // --- Mode markers ▲ for first N modes ---
    if(state.showModes){
      // compute a set of lowest modes using a small search window
      const pairs = [];
      const M = Math.max(5, Math.min(20, Math.max(state.mx, state.ny)));
      for(let m=1;m<=M;m++){
        for(let n=1;n<=M;n++){
          const f = fmn(m,n, state.Lx,state.Ly, state.E,state.nu,state.rho,state.h);
          if(f>=fmin && f<=fmax){
            pairs.push({m,n,f});
          }
        }
      }
      pairs.sort((a,b)=> a.f-b.f);
      const take = pairs.slice(0, state.markN);
      const xM = [], yM = [], txt = [];
      take.forEach(p=>{
        xM.push(p.f);
        // place marker on primary curve for readability
        const yv = interpY(freqs, primaryY, p.f);
        yM.push(yv);
        txt.push(`f_${p.m}${p.n} = ${p.f.toFixed(1)} Hz`);
      });
      if(xM.length){
        traces.push({x:xM, y:yM, type:"scatter", mode:"markers+text", name:"Modes (▲)",
          marker:{symbol:"triangle-up", size:9}, text:txt, textposition:"top center", hoverinfo:"text+x+y"});
      }
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
        title: state.db
          ? (state.showZ||!state.showY ? "Level [dB re 1 Ns/m]" : "Level [dB re 1 m/N·s]")
          : (state.showZ||!state.showY ? "Impedance |Z|" : "Mobility |Y|"),
        type: state.db ? "linear" : "log",
        gridcolor: "#e5e7eb", zerolinecolor: "#e5e7eb", rangemode: "tozero"
      },
      legend:{orientation:"h", x:0, y:1.06},
      plot_bgcolor:"#ffffff", paper_bgcolor:"#ffffff",
      margin:{l:80,r:30,t:60,b:60}
    };
    Plotly.react("plot", traces, layout, {responsive:true, displaylogo:false});
  }

  // ---------- Init ----------
  function init(){
    // inflate inputs from HTML defaults
    $("material").value = "concrete";
    applyPreset("concrete");
  }
  if(document.readyState === "loading"){
    document.addEventListener("DOMContentLoaded", ()=>{ attachEvents(); init(); });
  }else{
    attachEvents(); init();
  }
})();
