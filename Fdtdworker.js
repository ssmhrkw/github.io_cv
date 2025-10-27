// Fdtdworker.js — FDTD本体（単純双ラプラシアン、SS境界）
self.onmessage = function(ev){
  const msg = ev.data;
  if(msg.type!=='start') return;
  const cfg = msg.cfg;

  const E=cfg.E, nu=cfg.nu, rho=cfg.rho, h=cfg.h;
  const Lx=cfg.Lx, Ly=cfg.Ly, Ngrid=cfg.Ngrid;
  const dx=cfg.dx, dt=cfg.dt, Fs=cfg.Fs;
  const framesCount=cfg.framesToSave;

  const D = E * Math.pow(h,3) / (12*(1-nu*nu));
  const marea = rho*h;
  const Nx=Ngrid, Ny=Ngrid, size=Nx*Ny;

  let w=new Float32Array(size), wprev=new Float32Array(size), wnext=new Float32Array(size);
  function idx(i,j){ return j*Nx + i; }

  const sx = Math.max(1, Math.min(Nx-2, Math.round(cfg.x0 / Lx * (Nx-1))));
  const sy = Math.max(1, Math.min(Ny-2, Math.round(cfg.y0 / Ly * (Ny-1))));
  const srcIdx = idx(sx, sy);

  const forceSteps = Math.max(1, Math.floor(Fs * 0.002)); // 2ms矩形

  function lap(arr){
    const out=new Float32Array(size);
    for(let j=1;j<Ny-1;j++){
      for(let i=1;i<Nx-1;i++){
        const id=idx(i,j);
        out[id]=(arr[id-1]+arr[id+1]+arr[id-Nx]+arr[id+Nx]-4*arr[id])/(dx*dx);
      }
    }
    return out;
  }
  function enforceSS(a){
    for(let i=0;i<Nx;i++){ a[idx(i,0)]=0; a[idx(i,Ny-1)]=0; }
    for(let j=0;j<Ny;j++){ a[idx(0,j)]=0; a[idx(Nx-1,j)]=0; }
  }
  const coeff = (dt*dt)/marea;

  const frames=[];
  for(let step=0; step<framesCount; step++){
    const F=new Float32Array(size);
    if(step < forceSteps) F[srcIdx]=1.0;

    const lap1=lap(w);
    const bih = lap(lap1);

    for(let id=0; id<size; id++){
      wnext[id] = 2*w[id] - wprev[id] + coeff * (F[id] - D * bih[id]);
    }
    enforceSS(wnext);
    wprev.set(w);
    w.set(wnext);
    frames.push(w.slice(0));

    if(step % Math.max(1, Math.floor(framesCount/40)) === 0){
      self.postMessage({type:'progress', progress: Math.round(100*step/framesCount), msg:`step ${step}/${framesCount}`});
    }
  }

  const transferable=[]; const framesBuf=frames.map(f=>{ transferable.push(f.buffer); return f.buffer; });
  const result={
    framesCount: frames.length,
    Nx, Ny, dt, Fs,
    framesBuf,
    // 軽量な補助（メインで再計算も可）
    modes: cfg.modesLight,
    xs: Array.from({length:Nx}, (_,i)=> Lx*i/(Nx-1)),
    ys: Array.from({length:Ny}, (_,j)=> Ly*j/(Ny-1))
  };
  self.postMessage({type:'done', result}, transferable);
};