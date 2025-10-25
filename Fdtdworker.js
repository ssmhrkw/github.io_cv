// worker.js
self.onmessage = function(ev){
  const msg = ev.data;
  if(msg.type!=='start') return;
  const cfg = msg.cfg;
  const E = cfg.E, nu = cfg.nu, rho = cfg.rho, h = cfg.h;
  const Lx = cfg.Lx, Ly = cfg.Ly, Ngrid = cfg.Ngrid;
  const dx = cfg.dx, dt = cfg.dt;
  const framesCount = cfg.framesToSave;
  const Fs = cfg.Fs;
  const marea = rho*h;
  const D = E * Math.pow(h,3) / (12*(1-nu*nu));
  const Nx = Ngrid, Ny = Ngrid;
  const size = Nx*Ny;
  let w = new Float32Array(size), wprev = new Float32Array(size), wnext = new Float32Array(size);
  function idx(i,j){ return j*Nx + i; }
  const x0 = cfg.x0, y0 = cfg.y0;
  const sx = Math.max(1, Math.min(Nx-2, Math.round(x0 / cfg.Lx * (Nx-1))));
  const sy = Math.max(1, Math.min(Ny-2, Math.round(y0 / cfg.Ly * (Ny-1))));
  const srcIdx = idx(sx, sy);
  const forceSteps = Math.max(1, Math.floor(Fs * 0.002));
  const totalSteps = framesCount;
  function lap(arr){
    const out = new Float32Array(size);
    for(let j=1;j<Ny-1;j++){
      for(let i=1;i<Nx-1;i++){
        const id = idx(i,j);
        out[id] = (arr[id-1] + arr[id+1] + arr[id-Nx] + arr[id+Nx] - 4*arr[id])/(dx*dx);
      }
    }
    return out;
  }
  function enforceSS(a){
    for(let i=0;i<Nx;i++){ a[idx(i,0)] = 0; a[idx(i,Ny-1)]=0; }
    for(let j=0;j<Ny;j++){ a[idx(0,j)] = 0; a[idx(Nx-1,j)] = 0; }
  }
  const coeff = (dt*dt)/marea;
  const frames = [];
  for(let step=0; step<totalSteps; step++){
    const F = new Float32Array(size);
    if(step < forceSteps) F[srcIdx] = 1.0;
    const lap1 = lap(w);
    const bih = lap(lap1);
    for(let id=0; id<size; id++){
      wnext[id] = 2*w[id] - wprev[id] + coeff * (F[id] - D * bih[id]);
    }
    enforceSS(wnext);
    wprev.set(w);
    w.set(wnext);
    frames.push(w.slice(0));
    if(step % Math.max(1, Math.floor(totalSteps/40)) === 0){
      self.postMessage({type:'progress', progress: Math.round(100*step/totalSteps), msg: 'step '+step+'/'+totalSteps});
    }
  }
  const modes = cfg.modesLight;
  const Mmn = rho*h * (cfg.Lx * cfg.Ly) / 4.0;
  const omegaForced = 2*Math.PI * (modes[0].f || 1.0);
  const W = new Float32Array(size);
  for(const md of modes){
    const omega_m = md.omega;
    const phi0 = Math.sin(md.m*Math.PI*cfg.x0/cfg.Lx) * Math.sin(md.n*Math.PI*cfg.y0/cfg.Ly);
    const real = Mmn*(omega_m*omega_m - omegaForced*omegaForced);
    const imag = Mmn*(omegaForced * omega_m / 1e8);
    const denom = real*real + imag*imag;
    const er = (1.0*phi0*real)/denom;
    const ei = (-1.0*phi0*imag)/denom;
    for(let j=0;j<Ny;j++){
      const y = j/(Ny-1)*cfg.Ly;
      for(let i=0;i<Nx;i++){
        const x = i/(Nx-1)*cfg.Lx;
        const idx0 = idx(i,j);
        const cx = Math.sin(md.m * Math.PI * x / cfg.Lx);
        const cy = Math.sin(md.n * Math.PI * y / cfg.Ly);
        W[idx0] += Math.hypot(er*cx*cy, ei*cx*cy);
      }
    }
  }
  const transferable = [];
  const framesBuf = frames.map(f=>{ transferable.push(f.buffer); return f.buffer; });
  const result = {framesCount: frames.length, Nx: Nx, Ny: Ny, dt: dt, Fs: Fs, framesBuf: framesBuf, staticZ: Array.from(W), modes: modes, xs: Array.from({length:Nx}, (_,i)=>cfg.Lx*i/(Nx-1)), ys: Array.from({length:Ny}, (_,j)=>cfg.Ly*j/(Ny-1))};
  self.postMessage({type:'done', result: result}, transferable);
};
