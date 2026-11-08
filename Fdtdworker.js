// Fdtdworker.js — explicit FDTD for a simply supported Kirchhoff plate (SS on all edges).
// Rectangular grid Nx × Ny with spacings dx, dy, so the plate is exactly Lx × Ly.
// Force: 1 N rectangular pulse of 2 ms at (x0, y0), applied as a pressure 1/(dx·dy) on one cell.
// SS edges: w = 0 on the boundary and ∇²w = 0 there (the Laplacian buffer is never written on the boundary).
// Frames are stored every `saveEvery` steps; the result reports dt and Fs of the stored frames.
self.onmessage = function (ev) {
  const msg = ev.data; if (msg.type !== 'start') return;
  const c = msg.cfg;
  const E = c.E, nu = c.nu, rho = c.rho, h = c.h, Lx = c.Lx, Ly = c.Ly;
  const Nx = c.Nx, Ny = c.Ny, dx = c.dx, dy = c.dy, dt = c.dt;
  const steps = c.steps, saveEvery = Math.max(1, c.saveEvery | 0);

  const D = E * Math.pow(h, 3) / (12 * (1 - nu * nu));
  const m = rho * h;
  const size = Nx * Ny;
  let w = new Float64Array(size), wp = new Float64Array(size), wn = new Float64Array(size);
  const L1 = new Float64Array(size), L2 = new Float64Array(size);

  const sx = Math.max(1, Math.min(Nx - 2, Math.round(c.x0 / Lx * (Nx - 1))));
  const sy = Math.max(1, Math.min(Ny - 2, Math.round(c.y0 / Ly * (Ny - 1))));
  const src = sy * Nx + sx;
  const q = 1 / (dx * dy);                              // 1 N spread over one cell
  const forceSteps = Math.max(1, Math.round(0.002 / dt));
  const cx = 1 / (dx * dx), cy = 1 / (dy * dy);
  const coeff = dt * dt / m;

  function lap(a, out) {
    for (let j = 1; j < Ny - 1; j++) {
      const row = j * Nx;
      for (let i = 1; i < Nx - 1; i++) {
        const id = row + i, v2 = 2 * a[id];
        out[id] = (a[id - 1] - v2 + a[id + 1]) * cx + (a[id - Nx] - v2 + a[id + Nx]) * cy;
      }
    }
  }

  const frames = [];
  const report = Math.max(1, Math.floor(steps / 50));
  for (let step = 0; step < steps; step++) {
    lap(w, L1); lap(L1, L2);
    for (let j = 1; j < Ny - 1; j++) {
      const row = j * Nx;
      for (let i = 1; i < Nx - 1; i++) {
        const id = row + i;
        wn[id] = 2 * w[id] - wp[id] - coeff * D * L2[id];
      }
    }
    if (step < forceSteps) wn[src] += coeff * q;
    const t = wp; wp = w; w = wn; wn = t;             // rotate; boundary entries stay 0
    if ((step + 1) % saveEvery === 0) frames.push(Float32Array.from(w));
    if (step % report === 0) self.postMessage({type: 'progress', progress: Math.round(100 * step / steps), msg: `step ${step}/${steps}`});
  }

  const transferable = [];
  const framesBuf = frames.map(f => { transferable.push(f.buffer); return f.buffer; });
  const dtSave = dt * saveEvery;
  const result = {framesCount: frames.length, Nx, Ny, dt: dtSave, Fs: 1 / dtSave, dx, dy, Lx, Ly, framesBuf,
    xs: Array.from({length: Nx}, (_, i) => Lx * i / (Nx - 1)), ys: Array.from({length: Ny}, (_, j) => Ly * j / (Ny - 1))};
  self.postMessage({type: 'done', result}, transferable);
};
