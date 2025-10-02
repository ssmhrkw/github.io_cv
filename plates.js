// plates_worker.js
function clamp(x, lo, hi){ return Math.min(Math.max(x, lo), hi); }
function lossFactor(f){ const eta = 0.005 + 0.3/Math.sqrt(Math.max(1, f)); return clamp(eta, 1e-4, 0.2); }

function buildModalCache({Lx, Ly, h, E, rho, nu, px, py, P=12, Q=12}){
  const rho_s = rho*h;
  const D = (E*h**3)/(12*(1-nu**2));
  const S = Lx*Ly;
  const N = P*Q, psi2=new Float64Array(N), w2=new Float64Array(N);
  let k=0;
  for (let p=1;p<=P;p++){
    const sp = Math.sin(p*Math.PI*px/Lx), a=(p/Lx)**2;
    for (let q=1;q<=Q;q++){
      const sq = Math.sin(q*Math.PI*py/Ly), b=(q/Ly)**2;
      psi2[k]=(sp*sq)*(sp*sq);
      const f_pq=(Math.PI/2)*Math.sqrt(D/rho_s)*(a+b);
      w2[k]=(2*Math.PI*f_pq)**2; k++;
    }
  }
  const idx=[...Array(N).keys()].sort((i,j)=>w2[i]-w2[j]);
  const psi2s=new Float64Array(N), w2s=new Float64Array(N);
  for (let i=0;i<N;i++){ psi2s[i]=psi2[idx[i]]; w2s[i]=w2[idx[i]]; }
  return { psi2:psi2s, w2:w2s, rho_s, S };
}

function mobilityFromCache(omega, eta, cache, eps=1e-10){
  const { psi2, w2, rho_s, S } = cache;
  let Yr=0, Yi=0, accum=0;
  for (let k=0;k<psi2.length;k++){
    const w2k=w2[k], Ar=w2k-omega*omega, Ai=w2k*eta;
    const denom=Ar*Ar+Ai*Ai, scale=psi2[k]/denom;
    const termAbs = scale*Math.hypot(Ar,Ai);
    accum+=termAbs; Yr+=scale*Ar; Yi-=scale*Ai;
    if (termAbs<eps && accum>100*eps) break;
  }
  const C = 4/(rho_s*S);
  return { re: -omega*C*Yi, im: omega*C*Yr };
}

let cache=null;

self.onmessage = (e)=>{
  const { type, params, freqs } = e.data || {};
  if (type==='build'){
    cache = buildModalCache(params);
    self.postMessage({ type:'built' });
  } else if (type==='sweep' && cache){
    const out=[];
    for (const f of freqs){
      const omega=2*Math.PI*f, eta=lossFactor(f);
      const Y=mobilityFromCache(omega, eta, cache);
      const d = Y.re*Y.re + Y.im*Y.im;
      if (d===0){ out.push({f, Zr:null, Zi:null, ZdB:null, phaseDeg:null}); continue; }
      const Zr=Y.re/d, Zi=-Y.im/d;
      const ZdB = 20*Math.log10(Math.max(Math.hypot(Zr, Zi), 1e-12));
      const phaseDeg = Math.atan2(Zi, Zr)*180/Math.PI;
      out.push({ f, Zr, Zi, ZdB, phaseDeg });
    }
    self.postMessage({ type:'impedance', data: out });
  }
};
