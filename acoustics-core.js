/**
 * acoustics-core.js
 * Shared DSP and utilities for wav2Fmax, revtime_recorder, FilterDesign
 * ES2020 modules. Numbers for scalars, Float32Array for signals.
 * Coefficient shape per section: {b0,b1,b2,a1,a2}. First-order: b2=a2=0.
 * DF2T implementation. Normalize at fc => 0 dB.
 */
"use strict";

/** First-order via pre-warped bilinear */
export function d1_LP(fc, fs){
  const K = Math.tan(Math.PI*fc/fs);
  const a0 = 1 + K;
  return {b0:K/a0, b1:K/a0, b2:0, a1:(1-K)/a0, a2:0};
}
export function d1_HP(fc, fs){
  const K = Math.tan(Math.PI*fc/fs);
  const a0 = 1 + K;
  return {b0:1/a0, b1:-1/a0, b2:0, a1:(1-K)/a0, a2:0};
}

/** RBJ biquads with pre-warped w0 */
export function rbjLP(fc, fs, Q){
  const w0 = 2*Math.atan(Math.tan(Math.PI*fc/fs));
  const c = Math.cos(w0), s = Math.sin(w0), alpha = s/(2*Q);
  const a0 = 1+alpha, a1 = -2*c, a2 = 1-alpha;
  const b0 = (1-c)/2, b1 = 1-c, b2 = (1-c)/2;
  return {b0:b0/a0, b1:b1/a0, b2:b2/a0, a1:a1/a0, a2:a2/a0};
}
export function rbjHP(fc, fs, Q){
  const w0 = 2*Math.atan(Math.tan(Math.PI*fc/fs));
  const c = Math.cos(w0), s = Math.sin(w0), alpha = s/(2*Q);
  const a0 = 1+alpha, a1 = -2*c, a2 = 1-alpha;
  const b0 = (1+c)/2, b1 = -(1+c), b2 = (1+c)/2;
  return {b0:b0/a0, b1:b1/a0, b2:b2/a0, a1:a1/a0, a2:a2/a0};
}

/** Butterworth helpers */
export function butterQ(n,k){ return 1/(2*Math.sin((2*k+1)*Math.PI/(2*n))); }
export function butterLP_sos(fc, fs, order){
  if(order<=0) return [];
  const sos=[], p=Math.floor(order/2);
  for(let k=0;k<p;k++) sos.push(rbjLP(fc,fs,butterQ(order,k)));
  if(order%2===1) sos.push(d1_LP(fc,fs));
  return sos;
}
export function butterHP_sos(fc, fs, order){
  if(order<=0) return [];
  const sos=[], p=Math.floor(order/2);
  for(let k=0;k<p;k++) sos.push(rbjHP(fc,fs,butterQ(order,k)));
  if(order%2===1) sos.push(d1_HP(fc,fs));
  return sos;
}

/** SOS frequency response magnitude on fgrid */
export function sosFreqz(sections, fs, fgrid){
  const H = new Float64Array(fgrid.length);
  for(let i=0;i<fgrid.length;i++){
    const w = 2*Math.PI*fgrid[i]/fs;
    const c1 = Math.cos(-w), s1 = Math.sin(-w);
    const c2 = Math.cos(-2*w), s2 = Math.sin(-2*w);
    let Nr=1, Ni=0, Dr=1, Di=0;
    for(const s of sections){
      const nre = s.b0 + s.b1*c1 + s.b2*c2;
      const nim = s.b1*s1 + s.b2*s2;
      const dre = 1   + s.a1*c1 + s.a2*c2;
      const dim =       s.a1*s1 + s.a2*s2;
      const tNr = Nr*nre - Ni*nim, tNi = Nr*nim + Ni*nre; Nr=tNr; Ni=tNi;
      const tDr = Dr*dre - Di*dim, tDi = Dr*dim + Di*dre; Dr=tDr; Di=tDi;
    }
    const den = Dr*Dr + Di*Di;
    const Hr = (Nr*Dr + Ni*Di)/den;
    const Hi = (Ni*Dr - Nr*Di)/den;
    H[i] = Math.hypot(Hr, Hi);
  }
  return H;
}

/** Normalize cascade gain at fref to 0 dB (unity magnitude) */
export function normalizeAt(sections, fs, fref){
  const probe = [0.98*fref, fref, 1.02*fref];
  const mags = sosFreqz(sections, fs, probe);
  const Href = mags[1] || 1;
  if(!Number.isFinite(Href) || Href<=0) return;
  const K = Math.max(1, sections.length);
  const g = Math.pow(1/(Href), 1/K);
  for(const s of sections){ s.b0*=g; s.b1*=g; s.b2*=g; }
}

/** Butterworth band-pass via HP(f1)*LP(f2) split. N => split to ⌊N/2⌋ and ⌈N/2⌉. */
export function designButterworthBandpassSOS_N(fc, fs, N){
  const f1 = fc/Math.SQRT2, f2 = fc*Math.SQRT2;
  const M  = Math.max(2, 2*Math.max(1, Math.round(N)));
  const hpOrd = Math.floor(M/2), lpOrd = Math.ceil(M/2);
  const hp = butterHP_sos(f1, fs, hpOrd), lp = butterLP_sos(f2, fs, lpOrd);
  const sections = hp.concat(lp);
  normalizeAt(sections, fs, fc);
  return sections;
}

/** DF2T single biquad */
function _df2t_one(x, c){
  const N=x.length, y=new Float32Array(N);
  let s1=0, s2=0;
  const b0=c.b0, b1=c.b1, b2=c.b2, a1=c.a1, a2=c.a2;
  for(let n=0;n<N;n++){
    const w = x[n] - a1*s1 - a2*s2;
    const yn = b0*w + b1*s1 + b2*s2;
    s2 = s1; s1 = w;
    y[n] = yn;
  }
  return y;
}

/** Apply SOS cascade (DF2T). */
export function sosFilter(xFloat32, sosArray){
  let y = Float32Array.from(xFloat32);
  for(let i=0;i<sosArray.length;i++) y = _df2t_one(y, sosArray[i]);
  return y;
}

/** Band centers */
export function centersOct(){ const a=[]; let f=16; while(f<=8000){ a.push(f); f*=2; } return a; }
export function centersThird(){
  const a=[]; let f=16; const step=Math.pow(2,1/3);
  while(f<=8000){ a.push(Number(f.toFixed(6))); f*=step; }
  return a;
}
export function fullCenters(type){ return type==='oct' ? centersOct() : centersThird(); }

/** Analysis helpers */
const P0 = 2e-5;
export function percentile(arr, p){
  if(!arr.length) return NaN;
  const a = Float64Array.from(arr).sort();
  const pos = (a.length-1)*p;
  const lo = Math.floor(pos), hi = Math.ceil(pos), w = pos-lo;
  return (1-w)*a[lo] + w*a[hi];
}
export function ewmaSquaredToSPL(y, sr, tau){
  const out = new Float32Array(y.length);
  const alpha = 1 - Math.exp(-1/(sr*tau));
  let s=0;
  for(let i=0;i<y.length;i++){
    const v=y[i];
    s += alpha*(v*v - s);
    out[i] = 10*Math.log10( Math.max(1e-24, s/(P0*P0)) );
  }
  return out;
}
export function movingLeqSPL(y, sr, winSec){
  const win = Math.max(1, Math.round(winSec*sr));
  const out=new Float32Array(y.length);
  let acc=0;
  for(let i=0;i<y.length;i++){
    const v=y[i], vv=v*v;
    acc += vv;
    if(i>=win){ const old=y[i-win]; acc -= old*old; }
    const denom = P0*P0*Math.min(win, i+1)/sr;
    out[i] = 10*Math.log10( Math.max(1e-24, acc/denom) );
  }
  return out;
}
export function downsampleForPlot(x,y,maxPoints){
  if(y.length<=maxPoints) return {x,y};
  const step=y.length/maxPoints, xn=new Array(maxPoints), yn=new Array(maxPoints);
  for(let i=0;i<maxPoints;i++){ const k=Math.floor(i*step); xn[i]=x[k]; yn[i]=y[k]; }
  return {x:xn,y:yn};
}

/** Internal hooks for tests */
export const _internal = { _df2t_one };
