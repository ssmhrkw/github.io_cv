/* Shared DSP. Runs in the page, in a Worker (TL processing) and in Node (unit tests). No DOM access. */
(function (root) {
'use strict';
const P0 = 2e-5;
const G10 = Math.pow(10, 0.3);                       // IEC 61260-1:2014 base-ten octave ratio

/* ---------- IIR sections ---------- */
function d1LP(fc, fs){ const K = Math.tan(Math.PI * fc / fs), a0 = 1 + K; return { b0: K / a0, b1: K / a0, b2: 0, a1: (K - 1) / a0, a2: 0 }; }
function d1HP(fc, fs){ const K = Math.tan(Math.PI * fc / fs), a0 = 1 + K; return { b0: 1 / a0, b1: -1 / a0, b2: 0, a1: (K - 1) / a0, a2: 0 }; }
function rbjLP(fc, fs, Q){
  const w0 = 2 * Math.PI * fc / fs, c = Math.cos(w0), s = Math.sin(w0), al = s / (2 * Q), a0 = 1 + al;
  return { b0: (1 - c) / 2 / a0, b1: (1 - c) / a0, b2: (1 - c) / 2 / a0, a1: -2 * c / a0, a2: (1 - al) / a0 };
}
function rbjHP(fc, fs, Q){
  const w0 = 2 * Math.PI * fc / fs, c = Math.cos(w0), s = Math.sin(w0), al = s / (2 * Q), a0 = 1 + al;
  return { b0: (1 + c) / 2 / a0, b1: -(1 + c) / a0, b2: (1 + c) / 2 / a0, a1: -2 * c / a0, a2: (1 - al) / a0 };
}
function butterQ(n, k){ return 1 / (2 * Math.sin((2 * k + 1) * Math.PI / (2 * n))); }
function butterLP(fc, fs, n){ const s = []; for (let k = 0; k < Math.floor(n / 2); k++) s.push(rbjLP(fc, fs, butterQ(n, k))); if (n % 2) s.push(d1LP(fc, fs)); return s; }
function butterHP(fc, fs, n){ const s = []; for (let k = 0; k < Math.floor(n / 2); k++) s.push(rbjHP(fc, fs, butterQ(n, k))); if (n % 2) s.push(d1HP(fc, fs)); return s; }
function sosMag(sections, fs, f){
  const w = 2 * Math.PI * f / fs, c1 = Math.cos(w), s1 = -Math.sin(w), c2 = Math.cos(2 * w), s2 = -Math.sin(2 * w);
  let Nr = 1, Ni = 0, Dr = 1, Di = 0;
  for (const q of sections){
    const nr = q.b0 + q.b1 * c1 + q.b2 * c2, ni = q.b1 * s1 + q.b2 * s2;
    const dr = 1 + q.a1 * c1 + q.a2 * c2, di = q.a1 * s1 + q.a2 * s2;
    let t = Nr * nr - Ni * ni; Ni = Nr * ni + Ni * nr; Nr = t;
    t = Dr * dr - Di * di; Di = Dr * di + Di * dr; Dr = t;
  }
  return Math.hypot(Nr, Ni) / Math.hypot(Dr, Di);
}
function normalizeAt(sections, fs, f){
  const h = sosMag(sections, fs, f);
  if (!(h > 0)) return sections;
  const g = Math.pow(1 / h, 1 / sections.length);
  for (const q of sections){ q.b0 *= g; q.b1 *= g; q.b2 *= g; }
  return sections;
}
/** Octave band-pass: Butterworth HP(f1) x LP(f2), order n each; edges fc*G^(-1/2), fc*G^(1/2). */
const BP_ORDER = 36;
function bandpassSOS(fc, fs, n){
  n = n || BP_ORDER;
  const g = Math.sqrt(G10);
  return normalizeAt(butterHP(fc / g, fs, n).concat(butterLP(fc * g, fs, n)), fs, fc);
}
/** A-weighting (IEC 61672-1:2013 Annex E poles), pre-warped bilinear first-order sections, 0 dB at 1 kHz. */
function aWeightSOS(fs){
  const f1 = 20.598997, f2 = 107.65265, f3 = 737.86223, f4 = 12194.217;
  const s = [d1HP(f1, fs), d1HP(f1, fs), d1HP(f2, fs), d1HP(f3, fs), d1LP(Math.min(f4, 0.45 * fs), fs), d1LP(Math.min(f4, 0.45 * fs), fs)];
  return normalizeAt(s, fs, 1000);
}
function sosFilter(x, sections){
  const y = Float64Array.from(x);
  for (const q of sections){
    let s1 = 0, s2 = 0;
    const b0 = q.b0, b1 = q.b1, b2 = q.b2, a1 = q.a1, a2 = q.a2;
    for (let n = 0; n < y.length; n++){
      const w = y[n] - a1 * s1 - a2 * s2;
      y[n] = b0 * w + b1 * s1 + b2 * s2;
      s2 = s1; s1 = w;
    }
  }
  return y;
}
function reverseInPlace(a){ for (let i = 0, j = a.length - 1; i < j; i++, j--){ const t = a[i]; a[i] = a[j]; a[j] = t; } return a; }
function filtfilt(x, sections, pad){
  const y = new Float64Array(x.length + 2 * pad);
  y.set(x, pad);
  let z = sosFilter(y, sections);
  z = sosFilter(reverseInPlace(z), sections);
  reverseInPlace(z);
  return z.subarray(pad, pad + x.length).slice();
}
/** Maximum of F time-weighted level (tau = 0.125 s), dB re 20 uPa. */
function fastMaxDb(x, fs){
  const al = 1 - Math.exp(-1 / (fs * 0.125));
  let s = 0, m = 0;
  for (let n = 0; n < x.length; n++){ s += al * (x[n] * x[n] - s); if (s > m) m = s; }
  return 10 * Math.log10(Math.max(m, 1e-30) / (P0 * P0));
}
/** A-weighting gain (linear) from IEC 61672-1:2013 eq. E.6 (A1000 = -2.00 dB). */
function aLin(f){
  const f2 = f * f;
  const ra = (12194.217 * 12194.217 * f2 * f2) /
    ((f2 + 20.598997 * 20.598997) * Math.sqrt((f2 + 107.65265 * 107.65265) * (f2 + 737.86223 * 737.86223)) * (f2 + 12194.217 * 12194.217));
  return ra * Math.pow(10, 2.0 / 20);
}
function aDb(f){ return 20 * Math.log10(aLin(f)); }

/* ---------- FFT ---------- */
function fftInPlace(re, im, inverse){
  const n = re.length;
  for (let i = 1, j = 0; i < n; i++){
    let bit = n >> 1;
    for (; j & bit; bit >>= 1) j ^= bit;
    j ^= bit;
    if (i < j){ let t = re[i]; re[i] = re[j]; re[j] = t; t = im[i]; im[i] = im[j]; im[j] = t; }
  }
  for (let len = 2; len <= n; len <<= 1){
    const ang = 2 * Math.PI / len * (inverse ? 1 : -1), wr = Math.cos(ang), wi = Math.sin(ang), h = len >> 1;
    for (let i = 0; i < n; i += len){
      let cr = 1, ci = 0;
      for (let j = 0; j < h; j++){
        const p = i + j, q = p + h;
        const br = re[q] * cr - im[q] * ci, bi = re[q] * ci + im[q] * cr;
        re[q] = re[p] - br; im[q] = im[p] - bi;
        re[p] += br; im[p] += bi;
        const t = cr * wr - ci * wi; ci = cr * wi + ci * wr; cr = t;
      }
    }
  }
  if (inverse){ for (let i = 0; i < n; i++){ re[i] /= n; im[i] /= n; } }
}
function mulberry32(a){
  return function(){
    a |= 0; a = a + 0x6D2B79F5 | 0;
    let t = Math.imul(a ^ a >>> 15, 1 | a);
    t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t;
    return ((t ^ t >>> 14) >>> 0) / 4294967296;
  };
}
/** Log-frequency raised-cosine edge: 0 below lo, 1 above hi. */
function edgeUp(f, lo, hi){
  if (f <= lo) return 0; if (f >= hi) return 1;
  const u = Math.log(f / lo) / Math.log(hi / lo);
  return 0.5 - 0.5 * Math.cos(Math.PI * u);
}
/** Largest |sample|. A NaN or infinite sample yields Infinity, so every level guard that uses the peak refuses it. */
function peakOf(a){ let p = 0; for (let i = 0; i < a.length; i++){ const v = Math.abs(a[i]); if (!(v <= 1e300)) return Infinity; if (v > p) p = v; } return p; }
function allFinite(a){ for (let i = 0; i < a.length; i++) if (!Number.isFinite(a[i])) return false; return true; }
function rmsOf(a){ let s = 0; for (let i = 0; i < a.length; i++) s += a[i] * a[i]; return Math.sqrt(s / Math.max(1, a.length)); }
/** Peak of the sum of copies of x started at the given sample offsets (exact overlap peak). */
function overlapPeak(x, offsets){
  if (offsets.length <= 1) return peakOf(x);
  const last = offsets[offsets.length - 1], y = new Float64Array(last + x.length);
  for (const o of offsets) for (let n = 0; n < x.length; n++) y[o + n] += x[n];
  return peakOf(y);
}

/* ---------- airborne insulation (TL / Dr / Rr) ---------- */
const TL_FREQ = [16, 31.5, 63, 125, 250, 500, 1000, 2000, 4000, 8000];   // nominal octave centres
/** Frequencies used for the curve, interpolation and band edges: exact octave steps from 1 kHz (15.625 ... 8000 Hz),
    so that each octave below 125 Hz is exactly 7.5 dB lower. Labels on screen use the nominal values above. */
const TL_FX = TL_FREQ.map((f, i) => 1000 * Math.pow(2, i - 6));
/** Dr / Rr grade curve. Defined 125 Hz - 2 kHz (grade g: 125 Hz = g-15, 500 Hz = g, 2 kHz = g+10, log-f interpolation).
    Outside that range the same slopes are continued (low) / held flat (high): listening-only extrapolation. */
function refCurve(f, g){ return f <= 500 ? g + 7.5 * Math.log2(f / 500) : f <= 2000 ? g + 5 * Math.log2(f / 500) : g + 10; }
/** Piecewise-linear over log f through (x[i], y[i]); held flat outside the end points. */
function interpLog(f, x, y){
  if (f <= x[0]) return y[0];
  for (let i = 1; i < x.length; i++) if (f <= x[i]) return y[i - 1] + (y[i] - y[i - 1]) * Math.log(f / x[i - 1]) / Math.log(x[i] / x[i - 1]);
  return y[y.length - 1];
}
/** 1 inside the octave bands lo..hi (band edges lo/sqrt2 .. hi*sqrt2), raised-cosine to 0 a quarter octave outside. */
function bandMask(f, lo, hi){
  if (!(f > 0)) return 0;
  const e1 = lo / Math.SQRT2, e2 = hi * Math.SQRT2, q = Math.pow(2, 0.25);
  return edgeUp(f, e1 / q, e1) * (1 - edgeUp(f, e2, e2 * q));
}
/** Apply attenuation curves to a mono signal with a zero-phase FFT gain.
    x/y: band centres and attenuation (dB) for offset 0; offsets: list of dB shifts; attenuation is floored at 0 dB.
    Returns the band-limited dry signal, one wet signal per offset, the overall peak and level statistics. */
function tlProcess(input, sr, x, y, offsets, lo, hi){
  let n = 1; while (n < input.length) n *= 2;
  const re = new Float64Array(n), im = new Float64Array(n);
  re.set(input); fftInPlace(re, im, false);
  const half = n >> 1, mask = new Float64Array(half + 1), att = new Float64Array(half + 1), aw = new Float64Array(half + 1);
  for (let k = 0; k <= half; k++){
    const f = k * sr / n;
    mask[k] = bandMask(f, lo, hi);
    att[k] = f > 0 ? interpLog(f, x, y) : 0;
    aw[k] = f > 0 ? aLin(f) : 0;
  }
  const bandOf = new Int16Array(half + 1).fill(-1);          // octave band (index into x) that each bin belongs to
  for (let k = 1; k <= half; k++){ const f = k * sr / n; for (let b = 0; b < x.length; b++) if (f >= x[b] / Math.SQRT2 && f < x[b] * Math.SQRT2){ bandOf[k] = b; break; } }
  const nfade = Math.max(1, Math.round(sr * 0.02)), L = input.length;
  /** Shape the spectrum, return the played buffer (L samples, 20 ms fades at both ends). */
  const render = (gainOf) => {
    const r = new Float64Array(n), i = new Float64Array(n);
    for (let k = 0; k <= half; k++){
      const g = mask[k] * gainOf(k);
      if (g === 0) continue;
      r[k] = re[k] * g; i[k] = im[k] * g;
      if (k > 0 && k < half){ r[n - k] = r[k]; i[n - k] = -i[k]; }
    }
    fftInPlace(r, i, true);
    const out = new Float32Array(L);
    let peak = 0;
    for (let m = 0; m < L; m++){
      const fade = Math.min(1, m / nfade, (L - 1 - m) / nfade);
      out[m] = r[m] * fade;
      const v = Math.abs(out[m]); if (!(v <= 1e300)) peak = Infinity; else if (v > peak) peak = v;
    }
    return { out, peak, lev: measure(out) };
  };
  /** Levels of a played buffer itself: overall power (Z), A-weighted power and octave-band powers. Zero-phase shaping
      spreads a short sound beyond its original length; what is cut off or faded there is not played, so it must not be
      counted. Measuring the finished buffer keeps the displayed reduction equal to the reduction that is heard. */
  function measure(out){
    const r = new Float64Array(n), i = new Float64Array(n), bp = new Float64Array(x.length);
    r.set(out); fftInPlace(r, i, false);
    let z = 0, a = 0;
    for (let k = 1; k <= half; k++){
      const p = (r[k] * r[k] + i[k] * i[k]) * (k === half ? 1 : 2);
      z += p; a += p * aw[k] * aw[k];
      if (bandOf[k] >= 0) bp[bandOf[k]] += p;
    }
    return { z, a, bands: Array.from(bp, v => 10 * Math.log10(Math.max(v, 1e-300))) };
  }
  const dry = render(() => 1);
  const variants = {}, stats = {}, bands = { dry: dry.lev.bands, variants: {} };
  let peak = dry.peak;
  for (const off of offsets){
    const w = render(k => Math.pow(10, -Math.max(0, att[k] + off) / 20));
    variants[String(off)] = w.out; bands.variants[String(off)] = w.lev.bands;
    peak = Math.max(peak, w.peak);
    stats[String(off)] = {
      dZ: 10 * Math.log10(Math.max(dry.lev.z, 1e-300) / Math.max(w.lev.z, 1e-300)),
      dA: 10 * Math.log10(Math.max(dry.lev.a, 1e-300) / Math.max(w.lev.a, 1e-300)),
      peakDb: 20 * Math.log10(Math.max(w.peak, 1e-30))
    };
  }
  return { dry: dry.out, variants, peak, dryPeak: dry.peak, stats, bands, empty: !(dry.lev.z > 0) };
}

const DSP = { P0, G10, d1LP, d1HP, butterLP, butterHP, sosMag, normalizeAt, bandpassSOS, aWeightSOS, sosFilter, filtfilt, fastMaxDb,
  aLin, aDb, fftInPlace, mulberry32, edgeUp, peakOf, allFinite, rmsOf, overlapPeak, TL_FREQ, TL_FX, refCurve, interpLog, bandMask, tlProcess };
root.DSP = DSP;
if (typeof module !== 'undefined' && module.exports) module.exports = DSP;
})(typeof self !== 'undefined' ? self : globalThis);
