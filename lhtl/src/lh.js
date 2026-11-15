/* LH: heavy floor impact sound listening. */
App.register('lh', (function () {
'use strict';
const D = DSP, P0 = D.P0, G10 = D.G10, el = App.el, fmt = App.fmt, sgn = App.sgn, db = App.db;

/* ---------- constants ---------- */
const KEY = 'lh-listening-v5';
const OLD_KEYS = ['lh-listening-v4-range', 'lh-listening-v3-2khz', 'lh-listening-v2', 'lh-listening-v1'];
const BANDS = [
  { label: '16', x: -6, modes: 8, add: true }, { label: '31.5', x: -5, modes: 8, add: true },
  { label: '63', x: -4, modes: 8 }, { label: '125', x: -3, modes: 12 }, { label: '250', x: -2, modes: 20 }, { label: '500', x: -1, modes: 30 },
  { label: '1k', x: 0, modes: 40, add: true }, { label: '2k', x: 1, modes: 50, add: true }
];
BANDS.forEach((b, i) => { b.i = i; b.fc = 1000 * Math.pow(G10, b.x); });          // exact (base-ten) centre frequency
const ADD_LABELS = BANDS.filter(b => b.add).map(b => b.label);
/** LH-65 reference: C(f) = round(62 - A(f_exact)); reproduces 88/78/71/65/62 and extends it with the same rounding. */
const REF65 = BANDS.map(b => Math.round(62 - D.aDb(b.fc)));
const GRADES = [30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80];
const BG_LEVELS = [null, 20, 30, 40, 50, 60, 70];
const K_REF = 126;                 // uncalibrated: model dB (peak) that maps to 0 dBFS at output level 0 dB
const CAL_DBFS = -20;
const CEIL = Math.pow(10, -1 / 20); // digital peak ceiling, -1 dBFS
const LEVEL_MIN = -40, LEVEL_MAX = 20, LEVEL_STEP = 5;
const QUIET_DBFS = -90;            // below this rms a signal is at the 16-bit quantisation floor
const TRIPLE_OFFSETS = [0, 1, 2], RUN_INTERVAL = 1.9, RUN_LOOKAHEAD = 3, RUN_TICK_MS = 250;
const LEAD = 0.8, FADE = 0.05, GAP = 0.08, BG_FADE_IN = 0.15;
const PATTERNS = ['single', 'triple', 'run'];
const TYPING_MS = 500;            // rebuild delay while a number is being typed

/* ---------- settings ---------- */
const isNum = v => typeof v === 'number' && Number.isFinite(v);
function defaults(){
  const ref50 = {}, off = {};
  BANDS.forEach((b, i) => { ref50[b.label] = REF65[i] - 15; off[b.label] = 0; });
  return { v: 5, mode: 'grade', add: { '16': false, '31.5': false, '1k': false, '2k': false }, custom4: [88, 78, 71, 65], ref50, off,
    pattern: 'single', T: 0.7, slope: -5, bgBand: 'wide', nums: false, level: 0, cal: {} };
}
function sanitize(s){
  const o = defaults();
  if (!s || typeof s !== 'object') return o;
  if (s.mode === 'custom' || s.mode === 'grade') o.mode = s.mode;
  if (s.add && typeof s.add === 'object') for (const k of ADD_LABELS) o.add[k] = s.add[k] === true;
  if (Array.isArray(s.custom4) && s.custom4.length === 4 && s.custom4.every(v => isNum(v) && v >= 0 && v <= 140)) o.custom4 = s.custom4.slice();
  for (const b of BANDS){
    if (s.ref50 && isNum(s.ref50[b.label]) && s.ref50[b.label] >= 0 && s.ref50[b.label] <= 140) o.ref50[b.label] = s.ref50[b.label];
    if (s.off && isNum(s.off[b.label]) && s.off[b.label] <= 0 && s.off[b.label] >= -60) o.off[b.label] = s.off[b.label];
    if (s.cal && isNum(s.cal[b.label]) && s.cal[b.label] >= 0 && s.cal[b.label] <= 140) o.cal[b.label] = s.cal[b.label];
  }
  if (PATTERNS.includes(s.pattern)) o.pattern = s.pattern;
  if (isNum(s.T)) o.T = Math.min(2, Math.max(0.2, s.T));
  if ([-3, -5, -8].includes(s.slope)) o.slope = s.slope;
  if (s.bgBand === 'b4' || s.bgBand === 'wide') o.bgBand = s.bgBand;
  if (typeof s.nums === 'boolean') o.nums = s.nums;
  if (isNum(s.level)) o.level = Math.min(LEVEL_MAX, Math.max(LEVEL_MIN, Math.round(s.level / LEVEL_STEP) * LEVEL_STEP));
  if (o.mode === 'grade' && o.add['16'] && o.T < 0.5) o.T = 0.5;
  return o;
}
/** Convert settings saved by earlier versions (v1 - v4). Old keys are only read, never deleted. */
function migrate(old){
  const s = defaults();
  if (!old || typeof old !== 'object') return s;
  const range = old.range;
  s.mode = range === 'custom4' ? 'custom' : 'grade';
  if (range === 'custom4') s.bgBand = 'b4';
  s.add['1k'] = range !== 'custom4';                       // v1 - v4 always included 1 kHz in the grade modes
  s.add['31.5'] = range === 'extended' || range === 'wide';
  s.add['16'] = range === 'wide';
  if (Array.isArray(old.custom4)) s.custom4 = old.custom4;
  ['63', '125', '250', '500'].forEach((lab, i) => {
    if (Array.isArray(old.curve) && isNum(old.curve[i])) s.ref50[lab] = old.curve[i];
    if (Array.isArray(old.off) && isNum(old.off[i])) s.off[lab] = Math.min(0, old.off[i]);
  });
  if (Array.isArray(old.ext) && isNum(old.ext[0]) && !(old.ext[0] === -10 && old.ext[1] === -20)){
    s.ref50['1k'] = s.ref50['500'] + old.ext[0];
    if (isNum(old.ext[1])) s.ref50['2k'] = s.ref50['500'] + old.ext[1];
  }
  s.pattern = old.pattern; s.T = old.T; s.slope = old.slope; s.nums = old.nums;
  if (old.calMap && typeof old.calMap === 'object') s.cal = Object.assign({}, old.calMap);
  else if (Array.isArray(old.cal)) ['63', '125', '250', '500', '1k', '2k'].forEach((lab, i) => { if (isNum(old.cal[i])) s.cal[lab] = old.cal[i]; });
  return s;                                                // outputMode (compare / audition) no longer exists
}
function loadSettings(){
  const cur = App.store.get(KEY);
  if (cur) return sanitize(cur);
  for (const k of OLD_KEYS){ const old = App.store.get(k); if (old) return sanitize(migrate(old)); }
  return defaults();
}
let S = defaults();
function save(){ App.store.set(KEY, S); }

/* ---------- derived state ---------- */
let source = 'synth', wav = null;
let impact = null, bg = null, buildError = null, dirty = true;
let selected = { G: 50, LA: 30 };
let playing = null;                // { G, LA } while something sounds

function activeBands(){ return BANDS.filter(b => S.mode === 'custom' ? !b.add : (!b.add || S.add[b.label])); }
function calibrated(){ return BANDS.some(b => isNum(S.cal[b.label])); }
/** K (dB SPL at 0 dBFS) for all eight bands; measured bands exactly, others interpolated over band index, ends held. */
function calKAll(){
  const idx = BANDS.filter(b => isNum(S.cal[b.label])).map(b => b.i);
  if (!idx.length) return BANDS.map(() => K_REF);
  return BANDS.map(b => {
    const m = j => S.cal[BANDS[j].label] - CAL_DBFS;
    if (idx.includes(b.i)) return m(b.i);
    let lo = -1, hi = -1;
    for (const j of idx){ if (j < b.i) lo = j; if (j > b.i && hi < 0) hi = j; }
    if (lo < 0) return m(hi);
    if (hi < 0) return m(lo);
    return m(lo) + (m(hi) - m(lo)) * (b.i - lo) / (hi - lo);
  });
}
function kAtF(K, f){
  const x = Math.log(f / BANDS[0].fc) / Math.log(G10);
  if (x <= 0) return K[0];
  if (x >= K.length - 1) return K[K.length - 1];
  const i = Math.floor(x), u = x - i;
  return K[i] + (K[i + 1] - K[i]) * u;
}
function ref50(b){ return S.mode === 'custom' ? S.custom4[b.i - 2] - 15 : S.ref50[b.label]; }
function target50(b){ return S.mode === 'custom' ? ref50(b) : S.ref50[b.label] + S.off[b.label]; }
function gradeName(G){ return S.mode === 'custom' ? (G === 65 ? '指定値' : '指定値 ' + sgn(G - 65) + ' dB') : 'LH-' + G; }
function cellName(G, LA){
  return G === null ? '暗騒音のみ ' + LA + ' dB(A)' : gradeName(G) + (LA === null ? '（暗騒音なし）' : ' ＋ 暗騒音 ' + LA + ' dB(A)');
}
const keyOf = (G, LA) => (G === null ? 'bg' : G) + '|' + (LA === null ? 'none' : LA);

/* ---------- impact sound ---------- */
const rawCache = new Map(), sosCache = new Map();
function bandSOS(b, fs){
  const k = b.label + ':' + fs;
  if (!sosCache.has(k)) sosCache.set(k, D.bandpassSOS(b.fc, fs));
  return sosCache.get(k);
}
/** One band of the modal impact model. The random seed depends only on the band, so adding a band never changes the others. */
function synthBand(b, fs, T, N){
  const k = b.label + ':' + fs + ':' + T;
  if (rawCache.has(k)) return rawCache.get(k);
  const rnd = D.mulberry32(20261109 + 7919 * b.i), decay = Math.log(1000) / T;
  const r = Math.exp(-decay / fs), na = 0.004 * fs, nf = Math.round(0.05 * fs);
  let x = new Float64Array(N);
  for (let m = 0; m < b.modes; m++){
    const f = b.fc * Math.pow(G10, (0.1 + 0.8 * rnd()) - 0.5);
    const amp = 0.6 + 0.8 * rnd(), ph = 2 * Math.PI * rnd(), w = 2 * Math.PI * f / fs;
    const cr = r * Math.cos(w), ci = r * Math.sin(w);
    let re = amp * Math.cos(ph), im = amp * Math.sin(ph);
    for (let n = 0; n < N; n++){ x[n] += im; const t = re * cr - im * ci; im = re * ci + im * cr; re = t; }
  }
  for (let n = 0; n < N; n++){
    let e = 1 - Math.exp(-n / na);
    if (n > N - nf) e *= 0.5 - 0.5 * Math.cos(Math.PI * (N - n) / nf);
    x[n] *= e;
  }
  if (b.fc < 50){   // keep the very-low-frequency attack from spilling into higher bands
    const z = D.filtfilt(x, bandSOS(b, fs), Math.round(fs));
    const edge = Math.min(Math.round(fs * 0.2), Math.floor(N / 3));
    for (let n = 0; n < edge; n++){ const w = 0.5 - 0.5 * Math.cos(Math.PI * n / edge); z[n] *= w; z[N - 1 - n] *= w; }
    x = z;
  }
  if (rawCache.size > 24) rawCache.clear();
  rawCache.set(k, x);
  return x;
}
function synthBands(fs, T, bands){
  if (bands[0].fc < 20 && T < 0.5) throw new Error('16 Hz 帯を含む合成音では T を 0.5 秒以上にしてください。');
  const dur = Math.min(4, 0.15 + 2.2 * T), N = Math.round(dur * fs);
  return bands.map(b => synthBand(b, fs, T, N));
}
const tailCache = new Map();
/** Samples to analyse after the end of a signal so that the octave filter and the F time weighting have passed their
    maximum: the time at which the F-weighted impulse response of the band filter has fallen 10 dB below its own maximum.
    (A 16 Hz octave filter keeps rising for several tenths of a second after a short signal has already ended.) */
function bandTail(b, fs){
  const k = b.label + ':' + fs;
  if (tailCache.has(k)) return tailCache.get(k);
  const N = Math.round(fs * 6), x = new Float64Array(N);
  x[0] = 1;
  const y = D.sosFilter(x, bandSOS(b, fs)), al = 1 - Math.exp(-1 / (fs * 0.125)), lev = new Float64Array(N);
  let acc = 0, m = 0;
  for (let n = 0; n < N; n++){ acc += al * (y[n] * y[n] - acc); lev[n] = acc; if (acc > m) m = acc; }
  let last = N - 1;
  while (last > 0 && lev[last] < m / 10) last--;
  const t = last + Math.round(0.05 * fs);
  tailCache.set(k, t);
  return t;
}
/** Scale band components so that the 1/1-octave L_Fmax of their sum hits the targets (Pa domain).
    The analysis runs on the sum followed by silence, so the level is the true maximum, not the value at the last sample. */
function fitBands(raw, targets, fs, bands){
  const nb = raw.length, N = raw[0].length;
  const sos = bands.map(b => bandSOS(b, fs)), tails = bands.map(b => bandTail(b, fs)), sosA = D.aWeightSOS(fs);
  const M = N + Math.max(Math.round(0.5 * fs), ...tails);
  const g = raw.map((x, b) => Math.pow(10, (targets[b] - D.fastMaxDb(x, fs)) / 20));
  const sum = new Float64Array(M);
  const mix = () => { sum.fill(0); for (let b = 0; b < nb; b++){ const gb = g[b], x = raw[b]; for (let n = 0; n < N; n++) sum[n] += gb * x[n]; } };
  const analyze = () => sos.map((q, b) => D.fastMaxDb(D.sosFilter(sum.subarray(0, N + tails[b]), q), fs));
  for (let it = 0; it < 24; it++){
    mix();
    const band = analyze();
    let maxErr = 0;
    for (let b = 0; b < nb; b++){
      const e = targets[b] - band[b];
      maxErr = Math.max(maxErr, Math.abs(e));
      g[b] = Math.max(1e-9, g[b] * Math.pow(10, e / 20));
    }
    if (!(maxErr >= 0.02)) break;
  }
  mix();
  const band = analyze();
  const err = Math.max(...targets.map((t, b) => Math.abs(t - band[b])));
  return { gains: g, bandFmax: band, LAF: D.fastMaxDb(D.sosFilter(sum.subarray(0, N + Math.round(0.5 * fs)), sosA), fs), err: Number.isFinite(err) ? err : Infinity };
}
function buildImpact(){
  const fs = App.audio.fs, bands = activeBands();
  if (bands[bands.length - 1].fc * Math.sqrt(G10) >= fs / 2 * 0.98) throw new Error('この再生サンプリング周波数（' + fs + ' Hz）では選択した帯域の上端を扱えません。');
  const raw = (source === 'wav' && wav) ? wav.bands : synthBands(fs, S.T, bands);
  const fit = fitBands(raw, bands.map(target50), fs, bands);
  const K = calKAll(), N = raw[0].length, dig = new Float32Array(N);
  for (let j = 0; j < raw.length; j++){
    const gd = fit.gains[j] / (P0 * Math.pow(10, K[bands[j].i] / 20)), x = raw[j];
    for (let n = 0; n < N; n++) dig[n] += gd * x[n];
  }
  const runCount = Math.ceil(N / fs / RUN_INTERVAL) + 1;
  const peaks = {
    single: D.peakOf(dig),
    triple: D.overlapPeak(dig, TRIPLE_OFFSETS.map(o => Math.round(o * fs))),
    run: D.overlapPeak(dig, Array.from({ length: runCount }, (_, i) => Math.round(i * RUN_INTERVAL * fs)))
  };
  const LAsum = 10 * Math.log10(bands.reduce((a, b, j) => a + Math.pow(10, (fit.bandFmax[j] + D.aDb(b.fc)) / 10), 0));
  return { bands, bandFmax: fit.bandFmax, LAF: fit.LAF, LAsum, err: fit.err, dig, peaks, dur: N / fs, buffer: null, fs };
}

/* ---------- background noise ---------- */
function bgWindow(f){
  if (S.bgBand === 'b4'){ const lo = BANDS[2].fc, hi = BANDS[5].fc; return D.edgeUp(f, lo / Math.sqrt(G10), lo) * (1 - D.edgeUp(f, hi, hi * Math.sqrt(G10))); }
  return D.edgeUp(f, 18, 28) * (1 - D.edgeUp(f, 9000, 14000));
}
let bgSource = 'synth', bgWav = null, bgNote = '';
/** Background noise at model L_Aeq = 50 dB. Synthetic: random-phase noise with the chosen slope. Recorded: the loaded
    loop with its own spectrum. Either way the A-weighted mean level inside the background band is set to the column value
    by one overall gain; the shape of the spectrum is not altered (apart from the band limits and, when calibrated, K(f)). */
function buildBg(){
  const fs = App.audio.fs, bands = activeBands(), K = calKAll(), useWav = bgSource === 'wav' && !!bgWav;
  const N = useWav ? bgWav.data.length : 1 << Math.ceil(Math.log2(fs * 10)), half = N >> 1;
  const re = new Float64Array(N), im = new Float64Array(N), psd = new Float64Array(half);
  if (useWav){ re.set(bgWav.data); D.fftInPlace(re, im, false); }
  let msA = 0, allA = 0;
  for (let k = 1; k < half; k++){
    const f = k * fs / N, w = bgWindow(f);
    if (useWav){ const a0 = D.aLin(f); allA += (re[k] * re[k] + im[k] * im[k]) * a0 * a0; }
    if (w <= 0) continue;
    psd[k] = (useWav ? re[k] * re[k] + im[k] * im[k] : Math.pow(10, ((S.slope - 3) / 10) * Math.log2(f / 1000))) * w;
    const a = D.aLin(f);
    msA += psd[k] * a * a;
  }
  msA *= 2 / (N * N);
  // Less than 0.01 % of the recording's A-weighted energy lies in the selected band: scaling that residue up to the
  // column level would only amplify the recording's own noise floor.
  if (useWav && !(msA > allA * 2 / (N * N) * 1e-4)){ bgSource = 'synth'; bgNote = '録音のエネルギーが、選んだ暗騒音の帯域にほとんどないため（A 特性で全体の 0.01 % 未満）、合成ノイズに戻しました。'; return buildBg(); }
  const sc = Math.sqrt(P0 * P0 * 1e5 / msA);           // model LAeq = 50 dB
  const band = bands.map(b => {
    const lo = b.fc / Math.sqrt(G10), hi = b.fc * Math.sqrt(G10);
    let m = 0;
    for (let k = Math.ceil(lo * N / fs); k <= Math.floor(hi * N / fs) && k < half; k++) m += psd[k];
    return 10 * Math.log10(Math.max(1e-30, m * 2 / (N * N) * sc * sc) / (P0 * P0));
  });
  const rnd = D.mulberry32(7);
  re[0] = 0; im[0] = 0; re[half] = 0; im[half] = 0;
  for (let k = 1; k < half; k++){
    const ph = 2 * Math.PI * rnd();
    if (psd[k] <= 0){ re[k] = im[k] = re[N - k] = im[N - k] = 0; continue; }
    const kk = sc / (P0 * Math.pow(10, kAtF(K, k * fs / N) / 20));
    if (useWav){ const g = kk * Math.sqrt(bgWindow(k * fs / N)); re[k] *= g; im[k] *= g; }
    else { const amp = kk * Math.sqrt(psd[k]); re[k] = amp * Math.cos(ph); im[k] = amp * Math.sin(ph); }
    re[N - k] = re[k]; im[N - k] = -im[k];
  }
  D.fftInPlace(re, im, true);
  const dig = new Float32Array(N);
  for (let n = 0; n < N; n++) dig[n] = re[n];
  return { dig, peak: D.peakOf(dig), rms: D.rmsOf(dig), band, buffer: null, fs, key: bgKey(), src: useWav ? 'wav' : 'synth' };
}
function bgKey(){ return JSON.stringify([App.audio.fs, S.slope, S.bgBand, calKAll(), activeBands().map(b => b.label), bgSource === 'wav' && bgWav ? bgWav.id : 0]); }
/** Octave-band noise at CAL_DBFS (rms) for calibration, digital domain. */
function buildCalNoise(fs, fc){
  const N = 1 << Math.ceil(Math.log2(fs * 3)), half = N >> 1;
  const re = new Float64Array(N), im = new Float64Array(N);
  const lo = fc / Math.sqrt(G10), hi = fc * Math.sqrt(G10), t = Math.pow(2, 0.08);
  const rnd = D.mulberry32(Math.round(fc));
  let ms = 0;
  for (let k = 1; k < half; k++){
    const f = k * fs / N;
    const w = D.edgeUp(f, lo, lo * t) * (1 - D.edgeUp(f, hi / t, hi));
    if (w <= 0) continue;
    const m = Math.sqrt(w / f), ph = 2 * Math.PI * rnd();
    re[k] = m * Math.cos(ph); im[k] = m * Math.sin(ph); re[N - k] = re[k]; im[N - k] = -im[k];
    ms += 2 * m * m / (N * N);
  }
  D.fftInPlace(re, im, true);
  const sc = Math.pow(10, CAL_DBFS / 20) / Math.sqrt(ms);
  const out = new Float32Array(N);
  for (let n = 0; n < N; n++) out[n] = re[n] * sc;
  return out;
}

/* ---------- level design ----------
   digital = model pressure / (20 uPa * 10^(K/20)); the master gain is the user's output level and nothing else.
   Nothing here depends on which cell or which bands are selected, so comparisons stay valid across cells and band choices. */
const levelLin = (L = S.level) => Math.pow(10, L / 20);
const gradeGain = G => Math.pow(10, (G - 50) / 20);
function cellPeak(G, LA){      // upper bound of the digital peak at output level 0 dB
  if (!impact || !bg) return Infinity;
  return (G === null ? 0 : impact.peaks[S.pattern] * gradeGain(G)) + (LA === null ? 0 : bg.peak * gradeGain(LA));
}
function maxLevelFor(G, LA){ return db(CEIL / cellPeak(G, LA)); }
function playableAt(G, LA, L){ return L <= maxLevelFor(G, LA) + 1e-9; }
function playable(G, LA){ return playableAt(G, LA, S.level); }
/** Highest 5 dB step at which the cell fits, or null when it does not fit even at the minimum level. */
function neededLevel(G, LA){
  const L = Math.min(LEVEL_MAX, Math.floor((maxLevelFor(G, LA) + 1e-9) / LEVEL_STEP) * LEVEL_STEP);
  return L < LEVEL_MIN ? null : L;
}
function cellInfo(G, LA){
  const LAF = G === null ? null : impact.LAF + (G - 50);
  return { LAF, dLA: (LAF !== null && LA !== null) ? LAF - LA : null,
    impPeakDb: G === null ? null : db(impact.peaks[S.pattern] * gradeGain(G)) + S.level,
    bgRmsDb: LA === null ? null : db(bg.rms * gradeGain(LA)) + S.level,
    peakDb: db(cellPeak(G, LA)) + S.level, playable: playable(G, LA), needed: neededLevel(G, LA) };
}

/* ---------- audio engine ---------- */
let graph = null;
const V = { imps: [], bg: null, cal: null };
let playToken = 0, endTimer = 0, loopTimer = 0, loop = null, session = null;

function getGraph(){
  const c = App.audio.ctx;
  if (!c) return null;
  if (!graph || graph.ctx !== c){
    const master = c.createGain();
    master.gain.value = levelLin();
    master.connect(c.destination);
    graph = { ctx: c, master };
    quietAt = 0;
  }
  return graph;
}
/* Every voice gain is a known linear ramp (v0 at t0 -> v1 at t1), so its value at any time is computed here and never
   read back from the AudioParam (no reliance on AudioParam.value or cancelAndHoldAtTime). */
function gainAt(v, t){ return v.t1 > v.t0 ? v.v0 + (v.v1 - v.v0) * Math.min(1, Math.max(0, (t - v.t0) / (v.t1 - v.t0))) : v.v1; }
function rampTo(v, now, to, dur){
  const cur = gainAt(v, now);
  try { const g = v.gain.gain; g.cancelScheduledValues(now); g.setValueAtTime(cur, now); g.linearRampToValueAtTime(to, now + dur); } catch (e) { /* node gone */ }
  v.v0 = cur; v.v1 = to; v.t0 = now; v.t1 = now + dur;
}
let quietAt = 0;                   // context time after which no faded-out or left-over voice can still sound
function holdUntil(t){ if (t > quietAt) quietAt = t; }
function fadeOutVoice(v, now){ rampTo(v, now, 0, FADE); holdUntil(now + FADE + 0.01); try { v.src.stop(now + FADE + 0.01); } catch (e) { /* already stopped */ } }
function stopVoices(now){
  for (const v of V.imps) fadeOutVoice(v, now);
  V.imps = [];
  if (V.bg){ fadeOutVoice(V.bg, now); V.bg = null; }
  if (V.cal){ fadeOutVoice(V.cal, now); V.cal = null; }
}
function cancelTimers(){
  clearTimeout(endTimer); endTimer = 0;
  clearTimeout(loopTimer); loopTimer = 0; loop = null;
}
function stopAll(){
  playToken++;                      // invalidates every pending timer, loop tick and deferred "rebuild then play"
  cancelTimers(); session = null;
  const c = App.audio.ctx;
  if (c) stopVoices(c.currentTime);
  setPlaying(null);
}
function startImpactAt(g, t, gl){
  const src = g.ctx.createBufferSource(), gain = g.ctx.createGain();
  src.buffer = impact.buffer; gain.gain.value = gl;
  src.connect(gain); gain.connect(g.master);
  const v = { src, gain, v0: gl, v1: gl, t0: 0, t1: 0, start: t };
  V.imps.push(v);
  src.onended = () => {
    const i = V.imps.indexOf(v);
    if (i >= 0) V.imps.splice(i, 1);
    try { src.disconnect(); gain.disconnect(); } catch (e) { /* already disconnected */ }
  };
  src.start(t);
}
function startBg(g, LA, start, offset, fadeIn){
  const src = g.ctx.createBufferSource(), gain = g.ctx.createGain(), target = gradeGain(LA);
  src.buffer = bg.buffer; src.loop = true;
  gain.gain.value = 0;
  gain.gain.setValueAtTime(0, start); gain.gain.linearRampToValueAtTime(target, start + fadeIn);
  src.connect(gain); gain.connect(g.master);
  src.start(start, offset);
  V.bg = { src, gain, v0: 0, v1: target, t0: start, t1: start + fadeIn, start, offset, buf: bg.buffer };
}
function runTick(){
  const L = loop, g = graph;
  if (!L || !g || L.token !== playToken) return;
  const now = g.ctx.currentTime;
  // After a throttled/suspended tab, skip missed beats rather than burst-play them.
  if (L.next < now + 0.02) L.next += Math.ceil((now + 0.02 - L.next) / RUN_INTERVAL) * RUN_INTERVAL;
  while (L.next < now + RUN_LOOKAHEAD){ startImpactAt(g, L.next, L.gl); L.next += RUN_INTERVAL; }
  loopTimer = setTimeout(runTick, RUN_TICK_MS);
}
function armEnd(token, tEnd, now){
  clearTimeout(endTimer);
  endTimer = setTimeout(() => {
    if (token !== playToken) return;
    const c = App.audio.ctx;
    if (c && V.bg){ fadeOutVoice(V.bg, c.currentTime); V.bg = null; }
    session = null;
    setPlaying(null);
  }, Math.max(0, (tEnd - now) * 1000));
}
function play(G, LA){
  // The continuous pattern of this very cell is already running: pressing it again changes nothing.
  if (playing && session && S.pattern === 'run' && playing.G === G && playing.LA === LA && loop && loop.token === playToken && !dirty) return true;
  playToken++;
  const token = playToken;
  cancelTimers();
  if (!impact || !bg){ stopAll(); say('lh-over', buildError || '音源を計算できていません。'); return false; }
  if (G !== null && !(impact.err <= 0.25)){
    stopAll(); say('lh-over', '帯域レベルを目標に合わせ込めなかったため再生しません（誤差 ' + fmt(impact.err, 2) + ' dB）。帯域・減衰時間・録音を見直してください。'); return false;
  }
  if (!playable(G, LA)){ stopAll(); showOver(G, LA); return false; }
  if (!App.audio.ensureRunning()){
    stopAll(); say('lh-over', 'この表示画面では音声を使えません。HTML を保存して Safari または Chrome で開いてください。'); return false;
  }
  const g = getGraph(), ctx = g.ctx, now = ctx.currentTime;
  if (!impact.buffer) impact.buffer = App.audio.newBuffer(impact.dig);
  if (!bg.buffer) bg.buffer = App.audio.newBuffer(bg.dig);
  for (const v of V.imps) fadeOutVoice(v, now);      // old impacts are silent by now + FADE; new ones start at now + GAP or later
  V.imps = [];
  if (V.cal){ fadeOutVoice(V.cal, now); V.cal = null; }
  const start = now + GAP;
  let t0, tEnd;
  if (LA === null){
    if (V.bg){ fadeOutVoice(V.bg, now); V.bg = null; }
    t0 = start + 0.04;
  } else if (V.bg && V.bg.buf === bg.buffer){
    rampTo(V.bg, now, gradeGain(LA), FADE);           // same noise keeps running: only its level moves (no gap, same position)
    t0 = start + 0.15;
  } else {
    if (V.bg) fadeOutVoice(V.bg, now);
    startBg(g, LA, start, App.audio.offline ? 0 : Math.random() * bg.buffer.duration * 0.9, BG_FADE_IN);
    t0 = start + LEAD;
  }
  session = { G, LA, times: [] };
  if (G !== null){
    const gl = gradeGain(G);
    if (S.pattern === 'run'){ loop = { token, next: t0, gl }; runTick(); }
    else {
      session.times = (S.pattern === 'triple' ? TRIPLE_OFFSETS : [0]).map(o => t0 + o);
      for (const t of session.times) startImpactAt(g, t, gl);
      tEnd = session.times[session.times.length - 1] + impact.dur + (LA !== null ? 0.4 : 0);
    }
  } else if (S.pattern !== 'run') tEnd = start + 4;
  if (tEnd !== undefined){ session.tEnd = tEnd; armEnd(token, tEnd, now); }
  say('lh-over', '');
  setPlaying({ G, LA });
  return true;
}
/** Settings were edited while this cell was sounding and the sources have been rebuilt: carry the playback over
    without a restart. Strikes that have not started yet are re-issued with the new sound at the same times;
    the background cross-fades to the rebuilt noise at the same loop position. Never starts playback by itself. */
function liveUpdate(cell){
  const g = graph;
  if (!g || !playing || !session) return;
  const { G, LA } = cell;
  if (!impact || !bg || (G !== null && !(impact.err <= 0.25)) || !playable(G, LA)){
    stopAll();
    say('lh-over', '変更後の条件はいまの出力レベルでは再生できない（または計算できない）ため、再生を止めました。');
    return;
  }
  const ctx = g.ctx, now = ctx.currentTime, token = playToken;
  if (!impact.buffer) impact.buffer = App.audio.newBuffer(impact.dig);
  if (!bg.buffer) bg.buffer = App.audio.newBuffer(bg.dig);
  let bgChanged = false;
  if (LA !== null && V.bg && V.bg.buf !== bg.buffer){
    const old = V.bg, t = now + 0.02, d = bg.buffer.duration, pos = ((old.offset + Math.max(0, t - old.start)) % d + d) % d;
    try { const p = old.gain.gain, cur = gainAt(old, t); p.cancelScheduledValues(t); p.setValueAtTime(cur, t); p.linearRampToValueAtTime(0, t + FADE); old.src.stop(t + FADE + 0.01); } catch (e) { /* gone */ }
    holdUntil(t + FADE + 0.01);
    startBg(g, LA, t, pos % (old.buf.duration === d ? d : d), FADE);
    bgChanged = true;
  }
  if (G !== null){
    const gl = gradeGain(G), spacing = S.pattern === 'run' ? RUN_INTERVAL : 1, pending = [];
    for (const v of V.imps.slice()){
      if (v.start > now + 0.03){ pending.push(v.start); try { v.src.onended = null; v.src.stop(0); v.src.disconnect(); v.gain.disconnect(); } catch (e) { /* gone */ } V.imps.splice(V.imps.indexOf(v), 1); }
      else if (bgChanged || impact.dur > spacing) fadeOutVoice(v, now);   // otherwise the sounding strike simply rings out
      else holdUntil(v.start + (v.src.buffer ? v.src.buffer.duration : 0) + 0.01);   // it is the old sound: no level increase until it has ended
    }
    pending.sort((x, y) => x - y);
    if (S.pattern === 'run' && loop && loop.token === token){
      loop.gl = gl; if (pending.length) loop.next = pending[0];
      clearTimeout(loopTimer); runTick();
    } else {
      for (const t of pending) startImpactAt(g, t, gl);
      const last = session.times.length ? session.times[session.times.length - 1] : now;
      session.tEnd = last + impact.dur + (LA !== null ? 0.4 : 0);
      armEnd(token, session.tEnd, now);
    }
  }
  setPlaying({ G, LA });
}
function playCal(i){
  const was = V.cal && V.cal.i === i;
  stopAll();
  if (!was){
    if (!App.audio.ensureRunning()){ say('lh-cal-msg', 'この表示画面では音声を使えません。'); return; }
    const g = getGraph(), ctx = g.ctx, start = ctx.currentTime + GAP, fs = App.audio.fs;
    const src = ctx.createBufferSource(), gain = ctx.createGain();
    src.buffer = App.audio.newBuffer(buildCalNoise(fs, BANDS[i].fc)); src.loop = true;
    gain.gain.value = 0;
    gain.gain.setValueAtTime(0, start); gain.gain.linearRampToValueAtTime(1, start + 0.3);
    src.connect(gain); gain.connect(ctx.destination);       // not through the master: always -20 dBFS
    src.start(start);
    V.cal = { src, gain, v0: 0, v1: 1, t0: start, t1: start + 0.3, start, i };
  }
  updateCalButtons(); renderTransport();
}
function setLevel(L){
  L = Math.min(LEVEL_MAX, Math.max(LEVEL_MIN, Math.round(L / LEVEL_STEP) * LEVEL_STEP));
  if (L === S.level) return;
  let note = '';
  const up = L > S.level;
  if (playing && !playableAt(playing.G, playing.LA, L)){
    stopAll();
    note = '出力レベルを上げると再生中の条件が上限（−1 dBFS）を超えるため、再生を止めました。';
  }
  S.level = L; save();
  if (graph){
    // The shared gain may only rise once every voice that is fading out (stop, cell switch, live edit) is silent;
    // those voices belong to an earlier, possibly louder condition and are not covered by the check above.
    // A decrease is applied at once.
    // Voices created after this point are checked against the new level before they start.
    const p = graph.master.gain, now = graph.ctx.currentTime, t = up ? Math.max(now + 0.01, quietAt + 0.005) : now + 0.01;
    p.cancelScheduledValues(now);     // drops a change that has not begun yet; one already under way simply continues until t
    p.setTargetAtTime(levelLin(), t, 0.02);
  }
  renderMatrix(); renderReadout(); renderTransport(); renderRangeLine();
  say('lh-over', note);
}

/* ---------- rebuild scheduling: one job, always from the latest settings ----------
   live = true: an edit of the sound itself. If something is playing it keeps playing and is carried over to the new
   sound when the rebuild is done; if nothing is playing, nothing starts. */
let rebuildTimer = 0, afterPlay = null, keepLive = null;
function requestRebuild(delay, thenPlay, live){
  if (live && playing && session && !V.cal) keepLive = { G: playing.G, LA: playing.LA, token: playToken };
  else { keepLive = null; stopAll(); }
  dirty = true; App.busy(true);
  afterPlay = thenPlay ? { token: playToken, G: selected.G, LA: selected.LA } : null;
  renderTransport();
  clearTimeout(rebuildTimer);
  rebuildTimer = setTimeout(() => {
    rebuildTimer = 0;
    const k = keepLive, a = afterPlay; keepLive = null; afterPlay = null;
    rebuildNow();
    if (k && k.token === playToken && playing) liveUpdate(k);
    else if (a && a.token === playToken) play(a.G, a.LA);
  }, delay === undefined ? 60 : delay);
}
function rebuildNow(){
  clearTimeout(rebuildTimer); rebuildTimer = 0;
  try { impact = buildImpact(); if (!bg || bg.key !== bgKey()) bg = buildBg(); buildError = null; }
  catch (e) { impact = null; bg = null; buildError = (e && e.message) ? e.message : String(e); }
  dirty = false; App.busy(false);
  renderAll();
}

/* ---------- rendering ---------- */
function say(id, text){ el(id).textContent = text || ''; }
function setPlaying(p){
  playing = p;
  document.querySelectorAll('#lh-matrix .cell.playing').forEach(c => c.classList.remove('playing'));
  if (p){ const c = cellEl(p.G, p.LA); if (c) c.classList.add('playing'); }
  updateCalButtons(); renderTransport();
}
function cellEl(G, LA){ return document.querySelector('#lh-matrix .cell[data-key="' + keyOf(G, LA) + '"]'); }
const patternName = () => ({ single: '単発', triple: '3回', run: '連続' })[S.pattern];
function renderTransport(){
  const st = el('lh-state'), on = !!playing || !!V.cal;
  st.textContent = V.cal ? '校正信号を再生中（−20 dBFS）' : playing ? '再生中（' + (playing.G === null ? (S.pattern === 'run' ? '暗騒音のみ・停止まで' : '暗騒音のみ・4 秒') : patternName()) + '）' : dirty ? '計算中…' : '停止中';
  st.classList.toggle('on', on); el('lh-dot').classList.toggle('on', on);
  const p = playing || selected;
  el('lh-now').textContent = V.cal ? BANDS[V.cal.i].label + ' Hz 帯域ノイズ' : cellName(p.G, p.LA);
  el('lh-level').textContent = sgn(S.level) + ' dB';
  el('lh-level-down').disabled = S.level <= LEVEL_MIN; el('lh-level-up').disabled = S.level >= LEVEL_MAX;
}
function heatPct(info, G, LA){
  if (G === null) return 0;
  if (LA === null) return Math.min(1, Math.max(0, (info.LAF - 25) / 70)) * 35;
  return Math.min(1, Math.max(0, (info.dLA + 10) / 60)) * 55;
}
function renderMatrix(){
  const m = el('lh-matrix'), hadFocus = m.contains(document.activeElement);
  m.textContent = '';
  const add = (tag, cls, html) => { const e = document.createElement(tag); e.className = cls; if (html !== undefined) e.innerHTML = html; m.appendChild(e); return e; };
  add('div', 'mh', ''); add('div', 'mh span', '<span>暗騒音（モデル L<sub>Aeq</sub>、dB(A)）</span>');
  add('div', 'mh', S.mode === 'custom' ? '<span>指定値<br>からの差</span>' : '<span>等級</span>');
  for (const LA of BG_LEVELS) add('div', 'mh col', '<b>' + (LA === null ? 'なし' : LA) + '</b>');
  const ok = !!impact && !!bg;
  for (const G of GRADES.concat([null])){
    const sep = G === null ? ' sep' : '';
    add('div', 'mh' + sep, G === null ? '<b style="font-family:var(--font-ui)">暗騒音<br>のみ</b>' : '<b>' + (S.mode === 'custom' ? (G === 65 ? '指定値' : sgn(G - 65) + ' dB') : 'LH-' + G) + '</b>');
    for (const LA of BG_LEVELS){
      const c = add('button', 'cell' + sep);
      c.type = 'button'; c.tabIndex = -1;
      if (G === null && LA === null){ c.classList.add('blank'); c.disabled = true; c.setAttribute('aria-hidden', 'true'); continue; }
      c.dataset.key = keyOf(G, LA); c.dataset.g = G === null ? '' : String(G); c.dataset.la = LA === null ? '' : String(LA);
      const name = cellName(G, LA);
      if (!ok){ c.classList.add('over'); c.textContent = '—'; c.setAttribute('aria-disabled', 'true'); c.setAttribute('aria-label', name + '：計算できていません'); continue; }
      const info = cellInfo(G, LA);
      if (!info.playable){
        c.classList.add('over'); c.setAttribute('aria-disabled', 'true');
        c.textContent = info.needed === null ? '×' : '≤' + sgn(info.needed).replace('±', '');
        c.setAttribute('aria-label', name + '：いまの出力レベルでは再生できません' + (info.needed === null ? '' : '（出力レベル ' + sgn(info.needed) + ' dB 以下で再生可）'));
      } else {
        const pct = heatPct(info, G, LA);
        if (pct > 0) c.style.background = 'color-mix(in oklab, var(--fill) ' + pct.toFixed(0) + '%, var(--surface))';
        c.textContent = (S.nums && info.dLA !== null) ? 'Δ' + fmt(info.dLA, 0) : '▶';
        c.setAttribute('aria-label', name + ' を再生' + (info.dLA !== null ? '（ΔLA ' + sgn(info.dLA) + ' dB）' : ''));
      }
      if (selected.G === G && selected.LA === LA){ c.classList.add('sel'); c.tabIndex = 0; c.setAttribute('aria-pressed', 'true'); }
      else c.setAttribute('aria-pressed', 'false');
      if (playing && playing.G === G && playing.LA === LA) c.classList.add('playing');
    }
  }
  if (hadFocus){ const c = cellEl(selected.G, selected.LA); if (c) c.focus(); }
}
function bandText(b){ return b.label.replace('k', ' k') + ' Hz'; }
function renderBandLine(){
  const bands = activeBands(), extra = bands.filter(b => b.add);
  const imp = S.mode === 'custom'
    ? '63・125・250・500 Hz（指定値の4帯域のみ）'
    : '63・125・250・500 Hz' + (extra.length ? ' ＋ ' + extra.map(bandText).join('・') + '（参考外挿）' : '（評価範囲のみ）');
  const bgt = (S.bgBand === 'b4' ? '63〜500 Hz の4帯域内（約 45〜710 Hz）' : '広帯域 約 22 Hz〜11 kHz（衝撃音の帯域選択とは独立）') + (bgSource === 'wav' && bgWav ? '／音源：録音「' + bgWav.name + '」' : '／音源：合成ノイズ');
  el('lh-bandline').innerHTML = '<span><b>衝撃音の帯域：</b>' + imp + '</span><span><b>暗騒音の帯域：</b>' + bgt + '</span>';
  el('lh-grade-box').hidden = S.mode === 'custom'; el('lh-custom-box').hidden = S.mode !== 'custom';
}
function renderRangeLine(){
  if (!impact || !bg){ say('lh-rangeline', buildError ? '計算できません：' + buildError : ''); return; }
  const okG = GRADES.filter(G => playable(G, null)), top = GRADES[GRADES.length - 1], all = neededLevel(top, 70);
  let t = 'いまの出力レベル（' + sgn(S.level) + ' dB）で再生できる衝撃音：';
  t += okG.length === GRADES.length ? 'すべて（' + gradeName(GRADES[0]) + '〜' + gradeName(top) + '）。'
    : okG.length ? gradeName(okG[0]) + '〜' + gradeName(okG[okG.length - 1]) + '。' : 'なし。';
  if (okG.length < GRADES.length || !playable(top, 70)){
    t += all === null ? gradeName(top) + ' はこの条件では上限を超えるため再生できません。'
      : gradeName(top) + ' ＋暗騒音 70 まで全部を比べるには出力レベルを ' + sgn(all) + ' dB にします。全条件が同じだけ小さくなり、そのとき暗騒音 20 dB(A) のデジタル出力は約 ' + fmt(db(bg.rms * gradeGain(20)) + all, 0) + ' dBFS です。';
  }
  say('lh-rangeline', t);
}
function renderReadout(){
  const { G, LA } = selected;
  el('lh-ro-title').textContent = cellName(G, LA);
  const nums = el('lh-ro-nums'), warn = el('lh-ro-warn'), model = el('lh-model-nums');
  if (!impact || !bg){
    el('lh-ro-chart').innerHTML = ''; nums.innerHTML = ''; el('lh-ro-table').innerHTML = ''; model.innerHTML = '';
    warn.textContent = buildError ? '計算できません：' + buildError : ''; return;
  }
  const info = cellInfo(G, LA), bands = impact.bands;
  const imp = G === null ? null : impact.bandFmax.map(v => v + G - 50);
  const bgb = LA === null ? null : bg.band.map(v => v + LA - 50);
  drawChart(G, imp, bgb, bands);
  const cal = calibrated(), rows = [];
  rows.push(['再生', info.playable ? '可' : '不可（上限 −1 dBFS 超過）']);
  if (G !== null) rows.push(['衝撃音のデジタル出力（ピーク）', (info.impPeakDb > 0 ? '+' : '') + fmt(info.impPeakDb) + ' dBFS' + (info.impPeakDb > -1 ? '（上限超過）' : '')]);
  if (LA !== null) rows.push(['暗騒音のデジタル出力（実効値）', fmt(info.bgRmsDb) + ' dBFS']);
  rows.push(['出力レベル（全条件共通）', sgn(S.level) + ' dB']);
  rows.push(['音量の基準', cal ? (S.level === 0 ? '校正済み：モデル dB ≒ 受聴位置の dB' : '校正値から ' + sgn(S.level) + ' dB ずらして再生') : '未校正：モデル上の ' + (K_REF - S.level) + ' dB が 0 dBFS']);
  nums.innerHTML = rows.map(r => '<dt>' + r[0] + '</dt><dd>' + r[1] + '</dd>').join('');
  let h = '<thead><tr><th>帯域 Hz</th>' + bands.map(b => '<th>' + b.label + (b.add ? '*' : '') + '</th>').join('') + '</tr></thead><tbody>';
  if (imp) h += '<tr><td>衝撃音 L<sub>Fmax</sub></td>' + imp.map(v => '<td>' + fmt(v, 0) + '</td>').join('') + '</tr>';
  if (bgb) h += '<tr><td>暗騒音 L<sub>eq</sub></td>' + bgb.map(v => '<td>' + fmt(v, 0) + '</td>').join('') + '</tr>';
  if (imp && bgb) h += '<tr><td>差</td>' + imp.map((v, b) => '<td>' + sgn(v - bgb[b]) + '</td>').join('') + '</tr>';
  el('lh-ro-table').innerHTML = h + '</tbody>' + (bands.some(b => b.add) ? '<caption class="hint" style="caption-side:bottom;text-align:left">* 評価範囲外の追加帯域（参考外挿）</caption>' : '');
  let w = '';
  if (!info.playable) w = info.needed === null ? 'この条件はデジタル出力の上限を超えるため再生できません。' : 'この条件を再生するには出力レベルを ' + sgn(info.needed) + ' dB 以下にします（全条件が同じだけ小さくなります）。';
  else if (LA !== null && info.bgRmsDb < QUIET_DBFS) w = '暗騒音のデジタル出力が ' + fmt(info.bgRmsDb, 0) + ' dBFS です。16 bit 出力の量子化雑音（約 −96 dBFS）に近く、ほぼ無音になります。端末の音量を上げずに、出力レベルを上げる（大きい等級は再生不可になります）か、暗騒音の大きい列で比べてください。';
  warn.textContent = w;
  const mrows = [];
  if (G !== null){
    mrows.push(['① 波形から求めた参考 L<sub>A,Fmax</sub>（衝撃音 1 発）', fmt(info.LAF) + ' dB']);
    mrows.push(['② 帯域別 L<sub>Fmax</sub> ＋ A 補正のエネルギー合算（' + impact.bands.map(b => b.label).join('・') + ' Hz）', fmt(impact.LAsum + G - 50) + ' dB']);
    mrows.push(['②−①', sgn(impact.LAsum - impact.LAF, 1) + ' dB']);
  }
  if (LA !== null) mrows.push(['モデル暗騒音 L<sub>Aeq</sub>', fmt(LA) + ' dB']);
  if (info.dLA !== null) mrows.push(['Δ ＝ ① − 暗騒音 L<sub>Aeq</sub>（表のセルの Δ）', sgn(info.dLA, 1) + ' dB']);
  mrows.push(['帯域の合わせ込み誤差（最大）', fmt(impact.err, 2) + ' dB']);
  model.innerHTML = mrows.map(r => '<dt>' + r[0] + '</dt><dd>' + r[1] + '</dd>').join('');
}
function drawChart(G, imp, bgb, bands){
  const custom = S.mode === 'custom', fam = [70, 65, 60, 55, 50, 45, 40, 35, 30];
  if (G !== null && !fam.includes(G)) fam.push(G);
  const refAt = (g, b) => ref50(b) + g - 50, tgt = G === null ? null : bands.map(b => target50(b) + G - 50);
  const all = [].concat(...fam.map(g => bands.map(b => refAt(g, b))), imp || [], bgb || [], tgt || []);
  const lo = Math.max(0, Math.floor((Math.min(...all) - 4) / 10) * 10), hi = Math.ceil((Math.max(...all) + 4) / 10) * 10;
  const x0 = 34, x1 = 312, yT = 8, yB = 206, gw = (x1 - x0) / bands.length;
  const y = v => yB - (Math.max(lo, Math.min(hi, v)) - lo) / (hi - lo) * (yB - yT), cx = j => x0 + gw * (j + 0.5);
  let s = '';
  for (let v = lo; v <= hi; v += (hi - lo > 100 ? 20 : 10)){
    s += '<line class="grid" x1="' + x0 + '" x2="' + x1 + '" y1="' + y(v) + '" y2="' + y(v) + '"/><text x="' + (x0 - 5) + '" y="' + (y(v) + 3.5) + '" text-anchor="end">' + v + '</text>';
  }
  bands.forEach((b, j) => {
    if (bgb) s += '<rect class="bar-bg" x="' + (cx(j) - 5) + '" y="' + y(bgb[j]) + '" width="10" height="' + Math.max(0, yB - y(bgb[j])) + '"><title>暗騒音 ' + b.label + ' Hz L_eq ' + bgb[j].toFixed(1) + ' dB</title></rect>';
    s += '<text x="' + cx(j) + '" y="' + (yB + 14) + '" text-anchor="middle">' + b.label + (b.add ? '*' : '') + '</text>';
  });
  for (const g of fam){
    const on = g === G ? ' on' : '';
    bands.forEach((b, j) => { if (j > 0) s += '<line class="fam' + on + ((b.add || bands[j - 1].add) && !custom ? ' dash' : '') + '" x1="' + cx(j - 1) + '" y1="' + y(refAt(g, bands[j - 1])) + '" x2="' + cx(j) + '" y2="' + y(refAt(g, b)) + '"/>'; });
    s += '<text class="famlab' + on + '" x="' + (cx(bands.length - 1) + 7) + '" y="' + (y(refAt(g, bands[bands.length - 1])) + 3.5) + '">' + (custom ? (g === 65 ? '指定値' : sgn(g - 65)) : 'LH-' + g) + '</text>';
  }
  if (G !== null) bands.forEach((b, j) => {
    const r = refAt(G, b), t = tgt[j];
    if (Math.abs(r - t) > 0.05) s += '<line class="offl" x1="' + cx(j) + '" y1="' + y(r) + '" x2="' + cx(j) + '" y2="' + y(t) + '"><title>' + b.label + ' Hz オフセット ' + (t - r).toFixed(1) + ' dB</title></line>';
    s += '<rect class="tgt" x="' + (cx(j) - 5) + '" y="' + (y(t) - 5) + '" width="10" height="10" transform="rotate(45 ' + cx(j) + ' ' + y(t) + ')"><title>' + b.label + ' Hz 目標値 ' + t.toFixed(1) + ' dB</title></rect>';
    s += '<circle class="ana" cx="' + cx(j) + '" cy="' + y(imp[j]) + '" r="3.3"><title>' + b.label + ' Hz 生成波形の分析値 L_Fmax ' + imp[j].toFixed(1) + ' dB</title></circle>';
  });
  s += '<text x="4" y="' + (yT + 4) + '">dB</text><text x="' + (x1 + 44) + '" y="' + (yB + 14) + '" text-anchor="end">Hz</text>';
  const svg = el('lh-ro-chart');
  svg.innerHTML = s;
  svg.setAttribute('aria-label', (G === null ? '暗騒音のみ' : gradeName(G) + ' の帯域レベル：' + bands.map((b, j) => b.label + ' Hz 目標 ' + tgt[j].toFixed(0) + '、分析値 ' + imp[j].toFixed(1) + ' dB').join('；')) + (bgb ? '。暗騒音：' + bands.map((b, j) => b.label + ' Hz ' + bgb[j].toFixed(0) + ' dB').join('、') : ''));
}
/** Set a field's value unless the user is typing in it (never overwrite text that is being edited). */
function setVal(input, v){ if (document.activeElement !== input){ input.value = v; input.removeAttribute('aria-invalid'); } }
function buildCurveTable(){
  let h = '<thead><tr><th>帯域 Hz</th><th>LH-50 基準値 dB</th><th>オフセット dB</th><th>LH-65 換算</th></tr></thead><tbody>';
  for (const b of BANDS){
    h += '<tr><td>' + b.label + (b.add ? '（追加・参考外挿）' : '') + '</td>' +
      '<td><input type="number" id="lh-cv' + b.i + '" min="0" max="140" step="1" inputmode="decimal" aria-label="' + b.label + ' Hz の LH-50 基準値"></td>' +
      '<td><input type="number" id="lh-of' + b.i + '" min="-60" max="0" step="1" inputmode="decimal" aria-label="' + b.label + ' Hz のオフセット"></td>' +
      '<td id="lh-r65-' + b.i + '"></td></tr>';
  }
  el('lh-tbl-curve').innerHTML = h + '</tbody>';
  const msg = el('lh-curve-msg');
  for (const b of BANDS){
    App.bindNumber(el('lh-cv' + b.i), { min: 0, max: 140, label: b.label + ' Hz の基準値', msg, get: () => S.ref50[b.label],
      set: v => { S.ref50[b.label] = v; save(); el('lh-r65-' + b.i).textContent = fmt(v + 15, 0); requestRebuild(TYPING_MS, false, true); } });
    App.bindNumber(el('lh-of' + b.i), { min: -60, max: 0, label: b.label + ' Hz のオフセット', msg, get: () => S.off[b.label],
      set: v => { S.off[b.label] = v; save(); requestRebuild(TYPING_MS, false, true); } });
  }
}
function updateCurveTable(){
  for (const b of BANDS){
    const cv = el('lh-cv' + b.i), of = el('lh-of' + b.i);
    setVal(cv, S.ref50[b.label]); setVal(of, S.off[b.label]);
    cv.disabled = of.disabled = S.mode === 'custom';
    el('lh-r65-' + b.i).textContent = fmt(S.ref50[b.label] + 15, 0);
  }
  if (S.mode === 'custom') say('lh-curve-msg', '「指定値」の条件では、この表は使いません。');
  else if (el('lh-curve-msg').textContent.indexOf('この表は使いません') >= 0) say('lh-curve-msg', '');
}
function buildCalTable(){
  let h = '<thead><tr><th>帯域 Hz</th><th>信号</th><th>測定 L<sub>eq</sub> dB</th><th>K dB</th></tr></thead><tbody>';
  for (const b of BANDS){
    h += '<tr><td>' + b.label + '</td><td><button type="button" class="btn" data-cal="' + b.i + '" aria-pressed="false">再生</button></td>' +
      '<td><input type="number" step="0.1" min="0" max="140" inputmode="decimal" id="lh-cal' + b.i + '" aria-label="' + b.label + ' Hz の測定値"></td>' +
      '<td id="lh-k' + b.i + '"></td></tr>';
  }
  el('lh-tbl-cal').innerHTML = h + '</tbody>';
}
function updateCalButtons(){
  document.querySelectorAll('#lh-tbl-cal button[data-cal]').forEach(b => {
    const on = !!V.cal && V.cal.i === Number(b.dataset.cal);
    b.textContent = on ? '停止' : '再生'; b.setAttribute('aria-pressed', String(on));
  });
}
function updateCal(){
  const K = calKAll(), any = calibrated();
  for (const b of BANDS){
    const v = S.cal[b.label];
    setVal(el('lh-cal' + b.i), isNum(v) ? v : '');
    el('lh-k' + b.i).textContent = K[b.i].toFixed(1) + (!isNum(v) && any ? '（補間）' : '');
  }
  updateCalButtons();
  const n = BANDS.filter(b => isNum(S.cal[b.label])).length;
  el('lh-chip-cal').textContent = n ? '　校正値あり（' + n + ' 帯域）' : '　未校正';
}
function sourceStatus(){
  const st = el('lh-src-status');
  st.classList.remove('err');
  if (buildError){ st.classList.add('err'); st.textContent = buildError; return; }
  const bad = impact && !(impact.err <= 0.25);
  const e = impact ? '、合わせ込み誤差 ' + fmt(impact.err, 2) + ' dB' + (bad ? '：目標の帯域レベルに合わせられないため、衝撃音は再生しません' : '') : '';
  st.classList.toggle('err', !!bad);
  st.textContent = (source === 'wav' && wav) ? '録音：' + wav.name + '（' + wav.dur.toFixed(2) + ' s を使用' + e + '）' : '合成音（T = ' + S.T + ' s' + e + '）';
}
function syncInputs(){
  document.querySelectorAll('#lh-seg-mode button').forEach(b => b.setAttribute('aria-pressed', String(b.dataset.mode === S.mode)));
  document.querySelectorAll('#lh-add input[data-add]').forEach(c => { c.checked = !!S.add[c.dataset.add]; });
  for (let i = 0; i < 4; i++) setVal(el('lh-c' + i), S.custom4[i]);
  el('lh-bgband').value = S.bgBand;
  document.querySelectorAll('#lh-seg-pattern button').forEach(b => b.setAttribute('aria-pressed', String(b.dataset.p === S.pattern)));
  el('lh-nums').checked = S.nums;
  setVal(el('lh-T'), S.T); el('lh-T').disabled = source === 'wav';
  el('lh-slope').value = String(S.slope);
  el('lh-src-synth').checked = source === 'synth'; el('lh-src-wav').checked = source === 'wav'; el('lh-src-wav').disabled = !wav;
  el('lh-bgsrc-synth').checked = bgSource !== 'wav'; el('lh-bgsrc-wav').checked = bgSource === 'wav'; el('lh-bgsrc-wav').disabled = !bgWav; el('lh-slope').disabled = bgSource === 'wav';
}
function renderAll(){
  syncInputs(); renderBandLine(); renderMatrix(); renderReadout(); renderRangeLine(); updateCurveTable(); updateCal(); sourceStatus(); bgStatus(); renderTransport();
}

/* ---------- actions ---------- */
function showOver(G, LA){
  const box = el('lh-over'), need = neededLevel(G, LA);
  box.textContent = '';
  const t = document.createElement('span');
  t.textContent = cellName(G, LA) + ' は、いまの出力レベル（' + sgn(S.level) + ' dB）ではデジタル出力の上限を超えるため再生しません。';
  box.appendChild(t);
  if (need === null){ t.textContent += 'この条件は最小の出力レベルでも上限を超えます。'; return; }
  const b = document.createElement('button');
  b.type = 'button'; b.className = 'btn'; b.id = 'lh-over-btn';
  b.textContent = '出力レベルを ' + sgn(need) + ' dB に下げて再生（全条件が ' + (S.level - need) + ' dB 小さくなります）';
  b.addEventListener('click', () => { setLevel(need); play(G, LA); });
  box.appendChild(document.createTextNode(' ')); box.appendChild(b);
}
function select(G, LA){
  selected = { G, LA };
  document.querySelectorAll('#lh-matrix .cell').forEach(c => { const on = c.dataset.key === keyOf(G, LA); c.classList.toggle('sel', on); c.tabIndex = on ? 0 : -1; if (c.dataset.key) c.setAttribute('aria-pressed', String(on)); });
  say('lh-over', '');
  renderReadout(); renderTransport();
}
function playSelected(){
  App.audio.ensureRunning();        // resume inside the user gesture (needed on iOS when playback starts after a rebuild)
  if (dirty || !impact || !bg){ requestRebuild(0, true); return; }   // also retries after a failed build
  play(selected.G, selected.LA);
}
function setMode(mode){
  if (mode === S.mode) return;
  intentEpoch++;
  S.mode = mode;
  if (mode === 'custom'){ source = 'synth'; selected = { G: 65, LA: selected.LA }; }
  else if (S.add['16'] && S.T < 0.5) S.T = 0.5;
  if (!GRADES.includes(selected.G) && selected.G !== null) selected.G = 50;
  save(); const lost = reextractWav(); requestRebuild(60); renderBandLine(); syncInputs(); say('lh-over', lost);
}
/** Re-extract the loaded recording for the current bands. Returns a note when the recording had to be dropped. */
function reextractWav(){
  if (wav && wav.originalBuffer && S.mode !== 'custom'){
    try { wav = extractImpact(wav.originalBuffer, wav.name, gradeBands()); }
    catch (e) { wav = null; source = 'synth'; return '録音をこの帯域で使えないため、合成音に戻しました：' + (e && e.message ? e.message : ''); }
  }
  return '';
}
function gradeBands(){ return BANDS.filter(b => !b.add || S.add[b.label]); }
function extractImpact(ab, name, bands){
  const sr = ab.sampleRate, x = ab.getChannelData(0);
  bands = bands || gradeBands();
  if (!D.allFinite(x)) throw new Error('数値として読めないサンプル（NaN／∞）が含まれています。');
  if (sr !== App.audio.fs) throw new Error('サンプリング周波数が再生側と一致しません（' + sr + ' Hz）。');
  if (x.length < sr * 0.3) throw new Error('録音が短すぎます（0.3 s 以上必要）。');
  const fLo = bands[0].fc, fHi = bands[bands.length - 1].fc;
  const pre = D.butterHP(Math.max(8, fLo / Math.sqrt(G10)), sr, 2).concat(D.butterLP(Math.min(sr * 0.45, fHi * Math.sqrt(G10)), sr, 2));
  const y = D.sosFilter(x, pre);
  let best = 0, idx = 0;
  for (let n = 0; n < y.length; n++){ const a = Math.abs(y[n]); if (a > best){ best = a; idx = n; } }
  if (!(best > 0) || !Number.isFinite(best)) throw new Error('信号が入っていません。');
  const i0 = Math.max(0, idx - Math.round(0.08 * sr)), len = Math.min(x.length - i0, Math.round(3 * sr));
  const minLen = Math.max(0.3, 8 / bands[0].fc);            // at least 8 periods of the lowest band
  if (len < Math.round(minLen * sr)) throw new Error(bandText(bands[0]) + ' 帯を合わせ込むには、最大の衝撃の直前から ' + (Math.ceil(minLen * 100) / 100).toFixed(2) + ' s 以上の録音が必要です（この録音で使える長さは ' + (len / sr).toFixed(2) + ' s）。');
  const seg = new Float64Array(len);
  const nIn = Math.round(0.005 * sr), nOut = Math.min(Math.round(0.3 * sr), Math.floor(len / 3));
  for (let n = 0; n < len; n++){
    let e = 1;
    if (n < nIn) e = 0.5 - 0.5 * Math.cos(Math.PI * n / nIn);
    if (n > len - nOut) e *= 0.5 - 0.5 * Math.cos(Math.PI * (len - n) / nOut);
    seg[n] = x[i0 + n] * e;
  }
  const pad = Math.round(0.3 * sr);
  const out = bands.map(b => {
    const z = D.filtfilt(seg, bandSOS(b, sr), b.fc < 50 ? Math.round(sr) : pad);
    if (b.fc < 50){ const edge = Math.min(Math.round(sr * 0.2), Math.floor(z.length / 3)); for (let n = 0; n < edge; n++){ const w = 0.5 - 0.5 * Math.cos(Math.PI * n / edge); z[n] *= w; z[z.length - 1 - n] *= w; } }
    return z;
  });
  return { name, bands: out, dur: len / sr, originalBuffer: ab, channels: ab.numberOfChannels };
}
let wavToken = 0;
let intentEpoch = 0;                // bumped whenever the user changes the impact condition (mode, bands, source, values, new file)
async function loadWav(file){
  const st = el('lh-src-status'), token = ++wavToken, epoch = ++intentEpoch;
  st.classList.remove('err'); st.textContent = '読み込み中…';
  try {
    if (file.size > 200 * 1024 * 1024) throw new Error('200 MB 以下のファイルを選んでください。');
    const ctx = App.audio.create();
    if (!ctx) throw new Error('このブラウザは Web Audio に対応していません。');
    const arr = await file.arrayBuffer();
    const ab = await new Promise((res, rej) => { const p = ctx.decodeAudioData(arr, res, rej); if (p && typeof p.then === 'function') p.then(res, rej); });
    if (token !== wavToken) return;                         // a newer file was chosen meanwhile
    const next = extractImpact(ab, file.name, gradeBands());  // throws: mode, source and saved settings stay untouched
    wav = next;
    if (epoch !== intentEpoch && (S.mode === 'custom' || source !== 'wav')){
      // The user chose another condition while the file was loading: keep that choice, only make the recording available.
      syncInputs();
      st.textContent = '「' + file.name + '」を読み込みました。読み込み中に条件が変更されたため、音源は切り替えていません。「録音（WAV）」を選ぶと使えます。';
      return;
    }
    if (S.mode === 'custom'){ S.mode = 'grade'; save(); }
    source = 'wav';
    requestRebuild(0);
  } catch (err){
    if (token !== wavToken) return;
    st.classList.add('err');
    st.textContent = '読み込めませんでした（音源と条件は変更していません）：' + (err && err.message ? err.message : '音声ファイルとして解釈できません。');
  }
}

let bgWavToken = 0, bgWavId = 0;
/** Load a recording as the background. Nothing is changed unless the whole file was read and prepared successfully. */
async function loadBgWav(file){
  const st = el('lh-bg-status'), token = ++bgWavToken;
  st.classList.remove('err'); st.textContent = '読み込み中…';
  try {
    if (file.size > 200 * 1024 * 1024) throw new Error('200 MB 以下のファイルを選んでください。');
    const ctx = App.audio.create();
    if (!ctx) throw new Error('このブラウザは Web Audio に対応していません。');
    const arr = await file.arrayBuffer();
    const ab = await new Promise((res, rej) => { const p = ctx.decodeAudioData(arr, res, rej); if (p && typeof p.then === 'function') p.then(res, rej); });
    if (token !== bgWavToken) return;                       // a newer file was chosen meanwhile
    const fs = ab.sampleRate, len = Math.min(ab.length, Math.round(fs * 30));
    if (fs !== App.audio.fs) throw new Error('サンプリング周波数が再生側と一致しません（' + fs + ' Hz）。');
    if (len < fs * 2) throw new Error('暗騒音の録音は 2 秒以上必要です（この録音は ' + (ab.length / fs).toFixed(2) + ' s）。');
    const mono = new Float64Array(len);
    for (let c = 0; c < ab.numberOfChannels; c++){
      const d = ab.getChannelData(c).subarray(0, len);
      if (!D.allFinite(d)) throw new Error('数値として読めないサンプル（NaN／∞）が含まれています。');
      for (let i = 0; i < len; i++) mono[i] += d[i] / ab.numberOfChannels;
    }
    if (!(D.peakOf(mono) > 0)) throw new Error('信号が入っていません。');
    let N = 1; while (N * 2 <= len - Math.round(0.05 * fs)) N *= 2;      // loop length: a power of two, leaving a tail to blend
    const xf = Math.min(Math.round(0.5 * fs), len - N), loopData = new Float64Array(N);
    for (let i = 0; i < N; i++) loopData[i] = mono[i];
    for (let i = 0; i < xf; i++){ const u = (i + 0.5) / xf * Math.PI / 2; loopData[i] = mono[i] * Math.sin(u) + mono[N + i] * Math.cos(u); }   // seam: what follows the loop end fades into its start
    bgWav = { id: ++bgWavId, name: file.name, data: loopData, dur: N / fs,
      info: file.name + '：' + (N / fs).toFixed(1) + ' s をループ（全 ' + (ab.length / fs).toFixed(1) + ' s、継ぎ目は ' + (xf / fs).toFixed(2) + ' s のクロスフェード）／' + ab.numberOfChannels + ' ch → モノラル／' + fs + ' Hz／元の A 特性実効値 ' + fmt(aRmsDb(loopData, fs)) + ' dBFS' };
    bgSource = 'wav'; bgNote = '';
    requestRebuild(0, false, true);
  } catch (err){
    if (token !== bgWavToken) return;
    st.classList.add('err');
    st.textContent = '読み込めませんでした（暗騒音は変更していません）：' + (err && err.message ? err.message : '音声ファイルとして解釈できません。');
  }
}
function aRmsDb(x, fs){ const y = D.sosFilter(x, D.aWeightSOS(fs)); return 20 * Math.log10(Math.max(D.rmsOf(y), 1e-30)); }
function bgStatus(){
  const st = el('lh-bg-status');
  if (st.textContent.indexOf('読み込めませんでした') === 0 && !bgNote) return;       // keep the reason of a failed load visible
  st.classList.toggle('err', !!bgNote);
  st.textContent = bgNote || (bgSource === 'wav' && bgWav ? '録音：' + bgWav.info + '。列の dB(A) になるよう全体を一律に増減して再生します。' : '合成ノイズを使用中。');
}

/* ---------- events ---------- */
function bind(){
  el('lh-seg-mode').addEventListener('click', e => { const b = e.target.closest('button[data-mode]'); if (b) setMode(b.dataset.mode); });
  el('lh-add').addEventListener('change', e => {
    const c = e.target.closest('input[data-add]'); if (!c) return;
    S.add[c.dataset.add] = c.checked; intentEpoch++;
    let note = '';
    if (S.add['16'] && S.T < 0.5){ S.T = 0.5; note = '16 Hz を含むため、減衰時間 T を 0.5 s に変更しました。'; }
    save(); const lost = reextractWav(); renderBandLine(); requestRebuild(60, false, true); say('lh-over', [note, lost].filter(Boolean).join(' '));
  });
  for (let i = 0; i < 4; i++){
    App.bindNumber(el('lh-c' + i), { min: 0, max: 140, label: BANDS[i + 2].label + ' Hz の指定値', msg: el('lh-custom-msg'), get: () => S.custom4[i],
      set: v => { S.custom4[i] = v; intentEpoch++; save(); requestRebuild(TYPING_MS, false, true); } });
  }
  el('lh-bgband').addEventListener('change', e => { S.bgBand = e.target.value === 'b4' ? 'b4' : 'wide'; save(); renderBandLine(); requestRebuild(60, false, true); });
  el('lh-seg-pattern').addEventListener('click', e => {
    const b = e.target.closest('button[data-p]'); if (!b || !PATTERNS.includes(b.dataset.p)) return;
    stopAll(); S.pattern = b.dataset.p; save(); syncInputs(); renderMatrix(); renderReadout(); renderRangeLine(); renderTransport();
  });
  el('lh-nums').addEventListener('change', e => { S.nums = e.target.checked; save(); renderMatrix(); });
  const m = el('lh-matrix');
  m.addEventListener('click', e => {
    const c = e.target.closest('.cell'); if (!c || !c.dataset.key) return;
    select(c.dataset.g === '' ? null : Number(c.dataset.g), c.dataset.la === '' ? null : Number(c.dataset.la));
    playSelected();
  });
  m.addEventListener('keydown', e => {
    const d = { ArrowLeft: [0, -1], ArrowRight: [0, 1], ArrowUp: [-1, 0], ArrowDown: [1, 0] }[e.key];
    if (!d) return;
    const rows = GRADES.concat([null]);
    let r = rows.indexOf(selected.G), c = BG_LEVELS.indexOf(selected.LA);
    const nr = Math.min(rows.length - 1, Math.max(0, r + d[0])), nc = Math.min(BG_LEVELS.length - 1, Math.max(0, c + d[1]));
    if (rows[nr] === null && BG_LEVELS[nc] === null) return;
    e.preventDefault();
    select(rows[nr], BG_LEVELS[nc]);
    const t = cellEl(rows[nr], BG_LEVELS[nc]); if (t) t.focus();
  });
  el('lh-play').addEventListener('click', playSelected);
  el('lh-stop').addEventListener('click', stopAll);
  el('lh-level-down').addEventListener('click', () => setLevel(S.level - LEVEL_STEP));
  el('lh-level-up').addEventListener('click', () => setLevel(S.level + LEVEL_STEP));
  App.bindNumber(el('lh-T'), { min: 0.2, max: 2, label: '減衰時間 T', msg: el('lh-src-status'), get: () => S.T,
    set: v => { if (S.mode === 'grade' && S.add['16'] && v < 0.5){ el('lh-T').setAttribute('aria-invalid', 'true'); say('lh-src-status', '16 Hz を含むときは T を 0.5 s 以上にしてください（まだ反映していません）。'); return; } S.T = v; save(); requestRebuild(TYPING_MS, false, true); } });
  el('lh-slope').addEventListener('change', e => { const v = Number(e.target.value); if ([-3, -5, -8].includes(v)){ S.slope = v; save(); requestRebuild(60, false, true); } });
  el('lh-curve-reset').addEventListener('click', () => { const d = defaults(); S.ref50 = d.ref50; S.off = d.off; save(); requestRebuild(60, false, true); });
  const cal = el('lh-tbl-cal');
  cal.addEventListener('click', e => { const b = e.target.closest('button[data-cal]'); if (b) playCal(Number(b.dataset.cal)); });
  cal.addEventListener('input', e => {
    const mm = /^lh-cal(\d)$/.exec(e.target.id || ''); if (!mm) return;
    const b = BANDS[Number(mm[1])], raw = e.target.value.trim(), v = Number(raw);
    if (raw === ''){ delete S.cal[b.label]; e.target.removeAttribute('aria-invalid'); say('lh-cal-msg', ''); }
    else if (!Number.isFinite(v) || v < 0 || v > 140){ e.target.setAttribute('aria-invalid', 'true'); say('lh-cal-msg', b.label + ' Hz：0〜140 の数値を入力してください（まだ反映していません）。'); return; }
    else { S.cal[b.label] = v; e.target.removeAttribute('aria-invalid'); say('lh-cal-msg', ''); }
    save(); requestRebuild(TYPING_MS);
  });
  cal.addEventListener('change', e => {
    const mm = /^lh-cal(\d)$/.exec(e.target.id || ''); if (!mm || e.target.getAttribute('aria-invalid') !== 'true') return;
    const b = BANDS[Number(mm[1])];
    e.target.value = isNum(S.cal[b.label]) ? S.cal[b.label] : ''; e.target.removeAttribute('aria-invalid');
    say('lh-cal-msg', b.label + ' Hz：範囲外だったため、元の値に戻しました。');
  });
  el('lh-cal-clear').addEventListener('click', () => { S.cal = {}; save(); requestRebuild(60); });
  document.querySelectorAll('input[name="lh-src"]').forEach(r => r.addEventListener('change', e => {
    source = (e.target.value === 'wav' && wav) ? 'wav' : 'synth'; intentEpoch++;
    let lost = '';
    if (source === 'wav' && S.mode === 'custom'){ S.mode = 'grade'; save(); lost = reextractWav(); }
    requestRebuild(60); say('lh-over', lost);
  }));
  el('lh-wav').addEventListener('change', e => { const f = e.target.files && e.target.files[0]; if (f) loadWav(f); });
  el('lh-bgwav').addEventListener('change', e => { const f = e.target.files && e.target.files[0]; if (f) loadBgWav(f); });
  document.querySelectorAll('input[name="lh-bgsrc"]').forEach(r => r.addEventListener('change', e => {
    bgSource = (e.target.value === 'wav' && bgWav) ? 'wav' : 'synth'; bgNote = ''; el('lh-bg-status').textContent = '';
    renderBandLine(); requestRebuild(60, false, true);
  }));
}

let started = false;
return {
  init(){ S = loadSettings(); buildCurveTable(); buildCalTable(); bind(); syncInputs(); renderBandLine(); renderMatrix(); updateCurveTable(); updateCal(); renderTransport(); },
  show(){
    if (!started || dirty){ started = true; rebuildNow(); }
  },
  stop(){ stopAll(); },
  /* test / inspection hooks (real objects, no stand-ins) */
  get S(){ return S; }, get impact(){ return impact; }, get bg(){ return bg; }, get playing(){ return playing; }, get dirty(){ return dirty; },
  get selected(){ return selected; }, get bgSource(){ return bgSource; }, get bgWav(){ return bgWav; }, get voices(){ return { imps: V.imps.length, bg: !!V.bg, cal: !!V.cal }; }, get buildError(){ return buildError; }, get graph(){ return getGraph(); },
  BANDS, REF65, GRADES, BG_LEVELS, K_REF, CEIL, RUN_INTERVAL, TRIPLE_OFFSETS, KEY, FADE, GAP, liveUpdate,
  activeBands, calKAll, cellPeak, cellInfo, maxLevelFor, neededLevel, playable, play, playCal, stopAll, setLevel, select, rebuildNow, requestRebuild,
  runTick, sanitize, migrate, defaults, extractImpact, synthBands, fitBands, buildCalNoise
};
})());
