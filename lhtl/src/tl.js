/* TL: airborne sound insulation (Dr / Rr grade curve, band TL values) applied to a source. */
App.register('tl', (function () {
'use strict';
const D = DSP, el = App.el, fmt = App.fmt, sgn = App.sgn, db = App.db;

const KEY = 'tl-listening-v2', OLD_KEY = 'tl-listening-v1';
const FREQ = D.TL_FREQ, FX = D.TL_FX;      // nominal labels / exact octave-step frequencies used in all calculations
const FAMILIES = { dr: { grades: [30, 35, 40, 45, 50, 55, 60, 65, 70], def: 50, hi: 2000, range: '125,2000' }, t: { grades: [25, 30, 35, 40], def: 30, hi: 4000, range: '125,4000' } };
const BG_LEVELS = [0, 20, 30, 40, 50], SRC_LEVELS = [60, 70, 80, 90];
/** Grade curve of the selected family. dr: Dr / Rr (see DSP.refCurve). t: sash grade lines T-1 ... T-4, stored as their
    500 Hz value 25 / 30 / 35 / 40 dB: flat from 500 Hz to 4 kHz, 5 dB per octave lower below 500 Hz (values not checked against JIS A 4706). */
function refAt(f, g){ return S.family === 't' ? (f >= 500 ? g : g + 5 * Math.log2(f / 500)) : D.refCurve(f, g); }
function defHi(){ return FAMILIES[S.family].hi; }
function gName(g){ return S.family === 't' ? (FAMILIES.t.grades.includes(g) ? 'T-' + (g - 20) / 5 : '等級線 ' + g + ' dB') : '等級 ' + g; }
const OFFSETS = [-10, -5, 0, 5, 10];
const RANGES = ['125,2000', '125,4000', '63,8000', '31.5,8000', '16,8000'];
const DEF_LO = 125;                            // lower end of the range where the grade curves are defined
const NORM_PEAK = 0.25;                        // -12 dBFS, applied once when a source is loaded
const MAX_SECONDS = 30, CEIL = Math.pow(10, -1 / 20);
// Switching: old and new voices cross-fade over the same 10 ms with complementary raised-cosine gains.
// All conditions are zero-phase versions of one signal (fully coherent), so the sum moves monotonically from one
// to the other: no dip, no overshoot. 10 ms keeps the switching splatter > 60 dB below even the attenuated signal
// (tests/study-crossfade.mjs); an abrupt switch is +4 dB above it (audible click).
const XF = 0.010, LEAD = 0.010, START_FADE = 0.2, STOP_FADE = 0.02;
const VOL_MIN = -60, VOL_MAX = 0, VOL_STEP = 5;
const PADS = [0, -10, -20, -30];
const TYPING_MS = 500;

/* ---------- settings ---------- */
const isNum = v => typeof v === 'number' && Number.isFinite(v);
function defaults(){
  return { v: 2, mode: 'ref', grade: 50, range: '125,2000', values: FX.map(f => Math.max(0, Math.round(D.refCurve(f, 50) * 10) / 10)),
    room: { on: false, S: 10, A: 10 }, pad: 0, vol: -10, family: 'dr', bg: 0, srcLevel: 60 };
}
function sanitize(s){
  const o = defaults();
  if (!s || typeof s !== 'object') return o;
  if (s.mode === 'ref' || s.mode === 'custom') o.mode = s.mode;
  if (s.family === 't') o.family = 't';
  if (FAMILIES[o.family].grades.includes(Number(s.grade))) o.grade = Number(s.grade); else o.grade = FAMILIES[o.family].def;
  if (BG_LEVELS.includes(s.bg)) o.bg = s.bg;
  if (SRC_LEVELS.includes(s.srcLevel)) o.srcLevel = s.srcLevel;
  if (RANGES.includes(s.range)) o.range = s.range;
  if (Array.isArray(s.values) && s.values.length === FREQ.length && s.values.every(v => isNum(v) && v >= 0 && v <= 120)) o.values = s.values.slice();
  if (s.room && typeof s.room === 'object'){
    o.room.on = s.room.on === true;
    if (isNum(s.room.S) && s.room.S >= 1 && s.room.S <= 200) o.room.S = s.room.S;
    if (isNum(s.room.A) && s.room.A >= 1 && s.room.A <= 200) o.room.A = s.room.A;
  }
  if (PADS.includes(s.pad)) o.pad = s.pad;
  if (isNum(s.vol)) o.vol = Math.min(VOL_MAX, Math.max(VOL_MIN, Math.round(s.vol / VOL_STEP) * VOL_STEP));
  return o;
}
function loadSettings(){
  const cur = App.store.get(KEY);
  if (cur) return sanitize(cur);
  const old = App.store.get(OLD_KEY);          // v1: { theme, range, mode, grade, values }; read only, never deleted
  return old ? sanitize({ mode: old.mode, grade: old.grade, range: old.range, values: old.values }) : defaults();
}
let S = defaults();
function save(){ App.store.set(KEY, S); }

/* ---------- derived ---------- */
function rangeOf(){ const p = S.range.split(',').map(Number); return { lo: p[0], hi: p[1], fxLo: FX[FREQ.indexOf(p[0])], fxHi: FX[FREQ.indexOf(p[1])] }; }
function bandIdx(){ const r = rangeOf(); return FREQ.map((f, i) => i).filter(i => FREQ[i] >= r.lo && FREQ[i] <= r.hi); }
function roomDb(){ return S.room.on ? 10 * Math.log10(S.room.A / S.room.S) : 0; }
/** Attenuation (dB) at the in-range band centres for offset 0, before the 0 dB floor. */
function curve(){
  const idx = bandIdx(), r = roomDb();
  return { x: idx.map(i => FX[i]), y: idx.map(i => (S.mode === 'custom' ? S.values[i] : refAt(FX[i], S.grade)) + r), idx };
}
function condName(c){
  if (c === 'dry') return '原音（壁なし）';
  if (S.mode === 'ref') return c === 0 ? gName(S.grade) + ' の基準曲線' : '基準曲線 ' + sgn(c) + ' dB（' + gName(S.grade + c) + ' 相当）';
  return c === 0 ? '入力値どおり' : '入力値 ' + sgn(c) + ' dB';
}
const hz = f => (f >= 1000 ? (f / 1000) + ' k' : String(f)) + ' Hz';

/* ---------- source ---------- */
let source = null, srcInfo = '', srcToken = 0;
function normalize(a){
  const pk = D.peakOf(a);
  if (!(pk > 0) || !Number.isFinite(pk)) throw new Error('無音、または数値として読めない音声です。');
  const g = NORM_PEAK / pk;
  for (let i = 0; i < a.length; i++) a[i] *= g;
  return 20 * Math.log10(g);
}
function makeNoise(){
  const fs = App.audio.fs, a = new Float32Array(fs * 8), b = [0, 0, 0, 0, 0, 0, 0];
  let seed = 13421;
  for (let i = 0; i < a.length; i++){
    seed = (1664525 * seed + 1013904223) >>> 0;
    const w = seed / 2147483648 - 1;
    b[0] = .99886 * b[0] + w * .0555179; b[1] = .99332 * b[1] + w * .0750759; b[2] = .969 * b[2] + w * .153852;
    b[3] = .8665 * b[3] + w * .3104856; b[4] = .55 * b[4] + w * .5329522; b[5] = -.7616 * b[5] - w * .016898;
    a[i] = (b[0] + b[1] + b[2] + b[3] + b[4] + b[5] + b[6] + w * .5362) * .11;
    b[6] = w * .115926;
  }
  normalize(a);
  return a;
}
function setNoiseSource(traffic){
  source = makeNoise();
  if (traffic){ const a = source, k = 1 - Math.exp(-2 * Math.PI * 300 / App.audio.fs); let y = 0; for (let i = 0; i < a.length; i++){ y += k * (a[i] - y); a[i] = 0.35 * a[i] + y; } normalize(a); }
  srcInfo = (traffic ? '道路交通騒音の模擬（低音寄りの定常ノイズ。実測スペクトルではありません）' : 'ピンクノイズ') + '（8 秒・ループ、ピーク −12 dBFS、' + App.audio.fs + ' Hz）';
  el('tl-src-info').classList.remove('err'); say('tl-src-info', srcInfo);
}
function useNoise(traffic){
  srcToken++;
  stop();
  setNoiseSource(traffic === true);
  invalidate(false);
}
async function loadFile(f){
  const token = ++srcToken, info = el('tl-src-info');
  info.classList.remove('err'); info.textContent = '読み込み中…';
  try {
    if (f.size > 100 * 1024 * 1024) throw new Error('100 MB 以下の音声を選んでください。');
    const ctx = App.audio.create();
    if (!ctx) throw new Error('このブラウザは Web Audio に対応していません。');
    const arr = await f.arrayBuffer();
    const ab = await new Promise((res, rej) => { const p = ctx.decodeAudioData(arr, res, rej); if (p && typeof p.then === 'function') p.then(res, rej); });
    if (token !== srcToken) return;                          // a newer source was chosen meanwhile
    const rate = ab.sampleRate, n = Math.min(ab.length, Math.round(rate * MAX_SECONDS));
    if (n < rate * 0.2) throw new Error('音声が短すぎます（0.2 秒以上必要）。');
    const a = new Float32Array(n);
    for (let c = 0; c < ab.numberOfChannels; c++){
      const d = ab.getChannelData(c);
      if (!D.allFinite(d.subarray(0, n))) throw new Error('数値として読めないサンプル（NaN／∞）が含まれています。');
      for (let i = 0; i < n; i++) a[i] += d[i] / ab.numberOfChannels;
    }
    const gainDb = normalize(a);
    const prev = { source, srcInfo };
    stop();
    source = a;
    srcInfo = f.name + '：' + (n / rate).toFixed(1) + ' 秒を使用' + (ab.length > n ? '（全 ' + (ab.length / rate).toFixed(1) + ' 秒の先頭）' : '') +
      '／' + ab.numberOfChannels + ' ch → モノラル／ピークを −12 dBFS に正規化（' + sgn(gainDb, 1) + ' dB）／' + rate + ' Hz';
    resetProc();
    info.classList.remove('err'); info.textContent = srcInfo;
    const myGen = gen, ok = await prepare();
    if (!ok && myGen === gen && token === srcToken && prev.source){   // could not be processed: go back to the working source
      const why = el('tl-status').textContent;
      source = prev.source; srcInfo = prev.srcInfo;
      resetProc();
      info.classList.add('err');
      info.textContent = '「' + f.name + '」は処理できなかったため、元の音源に戻しました。' + why + '　現在の音源：' + srcInfo;
      prepare();
    }
  } catch (e){
    if (token !== srcToken) return;
    info.classList.add('err');
    info.textContent = '読み込めませんでした（音源は変更していません）：' + (e && e.message ? e.message : '音声ファイルとして解釈できません。') + (source ? '　現在の音源：' + srcInfo : '');
  }
}

/** Record from the microphone for a fixed time and use it as the source. Everything stays in this browser. */
let recBusy = false, recGen = 0, recStream = null;   // recGen: bumped by stop / tab switch, which abandons a running recording
async function record(seconds){
  const info = el('tl-src-info');
  if (recBusy) return;
  info.classList.remove('err');
  if (!navigator.mediaDevices || !navigator.mediaDevices.getUserMedia || typeof MediaRecorder === 'undefined'){
    info.classList.add('err'); info.textContent = 'この開き方ではマイクを使えません（音源は変更していません）。https で公開したページか、Chrome で開いてください。'; return;
  }
  stop();
  recBusy = true;
  const token = ++srcToken, rg = ++recGen, live = () => token === srcToken && rg === recGen;
  let stream = null;
  try {
    stream = await navigator.mediaDevices.getUserMedia({ audio: { echoCancellation: false, noiseSuppression: false, autoGainControl: false } });
    recStream = stream;
    if (!live()) return;                                      // stopped while the permission prompt was open
    const rec = new MediaRecorder(stream), chunks = [];
    rec.ondataavailable = e => { if (e.data && e.data.size) chunks.push(e.data); };
    const done = new Promise((res, rej) => { rec.onstop = res; rec.onerror = e => rej(e.error || new Error('録音エラー')); });
    rec.start();
    const end = Date.now() + seconds * 1000;
    while (Date.now() < end && live()){ info.textContent = '● 録音中… あと ' + Math.ceil((end - Date.now()) / 1000) + ' 秒（マイクに向かって話してください。■ 停止で中止）'; await new Promise(r => setTimeout(r, 100)); }
    if (rec.state !== 'inactive') rec.stop();
    await done;
    if (!live()){ info.textContent = '録音を中止しました（音源は変更していません）。' + (source ? '現在の音源：' + srcInfo : ''); return; }   // stop, tab switch or another source
    await loadFile(new File([new Blob(chunks, { type: rec.mimeType || 'audio/webm' })], 'マイク録音（' + seconds + ' 秒）'));
  } catch (e){
    info.classList.add('err');
    info.textContent = 'マイクで録音できませんでした（音源は変更していません）：' + (e && e.name === 'NotAllowedError' ? 'マイクの使用が許可されませんでした。' : (e && e.message) || String(e)) + (source ? '　現在の音源：' + srcInfo : '');
  } finally { if (stream) stream.getTracks().forEach(t => t.stop()); recStream = null; recBusy = false; }
}

/* ---------- processing (Worker; falls back to the main thread) ---------- */
let proc = null, gen = 0, worker = null, procTimer = 0, pendingPrep = null;
function killWorker(){ if (worker){ const w = worker; worker = null; try { w.terminate(); } catch (e) { /* gone */ } if (w.cancelJob) w.cancelJob(); } }
function runProcess(job){
  return new Promise((resolve, reject) => {
    let w = null;
    try {
      const code = el('dsp').textContent + '\nonmessage=e=>{try{const d=e.data,r=DSP.tlProcess(d.input,d.sr,d.x,d.y,d.offsets,d.lo,d.hi);postMessage(r);}catch(err){postMessage({error:String(err&&err.message||err)});}};';
      const url = URL.createObjectURL(new Blob([code], { type: 'text/javascript' }));
      w = new Worker(url); URL.revokeObjectURL(url);
    } catch (e) { w = null; }
    if (!w){                                                  // e.g. Workers blocked: compute here (UI pauses briefly)
      try { resolve(D.tlProcess(job.input, job.sr, job.x, job.y, job.offsets, job.lo, job.hi)); } catch (e) { reject(e); }
      return;
    }
    killWorker(); worker = w;
    w.cancelJob = () => reject(new Error('cancelled'));
    w.onmessage = e => { if (worker === w) worker = null; w.terminate(); if (e.data && e.data.error) reject(new Error(e.data.error)); else resolve(e.data); };
    w.onerror = e => { if (worker === w) worker = null; w.terminate(); reject(new Error(e.message || 'Worker error')); };
    w.postMessage(job);
  });
}
/** Build all conditions from the current source and curve. Resolves false when superseded or cancelled.
    Callers arriving while a job for the same generation is running share that job. */
function prepare(){
  if (proc) return Promise.resolve(true);
  if (!source) return Promise.resolve(false);
  if (pendingPrep && pendingPrep.gen === gen) return pendingPrep.promise;
  const promise = doPrepare(gen);
  pendingPrep = { gen, promise };
  const done = () => { if (pendingPrep && pendingPrep.promise === promise) pendingPrep = null; };   // a failed job is retried next time
  promise.then(done, done);
  return promise;
}
async function doPrepare(myGen){
  const c = curve(), r = rangeOf();
  App.busy(true); renderTransport('処理中…');
  try {
    const res = await runProcess({ input: source, sr: App.audio.fs, x: c.x, y: c.y, offsets: OFFSETS, lo: r.fxLo, hi: r.fxHi });
    if (myGen !== gen) return false;
    if (res.empty) throw new Error('この音源には、選んだ帯域の成分がありません。');
    if (!Number.isFinite(res.peak) || !(res.peak > 0) || !OFFSETS.every(o => Number.isFinite(res.stats[String(o)].dZ))) throw new Error('計算結果に数値として扱えない値が含まれています。');
    const wet = {};
    for (const o of OFFSETS) wet[String(o)] = App.audio.newBuffer(res.variants[String(o)]);
    proc = { dry: App.audio.newBuffer(res.dry), wet, stats: res.stats, peak: res.peak, dryPeak: res.dryPeak,
      safe: Math.min(1, CEIL / (res.peak + bgWorstPeak(res.dryMs.a))), dryMsA: res.dryMs.a, bands: res.bands, bandIdx: c.idx, dur: res.dry.length / App.audio.fs, data: { dry: res.dry, variants: res.variants } };
    App.busy(false); renderResult(); renderSpec(); renderTransport();
    return true;
  } catch (e){
    if (myGen !== gen) return false;
    App.busy(false); proc = null;
    say('tl-status', '処理できませんでした：' + (e && e.message ? e.message : String(e)));
    renderTransport();
    return false;
  }
}
/** Curve / range / source changed: drop the processed audio; optionally resume playback with the new condition. */
let resumePending = false;          // playback was interrupted by a condition change and should continue afterwards
/** Drop the processed audio and any running job (the source or the curve changed). */
function resetProc(){
  haltVoices();
  gen++; killWorker(); proc = null; pendingPrep = null; App.busy(false);
  clearTimeout(procTimer); procTimer = 0;
  renderChart(); renderSpec(); renderResult(); renderTransport();
}
function invalidate(resume, delay){
  resumePending = !!resume && (playing || resumePending);
  resetProc();
  if (!source) return;
  const req = playReq;
  procTimer = setTimeout(async () => {
    procTimer = 0;
    const ok = await prepare();
    if (ok && resumePending && req === playReq){ resumePending = false; startPlayback(); }
  }, delay === undefined ? 60 : delay);
}

/* ---------- background noise on the receiving side ----------
   The source is assumed to have the chosen A-weighted level (e.g. conversation 60 dB(A)); the original as played is that
   level. The background is pink noise set to its own A-weighted level relative to it. It is one independent looping
   voice: switching conditions never touches it. Its largest possible peak is reserved in the common gain. */
let bgNoise = null;
function ensureBgNoise(){
  if (bgNoise && bgNoise.fs === App.audio.fs) return bgNoise;
  const fs = App.audio.fs, a = new Float32Array(fs * 8), b = [0, 0, 0]; let seed = 99991;
  for (let i = 0; i < a.length; i++){ seed = (1664525 * seed + 1013904223) >>> 0; const w = seed / 2147483648 - 1; b[0] = .99765 * b[0] + w * .099046; b[1] = .963 * b[1] + w * .2965164; b[2] = .57 * b[2] + w * 1.0526913; a[i] = (b[0] + b[1] + b[2] + w * .1848) * .05; }
  const nf = Math.round(fs * 0.02); for (let i = 0; i < nf; i++){ a[i] *= i / nf; a[a.length - 1 - i] *= i / nf; }
  bgNoise = { fs, data: a, aRms: D.rmsOf(D.sosFilter(a, D.aWeightSOS(fs))), peak: D.peakOf(a), buffer: null };
  return bgNoise;
}
function bgRel(bg, src, dryMsA){ const n = ensureBgNoise(); return Math.sqrt(dryMsA) / n.aRms * Math.pow(10, (bg - src) / 20); }   // gain before the common "safe" factor
function bgWorstPeak(dryMsA){ return ensureBgNoise().peak * bgRel(Math.max(...BG_LEVELS), Math.min(...SRC_LEVELS), dryMsA); }
function bgGain(){ return proc && S.bg ? proc.safe * bgRel(S.bg, S.srcLevel, proc.dryMsA) : 0; }
function updateBg(){
  const g = graph && graph.ctx === App.audio.ctx ? graph : null;
  if (!g) return;
  const now = g.ctx.currentTime, want = playing ? bgGain() : 0;
  if (!g.bg){
    if (!(want > 0)) return;
    // Created once per context and never replaced: toggling can only move this one gain, so voices cannot pile up.
    const n = ensureBgNoise(), s = g.ctx.createBufferSource(), gn = g.ctx.createGain();
    if (!n.buffer) n.buffer = App.audio.newBuffer(n.data);
    s.buffer = n.buffer; s.loop = true; gn.gain.value = 0; s.connect(gn); gn.connect(g.master); s.start(now);
    g.bg = { s, gn };
  }
  const p = g.bg.gn.gain;
  p.cancelScheduledValues(now);       // drop changes that have not begun; the gain then heads for the latest request only
  p.setTargetAtTime(want, now + 0.003, want > 0 ? 0.03 : 0.006);
}

/* ---------- audio engine ---------- */
let graph = null, playing = false, cond = 0, epoch = 0, playReq = 0;
function getGraph(){
  const c = App.audio.ctx;
  if (!c) return null;
  if (!graph || graph.ctx !== c){
    const master = c.createGain();
    master.gain.value = Math.pow(10, S.vol / 20);
    master.connect(c.destination);
    graph = { ctx: c, master, mix: null, live: [], all: new Set(), killEnd: 0 };
  }
  return graph;
}
function fadeCurve(from, to){
  const n = 65, c = new Float32Array(n);
  for (let i = 0; i < n; i++) c[i] = from + (to - from) * (0.5 - 0.5 * Math.cos(Math.PI * i / (n - 1)));
  return c;
}
function fadeParam(param, from, to, t, dur){
  try { param.setValueCurveAtTime(fadeCurve(from, to), t, dur); }
  catch (e) { param.setValueAtTime(from, t); param.linearRampToValueAtTime(to, t + dur); }   // very old engines
}
/** Round a time up to a whole sample frame, so every voice starts sample-aligned with the shared loop position. */
function onFrame(t){ const fs = App.audio.fs; return Math.ceil(t * fs - 1e-6) / fs; }
function condTarget(c){ return proc.safe * (c === 'dry' ? Math.pow(10, S.pad / 20) : 1); }
/* Switching without queues. Every switch builds  mix_n = (1 - w)·mix_(n-1) + w·newVoice  with one raised-cosine w over
   [t, t + XF], t = now + LEAD. Each gain parameter receives exactly one automation and is never interrupted, so a new
   switch never has to wait for an earlier one, however fast they arrive. All conditions are coherent, and the weights
   of every stage add up to one, so the output stays between the signals involved (no dip, no overshoot).
   Stop inserts one more stage that fades the whole current mix, including switches that have not started yet. */
function swap(fade){
  if (!playing || !proc) return;
  const g = getGraph(), ctx = g.ctx, now = ctx.currentTime, t = onFrame(Math.max(now + LEAD, g.killEnd)), dur = proc.dur, T = fade || XF;
  const pos = (((t - epoch) % dur) + dur) % dur;
  const s = ctx.createBufferSource(), v = ctx.createGain(), mix = ctx.createGain();
  s.buffer = cond === 'dry' ? proc.dry : proc.wet[String(cond)]; s.loop = true;
  v.gain.value = 0; fadeParam(v.gain, 0, condTarget(cond), t, T);
  s.connect(v); v.connect(mix); mix.connect(g.master);
  if (g.mix){
    const prev = g.mix, p = ctx.createGain();
    fadeParam(p.gain, 1, 0, t, T);
    prev.connect(p); p.connect(mix);
    try { prev.disconnect(g.master); } catch (e) { /* not connected */ }
    for (const o of g.live) try { o.stop(t + T + 0.01); } catch (e) { /* ended */ }
  }
  s.onended = () => { g.all.delete(s); try { s.disconnect(); v.disconnect(); } catch (e) { /* gone */ } };
  g.mix = mix; g.live = [s]; g.all.add(s);
  s.startAt = t;
  s.start(t, pos);
}
function startPlayback(){
  if (!proc) return false;
  if (!App.audio.ensureRunning()){ say('tl-status', 'この表示画面では音声を使えません。HTML を保存して Safari または Chrome で開いてください。'); return false; }
  const g = getGraph();
  if (!playing){ playing = true; epoch = onFrame(Math.max(g.ctx.currentTime + LEAD, g.killEnd)); swap(cond === 'dry' ? START_FADE : XF); }
  else swap();
  updateBg();
  renderResult(); renderTransport();
  return true;
}
async function play(){
  const req = ++playReq;
  resumePending = false;
  if (playing) return;
  App.audio.ensureRunning();                                  // inside the user gesture
  if (!source) setNoiseSource();
  const ok = await prepare();
  if (!ok || req !== playReq) return;                         // stopped or superseded while processing
  startPlayback();
}
function stop(){               // user stop or tab switch: also cancels a pending resume and a running recording
  resumePending = false;
  if (recBusy){ recGen++; if (recStream) recStream.getTracks().forEach(t => t.stop()); }   // the microphone is released at once
  haltVoices();
}
function haltVoices(){
  playReq++;                                                  // cancels a play that is still waiting for processing
  const g = graph && graph.ctx === App.audio.ctx ? graph : null;
  if (g && g.mix){
    const ctx = g.ctx, t = ctx.currentTime + 0.003, k = ctx.createGain();
    fadeParam(k.gain, 1, 0, t, STOP_FADE);                    // fades everything at once, whatever is still scheduled inside
    g.mix.connect(k); k.connect(g.master);
    try { g.mix.disconnect(g.master); } catch (e) { /* not connected */ }
    // A switch that has not started yet is dropped altogether (stop before start: it never sounds), so pressing
    // "original" and then stop cannot let the loud original through during the fade.
    for (const o of g.all) try { o.stop(o.startAt > ctx.currentTime ? 0 : t + STOP_FADE + 0.01); } catch (e) { /* ended */ }
    g.killEnd = t + STOP_FADE + 0.001;                        // the next start waits these few ms, so old and new never add up
    g.mix = null; g.live = [];
  }
  playing = false;
  updateBg();
  renderTransport();
}
function setCond(c){
  if (c === cond) return;                                     // pressing the selected condition again changes nothing
  cond = c;
  document.querySelectorAll('#tl-cmp button').forEach(b => b.setAttribute('aria-pressed', String(b.dataset.cond === String(c))));
  renderChart(); renderSpec(); renderResult(); renderTransport();
  if (playing) swap();
}
function setVol(v){
  v = Math.min(VOL_MAX, Math.max(VOL_MIN, Math.round(v / VOL_STEP) * VOL_STEP));
  if (v === S.vol) return;
  S.vol = v; save();
  if (graph) graph.master.gain.setTargetAtTime(Math.pow(10, v / 20), graph.ctx.currentTime, 0.03);
  renderTransport();
}

/* ---------- rendering ---------- */
function say(id, text){ el(id).textContent = text || ''; }
function renderTransport(note){
  const st = el('tl-state');
  st.textContent = note || (playing ? '再生中（停止までループ）' : proc ? '停止中（準備完了）' : source ? '停止中' : '音源を用意してください');
  st.classList.toggle('on', playing); el('tl-dot').classList.toggle('on', playing);
  el('tl-now').textContent = condName(cond) + (cond === 'dry' && S.pad ? '　※原音を ' + (-S.pad) + ' dB 下げて再生' : '');
  el('tl-vol').textContent = sgn(S.vol) + ' dB';
  el('tl-vol-down').disabled = S.vol <= VOL_MIN; el('tl-vol-up').disabled = S.vol >= VOL_MAX;
}
function renderCmp(){
  document.querySelectorAll('#tl-cmp button').forEach(b => {
    const c = b.dataset.cond;
    b.setAttribute('aria-pressed', String(c === String(cond)));
    if (c === 'dry') return;
    const o = Number(c);
    if (S.mode === 'ref') b.innerHTML = (o === 0 ? '基準' : sgn(o) + ' dB') + '<small>' + gName(S.grade + o) + (o === 0 ? '' : ' 相当') + '</small>';
    else b.innerHTML = (o === 0 ? '入力値' : sgn(o) + ' dB') + '<small>' + (o === 0 ? 'そのまま' : '平行移動') + '</small>';
  });
}
function renderResult(){
  const res = el('tl-result'), tbl = el('tl-stats');
  if (!proc){ res.textContent = source ? '「▶ 再生」で、この条件を音に適用します。' : '音源を用意して再生してください。'; tbl.innerHTML = ''; return; }
  let h = '<thead><tr><th>条件</th><th>低下量 dB</th><th>A 特性 dB</th></tr></thead><tbody>';
  h += '<tr class="' + (cond === 'dry' ? 'on' : '') + '"><td>原音（壁なし）</td><td>0.0</td><td>0.0</td></tr>';
  for (const o of OFFSETS){
    const s = proc.stats[String(o)];
    h += '<tr class="' + (cond === o ? 'on' : '') + '"><td>' + (o === 0 ? (S.mode === 'ref' ? '基準（' + gName(S.grade) + '）' : '入力値') : sgn(o) + ' dB' + (S.mode === 'ref' ? '（' + gName(S.grade + o) + ' 相当）' : '')) + '</td><td>' + fmt(s.dZ) + '</td><td>' + fmt(s.dA) + '</td></tr>';
  }
  tbl.innerHTML = h + '</tbody>';
  let t, note = '';
  if (cond === 'dry') t = '原音です（選んだ帯域に制限した元の音）。';
  else { const s = proc.stats[String(cond)]; t = 'この条件の音は、原音より ' + fmt(s.dZ) + ' dB 小さい（A 特性では ' + fmt(s.dA) + ' dB）。'; }
  if (S.pad) note += '原音だけ ' + (-S.pad) + ' dB 下げて再生する設定です。聞こえる差は、実際より ' + (-S.pad) + ' dB 小さくなります。';
  if (proc.safe < 1) note += '全条件共通のピーク保護 ' + fmt(db(proc.safe)) + ' dB を掛けています（条件間の差は変わりません）。';
  res.textContent = t;
  say('tl-status', note);
  const lv = S.srcLevel - (cond === 'dry' ? 0 : proc.stats[String(cond)].dA);
  el('tl-est').textContent = '想定：音源が A 特性で ' + S.srcLevel + ' dB のとき、' + (cond === 'dry' ? '壁がなければ約 ' : 'この条件の向こう側では約 ') + fmt(lv, 0) + ' dB(A)。' + (S.bg ? '暗騒音 ' + S.bg + ' dB(A) より ' + fmt(Math.abs(lv - S.bg), 0) + ' dB ' + (lv >= S.bg ? '大きい（聞こえる目安）。' : '小さい（暗騒音に埋もれやすい目安）。') : '') + '音源の大きさは想定値で、実際の音圧ではありません。';
}
function renderBands(){
  const box = el('tl-bands'), c = curve(), custom = S.mode === 'custom', msg = el('tl-band-msg');
  box.textContent = '';
  for (const i of c.idx){
    const f = FREQ[i], l = document.createElement('label'), a = document.createElement('input');
    if (!custom && (f < DEF_LO || f > defHi())) l.className = 'outside';
    l.textContent = hz(f) + (!custom && (f < DEF_LO || f > defHi()) ? '（外挿）' : '');
    a.type = 'number'; a.min = '0'; a.max = '120'; a.step = '0.1'; a.inputMode = 'decimal'; a.id = 'tl-b' + i;
    a.value = custom ? S.values[i] : Math.max(0, refAt(FX[i], S.grade)).toFixed(1);
    a.disabled = !custom; a.setAttribute('aria-label', hz(f) + ' の値 dB');
    if (custom) App.bindNumber(a, { min: 0, max: 120, label: hz(f), msg, get: () => S.values[i], set: v => { S.values[i] = v; save(); invalidate(true, TYPING_MS); } });
    l.appendChild(a); box.appendChild(l);
  }
  el('tl-grade-row').hidden = custom; el('tl-fill-row').hidden = !custom;
  el('tl-fill').textContent = gName(S.grade) + ' の基準曲線の値を入力欄に写す';
  const ext = c.idx.some(i => FREQ[i] < DEF_LO || FREQ[i] > defHi());
  el('tl-curve-note').textContent = custom
    ? '入力した帯域の値だけを使います（Dr／Rr の等級とは表示しません）。帯域の間は対数周波数軸で補間します。負になる値は 0 dB に制限します。'
    : ext ? '点線の帯域は基準曲線の定義範囲（125 Hz〜' + hz(defHi()) + '）の外です。傾きの延長・水平延長による参考外挿で、実際の壁や窓の性能予測ではありません。'
      : '基準曲線の定義範囲（125 Hz〜' + hz(defHi()) + '）だけを鳴らします。外挿はしていません。';
  el('tl-grade-hint').textContent = S.family === 't' ? 'サッシの遮音等級線（T-1：500 Hz 以上 25 dB、T-2：30、T-3：35、T-4：40 dB。500 Hz 未満は 1 オクターブごとに 5 dB 低い）。数値は JIS A 4706 の原典と未照合です。' : '等級 50 の基準点：125 Hz = 35、500 Hz = 50、2 kHz = 60 dB（間は対数周波数軸で補間）';
}
function renderChart(){
  const c = curve(), off = cond === 'dry' ? 0 : cond, OW = App.CHART.OCT, DH = OW * App.CHART.RATIO, x0 = 36, yT = 14;
  const n = c.idx.length, x1 = x0 + n * OW, fx0 = FX[c.idx[0]], fxOf = f => { const i = FREQ.indexOf(f); return i >= 0 ? FX[i] : f; };
  const top = Math.max(60, Math.ceil((Math.max(...c.y.map(v => v + off), refAt(FX[c.idx[n - 1]], Math.max(...FAMILIES[S.family].grades))) + 2) / 10) * 10);
  const yB = yT + top / 10 * DH, W = x1 + 34, H = yB + 32;
  // fixed scale: one octave = OW, 10 dB = OW * 4/3; only the selected bands are drawn
  const X = f => x0 + (Math.log2(fxOf(f) / fx0) + 0.5) * OW, Y = v => yB - Math.max(0, Math.min(top, v)) / 10 * DH;
  const zx0 = Math.max(x0, X(DEF_LO / Math.SQRT2)), zx1 = Math.min(x1, X(defHi() * Math.SQRT2));
  let s = '<rect class="zone" x="' + zx0 + '" y="' + yT + '" width="' + Math.max(0, zx1 - zx0) + '" height="' + (yB - yT) + '"/>';
  for (let v = 0; v <= top; v += 10) s += '<line class="grid" x1="' + x0 + '" x2="' + x1 + '" y1="' + Y(v) + '" y2="' + Y(v) + '"/><text x="' + (x0 - 5) + '" y="' + (Y(v) + 3.5) + '" text-anchor="end">' + v + '</text>';
  for (const bi of c.idx) s += '<text x="' + X(FREQ[bi]) + '" y="' + (yB + 14) + '" text-anchor="middle">' + (FREQ[bi] >= 1000 ? FREQ[bi] / 1000 + 'k' : FREQ[bi]) + '</text>';
  s += '<text x="' + (zx0 + zx1) / 2 + '" y="' + (yT + 11) + '" text-anchor="middle">基準曲線の定義範囲</text><text x="' + (W - 4) + '" y="' + (yB + 28) + '" text-anchor="end">Hz</text><text x="4" y="10">dB</text>';
  const custom = S.mode === 'custom', y = c.y.map(v => Math.max(0, v + off));
  for (const g of FAMILIES[S.family].grades){                                  // the grade reference curves themselves (no offset, no room term)
    const on = !custom && cond !== 'dry' && g === S.grade + off && !S.room.on;
    c.idx.forEach((bi, j) => { if (j > 0) s += '<line class="fam' + (FREQ[bi] > defHi() || FREQ[c.idx[j - 1]] < DEF_LO ? ' dash' : '') + '" x1="' + X(FREQ[c.idx[j - 1]]) + '" y1="' + Y(Math.max(0, refAt(FX[c.idx[j - 1]], g))) + '" x2="' + X(FREQ[bi]) + '" y2="' + Y(Math.max(0, refAt(FX[bi], g))) + '"/>'; });
    if (!on) s += '<text class="famlab" x="' + (X(FREQ[c.idx[c.idx.length - 1]]) + 7) + '" y="' + (Y(Math.max(0, refAt(FX[c.idx[c.idx.length - 1]], g))) + 3.5) + '">' + (S.family === 't' ? 'T-' + (g - 20) / 5 : g) + '</text>';
  }
  c.idx.forEach((bi, j) => {
    const f = FREQ[bi];
    if (j > 0){
      const dash = !custom && (f > defHi() || FREQ[c.idx[j - 1]] < DEF_LO) ? ' dash' : '';
      s += '<line class="curve' + dash + '" style="stroke-width:2.5" x1="' + X(FREQ[c.idx[j - 1]]) + '" y1="' + Y(y[j - 1]) + '" x2="' + X(f) + '" y2="' + Y(y[j]) + '"/>';
    }
    s += '<circle class="dot" cx="' + X(f) + '" cy="' + Y(y[j]) + '" r="3"/><text x="' + X(f) + '" y="' + (Y(y[j]) - 7) + '" text-anchor="middle">' + y[j].toFixed(y[j] % 1 ? 1 : 0) + '</text>';
  });
  el('tl-chart').innerHTML = s;
  el('tl-chart').setAttribute('viewBox', '0 0 ' + W + ' ' + H); el('tl-chart').style.maxWidth = Math.round(W * 1.3) + 'px';
  el('tl-chart').setAttribute('aria-label', '適用する減衰量：' + c.idx.map((bi, j) => hz(FREQ[bi]) + ' ' + y[j].toFixed(1) + ' dB').join('、'));
}
/** Octave-band levels of the actual sound: source vs. what remains behind the selected condition (relative dB). */
function renderSpec(){
  const svg = el('tl-spec'), note = el('tl-spec-note');
  if (!proc || !proc.bands){ svg.innerHTML = ''; svg.setAttribute('viewBox', '0 0 10 1'); svg.setAttribute('aria-label', '音のスペクトル（未計算）'); note.textContent = '「▶ 再生」で音源を処理すると、原音と遮音後の帯域レベルをここに表示します。'; return; }
  const idx = proc.bandIdx, dry = proc.bands.dry, wet = cond === 'dry' ? null : proc.bands.variants[String(cond)], top = Math.max(...dry);
  const OW = App.CHART.OCT, DH = OW * App.CHART.RATIO, x0 = 36, yT = 14, gw = OW, x1 = x0 + idx.length * OW;
  const lo = Math.max(-120, Math.floor((Math.min(...(wet || dry)) - top - 2) / 10) * 10), yB = yT + -lo / 10 * DH, Y = v => yB - (Math.max(lo, Math.min(0, v)) - lo) / 10 * DH;
  let s = '';
  for (let v = 0; v >= lo; v -= 10) s += '<line class="grid" x1="' + x0 + '" x2="' + x1 + '" y1="' + Y(v) + '" y2="' + Y(v) + '"/><text x="' + (x0 - 5) + '" y="' + (Y(v) + 3.5) + '" text-anchor="end">' + v + '</text>';
  const parts = [];
  idx.forEach((bi, j) => {
    const cx = x0 + gw * (j + 0.5), d = dry[j] - top, w = wet ? wet[j] - top : null, bw = Math.min(14, gw / 2 - 3);
    s += '<rect class="bar-bg" x="' + (cx - bw - 1) + '" y="' + Y(d) + '" width="' + bw + '" height="' + (yB - Y(d)) + '"><title>原音 ' + hz(FREQ[bi]) + ' ' + d.toFixed(1) + ' dB</title></rect>';
    if (wet) s += '<rect class="bar-imp" x="' + (cx + 1) + '" y="' + Y(w) + '" width="' + bw + '" height="' + (yB - Y(w)) + '"><title>遮音後 ' + hz(FREQ[bi]) + ' ' + w.toFixed(1) + ' dB</title></rect><text x="' + cx + '" y="' + (yB + 26) + '" text-anchor="middle">−' + (d - w).toFixed(0) + '</text>';
    s += '<text x="' + cx + '" y="' + (yB + 13) + '" text-anchor="middle">' + (FREQ[bi] >= 1000 ? FREQ[bi] / 1000 + 'k' : FREQ[bi]) + '</text>';
    parts.push(hz(FREQ[bi]) + ' 原音 ' + d.toFixed(0) + (wet ? '、遮音後 ' + w.toFixed(0) : '') + ' dB');
  });
  s += '<text x="4" y="10">dB</text><text x="' + (x1 + 6) + '" y="' + (yB + 13) + '" text-anchor="end">Hz</text>' + (wet ? '<text x="' + (x0 - 5) + '" y="' + (yB + 26) + '" text-anchor="end">差</text>' : '');
  svg.setAttribute('viewBox', '0 0 ' + (x1 + 34) + ' ' + (yB + 34)); svg.style.maxWidth = Math.round((x1 + 34) * 1.3) + 'px';
  svg.innerHTML = s; svg.setAttribute('aria-label', '音の帯域レベル（原音の最大帯域を 0 dB とした相対値）：' + parts.join('；'));
  note.textContent = 'いまの音源の 1/1 オクターブ帯域レベルです（原音の最大帯域を 0 dB とした相対値。実際の音圧ではありません）。灰色＝原音、濃い色＝' + (wet ? condName(cond) + ' の遮音後。下段の数字は帯域ごとの差で、上のグラフの減衰量に対応します。' : '（原音を選択中のため遮音後は表示していません）。');
}
function renderRoom(){
  const r = roomDb();
  el('tl-room-S').disabled = el('tl-room-A').disabled = !S.room.on;
  say('tl-room-out', S.room.on ? '10 log10(A/S) = ' + sgn(r, 1) + ' dB を全帯域の値に加えて、室間音圧レベル差 D の目安として適用します（グラフに反映）。' : '換算しません（A/S = 1：入力値をそのまま減衰量として適用）。');
}
function setVal(input, v){ if (document.activeElement !== input){ input.value = v; input.removeAttribute('aria-invalid'); } }
function syncInputs(){
  document.querySelectorAll('#tl-seg-mode button').forEach(b => b.setAttribute('aria-pressed', String(b.dataset.mode === S.mode)));
  el('tl-family').value = S.family; el('tl-grade').value = String(S.grade); el('tl-range').value = S.range; el('tl-pad').value = String(S.pad);
  el('tl-bg').value = String(S.bg); el('tl-srclevel').value = String(S.srcLevel);
  el('tl-room-on').checked = S.room.on; setVal(el('tl-room-S'), S.room.S); setVal(el('tl-room-A'), S.room.A);
}
function renderAll(){ syncInputs(); renderBands(); renderCmp(); renderRoom(); renderChart(); renderSpec(); renderResult(); renderTransport(); }

/* ---------- events ---------- */
function fillGrades(){ const sel = el('tl-grade'); sel.textContent = ''; for (const g of FAMILIES[S.family].grades) sel.add(new Option(S.family === 't' ? 'T-' + (g - 20) / 5 : String(g), String(g))); sel.value = String(S.grade); }
function bind(){
  fillGrades();
  el('tl-family').addEventListener('change', e => { S.family = e.target.value === 't' ? 't' : 'dr'; S.grade = FAMILIES[S.family].def; S.range = FAMILIES[S.family].range; save(); fillGrades(); syncInputs(); renderBands(); renderCmp(); invalidate(true); });
  el('tl-traffic').addEventListener('click', () => useNoise(true));
  el('tl-rec5').addEventListener('click', () => record(5)); el('tl-rec10').addEventListener('click', () => record(10));
  el('tl-bg').addEventListener('change', e => { const v = Number(e.target.value); if (BG_LEVELS.includes(v)){ S.bg = v; save(); updateBg(); renderResult(); } });
  el('tl-srclevel').addEventListener('change', e => { const v = Number(e.target.value); if (SRC_LEVELS.includes(v)){ S.srcLevel = v; save(); updateBg(); renderResult(); } });
  el('tl-noise').addEventListener('click', useNoise);
  el('tl-file').addEventListener('change', e => { const f = e.target.files && e.target.files[0]; if (f) loadFile(f); });
  el('tl-seg-mode').addEventListener('click', e => {
    const b = e.target.closest('button[data-mode]'); if (!b || b.dataset.mode === S.mode) return;
    S.mode = b.dataset.mode; save(); syncInputs(); renderBands(); renderCmp(); invalidate(true);
  });
  el('tl-grade').addEventListener('change', e => { const g = Number(e.target.value); if (!FAMILIES[S.family].grades.includes(g)) return; S.grade = g; save(); renderBands(); renderCmp(); invalidate(true); });
  el('tl-fill').addEventListener('click', () => { S.values = FX.map(f => Math.max(0, Math.round(refAt(f, S.grade) * 10) / 10)); save(); renderBands(); invalidate(true); });
  el('tl-range').addEventListener('change', e => { if (!RANGES.includes(e.target.value)) return; S.range = e.target.value; save(); renderBands(); invalidate(true); });
  el('tl-room-on').addEventListener('change', e => { S.room.on = e.target.checked; save(); renderRoom(); invalidate(true); });
  App.bindNumber(el('tl-room-S'), { min: 1, max: 200, label: '壁の面積 S', msg: el('tl-room-out'), get: () => S.room.S, set: v => { S.room.S = v; save(); renderRoom(); invalidate(true, TYPING_MS); } });
  App.bindNumber(el('tl-room-A'), { min: 1, max: 200, label: '等価吸音面積 A', msg: el('tl-room-out'), get: () => S.room.A, set: v => { S.room.A = v; save(); renderRoom(); invalidate(true, TYPING_MS); } });
  el('tl-cmp').addEventListener('click', e => { const b = e.target.closest('button[data-cond]'); if (!b) return; setCond(b.dataset.cond === 'dry' ? 'dry' : Number(b.dataset.cond)); });
  el('tl-pad').addEventListener('change', e => { const v = Number(e.target.value); if (!PADS.includes(v)) return; S.pad = v; save(); renderResult(); renderTransport(); if (playing && cond === 'dry') swap(); });
  el('tl-play').addEventListener('click', play);
  el('tl-stop').addEventListener('click', stop);
  el('tl-vol-down').addEventListener('click', () => setVol(S.vol - VOL_STEP));
  el('tl-vol-up').addEventListener('click', () => setVol(S.vol + VOL_STEP));
}

return {
  init(){ S = loadSettings(); bind(); renderAll(); },
  record, get bgInfo(){ return bgNoise ? { gain: bgGain(), aRms: bgNoise.aRms, peak: bgNoise.peak, data: bgNoise.data } : null; },
  show(){ renderAll(); },
  stop,
  get S(){ return S; }, get proc(){ return proc; }, get playing(){ return playing; }, get cond(){ return cond; }, get source(){ return source; },
  get voices(){ return graph && graph.mix ? graph.live.length : 0; }, get sounding(){ return graph ? graph.all.size : 0; }, XF, LEAD, START_FADE, get srcInfo(){ return srcInfo; }, get graph(){ return getGraph(); },
  FREQ, OFFSETS, KEY, curve, prepare, play, startPlayback, setCond, setVol, useNoise, invalidate, sanitize, defaults
};
})());
