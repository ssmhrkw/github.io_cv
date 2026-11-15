/* TL: airborne sound insulation (Dr / Rr grade curve, band TL values) applied to a source. */
App.register('tl', (function () {
'use strict';
const D = DSP, el = App.el, fmt = App.fmt, sgn = App.sgn, db = App.db;

const KEY = 'tl-listening-v2', OLD_KEY = 'tl-listening-v1';
const FREQ = D.TL_FREQ, FX = D.TL_FX;      // nominal labels / exact octave-step frequencies used in all calculations
const GRADES = [30, 35, 40, 45, 50, 55, 60, 65, 70];
const OFFSETS = [-10, -5, 0, 5, 10];
const RANGES = ['125,2000', '125,4000', '63,8000', '31.5,8000', '16,8000'];
const DEF_LO = 125, DEF_HI = 2000;             // range where the grade curve is defined
const NORM_PEAK = 0.25;                        // -12 dBFS, applied once when a source is loaded
const MAX_SECONDS = 30, CEIL = Math.pow(10, -1 / 20);
const XFADE = 0.08, DRY_RAMP = 1.5, DRY_START_DB = -40, FADE = 0.05;
const VOL_MIN = -60, VOL_MAX = 0, VOL_STEP = 5;
const PADS = [0, -10, -20, -30];
const TYPING_MS = 500;

/* ---------- settings ---------- */
const isNum = v => typeof v === 'number' && Number.isFinite(v);
function defaults(){
  return { v: 2, mode: 'ref', grade: 50, range: '125,2000', values: FX.map(f => Math.max(0, Math.round(D.refCurve(f, 50) * 10) / 10)),
    room: { on: false, S: 10, A: 10 }, pad: 0, vol: -10 };
}
function sanitize(s){
  const o = defaults();
  if (!s || typeof s !== 'object') return o;
  if (s.mode === 'ref' || s.mode === 'custom') o.mode = s.mode;
  if (GRADES.includes(Number(s.grade))) o.grade = Number(s.grade);
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
  return { x: idx.map(i => FX[i]), y: idx.map(i => (S.mode === 'custom' ? S.values[i] : D.refCurve(FX[i], S.grade)) + r), idx };
}
function condName(c){
  if (c === 'dry') return '原音（壁なし）';
  if (S.mode === 'ref') return c === 0 ? '等級 ' + S.grade + ' の基準曲線' : '基準曲線 ' + sgn(c) + ' dB（等級 ' + (S.grade + c) + ' 相当）';
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
function setNoiseSource(){
  source = makeNoise();
  srcInfo = 'ピンクノイズ（8 秒・ループ、ピーク −12 dBFS、' + App.audio.fs + ' Hz）';
  el('tl-src-info').classList.remove('err'); say('tl-src-info', srcInfo);
}
function useNoise(){
  srcToken++;
  stop();
  setNoiseSource();
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
    for (let c = 0; c < ab.numberOfChannels; c++){ const d = ab.getChannelData(c); for (let i = 0; i < n; i++) a[i] += d[i] / ab.numberOfChannels; }
    const gainDb = normalize(a);
    stop();
    source = a;
    srcInfo = f.name + '：' + (n / rate).toFixed(1) + ' 秒を使用' + (ab.length > n ? '（全 ' + (ab.length / rate).toFixed(1) + ' 秒の先頭）' : '') +
      '／' + ab.numberOfChannels + ' ch → モノラル／ピークを −12 dBFS に正規化（' + sgn(gainDb, 1) + ' dB）／' + rate + ' Hz';
    invalidate(false);
    info.classList.remove('err'); info.textContent = srcInfo;
  } catch (e){
    if (token !== srcToken) return;
    info.classList.add('err');
    info.textContent = '読み込めませんでした（音源は変更していません）：' + (e && e.message ? e.message : '音声ファイルとして解釈できません。') + (source ? '　現在の音源：' + srcInfo : '');
  }
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
  return promise;
}
async function doPrepare(myGen){
  const c = curve(), r = rangeOf();
  App.busy(true); renderTransport('処理中…');
  try {
    const res = await runProcess({ input: source, sr: App.audio.fs, x: c.x, y: c.y, offsets: OFFSETS, lo: r.fxLo, hi: r.fxHi });
    if (myGen !== gen) return false;
    if (res.empty) throw new Error('この音源には、選んだ帯域の成分がありません。');
    const wet = {};
    for (const o of OFFSETS) wet[String(o)] = App.audio.newBuffer(res.variants[String(o)]);
    proc = { dry: App.audio.newBuffer(res.dry), wet, stats: res.stats, peak: res.peak, dryPeak: res.dryPeak,
      safe: Math.min(1, CEIL / res.peak), dur: res.dry.length / App.audio.fs, data: { dry: res.dry, variants: res.variants } };
    App.busy(false); renderResult(); renderTransport();
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
function invalidate(resume, delay){
  resumePending = !!resume && (playing || resumePending);
  haltVoices();
  gen++; killWorker(); proc = null; pendingPrep = null; App.busy(false);
  clearTimeout(procTimer);
  renderChart(); renderResult(); renderTransport();
  if (!source) return;
  const req = playReq;
  procTimer = setTimeout(async () => {
    procTimer = 0;
    const ok = await prepare();
    if (ok && resumePending && req === playReq){ resumePending = false; startPlayback(); }
  }, delay === undefined ? 60 : delay);
}

/* ---------- audio engine ---------- */
let graph = null, voices = [], playing = false, cond = 0, epoch = 0, playReq = 0;
function getGraph(){
  const c = App.audio.ctx;
  if (!c) return null;
  if (!graph || graph.ctx !== c){
    const master = c.createGain();
    master.gain.value = Math.pow(10, S.vol / 20);
    master.connect(c.destination);
    graph = { ctx: c, master };
  }
  return graph;
}
function rampValue(r, t){
  if (t < r.t0) return 0;
  if (t >= r.t1) return r.v1;
  const u = (t - r.t0) / (r.t1 - r.t0);
  return r.type === 'exp' ? r.v0 * Math.pow(r.v1 / r.v0, u) : r.v0 + (r.v1 - r.v0) * u;
}
function retire(v, now, end){
  try {
    const g = v.g.gain;
    g.cancelScheduledValues(now);
    g.setValueAtTime(rampValue(v.ramp, now), now);
    g.linearRampToValueAtTime(0, end);
    v.s.stop(end + 0.02);
  } catch (e) { /* already stopped */ }
}
function condTarget(c){ return proc.safe * (c === 'dry' ? Math.pow(10, S.pad / 20) : 1); }
/** Start the current condition at the shared loop position and cross-fade from whatever is sounding. */
function swap(){
  if (!playing || !proc) return;
  const g = getGraph(), ctx = g.ctx, now = ctx.currentTime, t = now + 0.02, dur = proc.dur;
  const pos = (((t - epoch) % dur) + dur) % dur;
  const s = ctx.createBufferSource(), gn = ctx.createGain(), target = condTarget(cond);
  s.buffer = cond === 'dry' ? proc.dry : proc.wet[String(cond)]; s.loop = true;
  s.connect(gn); gn.connect(g.master);
  gn.gain.value = 0; gn.gain.setValueAtTime(0, now);
  let ramp;
  if (cond === 'dry'){              // never jump to the loud original: rise from -40 dB over 1.5 s
    const v0 = target * Math.pow(10, DRY_START_DB / 20);
    gn.gain.setValueAtTime(v0, t); gn.gain.exponentialRampToValueAtTime(target, t + DRY_RAMP);
    ramp = { type: 'exp', t0: t, t1: t + DRY_RAMP, v0, v1: target };
  } else {
    gn.gain.linearRampToValueAtTime(target, t + XFADE);
    ramp = { type: 'lin', t0: now, t1: t + XFADE, v0: 0, v1: target };
  }
  s.start(t, pos);
  const old = voices;
  voices = [{ s, g: gn, ramp }];
  for (const v of old) retire(v, now, t + XFADE);
}
function startPlayback(){
  if (!proc) return false;
  if (!App.audio.ensureRunning()){ say('tl-status', 'この表示画面では音声を使えません。HTML を保存して Safari または Chrome で開いてください。'); return false; }
  const g = getGraph();
  if (!playing){ playing = true; epoch = g.ctx.currentTime + 0.02; }
  swap(); renderResult(); renderTransport();
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
function stop(){ resumePending = false; haltVoices(); }     // user stop: also cancels a pending resume
function haltVoices(){
  playReq++;                                                  // cancels a play that is still waiting for processing
  if (playing || voices.length){
    const c = App.audio.ctx;
    if (c){ const now = c.currentTime; for (const v of voices) retire(v, now, now + FADE); }
  }
  voices = []; playing = false;
  renderTransport();
}
function setCond(c){
  cond = c;
  document.querySelectorAll('#tl-cmp button').forEach(b => b.setAttribute('aria-pressed', String(b.dataset.cond === String(c))));
  renderChart(); renderResult(); renderTransport();
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
    if (S.mode === 'ref') b.innerHTML = (o === 0 ? '基準' : sgn(o) + ' dB') + '<small>等級 ' + (S.grade + o) + (o === 0 ? '' : ' 相当') + '</small>';
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
    h += '<tr class="' + (cond === o ? 'on' : '') + '"><td>' + (o === 0 ? (S.mode === 'ref' ? '基準（等級 ' + S.grade + '）' : '入力値') : sgn(o) + ' dB' + (S.mode === 'ref' ? '（等級 ' + (S.grade + o) + ' 相当）' : '')) + '</td><td>' + fmt(s.dZ) + '</td><td>' + fmt(s.dA) + '</td></tr>';
  }
  tbl.innerHTML = h + '</tbody>';
  let t, note = '';
  if (cond === 'dry') t = '原音です（選んだ帯域に制限した元の音）。';
  else { const s = proc.stats[String(cond)]; t = 'この条件の音は、原音より ' + fmt(s.dZ) + ' dB 小さい（A 特性では ' + fmt(s.dA) + ' dB）。'; }
  if (S.pad) note += '原音だけ ' + (-S.pad) + ' dB 下げて再生する設定です。聞こえる差は、実際より ' + (-S.pad) + ' dB 小さくなります。';
  if (proc.safe < 1) note += '全条件共通のピーク保護 ' + fmt(db(proc.safe)) + ' dB を掛けています（条件間の差は変わりません）。';
  res.textContent = t;
  say('tl-status', note);
}
function renderBands(){
  const box = el('tl-bands'), c = curve(), custom = S.mode === 'custom', msg = el('tl-band-msg');
  box.textContent = '';
  for (const i of c.idx){
    const f = FREQ[i], l = document.createElement('label'), a = document.createElement('input');
    if (!custom && (f < DEF_LO || f > DEF_HI)) l.className = 'outside';
    l.textContent = hz(f) + (!custom && (f < DEF_LO || f > DEF_HI) ? '（外挿）' : '');
    a.type = 'number'; a.min = '0'; a.max = '120'; a.step = '0.1'; a.inputMode = 'decimal'; a.id = 'tl-b' + i;
    a.value = custom ? S.values[i] : Math.max(0, D.refCurve(FX[i], S.grade)).toFixed(1);
    a.disabled = !custom; a.setAttribute('aria-label', hz(f) + ' の値 dB');
    if (custom) App.bindNumber(a, { min: 0, max: 120, label: hz(f), msg, get: () => S.values[i], set: v => { S.values[i] = v; save(); invalidate(true, TYPING_MS); } });
    l.appendChild(a); box.appendChild(l);
  }
  el('tl-grade-row').hidden = custom; el('tl-fill-row').hidden = !custom;
  el('tl-fill').textContent = '等級 ' + S.grade + ' の基準曲線の値を入力欄に写す';
  const ext = c.idx.some(i => FREQ[i] < DEF_LO || FREQ[i] > DEF_HI);
  el('tl-curve-note').textContent = custom
    ? '入力した帯域の値だけを使います（Dr／Rr の等級とは表示しません）。帯域の間は対数周波数軸で補間します。負になる値は 0 dB に制限します。'
    : ext ? '点線の帯域は基準曲線の定義範囲（125 Hz〜2 kHz）の外です。傾きの延長・水平延長による参考外挿で、実際の壁の性能予測ではありません。'
      : '基準曲線の定義範囲（125 Hz〜2 kHz）だけを鳴らします。外挿はしていません。';
}
function renderChart(){
  const c = curve(), off = cond === 'dry' ? 0 : cond, W = 520, x0 = 36, x1 = 508, yT = 14, yB = 186;
  const X = f => x0 + (x1 - x0) * Math.log(f / 11.3) / Math.log(11314 / 11.3), Y = v => yB - Math.max(0, Math.min(100, v)) / 100 * (yB - yT);
  let s = '<rect class="zone" x="' + X(DEF_LO / Math.SQRT2) + '" y="' + yT + '" width="' + (X(DEF_HI * Math.SQRT2) - X(DEF_LO / Math.SQRT2)) + '" height="' + (yB - yT) + '"/>';
  for (let v = 0; v <= 100; v += 20) s += '<line class="grid" x1="' + x0 + '" x2="' + x1 + '" y1="' + Y(v) + '" y2="' + Y(v) + '"/><text x="' + (x0 - 5) + '" y="' + (Y(v) + 3.5) + '" text-anchor="end">' + v + '</text>';
  for (const f of FREQ) s += '<text x="' + X(f) + '" y="' + (yB + 14) + '" text-anchor="middle">' + (f >= 1000 ? f / 1000 + 'k' : f) + '</text>';
  s += '<text x="' + X(500) + '" y="' + (yT + 11) + '" text-anchor="middle">基準曲線の定義範囲</text><text x="' + (W - 12) + '" y="' + (yB + 28) + '" text-anchor="end">Hz</text><text x="4" y="10">dB</text>';
  const custom = S.mode === 'custom', y = c.y.map(v => Math.max(0, v + off));
  c.idx.forEach((bi, j) => {
    const f = FREQ[bi];
    if (j > 0){
      const dash = !custom && (f > DEF_HI || FREQ[c.idx[j - 1]] < DEF_LO) ? ' dash' : '';
      s += '<line class="curve' + dash + '" style="stroke-width:2.5" x1="' + X(FREQ[c.idx[j - 1]]) + '" y1="' + Y(y[j - 1]) + '" x2="' + X(f) + '" y2="' + Y(y[j]) + '"/>';
    }
    s += '<circle class="dot" cx="' + X(f) + '" cy="' + Y(y[j]) + '" r="3"/><text x="' + X(f) + '" y="' + (Y(y[j]) - 7) + '" text-anchor="middle">' + y[j].toFixed(y[j] % 1 ? 1 : 0) + '</text>';
  });
  el('tl-chart').innerHTML = s;
  el('tl-chart').setAttribute('aria-label', '適用する減衰量：' + c.idx.map((bi, j) => hz(FREQ[bi]) + ' ' + y[j].toFixed(1) + ' dB').join('、'));
}
function renderRoom(){
  const r = roomDb();
  el('tl-room-S').disabled = el('tl-room-A').disabled = !S.room.on;
  say('tl-room-out', S.room.on ? '10 log10(A/S) = ' + sgn(r, 1) + ' dB を全帯域の値に加えて、室間音圧レベル差 D の目安として適用します（グラフに反映）。' : '換算しません（A/S = 1：入力値をそのまま減衰量として適用）。');
}
function setVal(input, v){ if (document.activeElement !== input){ input.value = v; input.removeAttribute('aria-invalid'); } }
function syncInputs(){
  document.querySelectorAll('#tl-seg-mode button').forEach(b => b.setAttribute('aria-pressed', String(b.dataset.mode === S.mode)));
  el('tl-grade').value = String(S.grade); el('tl-range').value = S.range; el('tl-pad').value = String(S.pad);
  el('tl-room-on').checked = S.room.on; setVal(el('tl-room-S'), S.room.S); setVal(el('tl-room-A'), S.room.A);
}
function renderAll(){ syncInputs(); renderBands(); renderCmp(); renderRoom(); renderChart(); renderResult(); renderTransport(); }

/* ---------- events ---------- */
function bind(){
  for (const g of GRADES) el('tl-grade').add(new Option(String(g), String(g)));
  el('tl-noise').addEventListener('click', useNoise);
  el('tl-file').addEventListener('change', e => { const f = e.target.files && e.target.files[0]; if (f) loadFile(f); });
  el('tl-seg-mode').addEventListener('click', e => {
    const b = e.target.closest('button[data-mode]'); if (!b || b.dataset.mode === S.mode) return;
    S.mode = b.dataset.mode; save(); syncInputs(); renderBands(); renderCmp(); invalidate(true);
  });
  el('tl-grade').addEventListener('change', e => { const g = Number(e.target.value); if (!GRADES.includes(g)) return; S.grade = g; save(); renderBands(); renderCmp(); invalidate(true); });
  el('tl-fill').addEventListener('click', () => { S.values = FX.map(f => Math.max(0, Math.round(D.refCurve(f, S.grade) * 10) / 10)); save(); renderBands(); invalidate(true); });
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
  show(){ renderAll(); },
  stop,
  get S(){ return S; }, get proc(){ return proc; }, get playing(){ return playing; }, get cond(){ return cond; }, get source(){ return source; },
  get voices(){ return voices.length; }, get srcInfo(){ return srcInfo; }, get graph(){ return getGraph(); },
  FREQ, OFFSETS, KEY, curve, prepare, play, startPlayback, setCond, setVol, useNoise, invalidate, sanitize, defaults
};
})());
