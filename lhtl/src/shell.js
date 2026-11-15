/* Shell: shared audio context, storage, theme, tabs. */
const App = (function () {
'use strict';
const VERSION = '__VERSION__';
const SHELL_KEY = 'listening-shell-v1';
const modules = [];
const ALL_TABS = ['lh', 'tl'];
/** Tabs that exist in this build (the single-purpose builds contain one panel and no tab bar). */
function tabsPresent(){ return ALL_TABS.filter(t => document.getElementById('panel-' + t)); }

/** Strict element lookup: a missing element is a bug and must be loud, never silently replaced. */
function el(id){
  const e = document.getElementById(id);
  if (!e) throw new Error('画面要素 #' + id + ' が見つかりません');
  return e;
}
const store = {
  get(key){ try { const v = localStorage.getItem(key); return v ? JSON.parse(v) : null; } catch (e) { return null; } },
  set(key, value){ try { localStorage.setItem(key, JSON.stringify(value)); return true; } catch (e) { return false; } }
};
const fmt = (v, d = 1) => (v === null || v === undefined || !Number.isFinite(v)) ? '—' : v.toFixed(d).replace('-', '−');
const sgn = (v, d = 0) => (v > 0 ? '+' : v < 0 ? '−' : '±') + Math.abs(v).toFixed(d);
const db = v => 20 * Math.log10(Math.max(v, 1e-30));
/** Text that did not originate in this page (file names, stored strings) must pass through here before it is put into HTML. */
const esc = v => String(v).replace(/[&<>"']/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' })[c]);

/** Numeric field: every valid keystroke is stored immediately; invalid text is never stored and is reverted on blur. */
function bindNumber(input, o){
  const say = t => { if (o.msg) o.msg.textContent = t; };
  input.addEventListener('input', () => {
    const raw = input.value.trim(), v = Number(raw);
    if (raw === '' || !Number.isFinite(v) || v < o.min || v > o.max){
      input.setAttribute('aria-invalid', 'true');
      say(o.label + '：' + o.min + '〜' + o.max + ' の数値を入力してください（まだ反映していません）。');
      return;
    }
    input.removeAttribute('aria-invalid'); say('');
    if (v !== o.get()){ o.set(v); }
  });
  input.addEventListener('change', () => {
    if (input.getAttribute('aria-invalid') === 'true'){
      input.value = String(o.get());
      input.removeAttribute('aria-invalid');
      say(o.label + '：範囲外だったため、元の値 ' + o.get() + ' に戻しました。');
    }
  });
}

const audio = {
  ctx: null, fs: 48000, offline: false,
  create(){
    if (this.ctx) return this.ctx;
    const AC = window.AudioContext || window.webkitAudioContext;
    if (!AC) return null;
    try { this.ctx = new AC(); this.fs = this.ctx.sampleRate; } catch (e) { this.ctx = null; }
    return this.ctx;
  },
  ensureRunning(){
    const c = this.create();
    if (!c) return false;
    if (this.offline) return true;
    try { if (navigator.audioSession) navigator.audioSession.type = 'playback'; } catch (e) { /* not supported */ }
    if (c.state !== 'running'){ const p = c.resume(); if (p && p.catch) p.catch(() => {}); }
    return true;
  },
  newBuffer(data){
    let b;
    try { b = new AudioBuffer({ length: data.length, numberOfChannels: 1, sampleRate: this.fs }); }
    catch (e) { b = this.ctx.createBuffer(1, data.length, this.fs); }
    if (b.copyToChannel) b.copyToChannel(data, 0); else b.getChannelData(0).set(data);
    return b;
  },
  /** Test hook: render through another real context (e.g. OfflineAudioContext) with the unchanged engine code. */
  useContext(c){
    stopAll();
    this.ctx = c; this.fs = c.sampleRate;
    this.offline = (typeof OfflineAudioContext !== 'undefined' && c instanceof OfflineAudioContext);
  }
};

let shell = { tab: 'lh', theme: 'light' };
function loadShell(){
  const s = store.get(SHELL_KEY) || {};
  const legacyLh = store.get('lh-listening-v4-range') || store.get('lh-listening-v3-2khz') || store.get('lh-listening-v2') || store.get('lh-listening-v1') || {};
  const legacyTl = store.get('tl-listening-v1') || {};
  const themes = ['light', 'dark', 'black'];
  let theme = [s.theme, legacyLh.theme, legacyTl.theme].find(t => themes.includes(t));
  if (!theme) theme = (window.matchMedia && window.matchMedia('(prefers-color-scheme: dark)').matches) ? 'dark' : 'light';
  shell = { tab: s.tab === 'tl' ? 'tl' : 'lh', theme };   // the saved tab is only used by the two-tab build
}
function applyTheme(){
  document.documentElement.dataset.theme = shell.theme;
  document.querySelectorAll('#seg-theme button').forEach(b => b.setAttribute('aria-pressed', String(b.dataset.theme === shell.theme)));
}
function stopAll(){ for (const m of modules) { try { m.mod.stop(); } catch (e) { /* keep stopping the others */ } } }
function showTab(name, focus){
  const tabs = tabsPresent();
  if (!tabs.includes(name)) name = tabs[0];
  stopAll();
  if (tabs.length > 1){ shell.tab = name; store.set(SHELL_KEY, shell); }
  for (const t of tabs){
    const on = t === name, tab = document.getElementById('tab-' + t);
    if (tab){ tab.setAttribute('aria-selected', String(on)); tab.tabIndex = on ? 0 : -1; if (on && focus) tab.focus(); }
    el('panel-' + t).hidden = !on;
  }
  const m = modules.find(x => x.name === name);
  if (m && m.ok && m.mod.show) { try { m.mod.show(); } catch (e) { showError(e); } }
}
function showError(e){
  const box = document.getElementById('app-error');
  const text = '起動または処理でエラーが発生しました：' + (e && e.message ? e.message : String(e));
  if (box){ box.hidden = false; box.textContent = text; }
  if (window.console) console.error(e);
}
function busy(on){ const b = document.getElementById('busy'); if (b) b.hidden = !on; }
function register(name, mod){ modules.push({ name, mod, ok: false }); api[name] = mod; }
function start(){
  loadShell(); applyTheme();
  audio.create();
  el('seg-theme').addEventListener('click', e => {
    const b = e.target.closest('button[data-theme]'); if (!b) return;
    shell.theme = b.dataset.theme; applyTheme(); store.set(SHELL_KEY, shell);
    for (const m of modules) if (m.ok && m.mod.themeChanged) m.mod.themeChanged();
  });
  const tabs = tabsPresent().filter(t => document.getElementById('tab-' + t));
  for (const t of tabs){
    el('tab-' + t).addEventListener('click', () => showTab(t, false));
    el('tab-' + t).addEventListener('keydown', e => {
      const i = tabs.indexOf(t);
      let j = -1;
      if (e.key === 'ArrowRight') j = (i + 1) % tabs.length; else if (e.key === 'ArrowLeft') j = (i + tabs.length - 1) % tabs.length;
      else if (e.key === 'Home') j = 0; else if (e.key === 'End') j = tabs.length - 1;
      if (j >= 0){ e.preventDefault(); showTab(tabs[j], true); }
    });
  }
  document.addEventListener('keydown', e => { if (e.key === 'Escape') stopAll(); });
  window.addEventListener('pagehide', stopAll);
  for (const m of modules){
    try { m.mod.init(); m.ok = true; } catch (e) { showError(e); }
  }
  const hash = (location.hash || '').replace('#', '');
  showTab(hash === 'tl' || hash === 'lh' ? hash : shell.tab, false);
}
/** Chart scale used for every level-versus-frequency chart: one octave = 15 mm and 10 dB = 20 mm on paper
    (the format used for building-acoustics test reports), i.e. 10 dB is 4/3 of an octave. OCT is the octave width in SVG units. */
const CHART = { OCT: 36, RATIO: 20 / 15 };
const api = { VERSION, el, store, fmt, sgn, db, esc, CHART, bindNumber, audio, register, start, busy, showError, stopAll, showTab, get shell(){ return shell; } };
return api;
})();
window.__app = App;
