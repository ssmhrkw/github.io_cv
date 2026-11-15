# LH / TL listening 改修 — 進捗と再開手順 (resume file)

作業場所: /home/claude/work
- orig/   : 受領した元HTML（読み取り専用の参照）
- src/    : 新実装のソース（build.mjs が dist/ に単一HTMLを生成）
- tests/  : 再現・回帰試験（Playwright + 実Chromium, NODE_PATH=/opt/npm-tools/node_modules）
- dist/   : 納品物
- shots/  : スクリーンショット

## 再開手順
1. このファイルの「状態」を読む。 2. `ls src tests dist` で実体を確認。 3. 未完了の最初の項目から続行。

## 注意
- _previous_attempt/ は前回の中断した試行の残り（未使用）。src/lh.js, src/tl.js は今回の内容で上書き済み。納品zipに含めない。
- 追加決定: build.mjs は3種を出力（lh-tl-listening.html=タブ付, lh-listening.html, tl-listening.html=単体）。tests/06-variants.mjs で単体版を確認。

## 確定した設計判断（変更時はここを更新）
- 納品: 1ファイルにLH/TLの2タブ（dist/lh-tl-listening.html）。外部フォント読込は削除（外部通信ゼロ）。
- LH 帯域: 基本 63/125/250/500。追加チェック4種 16 / 31.5 / 1k / 2k（独立）。
- LH 基準: C(f)=round(62−A(f)) (LH-65)。16:119, 31.5:101, 63:88,125:78,250:71,500:65,1k:62,2k:61。
  A(f)は厳密中心周波数 1000·10^(0.3x) で評価（IEC 61672-1 表値 16:−56.7, 31.5:−39.4, 2k:+1.2）、整数丸め。表示ラベルは公称。
- LH 音量設計: 「固定スケール + 明示的な出力レベル」。セル切替・帯域切替でスケールを変えない。
  digital = p_model /(P0·10^(K/20)) ; master gain = ユーザーの出力レベルのみ。
  クリップ上限を超えるセルは再生不可表示（必要な出力レベルを提示）。試聴モード(audition)は廃止。
  重なりピークは ×3 の仮定をやめ実波形の重ね合わせで算出。
- LH 3回=1.0 s間隔（元のまま）、連続=1.8 s間隔（仕様）。
- TL: 既定は外挿なし 125 Hz–4 kHz（基準曲線の定義範囲）。拡張（参考外挿）は選択式で残す。
  原音→透過後の実レベル差を保持。原音復帰は緩やかなランプ。原音パッド(明示)あり。D=R+10log10(A/S) 任意換算。
- 保存キー: 新 lh-listening-v5 / tl-listening-v2 / listening-shell-v1。旧キーは読むだけで消さない。

## 実行方法
`cd /home/claude/work && node build.mjs && PLAYWRIGHT_BROWSERS_PATH=/opt/pw-browsers node tests/<name>.mjs`
（node_modules は /opt/npm-tools/node_modules への symlink）

## 状態
- [x] 0. 元HTML読了、参考資料確認
- [x] 1. 元HTMLで16 Hz問題を実ブラウザ再現 (tests/repro-orig.mjs)。結果: K=118固定, 共通減衰 base −14.0 / 31.5Hz込 −21.0 / 16Hz込 −40.4 dB。
      impact peak@LH-80: base −2.1, +31.5Hz +5.3, +16Hz +24.9 dBFS(K=118)。bg50 rms −48.6, peak −36.5 dBFS。原因=LH-80×3重なり仮定+bg70を共通上限に。
      新設計の定数: K_REF=126 dB(未校正, 出力レベル0 dB時に0 dBFS), 上限 −1 dBFS, 出力レベル −40..+20 dB(5 dB刻み)。
- [x] 2. 新実装 src/ 作成（style.css, index.template.html, dsp.js, shell.js, lh.js, tl.js; `node build.mjs` → dist/lh-tl-listening.html）。smoke OK (tests/smoke.mjs)
- [x] 2b. 画面確認とUI調整（360px対応、transport圧縮済）
- [x] 3. 回帰試験 tests/ (bash tests/run-all.sh → 245/245 PASS, ログ tests/out/run.log) : helpers.mjs, 01-unit-static(26/26 OK), 02-lh-audio(87/90; 残3件は試験側の窓の取り方→修正済・要再実行), 03-lh-live(作成済・全PASS想定; WAV fixtureは広帯域バースト), 04-tl(39/39 OK; 途中でTLの再開抜けバグを発見→修正済), 05-ui-persist(未), screenshots(未)
      新設計の実測: base全セル level0で可 / 31.5Hz込も可 / 16Hz込は level0でLH-60まで, 全セルは −20 dB, そのとき bg20 = −106.6 dBFS
- [x] 4. スクリーンショット tests/screenshots.mjs → shots/before-*.png, after-*.png（最終ビルドで再取得済）。最終 245/245 PASS。TLは公称/厳密周波数を統一済(TL_FX)
- [ ] 5. 納品: dist/lh-tl-listening.html, REPORT.md, CHANGELOG.md, zip(src,tests,shots,build.mjs) を SendUserFile
