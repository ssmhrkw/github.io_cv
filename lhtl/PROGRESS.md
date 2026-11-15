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
- [x] 5. 納品: /home/claude/deliver/ に3HTML, REPORT.md, CHANGELOG.md, compare-*.png, project.zip。最終 252/252 PASS (tests/out/run.log)。SendUserFile 済みなら完了。

## 仕上げ作業（2回目の指示: 最後まで実装してZIP）
- [ ] 6a. Firefox / WebKit を導入できれば同じ試験を実行（helpers.mjs の BROWSER 環境変数）。不可なら未検証と明記
- [ ] 6b. Dr/Rr 基準曲線の値と定義範囲を検索で確認
- [ ] 6c. 細部修正 → version 2.0.1 → run-all → REPORT/CHANGELOG 更新 → /home/claude/deliver に再zip → SendUserFile

## ユーザー報告の不具合5件（v2.0.1 に含める）— 進行中
1 短いWAV+16Hzで帯域分析が信号長で打切り(fitBands)→ゼロ詰めの末尾込み分析(bandTail)+最小録音長。 2 WAV失敗時にS.modeを先に書換→成功後に一括確定。
3 読込中の操作を上書き→intentEpoch。 4 TL/LH NaN/∞サンプル拒否+peakOfを非有限でInfinity+失敗時は元音源へ復帰。 5 TL pendingPrepが失敗結果を保持→settle時にクリア。
試験: tests/08-reported-bugs.mjs（未作成なら作る）→ run-all → REPORT/CHANGELOG → zip → SendUserFile

## 完了状態（v2.0.1）
- 報告5件を修正、tests/08 追加。run-all 301/301 PASS。REPORT 12章・CHANGELOG 2.0.1 更新済み。
- 6a Firefox/WebKit: ダウンロード遮断で導入不可→未検証のまま。6b Dr評価範囲125Hz-2kHzはGBRC資料で確認。6c 完了。
- 納品: /home/claude/deliver/lh-tl-listening-v2.0.1-project.zip

# ===== v2.1.0 作業（3回目の依頼: attachment.txt の11章＋途中ZIP必須）=====
土台: v2.0.1（src/ 現状）。元ZIPは orig/lh-tl-listening-v2.0.0-project.zip と /home/claude/deliver/ に保持。
途中ZIP: `bash tools/pack.sh <name>` → /home/claude/deliver/wip/<name>.zip（STATUS.md を必ず更新してから）。
マイルストーン: M0 現状ZIP / M1 TL切替 / M2 LH 1.9s・T0.7・再生中編集 / M4 表とLAFmax表示 / M3 可視化 / M5 暗騒音WAV / M6 全試験・REPORT・最終ZIP
各M完了時: build → 対象試験 → STATUS/PROGRESS/CHANGELOG 更新 → pack。
## 進捗
- [x] M0 WIP-0 作成
- [x] M1 TL切替（XF=10ms cosine, LEAD=10ms, onFrame整列, graph.lastEndで直列化）tests/04 42/42, tests/study-crossfade.mjs
- [x] M2 LH: RUN_INTERVAL 1.9, T既定0.7, play()でbg継続(rampTo), 同セルrun再押下no-op, liveUpdate()+requestRebuild(delay,thenPlay,live), 表記(≤−n, Δn), LAsum。tests/09-v210.mjs A節 16/16。全体 320/320
- [ ] M3 可視化（LH: 等級曲線群+分析値, TL: 等級曲線群+帯域スペクトル→tlProcessに帯域レベル追加が必要）
- [ ] M5 暗騒音WAV（S.bgSrcは保存しない/ bgWav変数, buildBgでWAV分岐, ループ継ぎ目クロスフェード）
- [ ] M6 09に B(可視化) C(bg WAV) D(TL表示=出力) 追加 → run-all → REPORT.md を v2.1.0 に更新 → 最終zip
- [x] M3 可視化・M5 暗騒音WAV 実装済み（tests/09 B,C,D 合格。WIP-1には未収録）
- [ ] ユーザー報告3件（WIP-1に対して）: (1) LH 停止直後のレベル上げで旧音源が増幅 → fadeEnd まで上げを遅延 (2) TL 停止が切替予約の後ろに並ぶ → 入れ子クロスフェード(mixノード)+kill段で即停止 (3) TL 低下量が短い音でずれる → 実際の出力バッファから統計
  → 修正後 tests/10-reported-2.mjs、run-all、WIP-2 zip、REPORT 更新、最終zip
- [x] 報告3件修正（lh.js quietAt/holdUntil、tl.js 入れ子mix+kill段、dsp.js measure()）tests/10 29/29。全体 381/381。WIP-2 作成。
- [ ] 残り: REPORT.md を v2.1.0 に更新、screenshots 再取得（node tests/screenshots.mjs && node tests/montage.mjs）、最終zip（tools/pack.sh lh-tl-listening-v2.1.0-project）
- [x] REPORT 13章、スクリーンショット更新、最終zip /home/claude/deliver/lh-tl-listening-v2.1.0-project.zip
