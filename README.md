# 最重要

ここのファイルは実用的にこういうことをやりたいと思いながら作ってみたファイルであって、あっているかどうかの保証はできません。
あくまでも「コンセプト」です。何があっても何の責任もとれません。

# やりたかったこと
- 環境依存しないものがほしかった。スマホやiPadからうごくものがほしいなぁ。ずっとパソコン持ち歩くの面倒だなぁというのがモチベーション。
- ここ10年ぐらいMatlabとPythonで使ってきたものをChatGPTやGeminiに「アップロードした～.pyを.htmlとjavascriptベースで動作するようにして」⇒デバッグしたものを公開。デバッグの内容は例えば[これ](HTMLfilterCheck.md)
- 仕事上のコミュニケーションをとる上であったらいいな、とおもったものを対象にしている。
- 測定のあと同じルーティーン。やることもおなじ。だったらGUI化してしまえー速報値はさっさとみたい！

[index](https://ssmhrkw.github.io/github.io_cv/index.html)

## 🛠️ 公開ツール (HTML Files)

| ツール名 | ファイル名 (リンク) | 主な機能/概要 |
| :--- | :--- | :--- |
| **必要スラブ厚の逆算（かんたん入力）** | [slab-thickness-quick.html](https://ssmhrkw.github.io/github.io_cv/slab-thickness-quick.html) | 目標 LH・スラブ・下室寸法から必要な等価スラブ厚を求めます。梁条件 81 通りの分布つき。 |
| **必要スラブ厚の逆算（詳細入力）** | [slab-thickness.html](https://ssmhrkw.github.io/github.io_cv/slab-thickness.html) | 同上の詳細版。計算書（式・代入値）を表示します。 |
| **Floor Impact** | [Floorimpact.html](https://ssmhrkw.github.io/github.io_cv/Floorimpact.html) | 重量床衝撃音レベルと LH 等級の予測（インピーダンス法）。 |
| **NA28 FIC** | [NA28_FIC.html](https://ssmhrkw.github.io/github.io_cv/NA28_FIC.html) | NA-28 の RND を読み込み、床衝撃音遮断性能を等級曲線で評価します。 |
| **RND Viewer** | [rnd_viewer.html](https://ssmhrkw.github.io/github.io_cv/rnd_viewer.html) | NA-28 の RND をグラフと表で確認します。 |
| **Measurement Positions** | [MeasurementPositions.html](https://ssmhrkw.github.io/github.io_cv/MeasurementPositions.html) | 受音点・加振点の配置計画（JIS A 1418、ISO 16283-2、ISO 10140-3）。 |
| **CSV Segmenter** | [CSV_Segmenter.html](https://ssmhrkw.github.io/github.io_cv/CSV_Segmenter.html) | 計測 CSV から打撃区間を検出・切り出して集計します。 |
| **wav2Fmax** | [wav2Fmax.html](https://ssmhrkw.github.io/github.io_cv/wav2Fmax.html) | WAV/CSV の 1/1・1/3 オクターブ分析（Fast/Slow 最大値、Leq、LN）。DA-20 の校正値を自動適用。 |
| **WAV Inspector** | [wav_Inspector.html](https://ssmhrkw.github.io/github.io_cv/wav_Inspector.html) | WAV のヘッダ・DA-20 メタ情報（CSET など）を表示します。 |
| **Revtime Recorder** | [revtime_recorder.html](https://ssmhrkw.github.io/github.io_cv/revtime_recorder.html) | 手拍子から残響時間を推定します（Lundeby 法のノイズ補償）。端末によっては勝手に信号処理されます。 |
| **Mic Calibration** | [calcorr.html](https://ssmhrkw.github.io/github.io_cv/calcorr.html) | 校正器の録音からチャンネルごとの換算係数を求めます。 |
| **Filter Design** | [FilterDesign.html](https://ssmhrkw.github.io/github.io_cv/FilterDesign.html) | Butterworth 1/n オクターブフィルタの設計と IEC 61260:1995 許容範囲の確認。 |
| **Plates** | [plates.html](https://ssmhrkw.github.io/github.io_cv/plates.html) | 矩形板の駆動点インピーダンス（有限板・無限板）。 |
| **Standing Wave** | [standing-wave.html](https://ssmhrkw.github.io/github.io_cv/standing-wave.html) | 単純支持板の速度/加振力分布をモード別に可視化。 |
| **Slab FDTD** | [SlabFDTD.html](https://ssmhrkw.github.io/github.io_cv/SlabFDTD.html) | 単純支持板の FDTD 時間応答とモード同定（Fdtdworker.js を使用）。 |
| **Laminated Calc** | [Laminated-Calc.html](https://ssmhrkw.github.io/github.io_cv/Laminated-Calc.html) | 単板・積層・二重壁の透過損失の参考計算。 |
| **TSEA** | [TSEA.html](https://ssmhrkw.github.io/github.io_cv/TSEA.html) | 過渡 SEA（TSEA）の計算。 |
| **集合住宅の音 30秒CM** | [cm/](https://ssmhrkw.github.io/github.io_cv/cm/) | 普及啓発用の短い動画。 |

GitHub Pages はファイル名の大文字・小文字を区別します（例：`wav_Inspector.html`）。

## 📄 ドキュメント (Markdown Files)

| ドキュメント名 | ファイル名 (リンク) | 概要 |
| :--- | :--- | :--- |
| **HTML Filter Check** | [HTMLfilterCheck.md](HTMLfilterCheck.md) | ブラウザ版フィルタのデバッグ記録。 |
| **acoustics-core** | [acoustics-core.md](acoustics-core.md) | 共通 DSP ライブラリ acoustics-core.js の説明。 |
| **Floor Impact** | [Floorimpact.md](Floorimpact.md) | Floorimpact.html の入力・計算式・判定。 |
| **Filter Design** | [FilterDesign.md](FilterDesign.md) | FilterDesign.html の補足。 |
| **NA28 FIC** | [NA28_FIC.md](NA28_FIC.md) | NA28_FIC.html の補足。 |
| **Plates** | [plates.md](plates.md) | plates.html の補足。 |
| **Revtime Recorder** | [revtime_recorder.md](revtime_recorder.md) | revtime_recorder.html の補足。 |
| **RND Viewer** | [rnd_viewer.md](rnd_viewer.md) | rnd_viewer.html の補足。 |
| **wav2Fmax** | [wav2Fmax.md](wav2Fmax.md) | wav2Fmax.html の補足。 |
| **WAV Inspector** | [wav_Inspector.md](wav_Inspector.md) | wav_Inspector.html の補足。 |
| **CSV Segmenter** | [CSV_Segmenter.md](CSV_Segmenter.md) | CSV_Segmenter.html の補足。 |

---

## ⚙️ コアファイル

-   **`index.html`**: トップページ。
-   **`acoustics-core.js`**: wav2Fmax・revtime_recorder・calcorr・FilterDesign が使う共通 DSP（Butterworth 帯域フィルタ、帯域中心周波数、時間重み付け）。帯域は IEC 61260-1:2014 の基数10で計算します。
-   **`Fdtdworker.js`**: SlabFDTD.html の計算部分（Web Worker）。
-   **`slab-thickness*.html`**: 計算部分をページ内に含めているので、1 ファイルだけでも動きます。
-   **Google Analytics**: CSV_Segmenter・Floorimpact・plates・rnd_viewer・wav_Inspector の各ページに入っています（index には入っていません）。
-   **`README.md`**: このファイル。

---
