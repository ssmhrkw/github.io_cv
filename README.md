# 最重要
ここのファイルは実用的にこういうことをやりたいと思いながら作ってみたファイルであって、あっているかどうかの保証はできません。あくまでも「コンセプト」です。
何があっても何の責任もとれません。

# やりたかったこと
- 環境依存しないものがほしかった。スマホやiPadからうごくものがほしい。
- 測定のあと同じルーティーン。やることもおなじ。だったらGUI化してしまえ。
- ここ１０年ぐらいMatlabとPythonで使ってきたものをChatGPTやGeminiに～.pyを.htmlにしてjavascriptベースで動作するようにして⇒デバッグしたものを公開
- 仕事上のコミュニケーションをとる上であったらいいな、とおもったものを対象。

[index](https://ssmhrkw.github.io/github.io_cv/index.html).
## 🛠️ 公開ツール (HTML Files)

| ツール名 | ファイル名 (リンク) | 主な機能/概要 |
| :--- | :--- | :--- |
| **CSV Segmenter** | [CSV_Segmenter.html](CSV_Segmenter.html) | CSVファイルのセグメンテーション（区切り、分割）を行うツールです。 |
| **Filter Design** | [FilterDesign.html](FilterDesign.html) | フィルタ設計を行うツールです。（HTML選択オプションのフォーマット修正履歴あり） |
| **Floor Impact** | [Floorimpact.html](Floorimpact.html) | 床衝撃音の解析、表示に関するツールです。 |
| **Laminated Calc** | [Laminated-Calc.html](Laminated-Calc.html) | 積層構造に関する計算を行うツールです。 |
| **NA28 FIC** | [NA28_FIC.html](NA28_FIC.html) | NA28_FIC規格に関連する解析ツールです。 |
| **Plates** | [plates.html](plates.html) | プレート（板）の解析に関連するツールです。 |
| **Revtime Recorder** | [revtime_recorder.html](revtime_recorder.html) | 残響時間測定・記録を行うツールです。 |
| **RND Viewer** | [rnd_viewer.html](rnd_viewer.html) | RNDファイルの内容を表示・確認するビューアツールです。 |
| **wav2Fmax** | [wav2Fmax.html](wav2Fmax.html) | WAVファイルを解析し、Fmax（最大周波数）などの値を取得するツールです。 |
| **WAV Inspector** | [wav_inspector.html](wav_inspector.html) | WAVファイルの情報を詳細に検査するツールです。 |
| **Cal Corr** | [calcorr.html](calcorr.html) | キャリブレーション（校正）や補正に関する計算ツールです。 |

## 📄 ドキュメント (Markdown Files)

ツールやライブラリに関する説明文書です。

| ドキュメント名 | ファイル名 (リンク) | 概要 |
| :--- | :--- | :--- |
| **HTML Filter Check** | [HTMLFilterCheck.md](HTMLFilterCheck.md) | HTMLフィルタのチェックに関するドキュメントです。 |
| **NA28 FIC** | [NA28_FIC.md](NA28_FIC.md) | NA28_FICツールに関する補足ドキュメントです。 |
| **acoustics-core** | [acoustics-core.md](acoustics-core.md) | 音響解析の核となる`acoustics-core.js`ライブラリに関するドキュメントです。 |
| **Plates** | [plates.md](plates.md) | プレート解析ツールの補足ドキュメントです。 |
| **Revtime Recorder** | [revtime_recorder.md](revtime_recorder.md) | 残響時間レコーダーに関する補足ドキュメントです。 |
| **RND Viewer** | [rnd_viewer.md](rnd_viewer.md) | RNDビューアに関する補足ドキュメントです。 |
| **WAV Inspector** | [wav_inspector.md](wav_inspector.md) | WAVインスペクターに関する補足ドキュメントです。 |

---

## ⚙️ コアファイル

-   **`index.html`**: トップページ。Google Analyticsトラッキングスクリプトが統合されています。
-   **`acoustics-core.js`**: 各ツールで利用される音響計算のコアとなるJavaScriptライブラリです。
-   **`README.md`**: このファイル。フォーマット改善の履歴があります。

---
