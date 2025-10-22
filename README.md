# 最重要
ここのファイルは実用的にこういうことをやりたいと思いながら作ってみたファイルであって、あっているかどうかの保証はできません。あくまでも「コンセプト」です。
何があっても何の責任もとれません。

# やりたかったこと
- 環境依存しないものがほしかった。スマホやiPadからうごくものがほしい。
- 測定のあと同じルーティーン。やることもおなじ。だったらGUI化してしまえ。
- ここ１０年ぐらいMatlabとPythonで使ってきたものをChatGPTやGeminiに～.pyを.htmlにしてjavascriptベースで動作するようにして⇒デバッグしたものを公開
- 仕事上のコミュニケーションをとる上であったらいいな、とおもったものを対象。

[index](https://ssmhrkw.github.io/github.io_cv/index.html).

### Revtime_recorder
A[録音/ファイル読込] <br>
B[帯域通過 SoS BPF designButterworthBandpassSOS_N + sosFilter]
C[シュレーダー積分<br/>(BG減算オプション)]
D[回帰: EDT/T10/T20/T30]
E[SNR計算<br/>前0.5s vs 直後0.2s]
F[代表回帰選択<br/>(T30→T20→T10→EDT)]
G[帯域ごとOK/NG]
H[表とプロット更新]
I[全体判定]
