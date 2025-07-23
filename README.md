# FastICA
## 概要
C++11でFastICAを実行するプログラム

[こちらの記事](https://qiita.com/kikikitom/items/aa2d70a6d878aefc02ae)で書いたもののプログラムです

## 使用法
headerライブラリのため，src/fast_ica.cppをインクルードするだけで構いません.
計算ライブラリとしてEigenを使用しているため，コンパイルの際には-Iにeigenのインクルードパスを指定のこと．

具体的な使用方法はsrc/main.cppを見てください．

## デモ

gnuplotが必要です．

```bash
make
./main
plot plot.gpl
```
