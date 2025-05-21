/**
 * @file main.cpp
 * @brief FastICA implementation test
 * @author Based on Kikuchi Tomoo's implementation
 * @date 2025-05-21
 */

#include <Eigen/Dense>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>

#include "fast_ica.cpp"  // FastICAクラスのインクルード

// ファイル出力関数の定義
void output_file(const char* fpath, double datas[3][2048], double t[2048]) {
    std::ofstream outFile(fpath);
    if (!outFile.is_open()) {
        std::cerr << "Error: Could not open file " << fpath << std::endl;
        return;
    }

    // データをスペース区切りで出力（rad data0 data1 data2の形式）
    for (int i = 0; i < 2048; ++i) {
        outFile << t[i] << " " << datas[0][i] << " " << datas[1][i] << " "
                << datas[2][i] << "\n";
    }

    outFile.close();
}

int main() {
    // パラメータ設定
    const int numSamples = 2048;    // サンプル数
    const int rows = 3;             // 入力データのシグナル数
    const int numSignals = 3;       // 抽出したい独立成分の数
    const int maxIterations = 100;  // 最大反復回数
    const double tolerance = 1e-6;  // 収束判定の閾値

    // 乱数生成器の初期化（ホワイトノイズ用）
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> dist(0.0, 1.0);

    // 入力データの準備
    double t[2048];
    double original[3][2048];
    double input[3][2048];
    double recomp[3][2048];  // 復元された信号用

    // 0から6πまでのt値を生成
    for (int i = 0; i < 2048; ++i) {
        t[i] = 6.0 * M_PI * i / 2047.0;
    }

    // 元信号の生成
    for (int i = 0; i < 2048; ++i) {
        original[0][i] = std::sin(t[i]);         // sin(t)
        original[1][i] = std::sin(t[i] * 2.0);   // sin(2t)
        original[2][i] = dist(gen) * 2.0 - 1.0;  // ホワイトノイズ
    }

    // 混合信号の生成
    for (int i = 0; i < 2048; ++i) {
        input[0][i] = (original[0][i] * 0.2) + (original[1][i] * 0.3) +
                      (original[2][i] * 0.5);
        input[1][i] = (original[0][i] * 0.5) + (original[1][i] * 0.1) +
                      (original[2][i] * 0.4);
        input[2][i] = (original[0][i] * 0.1) + (original[1][i] * 0.6) +
                      (original[2][i] * 0.3);
    }

    // 混合前の信号と混合後の信号をファイルに出力
    output_file("original.dat", original, t);
    output_file("mixed.dat", input, t);

    // Eigen行列に変換
    Eigen::MatrixXd X(rows, numSamples);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < numSamples; ++j) {
            X(i, j) = input[i][j];
        }
    }

    // FastICAクラスを使用して独立成分分析を実行
    KT::FastICA fastICA;
    fastICA.set_max_iterations(maxIterations);
    fastICA.set_tolerance(tolerance);

    // 独立成分の抽出
    Eigen::MatrixXd S = fastICA.transform(X, numSignals);

    // 混合行列の取得
    Eigen::MatrixXd A = fastICA.get_mixing_matrix();

    // 結果の出力
    std::cout << "\n独立成分の先頭10サンプル:" << std::endl;
    std::cout << S.leftCols(10) << std::endl;

    std::cout << "\n推定された混合行列:" << std::endl;
    std::cout << A << std::endl;

    // 元の混合行列（理論値）
    std::cout << "\n理論上の混合行列:" << std::endl;
    std::cout << "0.2 0.3 0.5\n0.5 0.1 0.4\n0.1 0.6 0.3" << std::endl;

    // 相関係数の計算で検証
    std::cout << "\n抽出された独立成分と元信号の相関:" << std::endl;

    // 元信号をEigen形式に変換
    Eigen::MatrixXd Original(rows, numSamples);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < numSamples; ++j) {
            Original(i, j) = original[i][j];
        }
    }

    // 各独立成分と元信号の相関を計算
    for (int i = 0; i < numSignals; ++i) {
        std::cout << "Independent Component " << i + 1 << ":" << std::endl;
        for (int j = 0; j < rows; ++j) {
            // 相関係数の計算（符号が重要、大きさは信号のスケーリングにより異なる）
            Eigen::VectorXd ic = S.row(i);
            Eigen::VectorXd orig = Original.row(j);

            // 正規化
            ic = (ic.array() - ic.mean()) / ic.norm();
            orig = (orig.array() - orig.mean()) / orig.norm();

            double correlation = std::abs(ic.dot(orig));
            std::cout << "  with Original Signal " << j + 1 << ": "
                      << correlation << std::endl;
        }
    }

    // 復元信号を配列に変換
    for (int i = 0; i < numSignals; ++i) {
        for (int j = 0; j < numSamples; ++j) {
            recomp[i][j] = S(i, j);
        }
    }

    // 復元信号をファイルに出力
    output_file("recomp.dat", recomp, t);

    std::cout << "\n出力ファイル:" << std::endl;
    std::cout << "元信号: original.dat" << std::endl;
    std::cout << "混合信号: mixed.dat" << std::endl;
    std::cout << "復元信号: recomp.dat" << std::endl;

    return 0;
}