/**
 * @file fast_ica.cpp
 * @brief Fast Independent Component Analysis implementation
 * @author Kikuchi Tomoo
 * @date 2025-05-21
 * @details
 * @copyright Kikuchi Tomoo 2025 Copyright All right reserved.
 */

#ifndef __KIKUCHI_FAST_ICA
#define __KIKUCHI_FAST_ICA

#define DEFAULT_MAX_ITERATION (1000)
#define DEFAULT_TOLERANCE (1e-6)

#include <Eigen/Dense>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>

namespace KT {
class FastICA {
   private:
    int rows_;
    int cols_;
    int sigs_;
    int max_iter_;
    double tole_;

    Eigen::MatrixXd x_;
    Eigen::MatrixXd y_;

    Eigen::RowVectorXd mean_;
    Eigen::MatrixXd conv_;
    Eigen::VectorXd eigen_;
    Eigen::MatrixXd e_;
    Eigen::MatrixXd w_;
    Eigen::MatrixXd z_;
    Eigen::MatrixXd d_sq_inv_;
    Eigen::MatrixXd transform_;

    void __centering() {
        mean_ = x_.rowwise().mean();
        for (int i = 0; i < rows_; i++) {
            x_.row(i).array() -= mean_(i);
        }
    }

    void __whitening() {
       
        conv_ = (x_ * x_.transpose()) / (cols_ - 1);

        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen_solver(conv_);
        eigen_ = eigen_solver.eigenvalues();
        e_ = eigen_solver.eigenvectors();

        d_sq_inv_ = Eigen::MatrixXd::Zero(rows_, rows_);
        for (int i = 0; i < rows_; i++) {
            d_sq_inv_(i, i) = 1.0 / std::sqrt(eigen_(i));
        }

        // ホワイトニング
        Eigen::MatrixXd whitening_mat_ = d_sq_inv_ * e_.transpose();
        z_ = whitening_mat_ * x_;
    }

    void __ica() {
        w_ = Eigen::MatrixXd(sigs_, rows_);

        for (int c = 0; c < sigs_; c++) {
            Eigen::VectorXd w = Eigen::VectorXd::Random(rows_);
            w.normalize();

            for (int iter = 0; iter < max_iter_; iter++) {
                Eigen::VectorXd w_prv = w;

                // g(x) = tanh(x) を採用
                Eigen::VectorXd wx = w.transpose() * z_;
                Eigen::VectorXd g_wx(cols_);
                Eigen::VectorXd g_p_wx(cols_);

                for (int j = 0; j < cols_; j++) {
                    // tanh'(x) = 1 - tanh^2(x)
                    g_wx(j) = std::tanh(wx(j));
                    g_p_wx(j) = 1.0 - std::pow(std::tanh(wx(j)), 2);
                }

                // ウェイトベクトルの更新
                w = (z_ * g_wx) / cols_ - (g_p_wx.mean() * w);

                // ベクトルの直交化
                if (c > 0) {
                    Eigen::MatrixXd proj = Eigen::MatrixXd::Zero(rows_, rows_);
                    for (int i = 0; i < c; i++) {
                        proj += w_.row(i).transpose() * w_.row(i);
                    }
                    w = w - proj * w;
                }

                w.normalize();

                double d = std::abs(1.0 - std::abs(w.dot(w_prv)));
                if (d < tole_) {
                    break;
                }

                if (iter == max_iter_ - 1) {
                    // TODO: warn log
                    std::cerr << "Component " << c + 1 << " did not converge."
                              << std::endl;
                }
            }

            w_.row(c) = w.transpose();
        }
    }

    void __transform() {
        y_ = w_ * z_;

        Eigen::MatrixXd d_sq = Eigen::MatrixXd::Zero(rows_, rows_);
        for (int i = 0; i < rows_; i++) {
            d_sq(i, i) = std::sqrt(eigen_(i));
        }

        transform_ = e_ * d_sq * w_.inverse();
    }

   public:
    FastICA() {
        max_iter_ = DEFAULT_MAX_ITERATION;
        tole_ = DEFAULT_TOLERANCE;
    }

    ~FastICA() {}

    void set_max_iterations(int _max_iter) { max_iter_ = _max_iter; }

    void set_tolerance(double _tole) { tole_ = _tole; }

    Eigen::MatrixXd transform(Eigen::MatrixXd _input, int _sigs) {
        x_ = _input;

        rows_ = x_.rows();
        cols_ = x_.cols();

        sigs_ = _sigs;

        if (rows_ < sigs_) {
            std::cerr << "Error: Number of signals to extract cannot exceed "
                         "number of input rows."
                      << std::endl;
            return Eigen::MatrixXd::Zero(0, 0);  // 空の行列を返す
        }

        __centering();
        __whitening();
        __ica();
        __transform();

        return y_;
    }

    // 混合行列を取得するメソッドを追加
    Eigen::MatrixXd get_mixing_matrix() { return transform_; }
};
};  // namespace KT

#endif
