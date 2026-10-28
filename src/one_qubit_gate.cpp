//
// Created by adenilton on 10/3/25.
//

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "../include/one_qubit_gate.h"
#include <complex>

zyz_result OneQubit::zyz_decomposition(Eigen::Matrix2cd uMatrix) {
    /*
     * U = e^(i alpha) V, with V in SU(2):
     *   V = [ e^{-i(beta+delta)/2} cos(gamma/2)   -e^{-i(beta-delta)/2} sin(gamma/2) ]
     *       [ e^{ i(beta-delta)/2} sin(gamma/2)    e^{ i(beta+delta)/2} cos(gamma/2) ]
     * so (beta+delta)/2 = arg(V11) and (beta-delta)/2 = arg(V10).
     *
     * The previous implementation took the phases of U00, U01 and U10, which are
     * undefined when an entry is zero (identity, X, diagonal gates, ...).
     * Here, when V11 (or V10) is zero, std::arg returns 0 and the corresponding
     * combination of angles is irrelevant, so the result is still exact.
     */
    const std::complex<double> det = uMatrix.determinant();
    const double alpha = std::arg(det) / 2.0;
    const Eigen::Matrix2cd v = uMatrix * std::exp(std::complex<double>(0.0, -alpha));

    const double gamma = 2.0 * std::atan2(std::abs(v(1, 0)), std::abs(v(0, 0)));
    const double half_sum = std::arg(v(1, 1));   // (beta + delta) / 2
    const double half_diff = std::arg(v(1, 0));  // (beta - delta) / 2
    const double beta = half_sum + half_diff;
    const double delta = half_sum - half_diff;
    return {alpha, beta, gamma, delta};
}

Eigen::Matrix2cd OneQubit::h_matrix() {
    Eigen::Matrix2cd const h =
        (Eigen::Matrix2cd() <<  1.0/sqrt(2.0)+0.0i, 1.0/sqrt(2.0)+0.0i,
                                1.0/sqrt(2.0)+0.0i, -1.0/sqrt(2.0)+0.0i).finished();
    return h;
}

Eigen::Matrix2cd OneQubit::x_matrix() {
    Eigen::Matrix2cd const x =
        (Eigen::Matrix2cd() << 0.0 + 0.0i, 1.0 + 0.0i,
                                  1.0 + 0.0i, 0.0 + 0.0i).finished();
    return x;
}

Eigen::Matrix2cd OneQubit::ry_matrix(double theta) {
    Eigen::Matrix2cd const ry =
        (Eigen::Matrix2cd() << cos(theta/2)+0.i, -sin(theta/2)+0.i,
                               sin(theta/2)+0.i, cos(theta/2)+0.i).finished();
    return ry;
}

Eigen::Matrix2cd OneQubit::rz_matrix(double theta) {
    Eigen::Matrix2cd const rz =
        (Eigen::Matrix2cd() << exp(-theta*1.i/2.0), 0.i,
                               0.i,              exp(theta*1.i/2.0)).finished();
    return rz;
}
