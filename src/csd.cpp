//
// Created by adenilton on 10/10/25.
//

#include "../include/csd.h"
#include<iostream>

csd_result csd(const Eigen::MatrixXcd& unitary, Eigen::Index p, Eigen::Index q) {
    auto m = unitary.rows();
    auto n = unitary.cols();

    const auto X11 = eigen2lapack(unitary, 0, 0, p, p);
    const auto X21 = eigen2lapack(unitary, 0, q, p, n-p);
    const auto X12 = eigen2lapack(unitary, p, 0,m-p, q);
    const auto X22 = eigen2lapack(unitary, p, q, m-p, n-q);

    const int r = std::min({p, q, m-p, m-q});
    std::vector<double> theta(r);

    std::vector<lapack_complex_double> V1T(q * q);
    std::vector<lapack_complex_double> V2T((m - q) * (m - q));
    std::vector<lapack_complex_double> U1(p * p);
    std::vector<lapack_complex_double> U2((m - p) * (m - p));

    const int ldx11 = p;
    const int ldx12 = p;
    const int ldx21 = m-p;
    const int ldx22 = m-p;
    const int ldu1 = p;
    const int ldu2 = m-p;
    const int ldv1t = q;
    const int ldv2t = m-q;

    const int info = LAPACKE_zuncsd(
        LAPACK_COL_MAJOR, // matrix_layout
        'Y', // JOBU1
        'Y', // JOBU2
        'Y', // JOBV1T
        'Y', // JOBV2T
        'T', // TRANS
        'O', // SIGNS
        m, p, q, // M, P, Q
        X11, ldx11, // X11, ldx11
        X12, ldx12, // x12, ldx12
        X21, ldx21, // x21, ldx21
        X22, ldx22, // x11, ldx11
        theta.data(),
        U1.data(), ldu1,
        U2.data(), ldu2,
        V1T.data(), ldv1t,
        V2T.data(), ldv2t // x22, ldx22
    );

    if (info != 0) std::cerr << "LAPACKE_cuncsd: " << info << std::endl;

    const Eigen::Map<Eigen::MatrixXcd> out_U1(reinterpret_cast<std::complex<double>*> (U1.data()), p, p);
    const Eigen::Map<Eigen::MatrixXcd> out_U2(reinterpret_cast<std::complex<double>*> (U2.data()), m-p, m-p);
    const Eigen::Map<Eigen::MatrixXcd> out_V1T(reinterpret_cast<std::complex<double>*> (V1T.data()), p, p);
    const Eigen::Map<Eigen::MatrixXcd> out_V2T(reinterpret_cast<std::complex<double>*> (V2T.data()), m-p, m-p);
    const std::vector<double> out_theta(theta);

    csd_result out;
    out.U1 = out_U1;
    out.U2 = out_U2;
    out.V1T = out_V1T;
    out.V2T = out_V2T;
    out.theta = out_theta;
    delete X11;
    delete X12;
    delete X21;
    delete X22;

    return out;
}

lapack_complex_double* eigen2lapack(const Eigen::MatrixXcd& unitary, const int start_line, const int start_col, const int n_rows, const int n_cols) {
    lapack_complex_double* U1 = new lapack_complex_double[n_rows * n_cols];

    for (Eigen::Index j = 0; j < n_cols; ++j)
        for (Eigen::Index i = 0; i < n_rows; ++i) {
            const double a = unitary(i + start_line, j+start_col).real();
            const double b = unitary(i + start_line, j+start_col).imag();
            U1[n_rows*j +i] = complexD(a, b);
        }
    return U1;
}

bool verify(const Eigen::MatrixXcd& unitary, csd_result result) {

    const int n_rows = unitary.rows();
    const int n_cols = unitary.cols();
    const int p = result.U1.rows();
    const int q = result.V1T.rows();

    Eigen::MatrixXcd U = Eigen::MatrixXcd::Zero(n_rows, n_cols);
    U.block(0, 0, p, p) = result.U1;
    U.block(p, p, n_rows-p, n_rows-p) = result.U2;

    Eigen::MatrixXcd V = Eigen::MatrixXcd::Zero(n_rows, n_cols);
    V.block(0, 0, q, q) = result.V1T;
    V.block(q, q, n_rows-q, n_rows-q) = result.V2T;


    const size_t r = result.theta.size();
    Eigen::MatrixXcd S = Eigen::MatrixXcd::Zero(r, r);
    Eigen::MatrixXcd C = Eigen::MatrixXcd::Zero(r, r);
    Eigen::MatrixXcd sigma(2*r, 2*r);
    for (int i = 0; i < r; ++i) {
        S(i, i) = sin(result.theta[i]);
        C(i, i) = cos(result.theta[i]);
    }

    sigma << C, -S, S, C;

    Eigen::MatrixXcd out = V * sigma * U;
    if (out.isApprox(unitary)) return true;
    return false;
}