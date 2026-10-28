#include "../include/zxz.h"

#include <complex>

zxz_result blockZxz(const Eigen::MatrixXcd& u) {
    using Eigen::MatrixXcd;
    const Eigen::Index h = u.rows() / 2;
    const MatrixXcd x = u.topLeftCorner(h, h);
    const MatrixXcd y = u.topRightCorner(h, h);

    // X = V S W†
    const Eigen::BDCSVD<MatrixXcd> svd_x(x, Eigen::ComputeFullU | Eigen::ComputeFullV);
    const MatrixXcd& v = svd_x.matrixU();
    const MatrixXcd& w = svd_x.matrixV();

    // Unitary polar factor of Y†V = Q K. It is not unique when K is singular (for
    // example when Y = 0); any choice satisfies Y = V K Q†.
    const Eigen::BDCSVD<MatrixXcd> svd_m(y.adjoint() * v, Eigen::ComputeFullU | Eigen::ComputeFullV);
    const MatrixXcd q = svd_m.matrixU() * svd_m.matrixV().adjoint();

    const MatrixXcd c_adj = std::complex<double>(0.0, 1.0) * q * w.adjoint();

    zxz_result r;
    r.a1 = x + y * c_adj;
    r.a2 = u.bottomLeftCorner(h, h) + u.bottomRightCorner(h, h) * c_adj;
    r.b = 2.0 * r.a1.adjoint() * x - MatrixXcd::Identity(h, h);
    r.c = c_adj.adjoint();
    return r;
}
