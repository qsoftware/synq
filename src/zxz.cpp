#include "../include/zxz.h"

#include <complex>

#ifdef USE_MKL
    #include <mkl_lapacke.h>
#else
    #define lapack_complex_double std::complex<double>
    #include <lapacke.h>
#endif

namespace {

using Eigen::MatrixXcd;

/** Full SVD a = u diag(s) v†. */
struct svd_result {
    MatrixXcd u;
    MatrixXcd v;
};

svd_result svd(const MatrixXcd& a) {
    // LAPACK zgesdd (divide and conquer): several times faster than Eigen::BDCSVD.
    const auto n = static_cast<lapack_int>(a.rows());
    MatrixXcd work = a;
    MatrixXcd u(n, n);
    MatrixXcd vt(n, n);
    Eigen::VectorXd s(n);
    const auto info = LAPACKE_zgesdd(LAPACK_COL_MAJOR, 'A', n, n,
                                     reinterpret_cast<lapack_complex_double*>(work.data()), n, s.data(),
                                     reinterpret_cast<lapack_complex_double*>(u.data()), n,
                                     reinterpret_cast<lapack_complex_double*>(vt.data()), n);
    if (info == 0) {
        return {u, vt.adjoint()};
    }
    const Eigen::BDCSVD<MatrixXcd> eigen_svd(a, Eigen::ComputeFullU | Eigen::ComputeFullV);
    return {eigen_svd.matrixU(), eigen_svd.matrixV()};
}

}  // namespace

zxz_result blockZxz(const MatrixXcd& u) {
    const Eigen::Index h = u.rows() / 2;
    const MatrixXcd x = u.topLeftCorner(h, h);
    const MatrixXcd y = u.topRightCorner(h, h);

    // X = V S W†
    const svd_result svd_x = svd(x);
    const MatrixXcd& v = svd_x.u;
    const MatrixXcd& w = svd_x.v;

    // Unitary polar factor of Y†V = Q K. It is not unique when K is singular (for
    // example when Y = 0); any choice satisfies Y = V K Q†.
    const svd_result svd_m = svd(y.adjoint() * v);
    const MatrixXcd q = svd_m.u * svd_m.v.adjoint();

    const MatrixXcd c_adj = std::complex<double>(0.0, 1.0) * q * w.adjoint();

    zxz_result r;
    r.a1 = x + y * c_adj;
    r.a2 = u.bottomLeftCorner(h, h) + u.bottomRightCorner(h, h) * c_adj;
    r.b = 2.0 * r.a1.adjoint() * x - MatrixXcd::Identity(h, h);
    r.c = c_adj.adjoint();
    return r;
}
