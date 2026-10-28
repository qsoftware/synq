//
// Created by Adenilton on 08/02/26.
//

#include <complex>
#include "../include/qsd.h"

using namespace std;

qsd_result qsd(Eigen::MatrixXcd gate1, Eigen::MatrixXcd gate2) {
    /*
     * gate1 = V D W and gate2 = V D^dagger W, where V D^2 V^dagger = gate1 gate2^dagger.
     *
     * gate1 gate2^dagger is unitary (normal), so its complex Schur form is diagonal
     * and the Schur vectors are an orthonormal eigenbasis even when eigenvalues are
     * repeated. ComplexEigenSolver does not guarantee orthogonal eigenvectors in the
     * degenerate case (identity, permutations, diagonal gates, ...), which produced
     * non-unitary V and W.
     */
    const Eigen::MatrixXcd c = gate1 * gate2.adjoint();
    const Eigen::ComplexSchur<Eigen::MatrixXcd> schur(c);

    const Eigen::MatrixXcd& gate_v = schur.matrixU();
    const Eigen::VectorXcd d_square = schur.matrixT().diagonal();

    Eigen::VectorXcd d_vector(d_square.size());
    for (Eigen::Index i = 0; i < d_square.size(); ++i) {
        // unit-modulus square root (removes rounding drift of |d|)
        d_vector[i] = std::polar(1.0, std::arg(d_square[i]) / 2.0);
    }

    qsd_result b;
    b.gate_v = gate_v;
    b.list_d = d_vector;
    b.gate_w = d_vector.asDiagonal() * gate_v.adjoint() * gate2;
    return b;
}
