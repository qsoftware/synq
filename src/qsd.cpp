//
// Created by Adenilton on 08/02/26.
//

#include <complex>
#include "../include/qsd.h"

using namespace std;

qsd_result qsd(Eigen::MatrixXcf gate1, Eigen::MatrixXcf gate2) {

    Eigen::MatrixXcf c = gate1 * gate2.adjoint();
    Eigen::ComplexEigenSolver<Eigen::MatrixXcf> ces(c);

    auto gate_v = ces.eigenvectors();
    auto d_square = ces.eigenvalues();
    auto eig_list_d = d_square.array().sqrt();
    Eigen::VectorXcf d_vector(d_square.size());
    std::vector<complex<float>> list_d2(d_square.size());
    for (int i = 0; i < d_square.size(); ++i) {
        list_d2[i] = eig_list_d[i];
        d_vector[i] = eig_list_d[i];
    }
    qsd_result b;
    b.left_gate = gate_v;
    b.list_d = d_vector;

    b.right_gate = d_vector.asDiagonal() *  gate_v.adjoint() * gate2;
    return b;
}