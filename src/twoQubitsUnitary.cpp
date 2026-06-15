//
// Created by Alex on 06/05/26.
//

#include <complex>
#include "../include/twoQubitsUnitary.h"
#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include <numbers>

using namespace std;

using twoQubitsTensor = std::pair<Eigen::Matrix2cd, Eigen::Matrix2cd>;

complex<double> i(0,1);
double pi = std::numbers::pi;

Eigen::Matrix4cd E = (Eigen::Matrix4cd() << 
    1,  i,  0,  0,
    0,  0,  i,  1,
    0,  0,  i, -1,
    1, -i,  0,  0
).finished() * (1.0 / std::sqrt(2.0));

Eigen::Matrix4cd SWAP = (Eigen::Matrix4cd() << 
    1, 0, 0, 0,
    0, 0, 1, 0,
    0, 1, 0, 0,
    0, 0, 0, 1
).finished();

Eigen::Matrix2cd Rz(double theta) {
    Eigen::Matrix2cd mat;
    std::complex<double> exp_minus = std::exp(-i * (theta / 2.0));
    std::complex<double> exp_plus  = std::exp(i * (theta / 2.0));
    mat << exp_minus, 0.0,
           0.0, exp_plus;
    return mat;
}

twoQubitsTensor decomposeKronecker(const Eigen::Matrix4cd& M) {

    /** 
     * Kronecker product decomposition of a 4x4 matrix M into A ⊗ B, where A and B are 2x2 matrices.
     * Based on "Approximation with Kronecker Products" by Van Loan and Pitsianis (1993).
     **/
 
    Eigen::Matrix4cd R;
    
    R.row(0) << M(0,0), M(1,0), M(0,1), M(1,1); 
    R.row(1) << M(2,0), M(3,0), M(2,1), M(3,1);
    R.row(2) << M(0,2), M(1,2), M(0,3), M(1,3);
    R.row(3) << M(2,2), M(3,2), M(2,3), M(3,3);

    Eigen::JacobiSVD<Eigen::Matrix4cd> svd(R, Eigen::ComputeFullU | Eigen::ComputeFullV);
    
    double sigma0 = svd.singularValues()(0);
    Eigen::Vector4cd u0 = svd.matrixU().col(0);
    Eigen::Vector4cd v0 = svd.matrixV().col(0).conjugate();

    Eigen::Matrix2cd A, B;
    
    A << u0(0), u0(2),
         u0(1), u0(3);
         
    B << v0(0), v0(2),
         v0(1), v0(3);

    A *= std::sqrt(sigma0);
    B *= std::sqrt(sigma0);

    return {A, B};
}

std::pair<Eigen::Matrix4cd, Eigen::Matrix4cd> diagonalizes(const Eigen::Matrix4cd& SUm){

    /** Diagonalizes M = SUm * SUm^T while forcing the eigenvectors (op1) into SO(4, R).
     * A unitary in SU(4) can be decomposed in op1 * diag * op2, with op1 and op2 into SO(4, R).
     * based on "Nonlocal properties of two-qubit gates and mixed states and optimization of quantum computations"
     * by Makhlin (2000).
     **/

    Eigen::Matrix4cd M = SUm * SUm.transpose();

    Eigen::ComplexEigenSolver<Eigen::Matrix4cd> solver (M, false);
    Eigen::Vector4cd eigenvalues = solver.eigenvalues();
    Eigen::Matrix4d A = M.real();
    Eigen::Matrix4d B = M.imag();
    Eigen::Matrix4d I = Eigen::Matrix4d::Identity();

    std::vector<std::complex<double>> unique_eigenvalues;
    bool visited[4] = {false, false, false, false};

    //collects unique eigenvalues of M, which can be degenerate.
    for (int i = 0; i<4; i++){

        if (visited[i]) continue;

        visited[i] = true;
        std::complex<double> sigma = eigenvalues(i);
        unique_eigenvalues.push_back(sigma);
        for (int j = i+1; j<4; j++){
            if (std::abs(sigma - eigenvalues(j)) < 1e-10 && !visited[j]){
                visited[j] = true;
            }
        }
    }

    Eigen::Matrix4d real_op1;
    int col_idx = 0;

    // Hk * v = 0, where Hk = (A - ak * I)^2 + (B - bk * I)^2
    // and where ak and bk are the real and imaginary parts of the k-th eigenvalue of M
    // forces the eigenvectors v to be real.
    for (int k = 0; k < unique_eigenvalues.size(); k++) {
        double ak = unique_eigenvalues[k].real();
        double bk = unique_eigenvalues[k].imag();
        Eigen::Matrix4d Hk = (A - ak * I) * (A - ak * I) + (B - bk * I) * (B - bk * I);

        // resolves the linear system by using svd
        Eigen::JacobiSVD<Eigen::Matrix4d> svd(Hk, Eigen::ComputeFullV);
        Eigen::Vector4d singular_values = svd.singularValues();

        for (int i = 3; i >=0; i--) {
            if (singular_values(i) < 1e-10) {
                real_op1.col(col_idx) = svd.matrixV().col(i);
                col_idx++;
            }
        }
    }

    //guarantees op1 in SO(4)
    if (real_op1.determinant() < 0.0) {
        real_op1.col(0) *= -1.0;
    }

    Eigen::Matrix4cd op1 = real_op1.cast<std::complex<double>>();

    // W = op1^T * SUm and W^T * W = diag^2
    Eigen::Matrix4cd W = op1.transpose() * SUm;
    Eigen::Matrix4cd W_WT = W * W.transpose();
    Eigen::Vector4cd delta_diag(4);

    for (int k = 0; k < 4; k++) {
        delta_diag(k) = std::sqrt(W_WT(k, k));
        delta_diag(k) /= std::abs(delta_diag(k));
    }

    std::complex<double> prod = delta_diag(0) * delta_diag(1) * delta_diag(2) * delta_diag(3);
    
    // Delta with det = 1
    if (prod.real() < 0.0) {
        delta_diag(3) *= -1.0; 
    }

    Eigen::Matrix4cd diag = delta_diag.asDiagonal();
    return {op1, diag};
    
}

std::tuple<twoQubitsTensor, Eigen::Matrix4cd, twoQubitsTensor> su4Decomposition(const Eigen::Matrix4cd& specialUnitary) {

    // Puts specialUnitary in magic basis
    Eigen::Matrix4cd SUm = E.adjoint() * specialUnitary * E;

    auto [op1, diag] = diagonalizes(SUm);

    Eigen::Matrix4cd complexOp2 = diag.adjoint() * op1.transpose() * SUm;    
    Eigen::Matrix4d realOp2 = complexOp2.real();
    Eigen::Matrix4cd op2 = realOp2.cast<std::complex<double>>();
    
    // Puts M1 and M2 in canonical basis
    Eigen::Matrix4cd M1 = E * op1 * E.adjoint();
    Eigen::Matrix4cd M2 = E * op2 * E.adjoint();

    // Matrices in SO(4) can be decomposed into Local qubit operations by Makhlin (2000).
    twoQubitsTensor tensor1 = decomposeKronecker(M1);
    twoQubitsTensor tensor2 = decomposeKronecker(M2);
    
    return std::make_tuple(tensor1, diag, tensor2);
}

DecompTwoQubitsResult twoQubitsUnitaryDecomposition(const Eigen::Matrix4cd& U) {

    /** Decomposes a two-qubit unitary based on "Minimal universal two-qubit controlled CNOT-based circuits" 
     * by Shende, Markov, and Bullock (2004) **/

    std::complex<double> det = U.determinant();
    double extracted_phase = std::arg(det) / 4.0;
    Eigen::Matrix4cd SU4_U = U * std::exp(-i * extracted_phase);

    std::complex<double> exp_phase = std::exp(i * pi / 4.0);
    Eigen::Matrix4cd specialU = exp_phase * SWAP * SU4_U;

    auto [tensor1, diag, tensor2] = su4Decomposition(specialU);
    Eigen::Matrix2cd A = tensor1.first;
    Eigen::Matrix2cd B = tensor1.second * Rz(pi/2).adjoint();
    Eigen::Matrix2cd C = tensor2.first;
    Eigen::Matrix2cd D = Rz(pi/2) * tensor2.second;

    double theta1 = std::arg(diag.diagonal()[0]);
    double theta2 = std::arg(diag.diagonal()[1]);
    double theta3 = std::arg(diag.diagonal()[2]);

    double rz_theta1 = (theta1 + theta2);
    double rz_theta2 = (theta1 + theta3);
    double rz_theta3 = (theta2 + theta3);

    double phase = -std::arg(exp_phase) + extracted_phase;

    return DecompTwoQubitsResult{A, B, C, D, rz_theta1, rz_theta2, rz_theta3, phase};
}