//
// Created by Alex on 06/05/26.
//

#include <complex>
#include "../include/twoQubitsUnitary.h"
#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include <numbers>
#include <limits>
#include <algorithm>
#include <cmath>

using namespace std;

using twoQubitsTensor = std::pair<Eigen::Matrix2cd, Eigen::Matrix2cd>;

static const complex<double> i(0,1);
static const double pi = std::numbers::pi;

static const Eigen::Matrix4cd E = (Eigen::Matrix4cd() << 
    1,  i,  0,  0,
    0,  0,  i,  1,
    0,  0,  i, -1,
    1, -i,  0,  0
).finished() * (1.0 / std::sqrt(2.0));

static const Eigen::Matrix4cd SWAP = (Eigen::Matrix4cd() << 
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

// Real orthogonal P (det +1) such that P^T M P is diagonal, for M complex
// symmetric and unitary.
static Eigen::Matrix4d realEigenbasis(const Eigen::Matrix4cd& M) {
    // M is complex symmetric and unitary, so Re(M) and Im(M) are commuting real
    // symmetric matrices and share a real orthonormal eigenbasis. A generic real
    // combination x*Re(M) + y*Im(M) has that basis as its eigenbasis, so a real
    // symmetric eigensolver gives op1 directly, also when eigenvalues of M are
    // degenerate (the previous kernel-of-Hk approach depended on a 1e-10 grouping
    // threshold, lost precision and could leave columns of op1 uninitialized).
    // A few fixed coefficient pairs are tried in case one is accidentally degenerate.
    const Eigen::Matrix4d A = M.real();
    const Eigen::Matrix4d B = M.imag();
    static const double coefficients[][2] = {
        {1.0, 0.0}, {0.7853981633974483, 0.6180339887498949}, {0.3141592653589793, 1.4142135623730951},
        {1.7320508075688772, -0.5772156649015329}, {-0.4142135623730951, 1.2599210498948732},
        {0.9189385332046727, 2.718281828459045}};

    Eigen::Matrix4d real_op1 = Eigen::Matrix4d::Identity();
    double best_error = std::numeric_limits<double>::infinity();
    for (const auto& coef : coefficients) {
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix4d> solver(coef[0] * A + coef[1] * B);
        const Eigen::Matrix4d P = solver.eigenvectors();
        Eigen::Matrix4cd D = P.transpose().cast<std::complex<double>>() * M * P.cast<std::complex<double>>();
        const double error = (D - Eigen::Matrix4cd(D.diagonal().asDiagonal())).cwiseAbs().maxCoeff();
        if (error < best_error) {
            best_error = error;
            real_op1 = P;
        }
        if (error < 1e-13) break;
    }

    if (real_op1.determinant() < 0.0) {
        real_op1.col(0) *= -1.0;
    }
    return real_op1;
}

std::pair<Eigen::Matrix4cd, Eigen::Matrix4cd> diagonalizes(const Eigen::Matrix4cd& SUm){

    /** Diagonalizes M = SUm * SUm^T while forcing the eigenvectors (op1) into SO(4, R).
     * A unitary in SU(4) can be decomposed in op1 * diag * op2, with op1 and op2 into SO(4, R).
     * based on "Nonlocal properties of two-qubit gates and mixed states and optimization of quantum computations"
     * by Makhlin (2000).
     **/

    Eigen::Matrix4cd M = SUm * SUm.transpose();

    Eigen::Matrix4d real_op1 = realEigenbasis(M);

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

// ---------------------------------------------------------------------------
// Optimization A.2 (Shende, Bullock, Markov, "Synthesis of quantum logic
// circuits", 2006): a two-qubit unitary U is written as U = D * V, where D is
// diagonal and V needs only two CNOTs. D is then merged into the next
// two-qubit block of the decomposition.
// ---------------------------------------------------------------------------

namespace {

const Eigen::Matrix4cd YY = (Eigen::Matrix4cd() <<
     0, 0, 0, -1,
     0, 0, 1,  0,
     0, 1, 0,  0,
    -1, 0, 0,  0).finished();

// Eigenvalue signs of XX, ZZ and YY in the magic basis E (columns of E).
const double X_SIGNS[4] = {1.0, -1.0, 1.0, -1.0};
const double Z_SIGNS[4] = {1.0, 1.0, -1.0, -1.0};
const double Y_SIGNS[4] = {-1.0, 1.0, 1.0, -1.0};

// CNOT with control q0 and target q1 (little-endian order, as in Qiskit).
const Eigen::Matrix4cd CX01 = (Eigen::Matrix4cd() <<
    1, 0, 0, 0,
    0, 0, 0, 1,
    0, 0, 1, 0,
    0, 1, 0, 0).finished();

Eigen::Matrix2cd ry2(double theta) {
    Eigen::Matrix2cd m;
    m << std::cos(theta / 2), -std::sin(theta / 2),
         std::sin(theta / 2),  std::cos(theta / 2);
    return m;
}

Eigen::Matrix4cd kron2(const Eigen::Matrix2cd& q1, const Eigen::Matrix2cd& q0) {
    Eigen::Matrix4cd m;
    for (int r = 0; r < 2; ++r)
        for (int c = 0; c < 2; ++c)
            m.block<2, 2>(2 * r, 2 * c) = q1(r, c) * q0;
    return m;
}

std::complex<double> fourthRootOfDet(const Eigen::Matrix4cd& U) {
    return std::pow(U.determinant(), 0.25);
}

}  // namespace

Eigen::Matrix4cd diagonalForTwoCx(const Eigen::Matrix4cd& U) {
    /*
     * V = exp(-i theta ZZ) U needs two CNOTs iff tr(gamma(V)) is real, where
     * gamma(V) = V (Y x Y) V^T (Y x Y) and V is normalized to SU(4).
     * ZZ commutes with YY, so tr(gamma(V)) = cos(2 theta) tr(G) - i sin(2 theta) tr(ZZ G),
     * with G = gamma(U); its imaginary part vanishes for the theta below.
     * Returns D = exp(i theta ZZ), so that U = D V.
     */
    const Eigen::Matrix4cd Us = U / fourthRootOfDet(U);
    const Eigen::Matrix4cd G = Us * YY * Us.transpose() * YY;
    const std::complex<double> t0 = G.trace();
    std::complex<double> t1 = 0.0;
    for (int k = 0; k < 4; ++k) {
        const double zz = (k == 0 || k == 3) ? 1.0 : -1.0;
        t1 += zz * G(k, k);
    }
    const double theta = 0.5 * std::atan2(t0.imag(), t1.real());
    Eigen::Matrix4cd D = Eigen::Matrix4cd::Zero();
    D(0, 0) = std::exp(i * theta);
    D(1, 1) = std::exp(-i * theta);
    D(2, 2) = std::exp(-i * theta);
    D(3, 3) = std::exp(i * theta);
    return D;
}

DecompTwoCxResult twoCxDecomposition(const Eigen::Matrix4cd& V) {
    /*
     * Canonical (KAK) decomposition with one interaction coefficient equal to zero:
     *   V = e^{i phase} (A x B) exp(i (a XX + c ZZ)) (C x D),
     * and exp(i (a XX + c ZZ)) = CX01 (Rz(-2c) x Rx(-2a)) CX01, where
     * Rx(t) = Rz(-pi/2) Ry(t) Rz(pi/2) and Rz on the control commutes with CX01.
     *
     * In the magic basis V_B = O1 Delta O2 with O1, O2 in SO(4). The phases of Delta
     * (defined modulo pi through Delta^2) are matched to a*x_k + c*z_k + g by trying
     * the orderings and branches; the YY component must vanish and the phases must
     * sum to 0 mod 2pi so that O2 has determinant +1.
     */
    const Eigen::Matrix4cd Vs = V / fourthRootOfDet(V);
    const Eigen::Matrix4cd VB = E.adjoint() * Vs * E;
    const Eigen::Matrix4cd M = VB * VB.transpose();
    const Eigen::Matrix4d O1 = realEigenbasis(M);
    const Eigen::Matrix4cd Md = O1.transpose().cast<std::complex<double>>() * M * O1.cast<std::complex<double>>();

    double half[4];
    for (int k = 0; k < 4; ++k) half[k] = std::arg(Md(k, k)) / 2.0;

    // With the phases ordered as the X/Z sign patterns, lambda = (g+a+c, g-a+c, g+a-c, g-a-c):
    // the YY component vanishes iff lambda_0 + lambda_3 == lambda_1 + lambda_2. So the
    // four phases must split into two pairs with equal sums. Each half[k] is only known
    // modulo pi, so try the 3 pairings and the 16 branch choices of each (48 cases):
    // the pair sums must be exactly equal and the total must be 0 mod 2pi (so that O2
    // has determinant +1). All pairings are tried because, with repeated eigenvalues
    // (SWAP, CNOT, CZ, ...), several pairings can look equally good modulo pi while
    // only some of them have a valid choice of branches.
    static const int pairings[3][4] = {{0, 2, 3, 1}, {0, 1, 3, 2}, {0, 1, 2, 3}};  // order: p, r, s, q
    int best_perm[4] = {0, 1, 2, 3};
    double best_lambda[4] = {half[0], half[1], half[2], half[3]};
    double best_residual = std::numeric_limits<double>::infinity();
    for (const auto& order : pairings) {
        for (int t = 0; t < 16; ++t) {
            double lambda[4];
            double total = 0.0;
            for (int k = 0; k < 4; ++k) {
                lambda[k] = half[order[k]] + (((t >> k) & 1) ? pi : 0.0);
                total += lambda[k];
            }
            if (std::abs(std::remainder(total, 2 * pi)) > 1e-9) continue;
            double y = 0.0;
            for (int k = 0; k < 4; ++k) y += Y_SIGNS[k] * lambda[k];
            const double residual = std::abs(y);
            if (residual < best_residual) {
                best_residual = residual;
                std::copy(order, order + 4, best_perm);
                std::copy(lambda, lambda + 4, best_lambda);
            }
        }
    }

    double a = 0.0, c = 0.0, g = 0.0;
    for (int k = 0; k < 4; ++k) {
        a += X_SIGNS[k] * best_lambda[k] / 4.0;
        c += Z_SIGNS[k] * best_lambda[k] / 4.0;
        g += best_lambda[k] / 4.0;
    }

    Eigen::Matrix4d P = Eigen::Matrix4d::Zero();
    for (int k = 0; k < 4; ++k) P(k, best_perm[k]) = 1.0;
    Eigen::Matrix4d O1p = O1 * P.transpose();
    if (O1p.determinant() < 0.0) O1p.col(0) *= -1.0;

    Eigen::Vector4cd delta;
    for (int k = 0; k < 4; ++k) delta(k) = std::exp(i * (a * X_SIGNS[k] + c * Z_SIGNS[k] + g));
    const Eigen::Matrix4cd O2 = delta.conjugate().asDiagonal() * O1p.transpose().cast<std::complex<double>>() * VB;

    const Eigen::Matrix4cd K1 = E * O1p.cast<std::complex<double>>() * E.adjoint();
    const Eigen::Matrix4cd K2 = E * O2.real().cast<std::complex<double>>() * E.adjoint();
    const auto [A, B] = decomposeKronecker(K1);
    const auto [C, D] = decomposeKronecker(K2);

    DecompTwoCxResult result;
    result.A = A;                           // qubit 1, applied last
    result.B = B * Rz(-pi / 2);             // qubit 0, applied last
    result.C = C;                           // qubit 1, applied first
    result.D = Rz(pi / 2) * D;              // qubit 0, applied first
    result.angleRy = -2.0 * a;              // qubit 0, between the CNOTs
    result.angleRz = -2.0 * c;              // qubit 1, between the CNOTs

    const Eigen::Matrix4cd W = kron2(result.A, result.B) * CX01 *
                               kron2(Rz(result.angleRz), ry2(result.angleRy)) * CX01 *
                               kron2(result.C, result.D);
    const std::complex<double> overlap = (W.conjugate().cwiseProduct(V)).sum();
    result.phase = std::arg(overlap);
    result.error = (std::exp(i * result.phase) * W - V).cwiseAbs().maxCoeff();
    return result;
}
