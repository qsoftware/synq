//
// Created by Alex on 06/05/26.
//

#ifndef TWO_QUBITS_UNITARY_H
#define TWO_QUBITS_UNITARY_H

#include <Eigen/Dense>
#include <complex>
#include <tuple>
#include <utility>

using twoQubitsTensor = std::pair<Eigen::Matrix2cd, Eigen::Matrix2cd>;

struct DecompTwoQubitsResult {
    Eigen::Matrix2cd A, B, C, D;
    double angleRz, angleRy1, angleRy2, phase;
};

DecompTwoQubitsResult twoQubitsUnitaryDecomposition(const Eigen::Matrix4cd& U);

std::tuple<twoQubitsTensor, Eigen::Matrix4cd, twoQubitsTensor> su4Decomposition(const Eigen::Matrix4cd& specialUnitary);

twoQubitsTensor decomposeKronecker(const Eigen::Matrix4cd& M);

/**
 * Two-CNOT circuit, applied in time order:
 *   C (q1), D (q0); CX q0->q1; Rz(angleRz) (q1), Ry(angleRy) (q0); CX q0->q1; A (q1), B (q0);
 * times e^{i phase}. `error` is the max entry error of the reconstruction.
 */
struct DecompTwoCxResult {
    Eigen::Matrix2cd A, B, C, D;
    double angleRy, angleRz, phase, error;
};

/** Diagonal D = exp(i theta ZZ) such that D^dagger U can be implemented with two CNOTs. */
Eigen::Matrix4cd diagonalForTwoCx(const Eigen::Matrix4cd& U);

/** Two-CNOT decomposition of V; valid when V needs at most two CNOTs (check `error`). */
DecompTwoCxResult twoCxDecomposition(const Eigen::Matrix4cd& V);

#endif 