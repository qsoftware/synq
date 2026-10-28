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

#endif 