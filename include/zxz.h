/**
 * @file zxz.h
 * @brief Block-ZXZ decomposition of a unitary matrix (Krol and Al-Ars, 2024).
 */

#pragma once

#include <Eigen/Dense>

/**
 * @brief Blocks of the Block-ZXZ decomposition
 *
 *     U = (A1 ⊕ A2) (H ⊗ I) (I ⊕ B) (H ⊗ I) (I ⊕ C),
 *
 * where H acts on the most significant qubit and all blocks are unitary with half
 * the dimension of U.
 */
struct zxz_result {
    Eigen::MatrixXcd a1;
    Eigen::MatrixXcd a2;
    Eigen::MatrixXcd b;
    Eigen::MatrixXcd c;
};

/**
 * @brief Computes the Block-ZXZ decomposition of a unitary matrix of even dimension.
 *
 * With U = [[X, Y], [U21, U22]], the middle factor is ½[[I + B, I - B], [I - B, I + B]], so
 *
 *     X = A1 (I + B) / 2,   Y = A1 (I - B) C / 2,   A1 = X + Y C†,   A2 = U21 + U22 C†.
 *
 * C is chosen so that A1 is unitary: with the SVD X = V S W† and XX† + YY† = I,
 * Y = V K Q† with K = sqrt(I - S²) and Q the unitary polar factor of Y†V. Then
 * C† = i Q W† gives A1 = V (S + iK) W†, which is unitary, and B = 2 A1† X - I.
 */
zxz_result blockZxz(const Eigen::MatrixXcd& u);
