//
// Created by Adenilton on 08/02/26.
//

#include<iostream>
#include <gtest/gtest.h>
#include "../include/qsd.h"
#include "../include/one_qubit_gate.h"
#include "../include/randomUnitary.h"
#include <unsupported/Eigen/KroneckerProduct>



using namespace std;

TEST(QSD_TEST, h_x) {

    const auto id = MatrixXcf::Identity(2, 2);

    const auto M1 = OneQubit::h_matrix();
    const auto M2 = OneQubit::x_matrix();


    // Expected result
    const MatrixXcf Zero = MatrixXcf::Zero(M1.rows(), M1.cols());
    MatrixXcf expected_result = MatrixXcf::Zero(4, 4);
    expected_result << M1,  Zero, Zero, M2;

    // quantum Shannon decomposition
    auto result = qsd(M1, M2);

    // correct?
    const auto V = KroneckerProduct(id, result.left_gate);
    const auto W = KroneckerProduct(id, result.right_gate);
    MatrixXcf d_matrix = MatrixXcf::Zero(2*result.list_d.size(), 2*result.list_d.size());
    const MatrixXcf D = result.list_d.asDiagonal();

    d_matrix << D, Zero, Zero, D.adjoint();

    const auto result_matrix = V * d_matrix * W;
    
    EXPECT_TRUE(result_matrix.isApprox(expected_result));
}

TEST(QSD_TEST, rnd_rnd) {

    const auto id = Eigen::MatrixXcf::Identity(2, 2);

    const auto M1 = random_unitary_matrix(4);
    const auto M2 = random_unitary_matrix(4);


    // Expected result
    const MatrixXcf Zero = MatrixXcf::Zero(M1.rows(), M1.cols());
    MatrixXcf expected_result = MatrixXcf::Zero(2*M1.rows(), 2*M1.cols());
    expected_result << M1,  Zero, Zero, M2;

    // quantum Shannon decomposition
    auto result = qsd(M1, M2);

    // correct?
    const auto V = KroneckerProduct(id, result.left_gate);
    const auto W = KroneckerProduct(id, result.right_gate);
    MatrixXcf d_matrix = MatrixXcf::Zero(2*result.list_d.size(), 2*result.list_d.size());
    const MatrixXcf D = result.list_d.asDiagonal();

    d_matrix << D, Zero, Zero, D.adjoint();

    const auto result_matrix = V * d_matrix * W;

    EXPECT_TRUE(result_matrix.isApprox(expected_result));
}