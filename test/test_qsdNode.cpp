//
// Created by Adenilton on 15/02/26.
//
#include<iostream>
#include <gtest/gtest.h>
#include "../include/csdNode.h"
#include "../include/nodeVisitor.h"
#include "../include/qsd.h"
#include "../include/randomUnitary.h"
#include "unsupported/Eigen/KroneckerProduct"

TEST(QSDNODE_TEST, RND4) {
    Eigen::MatrixXcf A = Eigen::MatrixXcf::Zero(4,4);
    A(0, 0) = 1.0;
    A(1, 1) = 1.0;
    A(3, 2) = 1.0;
    A(2, 3) = 1.0;

    std::cout << A << std::endl;

    Eigen::MatrixXcf Id = Eigen::MatrixXcf::Identity(2,2);
    Eigen::MatrixXcf X = Eigen::MatrixXcf::Zero(2,2);
    X(1,0) = 1.0;
    X(0, 1) = 1.0;

    auto result = qsd(Id, X);

    const auto V = KroneckerProduct(Id, result.left_gate);
    const auto W = KroneckerProduct(Id, result.right_gate);

    std::cout << "left_gate\n" << result.left_gate << std::endl;
    std::cout << "right_gate\n" << result.right_gate << std::endl;

    std::cout << "list_d\n" << result.list_d << std::endl;


    MatrixXcf d_matrix = MatrixXcf::Zero(2*result.list_d.size(), 2*result.list_d.size());
    const MatrixXcf D = result.list_d.asDiagonal();
    const MatrixXcf Zero = MatrixXcf::Zero(2, 2);
    d_matrix << D, Zero, Zero, D.adjoint();

    const auto result_matrix = V * d_matrix * W;

    std::cout << "result: \n" <<  result_matrix << std::endl;

    auto tree = qsdNode(Id, X);
    auto visitor = qasmVisitor(2);
    tree.accept(visitor);

    std::cout << visitor.qasm_code;



}