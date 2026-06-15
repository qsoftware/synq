#include <gtest/gtest.h>
#include "../include/diagonalGateNode.h"
#include "../include/nodeVisitor.h"
#include <iostream>
#include <complex>
#include <Eigen/Dense>

using namespace std::complex_literals;

TEST(DiagonalGateNodeTests, Diagonal4x4Matrix) {

    Eigen::Matrix4cd diagonal;
    
    Eigen::Vector4d random_phases = Eigen::Vector4d::Random();
    Eigen::Vector4cd diag_elements = (1.0i * random_phases.cast<std::complex<double>>().array()).exp();

    diagonal = diag_elements.asDiagonal();

    auto diagonalNode = diagonalGateNode(diagonal);
    auto visitor = qasmVisitor(2);

    diagonalNode.accept(visitor);

    std::cout << visitor.qasm_code << std::endl;
    std::cout << "Expected diagonal elements: " << diag_elements.transpose() << std::endl;

    EXPECT_FALSE(visitor.qasm_code.empty());
    EXPECT_NE(visitor.qasm_code.find("OPENQASM"), std::string::npos);
}

TEST(DiagonalGateNodeTests, Diagonal8x8Matrix) {

    Eigen::MatrixXcd diagonal;
    
    Eigen::VectorXd random_phases = Eigen::VectorXd::Random(8);
    Eigen::VectorXcd diag_elements = (1.0i * random_phases.cast<std::complex<double>>().array()).exp();

    diagonal = diag_elements.asDiagonal();

    auto diagonalNode = diagonalGateNode(diagonal);
    auto visitor = qasmVisitor(3);

    diagonalNode.accept(visitor);

    std::cout << visitor.qasm_code << std::endl;
    std::cout << "Expected diagonal elements: " << diag_elements.transpose() << std::endl;

    EXPECT_FALSE(visitor.qasm_code.empty());
    EXPECT_NE(visitor.qasm_code.find("OPENQASM"), std::string::npos);
}


TEST(DiagonalGateNodeTests, Diagonal16x16Matrix) {

    Eigen::MatrixXcd diagonal;
    
    Eigen::VectorXd random_phases = Eigen::VectorXd::Random(16);
    Eigen::VectorXcd diag_elements = (1.0i * random_phases.cast<std::complex<double>>().array()).exp();

    diagonal = diag_elements.asDiagonal();

    auto diagonalNode = diagonalGateNode(diagonal);
    auto visitor = qasmVisitor(4);

    diagonalNode.accept(visitor);

    std::cout << visitor.qasm_code << std::endl;
    std::cout << "Expected diagonal elements: " << diag_elements.transpose() << std::endl;

    EXPECT_FALSE(visitor.qasm_code.empty());
    EXPECT_NE(visitor.qasm_code.find("OPENQASM"), std::string::npos);
}