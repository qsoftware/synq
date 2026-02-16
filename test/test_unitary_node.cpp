#include <gtest/gtest.h>
#include "../include/unitaryNode.h"
#include "../include/csd.h"
#include "../include/one_qubit_gate.h"
#include "../include/nodeVisitor.h"
#include "../include/randomUnitary.h"
#include <iostream>
#include <complex>
#include <unsupported/Eigen/KroneckerProduct>

TEST(UNITARY_NODE, SingleQubitGate) {
//TODO: fix the tests

    auto matrix = OneQubit::h_matrix();

    auto unitary = unitaryNode(matrix);
    auto visitor = qasmVisitor(1);

    unitary.accept(visitor);

    std::cout << visitor.qasm_code << std::endl;

    EXPECT_FALSE(visitor.qasm_code.empty());
    EXPECT_NE(visitor.qasm_code.find("OPENQASM"), std::string::npos);
}

TEST(UNITARY_NODE, TwoQubitGate) {
    //TODO: fix the tests

    auto matrix1 = OneQubit::h_matrix();
    auto matrix2 = OneQubit::x_matrix();
    Eigen::MatrixXcf matrix = Eigen::KroneckerProduct<Eigen::MatrixXcf, Eigen::MatrixXcf>(matrix1, matrix2);
    std::cout << "matrix\n" << matrix << std::endl;

    std::cout << "matrix1\n" << matrix1 << std::endl;

    auto result = csd(matrix, 2, 2);

    std::cout << "U1\n" << result.U1 << std::endl;
    std::cout << "U2\n" << result.U2 << std::endl;
    std::cout << "V1T\n" << result.V1T << std::endl;
    std::cout << "V2T\n" << result.V2T << std::endl;
    std::cout << "theta\n" << result.theta[0] <<", " << result.theta[1] << std::endl;


    // auto matrix = random_unitary_matrix(4);
    auto unitary = unitaryNode(matrix);
    auto visitor = qasmVisitor(2);

    unitary.accept(visitor);

    std::cout << visitor.qasm_code << std::endl;
}