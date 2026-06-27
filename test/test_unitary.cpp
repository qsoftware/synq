#include <gtest/gtest.h>
#include "test_utils.h"
#include "../include/unitaryGateNode.h"
#include "../include/nodeVisitor.h"
#include "../include/randomUnitary.h"
#include "../include/twoQubitsUnitary.h"
#include "../include/one_qubit_gate.h"
#include <iostream>
#include <complex>
#include <unsupported/Eigen/KroneckerProduct>


TEST(UnitaryGateNodeTests, Unitary2x2Matrix) {

    const double inv_sqrt2 = 1.0 / std::sqrt(2.0);

    Eigen::MatrixXcd matrix(2, 2);
    matrix << std::complex<double>(inv_sqrt2, 0.0),  std::complex<double>(0.0, -inv_sqrt2),
                std::complex<double>(0.0, -inv_sqrt2), std::complex<double>(inv_sqrt2, 0.0);

    auto unitaryNode = unitaryGateNode(matrix);
    auto visitor = qasmVisitor(1);

    unitaryNode.accept(visitor);

    std::cout << visitor.qasm_code << std::endl;
    Eigen::Matrix2cd converted_matrix = one_qubit_qasm_to_matrix(visitor.qasm_code);

    EXPECT_FALSE(visitor.qasm_code.empty());
    EXPECT_NE(visitor.qasm_code.find("OPENQASM"), std::string::npos);
    EXPECT_TRUE(converted_matrix.isApprox(matrix, 1e-6));
}

TEST(UnitaryGateNodeTests, IdentityGate) {
    Eigen::MatrixXcd eye(4, 4);
   
    eye << 1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1, 0,
            0, 0, 0, 1;

   
    auto unitaryNode = unitaryGateNode(eye);
    auto visitor = qasmVisitor(2);
    
    unitaryNode.accept(visitor);

    std::cout << visitor.qasm_code << std::endl;
    Eigen::MatrixXcd reconstructed = qasm_to_matrix(visitor.qasm_code, 2);
    
    EXPECT_FALSE(visitor.qasm_code.empty());
    EXPECT_NE(visitor.qasm_code.find("OPENQASM"), std::string::npos);
    EXPECT_TRUE(eye.isApprox(reconstructed, 1e-6));
}

TEST(UnitaryGateNodeTests, CNOTGate) {
    Eigen::MatrixXcd cnot(4, 4);
    cnot << 1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 0, 1,
            0, 0, 1, 0;

   
    auto unitaryNode = unitaryGateNode(cnot);
    auto visitor = qasmVisitor(2);
    
    unitaryNode.accept(visitor);

    std::cout << visitor.qasm_code << std::endl;
    Eigen::MatrixXcd reconstructed = qasm_to_matrix(visitor.qasm_code, 2);
    
    EXPECT_FALSE(visitor.qasm_code.empty());
    EXPECT_NE(visitor.qasm_code.find("OPENQASM"), std::string::npos);
    EXPECT_TRUE(cnot.isApprox(reconstructed, 1e-6));
}

TEST(UnitaryGateNodeTests, Random4x4Unitary) {
    Eigen::MatrixXcd T = random_unitary_matrix(4);
   
    auto unitaryNode = unitaryGateNode(T);
    auto visitor = qasmVisitor(2);
    
    unitaryNode.accept(visitor);

    std::cout << visitor.qasm_code << std::endl;
    std::cout << "Expected Unitary: \n" << std::endl;
    std::cout << T << std::endl;

    Eigen::MatrixXcd reconstructed = qasm_to_matrix(visitor.qasm_code, 2);
    
    EXPECT_FALSE(visitor.qasm_code.empty());
    EXPECT_TRUE(T.isApprox(reconstructed, 1e-6));
    EXPECT_NE(visitor.qasm_code.find("OPENQASM"), std::string::npos);
}

TEST(UnitaryGateNodeTests, Random8x8Unitary) {
    Eigen::MatrixXcd T = random_unitary_matrix(8);
   
    auto unitaryNode = unitaryGateNode(T);
    auto visitor = qasmVisitor(3);
    
    unitaryNode.accept(visitor);

    std::cout << visitor.qasm_code << std::endl;
    std::cout << "Expected Unitary: \n" << std::endl;
    std::cout << T << std::endl;

    Eigen::MatrixXcd reconstructed = qasm_to_matrix(visitor.qasm_code, 3);
    
    EXPECT_FALSE(visitor.qasm_code.empty());
    EXPECT_TRUE(T.isApprox(reconstructed, 1e-6));
    EXPECT_NE(visitor.qasm_code.find("OPENQASM"), std::string::npos);
}

TEST(UnitaryGateNodeTests, SU4Test1) {
    Eigen::Matrix4cd su4Matrix;
    su4Matrix <<
        1, 0, 0, 0,
        0, 0, 1, 0,
        0, 1, 0, 0,
        0, 0, 0, -1;

    auto [tensor1, diag, tensor2] = su4Decomposition(su4Matrix);

    Eigen::Matrix4cd E = (Eigen::Matrix4cd() << 
        1,  std::complex<double>(0,1),  0,  0,
        0,  0,  std::complex<double>(0,1),  1,
        0,  0,  std::complex<double>(0,1), -1,
        1, -std::complex<double>(0,1),  0,  0
    ).finished() * (1.0 / std::sqrt(2.0));

    Eigen::Matrix4cd delta = E * diag * E.adjoint();

    Eigen::Matrix4cd reconstructedMatrix = 
        (Eigen::kroneckerProduct(tensor1.first, tensor1.second)) * delta * 
        (Eigen::kroneckerProduct(tensor2.first, tensor2.second));

    EXPECT_TRUE(su4Matrix.isApprox(reconstructedMatrix, 1e-6)) << "ERROR.";
}

TEST(UnitaryGateNodeTests, SU4Test2) {
    Eigen::Matrix4cd su4Matrix;
    
    su4Matrix << 
        std::complex<double>( 0.5,  0.5), std::complex<double>( 0.0,  0.0), std::complex<double>(-0.5, -0.5), std::complex<double>( 0.0,  0.0),
        std::complex<double>( 0.0,  0.0), std::complex<double>( 0.5, -0.5), std::complex<double>( 0.0,  0.0), std::complex<double>( 0.5, -0.5),
        std::complex<double>( 0.5, -0.5), std::complex<double>( 0.0,  0.0), std::complex<double>( 0.5, -0.5), std::complex<double>( 0.0,  0.0),
        std::complex<double>( 0.0,  0.0), std::complex<double>(-0.5, -0.5), std::complex<double>( 0.0,  0.0), std::complex<double>( 0.5,  0.5);

    auto [tensor1, diag, tensor2] = su4Decomposition(su4Matrix);

    Eigen::Matrix4cd E = (Eigen::Matrix4cd() << 
        1,  std::complex<double>(0,1),  0,  0,
        0,  0,  std::complex<double>(0,1),  1,
        0,  0,  std::complex<double>(0,1), -1,
        1, -std::complex<double>(0,1),  0,  0
    ).finished() * (1.0 / std::sqrt(2.0));

    Eigen::Matrix4cd delta = E * diag * E.adjoint();

    Eigen::Matrix4cd reconstructedMatrix = 
        (Eigen::kroneckerProduct(tensor1.first, tensor1.second)) * delta * (Eigen::kroneckerProduct(tensor2.first, tensor2.second));

    EXPECT_TRUE(su4Matrix.isApprox(reconstructedMatrix, 1e-6)) << "ERROR.";
}