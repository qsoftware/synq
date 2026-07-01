#include<iostream>
using namespace std;
#include <gtest/gtest.h>
#include <vector>
#include <Eigen/Dense>
#include "../include/twoCtrlOperator.h"
#include "../include/one_qubit_gate.h"
#include "../include/nodeVisitor.h"
#include "test_utils.h"

TEST(CtrlOperators, oneCtrlOperator){
    Eigen::Matrix2cd X = OneQubit::x_matrix();
    vector<int> controls = {1}; // Control bit at position 1
    int target = 0;
    Eigen::MatrixXcd Identity = Eigen::MatrixXcd::Identity(4, 4);
    
    // Creating an instance of CtrlOperatorNode
    CtrlOperatorNode operation(controls, target, Identity, X);

    // Applies CNOT gate using the oneCtrlOperator method
    Eigen::MatrixXcd result = operation.oneCtrlOperator(X, controls[0], target);
    
    // The expected matrix is the identity matrix with the last two states swapped (|10> and |11>)
    Eigen::MatrixXcd expected = Eigen::MatrixXcd::Identity(4, 4);
    // Swap the rows and columns corresponding to states 2 (|10>) and 3 (|11>)
    expected(2, 2) = 0.0; expected(2, 3) = 1.0;
    expected(3, 2) = 1.0; expected(3, 3) = 0.0;

    printQuantumMatrix(result, "Result");
    printQuantumMatrix(expected, "Expected");
    
    // Verifies if the obtained matrix is numerically equal to the expected matrix (tolerance of 1e-6)
    EXPECT_TRUE(result.isApprox(expected, 1e-6));
}

TEST(CtrlOperators, twoCtrlOperator){
    Eigen::Matrix2cd X = OneQubit::x_matrix();
    vector<int> controls = {2, 1}; // Control bits at positions 1 and 2
    int target = 0;
    Eigen::MatrixXcd Identity = Eigen::MatrixXcd::Identity(8, 8);
    
    // Creating an instance of CtrlOperatorNode
    CtrlOperatorNode operation(controls, target, Identity, X);

    // Applies CNOT gate using the twoCtrlOperator method
    Eigen::MatrixXcd result = operation.twoCtrlOperator(X, controls[0], controls[1]);
    
    // The expected matrix is the identity matrix with the last two states swapped (|110> and |111>)
    Eigen::MatrixXcd expected = Eigen::MatrixXcd::Identity(8, 8);
    // Swap the rows and columns corresponding to states 6 (|110>) and 7 (|111>)
    expected(6, 6) = 0.0; expected(6, 7) = 1.0;
    expected(7, 6) = 1.0; expected(7, 7) = 0.0;

    printQuantumMatrix(result, "Result");
    printQuantumMatrix(expected, "Expected");

    // Verifies if the obtained matrix is numerically equal to the expected matrix (tolerance of 1e-6)
    EXPECT_TRUE(result.isApprox(expected, 1e-6));
}

TEST(CtrlOperators, oneCtrlOperatorVisitor){
    Eigen::Matrix2cd X = OneQubit::x_matrix();
    vector<int> controls = {1}; // Control bit at position 1
    int target = 0;
    Eigen::MatrixXcd Identity = Eigen::MatrixXcd::Identity(4, 4);
    
    // Creating an instance of CtrlOperatorNode
    CtrlOperatorNode operation(controls, target, Identity, X);

    // Creating an instance of qasmVisitor with 2 qubits
    auto visitor = qasmVisitor(2);
    
    // Accept the visitor
    operation.accept(visitor);

    std::cout << "Generated QASM code:\n" << visitor.qasm_code << "\n";
    Eigen::MatrixXcd reconstructed = qasm_to_matrix(visitor.qasm_code, 2);
    printQuantumMatrix(reconstructed, "Reconstructed");

    // The expected matrix is the identity matrix with the last two states swapped (|10> and |11>)
    Eigen::MatrixXcd expected = Eigen::MatrixXcd::Identity(4, 4);
   // Swap the rows and columns corresponding to states 2 (|10>) and 3 (|11>)
    expected(2, 2) = 0.0; expected(2, 3) = 1.0;
    expected(3, 2) = 1.0; expected(3, 3) = 0.0;

    // Verifies if the generated QASM code contains the expected CNOT instruction
    EXPECT_TRUE(reconstructed.isApprox(expected, 1e-6));
}

TEST(CtrlOperators, twoCtrlOperatorVisitor){
    Eigen::Matrix2cd X = OneQubit::x_matrix();
    vector<int> controls = {2, 1}; // Control bits at positions 1 and 2
    int target = 0;
    Eigen::MatrixXcd Identity = Eigen::MatrixXcd::Identity(8, 8);
    
    // Creating an instance of CtrlOperatorNode
    CtrlOperatorNode operation(controls, target, Identity, X);

    // Creating an instance of qasmVisitor with 3 qubits
    auto visitor = qasmVisitor(3);
    
    // Accept the visitor
    operation.accept(visitor);

    std::cout << "Generated QASM code:\n" << visitor.qasm_code << "\n";
    Eigen::MatrixXcd reconstructed = qasm_to_matrix(visitor.qasm_code, 3);
    printQuantumMatrix(reconstructed, "Reconstructed");

    // The expected matrix is the identity matrix with the last two states swapped (|110> and |111>)
    Eigen::MatrixXcd expected = Eigen::MatrixXcd::Identity(8, 8);
    // Swap the rows and columns corresponding to states 6 (|110>) and 7 (|111>)
    expected(6, 6) = 0.0; expected(6, 7) = 1.0;
    expected(7, 6) = 1.0; expected(7, 7) = 0.0;
    
    // Verifies if the generated QASM code contains the expected CNOT instruction
    EXPECT_TRUE(reconstructed.isApprox(expected, 1e-6));
}