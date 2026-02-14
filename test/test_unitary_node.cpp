#include <gtest/gtest.h>
#include "../include/unitaryNode.h"
#include "../include/one_qubit_gate.h"
#include "../include/nodeVisitor.h"
#include <iostream>
#include <complex>

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