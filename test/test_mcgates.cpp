//
// Created by adenilton on 10/4/25.
//

#include<iostream>
#include "../include/ucrNode.h"
#include "../include/nodeVisitor.h"
using namespace std;
#include <gtest/gtest.h>


TEST(MultiControlledRotations, OneControlUCRZ) {
    // Expect two strings not to be equal.
    auto visitor = qasmVisitor(2);
    vector<double> angles;
    angles.push_back(0.3);
    angles.push_back(0.2);
    auto a = firstUcrNode(&angles, ucrType::Z);
    a.accept(visitor);
    cout << visitor.qasm_code;

    EXPECT_STRNE("hello", "world");
    // Expect equality.
    EXPECT_EQ(7 * 6, 42);
}

TEST(MultiControlledRotations, TwoControlUCRZ) {
    // Expect two strings not to be equal.
    auto visitor = qasmVisitor(3);
    vector angles = {0.2, 0.1, 0.3, 0.4};

    auto a = firstUcrNode(&angles, ucrType::Z);
    a.accept(visitor);
    cout << visitor.qasm_code;

    EXPECT_STRNE("hello", "world");
    // Expect equality.
    EXPECT_EQ(7 * 6, 42);
}

TEST(MultiControlledRotations, ThreeControlUCRZ) {
    // Expect two strings not to be equal.
    ucrNode::index = 0;
    auto visitor = qasmVisitor(4);
    vector angles = {0.11, 0.22, 0.33, 0.44, 0.51, 0.678, 0.712, 0.821};

    auto a = firstUcrNode(&angles, ucrType::Z);
    a.accept(visitor);
    cout << visitor.qasm_code;

    EXPECT_STRNE("hello", "world");
    // Expect equality.
    EXPECT_EQ(7 * 6, 42);
}