//
// Created by adenilton on 10/11/25.
//
#include <gtest/gtest.h>
#include <Eigen/Dense>
#include "../csd/csd.h"
#include "../include/csdNode.h"
#include<complex.h>


TEST(HelloTest, qsd) {
    Eigen::MatrixXcf gate(4, 4);
    gate << 1.0f + 0.0fi,  0.0f + 0.0fi, 0.0f + 0.0fi,  0.0f + 0.0fi,
            0.0f + 0.0fi,  1.0f + 0.0fi, 0.0f + 0.0fi,  0.0f + 0.0fi,
            0.0f + 0.0fi,  0.0f + 0.0fi, 1.0f + 0.0fi,  0.0f + 0.0fi,
            0.0f + 0.0fi,  0.0f + 0.0fi, 0.0f + 0.0fi,  1.0f + 0.0fi;

    auto result = csdNode(gate);

}
