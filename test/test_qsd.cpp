//
// Created by adenilton on 10/11/25.
//
#include <gtest/gtest.h>
#include <Eigen/Dense>
#include "../csd/csd.h"
#include<complex.h>


TEST(HelloTest, qsd) {
    Eigen::MatrixXcf gate(2, 2);
    gate << 0.0f + 0.0fi,  1.0f + 0.0fi,
         1.0f + 0.0fi, 0.0f + 0.0fi;

    auto result = csd(gate, 1, 1);
    std::cout << result.theta[0] << std::endl;
}
