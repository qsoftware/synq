//
// Created by Adenilton on 03/02/26.
//

#include<iostream>
#include <gtest/gtest.h>
#include "../include/csd.h"


using namespace std;

TEST(CSD_TEST, Toffoli) {
    Eigen::MatrixXcf T(4, 4);
    T << 1.0f + 0.0if, 0.0f + 0.0if,  0.0f + 0.0if,  0.0f + 0.0if,
         0.0f + 0.0if, 1.0f + 0.0if,  0.0f + 0.0if,  0.0f + 0.0if,
         0.0f + 0.0if, 0.0f + 0.0if,  0.0f + 0.0if,  1.0f + 0.0if,
         0.0f + 0.0if, 0.0f + 0.0if,  1.0f + 0.0if,  0.0f + 0.0if;

    auto out = csd(T, 2, 2);
    auto test_ver = verify(T, out);

    if (test_ver) std::cout << "T gate ok" << std::endl;
    else std::cout << "no T gate ok" << std::endl;

    ASSERT_TRUE(test_ver);
}
