//
// Created by Adenilton on 03/02/26.
//

#include<iostream>
#include <gtest/gtest.h>
#include "../include/csd.h"
#include "../include/randomUnitary.h"
#include "../include/one_qubit_gate.h"
#include <unsupported/Eigen/KroneckerProduct>


using namespace std;

TEST(CSD_TEST, Toffoli) {
    Eigen::MatrixXcf T(4, 4);
    T << 1.0f + 0.0if, 0.0f + 0.0if,  0.0f + 0.0if,  0.0f + 0.0if,
         0.0f + 0.0if, 1.0f + 0.0if,  0.0f + 0.0if,  0.0f + 0.0if,
         0.0f + 0.0if, 0.0f + 0.0if,  0.0f + 0.0if,  1.0f + 0.0if,
         0.0f + 0.0if, 0.0f + 0.0if,  1.0f + 0.0if,  0.0f + 0.0if;

    auto out = csd(T, 2, 2);
    auto test_ver = verify(T, out);

    ASSERT_TRUE(test_ver);
}

TEST(CSD_TEST, UNITARY_GATE) {
    Eigen::MatrixXcf T(4, 4);
    T << (2.364882618857422347e-01f-4.822143073226670129e-01if),
    (6.829647732695640072e-02f+4.447897757628520299e-02if),
    (5.115007388835463020e-01f-9.529865837327188005e-02if),
    (5.520571865947265788e-01f-3.597469917942807127e-01if),
 (-1.114816726427168281e-01f+6.978032089401685623e-01if),
    (-2.494382786327934143e-01f+8.410701104345530421e-02if),
    (4.670592911721889728e-01f-4.319298125262461818e-01if),
    (1.622065409340211206e-01f-1.817602163908954571e-02if),
 (-2.235192632554647840e-01f-2.591835234472851801e-01if),
    (2.436499611612899796e-01f-5.375622959989125038e-02if),
    (-9.259829400673319499e-02f-5.381058778231553319e-01if),
    (2.726770843179722736e-01f+6.694197728029507299e-01if),
 (1.151995301476407552e-01f-2.859672891120846883e-01if),
    (-9.016492271986602525e-01f-2.209918615739195902e-01if),
    (-9.627758278964056171e-02f-1.310559524249576557e-01if),
    (-2.001520487551047139e-02f+1.276612783179109645e-01if);

    auto out = csd(T, 2, 2);
    auto test_ver = verify(T, out);

    ASSERT_TRUE(test_ver);
}

TEST(CSD_TEST, SEP) {
    auto matrix1 = OneQubit::h_matrix();
    auto matrix2 = OneQubit::x_matrix();
    Eigen::MatrixXcf T = Eigen::KroneckerProduct<Eigen::MatrixXcf, Eigen::MatrixXcf>(matrix1, matrix2);

    auto out = csd(T, 2, 2);
    auto test_ver = verify(T, out);

    ASSERT_TRUE(test_ver);
}

TEST(CSD_TEST, RND4) {
    auto T =random_unitary_matrix(4);

    auto out = csd(T, 2, 2);
    auto test_ver = verify(T, out);

    ASSERT_TRUE(test_ver);
}

TEST(CSD_TEST, RND8) {
    auto T =random_unitary_matrix(8);

    auto out = csd(T, 4, 4);
    auto test_ver = verify(T, out);

    ASSERT_TRUE(test_ver);
}

TEST(CSD_TEST, RND16) {
    auto T =random_unitary_matrix(16);

    auto out = csd(T, 8, 8);
    auto test_ver = verify(T, out);

    ASSERT_TRUE(test_ver);
}
