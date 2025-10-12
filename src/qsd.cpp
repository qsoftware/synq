//
// Created by adenilton on 10/11/25.
//

#include "../include/qsd.h"

qsd::qsd(Eigen::MatrixXcf unitary) {
    int p = unitary.cols() / 2;
    auto result = csd(unitary, p, p);
}
