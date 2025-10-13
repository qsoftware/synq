//
// Created by adenilton on 10/11/25.
//

#ifndef TEST_QSD_H
#define TEST_QSD_H

#include "../include/IASTnode.h"
#include <Eigen/Dense>

class qsd {
public:
    explicit qsd(Eigen::MatrixXcf first_unitary, Eigen::MatrixXcf second_unitary);
    IASTnode* first_unitary;
    IASTnode* mcrz;
    IASTnode* second_unitary;

};

#endif //TEST_QSD_H