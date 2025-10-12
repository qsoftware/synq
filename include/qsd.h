//
// Created by adenilton on 10/11/25.
//

#ifndef TEST_QSD_H
#define TEST_QSD_H

#include "../csd/csd.h"

class qsd {
public:
    explicit qsd(Eigen::MatrixXcf unitary);
    csd_result out;
};

#endif //TEST_QSD_H