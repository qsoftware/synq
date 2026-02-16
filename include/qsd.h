//
// Created by adenilton on 10/11/25.
//
#pragma once

#include<complex>
#include <Eigen/Dense>

struct qsd_result {
    Eigen::VectorXcf list_d;
    Eigen::MatrixXcf left_gate; // left_gate in matrix notation
    Eigen::MatrixXcf right_gate;
};

qsd_result qsd(Eigen::MatrixXcf gate1, Eigen::MatrixXcf gate2);

