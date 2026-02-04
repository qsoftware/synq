//
// Created by adenilton on 10/10/25.
//


#pragma once
#define EIGEN_USE_LAPACKE

#include<complex.h>
#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>

#include<Eigen/Dense>
#include <lapacke.h>


struct csd_result {
    Eigen::MatrixXcf U1;
    Eigen::MatrixXcf U2;
    std::vector<double> theta;
    Eigen::MatrixXcf V1T;
    Eigen::MatrixXcf V2T;
};
csd_result csd(Eigen::MatrixXcf& unitary, int p, int q);
lapack_complex_float* eigen2lapack(Eigen::MatrixXcf& unitary, int start_line, int start_col, int n_rows, int n_cols);

bool verify(Eigen::MatrixXcf& unitary, csd_result result);

