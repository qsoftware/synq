#pragma once

#ifdef USE_MKL
    #define EIGEN_USE_MKL_ALL
    #include<mkl.h>
    #include<Eigen/Dense>
    #include <complex>
    typedef MKL_Complex16 complexD;
#else
    #define EIGEN_USE_LAPACKE
    #include <complex>
    #define lapack_complex_double std::complex<double>
    #include<Eigen/Dense>
    #include <lapacke.h>
    typedef std::complex<double> complexD;
#endif






struct csd_result {
    Eigen::MatrixXcd U1;
    Eigen::MatrixXcd U2;
    std::vector<double> theta;
    Eigen::MatrixXcd V1T;
    Eigen::MatrixXcd V2T;
};
csd_result csd( const Eigen::MatrixXcd& unitary, Eigen::Index p, Eigen::Index q);
lapack_complex_double* eigen2lapack(const Eigen::MatrixXcd& unitary, int start_line, int start_col, int n_rows, int n_cols);

bool verify(const Eigen::MatrixXcd& unitary, csd_result result);

