//
// Created by adenilton on 10/3/25.
//

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "../include/one_qubit_gate.h"
#include <complex>

zyz_result OneQubit::zyz_decomposition(Eigen::Matrix2cd uMatrix) {

    double const x00 = arg(uMatrix(0, 0));
    double const x01 = arg(uMatrix(0, 1)) + M_PI;
    double const x10 = arg(uMatrix(1, 0));

    Eigen::Matrix3d A;
    A << 1.0, -1.0/2, -1.0/2,
         1.0, -1.0/2,  1.0/2,
         1.0,  1.0/2, -1.0/2;

    Eigen::Vector3d b;
    b << x00, x01, x10;

    Eigen::Vector3d angles = A.lu().solve(b);
    double const gamma =  acos(abs(uMatrix(0, 0))) * 2;
    double const alpha =  angles(0);
    double const beta =  angles(1);
    double const delta =  angles(2);
    return {alpha, beta, gamma, delta};

}

Eigen::Matrix2cd OneQubit::h_matrix() {
    Eigen::Matrix2cd const h =
        (Eigen::Matrix2cd() <<  1.0/sqrt(2.0)+0.0i, 1.0/sqrt(2.0)+0.0i,
                                1.0/sqrt(2.0)+0.0i, -1.0/sqrt(2.0)+0.0i).finished();
    return h;
}

Eigen::Matrix2cd OneQubit::x_matrix() {
    Eigen::Matrix2cd const x =
        (Eigen::Matrix2cd() << 0.0 + 0.0i, 1.0 + 0.0i,
                                  1.0 + 0.0i, 0.0 + 0.0i).finished();
    return x;
}

Eigen::Matrix2cd OneQubit::ry_matrix(double theta) {
    Eigen::Matrix2cd const ry =
        (Eigen::Matrix2cd() << cos(theta/2)+0.i, -sin(theta/2)+0.i,
                               sin(theta/2)+0.i, cos(theta/2)+0.i).finished();
    return ry;
}

Eigen::Matrix2cd OneQubit::rz_matrix(double theta) {
    Eigen::Matrix2cd const rz =
        (Eigen::Matrix2cd() << exp(-theta*1.i/2.0), 0.i,
                               0.i,              exp(theta*1.i/2.0)).finished();
    return rz;
}

std::complex<double> OneQubit::gphase(double theta) {
    return std::exp(std::complex<double>(0.0, theta));
}

Eigen::Matrix2cd OneQubit::p_matrix(double theta) {
    Eigen::Matrix2cd const p =
        (Eigen::Matrix2cd() << 1.0, 0.i,
                               0.i, OneQubit::gphase(theta)).finished();
    return p;
}
