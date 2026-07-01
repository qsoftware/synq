//
// Created by adenilton on 10/3/25.
//

#pragma once
#include<Eigen/Dense>
using namespace std::complex_literals;


struct zyz_result{
    double alpha, beta, gamma, delta;
};

class OneQubit{
public:
    /*
     * zyz decomposition of 2 x 2 unitary matriz
     * @param uMatrix 2 x 2 unitary matrix (quantum gate)
     * @return alpha, beta, gamma and delta where
     * e^(i alpha)Rz(beta)Ry(gamma)Rz(delta) = uMatrix
     */
    static zyz_result zyz_decomposition(Eigen::Matrix2cd uMatrix);

    /*
     * @return Hadamard matrix
     */
    static Eigen::Matrix2cd h_matrix();

    /*
     * @return X matrix
     */
    static Eigen::Matrix2cd x_matrix();

    /*
     * @return Ry(theta) matrix
     * @param theta rotation angle
     */
    static Eigen::Matrix2cd ry_matrix(double theta);

    /*
     * @param theta rotation angle
     * @return Rz(theta) matrix
     */
    static Eigen::Matrix2cd rz_matrix(double theta);

    /**
     * @brief Returns the global phase shift for a given angle
     * 
     * @param theta Angle in radians used to compute the global phase shift
     * @return std::complex<double> The complex number representing the global phase shift
     */
    static std::complex<double> gphase(double theta);

    /**
     * @brief Returns the relative phase shift matrix for a given angle
     * 
     * @param theta Angle in radians used to compute the relative phase shift
     * @return Eigen::Matrix2cd The 2x2 matrix representing the relative phase shift
     */
    static Eigen::Matrix2cd p_matrix(double theta);
};


