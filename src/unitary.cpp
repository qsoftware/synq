//
// Created by adenilton on 10/4/25.
//

#include "../include/unitary.h"

unitary::unitary(const Eigen::MatrixXcf &uMatrix) {
    data = uMatrix;
}

std::string unitary::unitary2qasm() {
    std::string qasm = "";
    if (data.cols() == 2) {
        auto results = OneQubit::zyz_decomposition(data);
        qasm += "rz(" + std::to_string(results.delta) + ") q[0];\n";
        qasm += "ry(" + std::to_string(results.gamma) + ") q[0];\n";
        qasm += "rz(" + std::to_string(results.beta) + ") q[0];\n";
        qasm += "gphase("+ std::to_string(results.alpha) +");";
    }
    return qasm;
}
