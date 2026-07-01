#include "test_utils.h"
#include "../include/one_qubit_gate.h"
#include <sstream>
#include <unsupported/Eigen/KroneckerProduct>
#include <iostream>
#include <iomanip>

Eigen::MatrixXcd expand_single(const Eigen::Matrix2cd& gate, int target, int n_qubits){
    Eigen::Matrix2cd I2 = Eigen::Matrix2cd::Identity();
    Eigen::MatrixXcd result(1, 1);
    result(0, 0) = 1.0;

    for (int q = n_qubits-1; q >= 0; q--){
        const Eigen::Matrix2cd& factor = (q == target) ? gate : I2;
        Eigen::MatrixXcd next = Eigen::kroneckerProduct(result, factor).eval();
        result = next;
    }
    return result;
}

Eigen::MatrixXcd cx_mat(int control, int target, int n_qubits){
    OneQubit gate = OneQubit();
    Eigen::Matrix2cd I, P0, P1, X;
    I << 1, 0, 0, 1;
    P0 << 1, 0, 0, 0;
    P1 << 0, 0, 0, 1;
    X = gate.x_matrix();

    Eigen::MatrixXcd term0 = Eigen::MatrixXcd::Identity(1, 1);
    Eigen::MatrixXcd term1 = Eigen::MatrixXcd::Identity(1, 1);

    for (int i = n_qubits-1; i >= 0; i--){
        Eigen::Matrix2cd op0 = I;
        Eigen::Matrix2cd op1 = I;

        if (i == control) {
            op0 = P0;
            op1 = P1;
        } else if (i == target){
            op1 = X;
        }

        term0 = Eigen::kroneckerProduct(term0, op0).eval();
        term1 = Eigen::kroneckerProduct(term1, op1).eval();
    }
    return term0 + term1;
}

Eigen::Matrix2cd one_qubit_qasm_to_matrix(const std::string& qasm_code) {
    OneQubit gate = OneQubit();
    Eigen::Matrix2cd U = Eigen::Matrix2cd::Identity();

    std::istringstream iss(qasm_code);
    std::string line;

    auto parse_angle = [](const std::string& l) -> double {
        size_t open  = l.find('(');
        size_t close = l.find(')');
        return std::stod(l.substr(open + 1, close - open - 1));
    };

    while (std::getline(iss, line)) {
        if (line.find("rz") != std::string::npos) {
            U = gate.rz_matrix(parse_angle(line)) * U;
        } else if (line.find("ry") != std::string::npos) {
            U = gate.ry_matrix(parse_angle(line)) * U;
        }
        if (line.find("gphase") != std::string::npos) {
            U *= gate.gphase(parse_angle(line));
        }
    }
    return U;
}

Eigen::MatrixXcd qasm_to_matrix(const std::string& qasm_code, int n_qubits){
    OneQubit gate = OneQubit();
    int dim = 1 << n_qubits;
    Eigen::MatrixXcd U = Eigen::MatrixXcd::Identity(dim, dim);

    std::istringstream iss(qasm_code);
    std::string line;

    auto parse_angle = [](const std::string& l) -> double {
        size_t open  = l.find('(');
        size_t close = l.find(')');
        return std::stod(l.substr(open + 1, close - open - 1));
    };

    auto parse_qubit = [](const std::string& l, size_t start) -> int{
        size_t open = l.find('[', start);
        size_t close = l.find(']', open);
        return std::stoi(l.substr(open + 1, close - open - 1));
    };

    while (std::getline(iss, line)) {
        if (line.find("rz") != std::string::npos){
            double angle = parse_angle(line);
            int qubit = parse_qubit(line, 0);
            U = expand_single(gate.rz_matrix(angle), qubit, n_qubits) * U;
        } else if (line.find("ry") != std::string::npos) {
            double angle = parse_angle(line);
            int qubit    = parse_qubit(line, 0);
            U = expand_single(gate.ry_matrix(angle), qubit, n_qubits) * U;
        } else if (line.find("cx") != std::string::npos) {
            int ctrl = parse_qubit(line, 0);
            size_t comma = line.find(',');
            int tgt = parse_qubit(line, comma);
            U = cx_mat(ctrl, tgt, n_qubits) * U;
        } else if (line.find("gphase") != std::string::npos){
            U *= gate.gphase(parse_angle(line));
        } else if (line.find("p") != std::string::npos){
            double angle = parse_angle(line);
            int qubit    = parse_qubit(line, 0);
            U = expand_single(gate.p_matrix(angle), qubit, n_qubits) * U;
        }
    }
    return U;
}

void printQuantumMatrix(const Eigen::MatrixXcd& matrix, const std::string& name) {
    std::cout << std::fixed << std::setprecision(2); 
    std::cout << name << " matrix:\n";
    for (int i = 0; i < matrix.rows(); ++i) {
        for (int j = 0; j < matrix.cols(); ++j) {
            std::complex<double> val = matrix(i, j);
            
            double real_part = (std::abs(val.real()) < 1e-6) ? 0.0 : val.real();
            double imag_part = (std::abs(val.imag()) < 1e-6) ? 0.0 : val.imag();

            std::cout << "(" << real_part << "," << imag_part << ")\t";
        }
        std::cout << "\n";
    }
}