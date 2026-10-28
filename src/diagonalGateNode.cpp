//
// Created by alex on 22/05/26.
//

#include "../include/diagonalGateNode.h"
#include "../include/nodeVisitor.h"
#include "../include/rzNode.h"
#include <stdexcept>
#include <cmath>
#include <complex>

using namespace std::complex_literals;

int diagonalGateNode::index = 0;

diagonalGateNode::diagonalGateNode(const Eigen::MatrixXcd& matrix) : data(matrix) {
    if (matrix.size() == 0 || matrix.rows() != matrix.cols()) {
        throw std::invalid_argument("Matrix must be square and non-empty");
    }

    if ((matrix.rows() & (matrix.rows() - 1)) != 0) {
        throw std::invalid_argument("Matrix dimension must be a power of 2 (e.g., 2, 4, 8, 16...)");
    }

    num_qubits = static_cast<int>(std::log2(static_cast<double>(matrix.rows())));
    name = index++;
    
    Eigen::VectorXcd diagonal = data.diagonal();
    Eigen::VectorXd theta = diagonal.array().arg();

    if (num_qubits == 1) {
        child_diagonal = oneQubitDiagonalGateNode::createDiagonalOneQubitGate(theta);
        return;
    }

    int half = theta.size() / 2;
    Eigen::VectorXd new_theta(half);
    Eigen::VectorXd rz_angles(half);

    for (int i = 0; i < theta.size()/2; i++) {
        new_theta(i) = (theta(i) + theta(i + half)) / 2.0;
        rz_angles(i) = (theta(i + half) - theta(i));
    }

    Eigen::MatrixXcd new_diagonal = (1.0i * new_theta.array()).exp().matrix().asDiagonal();
    std::vector<double> rz_vec(rz_angles.data(), rz_angles.data() + rz_angles.size());

    child_diagonal = createNextDiagonalGate(new_diagonal);
    ucrz = createUcrz(rz_vec);
}

std::unique_ptr<ucrzNode> diagonalGateNode::createUcrz(const std::vector<double> &angles) {
    return std::make_unique<firstUcrzNode>(&angles, false);
}

void diagonalGateNode::accept(nodeVisitor &visitor) {
    visitor.visit(*this);
}

return_type diagonalGateNode::get_data() {
    return data;
}

oneQubitDiagonalGateNode::oneQubitDiagonalGateNode(double global_phase, double p_phase) : global_phase(global_phase), p_phase(p_phase){
    num_qubits = 1;
    name = diagonalGateNode::index++;
}

std::unique_ptr<IASTnode> diagonalGateNode::createNextDiagonalGate(const Eigen::MatrixXcd& matrix) {
    return std::make_unique<diagonalGateNode>(matrix); 
}

std::unique_ptr<IASTnode> oneQubitDiagonalGateNode::createDiagonalOneQubitGate(const Eigen::VectorXd& angles) {
    double global_phase = angles(0);
    double p_phase = angles(1) - global_phase;
    return std::make_unique<oneQubitDiagonalGateNode>(global_phase, p_phase);
}

void oneQubitDiagonalGateNode::accept(nodeVisitor &visitor) {
    visitor.visit(*this);
}

return_type oneQubitDiagonalGateNode::get_data() {
    return std::vector<double>({global_phase, p_phase});
}