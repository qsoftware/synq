//
// Created by alex on 07/04/26.
//

#include "../include/csdNode.h"
#include "../include/csd.h"
#include "../include/nodeVisitor.h"
#include "../include/qsdNode.h"
#include "../include/ucryNode.h"    

csdNode::csdNode(const Eigen::MatrixXcd& matrix, bool defer_two_qubit)
    : m_defer_two_qubit(defer_two_qubit), m_matrix(matrix) {
    auto csdValues = csd(matrix, matrix.rows() / 2, matrix.cols() / 2);
    // Optimization A.1 (Shende, Bullock, Markov 2006): the RY multiplexor is built
    // with CZs and its last CZ is omitted. That CZ acts on the select qubit (MSB)
    // and on the MSB of the lower register, so it is absorbed into the next
    // multiplexor (V1T, V2T) as a Z on that qubit inside the V2T block, i.e. a sign
    // flip of the right half of V2T's columns.
    const Eigen::Index half = csdValues.V2T.cols() / 2;
    csdValues.V2T.rightCols(half) *= -1.0;
    qsd1 = createQSD(csdValues.V1T, csdValues.V2T);
    qsd2 = createQSD(csdValues.U1, csdValues.U2);
    Eigen::VectorXd theta_eigen = Eigen::Map<Eigen::VectorXd>(csdValues.theta.data(), csdValues.theta.size());
    Eigen::VectorXd angles_eigen = 2.0 * theta_eigen;
    std::vector<double> angles(angles_eigen.data(), angles_eigen.data() + angles_eigen.size());
    ucry = createUcry(angles);
    ucry->setUseCz(true);
    ucry->drop_last = true;
}

std::unique_ptr<qsdNode> csdNode::createQSD(const Eigen::MatrixXcd& gate1, const Eigen::MatrixXcd& gate2) {
    return std::make_unique<qsdNode>(gate1, gate2, m_defer_two_qubit);
}

std::unique_ptr<ucryNode> csdNode::createUcry(const std::vector<double>& angles) {
    return std::make_unique<firstUcryNode>(&angles, false);
}

void csdNode::accept(nodeVisitor &visitor) {
    visitor.visit(*this);
}

return_type csdNode::get_data() {
    return m_matrix;
}