//
// Created by alex on 07/04/26.
//>
#include "../include/qsdNode.h"
#include "../include/qsd.h"
#include "../include/nodeVisitor.h"
#include "../include/unitaryGateNode.h"
#include "../include/ucrzNode.h"


qsdNode::qsdNode(const Eigen::MatrixXcd& gate1, const Eigen::MatrixXcd& gate2, bool defer_two_qubit,
                 unitaryDecomposition method)
    : m_defer_two_qubit(defer_two_qubit), m_method(method), m_gate1(gate1), m_gate2(gate2) {
    qsd_result result = qsd(gate1, gate2);
    gate_v = createUnitary(result.gate_v);
    gate_w = createUnitary(result.gate_w);
    Eigen::VectorXd angles_list_d = -2.0 * result.list_d.array().arg();
    std::vector<double> angles(angles_list_d.data(), angles_list_d.data() + angles_list_d.size());
    ucrz = createUcrz(angles);
}

std::unique_ptr<unitaryGateNode> qsdNode::createUnitary(const Eigen::MatrixXcd& matrix) {
    return std::make_unique<unitaryGateNode>(matrix, false, m_defer_two_qubit, m_method);
}

std::unique_ptr<ucrzNode> qsdNode::createUcrz(const std::vector<double>& angles) {
    return std::make_unique<firstUcrzNode>(&angles, false);
}

void qsdNode::accept(nodeVisitor &visitor) {
    visitor.visit(*this);
}

return_type qsdNode::get_data() {
    return std::make_pair(m_gate1, m_gate2);
}