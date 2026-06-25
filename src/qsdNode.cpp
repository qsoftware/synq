//
// Created by alex on 07/04/26.
//>
#include "../include/qsdNode.h"
#include "../include/qsd.h"
#include "../include/nodeVisitor.h"
#include "../include/unitaryGateNode.h"
#include "../include/ucrzNode.h"
#include "../include/diagonalGateNode.h"


qsdNode::qsdNode(const Eigen::MatrixXcd& gate1, const Eigen::MatrixXcd& gate2) : m_gate1(gate1), m_gate2(gate2) {
    
    if (gate1.isDiagonal() && gate2.isDiagonal()){
        int size = 2 * gate1.rows();
        Eigen::MatrixXcd mux = Eigen::MatrixXcd::Zero(size, size);
        mux.topLeftCorner(size/2, size/2) = gate1;
        mux.bottomRightCorner(size/2, size/2) = gate2;

        gate_v = createDiagonal(mux);
    } else {
        qsd_result result = qsd(gate1, gate2);
        gate_v = createUnitary(result.gate_v);
        gate_w = createUnitary(result.gate_w);
        Eigen::VectorXd angles_list_d = -2.0 * result.list_d.array().arg();
        std::vector<double> angles(angles_list_d.data(), angles_list_d.data() + angles_list_d.size());
        ucrz = createUcrz(angles);
    }

}

std::unique_ptr<unitaryGateNode> qsdNode::createUnitary(const Eigen::MatrixXcd& matrix) {
    return std::make_unique<unitaryGateNode>(matrix);
}

std::unique_ptr<ucrzNode> qsdNode::createUcrz(const std::vector<double>& angles) {
    return std::make_unique<firstUcrzNode>(&angles, false);
}

std::unique_ptr<diagonalGateNode> qsdNode::createDiagonal(const Eigen::MatrixXcd& matrix) {
    return std::make_unique<diagonalGateNode>(matrix);
}
void qsdNode::accept(nodeVisitor &visitor) {
    visitor.visit(*this);
}

return_type qsdNode::get_data() {
    return std::make_pair(m_gate1, m_gate2);
}