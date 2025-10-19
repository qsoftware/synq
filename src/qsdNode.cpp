//
// Created by adenilton on 10/11/25.
//

#include "../include/qsdNode.h"

#include "../include/csdNode.h"

qsdNode::qsdNode(Eigen::MatrixXcf gate1, Eigen::MatrixXcf gate2) {
    num_qubits = static_cast<int>(log2(gate1.rows())) + 1;
    auto output = _compute_gate(gate1, gate2);
    first_unitary = new csdNode(output.gate_v);
    second_unitary = new csdNode(output.gate_w);
    mcrz = new firstUcrNode(&output.list_d, Z);

}

void qsdNode::accept(nodeVisitor &visitor) {
}

return_type qsdNode::get_data() {
    return_type a;
    return a;
}

qsd_result qsdNode::_compute_gate(Eigen::MatrixXcf gate1, Eigen::MatrixXcf gate2) {


    Eigen::MatrixXcf c = gate1 * gate2.adjoint();
    Eigen::ComplexEigenSolver<Eigen::MatrixXcf> ces(c);

    auto gate_v = ces.eigenvectors();
    auto d_square = ces.eigenvalues();
    auto eig_list_d = d_square.array().sqrt();
    std::vector<double> list_d2(d_square.size());
    for (int i = 0; i < d_square.size(); ++i) {
        list_d2[i] = std::sqrt(eig_list_d[i]).real();
    }
    qsd_result b;
    b.gate_v = gate_v;
    b.list_d = list_d2;

    b.gate_w = gate_v.adjoint() * gate2;
    return b;
}
