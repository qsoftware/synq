#include "../include/zxzNode.h"

#include <cmath>
#include <vector>
#include "../include/hNode.h"
#include "../include/nodeVisitor.h"
#include "../include/qsd.h"
#include "../include/zxz.h"

namespace {

std::vector<double> rzAngles(const Eigen::VectorXcd& d) {
    // D ⊕ D† = RZ(-2 arg d_j) on the most significant qubit, controlled by j
    const Eigen::VectorXd angles = -2.0 * d.array().arg();
    return {angles.data(), angles.data() + angles.size()};
}

}  // namespace

zxzNode::zxzNode(const Eigen::MatrixXcd& matrix, const bool defer_two_qubit) : m_matrix(matrix) {
    num_qubits = static_cast<int>(std::log2(static_cast<double>(matrix.rows())));
    const Eigen::Index half = matrix.rows() / 2;
    constexpr auto method = unitaryDecomposition::BlockZXZ;

    const zxz_result zxz = blockZxz(matrix);
    const qsd_result c = qsd(Eigen::MatrixXcd::Identity(half, half), zxz.c);
    const qsd_result a = qsd(zxz.a1, zxz.a2);

    // V_C and W_A are merged into the middle multiplexor, and the CZs left by the two
    // omitted CNOTs are a Z on the MSB of the lower register in its second block.
    const Eigen::MatrixXcd b1 = a.gate_w * c.gate_v;
    Eigen::MatrixXcd b2 = a.gate_w * zxz.b * c.gate_v;
    const Eigen::Index quarter = half / 2;
    b2.bottomRows(quarter) *= -1.0;
    b2.rightCols(quarter) *= -1.0;

    gate_w_c = std::make_unique<unitaryGateNode>(c.gate_w, false, defer_two_qubit, method);

    const std::vector<double> angles_c = rzAngles(c.list_d);
    ucrz_c = std::make_unique<firstUcrzNode>(&angles_c);
    ucrz_c->drop_last = true;

    h1 = hNode::createHNode(num_qubits - 1);
    middle = std::make_unique<qsdNode>(b1, b2, defer_two_qubit, method);
    h2 = hNode::createHNode(num_qubits - 1);

    const std::vector<double> angles_a = rzAngles(a.list_d);
    ucrz_a = std::make_unique<firstUcrzNode>(&angles_a);
    ucrz_a->mirror = true;
    ucrz_a->drop_first = true;

    gate_v_a = std::make_unique<unitaryGateNode>(a.gate_v, false, defer_two_qubit, method);
}

void zxzNode::accept(nodeVisitor &visitor) {
    visitor.visit(*this);
}

return_type zxzNode::get_data() {
    return m_matrix;
}
