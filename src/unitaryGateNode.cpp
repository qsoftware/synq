//
// Created by alex on 07/04/26.
//

#include "../include/unitaryGateNode.h"
#include "../include/one_qubit_gate.h"
#include "../include/csdNode.h"
#include "../include/nodeVisitor.h"
#include "../include/twoQubitsUnitary.h"
#include "../include/ryNode.h"
#include "../include/rzNode.h"
#include "../include/qsdNode.h"
#include <stdexcept>
#include <cmath>

int unitaryGateNode::index = 0;

unitaryGateNode::unitaryGateNode(const Eigen::MatrixXcd& matrix, bool optimize, bool defer_two_qubit) : data(matrix) {

    if (matrix.size() == 0 || matrix.rows() != matrix.cols()) {
        throw std::invalid_argument("Matrix must be square and non-empty");
    }

    if ((matrix.rows() & (matrix.rows() - 1)) != 0) {
        throw std::invalid_argument("Matrix dimension must be a power of 2 (e.g., 2, 4, 8, 16...)");
    }

    if ((matrix * matrix.adjoint()).isIdentity(1e-10) == false) {
        throw std::invalid_argument("Matrix must be unitary");
    }

    num_qubits = static_cast<int>(log2(static_cast<double>(matrix.rows())));
    name = index++;
    const bool apply_a2 = optimize && num_qubits > 2;
    // With A.2, every two-qubit block below this node is synthesized once, by applyA2().
    defer_children = apply_a2 || defer_two_qubit;
    if (num_qubits == 2 && defer_two_qubit) {
        return;  // synthesized by the root's applyA2()
    }
    init();
    if (apply_a2) {
        applyA2();
    }
}

void unitaryGateNode::collectTwoQubitLeaves(std::vector<unitaryGateNode*>& leaves) {
    // Same (time) order as qasmVisitor: csd -> qsd2, ucry, qsd1; qsd -> gate_w, ucrz, gate_v.
    if (num_qubits == 2) {
        leaves.push_back(this);
        return;
    }
    auto* csd_node = dynamic_cast<csdNode*>(csd.get());
    if (csd_node == nullptr) return;
    for (qsdNode* qsd_node : {csd_node->qsd2.get(), csd_node->qsd1.get()}) {
        qsd_node->gate_w->collectTwoQubitLeaves(leaves);
        qsd_node->gate_v->collectTwoQubitLeaves(leaves);
    }
}

void unitaryGateNode::applyA2() {
    /*
     * Optimization A.2 (Shende, Bullock, Markov 2006): every two-qubit block except
     * the last is implemented up to a diagonal D with two CNOTs, and D is moved
     * forward into the next block (in time order). Between two consecutive blocks
     * there are only multiplexed rotations, and the CZs of A.1, whose targets are
     * qubits >= 2. They are diagonal in qubits 0 and 1, so D commutes with them.
     */
    std::vector<unitaryGateNode*> leaves;
    collectTwoQubitLeaves(leaves);

    Eigen::Matrix4cd carry = Eigen::Matrix4cd::Identity();
    for (std::size_t k = 0; k < leaves.size(); ++k) {
        const Eigen::Matrix4cd matrix = Eigen::Matrix4cd(leaves[k]->data) * carry;
        const bool last = (k + 1 == leaves.size());
        auto node = std::make_unique<twoQubitGateNode>(matrix, !last);
        carry = node->diagonal;
        leaves[k]->csd = std::move(node);
    }
}

void unitaryGateNode::init(){
    if (num_qubits > 2) {
        csd = createCSD(data);
    } else if (num_qubits == 2) {
        csd = createTwoQubitUnitary(data);
    } else {
        csd = createBaseUnitary(data);
    }
}

unitaryOneQubitGateNode::unitaryOneQubitGateNode(const zyz_result& zyz_params, int position) : params(zyz_params), position(position) {
}

twoQubitGateNode::twoQubitGateNode(const Eigen::MatrixXcd& matrix, bool up_to_diagonal) {
    num_qubits = 2;
    if (up_to_diagonal) {
        const Eigen::Matrix4cd D = diagonalForTwoCx(matrix);
        const auto two = twoCxDecomposition(D.adjoint() * Eigen::Matrix4cd(matrix));
        if (two.error < 1e-9) {
            two_cx = true;
            diagonal = D;
            gate_a = unitaryGateNode::createBaseUnitary(two.A, 1);
            gate_b = unitaryGateNode::createBaseUnitary(two.B, 0);
            gate_c = unitaryGateNode::createBaseUnitary(two.C, 1);
            gate_d = unitaryGateNode::createBaseUnitary(two.D, 0);
            rz_theta = rzNode::createRzNode(two.angleRz, 1);
            ry1_theta = ryNode::createRyNode(two.angleRy, 0);
            phase = two.phase;
            return;
        }
    }
    auto decomp_result = twoQubitsUnitaryDecomposition(matrix);
    
    gate_a = unitaryGateNode::createBaseUnitary(decomp_result.A, 0);
    gate_b = unitaryGateNode::createBaseUnitary(decomp_result.B, 1);
    gate_c = unitaryGateNode::createBaseUnitary(decomp_result.C, 1);
    gate_d = unitaryGateNode::createBaseUnitary(decomp_result.D, 0);
    
    rz_theta = rzNode::createRzNode(-1 * decomp_result.angleRz, 1); 
    ry1_theta = ryNode::createRyNode(-1 *decomp_result.angleRy1, 0);
    ry2_theta = ryNode::createRyNode( decomp_result.angleRy2, 0);
    
    phase = decomp_result.phase;
}

void unitaryOneQubitGateNode::accept(nodeVisitor &visitor) {
    visitor.visit(*this);
}

return_type unitaryOneQubitGateNode::get_data() {
    return std::vector<double>{
        params.alpha, 
        params.beta, 
        params.gamma, 
        params.delta
    }; 
}

std::unique_ptr<IASTnode> unitaryGateNode::createCSD(const Eigen::MatrixXcd& matrix) {
    return std::make_unique<csdNode>(matrix, defer_children);
}

std::unique_ptr<IASTnode> unitaryGateNode::createBaseUnitary(const Eigen::MatrixXcd& matrix, int position) {
    auto zyz = OneQubit::zyz_decomposition(matrix);
    return std::make_unique<unitaryOneQubitGateNode>(zyz, position);
}

std::unique_ptr<IASTnode> unitaryGateNode::createTwoQubitUnitary(const Eigen::MatrixXcd& matrix) {
    return std::make_unique<twoQubitGateNode>(matrix);
}

void unitaryGateNode::accept(nodeVisitor &visitor) {
    visitor.visit(*this);
}

int unitaryGateNode::get_num_qubits() {
    return num_qubits;
}

return_type unitaryGateNode::get_data() {
    return data;
}

void twoQubitGateNode::accept(nodeVisitor &visitor) {
    visitor.visit(*this);
}

return_type twoQubitGateNode::get_data() {
    return phase; 
}