//
// Created by alex on 07/04/26.
//

#pragma once
#include "iast.h"
#include <Eigen/Dense>
#include <memory>
#include <vector>
#include "one_qubit_gate.h"
#include "synthesisOptions.h"

class csdNode;

class unitaryOneQubitGateNode final : public IASTnode {
public:
    explicit unitaryOneQubitGateNode(const zyz_result& zyz_params, int position = 0);
    void accept(nodeVisitor &visitor) override;
    return_type get_data() override;
    zyz_result params; 
    int position;
};

class twoQubitGateNode final : public IASTnode {
public:
    /**
     * @param matrix 4x4 unitary.
     * @param up_to_diagonal If true, implements D^dagger * matrix with two CNOTs
     *        (optimization A.2) and stores the diagonal D in `diagonal`, so that
     *        matrix = diagonal * (this circuit). Falls back to three CNOTs (and
     *        diagonal = identity) if the two-CNOT reconstruction is not exact.
     */
    explicit twoQubitGateNode(const Eigen::MatrixXcd& matrix, bool up_to_diagonal = false);
    void accept(nodeVisitor &visitor) override;
    return_type get_data() override;
    
    std::unique_ptr<IASTnode> gate_a;
    std::unique_ptr<IASTnode> gate_b;
    std::unique_ptr<IASTnode> gate_c;
    std::unique_ptr<IASTnode> gate_d;
    
    std::unique_ptr<IASTnode> rz_theta;
    std::unique_ptr<IASTnode> ry1_theta;
    std::unique_ptr<IASTnode> ry2_theta;
    
    double phase;

    /** True when the node uses the two-CNOT circuit (A.2). */
    bool two_cx = false;
    /** Diagonal left out of the circuit (identity unless two_cx). */
    Eigen::Matrix4cd diagonal = Eigen::Matrix4cd::Identity();
};

/*
 * unitary -> csd
 * unitary -> two_qubit_gate
 * unitary -> one_qubit_gate
  * @ params unitary matrix to be decomposed
 */
class unitaryGateNode final : public IASTnode {
public:
    /**
     * @param matrix Unitary matrix to be decomposed.
     * @param optimize Apply optimization A.2 to the whole tree. Only the root does
     *        it; the internal nodes of the decomposition are created with false.
     * @param defer_two_qubit Internal use: a two-qubit node only stores its matrix and
     *        is synthesized later by the root's applyA2(), so that every two-qubit
     *        block is decomposed exactly once. The root sets it for its descendants
     *        when it applies A.2; after the root's constructor the tree is complete.
     * @param method Decomposition of this node and its descendants with three or more qubits.
     * @param check_unitary Verify that `matrix` is unitary (O(d^3)). Internal use: the
     *        nodes of the decomposition pass false, because their matrices are unitary
     *        by construction; only the root checks its input.
     */
    explicit unitaryGateNode(const Eigen::MatrixXcd& matrix, bool optimize = true, bool defer_two_qubit = false,
                             unitaryDecomposition method = unitaryDecomposition::QSD, bool check_unitary = true);

    /**
     * @brief Decomposes a unitary matrix with the given options (this node is the root).
     *
     * Same as unitaryGateNode(matrix, options.a2, false, options.method).
     */
    unitaryGateNode(const Eigen::MatrixXcd& matrix, const synthesisOptions& options);
    void accept(nodeVisitor &visitor) override;
    return_type get_data() override;
    std::unique_ptr<IASTnode> csd;
    int get_num_qubits();
    static std::unique_ptr<IASTnode> createBaseUnitary(const Eigen::MatrixXcd& matrix, int position = 0);

private:    
    Eigen::MatrixXcd data;
    /** Children of this node defer their two-qubit blocks to the root's applyA2(). */
    bool defer_children = false;
    unitaryDecomposition method = unitaryDecomposition::QSD;
    static int index;
    int name;

    void init();
    void applyA2();
    void collectTwoQubitLeaves(std::vector<unitaryGateNode*>& leaves);
    std::unique_ptr<IASTnode> createCSD(const Eigen::MatrixXcd& matrix);
    std::unique_ptr<IASTnode> createZXZ(const Eigen::MatrixXcd& matrix);
    std::unique_ptr<IASTnode> createTwoQubitUnitary(const Eigen::MatrixXcd& matrix);
};