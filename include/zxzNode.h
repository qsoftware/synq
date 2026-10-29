/**
 * @file zxzNode.h
 * @brief Definition of the zxzNode class representing a Block-ZXZ decomposition in the AST.
 *
 * This file is part of the SynQ compiler frontend.
 *
 * @author SynQ team
 */

#pragma once
#include <memory>
#include <Eigen/Dense>
#include "iast.h"
#include "qsdNode.h"
#include "ucrzNode.h"
#include "unitaryGateNode.h"

/**
 * @class zxzNode
 * @brief AST node representing the Block-ZXZ decomposition (Krol and Al-Ars, 2024).
 *
 * U = (A1 ⊕ A2) (H ⊗ I) (I ⊕ B) (H ⊗ I) (I ⊕ C), H on the most significant qubit.
 * The three multiplexors are demultiplexed as in the QSD, mux(U1, U2) = V (D ⊕ D†) W,
 * where D ⊕ D† is an RZ multiplexor on the most significant qubit:
 *
 *     mux(I, C)  = V_C (D_C ⊕ D_C†) W_C
 *     mux(A1, A2) = V_A (D_A ⊕ D_A†) W_A
 *
 * V_C and W_A act on the lower qubits only, so they commute with the H gates and are
 * merged into the middle multiplexor. The last CNOT of the RZ multiplexor of C and
 * the first CNOT of the RZ multiplexor of A (which is emitted in mirrored order) act
 * on (lower MSB -> MSB); next to an H on the target they become CZs, which are
 * absorbed into the middle multiplexor as a Z on the MSB of the lower register:
 *
 *     middle = mux(W_A V_C, Z W_A B V_C Z).
 *
 * zxz --> unitary(W_C) ucrz(D_C, no last CX) h qsd(middle) h ucrz(D_A, mirrored,
 *         no first CX) unitary(V_A)
 *
 * Each step uses 3·2^(n-1) - 2 CNOTs besides the blocks on n - 1 qubits (no RY
 * multiplexor), one less than the QSD with optimization A.1 and two less than
 * without it.
 */
class zxzNode final : public IASTnode {
public:
    /**
     * @param matrix Unitary matrix with three or more qubits.
     * @param defer_two_qubit Passed to the unitary children (see unitaryGateNode).
     */
    explicit zxzNode(const Eigen::MatrixXcd& matrix, bool defer_two_qubit = false);

    void accept(nodeVisitor &visitor) override;

    /** @return The decomposed matrix. */
    return_type get_data() override;

    // Children in time order.
    /** @brief W_C, on the lower qubits. */
    std::unique_ptr<unitaryGateNode> gate_w_c;
    /** @brief RZ multiplexor D_C ⊕ D_C† without its last CNOT. */
    std::unique_ptr<ucrzNode> ucrz_c;
    /** @brief First H on the most significant qubit. */
    std::unique_ptr<IASTnode> h1;
    /** @brief Middle multiplexor mux(W_A V_C, Z W_A B V_C Z). */
    std::unique_ptr<qsdNode> middle;
    /** @brief Second H on the most significant qubit. */
    std::unique_ptr<IASTnode> h2;
    /** @brief RZ multiplexor D_A ⊕ D_A†, mirrored, without its first CNOT. */
    std::unique_ptr<ucrzNode> ucrz_a;
    /** @brief V_A, on the lower qubits. */
    std::unique_ptr<unitaryGateNode> gate_v_a;

private:
    Eigen::MatrixXcd m_matrix;
};
