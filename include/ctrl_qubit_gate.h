//
// Created by adenilton on 10/4/25.
//

/**
 * @file ctrl_qubit_gate.h
 * @brief Definition of the ctrl_qubit_gate class representing a controlled arbitrary unitary gate in the AST.
 *
 * This file is part of the SynQ compiler frontend and defines the abstract AST node
 * responsible for representing and decomposing a controlled arbitrary unitary gate.
 *
 * @author SynQ team
 */

#ifndef CTRL_QUBIT_GATE_H
#define CTRL_QUBIT_GATE_H

#include <Eigen/Dense>
#include "iast.h"
#include "one_qubit_gate.h"

/**
 * @class ctrl_qubit_gate
 * @brief Abstract AST Node representing a controlled arbitrary unitary gate.
 *
 * This node decomposes a controlled unitary into A, B, and C matrices.
 * Code generation is delegated to the nodeVisitor.
 */
class ctrl_qubit_gate : public IASTnode {
public:
    /** @brief Index of the control qubit. */
    int control;
    
    /** @brief Index of the target qubit. */
    int target;
    
    /** @brief Original 2x2 unitary matrix. */
    Eigen::Matrix2cd data;
    
    /** @brief Extracted Euler angles from ZYZ decomposition. */
    zyz_result angles;

    /** @brief Global phase alpha. */
    double phase_alpha;

    /** @brief Matrix A from the decomposition. */
    Eigen::Matrix2cd matrix_A;
    
    /** @brief Matrix B from the decomposition. */
    Eigen::Matrix2cd matrix_B;
    
    /** @brief Matrix C from the decomposition. */
    Eigen::Matrix2cd matrix_C;

    /**
     * @brief Construct a controlled unitary gate node.
     *
     * @param control_qubit Index of the control qubit.
     * @param target_qubit Index of the target qubit.
     * @param uMatrix The 2x2 unitary matrix to be decomposed.
     */
    explicit ctrl_qubit_gate(int control_qubit, int target_qubit, Eigen::Matrix2cd uMatrix);

    /**
     * @brief Accept a visitor according to the Visitor pattern.
     *
     * @param visitor Reference to a nodeVisitor instance.
     */
    void accept(nodeVisitor &visitor) override;

    /**
     * @brief Retrieve the node's stored unitary matrix.
     *
     * @return A return_type object containing the original matrix.
     */
    return_type get_data() override;

private:
    /**
     * @brief Performs the ZYZ decomposition and calculates A, B, and C matrices.
     */
    void decompose();
};

#endif //CTRL_QUBIT_GATE_H