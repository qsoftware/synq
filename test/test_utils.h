#ifndef TEST_UTILS_H
#define TEST_UTILS_H

/**
 * @file test_utils.h
 * @brief Utility functions for testint QASM generation and Matrix representations.
 * 
 * This file is part of the SynQ test suite. It provides helper functions
 * to convert generated QASM code back into Eigen Matrices, allowing verification of the 
 * unitary operations represented by the AST nodes.
 * 
 * @author SynQ team
 */

#pragma once

#include <string>
#include <Eigen/Dense>

/**
 * @brief Expands a single-qubit gate into an n-qubit Hilbert space.
 * 
 * Applies the Kronecker product to place the given 2x2 gate on the target qubit,
 * while applying the Identity matrix to the remaining qubits.
 * 
 * @param gate The 2x2 unitary matrix representing the single-qubit gate.
 * @param target The index of the target qubit (0-indexed)
 * @param n_qubits The total number of qubits in the system.
 * @return Eigen::MatrixXcd the expanded 2^n x 2^n matrix.
 */
Eigen::MatrixXcd expand_single(const Eigen::Matrix2cd& gate, int target, int n_qubits);

/**
 * @brief Constructs the matrix representation of a Controlled-X (CNOT) gate.
 * 
 * Builds the full 2^n x 2^n matrix for a CNOT operation between a control 
 * and a target qubit within an n-qubit system.
 * 
 * @param control The index of the control qubit.
 * @param target The index of the target qubit.
 * @param n_qubits The total number of qubits in the system.
 * @return Eigen::MatrixXcd The resulting 2^n x 2^n CNOT matrix.
 */
Eigen::MatrixXcd cx_mat(int control, int target, int n_qubits);

/**
 * @brief Converts a single-qubit QASM code string into its 2x2 matrix representation.
 * 
 * Parses basic QASM instructions (like rz, ry, gphase) and multiplies their
 * corresponding matrices to reconstruct the final unitary operation.
 * 
 * @param qasm_code The string containing the generated QASM instructions.
 * @return Eigen::Matrix2cd The reconstructed 2x2 unitary matrix.
 */
Eigen::Matrix2cd one_qubit_qasm_to_matrix(const std::string& qasm_code);

/**
 * @brief Converts a general QASM code string into its full matrix representation.
 * 
 * Parses multi-qubit QASM instructions (rz, ry, cx, gphase) and builds the full
 * 2^n x 2^n system matrix by multiplying the expanded individual operations.
 * 
 * @param qasm_code The string containing the generated QASM instructions.
 * @param n_qubits The total number of qubits in the system.
 * @return Eigen::MatrixXcd The reconstructed 2^n x 2^n unitary matrix.
 */
Eigen::MatrixXcd qasm_to_matrix(const std::string& qasm_code, int n_qubits);

#endif 