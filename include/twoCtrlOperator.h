/**
 * @file twoCtrlOperator.h
 * @brief Definition of the twoCtrlOperator function representing a operation controlled by 2 control bits.
 *
 * This file is part of the SynQ compiler frontend and defines the AST node
 * responsible for representing an operator controlled by two control bits applied to a state.
 *
 * @author SynQ team 
 */

#pragma once
#include <Eigen/Dense>
#include <vector>
#include "iast.h"

/**
 * @class CtrlOperatorNode
 * @brief AST node representing a controlled operator with one or two control qubits
 *
 * This node provides methods to apply a unitary operation controlled by one or two control qubits on a target qubit
 */
class CtrlOperatorNode :public IASTnode{
public:
/**
* @brief Construct an CtrlOperator node with a given number of qubits
*
* @param ctrl_bits Vector of integers representing the positions of the control bits
* @param target_bit Integer representing the position of the target bit
* @param matrix Eigen::MatrixXcd representing the unitary operation to be applied
*/
explicit CtrlOperatorNode(std::vector<int>& ctrl_bits, int target_bit, Eigen::MatrixXcd& matrix, Eigen::Matrix2cd& Op);

/** @brief Return the number of control bits */
int get_num_ctrl();

/** @brief Return the position of a control bit at a given index */
int get_ctrl(int index);

/** @brief Return the position of the target bit */
int get_target();

/**
* @brief Accept a visitor according to the Visitor pattern.
* 
* @param visitor Reference to a nodeVisitor instance.
*/
void accept(nodeVisitor &visitor) override;

/**
 * @brief Retrieve the node's stored data.
 *
 * @return A return_type object containing the number of qubits.
 */
return_type get_data() override;

/**
* @brief Factory method to create a unique pointer to an CtrlOperatorNode.
*
* @param num_qubits description
* @return std::unique_ptr<IASTnode> Newly created CtrlOperator node.
*/
static std::unique_ptr<IASTnode> createCtrlOperatorNode(std::vector<int>& controls, int target, Eigen::MatrixXcd& CtrlMatrix, Eigen::Matrix2cd& OperatorMatrix);

/**
* @brief Apply a operator controlled by one control bit on a target bit
*
* @param U Unitary matrix representing the applied operation
* @param ctrl Integer representing the position of the control bit
* @param target Integer representing the position of the target bit
*/
Eigen::MatrixXcd oneCtrlOperator(const Eigen::Matrix2cd& U, int ctrl, int target);

/**
* @brief Apply a operator controlled by two control bits on a target bit
*
* @param U Unitary matrix representing the applied operation
* @param ctrl1 Integer representing the position of the first control bit
* @param ctrl2 Integer representing the position of the second control bit
*/
Eigen::MatrixXcd twoCtrlOperator(const Eigen::Matrix2cd& U, int ctrl1, int ctrl2);

/** @brief Unitary matrix representing the operation to be applied */
Eigen::Matrix2cd OperatorMatrix;

private:
/** @brief Vector of integers representing the positions of the control bits. */
std::vector<int> controls;

/** @brief Position of the target bit. */
int target;

/** @brief Matrix representing the current state of the control matrix. */
Eigen::MatrixXcd CtrlMatrix;

};
