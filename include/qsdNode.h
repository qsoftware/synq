/**
 * @file qsdNode.h
 * @brief Definition of the qsdNode class representing a Quantum Shannon Decomposition in the AST.
 *
 * This file is part of the SynQ compiler frontend and defines the AST node
 * responsible for representing a QSD applied to a unitary matrix.
 *
 * @author SynQ team 
 */

#pragma once
#include <memory>
#include <vector>
#include "iast.h"
#include <Eigen/Dense>
#include "qsd.h"
#include "unitaryGateNode.h"
#include "ucrzNode.h"
#include "diagonalGateNode.h"

/**
 * @class qsdNode
 * @brief AST node representing a Quantum Shannon Decomposition.
 *
 * This node decomposes two unitary matrixes into two unitary gates nodes and an UCRZ node
 * qsd --> unitary ucrz unitary
  * @ params two unitary matrices corresponding to the mux gate in the CSD decomposition
 */
class qsdNode final : public IASTnode {
public:
    /**
     * @brief Construct an qsd node with two given unitary matrixes.
     *
     * @param gate1 First unitary matrix.
     * @param gate2 Second unitary matrix.
     */
    explicit qsdNode(const Eigen::MatrixXcd& gate1, const Eigen::MatrixXcd& gate2);
    
    /**
    * @brief Accept a visitor according to the Visitor pattern.
    * 
    * @param visitor Reference to a nodeVisitor instance.
    */
    void accept(nodeVisitor &visitor) override;

    /**
     * @brief Retrieve the node's stored data.
     *
     * @return A return_type object containing the two unitary gates.
     */
    return_type get_data() override;

    /**
    * @brief First unitary gate of the decomposed circuit
    */
    std::unique_ptr<IASTnode> gate_v;
    /**
    * @brief Second unitary gate of the decomposed circuit
    */
    std::unique_ptr<IASTnode> gate_w;
    /**
    * @brief Unitary controlled RZ of the decomposed circuit
    */
    std::unique_ptr<ucrzNode> ucrz;

private:
    /**
     * @brief Factory method to create a unique pointer to an unitary gate node.
     *
     * @param matrix An unitary matrix.
     * @return std::unique_ptr<IASTnode> Newly created unitaryGate node.
     */
    std::unique_ptr<unitaryGateNode> createUnitary(const Eigen::MatrixXcd& matrix);
    /**
     * @brief Factory method to create a unique pointer to an Ucrz node.
     *
     * @param angles Vector of rotation angles derived from the unitary gates values.
     * @return std::unique_ptr<IASTnode> Newly created Ucrz node encoding the controlled RZ rotations.
     */
    std::unique_ptr<ucrzNode> createUcrz(const std::vector<double>& angles);
    /**
     * @brief Factory method to create a unique pointer to a diagonal gate node.
     *
     * @param matrix A diagonal matrix representing a multi-qubit diagonal gate.
     * @return std::unique_ptr<diagonalGateNode> Newly created diagonalGateNode.
     */
    std::unique_ptr<diagonalGateNode> createDiagonal(const Eigen::MatrixXcd& matrix);

    /**
    * @brief First unitary base matrix for the decomposition
    */
    Eigen::MatrixXcd m_gate1;
    /**
    * @brief Second unitary base matrix for the decomposition
    */ 
    Eigen::MatrixXcd m_gate2;
};