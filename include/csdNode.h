/**
 * @file csdNode.h
 * @brief Definition of the csdNode class representing a Cosine-sine decomposition in the AST.
 * 
 * This file is part of the SynQ compiler frontend and defines the AST node 
 * responsible for representing an CSD applied to a matrix.
 * 
 * @author SynQ team
 */

#pragma once
#include <memory>
#include <vector>
#include <tuple>
#include "iast.h"
#include <Eigen/Dense>
#include "csd.h"
#include "qsdNode.h"
#include "ucryNode.h"

/**
 * @class csdNode   
 * @brief AST Node representing a Cosine-sine decomposition.
 * 
 * This node decomposes a matrix into one/two qubit gates using cosine-sine decomposition
 * using qsd (quantum shannon decomposition) nodes and an Ucry (Unitary controlled RY) node.
 * csd --> qsd ucry qsd
 * It is used in the AST to represent the decomposition of an arbitrary matrix in a quantum circuit.
 */

class csdNode final : public IASTnode {
public:
    /**
     * @brief Computes the CSD of the given matrix and initializes the child nodes (qsd1, qsd2, ucry).
     * 
     * @param matrix Unitary matrix that will be decomposed.
     */
    explicit csdNode(const Eigen::MatrixXcd& matrix);

    /**
     * @brief Accept a visitor according to the Visitor pattern.
     * 
     * @param visitor Reference to a nodeVisitor instance.
     */
    void accept(nodeVisitor &visitor) override;

    /**
     * @brief Retrieve the node's stored data.
     * 
     * @return A return_type object containing the og unitary matrix.
     */
    return_type get_data() override;

    /**
     * @brief first qsd of the decomposition circuit.
     */
    std::unique_ptr<qsdNode> qsd1;
    /**
     * @brief second qsd of the decomposition circuit.
     */
    std::unique_ptr<qsdNode> qsd2;
    /**
     * @brief Unitary controlled RY between both qsd nodes.
     */
    std::unique_ptr<ucryNode> ucry;

private:
    /**
     * @brief Factory method to create a unique pointer to an QSD node.
     * 
     * @param gate1 First unitary block matrix from the CSD output.
     * @param gate2 Second unitary block matrix from the CSD output.
     * @return std::unique_ptr<qsdNode> Newly crated QSD node wrapping both blocks.
     */
    std::unique_ptr<qsdNode> createQSD(const Eigen::MatrixXcd& gate1, const Eigen::MatrixXcd& gate2);
    /**
     * @brief Factory method to create a unique pointer to an Ucry node.
     * 
     * @param angles Vector of rotation angles (2*theta) derived from the CSD singular values.
     * @return std::unique_ptr<ucryNode> Newly created UCry node encoding the controlled RY rotations.
     */
    std::unique_ptr<ucryNode> createUcry(const std::vector<double>& angles);

    /**
     * @brief Check if a matrix is a generalized permutation matrix.
     * 
     * @param matrix The matrix to check.
     * @return bool True if each row and column has exactly one non-zero element.
     */
    bool is_generalized_permutation(const Eigen::MatrixXcd& matrix);

    /**
     * @brief Optimize synthesis when blocks are generalized permutation matrices by transforming them into diagonals.
     *
     * @param left_gate Vector with two left block matrices.
     * @param right_gate Vector with two right block matrices.
     * @param theta Vector of angles (theta) to transform.
     * @param side 'L' for left or 'R' for right optimization.
     * @return std::tuple<std::vector<Eigen::MatrixXcd>, std::vector<Eigen::MatrixXcd>, Eigen::VectorXd>
     *         New left/right gates and transformed theta.
     */
    std::tuple<std::vector<Eigen::MatrixXcd>, std::vector<Eigen::MatrixXcd>, Eigen::VectorXd>
    permutation_optimization(const std::vector<Eigen::MatrixXcd>& left_gate,
                             const std::vector<Eigen::MatrixXcd>& right_gate,
                             const Eigen::VectorXd& theta,
                             char side);

    /**
     * @brief Unitary base matrix for the decomposition
     */
    Eigen::MatrixXcd m_matrix;
};