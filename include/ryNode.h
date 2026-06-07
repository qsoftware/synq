/**
 * @file ryNode.h
 * @brief Definition of the ryNode class representing a quantum RY rotation gate in the AST.
 *
 * This file is part of the SynQ compiler frontend and defines the AST node
 * responsible for representing an RY rotation applied to a specific qubit.
 *
 * @author SynQ team 
 */

#pragma once
#include "iast.h"

/**
 * @class ryNode
 * @brief AST node representing a quantum RY rotation gate.
 *
 * This node stores a rotation angle (in radians) and an optional qubit position.
 * It is used in the AST to represent the application of an RY gate in a quantum circuit.
 */
class ryNode final :public IASTnode {
public:
    /**
     * @brief Construct an RY node with a given rotation angle and optional position.
     *
     * @param theta Rotation angle in radians.
     * @param position Qubit index where the gate is applied. Defaults to -1 (unspecified).
     */
    explicit ryNode(double theta, int position = -1);
    /**
     * @brief Default constructor.
     *
     * Initializes the angle to 0 and position to -1.
     */
    ryNode();

    /**
     * @brief Accept a visitor according to the Visitor pattern.
     *
     * @param visitor Reference to a nodeVisitor instance.
     */
    void accept(nodeVisitor &visitor) override;

    /**
     * @brief Retrieve the node's stored data.
     *
     * @return A return_type object containing the angle and position.
     */
    return_type get_data() override;

    /**
     * @brief Factory method to create a unique pointer to an RY node.
     *
     * @param angle Rotation angle in radians.
     * @param position Qubit index (default -1).
     * @return std::unique_ptr<IASTnode> Newly created RY node.
     */
    static std::unique_ptr<IASTnode> createRyNode(double angle, int position = -1);

    /** @brief Qubit index where the RY gate is applied. */
    int position;
private:
    /** @brief Rotation angle in radians. */
    double angle;
};