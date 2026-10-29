/**
 * @file hNode.h
 * @brief Definition of the hNode class representing a Hadamard gate in the AST.
 *
 * This file is part of the SynQ compiler frontend.
 *
 * @author SynQ team
 */

#pragma once
#include "iast.h"

/**
 * @class hNode
 * @brief AST node representing a Hadamard gate on a given qubit.
 */
class hNode final : public IASTnode {
public:
    /**
     * @brief Construct a Hadamard node.
     *
     * @param position Qubit index where the gate is applied.
     */
    explicit hNode(int position);

    /**
     * @brief Accept a visitor according to the Visitor pattern.
     *
     * @param visitor Reference to a nodeVisitor instance.
     */
    void accept(nodeVisitor &visitor) override;

    /**
     * @brief The gate has no parameters.
     *
     * @return An empty return_type (std::monostate).
     */
    return_type get_data() override;

    /**
     * @brief Factory method to create a unique pointer to a Hadamard node.
     *
     * @param position Qubit index.
     * @return std::unique_ptr<IASTnode> Newly created Hadamard node.
     */
    static std::unique_ptr<IASTnode> createHNode(int position);

    /** @brief Qubit index where the gate is applied. */
    int position;
};
