/**
 * @file ucRotationNode.h
 * @brief Definition of the UCRotationNode class representing a uniformly controlled rotation in the AST.
 *
 * This file is part of the SynQ compiler frontend and defines the abstract AST node
 * responsible for representing uniformly controlled rotations (UCR).
 *
 * @author SynQ team
 */

#pragma once

#include <vector>
#include "iast.h"

/**
 * @class UCRotationNode
 * @brief Abstract AST Node representing a Uniformly Controlled Rotation.
 *
 * This node decomposes a UCR into smaller UCR sub-nodes or base rotations 
 * using its init() method.
 */
class UCRotationNode : public IASTnode {
public:
    /** @brief Flag indicating if the node uses the first decomposition rule. */
    bool first = false;
    
    /** @brief Flag indicating if the inverse of the gate should be applied. */
    bool inverse = false;

    /**
     * @brief Use CZ instead of CNOT for every two-qubit gate of this node and its
     * sub-nodes (valid for RY multiplexors: Z Ry(t) Z = Ry(-t), and CZs commute).
     */
    bool use_cz = false;

    /**
     * @brief Omit the final two-qubit gate (optimization A.1). Only meaningful for
     * a "first" node with use_cz; the caller absorbs the omitted CZ elsewhere.
     */
    bool drop_last = false;

    /**
     * @brief Emit the gates of this RZ multiplexor in reverse order, so that it starts
     * with its outer CNOT. A multiplexed RZ is diagonal and its circuit has only RZ and
     * CNOT gates, which are symmetric matrices, so the reversed circuit (the transpose)
     * implements the same operator. Only meaningful for a "first" RZ node.
     */
    bool mirror = false;

    /**
     * @brief Omit the first two-qubit gate of a mirrored node (Block-ZXZ); the caller
     * absorbs the omitted CNOT elsewhere.
     */
    bool drop_first = false;

    /** @brief Sets use_cz on this node and all its UCR sub-nodes. */
    void setUseCz(bool value);
    
    /** @brief Unique index assigned to the node. */
    int name;
    
    /** @brief Vector of rotation angles. */
    std::vector<double> angles;
    
    /** @brief First child node from the decomposition. */
    std::unique_ptr<IASTnode> gate1;
    
    /** @brief Second child node from the decomposition. */
    std::unique_ptr<IASTnode> gate2;

    /**
     * @brief Construct a UCR node with given angles and flags.
     *
     * @param angles Pointer to a vector of rotation angles.
     * @param _first Flag for the first decomposition rule.
     * @param _reverse Flag to reverse the sub-node creation order.
     * @param _inverse Flag to apply the inverse of the gate.
     */
    explicit UCRotationNode(const std::vector<double>* angles, bool _first, bool _reverse, bool _inverse);

    /**
     * @brief Accept a visitor according to the Visitor pattern.
     *
     * @param visitor Reference to a nodeVisitor instance.
     */
    void accept(nodeVisitor &visitor) override;

    /**
     * @brief Retrieve the node's stored angles.
     *
     * @return A return_type object containing the rotation angles.
     */
    return_type get_data() override;

    /**
     * @brief Create a base single-qubit rotation node.
     *
     * @param angle The rotation angle.
     * @return std::unique_ptr<IASTnode> Newly created base rotation node.
     */
    virtual std::unique_ptr<IASTnode> createBaseRotation(double angle) = 0;

    /**
     * @brief Create a sub-node for recursive decomposition.
     *
     * @param subAngles Pointer to the subset of angles.
     * @param reverse Flag to reverse the order in the sub-node.
     * @return std::unique_ptr<UCRotationNode> Newly created child node.
     */
    virtual std::unique_ptr<UCRotationNode> createSubNode(const std::vector<double>* subAngles, bool reverse) = 0;

    /**
     * @brief Initializes the node by splitting angles and creating child nodes.
     */
    void init();
};
