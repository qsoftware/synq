//
// Created by Adenilton on 14/02/26.
//

#include "../include/unitaryNode.h"

#include "../include/nodeVisitor.h"

/**
 * unitaryNode --> unitary_matrix
 * unitaryNode --> csdNode
 * @param unitary_matrix
 */
unitaryNode::unitaryNode(Eigen::Matrix2cf unitary_matrix) {
    {
        if (unitary_matrix.rows() == 2) {
            data = unitary_matrix;
            num_qubits = 1;
        }
        else {
            data = Eigen::Matrix2cf::Zero(0, 0);
        }
    }
}

void unitaryNode::accept(nodeVisitor &visitor) {
    visitor.visit(*this);
}

return_type unitaryNode::get_data() {
    return data;
}
