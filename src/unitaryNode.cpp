//
// Created by Adenilton on 14/02/26.
//

#include "../include/unitaryNode.h"

#include "../include/csdNode.h"
#include "../include/nodeVisitor.h"

/**
 * unitaryNode --> unitary_matrix
 * unitaryNode --> csdNode
 * @param unitary_matrix
 */
unitaryNode::unitaryNode(Eigen::MatrixXcf unitary_matrix) {
    {
        if (unitary_matrix.rows() == 2) {
            data = unitary_matrix;
            num_qubits = 1;
        }
        else {
            num_qubits = static_cast<int>(std::log2(unitary_matrix.rows()));
            decomposition = std::make_unique<csdNode>(unitary_matrix);

        }
    }
}

void unitaryNode::accept(nodeVisitor &visitor) {
    visitor.visit(*this);
}

return_type unitaryNode::get_data() {
    return data;
}
