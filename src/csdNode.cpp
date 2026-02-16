//
// Created by Adenilton on 14/02/26.
//

#include "../include/csdNode.h"
#include "../include/csd.h"
#include "../include/qsd.h"
#include "../include/unitaryNode.h"
#include "../include/ucrzNode.h"
#include "../include/ucryNode.h"
#include "../include/nodeVisitor.h"
#include<iostream>


csdNode::csdNode(Eigen::MatrixXcf& unitary_matrix) {
    long p = unitary_matrix.rows() / 2;
    auto csd_result = csd(unitary_matrix, p, p);
    left_ucg = std::make_unique<qsdNode>(csd_result.U1, csd_result.U2);
    std::cout << csd_result.U1 << "\n\n" << csd_result.U2 << std::endl;
    right_ucg = std::make_unique<qsdNode>(csd_result.V1T, csd_result.V2T);

    left_ucg = std::make_unique<qsdNode>(csd_result.V1T, csd_result.V2T);
    mcry = std::make_unique<firstUcryNode>(&csd_result.theta);

}

void csdNode::accept(nodeVisitor &visitor) {
    visitor.visit(*this);
}

return_type csdNode::get_data() {
    return 0.0;
}

qsdNode::qsdNode(Eigen::MatrixXcf unitary_matrix1, Eigen::MatrixXcf unitary_matrix2) {
    auto result = qsd(unitary_matrix1, unitary_matrix2);
    left_unitary = std::make_unique<unitaryNode>(result.left_gate);
    right_unitary = std::make_unique<unitaryNode>(result.right_gate);

    std::vector<double> angles(result.list_d.size());
    for (size_t i = 0; i < result.list_d.size(); ++i) {
        angles[i] =  -2 * std::arg(result.list_d[i]);
    }
    mcrz = std::make_unique<firstUcrzNode>(&angles);

}

void qsdNode::accept(nodeVisitor &visitor) {
    visitor.visit(*this);
}

return_type qsdNode::get_data() {
    return 0.0;
}
