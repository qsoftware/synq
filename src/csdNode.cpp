//
// Created by adenilton on 10/12/25.
//

#include "../include/csdNode.h"

#include "../csd/csd.h"
#include "../include/qsdNode.h"


csdNode::csdNode(Eigen::MatrixXcf unitary) {
    const int rows = static_cast<int>(unitary.rows());
    if (rows == 2) {
        first_ucg1control = nullptr;
        second_ucg1control = nullptr;
        mcry = nullptr;
        matrix = unitary;
        return;
    }

    int p = rows / 2;
    auto csd_result = csd(unitary, p, p);

    // TODO: implementar mcry
    mcry = std::make_unique<firstUcrzNode>(&(csd_result.theta));

    // TODO implementar qsdNode
    first_ucg1control = std::make_unique<qsdNode>(csd_result.U1, csd_result.U2);
    second_ucg1control = std::make_unique<qsdNode>(csd_result.V1T, csd_result.V2T);


}
