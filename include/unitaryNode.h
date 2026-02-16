//
// Created by Adenilton on 14/02/26.
//

#pragma once
#include "iast.h"

class unitaryNode: public IASTnode {
public:
    unitaryNode(Eigen::MatrixXcf unitary_matrix);
    void accept(nodeVisitor &visitor) override;
    return_type get_data() override;
    std::unique_ptr<IASTnode> decomposition;
private:
    Eigen::Matrix2cf data;
};
