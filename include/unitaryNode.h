//
// Created by Adenilton on 14/02/26.
//

#pragma once
#include "iast.h"

class unitaryNode: public IASTnode {
public:
    unitaryNode(Eigen::Matrix2cf unitary_matrix);
    void accept(nodeVisitor &visitor) override;
    return_type get_data() override;
private:
    Eigen::Matrix2cf data;
};
