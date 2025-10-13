//
// Created by adenilton on 10/12/25.
//

#pragma once

#include <Eigen/Dense>

#include<memory>
#include "../include/IASTnode.h"
#include "../include/ucrNode.h"

class csdNode final :public IASTnode {
public:
    bool first=false;

    int name;
    explicit csdNode(Eigen::MatrixXcf unitary);
    void accept(nodeVisitor &visitor) override;
    return_type get_data() override;
    Eigen::MatrixXcf matrix;
    std::unique_ptr<IASTnode> second_ucg1control;
    std::unique_ptr<IASTnode> mcry;
    std::unique_ptr<IASTnode> first_ucg1control;
};
