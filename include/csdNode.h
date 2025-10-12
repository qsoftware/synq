//
// Created by adenilton on 10/12/25.
//

#ifndef TEST_CSDNODE_H
#define TEST_CSDNODE_H

#include "../include/ucrzNode.h"

class csdNode final :public IASTnode {
public:
    bool first=false;

    int name;
    csdNode(const std::vector<double>* angles, bool _first, bool _reverse);
    void accept(nodeVisitor &visitor) override;
    return_type get_data() override;
    std::vector<double> angles;
    std::unique_ptr<IASTnode> second_ucg1control;
    std::unique_ptr<IASTnode> mcry;
    std::unique_ptr<IASTnode> first_ucg1control;
};


#endif //TEST_CSDNODE_H