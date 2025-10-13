//
// Created by adenilton on 10/4/25.
//

#pragma once
#include <memory>
#include "IASTnode.h"

class uGate:public IASTnode {

};

class rzNode final :public IASTnode {
public:
    explicit rzNode(double theta);
    rzNode();
    void accept(nodeVisitor &visitor) override;
    return_type get_data() override;
private:
    double angle;

};

/*
 * ucrz --> ucrz cx ucrz cx
 * ucrz --> rz cx ucrz rz
  * @ param angles vector with two double representing angles of uniformly controlled rz gate
 */
class ucrzNode:public IASTnode {
public:
    bool first=false;

    int name;
    explicit ucrzNode(const std::vector<double>* angles, bool _first, bool _reverse);
    void accept(nodeVisitor &visitor) override;
    return_type get_data() override;
    std::vector<double> angles;
    std::unique_ptr<IASTnode> gate1;
    std::unique_ptr<IASTnode> gate2;
};

/*
 * first_ucrz --> ucrz cx ucrz cx
 *              | rz cx rz cx
 */
class firstUcrzNode:public ucrzNode {
public:
    explicit firstUcrzNode(const std::vector<double>* angles);
    void accept(nodeVisitor &visitor) override;
};