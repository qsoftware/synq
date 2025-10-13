//
// Created by adenilton on 10/4/25.
//

#pragma once
#include <memory>
#include "IASTnode.h"

enum ucrType {Z, Y};
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

class ryNode final :public IASTnode {
public:
    explicit ryNode(double theta);
    ryNode();
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
class ucrNode:public IASTnode {
public:
    bool first=false;

    int name;
    explicit ucrNode(const std::vector<double>* angles, ucrType type, bool _first, bool _reverse);
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
class firstUcrNode:public ucrNode {
public:
    explicit firstUcrNode(const std::vector<double>* angles, ucrType type);
    void accept(nodeVisitor &visitor) override;
};