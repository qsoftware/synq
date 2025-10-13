//
// Created by adenilton on 10/4/25.
//

#pragma once

#include "ucrNode.h"


class nodeVisitor {
public:
    virtual ~nodeVisitor() = default;
    virtual void visit(rzNode &node) = 0;
    virtual void visit(ryNode &node) = 0;
    virtual void visit(ucrNode &node) = 0;
    virtual void visit(firstUcrNode &node) = 0;
};

struct return_type_visitor {
  std::string operator()(const double a) {
      return std::to_string(a);
  }
    std::string operator()(const std::vector<double> a) {
      return "";
  }
};

class qasmVisitor : public nodeVisitor {
public:
    std::string qasm_code;
    int _num_qubits;
    explicit qasmVisitor(int num_qubits);
    void visit(rzNode &node) override;
    void visit(ryNode &node) override;
    void visit(firstUcrNode &node) override;
    void visit(ucrNode &node) override;
};
