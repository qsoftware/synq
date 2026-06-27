//
// Created by adenilton on 10/4/25.
//

#pragma once

#include <Eigen/Dense>
#include <variant>
#include <vector>
#include <complex>
#include <string>
#include <iomanip>
#include <sstream>

#include "ucrzNode.h"
#include "ucryNode.h"
#include "qspUcrNode.h"
#include "csdNode.h"
#include "qsdNode.h"
#include "unitaryGateNode.h"
class diagonalGateNode;
class oneQubitDiagonalGateNode;

class nodeVisitor {
public:
    virtual ~nodeVisitor() = default;
    virtual void visit(rzNode &node) = 0;
    virtual void visit(ucrzNode &node) = 0;
    virtual void visit(firstUcrzNode &node) = 0;
    virtual void visit(ryNode &node) = 0;
    virtual void visit(ucryNode &node) = 0;
    virtual void visit(firstUcryNode &node) = 0;
    virtual void visit(UCRotationNode &node) = 0;
    virtual void visit(qspUcrNode &node) = 0;
    virtual void visit(csdNode &node) = 0;
    virtual void visit(qsdNode &node) = 0;
    virtual void visit(unitaryGateNode &node) = 0;
    virtual void visit(twoQubitGateNode &node) = 0;
    virtual void visit(unitaryOneQubitGateNode &node) = 0;
    virtual void visit(diagonalGateNode &node) = 0;
    virtual void visit(oneQubitDiagonalGateNode &node) = 0;
};

struct return_type_visitor {
    std::string operator()(std::monostate) {
        return "";
    }

    std::string operator()(const double a) {
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(15) << a;
        return oss.str();
    }
    
    std::string operator()(const std::vector<double> a) { 
        return "";
    }
    
    std::string operator()(const std::vector<std::complex<double>>& a) {
        return ""; 
    }

    std::string operator()(const Eigen::MatrixXcd& a) {
        return "";
    }

    std::string operator()(const std::pair<Eigen::MatrixXcd, Eigen::MatrixXcd>& a) {
        return "";
    }
};

class qasmVisitor : public nodeVisitor {
public:
    std::string qasm_code;
    int _num_qubits;
    int current_msb;
    int active_target = -1;
    bool gate_synthesis = false;
    double global_phase = 0;
    explicit qasmVisitor(int num_qubits);
    void visit(rzNode &node) override;
    void visit(firstUcrzNode &node) override;
    void visit(ucrzNode &node) override;
    void visit(ryNode &node) override;
    void visit(firstUcryNode &node) override;
    void visit(ucryNode &node) override;
    void visit(UCRotationNode &node) override;
    void visit(qspUcrNode &node) override;
    void visit(csdNode &node) override;
    void visit(qsdNode &node) override;
    void visit(unitaryGateNode &node) override;
    void visit(twoQubitGateNode &node) override;
    void visit(unitaryOneQubitGateNode &node) override;
    void visit(diagonalGateNode &node) override;
    void visit(oneQubitDiagonalGateNode &node) override;

};
