//
// Created by alex on 22/05/26.
//

#pragma once

#include "iast.h"
#include <Eigen/Dense>
#include <memory>
#include <vector>

class ucrzNode;

/*
 * diagonal -> diagonal ucrz
 * diagonal -> one_qubit_diagonal
  * @ params diagonal unitary matrix to be decomposed
 */
class diagonalGateNode final : public IASTnode {
public:
    explicit diagonalGateNode(const Eigen::MatrixXcd& matrix);
    void accept(nodeVisitor &visitor) override;
    return_type get_data() override;

    static int index;
    std::unique_ptr<IASTnode> child_diagonal;
    std::unique_ptr<ucrzNode> ucrz;

    std::unique_ptr<IASTnode> createNextDiagonalGate(const Eigen::MatrixXcd& matrix);
    std::unique_ptr<ucrzNode> createUcrz(const std::vector<double>& angles);

private:
    Eigen::MatrixXcd data;
};

class oneQubitDiagonalGateNode final : public IASTnode {
public:
    explicit oneQubitDiagonalGateNode(double global_phase, double p_phase);
    void accept(nodeVisitor &visitor) override;
    return_type get_data() override;
    static std::unique_ptr<IASTnode> createDiagonalOneQubitGate(const Eigen::VectorXd& angles);

    double global_phase;
    double p_phase;
};
