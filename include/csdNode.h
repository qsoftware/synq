#pragma once
#include "iast.h"


class csdNode: public IASTnode {
public:
    std::unique_ptr<IASTnode> left_ucg;
    std::unique_ptr<IASTnode> mcry;
    std::unique_ptr<IASTnode> right_ucg;
    csdNode(Eigen::MatrixXcf& unitary_matrix);
    void accept(nodeVisitor &visitor) override;
    return_type get_data() override;

};

class qsdNode: public IASTnode {
public:
    std::unique_ptr<IASTnode> left_unitary;
    std::unique_ptr<IASTnode> mcrz;
    std::unique_ptr<IASTnode> right_unitary;
    qsdNode(Eigen::MatrixXcf unitary_matrix1, Eigen::MatrixXcf unitary_matrix2);
    void accept(nodeVisitor &visitor) override;
    return_type get_data() override;
};