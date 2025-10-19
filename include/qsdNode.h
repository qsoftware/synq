//
// Created by adenilton on 10/11/25.
//

#ifndef TEST_QSD_H
#define TEST_QSD_H

#include "../include/IASTnode.h"
#include <Eigen/Dense>

struct qsd_result {
    std::vector<double> list_d;
    Eigen::MatrixXcf gate_v;
    Eigen::MatrixXcf gate_w;
};

class qsdNode: public IASTnode {
public:

    explicit qsdNode(Eigen::MatrixXcf gate1, Eigen::MatrixXcf gate2);
    void accept(nodeVisitor &visitor) override;
    return_type get_data() override;
    IASTnode* first_unitary;
    IASTnode* mcrz;
    IASTnode* second_unitary;
private:
    qsd_result _compute_gate(Eigen::MatrixXcf matrix, Eigen::MatrixXcf matrix_xcf);

};

#endif //TEST_QSD_H