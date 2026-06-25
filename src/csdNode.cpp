//
// Created by alex on 07/04/26.
//

#include "../include/csdNode.h"
#include "../include/csd.h"
#include "../include/nodeVisitor.h"
#include "../include/qsdNode.h"
#include "../include/ucryNode.h"    

csdNode::csdNode(const Eigen::MatrixXcd& matrix) : m_matrix(matrix) {
    auto csdValues = csd(matrix, matrix.rows() / 2, matrix.cols() / 2);

    auto left_gate = std::vector<Eigen::MatrixXcd>{csdValues.V1T, csdValues.V2T};
    auto right_gate = std::vector<Eigen::MatrixXcd>{csdValues.U1, csdValues.U2};

    Eigen::VectorXd theta_eigen = Eigen::Map<Eigen::VectorXd>(csdValues.theta.data(), csdValues.theta.size());
    Eigen::VectorXd angles_eigen = 2.0 * theta_eigen;


    if (is_generalized_permutation(left_gate[0]) && is_generalized_permutation(left_gate[1])) {
        std::tie(left_gate, right_gate, angles_eigen) = permutation_optimization(left_gate, right_gate, angles_eigen, 'L');        
    }
    
    if (is_generalized_permutation(right_gate[0]) && is_generalized_permutation(right_gate[1])) {
        std::tie(left_gate, right_gate, angles_eigen) = permutation_optimization(left_gate, right_gate, angles_eigen, 'R');        
    }

    qsd1 = createQSD(left_gate[0], left_gate[1]);
    qsd2 = createQSD(right_gate[0], right_gate[1]);

    std::vector<double> angles(angles_eigen.data(), angles_eigen.data() + angles_eigen.size());
    ucry = createUcry(angles);
}

std::tuple<std::vector<Eigen::MatrixXcd>, std::vector<Eigen::MatrixXcd>, Eigen::VectorXd>
csdNode::permutation_optimization(const std::vector<Eigen::MatrixXcd>& left_gate, const std::vector<Eigen::MatrixXcd>& right_gate, const Eigen::VectorXd& theta, char side) {
    
    std::vector<Eigen::MatrixXcd> new_left_gate = left_gate;
    std::vector<Eigen::MatrixXcd> new_right_gate = right_gate;
    Eigen::VectorXd new_theta = theta;

    double tol = 1e-10;

    if (side == 'R') {

        Eigen::MatrixXd permutation = (right_gate[0].array().abs() > tol).cast<double>().matrix();
        
        new_right_gate[0] = permutation.transpose() * right_gate[0];
        new_right_gate[1] = permutation.transpose() * right_gate[1];

        new_theta = permutation.transpose() * theta;

        new_left_gate[0] = left_gate[0] * permutation;
        new_left_gate[1] = left_gate[1] * permutation;

    } else if (side == 'L') {
        
        Eigen::MatrixXd permutation = (left_gate[0].array().abs() > tol).cast<double>().matrix();
        
        new_left_gate[0] = left_gate[0] * permutation.transpose();
        new_left_gate[1] = left_gate[1] * permutation.transpose();

        new_theta = permutation * theta;

        new_right_gate[0] = permutation * right_gate[0];
        new_right_gate[1] = permutation * right_gate[1];
    }

    return std::make_tuple(new_left_gate, new_right_gate, new_theta);
}

bool csdNode::is_generalized_permutation(const Eigen::MatrixXcd& matrix) {
    double tol = 1e-10;
    Eigen::VectorXi row_count = (matrix.array().abs() > tol).cast<int>().rowwise().sum();
    Eigen::RowVectorXi col_count = (matrix.array().abs() > tol).cast<int>().colwise().sum();
    return (row_count.array() == 1).all() && (col_count.array() == 1).all();
}

std::unique_ptr<qsdNode> csdNode::createQSD(const Eigen::MatrixXcd& gate1, const Eigen::MatrixXcd& gate2) {
    return std::make_unique<qsdNode>(gate1, gate2);
}

std::unique_ptr<ucryNode> csdNode::createUcry(const std::vector<double>& angles) {
    return std::make_unique<firstUcryNode>(&angles, false);
}

void csdNode::accept(nodeVisitor &visitor) {
    visitor.visit(*this);
}

return_type csdNode::get_data() {
    return m_matrix;
}