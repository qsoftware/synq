//
// Created by adenilton on 10/4/25.
//

#include "../include/ctrl_qubit_gate.h"
#include "../include/nodeVisitor.h"

ctrl_qubit_gate::ctrl_qubit_gate(int control_qubit, int target_qubit, Eigen::Matrix2cd uMatrix) {
    control = control_qubit;
    target = target_qubit;
    data = uMatrix;
    num_qubits = 2; 
    
    decompose();
}

void ctrl_qubit_gate::decompose() {
   
    angles = OneQubit::zyz_decomposition(data);
    phase_alpha = angles.alpha;

    // Matriz C
    matrix_C = OneQubit::rz_matrix((angles.delta - angles.beta) / 2.0);

    // Matriz B 
    Eigen::Matrix2cd ry_B = OneQubit::ry_matrix(-angles.gamma / 2.0);
    Eigen::Matrix2cd rz_B = OneQubit::rz_matrix(-(angles.delta + angles.beta) / 2.0);
    matrix_B = ry_B * rz_B;

    // Matriz A     
    Eigen::Matrix2cd rz_A = OneQubit::rz_matrix(angles.beta);
    Eigen::Matrix2cd ry_A = OneQubit::ry_matrix(angles.gamma / 2.0);
    matrix_A = rz_A * ry_A;
}

void ctrl_qubit_gate::accept(nodeVisitor &visitor) {
    visitor.visit(*this);
}

return_type ctrl_qubit_gate::get_data() {
    return data;
}
