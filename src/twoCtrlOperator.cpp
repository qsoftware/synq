#include "../include/twoCtrlOperator.h"
#include "../include/one_qubit_gate.h"
#include <cmath>
#include <vector>
#include <iostream>
#include <unsupported/Eigen/MatrixFunctions>
#include "../include/nodeVisitor.h"
using namespace std;

// Create an instance of CtrlOperatorNode
CtrlOperatorNode::CtrlOperatorNode(std::vector<int>& ctrl_bits, int target_bit, Eigen::MatrixXcd& matrix, Eigen::Matrix2cd& Op) : controls(ctrl_bits), target(target_bit), CtrlMatrix(matrix), OperatorMatrix(Op) {}

// Returns the position of the specified control bit
int CtrlOperatorNode::get_ctrl(int index) {return controls.at(index);}

// Returns the number of control bits
int CtrlOperatorNode::get_num_ctrl() {return controls.size();}

// Returns the position of the target bit
int CtrlOperatorNode::get_target() {return target;}

// Makes visitor accept this node
void CtrlOperatorNode::accept(nodeVisitor &visitor) {
    visitor.visit(*this);
}

// Returns the control matrix stored in the node
return_type CtrlOperatorNode::get_data() {
    return CtrlMatrix;
}

// Factory method to create a unique pointer to an CtrlOperatorNode
static std::unique_ptr<IASTnode> createCtrlOperatorNode(std::vector<int>& controls, int target, Eigen::MatrixXcd& CtrlMatrix, Eigen::Matrix2cd& OperatorMatrix) {
    return std::make_unique<CtrlOperatorNode>(controls, target, CtrlMatrix, OperatorMatrix);
}

// Applies a operator controlled by one qubit on a target bit
Eigen::MatrixXcd CtrlOperatorNode::oneCtrlOperator(const Eigen::Matrix2cd& U, int ctrl, int targ) {

    int dim = CtrlMatrix.rows(); // 2^num_qubits
    Eigen::MatrixXcd result = Eigen::MatrixXcd::Identity(dim, dim); // Initialize result as the identity matrix of size 2^num_qubits

    // Sweep through all possible states of the system 
    for (int col = 0; col < dim; col++)
    {
        // Extract the bit representing the control qubit from the current column index
        // Example: col = 1101 and ctrl = 2 (from left to right), then (col >> ctrl) = 0011, and (col >> ctrl) & 0001 = 1
        bool c = (col >> ctrl) & 1;

        // If the control qubit is not set, skip the column 
        if (!c) continue;

        // Identify the index of the states where the target qubit is |0> and |1>
        int base0 = col & ~(1 << targ); //base0 is the current col with the target bit set to 0
        int base1 = base0 | (1 << targ); //base1 is the current col with the target bit set to 1

        // Map the 2x2 matrix U into the global matrix
        for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
        {
            int row = (i == 0) ? base0 : base1;
            int c_idx = (j == 0) ? base0 : base1;

            result(row, c_idx) = U(i, j);
        }
    }

    if(get_num_ctrl() == 1){
    // Updates the control matrix with the result of the operation if there is only one control qubit
    CtrlMatrix = result * CtrlMatrix;
    return CtrlMatrix;}
    return result;
}

/*
Applies a operator controlled by two qubits on a target bit.
The function uses the decomposition of a two-controlled gate into a sequence of one-controlled gates and CNOTs, as follows:

  q1 (ctrl1)    ───●───       ───────────●───────────●───●───
                   │                     │           │   │
  q2 (ctrl2)    ───●───   =   ───●───────⊕──●───────⊕──│───
                   │             │           │           │
  q3 (target)   ──[U]──       ──[V]─────────[V†]────────[V]──
*/
Eigen::MatrixXcd CtrlOperatorNode::twoCtrlOperator(const Eigen::Matrix2cd& U, int ctrl1, int ctrl2) {
    // V is the square root of U, such that V² = U
    Eigen::Matrix2cd V = U.sqrt();
    // V† is the adjoint of V
    Eigen::Matrix2cd V_dagger = V.adjoint();
    // X is the Pauli-X gate (CNOT)
    Eigen::Matrix2cd X = OneQubit::x_matrix();

    // Applies the V operator on target controlled by ctrl2
    Eigen::MatrixXcd op1 = oneCtrlOperator(V, ctrl2, target);
    
    // Applies the CNOT gate on ctrl2 controlled by ctrl1
    Eigen::MatrixXcd op2 = oneCtrlOperator(X, ctrl1, ctrl2);
    
    // Applies V† controlled by ctrl2 on target
    Eigen::MatrixXcd op3 = oneCtrlOperator(V_dagger, ctrl2, target);
    
    // Applies the CNOT gate on ctrl2 controlled by ctrl1
    Eigen::MatrixXcd op4 = oneCtrlOperator(X, ctrl1, ctrl2);
    
    // Applies V controlled by ctrl1 on target
    Eigen::MatrixXcd op5 = oneCtrlOperator(V, ctrl1, target);

    Eigen::MatrixXcd result = op5 * op4 * op3 * op2 * op1;

    // Updates the control matrix with the result of the operation
    CtrlMatrix = result * CtrlMatrix;
    return CtrlMatrix;
}
