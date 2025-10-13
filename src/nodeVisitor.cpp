//
// Created by adenilton on 10/4/25.
//

#include "../include/nodeVisitor.h"
#include <utility>

qasmVisitor::qasmVisitor(int num_qubits) {
    _num_qubits = num_qubits;
    qasm_code += "OPENQASM 3.0;\n";
    qasm_code += "include \"stdgates.inc\";\n\n";

    qasm_code += "qubit[" + std::to_string(_num_qubits) + "] q; \n";
}

void qasmVisitor::visit(rzNode &node)
{
    int target = _num_qubits - 1;
    qasm_code += "rz(" + std::visit(return_type_visitor{}, node.get_data()) + ") q[" + std::to_string(target) + "];\n";
    // qasm_code += "rz(" + std::visit(return_type_visitor{}, node.get_data()) + ") q[0];\n";
}

void qasmVisitor::visit(ryNode &node)
{
    int target = _num_qubits - 1;
    qasm_code += "ry(" + std::visit(return_type_visitor{}, node.get_data()) + ") q[" + std::to_string(target) + "];\n";
    // qasm_code += "rz(" + std::visit(return_type_visitor{}, node.get_data()) + ") q[0];\n";
}

void qasmVisitor::visit(firstUcrNode &node) {
    int target = _num_qubits - 1;
    // int control = _num_qubits - node.get_num_qubits();
    int control = target - node.get_num_qubits() + 1;
    // int control = node.get_num_qubits() - 2;

    std::string ctrl_not = "cx q[" + std::to_string(control) +  "], q[" + std::to_string(target) + "];\n";
    // string ctrl_not = "cx q[" + std::to_string(node.get_num_qubits()-1) +  "], q[0];\n";
    node.gate1->accept(*this);
    qasm_code += ctrl_not;
    node.gate2->accept(*this);
    qasm_code += ctrl_not;
}


void qasmVisitor::visit(ucrNode &node) {
    int target = _num_qubits - 1;
    // int control = _num_qubits - node.get_num_qubits();
    int control = target - node.get_num_qubits() + 1;
    // int control = node.get_num_qubits() - 2;

    std::string ctrl_not = "cx q[" + std::to_string(control) +  "], q[" + std::to_string(target) + "];\n";
    // string ctrl_not = "cx q[" + std::to_string(node.get_num_qubits()-1) +  "], q[0];\n";
    if (node.reverse_gate) {
        node.gate1->reverse_gate = !node.gate1->reverse_gate;
        node.gate2->reverse_gate = !node.gate2->reverse_gate;
        std::swap(node.gate1, node.gate2);
    }
    node.gate1->accept(*this);
    qasm_code += ctrl_not;
    node.gate2->accept(*this);
    // qasm_code += ctrl_not; // TODO: Remover este cnot
}
