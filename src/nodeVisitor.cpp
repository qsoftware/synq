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
    int target;
    if (active_target != -1) {
        target = active_target;
    } else {
        target = _num_qubits - 1;
    }
   
    qasm_code += "rz(" + std::visit(return_type_visitor{}, node.get_data()) + ") q[" + std::to_string(target) + "];\n";
    // qasm_code += "rz(" + std::visit(return_type_visitor{}, node.get_data()) + ") q[0];\n";
}

void qasmVisitor::visit(firstUcrzNode &node) {

    int target;
    int control;

    if(node.inverse){
        target = active_target;
        control = target + node.get_num_qubits() - 1;
    }else {
        target = _num_qubits - 1;
        // int control = _num_qubits - node.get_num_qubits();
        control = target - node.get_num_qubits() + 1;
        // int control = node.get_num_qubits() - 2;
    }

    std::string ctrl_not = "cx q[" + std::to_string(control) +  "], q[" + std::to_string(target) + "];\n";
    // string ctrl_not = "cx q[" + std::to_string(node.get_num_qubits()-1) +  "], q[0];\n";
    node.gate1->accept(*this);
    qasm_code += ctrl_not;
    node.gate2->accept(*this);
    qasm_code += ctrl_not;
}


void qasmVisitor::visit(ucrzNode &node) {

    int target;
    int control;

    if(node.inverse){
        target = active_target;
        control = target + node.get_num_qubits() - 1;
    }else {
        target = _num_qubits - 1;
        // int control = _num_qubits - node.get_num_qubits();
        control = target - node.get_num_qubits() + 1;
        // int control = node.get_num_qubits() - 2;
    }

   

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

void qasmVisitor::visit(ryNode &node)
{
    int target;
    if (active_target != -1) {
        target = active_target;
    } else {
        target = _num_qubits - 1;
    }
   
    qasm_code += "ry(" + std::visit(return_type_visitor{}, node.get_data()) + ") q[" + std::to_string(target) + "];\n";
}

void qasmVisitor::visit(firstUcryNode &node) {
    int control;
    int target;

    if(node.inverse){

        target = _num_qubits - node.get_num_qubits();

        if (active_target < target && active_target != -1) target = active_target;
        else active_target = target;
        control = target + node.get_num_qubits() - 1;
        
    }else {
        target = _num_qubits - 1;
        control = target - node.get_num_qubits() + 1;
    }


    std::string ctrl_not = "cx q[" + std::to_string(control) +  "], q[" + std::to_string(target) + "];\n";
    node.gate1->accept(*this);
    qasm_code += ctrl_not;
    node.gate2->accept(*this);
    qasm_code += ctrl_not;
}

void qasmVisitor::visit(ucryNode &node) {

    int control;
    int target;

    if(node.inverse){

        target = _num_qubits - node.get_num_qubits();

        if (active_target < target && active_target != -1) target = active_target;
        else active_target = target;
        control = target + node.get_num_qubits() - 1;
        
    }else {
        target = _num_qubits - 1;
        control = target - node.get_num_qubits() + 1;
    }

    std::string ctrl_not = "cx q[" + std::to_string(control) +  "], q[" + std::to_string(target) + "];\n";
    if (node.reverse_gate) {
        node.gate1->reverse_gate = !node.gate1->reverse_gate;
        node.gate2->reverse_gate = !node.gate2->reverse_gate;
        std::swap(node.gate1, node.gate2);
    }
    node.gate1->accept(*this);
    qasm_code += ctrl_not;
    node.gate2->accept(*this);
}

void qasmVisitor::visit(UCRotationNode &node) {
}

void qasmVisitor::visit(qspUcrNode &node) {
    if (std::abs(node.global_phase) > 1e-12) {
        qasm_code += "gphase(" + std::to_string(node.global_phase) + ");\n";
    }
    if (node.base_ry) {
        node.base_ry->accept(*this);
    }
    if (node.base_rz) {
        node.base_rz->accept(*this);
    }
    if (node.next_qsp) {
        node.next_qsp->accept(*this);
    }
    if (node.ucry) {
        node.ucry->accept(*this);
    }
    if (node.ucrz) {
        node.ucrz->accept(*this);
    }
}

void qasmVisitor::visit(unitaryNode &node) {
    if (node.num_qubits == 1) {
        qasm_code += std::visit(return_type_visitor{}, node.get_data());
    } else
        node.decomposition->accept(*this);
}

void qasmVisitor::visit(csdNode &node) {
    if (node.left_ucg)
        node.left_ucg->accept(*this);

    if (node.mcry)
        node.mcry->accept(*this);

    if (node.right_ucg)
        node.right_ucg->accept(*this);
}

void qasmVisitor::visit(qsdNode &node) {

    if (node.right_unitary)
        node.right_unitary->accept(*this);

    if (node.mcrz)
        node.mcrz->accept(*this);

    if (node.left_unitary)
        node.left_unitary->accept(*this);
}
