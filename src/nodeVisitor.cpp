//
// Created by adenilton on 10/4/25.
//

#include "../include/nodeVisitor.h"
#include "../include/diagonalGateNode.h"
#include <utility>
#include<numbers>

qasmVisitor::qasmVisitor(int num_qubits) {
    _num_qubits = num_qubits;
    current_msb = _num_qubits - 1;
    qasm_code += "OPENQASM 3.0;\n";
    qasm_code += "include \"stdgates.inc\";\n\n";

    qasm_code += "qubit[" + std::to_string(_num_qubits) + "] q; \n";
}

void qasmVisitor::visit(rzNode &node)
{
    int target;

    if (node.position != -1) {
        target = node.position;
    } else {
        if (active_target != -1) {
            target = active_target;
        } else {
            target = current_msb;
        }
    }

    qasm_code += "rz(" + std::visit(return_type_visitor{}, node.get_data()) + ") q[" + std::to_string(target) + "];\n";
}

void qasmVisitor::visit(firstUcrzNode &node) {

    int target;
    int control;

    if(node.inverse){
        target = active_target;
        control = target + node.get_num_qubits() - 1;
    }else {
        target = current_msb;
        control = gate_synthesis ? (node.get_num_qubits() - 2) : (target - node.get_num_qubits() + 1);
    }

    std::string ctrl_not = "cx q[" + std::to_string(control) +  "], q[" + std::to_string(target) + "];\n";
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
        target = current_msb;
        control = gate_synthesis ? (node.get_num_qubits() - 2) : (target - node.get_num_qubits() + 1);
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

void qasmVisitor::visit(ryNode &node)
{
    int target;

    if (node.position != -1) {
        target = node.position;
    } else {
        if (active_target != -1) {
            target = active_target;
        } else {
            target = current_msb;
        }
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
        target = current_msb;
        control = gate_synthesis ? (node.get_num_qubits() - 2) : (target - node.get_num_qubits() + 1);
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
        target = current_msb;
        control = gate_synthesis ? (node.get_num_qubits() - 2) : (target - node.get_num_qubits() + 1);
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

void qasmVisitor::visit(unitaryGateNode &node) {

    gate_synthesis = true;

    int previous_msb = current_msb;

    current_msb = node.get_num_qubits() - 1;

    node.csd->accept(*this);

    if (node.get_num_qubits() == _num_qubits) {
        global_phase = std::fmod(global_phase, 2 * std::numbers::pi);
        qasm_code += "gphase(" + std::to_string(global_phase) + ");\n";
    }
    
    current_msb = previous_msb;
}

void qasmVisitor::visit(qsdNode &node) {
    if (node.gate_w) {
        node.gate_w->accept(*this);
    }
    if (node.ucrz) {
        node.ucrz->accept(*this);
    }
    if (node.gate_v) {
        node.gate_v->accept(*this);
    }
}

void qasmVisitor::visit(diagonalGateNode &node) {
    
    bool prev_gate_synthesis = gate_synthesis;
    gate_synthesis = true; 

    int previous_msb = current_msb;

    if (node.child_diagonal) {
        current_msb = previous_msb - 1; 
        node.child_diagonal->accept(*this);
    }
    
    if (node.ucrz) {
        current_msb = previous_msb; 
        node.ucrz->accept(*this);
    }
    
    current_msb = previous_msb;
    gate_synthesis = prev_gate_synthesis;
}

void qasmVisitor::visit(oneQubitDiagonalGateNode &node) {
    if (std::abs(node.global_phase) > 1e-12) {
        qasm_code += "gphase(" + std::to_string(node.global_phase) + ");\n";
    }
    qasm_code += "p(" + std::to_string(node.p_phase) + ") q[0];\n";
}

void qasmVisitor::visit(csdNode &node) {
    if (node.qsd2) {
        node.qsd2->accept(*this);
    }
    if (node.ucry) {
        node.ucry->accept(*this);
    }
    if (node.qsd1) {
        node.qsd1->accept(*this);
    }
}

void qasmVisitor::visit(unitaryOneQubitGateNode &node){
    global_phase += node.params.alpha;
    qasm_code += "rz(" + std::to_string(node.params.delta) + ") q[" + std::to_string(node.position) + "];\n";
    qasm_code += "ry(" + std::to_string(node.params.gamma) + ") q[" + std::to_string(node.position) + "];\n";
    qasm_code += "rz(" + std::to_string(node.params.beta) + ") q[" + std::to_string(node.position) + "];\n";
}

void qasmVisitor::visit(twoQubitGateNode &node) {
    if (std::abs(node.phase) > 1e-12) {
        global_phase += node.phase;
    }

    node.gate_c->accept(*this);
    node.gate_d->accept(*this);

    qasm_code += "cx q[0], q[1];\n";

    node.rz_theta->accept(*this);
    node.ry1_theta->accept(*this);

    qasm_code += "cx q[1], q[0];\n";

    node.ry2_theta->accept(*this);

    qasm_code += "cx q[0], q[1];\n";

    node.gate_a->accept(*this);
    node.gate_b->accept(*this);

}