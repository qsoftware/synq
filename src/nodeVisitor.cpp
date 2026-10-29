//
// Created by adenilton on 10/4/25.
//

#include "../include/nodeVisitor.h"
#include "../include/diagonalGateNode.h"
#include <utility>
#include<numbers>

qasmVisitor::qasmVisitor(int num_qubits, bool emit_qasm) : emit_qasm(emit_qasm) {
    _num_qubits = num_qubits;
    current_msb = _num_qubits - 1;
    if (emit_qasm) {
        qasm_code += "OPENQASM 3.0;\n";
        qasm_code += "include \"stdgates.inc\";\n\n";
        qasm_code += "qubit[" + std::to_string(_num_qubits) + "] q; \n";
    }
}

void qasmVisitor::emit_rotation(const gateType type, const int qubit, const double angle) {
    if (!emit_qasm) {
        gates.push_back({type, qubit, -1, angle});
        return;
    }
    static const char* names[] = {"rz", "ry", "p"};
    qasm_code += std::string(names[type]) + "(" + format_angle(angle) + ") q[" + std::to_string(qubit) + "];\n";
}

void qasmVisitor::emit_two_qubit(const gateType type, const int control, const int target) {
    if (!emit_qasm) {
        gates.push_back({type, control, target, 0.0});
        return;
    }
    qasm_code += std::string(type == CZ ? "cz" : "cx") + " q[" + std::to_string(control) + "], q[" + std::to_string(target) + "];\n";
}

void qasmVisitor::emit_h(const int qubit) {
    if (!emit_qasm) {
        gates.push_back({H, qubit, -1, 0.0});
        return;
    }
    qasm_code += "h q[" + std::to_string(qubit) + "];\n";
}

void qasmVisitor::visit(hNode &node) {
    emit_h(node.position);
}

void qasmVisitor::emit_gphase(const double angle) {
    if (!emit_qasm) {
        total_phase += angle;
        return;
    }
    qasm_code += "gphase(" + format_angle(angle) + ");\n";
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

    emit_rotation(RZ, target, std::get<double>(node.get_data()));
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

    // A "first" node starts a multiplexor: its children have no inherited reversal.
    const bool saved_context = reverse_context;
    const bool saved_mirror = mirror_context;
    reverse_context = false;
    mirror_context = node.mirror;
    const gateType two_qubit = node.use_cz ? CZ : CX;
    if (node.mirror) {
        // the same gates in reverse order (see UCRotationNode::mirror)
        if (!node.drop_first) {  // Block-ZXZ: the caller absorbed this gate
            emit_two_qubit(two_qubit, control, target);
        }
        node.gate2->accept(*this);
        emit_two_qubit(two_qubit, control, target);
        node.gate1->accept(*this);
    } else {
        node.gate1->accept(*this);
        emit_two_qubit(two_qubit, control, target);
        node.gate2->accept(*this);
        if (!node.drop_last) {  // Block-ZXZ: the caller absorbed this gate
            emit_two_qubit(two_qubit, control, target);
        }
    }
    reverse_context = saved_context;
    mirror_context = saved_mirror;
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

   
    // The node is reversed if its flag differs from the reversal inherited from its
    // parent; a reversed node visits gate2 before gate1. This used to be done by
    // swapping gate1/gate2 and toggling the children's flags in the tree, so a second
    // visit of the same tree produced a different (wrong) circuit. The reversal is now
    // visitor state passed down the recursion, and the tree is not modified.
    const bool reversed = node.reverse_gate != reverse_context;
    IASTnode& first = reversed ? *node.gate2 : *node.gate1;
    IASTnode& second = reversed ? *node.gate1 : *node.gate2;
    const bool saved_context = reverse_context;
    reverse_context = reversed;  // the children's flags are relative to this node
    if (mirror_context) {  // inside a mirrored multiplexor: the same gates in reverse order
        second.accept(*this);
        emit_two_qubit(node.use_cz ? CZ : CX, control, target);
        first.accept(*this);
    } else {
        first.accept(*this);
        emit_two_qubit(node.use_cz ? CZ : CX, control, target);
        second.accept(*this);
    }
    reverse_context = saved_context;
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
   
    emit_rotation(RY, target, std::get<double>(node.get_data()));
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


    // A "first" node starts a multiplexor: its children have no inherited reversal.
    const bool saved_context = reverse_context;
    reverse_context = false;
    node.gate1->accept(*this);
    emit_two_qubit(node.use_cz ? CZ : CX, control, target);
    node.gate2->accept(*this);
    if (!node.drop_last) {  // A.1: the caller absorbed this gate
        emit_two_qubit(node.use_cz ? CZ : CX, control, target);
    }
    reverse_context = saved_context;
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

    // The node is reversed if its flag differs from the reversal inherited from its
    // parent; a reversed node visits gate2 before gate1. This used to be done by
    // swapping gate1/gate2 and toggling the children's flags in the tree, so a second
    // visit of the same tree produced a different (wrong) circuit. The reversal is now
    // visitor state passed down the recursion, and the tree is not modified.
    const bool reversed = node.reverse_gate != reverse_context;
    IASTnode& first = reversed ? *node.gate2 : *node.gate1;
    IASTnode& second = reversed ? *node.gate1 : *node.gate2;
    const bool saved_context = reverse_context;
    reverse_context = reversed;  // the children's flags are relative to this node
    first.accept(*this);
    emit_two_qubit(node.use_cz ? CZ : CX, control, target);
    second.accept(*this);
    reverse_context = saved_context;
}

void qasmVisitor::visit(UCRotationNode &node) {
}

void qasmVisitor::visit(qspUcrNode &node) {
    if (std::abs(node.global_phase) > 1e-12) {
        emit_gphase(node.global_phase);
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
        emit_gphase(global_phase);
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

void qasmVisitor::visit(zxzNode &node) {
    node.gate_w_c->accept(*this);
    node.ucrz_c->accept(*this);
    node.h1->accept(*this);
    node.middle->accept(*this);
    node.h2->accept(*this);
    node.ucrz_a->accept(*this);
    node.gate_v_a->accept(*this);
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
        emit_gphase(node.global_phase);
    }
    emit_rotation(P, 0, node.p_phase);
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
    emit_rotation(RZ, node.position, node.params.delta);
    emit_rotation(RY, node.position, node.params.gamma);
    emit_rotation(RZ, node.position, node.params.beta);
}

void qasmVisitor::visit(twoQubitGateNode &node) {
    if (std::abs(node.phase) > 1e-12) {
        global_phase += node.phase;
    }

    node.gate_c->accept(*this);
    node.gate_d->accept(*this);

    if (node.two_cx) {
        // A.2: CX q0->q1; Rz (q1), Ry (q0); CX q0->q1
        emit_two_qubit(CX, 0, 1);
        node.rz_theta->accept(*this);
        node.ry1_theta->accept(*this);
        emit_two_qubit(CX, 0, 1);
    } else {
        emit_two_qubit(CX, 0, 1);
        node.rz_theta->accept(*this);
        node.ry1_theta->accept(*this);
        emit_two_qubit(CX, 1, 0);
        node.ry2_theta->accept(*this);
        emit_two_qubit(CX, 0, 1);
    }

    node.gate_a->accept(*this);
    node.gate_b->accept(*this);
}
