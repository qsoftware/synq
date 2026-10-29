//
// Created by adenilton on 10/4/25.
//

#pragma once

#include <Eigen/Dense>
#include <variant>
#include <vector>
#include <complex>
#include <string>
#include <sstream>
#include <iomanip>
#include <limits>
#include <locale>
#include <cstdint>

#include "ucrzNode.h"
#include "ucryNode.h"
#include "qspUcrNode.h"
#include "csdNode.h"
#include "qsdNode.h"
#include "unitaryGateNode.h"
#include "hNode.h"
#include "zxzNode.h"
class diagonalGateNode;
class oneQubitDiagonalGateNode;

class nodeVisitor {
public:
    virtual ~nodeVisitor() = default;
    virtual void visit(rzNode &node) = 0;
    virtual void visit(ucrzNode &node) = 0;
    virtual void visit(firstUcrzNode &node) = 0;
    virtual void visit(ryNode &node) = 0;
    virtual void visit(hNode &node) = 0;
    virtual void visit(ucryNode &node) = 0;
    virtual void visit(firstUcryNode &node) = 0;
    virtual void visit(UCRotationNode &node) = 0;
    virtual void visit(qspUcrNode &node) = 0;
    virtual void visit(csdNode &node) = 0;
    virtual void visit(qsdNode &node) = 0;
    virtual void visit(zxzNode &node) = 0;
    virtual void visit(unitaryGateNode &node) = 0;
    virtual void visit(twoQubitGateNode &node) = 0;
    virtual void visit(unitaryOneQubitGateNode &node) = 0;
    virtual void visit(diagonalGateNode &node) = 0;
    virtual void visit(oneQubitDiagonalGateNode &node) = 0;
};

/**
 * @brief Formats an angle with full double precision (round-trip safe).
 *
 * std::to_string only keeps 6 decimal places, which limits the accuracy of the
 * generated circuits to ~1e-6.
 */
inline std::string format_angle(const double value) {
    std::ostringstream stream;
    stream.imbue(std::locale::classic());
    stream << std::setprecision(std::numeric_limits<double>::max_digits10) << value;
    return stream.str();
}

struct return_type_visitor {
    std::string operator()(std::monostate) {
        return "";
    }

    std::string operator()(const double a) {
        return format_angle(a);
    }
    
    std::string operator()(const std::vector<double> a) { 
        return "";
    }
    
    std::string operator()(const std::vector<std::complex<double>>& a) {
        return ""; 
    }

    std::string operator()(const Eigen::MatrixXcd& a) {
        return "";
    }

    std::string operator()(const std::pair<Eigen::MatrixXcd, Eigen::MatrixXcd>& a) {
        return "";
    }
};

/**
 * @brief A gate of the synthesized circuit (used when the visitor records gates).
 *
 * type: 0 = rz, 1 = ry, 2 = p (q0 = target, q1 = -1, angle),
 *       3 = cx, 4 = cz (q0 = control, q1 = target),
 *       5 = h (q0 = target, q1 = -1, angle = 0).
 */
struct circuitGate {
    std::uint8_t type;
    int q0;
    int q1;
    double angle;
};

/**
 * @brief Emits the circuit of an AST, as OpenQASM 3 text (default) or as a list of gates.
 *
 * Both outputs describe the same circuit: every gate goes through emit_rotation,
 * emit_two_qubit or emit_gphase, which either append QASM text or record the gate.
 */
class qasmVisitor : public nodeVisitor {
public:
    enum gateType : std::uint8_t { RZ = 0, RY = 1, P = 2, CX = 3, CZ = 4, H = 5 };

    std::string qasm_code;
    int _num_qubits;
    int current_msb;
    int active_target = -1;
    bool gate_synthesis = false;
    double global_phase = 0;

    /** If false, gates are recorded in `gates` (and `total_phase`) instead of `qasm_code`. */
    bool emit_qasm = true;
    /** Gates in time order (only when emit_qasm is false). */
    std::vector<circuitGate> gates;
    /** Sum of all global phases (only when emit_qasm is false). */
    double total_phase = 0.0;

    explicit qasmVisitor(int num_qubits, bool emit_qasm = true);
    void visit(rzNode &node) override;
    void visit(firstUcrzNode &node) override;
    void visit(ucrzNode &node) override;
    void visit(ryNode &node) override;
    void visit(hNode &node) override;
    void visit(firstUcryNode &node) override;
    void visit(ucryNode &node) override;
    void visit(UCRotationNode &node) override;
    void visit(qspUcrNode &node) override;
    void visit(csdNode &node) override;
    void visit(qsdNode &node) override;
    void visit(zxzNode &node) override;
    void visit(unitaryGateNode &node) override;
    void visit(twoQubitGateNode &node) override;
    void visit(unitaryOneQubitGateNode &node) override;
    void visit(diagonalGateNode &node) override;
    void visit(oneQubitDiagonalGateNode &node) override;

private:
    /** Reversal inherited from the parent UCR node (replaces modifying the tree). */
    bool reverse_context = false;
    /** Inside a mirrored RZ multiplexor: every RZ node emits its gates in reverse order. */
    bool mirror_context = false;

    void emit_rotation(gateType type, int qubit, double angle);
    void emit_two_qubit(gateType type, int control, int target);
    void emit_h(int qubit);
    void emit_gphase(double angle);
};
