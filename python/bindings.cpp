/**
 * @file bindings.cpp
 * @brief Python bindings for SynQ (module ``synq._synq``).
 *
 * Each function builds the corresponding AST node and runs the qasmVisitor.
 * state_preparation / unitary / diagonal return the OpenQASM 3 program as a
 * string; the *_gates variants return the same circuit as numpy arrays.
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include <cmath>
#include <complex>
#include <stdexcept>
#include <string>
#include <vector>
#include <cstdint>

#include "../include/diagonalGateNode.h"
#include "../include/nodeVisitor.h"
#include "../include/qspUcrNode.h"
#include "../include/unitaryGateNode.h"

namespace py = pybind11;

namespace {

int num_qubits_from_size(std::size_t size, const char* what) {
    if (size < 2 || (size & (size - 1)) != 0) {
        throw std::invalid_argument(std::string(what) + " length must be a power of 2 (>= 2)");
    }
    return static_cast<int>(std::lround(std::log2(static_cast<double>(size))));
}

/**
 * Converts the gates recorded by a qasmVisitor into numpy arrays:
 * (types uint8[N], qubits int32[N, 2], angles float64[N], global_phase float).
 * Gate types: 0 rz, 1 ry, 2 p, 3 cx. Single-qubit gates have qubits[k, 1] == -1;
 * two-qubit gates store (control, target).
 */
py::tuple to_arrays(const qasmVisitor& visitor) {
    const auto n = static_cast<py::ssize_t>(visitor.gates.size());
    py::array_t<std::uint8_t> types(n);
    py::array_t<std::int32_t> qubits({n, static_cast<py::ssize_t>(2)});
    py::array_t<double> angles(n);
    auto t = types.mutable_unchecked<1>();
    auto q = qubits.mutable_unchecked<2>();
    auto a = angles.mutable_unchecked<1>();
    for (py::ssize_t k = 0; k < n; ++k) {
        const auto& gate = visitor.gates[static_cast<std::size_t>(k)];
        t(k) = gate.type;
        q(k, 0) = gate.q0;
        q(k, 1) = gate.q1;
        a(k) = gate.angle;
    }
    return py::make_tuple(types, qubits, angles, visitor.total_phase);
}

int check_state(const std::vector<std::complex<double>>& state) {
    const int n = num_qubits_from_size(state.size(), "state");
    double norm2 = 0.0;
    for (const auto& a : state) norm2 += std::norm(a);
    if (std::abs(norm2 - 1.0) > 1e-8) {
        throw std::invalid_argument("state must be normalized");
    }
    return n;
}

Eigen::MatrixXcd diagonal_matrix(const Eigen::VectorXcd& diag) {
    num_qubits_from_size(static_cast<std::size_t>(diag.size()), "diagonal");
    for (Eigen::Index i = 0; i < diag.size(); ++i) {
        if (std::abs(std::abs(diag(i)) - 1.0) > 1e-8) {
            throw std::invalid_argument("diagonal entries must have modulus 1");
        }
    }
    return diag.asDiagonal();
}

std::string state_preparation(const std::vector<std::complex<double>>& state) {
    const int n = check_state(state);
    qspUcrNode node(&state);
    qasmVisitor visitor(n);
    node.accept(visitor);
    return visitor.qasm_code;
}

py::tuple state_preparation_gates(const std::vector<std::complex<double>>& state) {
    const int n = check_state(state);
    qspUcrNode node(&state);
    qasmVisitor visitor(n, false);
    node.accept(visitor);
    return to_arrays(visitor);
}

std::string unitary(const Eigen::MatrixXcd& matrix) {
    const int n = num_qubits_from_size(static_cast<std::size_t>(matrix.rows()), "matrix");
    unitaryGateNode node(matrix);  // validates shape and unitarity
    qasmVisitor visitor(n);
    node.accept(visitor);
    return visitor.qasm_code;
}

py::tuple unitary_gates(const Eigen::MatrixXcd& matrix) {
    const int n = num_qubits_from_size(static_cast<std::size_t>(matrix.rows()), "matrix");
    unitaryGateNode node(matrix);  // validates shape and unitarity
    qasmVisitor visitor(n, false);
    node.accept(visitor);
    return to_arrays(visitor);
}

std::string diagonal(const Eigen::VectorXcd& diag) {
    const Eigen::MatrixXcd matrix = diagonal_matrix(diag);
    const int n = num_qubits_from_size(static_cast<std::size_t>(diag.size()), "diagonal");
    diagonalGateNode node(matrix);
    qasmVisitor visitor(n);
    node.accept(visitor);
    return visitor.qasm_code;
}

py::tuple diagonal_gates(const Eigen::VectorXcd& diag) {
    const Eigen::MatrixXcd matrix = diagonal_matrix(diag);
    const int n = num_qubits_from_size(static_cast<std::size_t>(diag.size()), "diagonal");
    diagonalGateNode node(matrix);
    qasmVisitor visitor(n, false);
    node.accept(visitor);
    return to_arrays(visitor);
}

}  // namespace

PYBIND11_MODULE(_synq, m) {
    m.doc() = "SynQ: quantum circuit synthesis (C++ core). Functions return OpenQASM 3.";

    m.def("state_preparation", &state_preparation, py::arg("state"),
          "Circuit (OpenQASM 3) that prepares the normalized state vector `state` from |0...0>.");
    m.def("unitary", &unitary, py::arg("matrix"),
          "Circuit (OpenQASM 3) implementing the 2^n x 2^n unitary `matrix` (CSD/QSD).");
    m.def("diagonal", &diagonal, py::arg("diagonal"),
          "Circuit (OpenQASM 3) implementing diag(`diagonal`); entries must have modulus 1.");

    const char* gates_doc =
        "Same as the QASM function, but returns the gates as numpy arrays: "
        "(types uint8[N], qubits int32[N, 2], angles float64[N], global_phase). "
        "Types: 0 rz, 1 ry, 2 p, 3 cx; single-qubit gates have qubits[k, 1] == -1, "
        "two-qubit gates store (control, target).";
    m.def("state_preparation_gates", &state_preparation_gates, py::arg("state"), gates_doc);
    m.def("unitary_gates", &unitary_gates, py::arg("matrix"), gates_doc);
    m.def("diagonal_gates", &diagonal_gates, py::arg("diagonal"), gates_doc);

#ifdef SYNQ_VERSION
    m.attr("__version__") = SYNQ_VERSION;
#endif
}
