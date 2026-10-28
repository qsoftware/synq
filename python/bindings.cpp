/**
 * @file bindings.cpp
 * @brief Python bindings for SynQ (module ``synq._synq``).
 *
 * Each function builds the corresponding AST node, runs the qasmVisitor and
 * returns the generated OpenQASM 3 program as a string.
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/eigen.h>

#include <cmath>
#include <complex>
#include <stdexcept>
#include <string>
#include <vector>

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

std::string state_preparation(const std::vector<std::complex<double>>& state) {
    const int n = num_qubits_from_size(state.size(), "state");
    double norm2 = 0.0;
    for (const auto& a : state) norm2 += std::norm(a);
    if (std::abs(norm2 - 1.0) > 1e-8) {
        throw std::invalid_argument("state must be normalized");
    }
    qspUcrNode node(&state);
    qasmVisitor visitor(n);
    node.accept(visitor);
    return visitor.qasm_code;
}

std::string unitary(const Eigen::MatrixXcd& matrix) {
    const int n = num_qubits_from_size(static_cast<std::size_t>(matrix.rows()), "matrix");
    unitaryGateNode node(matrix);  // validates shape and unitarity
    qasmVisitor visitor(n);
    node.accept(visitor);
    return visitor.qasm_code;
}

std::string diagonal(const Eigen::VectorXcd& diag) {
    const int n = num_qubits_from_size(static_cast<std::size_t>(diag.size()), "diagonal");
    for (Eigen::Index i = 0; i < diag.size(); ++i) {
        if (std::abs(std::abs(diag(i)) - 1.0) > 1e-8) {
            throw std::invalid_argument("diagonal entries must have modulus 1");
        }
    }
    Eigen::MatrixXcd matrix = diag.asDiagonal();
    diagonalGateNode node(matrix);
    qasmVisitor visitor(n);
    node.accept(visitor);
    return visitor.qasm_code;
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

#ifdef SYNQ_VERSION
    m.attr("__version__") = SYNQ_VERSION;
#endif
}
