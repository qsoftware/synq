//
// Correctness tests: the generated OpenQASM is simulated and compared with the
// target state / unitary, including the global phase. Qubit order is
// little-endian (qubit 0 is the least significant bit), as in Qiskit.
//

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <numbers>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <Eigen/Dense>

#include "../include/diagonalGateNode.h"
#include "../include/nodeVisitor.h"
#include "../include/qspUcrNode.h"
#include "../include/unitaryGateNode.h"

namespace {

using cd = std::complex<double>;
constexpr double pi = std::numbers::pi;

// ---------------------------------------------------------------------------
// Minimal simulator for the OpenQASM subset emitted by SynQ:
// rz, ry, p, cx, cz and gphase on a single register q.
// ---------------------------------------------------------------------------

void apply_1q(Eigen::VectorXcd& state, const Eigen::Matrix2cd& gate, int qubit) {
    const Eigen::Index stride = Eigen::Index(1) << qubit;
    for (Eigen::Index i = 0; i < state.size(); ++i) {
        if (i & stride) continue;
        const cd a = state(i), b = state(i | stride);
        state(i) = gate(0, 0) * a + gate(0, 1) * b;
        state(i | stride) = gate(1, 0) * a + gate(1, 1) * b;
    }
}

void apply_2q(Eigen::VectorXcd& state, bool is_cz, int control, int target) {
    const Eigen::Index c = Eigen::Index(1) << control, t = Eigen::Index(1) << target;
    for (Eigen::Index i = 0; i < state.size(); ++i) {
        if (!(i & c)) continue;
        if (is_cz) {
            if (i & t) state(i) = -state(i);
        } else if (!(i & t)) {
            std::swap(state(i), state(i | t));
        }
    }
}

int parse_index(const std::string& text, std::size_t from) {
    const std::size_t open = text.find('[', from);
    return std::stoi(text.substr(open + 1, text.find(']', open) - open - 1));
}

double parse_angle(const std::string& line) {
    const std::size_t open = line.find('(');
    return std::stod(line.substr(open + 1, line.find(')') - open - 1));
}

Eigen::VectorXcd simulate(const std::string& qasm, Eigen::VectorXcd state) {
    std::istringstream stream(qasm);
    std::string line;
    while (std::getline(stream, line)) {
        if (line.empty() || line.rfind("OPENQASM", 0) == 0 || line.rfind("include", 0) == 0 ||
            line.rfind("qubit", 0) == 0) {
            continue;
        }
        if (line.rfind("gphase", 0) == 0) {
            state *= std::exp(cd(0.0, parse_angle(line)));
        } else if (line.rfind("rz", 0) == 0 || line.rfind("ry", 0) == 0 || line.rfind("p(", 0) == 0) {
            const double t = parse_angle(line);
            const int q = parse_index(line, line.find(')'));
            Eigen::Matrix2cd g;
            if (line[0] == 'p') {
                g << 1.0, 0.0, 0.0, std::exp(cd(0.0, t));
            } else if (line[1] == 'z') {
                g << std::exp(cd(0.0, -t / 2)), 0.0, 0.0, std::exp(cd(0.0, t / 2));
            } else {
                g << std::cos(t / 2), -std::sin(t / 2), std::sin(t / 2), std::cos(t / 2);
            }
            apply_1q(state, g, q);
        } else if (line.rfind("cx", 0) == 0 || line.rfind("cz", 0) == 0) {
            const int control = parse_index(line, 0);
            const int target = parse_index(line, line.find(','));
            apply_2q(state, line[1] == 'z', control, target);
        } else {
            throw std::runtime_error("unexpected QASM line: " + line);
        }
    }
    return state;
}

Eigen::MatrixXcd operator_of(const std::string& qasm, Eigen::Index dim) {
    Eigen::MatrixXcd result(dim, dim);
    for (Eigen::Index k = 0; k < dim; ++k) {
        result.col(k) = simulate(qasm, Eigen::VectorXcd::Unit(dim, k));
    }
    return result;
}

// ---------------------------------------------------------------------------
// Synthesis helpers and deterministic test data
// ---------------------------------------------------------------------------

int qubits_of(Eigen::Index dim) { return static_cast<int>(std::lround(std::log2(double(dim)))); }

std::string synthesize_state(const Eigen::VectorXcd& v) {
    std::vector<cd> state(v.data(), v.data() + v.size());
    qspUcrNode node(&state);
    qasmVisitor visitor(qubits_of(v.size()));
    node.accept(visitor);
    return visitor.qasm_code;
}

std::string synthesize_unitary(const Eigen::MatrixXcd& u) {
    unitaryGateNode node(u);
    qasmVisitor visitor(qubits_of(u.rows()));
    node.accept(visitor);
    return visitor.qasm_code;
}

std::string synthesize_diagonal(const Eigen::VectorXcd& d) {
    diagonalGateNode node(d.asDiagonal().toDenseMatrix());
    qasmVisitor visitor(qubits_of(d.size()));
    node.accept(visitor);
    return visitor.qasm_code;
}

std::mt19937& rng() {
    static std::mt19937 engine(20261028);
    return engine;
}

Eigen::VectorXcd random_state(Eigen::Index dim, bool real = false, double sparsity = 0.0) {
    std::normal_distribution<double> gauss;
    std::uniform_real_distribution<double> uniform;
    Eigen::VectorXcd v(dim);
    for (Eigen::Index i = 0; i < dim; ++i) {
        v(i) = (uniform(rng()) < sparsity) ? cd(0.0) : cd(gauss(rng()), real ? 0.0 : gauss(rng()));
    }
    if (v.norm() == 0.0) v(0) = 1.0;
    return v.normalized();
}

Eigen::MatrixXcd random_unitary(Eigen::Index dim) {
    std::normal_distribution<double> gauss;
    Eigen::MatrixXcd z(dim, dim);
    for (Eigen::Index i = 0; i < dim; ++i)
        for (Eigen::Index j = 0; j < dim; ++j) z(i, j) = cd(gauss(rng()), gauss(rng()));
    Eigen::HouseholderQR<Eigen::MatrixXcd> qr(z);
    return qr.householderQ() * Eigen::MatrixXcd::Identity(dim, dim);
}

Eigen::VectorXcd random_phases(Eigen::Index dim) {
    std::uniform_real_distribution<double> angle(-pi, pi);
    Eigen::VectorXcd d(dim);
    for (Eigen::Index i = 0; i < dim; ++i) d(i) = std::exp(cd(0.0, angle(rng())));
    return d;
}

Eigen::MatrixXcd permutation(Eigen::Index dim) {
    std::vector<Eigen::Index> p(static_cast<std::size_t>(dim));
    for (Eigen::Index i = 0; i < dim; ++i) p[static_cast<std::size_t>(i)] = i;
    std::shuffle(p.begin(), p.end(), rng());
    Eigen::MatrixXcd m = Eigen::MatrixXcd::Zero(dim, dim);
    for (Eigen::Index i = 0; i < dim; ++i) m(p[static_cast<std::size_t>(i)], i) = 1.0;
    return m;
}

Eigen::MatrixXcd qft(Eigen::Index dim) {
    Eigen::MatrixXcd m(dim, dim);
    for (Eigen::Index j = 0; j < dim; ++j)
        for (Eigen::Index k = 0; k < dim; ++k)
            m(j, k) = std::exp(cd(0.0, 2 * pi * double(j * k) / double(dim))) / std::sqrt(double(dim));
    return m;
}

Eigen::MatrixXcd kron(const Eigen::MatrixXcd& a, const Eigen::MatrixXcd& b) {
    Eigen::MatrixXcd m(a.rows() * b.rows(), a.cols() * b.cols());
    for (Eigen::Index i = 0; i < a.rows(); ++i)
        for (Eigen::Index j = 0; j < a.cols(); ++j) m.block(i * b.rows(), j * b.cols(), b.rows(), b.cols()) = a(i, j) * b;
    return m;
}

constexpr double kTolerance = 1e-9;

void expect_unitary(const Eigen::MatrixXcd& u, const std::string& name) {
    std::string qasm;
    ASSERT_NO_THROW(qasm = synthesize_unitary(u)) << name;
    const double error = (operator_of(qasm, u.rows()) - u).cwiseAbs().maxCoeff();
    EXPECT_LT(error, kTolerance) << name << " (" << qubits_of(u.rows()) << " qubits)";
}

}  // namespace

// ---------------------------------------------------------------------------
// State preparation
// ---------------------------------------------------------------------------

TEST(CorrectnessStatePreparation, RandomComplexStates) {
    for (Eigen::Index dim = 2; dim <= 64; dim *= 2) {
        for (int trial = 0; trial < 3; ++trial) {
            const Eigen::VectorXcd v = random_state(dim);
            const Eigen::VectorXcd out = simulate(synthesize_state(v), Eigen::VectorXcd::Unit(dim, 0));
            EXPECT_LT((out - v).cwiseAbs().maxCoeff(), kTolerance) << qubits_of(dim) << " qubits";
        }
    }
}

TEST(CorrectnessStatePreparation, RealAndSparseStates) {
    for (Eigen::Index dim = 2; dim <= 64; dim *= 2) {
        for (const auto& v : {random_state(dim, true), random_state(dim, false, 0.7)}) {
            const Eigen::VectorXcd out = simulate(synthesize_state(v), Eigen::VectorXcd::Unit(dim, 0));
            EXPECT_LT((out - v).cwiseAbs().maxCoeff(), kTolerance) << qubits_of(dim) << " qubits";
        }
    }
}

TEST(CorrectnessStatePreparation, BasisStatesWithPhase) {
    for (Eigen::Index k = 0; k < 8; ++k) {
        Eigen::VectorXcd v = Eigen::VectorXcd::Zero(8);
        v(k) = std::exp(cd(0.0, 0.7 * double(k + 1)));
        const Eigen::VectorXcd out = simulate(synthesize_state(v), Eigen::VectorXcd::Unit(8, 0));
        EXPECT_LT((out - v).cwiseAbs().maxCoeff(), kTolerance) << "|" << k << ">";
    }
}

// ---------------------------------------------------------------------------
// Unitary synthesis (exact, including the global phase)
// ---------------------------------------------------------------------------

TEST(CorrectnessUnitary, RandomUnitaries) {
    for (Eigen::Index dim = 2; dim <= 16; dim *= 2) {
        expect_unitary(random_unitary(dim), "random");
    }
}

TEST(CorrectnessUnitary, IdentityAndMinusIdentity) {
    for (Eigen::Index dim = 2; dim <= 16; dim *= 2) {
        expect_unitary(Eigen::MatrixXcd::Identity(dim, dim), "identity");
        expect_unitary(-Eigen::MatrixXcd::Identity(dim, dim), "-identity");
    }
}

TEST(CorrectnessUnitary, PermutationsAndDiagonals) {
    for (Eigen::Index dim = 2; dim <= 16; dim *= 2) {
        expect_unitary(permutation(dim), "permutation");
        expect_unitary(random_phases(dim).asDiagonal().toDenseMatrix(), "diagonal");
    }
}

TEST(CorrectnessUnitary, QftAndTensorProducts) {
    for (Eigen::Index dim = 4; dim <= 16; dim *= 2) {
        expect_unitary(qft(dim), "QFT");
        expect_unitary(kron(random_unitary(2), random_unitary(dim / 2)), "tensor product");
    }
}

TEST(CorrectnessUnitary, Toffoli) {
    Eigen::MatrixXcd toffoli = Eigen::MatrixXcd::Identity(8, 8);
    toffoli.row(6).swap(toffoli.row(7));
    expect_unitary(toffoli, "Toffoli");
}

// ---------------------------------------------------------------------------
// Diagonal gates
// ---------------------------------------------------------------------------

TEST(CorrectnessDiagonal, RandomPhases) {
    for (Eigen::Index dim = 2; dim <= 16; dim *= 2) {
        const Eigen::VectorXcd d = random_phases(dim);
        const double error = (operator_of(synthesize_diagonal(d), dim) - Eigen::MatrixXcd(d.asDiagonal())).cwiseAbs().maxCoeff();
        EXPECT_LT(error, kTolerance) << qubits_of(dim) << " qubits";
    }
}

// ---------------------------------------------------------------------------
// Output precision
// ---------------------------------------------------------------------------

TEST(CorrectnessPrecision, AnglesKeepFullPrecision) {
    // ry angle 2 * 0.1234567890123 needs more than the 6 decimals of std::to_string.
    const double half = 0.1234567890123;
    Eigen::VectorXcd v(2);
    v << std::cos(half), std::sin(half);
    const std::string qasm = synthesize_state(v);
    const std::size_t pos = qasm.find("ry(");
    ASSERT_NE(pos, std::string::npos);
    EXPECT_NEAR(parse_angle(qasm.substr(pos)), 2 * half, 1e-14);
}
