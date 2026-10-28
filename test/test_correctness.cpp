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
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <Eigen/Dense>

#include "../include/diagonalGateNode.h"
#include "../include/nodeVisitor.h"
#include "../include/qspUcrNode.h"
#include "../include/unitaryGateNode.h"
#include "../include/twoQubitsUnitary.h"
#include "../include/synthesisOptions.h"
#include "../include/hNode.h"
#include "../include/ucrzNode.h"
#include "../include/zxz.h"

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
        } else if (line.rfind("h ", 0) == 0) {
            Eigen::Matrix2cd h;
            h << 1.0, 1.0, 1.0, -1.0;
            apply_1q(state, h / std::sqrt(2.0), parse_index(line, 0));
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
// Gate-list output: same circuit as the QASM text, gate by gate
// ---------------------------------------------------------------------------

namespace {

// Parses the QASM text into the representation recorded by qasmVisitor(n, false).
std::pair<std::vector<circuitGate>, double> qasm_to_gates(const std::string& qasm) {
    std::vector<circuitGate> gates;
    double phase = 0.0;
    std::istringstream stream(qasm);
    std::string line;
    while (std::getline(stream, line)) {
        if (line.rfind("gphase", 0) == 0) {
            phase += parse_angle(line);
        } else if (line.rfind("rz", 0) == 0 || line.rfind("ry", 0) == 0 || line.rfind("p(", 0) == 0) {
            const auto type = line[0] == 'p' ? qasmVisitor::P : (line[1] == 'z' ? qasmVisitor::RZ : qasmVisitor::RY);
            gates.push_back({type, parse_index(line, line.find(')')), -1, parse_angle(line)});
        } else if (line.rfind("cx", 0) == 0 || line.rfind("cz", 0) == 0) {
            const auto type = line[1] == 'z' ? qasmVisitor::CZ : qasmVisitor::CX;
            gates.push_back({type, parse_index(line, 0), parse_index(line, line.find(',')), 0.0});
        } else if (line.rfind("h ", 0) == 0) {
            gates.push_back({qasmVisitor::H, parse_index(line, 0), -1, 0.0});
        }
    }
    return {gates, phase};
}

template <typename MakeNode>
void expect_same_circuit(MakeNode make_node, int n, const std::string& name) {
    auto text_node = make_node();
    qasmVisitor text(n);
    text_node->accept(text);

    auto list_node = make_node();  // the visitor modifies the tree: use a fresh one
    qasmVisitor list(n, false);
    list_node->accept(list);

    const auto [expected, phase] = qasm_to_gates(text.qasm_code);
    ASSERT_EQ(list.gates.size(), expected.size()) << name;
    for (std::size_t k = 0; k < expected.size(); ++k) {
        EXPECT_EQ(list.gates[k].type, expected[k].type) << name << " gate " << k;
        EXPECT_EQ(list.gates[k].q0, expected[k].q0) << name << " gate " << k;
        EXPECT_EQ(list.gates[k].q1, expected[k].q1) << name << " gate " << k;
        // exact: the QASM text prints angles with max_digits10 digits
        EXPECT_EQ(list.gates[k].angle, expected[k].angle) << name << " gate " << k;
    }
    EXPECT_NEAR(list.total_phase, phase, 1e-15) << name;
}

}  // namespace

TEST(CorrectnessGateList, SameCircuitAsQasm) {
    for (Eigen::Index dim = 2; dim <= 16; dim *= 2) {
        const int n = qubits_of(dim);
        const Eigen::MatrixXcd u = random_unitary(dim);
        expect_same_circuit([&] { return std::make_unique<unitaryGateNode>(u); }, n, "unitary");

        const Eigen::MatrixXcd d = random_phases(dim).asDiagonal().toDenseMatrix();
        expect_same_circuit([&] { return std::make_unique<diagonalGateNode>(d); }, n, "diagonal");

        const Eigen::VectorXcd v = random_state(dim);
        std::vector<cd> state(v.data(), v.data() + v.size());
        expect_same_circuit([&] { return std::make_unique<qspUcrNode>(&state); }, n, "state");
    }
}

// ---------------------------------------------------------------------------
// Two-qubit gate count (optimizations A.1 and A.2)
// ---------------------------------------------------------------------------

TEST(CorrectnessUnitary, TwoQubitGateCountWithA1A2) {
    // (23/48) 4^n - (3/2) 2^n + 4/3 two-qubit gates (Shende, Bullock, Markov 2006):
    // A.1 saves one CNOT per CSD step, A.2 one per two-qubit block except the last.
    // Structured matrices must give the same count: otherwise some two-qubit block
    // silently fell back to three CNOTs.
    const auto count_two_qubit = [](const Eigen::MatrixXcd& u) {
        unitaryGateNode node(u);
        qasmVisitor visitor(qubits_of(u.rows()), false);
        node.accept(visitor);
        std::size_t two_qubit = 0;
        for (const auto& gate : visitor.gates) two_qubit += (gate.type >= qasmVisitor::CX);
        return two_qubit;
    };
    const std::pair<int, std::size_t> expected[] = {{3, 20}, {4, 100}, {5, 444}};
    for (const auto& [n, count] : expected) {
        const Eigen::Index dim = Eigen::Index(1) << n;
        EXPECT_EQ(count_two_qubit(random_unitary(dim)), count) << n << " qubits, random";
        if (n > 4) continue;
        EXPECT_EQ(count_two_qubit(Eigen::MatrixXcd::Identity(dim, dim)), count) << n << " qubits, identity";
        EXPECT_EQ(count_two_qubit(permutation(dim)), count) << n << " qubits, permutation";
        EXPECT_EQ(count_two_qubit(qft(dim)), count) << n << " qubits, QFT";
        EXPECT_EQ(count_two_qubit(random_phases(dim).asDiagonal().toDenseMatrix()), count) << n << " qubits, diagonal";
    }
}

TEST(CorrectnessUnitary, WithoutA2EveryBlockUsesThreeCnots) {
    // optimize = false: no A.2 pass, every two-qubit block is synthesized with three
    // CNOTs when the tree is built (only A.1 applies): 23 and 115 two-qubit gates.
    const std::pair<int, std::size_t> expected[] = {{3, 23}, {4, 115}};
    for (const auto& [n, count] : expected) {
        const Eigen::MatrixXcd u = random_unitary(Eigen::Index(1) << n);
        unitaryGateNode node(u, false);
        qasmVisitor visitor(n, false);
        node.accept(visitor);
        std::size_t two_qubit = 0;
        for (const auto& gate : visitor.gates) two_qubit += (gate.type >= qasmVisitor::CX);
        EXPECT_EQ(two_qubit, count) << n << " qubits";
        qasmVisitor text(n);
        node.accept(text);
        EXPECT_LT((operator_of(text.qasm_code, u.rows()) - u).cwiseAbs().maxCoeff(), kTolerance) << n << " qubits";
    }
}

// ---------------------------------------------------------------------------
// Hadamard gate (used by the Block-ZXZ decomposition)
// ---------------------------------------------------------------------------

TEST(CorrectnessGates, Hadamard) {
    for (int qubit = 0; qubit < 3; ++qubit) {
        hNode node(qubit);
        qasmVisitor text(3);
        node.accept(text);
        EXPECT_NE(text.qasm_code.find("h q[" + std::to_string(qubit) + "];"), std::string::npos);

        qasmVisitor list(3, false);
        node.accept(list);
        ASSERT_EQ(list.gates.size(), 1u);
        EXPECT_EQ(list.gates[0].type, qasmVisitor::H);
        EXPECT_EQ(list.gates[0].q0, qubit);
        EXPECT_EQ(list.gates[0].q1, -1);

        // the simulated QASM is H on that qubit (little-endian), I elsewhere
        Eigen::Matrix2cd h;
        h << 1.0, 1.0, 1.0, -1.0;
        h /= std::sqrt(2.0);
        const Eigen::Matrix2cd i2 = Eigen::Matrix2cd::Identity();
        const Eigen::MatrixXcd expected = kron(qubit == 2 ? h : i2, kron(qubit == 1 ? h : i2, qubit == 0 ? h : i2));
        EXPECT_LT((operator_of(text.qasm_code, 8) - expected).cwiseAbs().maxCoeff(), 1e-15);
    }
}

// ---------------------------------------------------------------------------
// Synthesis options
// ---------------------------------------------------------------------------

TEST(CorrectnessUnitary, SynthesisOptions) {
    const Eigen::MatrixXcd u = random_unitary(16);
    const auto qasm_of = [](unitaryGateNode& node) {
        qasmVisitor visitor(4);
        node.accept(visitor);
        return visitor.qasm_code;
    };
    unitaryGateNode plain(u);
    unitaryGateNode defaults(u, synthesisOptions{});
    unitaryGateNode explicit_qsd(u, synthesisOptions{.method = unitaryDecomposition::QSD, .a2 = true});
    EXPECT_EQ(qasm_of(plain), qasm_of(defaults));
    EXPECT_EQ(qasm_of(plain), qasm_of(explicit_qsd));

    unitaryGateNode no_a2(u, synthesisOptions{.a2 = false});
    unitaryGateNode no_a2_bool(u, false);
    EXPECT_EQ(qasm_of(no_a2), qasm_of(no_a2_bool));
    EXPECT_NE(qasm_of(no_a2), qasm_of(plain));
    EXPECT_LT((operator_of(qasm_of(no_a2), 16) - u).cwiseAbs().maxCoeff(), kTolerance);
}

// ---------------------------------------------------------------------------
// A.2 building block: U = D * V, D diagonal, V with two CNOTs
// ---------------------------------------------------------------------------

namespace {

void expect_two_cx_up_to_diagonal(const Eigen::Matrix4cd& u, const std::string& name) {
    const Eigen::Matrix4cd d = diagonalForTwoCx(u);
    EXPECT_LT((d - Eigen::Matrix4cd(d.diagonal().asDiagonal())).cwiseAbs().maxCoeff(), 1e-15) << name;
    EXPECT_LT((d * d.adjoint() - Eigen::Matrix4cd::Identity()).cwiseAbs().maxCoeff(), 1e-14) << name;
    // twoCxDecomposition reconstructs its circuit and reports the max entry error
    EXPECT_LT(twoCxDecomposition(d.adjoint() * u).error, 1e-12) << name;
}

Eigen::Matrix4cd controlled(const Eigen::Matrix2cd& g) {  // control q0, target q1
    Eigen::Matrix4cd m = Eigen::Matrix4cd::Identity();
    m(1, 1) = g(0, 0); m(1, 3) = g(0, 1); m(3, 1) = g(1, 0); m(3, 3) = g(1, 1);
    return m;
}

}  // namespace

TEST(CorrectnessTwoQubit, TwoCxUpToDiagonal) {
    Eigen::Matrix4cd swap = Eigen::Matrix4cd::Zero();
    swap(0, 0) = swap(1, 2) = swap(2, 1) = swap(3, 3) = 1.0;
    Eigen::Matrix4cd iswap = swap;
    iswap(1, 2) = iswap(2, 1) = cd(0.0, 1.0);
    Eigen::Matrix2cd x, h, rz;
    x << 0, 1, 1, 0;
    h << 1, 1, 1, -1;
    h /= std::sqrt(2.0);
    rz << std::exp(cd(0, -0.2)), 0, 0, std::exp(cd(0, 0.2));

    expect_two_cx_up_to_diagonal(Eigen::Matrix4cd::Identity(), "identity");
    expect_two_cx_up_to_diagonal(swap, "SWAP");
    expect_two_cx_up_to_diagonal(iswap, "iSWAP");
    expect_two_cx_up_to_diagonal(controlled(x), "CNOT");
    expect_two_cx_up_to_diagonal(controlled(h), "CH");
    expect_two_cx_up_to_diagonal(controlled(rz), "CRZ");
    expect_two_cx_up_to_diagonal(Eigen::Matrix4cd(Eigen::Vector4cd(1, 1, 1, -1).asDiagonal()), "CZ");
    expect_two_cx_up_to_diagonal(kron(random_unitary(2), random_unitary(2)), "local");
    for (int k = 0; k < 300; ++k) {
        expect_two_cx_up_to_diagonal(random_unitary(4), "random " + std::to_string(k));
    }
}

// ---------------------------------------------------------------------------
// The visitor must not modify the tree: visiting it again gives the same circuit
// ---------------------------------------------------------------------------

namespace {

std::string visit_qasm(IASTnode& node, int n) {
    qasmVisitor visitor(n);
    node.accept(visitor);
    return visitor.qasm_code;
}

}  // namespace

TEST(CorrectnessVisitor, VisitingTwiceGivesTheSameCorrectCircuit) {
    for (Eigen::Index dim = 4; dim <= 32; dim *= 2) {
        const int n = qubits_of(dim);

        const Eigen::MatrixXcd u = random_unitary(dim);
        unitaryGateNode unitary_node(u);
        const std::string first = visit_qasm(unitary_node, n);
        const std::string second = visit_qasm(unitary_node, n);
        EXPECT_EQ(first, second) << "unitary, " << n << " qubits";
        EXPECT_LT((operator_of(second, dim) - u).cwiseAbs().maxCoeff(), kTolerance) << "unitary, " << n << " qubits";

        const Eigen::VectorXcd d = random_phases(dim);
        diagonalGateNode diagonal_node(d.asDiagonal().toDenseMatrix());
        EXPECT_EQ(visit_qasm(diagonal_node, n), visit_qasm(diagonal_node, n)) << "diagonal, " << n << " qubits";

        const Eigen::VectorXcd v = random_state(dim);
        std::vector<cd> state(v.data(), v.data() + v.size());
        qspUcrNode state_node(&state);
        const std::string first_state = visit_qasm(state_node, n);
        const std::string second_state = visit_qasm(state_node, n);
        EXPECT_EQ(first_state, second_state) << "state, " << n << " qubits";
        EXPECT_LT((simulate(second_state, Eigen::VectorXcd::Unit(dim, 0)) - v).cwiseAbs().maxCoeff(), kTolerance)
            << "state, " << n << " qubits";
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

// ---------------------------------------------------------------------------
// RZ multiplexor variants used by the Block-ZXZ: mirrored order, outer CNOT omitted
// ---------------------------------------------------------------------------

namespace {

std::vector<circuitGate> ucrz_gates(const std::vector<double>& angles, bool mirror, bool drop_first, bool drop_last) {
    firstUcrzNode node(&angles);
    node.mirror = mirror;
    node.drop_first = drop_first;
    node.drop_last = drop_last;
    qasmVisitor visitor(qubits_of(Eigen::Index(2 * angles.size())), false);
    visitor.gate_synthesis = true;  // as inside a unitary: controls are the lower qubits
    node.accept(visitor);
    return visitor.gates;
}

std::string ucrz_qasm(const std::vector<double>& angles, bool mirror) {
    firstUcrzNode node(&angles);
    node.mirror = mirror;
    const int n = qubits_of(Eigen::Index(2 * angles.size()));
    qasmVisitor visitor(n);
    visitor.gate_synthesis = true;
    node.accept(visitor);
    return visitor.qasm_code;
}

}  // namespace

TEST(CorrectnessUcrz, MirroredIsTheReversedCircuitWithTheSameOperator) {
    std::uniform_real_distribution<double> angle(-3.0, 3.0);
    for (std::size_t size = 2; size <= 16; size *= 2) {
        std::vector<double> angles(size);
        for (auto& a : angles) a = angle(rng());
        const int n = qubits_of(Eigen::Index(2 * size));
        const Eigen::Index dim = Eigen::Index(1) << n;

        // RZ multiplexor on the MSB: diag(e^{-i a_j / 2}) for MSB = 0, conjugate for MSB = 1
        Eigen::VectorXcd expected(dim);
        for (Eigen::Index j = 0; j < dim / 2; ++j) {
            expected(j) = std::polar(1.0, -angles[j] / 2);
            expected(j + dim / 2) = std::polar(1.0, angles[j] / 2);
        }
        const Eigen::MatrixXcd target = expected.asDiagonal().toDenseMatrix();
        EXPECT_LT((operator_of(ucrz_qasm(angles, false), dim) - target).cwiseAbs().maxCoeff(), 1e-14) << n;
        EXPECT_LT((operator_of(ucrz_qasm(angles, true), dim) - target).cwiseAbs().maxCoeff(), 1e-14) << n;

        auto normal = ucrz_gates(angles, false, false, false);
        const auto mirrored = ucrz_gates(angles, true, false, false);
        std::reverse(normal.begin(), normal.end());
        ASSERT_EQ(normal.size(), mirrored.size()) << n;
        for (std::size_t k = 0; k < normal.size(); ++k) {
            EXPECT_EQ(normal[k].type, mirrored[k].type) << n << " gate " << k;
            EXPECT_EQ(normal[k].q0, mirrored[k].q0) << n << " gate " << k;
            EXPECT_EQ(normal[k].q1, mirrored[k].q1) << n << " gate " << k;
            EXPECT_EQ(normal[k].angle, mirrored[k].angle) << n << " gate " << k;
        }
    }
}

TEST(CorrectnessUcrz, DroppedOuterCnot) {
    const std::vector<double> angles = {0.1, -0.7, 1.3, 2.2};  // three qubits, four CNOTs
    const auto full = ucrz_gates(angles, false, false, false);
    ASSERT_EQ(full.back().type, qasmVisitor::CX);
    EXPECT_EQ(full.back().q0, 1);
    EXPECT_EQ(full.back().q1, 2);

    // drop_last: the same circuit without the last CNOT (1 -> 2)
    const auto no_last = ucrz_gates(angles, false, false, true);
    ASSERT_EQ(no_last.size() + 1, full.size());
    EXPECT_TRUE(std::equal(no_last.begin(), no_last.end(), full.begin(), [](const auto& a, const auto& b) {
        return a.type == b.type && a.q0 == b.q0 && a.q1 == b.q1 && a.angle == b.angle;
    }));

    // mirrored with drop_first: the mirrored circuit without its first CNOT (1 -> 2)
    const auto mirrored = ucrz_gates(angles, true, false, false);
    const auto no_first = ucrz_gates(angles, true, true, false);
    ASSERT_EQ(mirrored.front().type, qasmVisitor::CX);
    EXPECT_EQ(mirrored.front().q0, 1);
    ASSERT_EQ(no_first.size() + 1, mirrored.size());
    EXPECT_TRUE(std::equal(no_first.begin(), no_first.end(), mirrored.begin() + 1, [](const auto& a, const auto& b) {
        return a.type == b.type && a.q0 == b.q0 && a.q1 == b.q1 && a.angle == b.angle;
    }));

    // drop_first only applies to a mirrored node, drop_last only to a normal one
    EXPECT_EQ(ucrz_gates(angles, false, true, false).size(), full.size());
    EXPECT_EQ(ucrz_gates(angles, true, false, true).size(), full.size());
}

// ---------------------------------------------------------------------------
// Block-ZXZ decomposition (Krol and Al-Ars, 2024)
// ---------------------------------------------------------------------------

namespace {

const synthesisOptions kZxz{.method = unitaryDecomposition::BlockZXZ};

std::string synthesize_zxz(const Eigen::MatrixXcd& u, bool a2 = true) {
    unitaryGateNode node(u, synthesisOptions{.method = unitaryDecomposition::BlockZXZ, .a2 = a2});
    qasmVisitor visitor(qubits_of(u.rows()));
    node.accept(visitor);
    return visitor.qasm_code;
}

void expect_zxz_unitary(const Eigen::MatrixXcd& u, const std::string& name) {
    std::string qasm;
    ASSERT_NO_THROW(qasm = synthesize_zxz(u)) << name;
    const double error = (operator_of(qasm, u.rows()) - u).cwiseAbs().maxCoeff();
    EXPECT_LT(error, kTolerance) << name << " (" << qubits_of(u.rows()) << " qubits)";
}

std::size_t zxz_two_qubit_count(const Eigen::MatrixXcd& u, bool a2 = true) {
    unitaryGateNode node(u, synthesisOptions{.method = unitaryDecomposition::BlockZXZ, .a2 = a2});
    qasmVisitor visitor(qubits_of(u.rows()), false);
    node.accept(visitor);
    std::size_t two_qubit = 0;
    for (const auto& gate : visitor.gates) two_qubit += (gate.type == qasmVisitor::CX || gate.type == qasmVisitor::CZ);
    return two_qubit;
}

std::vector<std::pair<std::string, Eigen::MatrixXcd>> structured_unitaries(Eigen::Index dim) {
    Eigen::MatrixXcd msb_x = Eigen::MatrixXcd::Zero(dim, dim);  // X on the MSB: U11 = 0
    msb_x.topRightCorner(dim / 2, dim / 2).setIdentity();
    msb_x.bottomLeftCorner(dim / 2, dim / 2).setIdentity();
    return {
        {"identity", Eigen::MatrixXcd::Identity(dim, dim)},
        {"-identity", -Eigen::MatrixXcd::Identity(dim, dim)},
        {"permutation", permutation(dim)},
        {"diagonal", random_phases(dim).asDiagonal().toDenseMatrix()},
        {"QFT", qft(dim)},
        {"tensor product", kron(random_unitary(2), random_unitary(dim / 2))},
        {"X on the MSB", msb_x},
    };
}

}  // namespace

TEST(CorrectnessZxz, BlockDecomposition) {
    for (Eigen::Index dim = 2; dim <= 64; dim *= 2) {
        auto cases = structured_unitaries(dim);
        cases.emplace_back("random", random_unitary(dim));
        for (const auto& [name, u] : cases) {
            const zxz_result r = blockZxz(u);
            const Eigen::Index h = dim / 2;
            const Eigen::MatrixXcd id = Eigen::MatrixXcd::Identity(h, h);
            for (const auto* block : {&r.a1, &r.a2, &r.b, &r.c}) {
                EXPECT_LT((*block * block->adjoint() - id).cwiseAbs().maxCoeff(), 1e-12) << name << ", dim " << dim;
            }
            // U = (A1 ⊕ A2) ½[[I + B, I - B], [I - B, I + B]] (I ⊕ C)
            Eigen::MatrixXcd rebuilt(dim, dim);
            rebuilt << r.a1 * (id + r.b) / 2, r.a1 * (id - r.b) * r.c / 2,
                       r.a2 * (id - r.b) / 2, r.a2 * (id + r.b) * r.c / 2;
            EXPECT_LT((rebuilt - u).cwiseAbs().maxCoeff(), 1e-12) << name << ", dim " << dim;
        }
    }
}

TEST(CorrectnessZxz, RandomUnitaries) {
    for (Eigen::Index dim = 2; dim <= 64; dim *= 2) {
        expect_zxz_unitary(random_unitary(dim), "random");
    }
}

TEST(CorrectnessZxz, StructuredUnitaries) {
    for (Eigen::Index dim = 8; dim <= 32; dim *= 2) {
        for (const auto& [name, u] : structured_unitaries(dim)) {
            expect_zxz_unitary(u, name);
        }
    }
    Eigen::MatrixXcd toffoli = Eigen::MatrixXcd::Identity(8, 8);
    toffoli.row(6).swap(toffoli.row(7));
    expect_zxz_unitary(toffoli, "Toffoli");
}

TEST(CorrectnessZxz, TwoQubitGateCountWithA2) {
    // One CNOT less per decomposition step than the QSD (20, 100, 444, 1868), the
    // same counts as Krol and Al-Ars. Structured matrices must give the same count:
    // otherwise some two-qubit block silently fell back to three CNOTs.
    const std::pair<int, std::size_t> expected[] = {{3, 19}, {4, 95}, {5, 423}, {6, 1783}};
    for (const auto& [n, count] : expected) {
        const Eigen::Index dim = Eigen::Index(1) << n;
        EXPECT_EQ(zxz_two_qubit_count(random_unitary(dim)), count) << n << " qubits, random";
        if (n > 4) continue;
        for (const auto& [name, u] : structured_unitaries(dim)) {
            EXPECT_EQ(zxz_two_qubit_count(u), count) << n << " qubits, " << name;
        }
    }
}

TEST(CorrectnessZxz, WithoutA2EveryBlockUsesThreeCnots) {
    // 3·2^(n-1) - 2 CNOTs per step and three per two-qubit block: 22 and 110
    const std::pair<int, std::size_t> expected[] = {{3, 22}, {4, 110}};
    for (const auto& [n, count] : expected) {
        const Eigen::MatrixXcd u = random_unitary(Eigen::Index(1) << n);
        EXPECT_EQ(zxz_two_qubit_count(u, false), count) << n << " qubits";
        EXPECT_LT((operator_of(synthesize_zxz(u, false), u.rows()) - u).cwiseAbs().maxCoeff(), kTolerance) << n;
    }
}

TEST(CorrectnessZxz, OnlyRzRyCxAndH) {
    // no RY multiplexor and no CZ: the two-qubit gates are CNOTs
    unitaryGateNode node(random_unitary(16), kZxz);
    qasmVisitor visitor(4, false);
    node.accept(visitor);
    std::size_t h = 0;
    for (const auto& gate : visitor.gates) {
        EXPECT_NE(gate.type, qasmVisitor::CZ);
        EXPECT_NE(gate.type, qasmVisitor::P);
        h += gate.type == qasmVisitor::H;
    }
    EXPECT_EQ(h, 2u + 4u * 2u);  // two per decomposition step (4 qubits, then four of 3 qubits)
}

TEST(CorrectnessZxz, GateListAndRepeatedVisits) {
    for (Eigen::Index dim = 8; dim <= 32; dim *= 2) {
        const int n = qubits_of(dim);
        const Eigen::MatrixXcd u = random_unitary(dim);
        expect_same_circuit([&] { return std::make_unique<unitaryGateNode>(u, kZxz); }, n, "zxz");

        unitaryGateNode node(u, kZxz);
        const std::string first = visit_qasm(node, n);
        EXPECT_EQ(first, visit_qasm(node, n)) << n << " qubits";
        EXPECT_LT((operator_of(first, dim) - u).cwiseAbs().maxCoeff(), kTolerance) << n << " qubits";
    }
}

TEST(CorrectnessZxz, QsdIsStillTheDefault) {
    const Eigen::MatrixXcd u = random_unitary(16);
    unitaryGateNode plain(u);
    unitaryGateNode zxz(u, kZxz);
    EXPECT_EQ(visit_qasm(plain, 4), synthesize_unitary(u));
    EXPECT_NE(visit_qasm(plain, 4), visit_qasm(zxz, 4));
    EXPECT_EQ(visit_qasm(plain, 4).find("h q"), std::string::npos);
}
