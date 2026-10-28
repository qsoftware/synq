"""
Correctness tests for the SynQ Python module.

The generated OpenQASM is simulated with a small numpy simulator (no Qiskit
needed) and compared with the target, including the global phase. Qubit order
is little-endian (qubit 0 is the least significant bit), as in Qiskit.

Run with:  python -m pytest python/tests
"""

import re
import numpy as np
import pytest
import synq

_ROT = re.compile(r"^(rz|ry|p)\(([^)]+)\) q\[(\d+)\];$")
_CX = re.compile(r"^(cx|cz) q\[(\d+)\], q\[(\d+)\];$")
_GPHASE = re.compile(r"^gphase\(([^)]+)\);$")


def _one_qubit(name, theta):
    if name == "rz":
        return np.diag([np.exp(-0.5j * theta), np.exp(0.5j * theta)])
    if name == "ry":
        c, s = np.cos(theta / 2), np.sin(theta / 2)
        return np.array([[c, -s], [s, c]], dtype=complex)
    return np.diag([1, np.exp(1j * theta)])


def _apply_1q(state, gate, qubit, n):
    psi = state.reshape([2] * n)  # axis k <-> qubit n-1-k
    axis = n - 1 - qubit
    psi = np.moveaxis(np.tensordot(gate, psi, axes=([1], [axis])), 0, axis)
    return psi.reshape(-1)


def _apply_cx(state, control, target):
    idx = np.arange(len(state))
    flipped = np.where((idx >> control) & 1, idx ^ (1 << target), idx)
    return state[flipped]


def _apply_cz(state, control, target):
    idx = np.arange(len(state))
    sign = np.where(((idx >> control) & 1) & ((idx >> target) & 1), -1.0, 1.0)
    return state * sign


def simulate(qasm, state):
    n = int(re.search(r"qubit\[(\d+)\]", qasm).group(1))
    state = np.array(state, dtype=complex)
    for line in map(str.strip, qasm.splitlines()):
        if m := _ROT.match(line):
            state = _apply_1q(state, _one_qubit(m[1], float(m[2])), int(m[3]), n)
        elif m := _CX.match(line):
            apply = _apply_cx if m[1] == "cx" else _apply_cz
            state = apply(state, int(m[2]), int(m[3]))
        elif m := _GPHASE.match(line):
            state = state * np.exp(1j * float(m[1]))
        elif line and not line.startswith(("OPENQASM", "include", "qubit")):
            raise ValueError(f"unexpected line {line!r}")
    return state


def operator(qasm):
    n = int(re.search(r"qubit\[(\d+)\]", qasm).group(1))
    return np.column_stack([simulate(qasm, col) for col in np.eye(2**n)])


def random_unitary(dim, rng):
    z = rng.normal(size=(dim, dim)) + 1j * rng.normal(size=(dim, dim))
    q, r = np.linalg.qr(z)
    return q * (np.diag(r) / np.abs(np.diag(r)))


rng = np.random.default_rng(1234)


@pytest.mark.parametrize("n", range(1, 8))
@pytest.mark.parametrize("kind", ["complex", "real", "sparse", "basis"])
def test_state_preparation(n, kind):
    v = rng.normal(size=2**n) + 1j * rng.normal(size=2**n)
    if kind == "real":
        v = v.real
    elif kind == "sparse":
        v[rng.random(2**n) < 0.7] = 0
        v[0] += 1
    elif kind == "basis":
        v = np.zeros(2**n, dtype=complex)
        v[rng.integers(2**n)] = np.exp(1j * rng.uniform(0, 2 * np.pi))
    v = v / np.linalg.norm(v)
    zero = np.zeros(2**n, dtype=complex)
    zero[0] = 1
    assert np.allclose(simulate(synq.state_preparation(v), zero), v, atol=1e-12)


def _structured(n):
    dim = 2**n
    yield np.eye(dim)
    yield -np.eye(dim)
    yield np.eye(dim)[rng.permutation(dim)]
    yield np.diag(np.exp(1j * rng.uniform(0, 2 * np.pi, dim)))
    yield random_unitary(dim, rng)
    if n > 1:
        yield np.kron(random_unitary(2, rng), random_unitary(dim // 2, rng))
        omega = np.exp(2j * np.pi / dim)
        yield np.array([[omega ** (j * k) for k in range(dim)] for j in range(dim)]) / np.sqrt(dim)


@pytest.mark.parametrize("n", range(1, 5))
def test_unitary(n):
    for u in _structured(n):
        assert np.allclose(operator(synq.unitary(u)), u, atol=1e-10)


def qasm_to_gates(qasm):
    """Parses the QASM text into the same representation as the *_gates functions."""
    names = {"rz": synq.RZ, "ry": synq.RY, "p": synq.P}
    types, qubits, angles, phase = [], [], [], 0.0
    for line in map(str.strip, qasm.splitlines()):
        if m := _ROT.match(line):
            types.append(names[m[1]]); qubits.append((int(m[3]), -1)); angles.append(float(m[2]))
        elif m := _CX.match(line):
            types.append(synq.CX if m[1] == "cx" else synq.CZ)
            qubits.append((int(m[2]), int(m[3]))); angles.append(0.0)
        elif m := _GPHASE.match(line):
            phase += float(m[1])
    return types, qubits, angles, phase


def assert_same_circuit(qasm, gates):
    """The gate list and the QASM text describe exactly the same circuit."""
    types, qubits, angles, phase = gates
    q_types, q_qubits, q_angles, q_phase = qasm_to_gates(qasm)
    assert types.dtype == np.uint8 and qubits.dtype == np.int32 and angles.dtype == np.float64
    assert types.tolist() == q_types
    assert [tuple(q) for q in qubits.tolist()] == q_qubits
    assert angles.tolist() == q_angles  # exact: QASM angles are printed with 17 digits
    assert phase == pytest.approx(q_phase, abs=1e-15)


@pytest.mark.parametrize("n", range(1, 6))
def test_gate_lists_match_qasm(n):
    u = random_unitary(2**n, rng)
    assert_same_circuit(synq.unitary(u), synq.unitary_gates(u))
    d = np.exp(1j * rng.uniform(-np.pi, np.pi, 2**n))
    assert_same_circuit(synq.diagonal(d), synq.diagonal_gates(d))
    v = rng.normal(size=2**n) + 1j * rng.normal(size=2**n)
    v /= np.linalg.norm(v)
    assert_same_circuit(synq.state_preparation(v), synq.state_preparation_gates(v))


def simulate_gates(gates, state):
    types, qubits, angles, phase = gates
    n = int(np.log2(len(state)))
    state = np.array(state, dtype=complex) * np.exp(1j * phase)
    names = {synq.RZ: "rz", synq.RY: "ry", synq.P: "p"}
    for gate, (q0, q1), angle in zip(types, qubits, angles):
        if gate in names:
            state = _apply_1q(state, _one_qubit(names[gate], angle), q0, n)
        elif gate == synq.CX:
            state = _apply_cx(state, q0, q1)
        else:
            state = _apply_cz(state, q0, q1)
    return state


@pytest.mark.parametrize("n", range(1, 5))
def test_gate_lists_are_correct(n):
    u = random_unitary(2**n, rng)
    out = np.column_stack([simulate_gates(synq.unitary_gates(u), c) for c in np.eye(2**n)])
    assert np.allclose(out, u, atol=1e-10)


def test_gate_lists_validate_inputs():
    with pytest.raises(ValueError):
        synq.state_preparation_gates([1.0, 1.0])
    with pytest.raises(ValueError):
        synq.unitary_gates(np.ones((4, 4)))
    with pytest.raises(ValueError):
        synq.diagonal_gates([1.0, 0.5])


@pytest.mark.parametrize("n,expected", [(2, 3), (3, 20), (4, 100), (5, 444)])
def test_unitary_two_qubit_gate_count(n, expected):
    """Optimizations A.1 and A.2: (23/48) 4^n - (3/2) 2^n + 4/3 two-qubit gates for n >= 3."""
    types, _, _, _ = synq.unitary_gates(random_unitary(2**n, rng))
    assert np.count_nonzero(types >= synq.CX) == expected


def test_toffoli():
    u = np.eye(8)
    u[[6, 7]] = u[[7, 6]]
    assert np.allclose(operator(synq.unitary(u)), u, atol=1e-10)


@pytest.mark.parametrize("n", range(1, 6))
def test_diagonal(n):
    d = np.exp(1j * rng.uniform(-np.pi, np.pi, 2**n))
    assert np.allclose(operator(synq.diagonal(d)), np.diag(d), atol=1e-12)


def test_full_precision_angles():
    qasm = synq.state_preparation([np.cos(0.1), np.sin(0.1)])
    angle = float(re.search(r"ry\(([^)]+)\)", qasm).group(1))
    assert angle == pytest.approx(0.2, abs=1e-15)


def test_invalid_inputs():
    with pytest.raises(ValueError):
        synq.state_preparation([1.0, 1.0])
    with pytest.raises(ValueError):
        synq.state_preparation([1.0, 0.0, 0.0])
    with pytest.raises(ValueError):
        synq.unitary(np.ones((4, 4)))
    with pytest.raises(ValueError):
        synq.diagonal([1.0, 0.5])
