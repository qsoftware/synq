# SynQ

Quantum circuit synthesis in C++ (Eigen + LAPACK), with Python bindings.

SynQ builds an abstract syntax tree for the decomposition and emits OpenQASM 3 that uses only `rz`, `ry`, `p`, `cx` and `gphase`:

| Function | Method |
|---|---|
| state preparation | uniformly controlled RY/RZ rotations (Möttönen et al.) |
| unitary | cosine-sine decomposition + quantum Shannon decomposition; 2-qubit blocks use the 3-CNOT decomposition |
| diagonal gate | recursive uniformly controlled RZ rotations |

The circuits use the Qiskit (little-endian) qubit order: qubit 0 is the least significant bit. They implement the target exactly, including the global phase.

## C++ build

Requirements: CMake ≥ 3.14, a C++20 compiler, Eigen 3, and LAPACKE (Linux/macOS) or Intel MKL (Windows, e.g. through vcpkg). The tests also need GoogleTest.

```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build
```

Options: `SYNQ_BUILD_TESTS` (default `ON`) and `SYNQ_BUILD_PYTHON` (default `OFF`; needs pybind11).

## Python package

```
pip install .
```

`pip` builds the extension with scikit-build-core and pybind11. The C++ tests are not built in this mode, so GoogleTest is not needed. On Windows, point CMake at vcpkg so that Eigen and MKL are found, for example:

```
set CMAKE_ARGS=-DCMAKE_TOOLCHAIN_FILE=C:/path/to/vcpkg/scripts/buildsystems/vcpkg.cmake
pip install .
```

The build records the vcpkg `bin` folder, and `import synq` registers it with `os.add_dll_directory` so the MKL/OpenMP DLLs are found. If they are somewhere else, set `SYNQ_DLL_DIRS` (folders separated by `;`).

```python
import numpy as np
import synq

qasm = synq.state_preparation(np.array([1, 1j, -1, 0]) / np.sqrt(3))
qasm = synq.unitary(matrix)         # 2^n x 2^n unitary
qasm = synq.diagonal(phases)        # entries with modulus 1

types, qubits, angles, phase = synq.unitary_gates(matrix)   # same circuit as arrays
```

`state_preparation`, `unitary` and `diagonal` return an OpenQASM 3 string. The `*_gates` variants (`state_preparation_gates`, `unitary_gates`, `diagonal_gates`) return the same circuit as `(types, qubits, angles, global_phase)` numpy arrays: types are `synq.RZ`, `RY`, `P`, `CX` (0–3), `qubits[k]` is `(target, -1)` for rotations and `(control, target)` for `cx`. They avoid formatting and parsing text, so they are the faster way to consume the circuit. [qclib](https://github.com/qclib/qclib) converts it into a Qiskit circuit (`qclib.synq`, `SynqInitialize`, `unitary(..., decomposition="synq")`).

Python tests (they simulate the output with numpy and check it against the target):

```
python -m pytest python/tests
```
