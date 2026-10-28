"""
SynQ: quantum circuit synthesis.

Python interface to the SynQ C++ library. ``state_preparation``, ``unitary``
and ``diagonal`` return an OpenQASM 3 program (``str``) using only the gates
``rz``, ``ry``, ``p``, ``cx``, ``cz``, ``h`` and ``gphase``. The ``*_gates`` variants return
the same circuit as numpy arrays, which avoids formatting and parsing text.

Functions
---------
state_preparation(state)
    Prepares the normalized vector ``state`` from ``|0...0>``.
unitary(matrix, method="qsd")
    Implements a ``2^n x 2^n`` unitary matrix with the cosine-sine / quantum Shannon
    decomposition (``method="qsd"``) or the Block-ZXZ decomposition (``method="zxz"``,
    one CNOT less per decomposition step: 19, 95, 423 instead of 20, 100, 444 for
    3, 4, 5 qubits).
diagonal(diagonal)
    Implements ``diag(diagonal)`` (entries with modulus 1).
state_preparation_gates, unitary_gates, diagonal_gates
    Same circuits as ``(types, qubits, angles, global_phase)`` numpy arrays.
"""

import os as _os
import sys as _sys

import numpy as _np


def _add_dll_directories():
    """
    On Windows, Python >= 3.8 does not use PATH to find the DLLs an extension
    depends on (MKL, OpenMP). Register the directories recorded at build time
    and those in the SYNQ_DLL_DIRS environment variable (separated by ';').
    """
    if _sys.platform != "win32" or not hasattr(_os, "add_dll_directory"):
        return
    try:
        from ._build_config import DLL_DIRS  # pylint: disable=import-outside-toplevel
    except ImportError:
        DLL_DIRS = []  # pylint: disable=invalid-name
    extra = [d for d in _os.environ.get("SYNQ_DLL_DIRS", "").split(";") if d]
    for directory in dict.fromkeys(list(DLL_DIRS) + extra):
        if _os.path.isdir(directory):
            _os.add_dll_directory(directory)


_add_dll_directories()

try:
    from . import _synq
except ImportError as _error:
    raise ImportError(
        "Could not load the SynQ extension module. On Windows, make sure the "
        "MKL/OpenMP DLLs can be found: set SYNQ_DLL_DIRS to the folder that "
        "contains them (e.g. <vcpkg>/installed/x64-windows/bin)."
    ) from _error

__version__ = getattr(_synq, "__version__", "0.0.1")


def state_preparation(state):
    """OpenQASM 3 circuit that prepares ``state`` (normalized, length ``2^n``)."""
    return _synq.state_preparation(_np.asarray(state, dtype=complex).ravel())


def unitary(matrix, method="qsd"):
    """
    OpenQASM 3 circuit that implements the unitary ``matrix`` (``2^n x 2^n``).

    ``method`` is ``"qsd"`` (quantum Shannon decomposition, default) or ``"zxz"``
    (Block-ZXZ decomposition).
    """
    return _synq.unitary(_np.asarray(matrix, dtype=complex), method)


def diagonal(diagonal):  # pylint: disable=redefined-outer-name
    """OpenQASM 3 circuit that implements ``diag(diagonal)``."""
    return _synq.diagonal(_np.asarray(diagonal, dtype=complex).ravel())


# Gate-list variants: return (types, qubits, angles, global_phase) as numpy arrays.
# types: 0 rz, 1 ry, 2 p, 3 cx, 4 cz, 5 h; qubits[k] = (target, -1) or (control, target).
RZ, RY, P, CX, CZ, H = 0, 1, 2, 3, 4, 5


def state_preparation_gates(state):
    """Gate list (types, qubits, angles, global_phase) that prepares ``state``."""
    return _synq.state_preparation_gates(_np.asarray(state, dtype=complex).ravel())


def unitary_gates(matrix, method="qsd"):
    """Gate list (types, qubits, angles, global_phase) that implements ``matrix`` (see ``unitary``)."""
    return _synq.unitary_gates(_np.asarray(matrix, dtype=complex), method)


def diagonal_gates(diagonal):  # pylint: disable=redefined-outer-name
    """Gate list (types, qubits, angles, global_phase) that implements ``diag(diagonal)``."""
    return _synq.diagonal_gates(_np.asarray(diagonal, dtype=complex).ravel())


__all__ = [
    "state_preparation", "unitary", "diagonal",
    "state_preparation_gates", "unitary_gates", "diagonal_gates",
    "RZ", "RY", "P", "CX", "CZ", "H", "__version__",
]
