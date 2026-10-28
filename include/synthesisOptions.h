/**
 * @file synthesisOptions.h
 * @brief Options of the unitary synthesis.
 */

#pragma once

/** Decomposition used for unitaries with three or more qubits. */
enum class unitaryDecomposition {
    /** Cosine-sine decomposition + quantum Shannon decomposition (Shende, Bullock, Markov 2006). */
    QSD,
};

/**
 * @brief Options of the unitary synthesis.
 *
 * Example: unitaryGateNode node(matrix, synthesisOptions{.a2 = false});
 */
struct synthesisOptions {
    /** Decomposition used for unitaries with three or more qubits. */
    unitaryDecomposition method = unitaryDecomposition::QSD;
    /**
     * Optimization A.2: two-qubit blocks are implemented up to a diagonal with two
     * CNOTs and the diagonal is moved into the next block.
     */
    bool a2 = true;
};
