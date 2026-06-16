#pragma once

#include "DDpackage.h"
#include "InteractionGraph.h"
#include "operations/StandardOperation.hpp"
#include <array>
#include <cassert>
#include <cmath>
#include <vector>

namespace dm {

// ---------------------------------------------------------------------------
// Density matrix construction
// ---------------------------------------------------------------------------

/// Build rho = |psi><psi| from a state-vector DD (mode=Vector).
/// Returns a matrix DD with incRef; caller owns the ref.
inline dd::Edge densityMatrixFromState(std::unique_ptr<dd::Package>& dd,
                                       dd::Edge psi) {
    dd->incRef(psi);
    dd->setMode(dd::Matrix);
    dd::Edge psiDag = dd->conjugateTranspose(psi);
    dd->incRef(psiDag);

    dd::Edge rho = dd->multiply(psi, psiDag);
    dd->incRef(rho);

    dd->decRef(psiDag);
    dd->decRef(psi);
    dd->garbageCollect();

    return rho;
}

// ---------------------------------------------------------------------------
// Kraus operator
// ---------------------------------------------------------------------------

/// A single-qubit Kraus operator: flat [K[0][0], K[0][1], K[1][0], K[1][1]]
using KrausOp = qc::GateMatrix;  // std::array<dd::ComplexValue, 4>

/// Apply single-qubit Kraus channel {K_m} to qubit `target` of density matrix `rho`.
///   rho' = sum_m  K_m * rho * K_m†
/// Returns new DD with incRef; caller owns the ref.
inline dd::Edge applyKrausChannel(std::unique_ptr<dd::Package>& dd,
                                  dd::Edge rho,
                                  const std::vector<KrausOp>& kraus,
                                  unsigned short target,
                                  unsigned short nqubits) {
    dd->setMode(dd::Matrix);
    std::array<short, dd::MAXN> line;
    line.fill(qc::LINE_DEFAULT);

    dd::Edge result = dd->DDzero;
    dd->incRef(result);

    for (const auto& K : kraus) {
        line[target] = 2;
        dd::Edge Kdd = dd->makeGateDD(K, nqubits, line);
        dd->incRef(Kdd);
        line[target] = qc::LINE_DEFAULT;

        dd::Edge Kdag = dd->conjugateTranspose(Kdd);
        dd->incRef(Kdag);

        dd::Edge tmp  = dd->multiply(Kdd, rho);
        dd->incRef(tmp);
        dd::Edge term = dd->multiply(tmp, Kdag);
        dd->incRef(term);

        dd::Edge newResult = dd->add(result, term);
        dd->incRef(newResult);

        dd->decRef(result);
        dd->decRef(tmp);
        dd->decRef(term);
        dd->decRef(Kdag);
        dd->decRef(Kdd);
        dd->garbageCollect();

        result = newResult;
    }

    return result;
}

// ---------------------------------------------------------------------------
// Standard noise channels
// ---------------------------------------------------------------------------

/// Depolarizing: E(rho) = (1-p)*rho + (p/3)*(X rho X + Y rho Y + Z rho Z)
inline std::vector<KrausOp> depolarizingKraus(double p) {
    double a = std::sqrt(1.0 - p);
    double b = std::sqrt(p / 3.0);
    return {
        KrausOp({dd::ComplexValue{a,0},  dd::ComplexValue{0,0},  dd::ComplexValue{0,0},  dd::ComplexValue{a,0}}),
        KrausOp({dd::ComplexValue{0,0},  dd::ComplexValue{b,0},  dd::ComplexValue{b,0},  dd::ComplexValue{0,0}}),
        KrausOp({dd::ComplexValue{0,0},  dd::ComplexValue{0,-b}, dd::ComplexValue{0,b},  dd::ComplexValue{0,0}}),
        KrausOp({dd::ComplexValue{b,0},  dd::ComplexValue{0,0},  dd::ComplexValue{0,0},  dd::ComplexValue{-b,0}})
    };
}

/// Amplitude damping (T1): gamma = decay probability |1> -> |0>
inline std::vector<KrausOp> amplitudeDampingKraus(double gamma) {
    double sqG  = std::sqrt(gamma);
    double sq1G = std::sqrt(1.0 - gamma);
    return {
        KrausOp({dd::ComplexValue{1,0}, dd::ComplexValue{0,0},   dd::ComplexValue{0,0}, dd::ComplexValue{sq1G,0}}),
        KrausOp({dd::ComplexValue{0,0}, dd::ComplexValue{sqG,0}, dd::ComplexValue{0,0}, dd::ComplexValue{0,0}})
    };
}

/// Dephasing (T2/phase damping): p = dephasing probability
inline std::vector<KrausOp> dephasingKraus(double p) {
    double sqP  = std::sqrt(p);
    double sq1P = std::sqrt(1.0 - p);
    return {
        KrausOp({dd::ComplexValue{1,0}, dd::ComplexValue{0,0}, dd::ComplexValue{0,0}, dd::ComplexValue{sq1P,0}}),
        KrausOp({dd::ComplexValue{0,0}, dd::ComplexValue{0,0}, dd::ComplexValue{0,0}, dd::ComplexValue{sqP,0}})
    };
}

// ---------------------------------------------------------------------------
// Step 2: igGroupSifting on density matrix
// ---------------------------------------------------------------------------

/// Apply igGroupSifting to compress a density matrix DD.
/// varMap maps logical qubit index -> current level position.
/// Returns the reordered DD (same root, updated internally); caller keeps existing ref.
/// Precondition: `rho` must be held with incRef by caller.
inline dd::Edge applyIGGroupSifting(std::unique_ptr<dd::Package>& dd,
                                    dd::Edge rho,
                                    qc::permutationMap& varMap) {
    dd->setMode(dd::Matrix);
    dd::InteractionGraph emptyIG;
    return dd->dynamicReorder(rho, varMap, dd::IGGroupSifting);
}

/// Apply plain Sifting to a density matrix DD (fallback when IG unavailable).
inline dd::Edge applySifting(std::unique_ptr<dd::Package>& dd,
                             dd::Edge rho,
                             qc::permutationMap& varMap) {
    dd->setMode(dd::Matrix);
    return dd->dynamicReorder(rho, varMap, dd::Sifting);
}

// ---------------------------------------------------------------------------
// Observables
// ---------------------------------------------------------------------------

/// Purity = Tr(rho^2).  Pure: 1.0.  Maximally mixed n-qubit: 1/2^n.
inline double purity(std::unique_ptr<dd::Package>& dd, dd::Edge rho) {
    dd->setMode(dd::Matrix);
    dd::Edge rho2 = dd->multiply(rho, rho);
    dd->incRef(rho2);
    dd::ComplexValue tr = dd->trace(rho2);
    dd->decRef(rho2);
    dd->garbageCollect();
    return tr.r;
}

// ---------------------------------------------------------------------------
// Step 3: VQE Hamiltonian expectation value  <H> = Tr(rho * H)
// ---------------------------------------------------------------------------

/// A Pauli term: coefficient * tensor product of single-qubit Paulis.
/// `paulis[q]` is 0=I, 1=X, 2=Y, 3=Z for qubit q.
struct PauliTerm {
    double coeff;
    std::vector<int> paulis;  // length == nqubits
};

namespace {
    // Single-qubit Pauli matrices as GateMatrix (flat row-major)
    inline KrausOp pauliGate(int p) {
        switch (p) {
            case 1: return qc::Xmat;
            case 2: return qc::Ymat;
            case 3: return qc::Zmat;
            default: return qc::Imat;
        }
    }
}

/// Build the DD for a single Pauli tensor product term.
/// Returns a matrix DD; caller must incRef/decRef.
inline dd::Edge buildPauliDD(std::unique_ptr<dd::Package>& dd,
                             const PauliTerm& term,
                             unsigned short nqubits) {
    dd->setMode(dd::Matrix);
    std::array<short, dd::MAXN> line;
    line.fill(qc::LINE_DEFAULT);

    dd::Edge result = dd->makeIdent(0, (short)(nqubits - 1));
    dd->incRef(result);

    for (unsigned short q = 0; q < nqubits; ++q) {
        if (term.paulis[q] == 0) continue;  // identity, skip
        line[q] = 2;
        dd::Edge gate = dd->makeGateDD(pauliGate(term.paulis[q]), nqubits, line);
        dd->incRef(gate);
        line[q] = qc::LINE_DEFAULT;

        dd::Edge newResult = dd->multiply(result, gate);
        dd->incRef(newResult);
        dd->decRef(result);
        dd->decRef(gate);
        dd->garbageCollect();
        result = newResult;
    }
    return result;
}

/// Compute expectation value <H> = Tr(rho * H) for a Hamiltonian
/// given as a sum of Pauli terms.
///
/// H = sum_i coeff_i * P_i   where P_i is an n-qubit Pauli tensor product.
///
/// Returns the real part of Tr(rho * H) (imaginary part is 0 for Hermitian H).
inline double expectationValue(std::unique_ptr<dd::Package>& dd,
                               dd::Edge rho,
                               const std::vector<PauliTerm>& hamiltonian,
                               unsigned short nqubits) {
    dd->setMode(dd::Matrix);
    double result = 0.0;

    for (const auto& term : hamiltonian) {
        // Build Pauli DD
        dd::Edge pauliDD = buildPauliDD(dd, term, nqubits);

        // rhoH = rho * pauliDD
        dd::Edge rhoH = dd->multiply(rho, pauliDD);
        dd->incRef(rhoH);

        // Tr(rhoH)
        dd::ComplexValue tr = dd->trace(rhoH);

        dd->decRef(rhoH);
        dd->decRef(pauliDD);
        dd->garbageCollect();

        result += term.coeff * tr.r;
    }
    return result;
}

}  // namespace dm
