#include "SteaneCode.hpp"
#include "operations/StandardOperation.hpp"
#include "operations/NonUnitaryOperation.hpp"

namespace qec {

SteaneCode::SteaneCode() {
    initStabilizers();
}

void SteaneCode::initStabilizers() {
    // Steane [[7,1,3]] has 6 stabilizer generators (3 X-type, 3 Z-type)
    // Qubit numbering: 0..6 = data qubits
    // X-stabilizers (from classical Hamming [7,4,3] code):
    stabilizers_ = {
        // S0 (X-type): X I X I X I X  — qubits {0,2,4,6}
        { {0,2,4,6}, {Pauli::X, Pauli::X, Pauli::X, Pauli::X} },
        // S1 (X-type): I X X I I X X  — qubits {1,2,5,6}
        { {1,2,5,6}, {Pauli::X, Pauli::X, Pauli::X, Pauli::X} },
        // S2 (X-type): I I I X X X X  — qubits {3,4,5,6}
        { {3,4,5,6}, {Pauli::X, Pauli::X, Pauli::X, Pauli::X} },
        // S3 (Z-type): Z I Z I Z I Z  — qubits {0,2,4,6}
        { {0,2,4,6}, {Pauli::Z, Pauli::Z, Pauli::Z, Pauli::Z} },
        // S4 (Z-type): I Z Z I I Z Z  — qubits {1,2,5,6}
        { {1,2,5,6}, {Pauli::Z, Pauli::Z, Pauli::Z, Pauli::Z} },
        // S5 (Z-type): I I I Z Z Z Z  — qubits {3,4,5,6}
        { {3,4,5,6}, {Pauli::Z, Pauli::Z, Pauli::Z, Pauli::Z} },
    };
}

std::vector<Pauli> SteaneCode::logicalX(int) const {
    // X_L = X X X X X X X (transversal X)
    return std::vector<Pauli>(7, Pauli::X);
}

std::vector<Pauli> SteaneCode::logicalZ(int) const {
    // Z_L = Z Z Z Z Z Z Z (transversal Z)
    return std::vector<Pauli>(7, Pauli::Z);
}

std::unique_ptr<qc::QuantumComputation> SteaneCode::generateEncodingCircuit() const {
    // Encode |0⟩_L for Steane code
    // Based on classical Hamming [7,4,3] encoding with Hadamard for sign
    unsigned short nq = 7 + 6; // 7 data + 6 ancilla
    auto qc = std::make_unique<qc::QuantumComputation>(nq);
    // name: steane_encode

    // Step 1: Prepare ancilla |0⟩ (already done by initialization)
    // Step 2: Encode logical |0⟩ state
    // Using method from Gottesman thesis:
    //   H on qubits 0,1,3; then CNOTs to create superposition

    addGate(*qc, nq, qc::H, 0);
    addGate(*qc, nq, qc::H, 1);
    addGate(*qc, nq, qc::H, 3);

    // CNOT gates for classical Hamming [7,4] encoder
    addCNOT(*qc, nq, 0, 2);
    addCNOT(*qc, nq, 0, 4);
    addCNOT(*qc, nq, 1, 2);
    addCNOT(*qc, nq, 1, 5);
    addCNOT(*qc, nq, 3, 4);
    addCNOT(*qc, nq, 3, 5);
    addCNOT(*qc, nq, 2, 6);
    addCNOT(*qc, nq, 4, 6);
    addCNOT(*qc, nq, 5, 6);

    return qc;
}

void SteaneCode::generateEncodingCircuitInto(qc::QuantumComputation& qc) const {
    unsigned short nq = nPhysical() + nAncilla();

    addGate(qc, nq, qc::H, 0);
    addGate(qc, nq, qc::H, 1);
    addGate(qc, nq, qc::H, 3);

    addCNOT(qc, nq, 0, 2);
    addCNOT(qc, nq, 0, 4);
    addCNOT(qc, nq, 1, 2);
    addCNOT(qc, nq, 1, 5);
    addCNOT(qc, nq, 3, 4);
    addCNOT(qc, nq, 3, 5);
    addCNOT(qc, nq, 2, 6);
    addCNOT(qc, nq, 4, 6);
    addCNOT(qc, nq, 5, 6);
}

std::unique_ptr<qc::QuantumComputation> SteaneCode::generateSyndromeExtraction() const {
    unsigned short nq = 7 + 6; // 7 data + 6 ancilla
    auto qc = std::make_unique<qc::QuantumComputation>(nq);
    // name: steane_syndrome

    // For each of the 6 stabilizer generators, measure the stabilizer
    // using a dedicated ancilla qubit (data qubits 0-6, ancilla 7-12)

    for (int s = 0; s < 6; ++s) {
        const auto& stab = stabilizers_[s];
        unsigned short ancilla = 7 + s;

        if (stab.hasX()) {
            // X-type measurement: |+⟩ → CNOTs → H → measure
            addGate(*qc, nq, qc::H, ancilla);
            for (size_t i = 0; i < stab.qubits.size(); ++i) {
                if (stab.paulis[i] == Pauli::X || stab.paulis[i] == Pauli::Y) {
                    addCNOT(*qc, nq, ancilla, stab.qubits[i]);
                }
            }
            addGate(*qc, nq, qc::H, ancilla);
        }

        if (stab.hasZ()) {
            // Z-type measurement: CNOTs → measure
            // (if both X and Z, we use separate halves)
            for (size_t i = 0; i < stab.qubits.size(); ++i) {
                if (stab.paulis[i] == Pauli::Z || stab.paulis[i] == Pauli::Y) {
                    addCNOT(*qc, nq, stab.qubits[i], ancilla);
                }
            }
        }

        // Measure the ancilla
        // Use NonUnitaryOperation (Measure) on the ancilla
        // Note: QEC benchmark simulates this via stochastic unraveling,
        // so we mark measurement with NonUnitaryOperation
        qc->emplace_back<qc::NonUnitaryOperation>(nq, ancilla, qc::Measure);
    }

    return qc;
}

void SteaneCode::generateSyndromeExtractionInto(qc::QuantumComputation& qc) const {
    unsigned short nq = nPhysical() + nAncilla();

    for (int s = 0; s < 6; ++s) {
        const auto& stab = stabilizers_[s];
        unsigned short ancilla = 7 + s;

        if (stab.hasX()) {
            addGate(qc, nq, qc::H, ancilla);
            for (size_t i = 0; i < stab.qubits.size(); ++i) {
                if (stab.paulis[i] == Pauli::X || stab.paulis[i] == Pauli::Y) {
                    addCNOT(qc, nq, ancilla, stab.qubits[i]);
                }
            }
            addGate(qc, nq, qc::H, ancilla);
        }

        if (stab.hasZ()) {
            for (size_t i = 0; i < stab.qubits.size(); ++i) {
                if (stab.paulis[i] == Pauli::Z || stab.paulis[i] == Pauli::Y) {
                    addCNOT(qc, nq, stab.qubits[i], ancilla);
                }
            }
        }

        qc.emplace_back<qc::NonUnitaryOperation>(nq, ancilla, qc::Measure);
    }
}

std::unique_ptr<qc::QuantumComputation> SteaneCode::generateCorrection(
    const std::vector<int>& syndrome) const
{
    unsigned short nq = 7 + 6;
    auto qc = std::make_unique<qc::QuantumComputation>(nq);
    // name: steane_correction

    if (syndrome.size() != 6) return qc;

    // Decode which data qubit has an X error (from syndrome bits 0-2)
    // Decode which data qubit has a Z error (from syndrome bits 3-5)
    // Hamming code: syndrome = column of parity check matrix

    // X error location decoding
    int xErrBit = syndrome[0] * 1 + syndrome[1] * 2 + syndrome[2] * 4;
    if (xErrBit >= 1 && xErrBit <= 7) {
        int qubit = xErrBit - 1; // 0-indexed
        addGate(*qc, nq, qc::X, qubit);
    }

    // Z error location decoding
    int zErrBit = syndrome[3] * 1 + syndrome[4] * 2 + syndrome[5] * 4;
    if (zErrBit >= 1 && zErrBit <= 7) {
        int qubit = zErrBit - 1;
        addGate(*qc, nq, qc::Z, qubit);
    }

    return qc;
}

void SteaneCode::encodeSingleQubitGate(
    qc::QuantumComputation& circuit,
    qc::OpType gate,
    int logicalQubit,
    const std::vector<fp>& params) const
{
    (void)logicalQubit; // k=1, always qubit 0
    unsigned short nq = 7 + 6;

    switch (gate) {
        case qc::X:
            for (int i = 0; i < 7; ++i) addGate(circuit, nq, qc::X, i);
            break;
        case qc::Z:
            for (int i = 0; i < 7; ++i) addGate(circuit, nq, qc::Z, i);
            break;
        case qc::H:
            for (int i = 0; i < 7; ++i) addGate(circuit, nq, qc::H, i);
            break;
        case qc::S:
            // Sdag is transversal for Steane code (not S)
            for (int i = 0; i < 7; ++i) addGate(circuit, nq, qc::Sdag, i);
            break;
        case qc::Sdag:
            for (int i = 0; i < 7; ++i) addGate(circuit, nq, qc::S, i);
            break;
        case qc::T:
        case qc::Tdag:
            // T gate: NOT transversal, would need magic state injection
            // Placeholder: skip for now (Clifford-only circuit)
            break;
        default:
            break;
    }
}

void SteaneCode::encodeTwoQubitGate(
    qc::QuantumComputation& circuit,
    qc::OpType gate,
    int logicalControl,
    int logicalTarget,
    const std::vector<fp>& params) const
{
    (void)params;
    unsigned short nq = 7 + 6;
    (void)logicalControl; (void)logicalTarget;

    if (gate == qc::X) {
        // Transversal CNOT: CNOT(phys_control[i], phys_target[i]) for i=0..6
        // For k=1: control block = qubits 0-6, target block = qubits 0-6 (same block)
        // When doing logical CNOT between two encoded qubits:
        //    CNOT(phys[i], phys[i+7]) for i=0..6
        // Here we only have k=1, so we can't do logical 2-qubit gates between
        // different logical qubits in a single code block.
        // For a 2-logical-qubit setup, you'd need two code blocks.
    }
}

std::unique_ptr<qc::QuantumComputation> SteaneCode::decodeLogicalState(bool zBasis) const {
    unsigned short nq = 7 + 6;
    auto qc = std::make_unique<qc::QuantumComputation>(nq);
    // name: steane_decode

    if (zBasis) {
        // Measure Z_L: Just measure each data qubit in Z, majority vote
        for (int i = 0; i < 7; ++i) {
            qc->emplace_back<qc::NonUnitaryOperation>(nq, i, qc::Measure);
        }
    } else {
        // Measure X_L: H then measure Z, majority vote
        for (int i = 0; i < 7; ++i) {
            addGate(*qc, nq, qc::H, i);
        }
        for (int i = 0; i < 7; ++i) {
            qc->emplace_back<qc::NonUnitaryOperation>(nq, i, qc::Measure);
        }
    }

    return qc;
}

int SteaneCode::decodeMeasurement(const std::string& physicalResult, bool zBasis) const {
    (void)zBasis;
    // Majority vote on 7 data qubit measurements
    // physicalResult is a bitstring, first 7 bits are data qubits
    int count1 = 0;
    for (int i = 0; i < 7 && i < (int)physicalResult.size(); ++i) {
        if (physicalResult[i] == '1') count1++;
    }
    return (count1 >= 4) ? 1 : 0;
}

std::vector<int> SteaneCode::extractSyndrome(const std::string& measurementResult) const {
    // Ancilla qubits 7-12 (indices 7..12 in the bitstring)
    // But measurement results are stored per-qubit in order
    // We expect the result format to match qubit ordering
    std::vector<int> syndrome(6, 0);
    for (int i = 0; i < 6; ++i) {
        int idx = 7 + i;
        if (idx < (int)measurementResult.size()) {
            syndrome[i] = (measurementResult[idx] == '1') ? 1 : 0;
        }
    }
    return syndrome;
}

} // namespace qec