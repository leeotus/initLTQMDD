#ifndef QEC_QECCODE_HPP
#define QEC_QECCODE_HPP

#include <vector>
#include <string>
#include <memory>
#include "QuantumComputation.hpp"
#include "operations/StandardOperation.hpp"
#include "operations/NonUnitaryOperation.hpp"
#include "operations/ClassicControlledOperation.hpp"

namespace qec {

    enum class Pauli { I, X, Y, Z };

    struct StabilizerGenerator {
        std::vector<unsigned short> qubits;   // physical qubits involved
        std::vector<Pauli>         paulis;     // Pauli operator per qubit
        bool hasX() const {
            for (auto p : paulis) if (p == Pauli::X || p == Pauli::Y) return true;
            return false;
        }
        bool hasZ() const {
            for (auto p : paulis) if (p == Pauli::Z || p == Pauli::Y) return true;
            return false;
        }
    };

    class QECCode {
    public:
        virtual ~QECCode() = default;

        // Code parameters
        [[nodiscard]] virtual int nPhysical() const = 0;   // n: physical data qubits
        [[nodiscard]] virtual int kLogical() const = 0;     // k: logical qubits encoded
        [[nodiscard]] virtual int distance() const = 0;     // d: code distance
        [[nodiscard]] virtual int nAncilla() const = 0;     // m = n-k: syndrome ancillae

        // Stabilizer information
        [[nodiscard]] virtual const std::vector<StabilizerGenerator>& stabilizers() const = 0;

        // Logical operators (Pauli strings on physical qubits for each logical)
        [[nodiscard]] virtual std::vector<Pauli> logicalX(int logicalIdx = 0) const = 0;
        [[nodiscard]] virtual std::vector<Pauli> logicalZ(int logicalIdx = 0) const = 0;

        // ===== Core methods =====

        // Generate encoding circuit: |0⟩^(n-k) ancilla ⊗ |0⟩^k logical → encoded state
        // Returns a QuantumComputation that prepares the logical |0⟩ state on physical qubits
        virtual std::unique_ptr<qc::QuantumComputation> generateEncodingCircuit() const = 0;

        // Generate ONE round of syndrome extraction
        // Returns (circuit, mapping: ancilla qubit index → stabilizer index)
        virtual std::unique_ptr<qc::QuantumComputation> generateSyndromeExtraction() const = 0;

        // Generate correction operations for a given syndrome pattern
        // syndrome[i] = measurement result of stabilizer i (0 or 1)
        virtual std::unique_ptr<qc::QuantumComputation> generateCorrection(
            const std::vector<int>& syndrome) const = 0;

        // Encode a single logical gate into physical gate sequence
        virtual void encodeSingleQubitGate(qc::QuantumComputation& circuit,
                                           qc::OpType gate,
                                           int logicalQubit,
                                           const std::vector<fp>& params = {}) const = 0;

        virtual void encodeTwoQubitGate(qc::QuantumComputation& circuit,
                                        qc::OpType gate,
                                        int logicalControl,
                                        int logicalTarget,
                                        const std::vector<fp>& params = {}) const = 0;

        // ===== Higher-level methods =====

        // Encode a full logical circuit into physical circuit with QEC interleaving
        // logical: the logical circuit to encode
        // numQecRounds: number of syndrome extraction rounds to insert
        // interleaveGap: insert a QEC round every N logical gates (0 = only at end)
        virtual std::unique_ptr<qc::QuantumComputation> encodeLogicalCircuit(
            const qc::QuantumComputation& logical,
            int numQecRounds,
            int interleaveGap = 0) const;

        // Decode: measure logical qubits in X or Z basis
        // Returns a circuit that measures the logical state
        virtual std::unique_ptr<qc::QuantumComputation> decodeLogicalState(
            bool zBasis = true) const = 0;

        // Get logical value from a measurement result string
        // physicalResult: bitstring of physical qubit measurements
        // returns: logical bit value (0 or 1)
        virtual int decodeMeasurement(const std::string& physicalResult,
                                      bool zBasis = true) const = 0;

        // Get syndrome bits from measurement result
        virtual std::vector<int> extractSyndrome(const std::string& measurementResult) const = 0;

    protected:
        // Helper: add a single-qubit gate
        static void addGate(qc::QuantumComputation& qc, unsigned short nq,
                           qc::OpType gate, unsigned short target,
                           const std::vector<fp>& params = {});

        // Helper: add a controlled gate (multiple controls, single target)
        static void addControlledGate(qc::QuantumComputation& qc, unsigned short nq,
                                      const std::vector<unsigned short>& controls,
                                      unsigned short target,
                                      qc::OpType gate = qc::X);

        // Helper: add a CNOT with single control
        static void addCNOT(qc::QuantumComputation& qc, unsigned short nq,
                           unsigned short control, unsigned short target);
    };

} // namespace qec

#endif // QEC_QECCODE_HPP