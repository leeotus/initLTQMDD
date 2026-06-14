#ifndef QEC_STEANECODE_HPP
#define QEC_STEANECODE_HPP

#include "QECCode.hpp"

namespace qec {

    // Steane [[7,1,3]] code
    // 7 physical qubits encode 1 logical qubit, distance 3
    // Clifford gates are transversal: H, S, X, Z, CNOT
    // Non-Clifford T gate requires magic state injection (not implemented here)
    class SteaneCode : public QECCode {
    public:
        SteaneCode();

        [[nodiscard]] int nPhysical() const override { return 7; }
        [[nodiscard]] int kLogical() const override { return 1; }
        [[nodiscard]] int distance() const override { return 3; }
        [[nodiscard]] int nAncilla() const override { return 6; }

        [[nodiscard]] const std::vector<StabilizerGenerator>& stabilizers() const override {
            return stabilizers_;
        }

        [[nodiscard]] std::vector<Pauli> logicalX(int logicalIdx = 0) const override;
        [[nodiscard]] std::vector<Pauli> logicalZ(int logicalIdx = 0) const override;

        std::unique_ptr<qc::QuantumComputation> generateEncodingCircuit() const override;
        std::unique_ptr<qc::QuantumComputation> generateSyndromeExtraction() const override;

        // In-place versions that build directly into an existing QuantumComputation
        void generateEncodingCircuitInto(qc::QuantumComputation& qc) const;
        void generateSyndromeExtractionInto(qc::QuantumComputation& qc) const;
        std::unique_ptr<qc::QuantumComputation> generateCorrection(
            const std::vector<int>& syndrome) const override;

        void encodeSingleQubitGate(qc::QuantumComputation& circuit,
                                   qc::OpType gate,
                                   int logicalQubit,
                                   const std::vector<fp>& params = {}) const override;

        void encodeTwoQubitGate(qc::QuantumComputation& circuit,
                                qc::OpType gate,
                                int logicalControl,
                                int logicalTarget,
                                const std::vector<fp>& params = {}) const override;

        std::unique_ptr<qc::QuantumComputation> decodeLogicalState(
            bool zBasis = true) const override;

        int decodeMeasurement(const std::string& physicalResult,
                              bool zBasis = true) const override;

        std::vector<int> extractSyndrome(const std::string& measurementResult) const override;

    private:
        std::vector<StabilizerGenerator> stabilizers_;

        void initStabilizers();
        void addSyndromeExtractionRound(qc::QuantumComputation& qc) const;
    };

} // namespace qec

#endif // QEC_STEANECODE_HPP