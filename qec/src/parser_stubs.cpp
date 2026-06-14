// Stub implementations for parser import functions
// These are needed to link QuantumComputation but are never called by QEC benchmark
#include "QuantumComputation.hpp"

namespace qc {
    // These are declared in QuantumComputation.hpp but defined in parser .cpp files.
    // We provide empty stubs since QEC benchmark never calls import().
    void QuantumComputation::importReal(std::istream&) { }
    void QuantumComputation::importOpenQASM(std::istream&) { }
    void QuantumComputation::importGRCS(std::istream&) { }
    void QuantumComputation::importTFC(std::istream&) { }
    void QuantumComputation::importQC(std::istream&) { }
}