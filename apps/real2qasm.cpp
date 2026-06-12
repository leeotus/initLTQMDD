#include "QuantumComputation.hpp"
#include <iostream>

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " input.real output.qasm" << std::endl;
        return 1;
    }
    try {
        qc::QuantumComputation qc(argv[1]);
        qc.dump(argv[2], qc::OpenQASM);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
