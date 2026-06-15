// Benchmark: StochSimulate with different sifting strategies
// Tests the "one-shot IG detection, N-run reuse" optimization.
// No Boost dependency.
//
// Usage: ./bench_stoch_sift <circuit.qasm> [stoch_runs=100] [noise_prob=0.01] [threshold=50]

#include <iostream>
#include <memory>
#include <chrono>
#include <iomanip>
#include <string>
#include <cstdlib>

#include "QuantumComputation.hpp"
#include "QFRSimulator.hpp"

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <circuit_file> [stoch_runs=100] [noise_prob=0.01] [threshold=50]\n";
        std::cerr << "\nCompares StochSimulate performance with different dynamic_reorder strategies:\n";
        std::cerr << "  0 = none (baseline)\n";
        std::cerr << "  1 = plain sifting\n";
        std::cerr << "  4 = group sifting\n";
        std::cerr << "  5 = ig group sifting\n";
        return 1;
    }

    std::string filename = argv[1];
    long stoch_runs = (argc > 2) ? std::atol(argv[2]) : 100;
    double noise_prob = (argc > 3) ? std::atof(argv[3]) : 0.01;
    unsigned int threshold = (argc > 4) ? std::atoi(argv[4]) : 50;

    std::string name = filename.substr(filename.find_last_of("/\\") + 1);

    // noise_effects: "A" = amplitude damping
    std::string noise_effects = "APD";
    // recorded_properties: measure all qubits + fidelity
    std::string recorded_props = "0-" + std::to_string(0) + ";fidelity";

    // First, get qubit count
    {
        auto qc_tmp = std::make_unique<qc::QuantumComputation>(filename);
        unsigned short nq = qc_tmp->getNqubits();
        recorded_props = "0-" + std::to_string(nq - 1) + ";fidelity";
        std::cout << "Circuit: " << name << ", qubits=" << nq << ", gates=" << qc_tmp->getNops() << "\n";
        std::cout << "Noise: " << noise_effects << ", prob=" << noise_prob
                  << ", runs=" << stoch_runs << ", threshold=" << threshold << "\n";
        std::cout << std::string(80, '-') << "\n";
    }

    std::vector<int> strategies = {0, 1, 5};
    const char* strat_names[] = {"none", "sifting", "move_top", "move_bot", "group", "iggroup"};

    std::cout << std::left
              << std::setw(12) << "Strategy"
              << std::setw(12) << "Time(s)"
              << std::setw(12) << "Fidelity"
              << std::setw(12) << "Reorders"
              << "\n";
    std::cout << std::string(48, '-') << "\n";

    for (int strat : strategies) {
        try {
            auto qc = std::make_unique<qc::QuantumComputation>(filename);
            auto sim = std::make_unique<QFRSimulator>(
                qc, noise_effects, noise_prob, stoch_runs,
                1, 1.0, recorded_props,
                0, strat, 0, threshold
            );

            auto start = std::chrono::high_resolution_clock::now();
            auto results = sim->StochSimulate();
            auto end = std::chrono::high_resolution_clock::now();

            double elapsed = std::chrono::duration<double>(end - start).count();
            double fidelity = 0.0;
            if (results.count("fidelity")) {
                fidelity = results["fidelity"];
            }

            auto stats = sim->AdditionalStatistics();

            std::cout << std::left
                      << std::setw(12) << strat_names[strat]
                      << std::setw(12) << std::fixed << std::setprecision(3) << elapsed
                      << std::setw(12) << std::fixed << std::setprecision(6) << fidelity
                      << std::setw(12) << stats["reorder_count"]
                      << "\n";
        } catch (std::exception& e) {
            std::cout << std::left
                      << std::setw(12) << strat_names[strat]
                      << "FAILED: " << e.what() << "\n";
        }
    }

    return 0;
}
