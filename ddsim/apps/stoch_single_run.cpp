// Single stochastic simulation run for benchmarking.
// Outputs one CSV line: circuit,strategy,stoch_runs,noise_prob,time_s,reorder_count,status
// Usage: ./stoch_single_run <circuit> <strategy:0-5> <stoch_runs> <noise_prob> <threshold>

#include <iostream>
#include <memory>
#include <chrono>
#include <string>
#include <cstdlib>

#include "QuantumComputation.hpp"
#include "QFRSimulator.hpp"

int main(int argc, char** argv) {
    if (argc < 6) {
        std::cerr << "Usage: " << argv[0] << " <circuit> <strategy> <stoch_runs> <noise_prob> <threshold>\n";
        return 1;
    }

    std::string filename = argv[1];
    int strategy = std::atoi(argv[2]);
    long stoch_runs = std::atol(argv[3]);
    double noise_prob = std::atof(argv[4]);
    unsigned int threshold = std::atoi(argv[5]);

    std::string name = filename.substr(filename.find_last_of("/\\") + 1);

    try {
        auto qc = std::make_unique<qc::QuantumComputation>(filename);
        unsigned short nq = qc->getNqubits();

        std::string noise_effects = "APD";
        std::string recorded_props = "0-" + std::to_string(nq - 1) + ";fidelity";

        auto sim = std::make_unique<QFRSimulator>(
            qc, noise_effects, noise_prob, stoch_runs,
            1, 1.0, recorded_props,
            0, strategy, 0, threshold
        );

        auto start = std::chrono::high_resolution_clock::now();
        auto results = sim->StochSimulate();
        auto end = std::chrono::high_resolution_clock::now();

        double elapsed = std::chrono::duration<double>(end - start).count();
        auto stats = sim->AdditionalStatistics();

        // CSV output: circuit,qubits,gates,strategy,stoch_runs,noise_prob,threshold,time_s,reorder_count,exchanges,status
        std::cout << name << ","
                  << nq << ","
                  << qc->getNops() << ","
                  << strategy << ","
                  << stoch_runs << ","
                  << noise_prob << ","
                  << threshold << ","
                  << elapsed << ","
                  << stats["reorder_count"] << ","
                  << stats["exchange_count"] << ","
                  << "ok" << std::endl;
    } catch (std::exception& e) {
        std::cout << name << ",,,," << strategy << "," << stoch_runs << "," << noise_prob << "," << threshold << ",,,,error:" << e.what() << std::endl;
        return 1;
    }
    return 0;
}
