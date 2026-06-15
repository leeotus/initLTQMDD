// Benchmark: compare ddsim simulation with different sifting strategies.
// Generates QFT/Entanglement circuits internally (no Boost dependency).
// Usage: ./bench_sift_strategies [circuit_file] [threshold]

#include <iostream>
#include <memory>
#include <chrono>
#include <iomanip>
#include <string>

#include "QuantumComputation.hpp"
#include "QFRSimulator.hpp"
#include <algorithms/QFT.hpp>
#include <algorithms/Entanglement.hpp>
#include <algorithms/Grover.hpp>

struct BenchResult {
    int strategy;
    double time_s;
    unsigned long final_dd_size;
    unsigned int reorder_count;
    unsigned long exchange_count;
};

BenchResult run_sim(std::unique_ptr<qc::QuantumComputation>& qc, int dyn_reorder) {
    // Re-create simulator each time for clean state
    auto sim = std::make_unique<QFRSimulator>(qc, 1, 1.0, 0, dyn_reorder, 0);

    auto start = std::chrono::high_resolution_clock::now();
    sim->Simulate(0);
    auto end = std::chrono::high_resolution_clock::now();

    double elapsed = std::chrono::duration<double>(end - start).count();
    auto stats = sim->AdditionalStatistics();

    return {
        dyn_reorder,
        elapsed,
        sim->getNodeCount(),
        static_cast<unsigned int>(std::stoul(stats["reorder_count"])),
        static_cast<unsigned long>(std::stoul(stats["exchange_count"]))
    };
}

void print_header() {
    std::cout << std::left
              << std::setw(30) << "Circuit"
              << std::setw(10) << "Strategy"
              << std::setw(12) << "Time(s)"
              << std::setw(12) << "FinalDD"
              << std::setw(10) << "Reorders"
              << std::setw(12) << "Exchanges"
              << "\n";
    std::cout << std::string(86, '-') << "\n";
}

void print_result(const std::string& name, const BenchResult& r) {
    const char* strat_names[] = {"none", "sifting", "move_top", "move_bot", "group", "iggroup"};
    std::cout << std::left
              << std::setw(30) << name
              << std::setw(10) << strat_names[r.strategy]
              << std::setw(12) << std::fixed << std::setprecision(4) << r.time_s
              << std::setw(12) << r.final_dd_size
              << std::setw(10) << r.reorder_count
              << std::setw(12) << r.exchange_count
              << "\n";
}

int main(int argc, char** argv) {
    std::vector<int> strategies = {0, 1, 4, 5};  // none, sifting, group, iggroup

    print_header();

    // --- Built-in algorithm circuits ---
    std::vector<unsigned int> qft_sizes = {10, 14, 18, 20};
    for (auto nq : qft_sizes) {
        std::string name = "QFT_" + std::to_string(nq) + "q";
        for (int strat : strategies) {
            try {
                auto qc = std::make_unique<qc::QFT>(nq);
                std::unique_ptr<qc::QuantumComputation> qc_base = std::move(qc);
                auto r = run_sim(qc_base, strat);
                print_result(name, r);
            } catch (std::exception& e) {
                std::cerr << name << " strat=" << strat << " FAILED: " << e.what() << "\n";
            }
        }
        std::cout << "\n";
    }

    // Entanglement circuits
    std::vector<unsigned int> ent_sizes = {14, 18, 22};
    for (auto nq : ent_sizes) {
        std::string name = "Entangle_" + std::to_string(nq) + "q";
        for (int strat : strategies) {
            try {
                auto qc = std::make_unique<qc::Entanglement>(nq);
                std::unique_ptr<qc::QuantumComputation> qc_base = std::move(qc);
                auto r = run_sim(qc_base, strat);
                print_result(name, r);
            } catch (std::exception& e) {
                std::cerr << name << " strat=" << strat << " FAILED: " << e.what() << "\n";
            }
        }
        std::cout << "\n";
    }

    // --- File-based circuits (if provided) ---
    if (argc > 1) {
        std::string filename = argv[1];
        std::string name = filename.substr(filename.find_last_of("/\\") + 1);
        for (int strat : strategies) {
            try {
                auto qc = std::make_unique<qc::QuantumComputation>(filename);
                std::unique_ptr<qc::QuantumComputation> qc_base = std::move(qc);
                auto r = run_sim(qc_base, strat);
                print_result(name, r);
            } catch (std::exception& e) {
                std::cerr << name << " strat=" << strat << " FAILED: " << e.what() << "\n";
            }
        }
    }

    return 0;
}
