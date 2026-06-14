#ifndef QEC_QECSIMULATOR_HPP
#define QEC_QECSIMULATOR_HPP

#include <memory>
#include <string>
#include <vector>
#include <map>
#include <functional>

#include "QECCode.hpp"
#include "SteaneCode.hpp"
#include "QFRSimulator.hpp"

namespace qec {

    struct NoiseModel {
        std::string type;             // "A"=amplitude, "D"=depolarization, "P"=phase, "APD"=combined
        double prob = 0.0;            // per-gate noise probability
        double idle_prob = 0.0;       // probability for idle qubits (default = prob)
    };

    struct QECSimResult {
        // Main metrics
        double logical_error_rate = 0.0;
        double output_fidelity = 0.0;
        double pseudo_threshold_est = 0.0;

        // Per-round fidelity tracking
        std::vector<double> round_fidelities;

        // DD statistics
        unsigned long max_dd_size = 0;
        unsigned long final_dd_size = 0;
        unsigned long peak_active_nodes = 0;

        // Simulation timing
        double simulation_time_sec = 0.0;
        double build_time_sec = 0.0;

        // Raw data
        int total_shots = 0;
        int logical_errors = 0;
        std::vector<double> per_shot_fidelities;
    };

    struct QECExperimentConfig {
        std::string logical_circuit_path;
        std::string output_prefix;

        // Noise sweep
        NoiseModel base_noise;
        int num_noise_points = 10;
        double noise_min = 1e-5;
        double noise_max = 1e-1;
        bool log_scale = true;

        // QEC rounds sweep
        std::vector<int> qec_rounds = {1, 2, 5, 10};

        // Simulation
        int shots = 1000;
        bool verbose = false;

        // Reordering strategy
        int dynamic_reorder = 1;      // 0=none, 1=sifting, 2=move2top
        int initial_reorder = 0;
    };

    class QECSimulator {
    public:
        explicit QECSimulator(std::unique_ptr<QECCode> code);

        // Run a single simulation
        QECSimResult runSingle(
            const qc::QuantumComputation& logicalCircuit,
            const NoiseModel& noise,
            int numQecRounds,
            int shots,
            int interleaveGap = 0);

        // Run threshold experiment: sweep physical error rate
        std::vector<QECSimResult> runThresholdSweep(
            const qc::QuantumComputation& logicalCircuit,
            const QECExperimentConfig& config);

        // Run robustness experiment: sweep QEC rounds
        std::vector<QECSimResult> runRoundSweep(
            const qc::QuantumComputation& logicalCircuit,
            const QECExperimentConfig& config);

        // Batch run from config
        void runExperiment(const QECExperimentConfig& config);

        // ===== Static analysis (no simulation) =====
        static void analyzeIGSymmetry(const qc::QuantumComputation& circuit);
        static unsigned long predictDDSize(const qc::QuantumComputation& circuit,
                                           int reorderStrategy = 6); // 6 = IGGroupSifting

        // ===== Utilities =====
        static NoiseModel makeDepolarizingNoise(double p);
        static NoiseModel makeAPDNoise(double p);
        static NoiseModel makeCombinedNoise(double p);
        static std::string noiseTypeToString(const NoiseModel& noise);

        // Precomputed symmetries for known codes
        struct CodeSymmetry {
            int nDataQubits;
            std::vector<std::vector<short>> dataGroups;
            std::vector<std::vector<short>> ancillaGroups;
        };
        static CodeSymmetry getSteaneSymmetry();

    private:
        std::unique_ptr<QECCode> code_;

        // Build the physical circuit with noise
        std::unique_ptr<qc::QuantumComputation> buildPhysicalCircuit(
            const qc::QuantumComputation& logicalCircuit,
            int numQecRounds,
            int interleaveGap) const;

        // Create a simple logical circuit for testing
        static std::unique_ptr<qc::QuantumComputation> makeIdentityCircuit(int numLogicalQubits);
    };

} // namespace qec

#endif // QEC_QECSIMULATOR_HPP