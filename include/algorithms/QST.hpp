#ifndef QST_HPP
#define QST_HPP

#include "QuantumComputation.hpp"
#include "InteractionGraph.h"

#include <vector>
#include <complex>
#include <random>
#include <cmath>
#include <functional>
#include <chrono>
#include <iostream>

namespace qc {

/**
 * @brief Quantum State Tomography (QST) using DD-based representation.
 * 
 * This class implements MLE (Maximum Likelihood Estimation) based QST
 * using the QMDD decision diagram package. The key idea is that the
 * density matrix ρ is represented as a DD, and measurement outcomes
 * are computed via DD operations (trace, partial trace, multiply).
 * 
 * The DD-based approach exploits sparsity in the density matrix,
 * achieving much lower memory usage compared to dense matrix methods
 * like Qiskit Aer or QuTiP.
 */
class QST {
public:
    struct MeasurementBasis {
        enum Basis { Z, X, Y };
        std::vector<Basis> bases;  // one per qubit
    };
    
    struct MeasurementOutcome {
        std::vector<int> bits;     // 0 or 1 per qubit
        unsigned long long count;  // number of times observed
    };
    
    struct QSTConfig {
        unsigned int nQubits;
        unsigned int shotsPerBasis;       // measurements per basis setting
        unsigned int nIterations;         // MLE iterations
        double convergenceTol;            // convergence tolerance
        double learningRate;              // gradient descent step size
        bool usePauliBasis;               // use full Pauli basis (3^n measurements)
        bool useReducedBasis;             // use reduced basis for compressed sensing
        unsigned int reducedBasisSize;    // number of random bases for compressed sensing
        unsigned int seed;                // random seed
        
        QSTConfig(unsigned int n = 2)
            : nQubits(n)
            , shotsPerBasis(1000)
            , nIterations(100)
            , convergenceTol(1e-6)
            , learningRate(0.1)
            , usePauliBasis(false)
            , useReducedBasis(true)
            , reducedBasisSize(0)
            , seed(42) {}
    };
    
    struct QSTResult {
        dd::Edge rho;                           // reconstructed density matrix DD
        dd::Edge pureStateRho;                  // true density matrix (if known)
        double fidelity;                        // fidelity between true and reconstructed
        double traceDistance;                   // trace distance
        double logLikelihood;                   // final log-likelihood
        double elapsedTimeMs;                   // execution time
        unsigned int ddSize;                    // DD node count
        unsigned int peakNodeCount;             // peak node count during reconstruction
        unsigned int convergedIteration;        // iteration at convergence
        bool converged;                         // whether MLE converged
        
        QSTResult() 
            : rho({nullptr, dd::ComplexNumbers::ZERO})
            , pureStateRho({nullptr, dd::ComplexNumbers::ZERO})
            , fidelity(0.0), traceDistance(0.0), logLikelihood(0.0)
            , elapsedTimeMs(0.0), ddSize(0), peakNodeCount(0)
            , convergedIteration(0), converged(false) {}
    };
    
    QST() = default;
    
    /**
     * @brief Generate Pauli basis measurement settings.
     * @param nQubits Number of qubits
     * @param usePauliBasis If true, generate all 3^n bases. If false, generate random subset.
     * @param nRandomBases Number of random bases (only used if usePauliBasis is false)
     * @param seed Random seed
     * @return Vector of measurement bases
     */
    static std::vector<MeasurementBasis> generateMeasurementBases(
        unsigned int nQubits, bool usePauliBasis, 
        unsigned int nRandomBases = 0, unsigned int seed = 42);
    
    /**
     * @brief Simulate measurement outcomes from a true density matrix.
     * @param dd DD package
     * @param rhoTrue True density matrix DD
     * @param config QST configuration
     * @return Vector of measurement outcomes
     */
    static std::vector<MeasurementOutcome> simulateMeasurements(
        std::unique_ptr<dd::Package>& dd, const dd::Edge& rhoTrue,
        const QSTConfig& config);
    
    /**
     * @brief Perform Maximum Likelihood Estimation (MLE) for QST.
     * @param dd DD package
     * @param measurements Measurement outcomes
     * @param bases Measurement bases used
     * @param config QST configuration
     * @return QST result with reconstructed density matrix
     */
    QSTResult performMLE(
        std::unique_ptr<dd::Package>& dd,
        const std::vector<MeasurementOutcome>& measurements,
        const std::vector<MeasurementBasis>& bases,
        const QSTConfig& config);
    
    /**
     * @brief Perform MLE with sifting optimization (Dynamic Reordering).
     * The DD is reordered after each iteration to keep it compact.
     */
    QSTResult performMLEWithReordering(
        std::unique_ptr<dd::Package>& dd,
        const std::vector<MeasurementOutcome>& measurements,
        const std::vector<MeasurementBasis>& bases,
        const QSTConfig& config,
        dd::DynamicReorderingStrategy strat,
        dd::InteractionGraph* ig = nullptr);
    
    /**
     * @brief Compute fidelity between two density matrices.
     */
    static double computeFidelity(
        std::unique_ptr<dd::Package>& dd, 
        const dd::Edge& rho1, const dd::Edge& rho2);
    
    /**
     * @brief Compute trace distance between two density matrices.
     */
    static double computeTraceDistance(
        std::unique_ptr<dd::Package>& dd,
        const dd::Edge& rho1, const dd::Edge& rho2);
    
    /**
     * @brief Create a pure state density matrix from a state vector DD.
     * ρ = |ψ⟩⟨ψ|
     */
    static dd::Edge stateToDensityMatrix(
        std::unique_ptr<dd::Package>& dd, const dd::Edge& stateVector);
    
    /**
     * @brief Create a random pure state density matrix.
     */
    static dd::Edge createRandomPureState(
        std::unique_ptr<dd::Package>& dd, unsigned int nQubits, unsigned int seed = 42);
    
    /**
     * @brief Create a random mixed state density matrix.
     */
    static dd::Edge createRandomMixedState(
        std::unique_ptr<dd::Package>& dd, unsigned int nQubits, 
        double purity = 0.8, unsigned int seed = 42);
    
    /**
     * @brief Create a GHZ state density matrix.
     */
    static dd::Edge createGHZState(
        std::unique_ptr<dd::Package>& dd, unsigned int nQubits);
    
    /**
     * @brief Create a W state density matrix.
     */
    static dd::Edge createWState(
        std::unique_ptr<dd::Package>& dd, unsigned int nQubits);
    
    /**
     * @brief Create a Pauli measurement operator for a single qubit.
     * This is needed by the full QST pipeline for building projectors.
     */
    static dd::Edge makePauliMeasurementDD(
        std::unique_ptr<dd::Package>& dd, unsigned int nQubits,
        unsigned int qubit, MeasurementBasis::Basis basis, int outcome);

    /**
     * @brief Project a density matrix onto a measurement basis and outcome.
     * Returns the probability of observing that outcome.
     */
    static double measureProbability(
        std::unique_ptr<dd::Package>& dd, const dd::Edge& rho,
        const MeasurementBasis& basis, const std::vector<int>& outcome);
    
private:
    /**
     * @brief Compute the gradient of the log-likelihood for MLE.
     */
    static dd::Edge computeGradient(
        std::unique_ptr<dd::Package>& dd, const dd::Edge& rho,
        const std::vector<MeasurementOutcome>& measurements,
        const std::vector<MeasurementBasis>& bases,
        const QSTConfig& config);
    
    /**
     * @brief Apply the RρR iteration for MLE (dilute MLE).
     * R = sum_k (f_k / p_k(rho)) * Π_k
     * ρ_{t+1} = R ρ_t R / Tr(R ρ_t R)
     */
    static dd::Edge diluteMLEIteration(
        std::unique_ptr<dd::Package>& dd, const dd::Edge& rho,
        const std::vector<MeasurementOutcome>& measurements,
        const std::vector<MeasurementBasis>& bases,
        const QSTConfig& config);
};

/**
 * @brief Benchmark runner for QST with different sifting strategies.
 */
class QSTBenchmark {
public:
    struct BenchmarkResult {
        std::string strategyName;
        double avgFidelity;
        double avgTraceDistance;
        double avgLogLikelihood;
        double avgTimeMs;
        double avgDDSize;
        double avgPeakNodes;
        int convergedCount;
        int totalTrials;
        
        void print(std::ostream& os = std::cout) const {
            os << strategyName << ": "
               << "fid=" << avgFidelity 
               << " trDist=" << avgTraceDistance
               << " logL=" << avgLogLikelihood
               << " time=" << avgTimeMs << "ms"
               << " ddSize=" << avgDDSize
               << " peakNodes=" << avgPeakNodes
               << " converged=" << convergedCount << "/" << totalTrials;
        }
        
        std::string toCSV() const {
            std::stringstream ss;
            ss << strategyName << ","
               << avgFidelity << ","
               << avgTraceDistance << ","
               << avgLogLikelihood << ","
               << avgTimeMs << ","
               << avgDDSize << ","
               << avgPeakNodes << ","
               << convergedCount;
            return ss.str();
        }
    };
    
    /**
     * @brief Run QST benchmark comparing multiple strategies.
     */
    static std::vector<BenchmarkResult> runBenchmark(
        const QST::QSTConfig& config,
        unsigned int nTrials = 5,
        bool useDynamic = true);
    
    /**
     * @brief Run QST with a specific density matrix as ground truth.
     */
    static std::vector<BenchmarkResult> runBenchmarkWithState(
        std::unique_ptr<dd::Package>& dd,
        const dd::Edge& rhoTrue,
        const QST::QSTConfig& config,
        unsigned int nTrials = 5);
};

} // namespace qc

#endif // QST_HPP