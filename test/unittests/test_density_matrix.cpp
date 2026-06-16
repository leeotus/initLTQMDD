#include "gtest/gtest.h"
#include "QuantumComputation.hpp"
#include "DensityMatrix.hpp"
#include <cmath>
#include <memory>

// Tolerance for floating-point comparisons
static constexpr double TOL = 1e-6;

class DensityMatrixTest : public ::testing::Test {
protected:
    void SetUp() override {
        dd = std::make_unique<dd::Package>();
    }
    void TearDown() override {
        dd->garbageCollect(true);
    }
    std::unique_ptr<dd::Package> dd;
};

// ---------------------------------------------------------------------------
// Test 1: rho = |0><0|  (pure state, 1 qubit)
//   Tr(rho) = 1, Purity = 1
// ---------------------------------------------------------------------------
TEST_F(DensityMatrixTest, PureZeroState1Qubit) {
    unsigned short n = 1;
    dd->setMode(dd::Vector);
    dd::Edge psi = dd->makeZeroState(n);  // |0>
    dd->incRef(psi);

    dd->setMode(dd::Matrix);
    dd::Edge rho = dm::densityMatrixFromState(dd, psi);

    // Tr(rho) should be 1
    dd::ComplexValue tr = dd->trace(rho);
    EXPECT_NEAR(tr.r, 1.0, TOL);
    EXPECT_NEAR(tr.i, 0.0, TOL);

    // Purity Tr(rho^2) = 1 for pure state
    double pur = dm::purity(dd, rho);
    EXPECT_NEAR(pur, 1.0, TOL);

    dd->decRef(rho);
}

// ---------------------------------------------------------------------------
// Test 2: |+> = (|0>+|1>)/sqrt(2), rho = |+><+|
//   Tr(rho) = 1, Purity = 1
// ---------------------------------------------------------------------------
TEST_F(DensityMatrixTest, PurePlusState1Qubit) {
    unsigned short n = 1;
    dd->setMode(dd::Vector);
    // Build |+> = H|0>
    dd::Edge zero = dd->makeZeroState(n);
    dd->incRef(zero);

    std::array<short, dd::MAXN> line;
    line.fill(qc::LINE_DEFAULT);
    line[0] = 2;
    dd->setMode(dd::Matrix);
    dd::Edge H = dd->makeGateDD(qc::Hmat, n, line);
    dd->incRef(H);
    line[0] = qc::LINE_DEFAULT;

    dd->setMode(dd::Vector);
    dd::Edge plus = dd->multiply(H, zero);
    dd->incRef(plus);
    dd->decRef(H);
    dd->decRef(zero);

    dd->setMode(dd::Matrix);
    dd::Edge rho = dm::densityMatrixFromState(dd, plus);

    dd::ComplexValue tr = dd->trace(rho);
    EXPECT_NEAR(tr.r, 1.0, TOL);
    EXPECT_NEAR(tr.i, 0.0, TOL);

    double pur = dm::purity(dd, rho);
    EXPECT_NEAR(pur, 1.0, TOL);

    dd->decRef(rho);
    dd->decRef(plus);
    dd->garbageCollect();
}

// ---------------------------------------------------------------------------
// Test 3: Depolarizing channel on |0><0| with p=0
//   Should leave rho unchanged: Tr=1, Purity=1
// ---------------------------------------------------------------------------
TEST_F(DensityMatrixTest, DepolarizingZeroNoise) {
    unsigned short n = 1;
    dd->setMode(dd::Vector);
    dd::Edge psi = dd->makeZeroState(n);
    dd->incRef(psi);

    dd->setMode(dd::Matrix);
    dd::Edge rho = dm::densityMatrixFromState(dd, psi);
    dd->incRef(rho);

    auto kraus = dm::depolarizingKraus(0.0);  // p=0 -> identity channel
    dd::Edge rho2 = dm::applyKrausChannel(dd, rho, kraus, 0, n);

    dd::ComplexValue tr = dd->trace(rho2);
    EXPECT_NEAR(tr.r, 1.0, TOL);
    EXPECT_NEAR(tr.i, 0.0, TOL);

    double pur = dm::purity(dd, rho2);
    EXPECT_NEAR(pur, 1.0, TOL);

    dd->decRef(rho2);
    dd->decRef(rho);
    dd->garbageCollect();
}

// ---------------------------------------------------------------------------
// Test 4: Depolarizing channel on |0><0| with p=3/4 -> maximally mixed I/2
//   Purity = 1/2
// ---------------------------------------------------------------------------
TEST_F(DensityMatrixTest, DepolarizingMaxMixedState) {
    unsigned short n = 1;
    dd->setMode(dd::Vector);
    dd::Edge psi = dd->makeZeroState(n);
    dd->incRef(psi);

    dd->setMode(dd::Matrix);
    dd::Edge rho = dm::densityMatrixFromState(dd, psi);
    dd->incRef(rho);

    // p=3/4: E(rho) = (1-3/4)*rho + (3/4)/3*(XrhoX+YrhoY+ZrhoZ)
    //               = rho/4 + (rho_x + rho_y + rho_z)/4 = I/2
    auto kraus = dm::depolarizingKraus(0.75);
    dd::Edge rho2 = dm::applyKrausChannel(dd, rho, kraus, 0, n);

    dd::ComplexValue tr = dd->trace(rho2);
    EXPECT_NEAR(tr.r, 1.0, TOL);

    double pur = dm::purity(dd, rho2);
    EXPECT_NEAR(pur, 0.5, TOL);  // Tr((I/2)^2) = Tr(I/4) = 1/2

    dd->decRef(rho2);
    dd->decRef(rho);
    dd->garbageCollect();
}

// ---------------------------------------------------------------------------
// Test 5: Amplitude damping on |1><1| with gamma=1
//   |1> -> |0> with probability 1 -> rho' = |0><0|, Purity=1
// ---------------------------------------------------------------------------
TEST_F(DensityMatrixTest, AmplitudeDampingFullDecay) {
    unsigned short n = 1;
    dd->setMode(dd::Vector);
    std::bitset<dd::MAXN> bits;
    bits.set(0);
    dd::Edge psi = dd->makeBasisState(n, bits);  // |1>
    dd->incRef(psi);

    dd->setMode(dd::Matrix);
    dd::Edge rho = dm::densityMatrixFromState(dd, psi);
    dd->incRef(rho);

    auto kraus = dm::amplitudeDampingKraus(1.0);
    dd::Edge rho2 = dm::applyKrausChannel(dd, rho, kraus, 0, n);

    dd::ComplexValue tr = dd->trace(rho2);
    EXPECT_NEAR(tr.r, 1.0, TOL);

    double pur = dm::purity(dd, rho2);
    EXPECT_NEAR(pur, 1.0, TOL);  // |0><0| is pure

    // Size should be minimal (pure state)
    unsigned int sz = dd->size(rho2);
    EXPECT_LE(sz, 4u);

    dd->decRef(rho2);
    dd->decRef(rho);
    dd->garbageCollect();
}

// ---------------------------------------------------------------------------
// Test 6: Multi-qubit — 2-qubit |00>, depolarizing on qubit 0
//   Tr(rho)=1 still holds after noise
// ---------------------------------------------------------------------------
TEST_F(DensityMatrixTest, TwoQubitDepolarizing) {
    unsigned short n = 2;
    dd->setMode(dd::Vector);
    dd::Edge psi = dd->makeZeroState(n);
    dd->incRef(psi);

    dd->setMode(dd::Matrix);
    dd::Edge rho = dm::densityMatrixFromState(dd, psi);
    dd->incRef(rho);

    auto kraus = dm::depolarizingKraus(0.1);
    dd::Edge rho2 = dm::applyKrausChannel(dd, rho, kraus, 0, n);

    dd::ComplexValue tr = dd->trace(rho2);
    EXPECT_NEAR(tr.r, 1.0, TOL);
    EXPECT_NEAR(tr.i, 0.0, TOL);

    dd->decRef(rho2);
    dd->decRef(rho);
    dd->garbageCollect();
}

// ---------------------------------------------------------------------------
// Test 7: DD size sanity — 4-qubit noisy density matrix
//   Build rho for |++++>, apply depolarizing on all qubits, verify Tr=1
//   and report DD node count (baseline for future igGroupSifting experiments).
// ---------------------------------------------------------------------------
TEST_F(DensityMatrixTest, NoisyDensityMatrixSizeReport) {
    unsigned short n = 4;
    dd->setMode(dd::Vector);

    dd::Edge state = dd->makeZeroState(n);
    dd->incRef(state);

    std::array<short, dd::MAXN> line;
    line.fill(qc::LINE_DEFAULT);
    for (unsigned short q = 0; q < n; ++q) {
        line[q] = 2;
        dd->setMode(dd::Matrix);
        dd::Edge H = dd->makeGateDD(qc::Hmat, n, line);
        dd->incRef(H);
        line[q] = qc::LINE_DEFAULT;
        dd->setMode(dd::Vector);
        dd::Edge newState = dd->multiply(H, state);
        dd->incRef(newState);
        dd->decRef(H);
        dd->decRef(state);
        dd->garbageCollect();
        state = newState;
    }

    dd->setMode(dd::Matrix);
    dd::Edge rho = dm::densityMatrixFromState(dd, state);
    dd->incRef(rho);

    unsigned int sizePure = dd->size(rho);

    auto kraus = dm::depolarizingKraus(0.05);
    for (unsigned short q = 0; q < n; ++q) {
        dd::Edge rhoNew = dm::applyKrausChannel(dd, rho, kraus, q, n);
        dd->incRef(rhoNew);
        dd->decRef(rho);
        dd->garbageCollect();
        rho = rhoNew;
    }

    unsigned int sizeNoisy = dd->size(rho);
    dd::ComplexValue tr = dd->trace(rho);

    std::cout << "[NoisyDensityMatrixSizeReport] n=" << n
              << " pure_size=" << sizePure
              << " noisy_size=" << sizeNoisy << "\n";

    EXPECT_NEAR(tr.r, 1.0, TOL);
    EXPECT_NEAR(tr.i, 0.0, TOL);

    dd->decRef(rho);
    dd->garbageCollect();
}

// ===========================================================================
// Step 1: 8-qubit scale-up + memory comparison baseline
// ===========================================================================

// Build helper: apply H to all qubits of state
static dd::Edge applyHAll(std::unique_ptr<dd::Package>& ddpkg,
                           dd::Edge state, unsigned short n) {
    std::array<short, dd::MAXN> line;
    line.fill(qc::LINE_DEFAULT);
    for (unsigned short q = 0; q < n; ++q) {
        line[q] = 2;
        ddpkg->setMode(dd::Matrix);
        dd::Edge H = ddpkg->makeGateDD(qc::Hmat, n, line);
        ddpkg->incRef(H);
        line[q] = qc::LINE_DEFAULT;
        ddpkg->setMode(dd::Vector);
        dd::Edge newState = ddpkg->multiply(H, state);
        ddpkg->incRef(newState);
        ddpkg->decRef(H);
        ddpkg->decRef(state);
        ddpkg->garbageCollect();
        state = newState;
    }
    return state;
}

TEST_F(DensityMatrixTest, EightQubitScaleUp) {
    unsigned short n = 8;
    dd->setMode(dd::Vector);
    dd::Edge state = dd->makeZeroState(n);
    dd->incRef(state);
    state = applyHAll(dd, state, n);

    dd->setMode(dd::Matrix);
    dd::Edge rho = dm::densityMatrixFromState(dd, state);
    dd->incRef(rho);
    unsigned int sizePure = dd->size(rho);

    auto kraus = dm::depolarizingKraus(0.02);
    for (unsigned short q = 0; q < n; ++q) {
        dd::Edge rhoNew = dm::applyKrausChannel(dd, rho, kraus, q, n);
        dd->incRef(rhoNew);
        dd->decRef(rho);
        dd->garbageCollect();
        rho = rhoNew;
    }
    unsigned int sizeNoisy = dd->size(rho);

    // Dense comparison: 4^n complex = 2 * 4^n * 8 bytes
    double denseMB = 2.0 * (1ULL << (2*n)) * 8.0 / (1024.0 * 1024.0);
    std::cout << "[EightQubitScaleUp] n=" << n
              << " pure_size=" << sizePure
              << " noisy_size=" << sizeNoisy
              << " dense_rho_MB=" << denseMB << "\n";

    dd::ComplexValue tr = dd->trace(rho);
    EXPECT_NEAR(tr.r, 1.0, TOL);
    EXPECT_GT(denseMB, 0.5);  // dense should be >0.5MB

    dd->decRef(rho);
    dd->garbageCollect();
}

// ===========================================================================
// Step 2: igGroupSifting compression on density matrix
// ===========================================================================

TEST_F(DensityMatrixTest, IGGroupSiftingCompression) {
    unsigned short n = 6;
    dd->setMode(dd::Vector);
    dd::Edge state = dd->makeZeroState(n);
    dd->incRef(state);
    state = applyHAll(dd, state, n);

    dd->setMode(dd::Matrix);
    dd::Edge rho = dm::densityMatrixFromState(dd, state);
    dd->incRef(rho);

    auto kraus = dm::depolarizingKraus(0.05);
    for (unsigned short q = 0; q < n; ++q) {
        dd::Edge rhoNew = dm::applyKrausChannel(dd, rho, kraus, q, n);
        dd->incRef(rhoNew);
        dd->decRef(rho);
        dd->garbageCollect();
        rho = rhoNew;
    }

    unsigned int sizeBefore = dd->size(rho);
    dd::ComplexValue trBefore = dd->trace(rho);
    EXPECT_NEAR(trBefore.r, 1.0, TOL);

    // Apply igGroupSifting
    qc::permutationMap varMap;
    for (unsigned short q = 0; q < n; ++q) varMap[q] = q;
    rho = dm::applyIGGroupSifting(dd, rho, varMap);

    unsigned int sizeAfter = dd->size(rho);
    dd::ComplexValue trAfter = dd->trace(rho);

    std::cout << "[IGGroupSiftingCompression] n=" << n
              << " before=" << sizeBefore
              << " after=" << sizeAfter
              << " ratio=" << (sizeBefore > 0 ? (double)sizeAfter/sizeBefore : 1.0) << "\n";

    EXPECT_NEAR(trAfter.r, 1.0, TOL);

    dd->decRef(rho);
    dd->garbageCollect();
}

// ===========================================================================
// Step 3: VQE expectation value  <H> = Tr(rho * H)
// ===========================================================================

// Test 3a: <Z> for |0><0|  = +1
TEST_F(DensityMatrixTest, ExpValZOnZeroState) {
    unsigned short n = 1;
    dd->setMode(dd::Vector);
    dd::Edge psi = dd->makeZeroState(n);
    dd->incRef(psi);
    dd->setMode(dd::Matrix);
    dd::Edge rho = dm::densityMatrixFromState(dd, psi);
    dd->incRef(rho);

    // H = Z on qubit 0
    std::vector<dm::PauliTerm> H = { {1.0, {3}} };  // 1.0 * Z
    double ev = dm::expectationValue(dd, rho, H, n);
    EXPECT_NEAR(ev, 1.0, TOL);  // <0|Z|0> = 1

    dd->decRef(rho);
    dd->garbageCollect();
}

// Test 3b: <Z> for |1><1|  = -1
TEST_F(DensityMatrixTest, ExpValZOnOneState) {
    unsigned short n = 1;
    dd->setMode(dd::Vector);
    std::bitset<dd::MAXN> bits;
    bits.set(0);
    dd::Edge psi = dd->makeBasisState(n, bits);
    dd->incRef(psi);
    dd->setMode(dd::Matrix);
    dd::Edge rho = dm::densityMatrixFromState(dd, psi);
    dd->incRef(rho);

    std::vector<dm::PauliTerm> H = { {1.0, {3}} };  // Z
    double ev = dm::expectationValue(dd, rho, H, n);
    EXPECT_NEAR(ev, -1.0, TOL);  // <1|Z|1> = -1

    dd->decRef(rho);
    dd->garbageCollect();
}

// Test 3c: <X> for |+><+|  = +1
TEST_F(DensityMatrixTest, ExpValXOnPlusState) {
    unsigned short n = 1;
    dd->setMode(dd::Vector);
    dd::Edge zero = dd->makeZeroState(n);
    dd->incRef(zero);
    std::array<short, dd::MAXN> line;
    line.fill(qc::LINE_DEFAULT);
    line[0] = 2;
    dd->setMode(dd::Matrix);
    dd::Edge H = dd->makeGateDD(qc::Hmat, n, line);
    dd->incRef(H);
    line[0] = qc::LINE_DEFAULT;
    dd->setMode(dd::Vector);
    dd::Edge plus = dd->multiply(H, zero);
    dd->incRef(plus);
    dd->decRef(H);
    dd->decRef(zero);
    dd->setMode(dd::Matrix);
    dd::Edge rho = dm::densityMatrixFromState(dd, plus);
    dd->incRef(rho);
    dd->decRef(plus);
    dd->garbageCollect();

    std::vector<dm::PauliTerm> Hx = { {1.0, {1}} };  // X
    double ev = dm::expectationValue(dd, rho, Hx, n);
    EXPECT_NEAR(ev, 1.0, TOL);  // <+|X|+> = 1

    dd->decRef(rho);
    dd->garbageCollect();
}

// Test 3d: 2-qubit Ising Hamiltonian H = Z⊗Z + 0.5*X⊗I + 0.5*I⊗X
//   on |00>: <Z⊗Z>=1, <X⊗I>=0, <I⊗X>=0  ->  <H>=1
TEST_F(DensityMatrixTest, ExpValIsingHamiltonian) {
    unsigned short n = 2;
    dd->setMode(dd::Vector);
    dd::Edge psi = dd->makeZeroState(n);
    dd->incRef(psi);
    dd->setMode(dd::Matrix);
    dd::Edge rho = dm::densityMatrixFromState(dd, psi);
    dd->incRef(rho);

    std::vector<dm::PauliTerm> Hising = {
        {1.0,  {3, 3}},   // Z⊗Z
        {0.5,  {1, 0}},   // X⊗I
        {0.5,  {0, 1}}    // I⊗X
    };
    double ev = dm::expectationValue(dd, rho, Hising, n);
    // |00>: Z0Z1=1, X0=0, X1=0  ->  <H>=1
    EXPECT_NEAR(ev, 1.0, TOL);

    dd->decRef(rho);
    dd->garbageCollect();
}

// Test 3e: noisy state -> <Z> decreases with depolarizing noise p
//   E_p(|0><0|) = (1-2p/3)|0><0| + (2p/3)|1><1|
//   <Z> = 1 - 4p/3
TEST_F(DensityMatrixTest, ExpValDecaysWithNoise) {
    unsigned short n = 1;
    double p = 0.3;

    dd->setMode(dd::Vector);
    dd::Edge psi = dd->makeZeroState(n);
    dd->incRef(psi);
    dd->setMode(dd::Matrix);
    dd::Edge rho = dm::densityMatrixFromState(dd, psi);
    dd->incRef(rho);

    auto kraus = dm::depolarizingKraus(p);
    dd::Edge rhoNoisy = dm::applyKrausChannel(dd, rho, kraus, 0, n);
    dd->incRef(rhoNoisy);
    dd->decRef(rho);
    dd->garbageCollect();

    std::vector<dm::PauliTerm> Hz = { {1.0, {3}} };
    double ev = dm::expectationValue(dd, rhoNoisy, Hz, n);
    double expected = 1.0 - 4.0 * p / 3.0;
    EXPECT_NEAR(ev, expected, 1e-5);

    dd->decRef(rhoNoisy);
    dd->garbageCollect();
}
