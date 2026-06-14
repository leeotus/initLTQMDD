#include "QECCode.hpp"

namespace qec {

void QECCode::addGate(qc::QuantumComputation& qc, unsigned short nq,
                      qc::OpType gate, unsigned short target,
                      const std::vector<fp>& params) {
    fp p0 = params.size() > 0 ? params[0] : 0;
    fp p1 = params.size() > 1 ? params[1] : 0;
    fp p2 = params.size() > 2 ? params[2] : 0;
    qc.emplace_back<qc::StandardOperation>(nq, target, gate, p0, p1, p2);
}

void QECCode::addControlledGate(qc::QuantumComputation& qc, unsigned short nq,
                                const std::vector<unsigned short>& controls,
                                unsigned short target,
                                qc::OpType gate) {
    std::vector<qc::Control> ctrls;
    for (auto c : controls) {
        ctrls.emplace_back(c);
    }
    qc.emplace_back<qc::StandardOperation>(nq, ctrls, target, gate);
}

void QECCode::addCNOT(qc::QuantumComputation& qc, unsigned short nq,
                      unsigned short control, unsigned short target) {
    addControlledGate(qc, nq, {control}, target, qc::X);
}

std::unique_ptr<qc::QuantumComputation> QECCode::encodeLogicalCircuit(
    const qc::QuantumComputation& logical,
    int numQecRounds,
    int interleaveGap) const
{
    unsigned short nq = nPhysical() + nAncilla();
    auto physical = std::make_unique<qc::QuantumComputation>(nq);
    // encoded circuit from logical circuit

    // Step 1: Build physical circuit by generating encoding, QEC rounds, and decoding
    // We rebuild using the helper functions directly rather than copying Operations

    {
        auto enc = generateEncodingCircuit();
        for (auto& op : *enc) {
            if (op->isStandardOperation()) {
                auto* sop = dynamic_cast<qc::StandardOperation*>(op.get());
                if (sop) {
                    if (sop->getNcontrols() == 0) {
                        physical->emplace_back<qc::StandardOperation>(
                            nq, sop->getTargets()[0], sop->getType(),
                            sop->getParameter()[0], sop->getParameter()[1], sop->getParameter()[2]);
                    } else if (sop->getNcontrols() == 1 && sop->getNtargets() == 1) {
                        auto ctrls = sop->getControls();
                        physical->emplace_back<qc::StandardOperation>(
                            nq, ctrls, sop->getTargets()[0], sop->getType());
                    }
                }
            } else if (op->isNonUnitaryOperation()) {
                auto* nu = dynamic_cast<qc::NonUnitaryOperation*>(op.get());
                if (nu && nu->getType() == qc::Measure) {
                    physical->emplace_back<qc::NonUnitaryOperation>(
                        nq, nu->getTargets()[0], qc::Measure);
                }
            }
        }
    }

    // Step 2: Encode logical gates with QEC interleaving
    int gateCount = 0;
    int roundsInserted = 0;
    for (auto& op : logical) {
        if (op->isStandardOperation()) {
            auto* sop = dynamic_cast<qc::StandardOperation*>(op.get());
            if (sop && sop->getNtargets() == 1 && sop->getNcontrols() == 0) {
                encodeSingleQubitGate(*physical, sop->getType(), 0);
            }
        }

        gateCount++;

        if (interleaveGap > 0 && gateCount % interleaveGap == 0
            && roundsInserted < numQecRounds) {
            auto synd = generateSyndromeExtraction();
            for (auto& op : *synd) {
                if (op->isStandardOperation()) {
                    auto* sop = dynamic_cast<qc::StandardOperation*>(op.get());
                    if (sop) {
                        if (sop->getNcontrols() == 0) {
                            physical->emplace_back<qc::StandardOperation>(
                                nq, sop->getTargets()[0], sop->getType(),
                                sop->getParameter()[0], sop->getParameter()[1], sop->getParameter()[2]);
                        } else if (sop->getNcontrols() == 1 && sop->getNtargets() == 1) {
                            auto ctrls = sop->getControls();
                            physical->emplace_back<qc::StandardOperation>(
                                nq, ctrls, sop->getTargets()[0], sop->getType());
                        }
                    }
                } else if (op->isNonUnitaryOperation()) {
                    auto* nu = dynamic_cast<qc::NonUnitaryOperation*>(op.get());
                    if (nu && nu->getType() == qc::Measure) {
                        physical->emplace_back<qc::NonUnitaryOperation>(
                            nq, nu->getTargets()[0], qc::Measure);
                    }
                }
            }
            roundsInserted++;
        }
    }

    // Step 3: Remaining QEC rounds
    while (roundsInserted < numQecRounds) {
        auto synd = generateSyndromeExtraction();
        for (auto& op : *synd) {
            if (op->isStandardOperation()) {
                auto* sop = dynamic_cast<qc::StandardOperation*>(op.get());
                if (sop) {
                    if (sop->getNcontrols() == 0) {
                        physical->emplace_back<qc::StandardOperation>(
                            nq, sop->getTargets()[0], sop->getType(),
                            sop->getParameter()[0], sop->getParameter()[1], sop->getParameter()[2]);
                    } else if (sop->getNcontrols() == 1 && sop->getNtargets() == 1) {
                        auto ctrls = sop->getControls();
                        physical->emplace_back<qc::StandardOperation>(
                            nq, ctrls, sop->getTargets()[0], sop->getType());
                    }
                }
            } else if (op->isNonUnitaryOperation()) {
                auto* nu = dynamic_cast<qc::NonUnitaryOperation*>(op.get());
                if (nu && nu->getType() == qc::Measure) {
                    physical->emplace_back<qc::NonUnitaryOperation>(
                        nq, nu->getTargets()[0], qc::Measure);
                }
            }
        }
        roundsInserted++;
    }

    // Step 4: Append decoding
    {
        auto dec = decodeLogicalState(true);
        for (auto& op : *dec) {
            if (op->isStandardOperation()) {
                auto* sop = dynamic_cast<qc::StandardOperation*>(op.get());
                if (sop && sop->getNcontrols() == 0) {
                    physical->emplace_back<qc::StandardOperation>(
                        nq, sop->getTargets()[0], sop->getType(),
                        sop->getParameter()[0], sop->getParameter()[1], sop->getParameter()[2]);
                }
            } else if (op->isNonUnitaryOperation()) {
                auto* nu = dynamic_cast<qc::NonUnitaryOperation*>(op.get());
                if (nu && nu->getType() == qc::Measure) {
                    physical->emplace_back<qc::NonUnitaryOperation>(
                        nq, nu->getTargets()[0], qc::Measure);
                }
            }
        }
    }

    return physical;
}

} // namespace qec