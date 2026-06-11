#ifndef INTERACTION_GRAPH_H
#define INTERACTION_GRAPH_H

#include <vector>
#include <map>
#include <algorithm>
#include <numeric>

namespace dd {

    struct InteractionGraph {
        int n = 0;
        std::vector<std::vector<int>> weight;
        std::vector<int> degree;

        InteractionGraph() = default;

        void initForNqubits(int nqubits) {
            n = nqubits;
            weight.assign(n, std::vector<int>(n, 0));
            degree.assign(n, 0);
        }

        template<typename Operation>
        void addGate(const Operation& op) {
            std::vector<unsigned short> involved;
            for (auto t : op->getTargets())
                involved.push_back(t);
            for (const auto& c : op->getControls())
                involved.push_back(c.qubit);

            if (involved.size() < 2) return;

            for (size_t a = 0; a < involved.size(); ++a) {
                for (size_t b = a + 1; b < involved.size(); ++b) {
                    unsigned short qa = involved[a];
                    unsigned short qb = involved[b];
                    if (qa < (unsigned short)n && qb < (unsigned short)n) {
                        weight[qa][qb]++;
                        weight[qb][qa]++;
                        degree[qa]++;
                        degree[qb]++;
                    }
                }
            }
        }

        template<typename QuantumCircuit>
        void build(const QuantumCircuit& qc) {
            initForNqubits(qc.getNqubits());

            for (const auto& op : qc) {
                addGate(op);
            }
        }

        std::vector<short> getSiftOrder() const {
            std::vector<short> order(n);
            std::iota(order.begin(), order.end(), 0);
            auto& deg = degree;
            std::sort(order.begin(), order.end(), [&deg](short a, short b) {
                return deg[a] > deg[b];
            });
            return order;
        }

        std::vector<short> getHybridSiftOrder(
            const std::vector<unsigned long>& activeNodes,
            const std::map<unsigned short, unsigned short>& varMap,
            double alpha = 0.6) const
        {
            unsigned long maxActive = 1;
            for (auto v : activeNodes)
                if (v > maxActive) maxActive = v;
            int maxDeg = 1;
            for (auto d : degree)
                if (d > maxDeg) maxDeg = d;

            std::vector<short> order(n);
            std::iota(order.begin(), order.end(), 0);

            const auto& deg = degree;
            std::sort(order.begin(), order.end(),
                [&activeNodes, &varMap, &deg, alpha, maxActive, maxDeg](short a, short b) {
                    auto itA = varMap.find(a);
                    auto itB = varMap.find(b);
                    unsigned long activeA = (itA != varMap.end() && itA->second < activeNodes.size())
                        ? activeNodes[itA->second] : 0;
                    unsigned long activeB = (itB != varMap.end() && itB->second < activeNodes.size())
                        ? activeNodes[itB->second] : 0;

                    double scoreA = alpha * (double)activeA / (double)maxActive
                                  + (1.0 - alpha) * (double)deg[a] / (double)maxDeg;
                    double scoreB = alpha * (double)activeB / (double)maxActive
                                  + (1.0 - alpha) * (double)deg[b] / (double)maxDeg;
                    return scoreA > scoreB;
                });
            return order;
        }

        bool shouldSiftUpFirst(short var, short pos, short numVars,
                               const std::map<unsigned short, unsigned short>& invVarMap) const
        {
            double gravityUp = 0, gravityDown = 0;
            for (short j = pos + 1; j < numVars; ++j) {
                auto it = invVarMap.find((unsigned short)j);
                if (it != invVarMap.end() && it->second < (unsigned short)this->n)
                    gravityUp += weight[var][it->second];
            }
            for (short j = 0; j < pos; ++j) {
                auto it = invVarMap.find((unsigned short)j);
                if (it != invVarMap.end() && it->second < (unsigned short)this->n)
                    gravityDown += weight[var][it->second];
            }
            return gravityUp > gravityDown;
        }
    };

} // namespace dd

#endif // INTERACTION_GRAPH_H
