#ifndef INTERACTION_GRAPH_H
#define INTERACTION_GRAPH_H

#include <vector>
#include <map>
#include <algorithm>
#include <numeric>
#include <unordered_map>
#include <functional>

namespace dd {

    struct InteractionGraph {
        int n = 0;
        std::vector<std::vector<int>> weight;
        std::vector<int> degree;

        // Symmetry detection results
        std::vector<std::vector<short>> symmetricGroups;  // groups of symmetric qubits
        std::vector<short> groupId;                        // qubit -> group index (-1 = singleton)

        InteractionGraph() = default;

        void initForNqubits(int nqubits) {
            n = nqubits;
            weight.assign(n, std::vector<int>(n, 0));
            degree.assign(n, 0);
            symmetricGroups.clear();
            groupId.assign(n, -1);
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

        // Detect symmetric qubit pairs based on IG structure.
        // Two qubits q_i and q_j are IG-symmetric if:
        //   for all k != i,j: w(q_i, k) == w(q_j, k)
        // This means they have identical interaction profiles with all other qubits.
        void detectSymmetry() {
            symmetricGroups.clear();
            groupId.assign(n, -1);

            // Compute interaction profile hash for each qubit
            // Profile of q_i = sorted multiset of (neighbor, weight) excluding self-interactions
            std::unordered_map<size_t, std::vector<short>> profileMap;

            for (int i = 0; i < n; ++i) {
                size_t h = hashProfile(i);
                profileMap[h].push_back(static_cast<short>(i));
            }

            // For each hash bucket, verify exact equality (handle hash collisions)
            for (auto& [hash, candidates] : profileMap) {
                if (candidates.size() < 2) continue;

                std::vector<bool> assigned(candidates.size(), false);
                for (size_t a = 0; a < candidates.size(); ++a) {
                    if (assigned[a]) continue;

                    std::vector<short> group;
                    group.push_back(candidates[a]);
                    assigned[a] = true;

                    for (size_t b = a + 1; b < candidates.size(); ++b) {
                        if (assigned[b]) continue;
                        if (areSymmetric(candidates[a], candidates[b])) {
                            group.push_back(candidates[b]);
                            assigned[b] = true;
                        }
                    }

                    if (group.size() >= 2) {
                        short gid = static_cast<short>(symmetricGroups.size());
                        for (short q : group) {
                            groupId[q] = gid;
                        }
                        symmetricGroups.push_back(std::move(group));
                    }
                }
            }
        }

        // Detect approximate symmetry: qubits whose profiles differ by at most 'tolerance'
        void detectApproximateSymmetry(int tolerance) {
            symmetricGroups.clear();
            groupId.assign(n, -1);

            std::vector<bool> assigned(n, false);
            for (int i = 0; i < n; ++i) {
                if (assigned[i]) continue;

                std::vector<short> group;
                group.push_back(static_cast<short>(i));
                assigned[i] = true;

                for (int j = i + 1; j < n; ++j) {
                    if (assigned[j]) continue;
                    if (profileDistance(i, j) <= tolerance) {
                        group.push_back(static_cast<short>(j));
                        assigned[j] = true;
                    }
                }

                if (group.size() >= 2) {
                    short gid = static_cast<short>(symmetricGroups.size());
                    for (short q : group) {
                        groupId[q] = gid;
                    }
                    symmetricGroups.push_back(std::move(group));
                }
            }
        }

        // Get the representative of each group (highest degree member)
        std::vector<short> getGroupRepresentatives() const {
            std::vector<short> reps;
            for (const auto& group : symmetricGroups) {
                short best = group[0];
                for (short q : group) {
                    if (degree[q] > degree[best]) best = q;
                }
                reps.push_back(best);
            }
            return reps;
        }

        // Check if a qubit is the representative of its group
        bool isRepresentative(short q) const {
            if (q < 0 || q >= n) return true;
            short gid = groupId[q];
            if (gid < 0) return true;  // singleton, always representative
            const auto& group = symmetricGroups[gid];
            short best = group[0];
            for (short member : group) {
                if (degree[member] > degree[best]) best = member;
            }
            return q == best;
        }

        // Get other members in same group (excluding self)
        std::vector<short> getGroupMembers(short q) const {
            std::vector<short> members;
            if (q < 0 || q >= n) return members;
            short gid = groupId[q];
            if (gid < 0) return members;
            for (short member : symmetricGroups[gid]) {
                if (member != q) members.push_back(member);
            }
            return members;
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

    private:
        // Hash a qubit's interaction profile (excluding self-edge)
        size_t hashProfile(int qi) const {
            size_t h = 0;
            for (int k = 0; k < n; ++k) {
                if (k == qi) continue;
                // Use position-independent hashing: hash the sorted (neighbor_degree, weight) pairs
                size_t val = static_cast<size_t>(weight[qi][k]);
                h ^= val * 2654435761ULL + static_cast<size_t>(degree[k]);
            }
            // Also include self-degree as discriminator
            h ^= static_cast<size_t>(degree[qi]) * 31;
            return h;
        }

        // Check exact symmetry between two qubits
        bool areSymmetric(short qi, short qj) const {
            if (degree[qi] != degree[qj]) return false;

            for (int k = 0; k < n; ++k) {
                if (k == qi || k == qj) continue;
                if (weight[qi][k] != weight[qj][k]) return false;
            }
            // Also check: w(qi, qj) should equal w(qj, qi) (always true for undirected)
            // and both should have same self-interaction pattern
            return true;
        }

        // Manhattan distance between two qubit profiles
        int profileDistance(int qi, int qj) const {
            int dist = 0;
            for (int k = 0; k < n; ++k) {
                if (k == qi || k == qj) continue;
                dist += std::abs(weight[qi][k] - weight[qj][k]);
            }
            return dist;
        }
    };

} // namespace dd

#endif // INTERACTION_GRAPH_H
