/*
 * Group Sifting for QMDD variable reordering based on IG symmetry detection.
 *
 * Core idea: qubits with identical interaction profiles (IG-symmetric) behave
 * similarly during sifting. We detect symmetric groups, sift only the group
 * representative, then place other members adjacent to it.
 *
 * Two variants:
 *   groupSifting:   symmetry detection + group-aware sifting (no IG direction)
 *   igGroupSifting: symmetry detection + IG-guided direction + LB pruning
 */

#include "DDpackage.h"
#include <algorithm>
#include <numeric>
#include <cassert>

namespace dd {

    // Move variable at position 'from' to position 'to' by repeated swaps
    static void moveVariable(Package& pkg, short from, short to, Edge& in) {
        if (from == to) return;
        if (from > to) {
            while (from > to) {
                pkg.exchangeBaseCase(from, in);
                --from;
            }
        } else {
            while (from < to) {
                pkg.exchangeBaseCase(from + 1, in);
                ++from;
            }
        }
    }

    // Find current position of circuit variable 'var' in varMap
    static short findPosition(short var, const std::map<unsigned short, unsigned short>& varMap) {
        auto it = varMap.find(static_cast<unsigned short>(var));
        if (it != varMap.end()) {
            // varMap maps circuit var -> DD level. We need to find position.
            // Actually varMap[circuitVar] = ddLevel, so position = ddLevel
            return static_cast<short>(it->second);
        }
        return -1;
    }

    // ================== groupSifting ==================
    // 1. Detect symmetric groups via IG
    // 2. For each group, sift the representative through all positions
    // 3. After representative finds optimal, place group members adjacent

    std::tuple<Edge, unsigned int, unsigned int> Package::groupSifting(
        Edge in, std::map<unsigned short, unsigned short>& varMap,
        InteractionGraph& ig)
    {
        if (ig.n <= 0) return sifting(in, varMap);

        ig.detectSymmetry();

        const auto n = static_cast<short>(in.p->v);
        std::vector<bool> free(n, true);

        computeMatrixProperties = Disabled;
        unsigned int total_max = size(in);
        unsigned int total_min = total_max;

        std::vector<bool> placed(ig.n, false);

        for (int iter = 0; iter < n; ++iter) {
            unsigned long max = 0;
            short pos = -1;

            for (short j = 0; j < n; j++) {
                if (free.at(varMap[j]) && active.at(varMap[j]) > max) {
                    max = active.at(varMap[j]);
                    pos = j;
                }
            }
            if (pos < 0) break;

            if (pos < ig.n && placed[pos]) {
                free.at(varMap[pos]) = false;
                continue;
            }

            free.at(varMap[pos]) = false;
            short optimalPos = pos;
            unsigned long min = size(in);

            short currentPos = pos;
            if (currentPos < n / 2) {
                while (currentPos > 0) {
                    exchangeBaseCase(currentPos, in, varMap);
                    auto s = size(in);
                    total_min = std::min(total_min, s);
                    total_max = std::max(total_max, s);
                    --currentPos;
                    if (s < min) { min = s; optimalPos = currentPos; }
                }
                while (currentPos < n - 1) {
                    exchangeBaseCase(currentPos + 1, in, varMap);
                    auto s = size(in);
                    total_min = std::min(total_min, s);
                    total_max = std::max(total_max, s);
                    ++currentPos;
                    if (s < min) { min = s; optimalPos = currentPos; }
                }
            } else {
                while (currentPos < n - 1) {
                    exchangeBaseCase(currentPos + 1, in, varMap);
                    auto s = size(in);
                    total_min = std::min(total_min, s);
                    total_max = std::max(total_max, s);
                    ++currentPos;
                    if (s < min) { min = s; optimalPos = currentPos; }
                }
                while (currentPos > 0) {
                    exchangeBaseCase(currentPos, in, varMap);
                    auto s = size(in);
                    total_min = std::min(total_min, s);
                    total_max = std::max(total_max, s);
                    --currentPos;
                    if (s < min) { min = s; optimalPos = currentPos; }
                }
            }

            while (currentPos > optimalPos) {
                exchangeBaseCase(currentPos, in, varMap); --currentPos;
            }
            while (currentPos < optimalPos) {
                exchangeBaseCase(currentPos + 1, in, varMap); ++currentPos;
            }

            // Group placement
            if (pos < ig.n && ig.groupId[pos] >= 0) {
                auto members = ig.getGroupMembers(pos);
                unsigned long sizeBeforeGroup = size(in);

                for (short member : members) {
                    if (member >= n || placed[member]) continue;

                    auto memberIt = varMap.find((unsigned short)member);
                    if (memberIt == varMap.end()) continue;
                    short memberPos = static_cast<short>(memberIt->second);

                    short repPos = optimalPos;
                    unsigned long bestSize = size(in);
                    short bestTarget = memberPos;

                    // Try above
                    if (repPos + 1 <= n - 1) {
                        short target = repPos + 1;
                        short tmpPos = memberPos;
                        while (tmpPos < target) { exchangeBaseCase(tmpPos + 1, in, varMap); ++tmpPos; }
                        while (tmpPos > target) { exchangeBaseCase(tmpPos, in, varMap); --tmpPos; }
                        auto s = size(in);
                        if (s < bestSize) { bestSize = s; bestTarget = target; }
                        // Undo
                        while (tmpPos < memberPos) { exchangeBaseCase(tmpPos + 1, in, varMap); ++tmpPos; }
                        while (tmpPos > memberPos) { exchangeBaseCase(tmpPos, in, varMap); --tmpPos; }
                    }

                    // Try below
                    if (repPos - 1 >= 0) {
                        short target = repPos - 1;
                        short tmpPos = memberPos;
                        while (tmpPos < target) { exchangeBaseCase(tmpPos + 1, in, varMap); ++tmpPos; }
                        while (tmpPos > target) { exchangeBaseCase(tmpPos, in, varMap); --tmpPos; }
                        auto s = size(in);
                        if (s < bestSize) { bestSize = s; bestTarget = target; }
                        // Undo
                        while (tmpPos < memberPos) { exchangeBaseCase(tmpPos + 1, in, varMap); ++tmpPos; }
                        while (tmpPos > memberPos) { exchangeBaseCase(tmpPos, in, varMap); --tmpPos; }
                    }

                    // Commit only if improvement
                    if (bestTarget != memberPos && bestSize < size(in)) {
                        short tmpPos = memberPos;
                        while (tmpPos < bestTarget) { exchangeBaseCase(tmpPos + 1, in, varMap); ++tmpPos; }
                        while (tmpPos > bestTarget) { exchangeBaseCase(tmpPos, in, varMap); --tmpPos; }
                    }

                    total_min = std::min(total_min, size(in));
                    total_max = std::max(total_max, size(in));
                    placed[member] = true;
                    free.at(varMap[member]) = false;
                }

                placed[pos] = true;
            }
        }

        return {in, total_min, total_max};
    }

    // ================== igGroupSifting ==================
    // Group Sifting with IG-guided direction + LB pruning

    std::tuple<Edge, unsigned int, unsigned int> Package::igGroupSifting(
        Edge in, std::map<unsigned short, unsigned short>& varMap,
        InteractionGraph& ig)
    {
        if (ig.n <= 0) return sifting(in, varMap);

        ig.detectSymmetry();

        const auto n = static_cast<short>(in.p->v);
        std::vector<bool> free(n, true);
        std::map<unsigned short, unsigned short> invVarMap{};
        for (const auto& i : varMap) invVarMap[i.second] = i.first;

        computeMatrixProperties = Disabled;
        unsigned int total_max = size(in);
        unsigned int total_min = total_max;

        std::vector<bool> placed(ig.n, false);

        for (int iter = 0; iter < n; ++iter) {
            // IG-weighted variable selection
            short pos = -1;
            int maxDeg = 1;
            for (int d = 0; d < std::min((int)n, ig.n); ++d)
                if (ig.degree[d] > maxDeg) maxDeg = ig.degree[d];

            unsigned long maxActive = 1;
            for (short j = 0; j < n; ++j) {
                auto it = varMap.find(j);
                if (it != varMap.end()) {
                    auto act = (unsigned long)active.at(it->second);
                    if (act > maxActive) maxActive = act;
                }
            }

            double bestScore = -1.0;
            constexpr double alpha = 0.85;
            for (short j = 0; j < n; j++) {
                auto it = varMap.find(j);
                if (it == varMap.end()) continue;
                if (!free.at(it->second)) continue;
                if (j < ig.n && placed[j]) continue;

                double actNorm = (double)active.at(it->second) / (double)maxActive;
                double degNorm = (j < ig.n) ? (double)ig.degree[j] / (double)maxDeg : 0.0;
                double score = alpha * actNorm + (1.0 - alpha) * degNorm;
                if (score > bestScore) {
                    bestScore = score;
                    pos = j;
                }
            }
            if (pos < 0) break;

            if (pos < ig.n && placed[pos]) {
                free.at(varMap[pos]) = false;
                continue;
            }

            free.at(varMap[pos]) = false;
            short optimalPos = pos;
            unsigned long min = size(in);
            short currentPos = pos;

            // IG-guided direction
            bool upFirst = false;
            if (pos < ig.n) {
                upFirst = ig.shouldSiftUpFirst(pos, currentPos, n, invVarMap);
            } else {
                upFirst = (currentPos >= n / 2);
            }

            // Sift with LB pruning
            if (!upFirst) {
                while (currentPos > 0) {
                    uint64_t lb = computeLowerBoundDown(varMap, currentPos);
                    if (lb > min) break;
                    exchangeBaseCase(currentPos, in, varMap);
                    auto s = size(in);
                    total_min = std::min(total_min, s);
                    total_max = std::max(total_max, s);
                    --currentPos;
                    if (s < min) { min = s; optimalPos = currentPos; }
                }
                while (currentPos < n - 1) {
                    uint64_t lb = computeLowerBoundUp(varMap, currentPos);
                    if (lb > min) break;
                    exchangeBaseCase(currentPos + 1, in, varMap);
                    auto s = size(in);
                    total_min = std::min(total_min, s);
                    total_max = std::max(total_max, s);
                    ++currentPos;
                    if (s < min) { min = s; optimalPos = currentPos; }
                }
            } else {
                while (currentPos < n - 1) {
                    uint64_t lb = computeLowerBoundUp(varMap, currentPos);
                    if (lb > min) break;
                    exchangeBaseCase(currentPos + 1, in, varMap);
                    auto s = size(in);
                    total_min = std::min(total_min, s);
                    total_max = std::max(total_max, s);
                    ++currentPos;
                    if (s < min) { min = s; optimalPos = currentPos; }
                }
                while (currentPos > 0) {
                    uint64_t lb = computeLowerBoundDown(varMap, currentPos);
                    if (lb > min) break;
                    exchangeBaseCase(currentPos, in, varMap);
                    auto s = size(in);
                    total_min = std::min(total_min, s);
                    total_max = std::max(total_max, s);
                    --currentPos;
                    if (s < min) { min = s; optimalPos = currentPos; }
                }
            }

            // Move back to optimal
            while (currentPos > optimalPos) {
                exchangeBaseCase(currentPos, in, varMap);
                --currentPos;
            }
            while (currentPos < optimalPos) {
                exchangeBaseCase(currentPos + 1, in, varMap);
                ++currentPos;
            }

            // Group placement for symmetric members
            if (pos < ig.n && ig.groupId[pos] >= 0) {
                auto members = ig.getGroupMembers(pos);

                for (short member : members) {
                    if (member >= n || placed[member]) continue;

                    auto memberIt = varMap.find((unsigned short)member);
                    if (memberIt == varMap.end()) continue;
                    short memberPos = static_cast<short>(memberIt->second);

                    // Determine target: adjacent to representative's optimal position
                    short repPos = optimalPos;
                    unsigned long curSize = size(in);
                    short bestTarget = memberPos;
                    unsigned long bestSize = curSize;

                    // Try position repPos+1 (above)
                    if (repPos + 1 <= n - 1) {
                        short target = repPos + 1;
                        short tmpPos = memberPos;
                        if (tmpPos < target) {
                            while (tmpPos < target) {
                                exchangeBaseCase(tmpPos + 1, in, varMap);
                                ++tmpPos;
                            }
                        } else if (tmpPos > target) {
                            while (tmpPos > target) {
                                exchangeBaseCase(tmpPos, in, varMap);
                                --tmpPos;
                            }
                        }
                        auto s = size(in);
                        if (s <= bestSize) {
                            bestSize = s;
                            bestTarget = target;
                        }
                        // Undo
                        while (tmpPos < memberPos) {
                            exchangeBaseCase(tmpPos + 1, in, varMap); ++tmpPos;
                        }
                        while (tmpPos > memberPos) {
                            exchangeBaseCase(tmpPos, in, varMap); --tmpPos;
                        }
                    }

                    // Try position repPos-1 (below)
                    if (repPos - 1 >= 0) {
                        short target = repPos - 1;
                        short tmpPos = memberPos;
                        if (tmpPos < target) {
                            while (tmpPos < target) {
                                exchangeBaseCase(tmpPos + 1, in, varMap);
                                ++tmpPos;
                            }
                        } else if (tmpPos > target) {
                            while (tmpPos > target) {
                                exchangeBaseCase(tmpPos, in, varMap);
                                --tmpPos;
                            }
                        }
                        auto s = size(in);
                        if (s < bestSize) {
                            bestSize = s;
                            bestTarget = target;
                        }
                        // Undo
                        while (tmpPos < memberPos) {
                            exchangeBaseCase(tmpPos + 1, in, varMap); ++tmpPos;
                        }
                        while (tmpPos > memberPos) {
                            exchangeBaseCase(tmpPos, in, varMap); --tmpPos;
                        }
                    }

                    // Commit to best target
                    if (bestTarget != memberPos) {
                        if (memberPos < bestTarget) {
                            while (memberPos < bestTarget) {
                                exchangeBaseCase(memberPos + 1, in, varMap);
                                ++memberPos;
                            }
                        } else {
                            while (memberPos > bestTarget) {
                                exchangeBaseCase(memberPos, in, varMap);
                                --memberPos;
                            }
                        }
                    }

                    total_min = std::min(total_min, size(in));
                    total_max = std::max(total_max, size(in));
                    placed[member] = true;
                    free.at(varMap[member]) = false;
                }

                placed[pos] = true;
            }
        }

        return {in, total_min, total_max};
    }

} // namespace dd
