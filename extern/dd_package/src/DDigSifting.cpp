/*
 * Interaction Graph driven Sifting for QMDD variable reordering.
 *
 * igSifting:   IG variable order + IG direction + no pruning (same result as sifting)
 * igLbSifting: IG variable order + IG direction + LB pruning (key improvement)
 *
 * The value of IG is in guiding LB pruning direction. Standard LB Sifting uses
 * pos < n/2 to decide direction, which can cause the pruning to skip the better
 * direction. IG gravity ensures we explore the direction with more interaction
 * partners first, finding better positions before LB prunes.
 */

#include "DDpackage.h"
#include <algorithm>
#include <numeric>

namespace dd {

    // IG-driven variable selection: hybrid score
    static short selectNextVariable(
        short n, const std::vector<bool>& free,
        const std::map<unsigned short, unsigned short>& varMap,
        const int* activeArr,
        const InteractionGraph& ig)
    {
        short pos = -1;

        if (ig.n > 0) {
            int maxDeg = 1;
            for (int d = 0; d < std::min((int)n, ig.n); ++d)
                if (ig.degree[d] > maxDeg) maxDeg = ig.degree[d];

            unsigned long maxActive = 1;
            for (short j = 0; j < n; ++j) {
                auto it = varMap.find(j);
                if (it != varMap.end()) {
                    auto act = (unsigned long)activeArr[it->second];
                    if (act > maxActive) maxActive = act;
                }
            }

            double bestScore = -1.0;
            for (short j = 0; j < n; j++) {
                auto it = varMap.find(j);
                if (it == varMap.end()) continue;
                if (free.at(it->second) && activeArr[it->second] > 0) {
                    double actScore = (double)activeArr[it->second] / (double)maxActive;
                    double degScore = (j < ig.n) ? (double)ig.degree[j] / (double)maxDeg : 0.0;
                    double score = 0.85 * actScore + 0.15 * degScore;
                    if (score > bestScore) {
                        bestScore = score;
                        pos = j;
                    }
                }
            }
        }

        // Fallback: max active
        if (pos < 0) {
            unsigned long max = 0;
            for (short j = 0; j < n; j++) {
                auto it = varMap.find(j);
                if (it == varMap.end()) continue;
                if (free.at(it->second) && activeArr[it->second] > max) {
                    max = activeArr[it->second];
                    pos = j;
                }
            }
        }
        return pos;
    }

    // IG gravity: should we sift up first?
    static bool shouldSiftUpFirst(
        short circuitVar, short pos, short n,
        const std::map<unsigned short, unsigned short>& invVarMap,
        const InteractionGraph& ig)
    {
        if (ig.n <= 0 || circuitVar >= ig.n) {
            return pos >= n / 2;  // standard heuristic
        }

        double gravityUp = 0, gravityDown = 0;
        for (short j = pos + 1; j < n; ++j) {
            auto it = invVarMap.find(j);
            if (it != invVarMap.end() && it->second < (unsigned short)ig.n)
                gravityUp += ig.weight[circuitVar][it->second];
        }
        for (short j = 0; j < pos; ++j) {
            auto it = invVarMap.find(j);
            if (it != invVarMap.end() && it->second < (unsigned short)ig.n)
                gravityDown += ig.weight[circuitVar][it->second];
        }

        if (gravityUp != gravityDown)
            return gravityUp > gravityDown;
        return pos >= n / 2;  // tie-break: standard heuristic
    }

    // ================== igSifting (no pruning) ==================

    std::tuple<Edge, unsigned int, unsigned int> Package::igSifting(
        Edge in, std::map<unsigned short, unsigned short>& varMap,
        const InteractionGraph& ig)
    {
        if (ig.n <= 0) return sifting(in, varMap);

        const auto n = static_cast<short>(in.p->v);
        std::vector<bool> free(n, true);
        std::map<unsigned short, unsigned short> invVarMap{};
        for (const auto& i : varMap) invVarMap[i.second] = i.first;

        computeMatrixProperties = Disabled;
        Edge root{in};
        unsigned int total_max = size(in);
        unsigned int total_min = total_max;

        short pos = -1;
        for (int i = 0; i < n; ++i) {
            assert(is_globally_consistent_dd(in));
            unsigned long min = size(in);

            pos = selectNextVariable(n, free, varMap, active.data(), ig);
            if (pos < 0) break;

            free.at(varMap[pos]) = false;
            short optimalPos = pos;
            short originalPos = pos;

            bool upFirst = shouldSiftUpFirst(pos, pos, n, invVarMap, ig);

            if (!upFirst) {
                while (pos > 0) {
                    exchangeBaseCase(pos, in);
                    auto s = size(in); total_min = std::min(total_min, s); total_max = std::max(total_max, s);
                    assert(is_locally_consistent_dd(in)); --pos;
                    if (s < min) { min = s; optimalPos = pos; }
                }
                while (pos < n - 1) {
                    exchangeBaseCase(pos + 1, in);
                    auto s = size(in); total_min = std::min(total_min, s); total_max = std::max(total_max, s);
                    assert(is_locally_consistent_dd(in)); ++pos;
                    if (s < min) { min = s; optimalPos = pos; }
                }
                while (pos > optimalPos) {
                    exchangeBaseCase(pos, in);
                    auto s = size(in); total_min = std::min(total_min, s); total_max = std::max(total_max, s);
                    assert(is_locally_consistent_dd(in)); --pos;
                }
            } else {
                while (pos < n - 1) {
                    exchangeBaseCase(pos + 1, in);
                    auto s = size(in); total_min = std::min(total_min, s); total_max = std::max(total_max, s);
                    assert(is_locally_consistent_dd(in)); ++pos;
                    if (s < min) { min = s; optimalPos = pos; }
                }
                while (pos > 0) {
                    exchangeBaseCase(pos, in); assert(is_locally_consistent_dd(in));
                    auto s = size(in); total_min = std::min(total_min, s); total_max = std::max(total_max, s);
                    --pos;
                    if (s < min) { min = s; optimalPos = pos; }
                }
                while (pos < optimalPos) {
                    exchangeBaseCase(pos + 1, in);
                    auto s = size(in); total_min = std::min(total_min, s); total_max = std::max(total_max, s);
                    assert(is_locally_consistent_dd(in)); ++pos;
                }
            }

            initComputeTable();
            if (unnormalizedNodes > 0) {
                auto oldroot = root; root = renormalize(root);
                decRef(oldroot); incRef(root); in.p = root.p; in.w = root.w;
            }
            computeMatrixProperties = Enabled;
            markForMatrixPropertyRecomputation(root);
            recomputeMatrixProperties(root);

            if (optimalPos > originalPos) {
                auto t = invVarMap[originalPos];
                for (int j = originalPos; j < optimalPos; ++j) { invVarMap[j] = invVarMap[j+1]; varMap[invVarMap[j]] = j; }
                invVarMap[optimalPos] = t; varMap[t] = optimalPos;
            } else if (optimalPos < originalPos) {
                auto t = invVarMap[originalPos];
                for (int j = originalPos; j > optimalPos; --j) { invVarMap[j] = invVarMap[j-1]; varMap[invVarMap[j]] = j; }
                invVarMap[optimalPos] = t; varMap[t] = optimalPos;
            }
        }
        return {in, total_min, total_max};
    }

    // ================== igLbSifting (LB pruning + IG guidance) ==================
    // Relaxed LB: use a factor to avoid premature pruning on complex circuits.
    // The factor trades off pruning aggressiveness vs. exploration depth.
    static constexpr double LB_RELAX_FACTOR = 1.1;  // was 1.0 (original: lb > min)

    std::tuple<Edge, unsigned int, unsigned int> Package::igLbSifting(
        Edge in, std::map<unsigned short, unsigned short>& varMap,
        const InteractionGraph& ig)
    {
        const auto n = static_cast<short>(in.p->v);

        std::vector<bool> free(n, true);
        std::map<unsigned short, unsigned short> invVarMap{};
        for (const auto& i : varMap) invVarMap[i.second] = i.first;

        computeMatrixProperties = Disabled;
        Edge root{in};
        unsigned int total_max = size(in);
        unsigned int total_min = total_max;

        short pos = -1;
        for (int i = 0; i < n; ++i) {
            assert(is_globally_consistent_dd(in));
            unsigned long min = size(in);

            // IG variable selection
            pos = selectNextVariable(n, free, varMap, active.data(), ig);
            if (pos < 0) break;

            free.at(varMap[pos]) = false;
            short optimalPos = pos;
            short originalPos = pos;

            // IG direction
            bool upFirst = shouldSiftUpFirst(pos, pos, n, invVarMap, ig);

            if (!upFirst) {
                // sifting down with relaxed LB pruning
                while (pos > 0) {
                    uint64_t lb = computeLowerBoundDown(varMap, pos);
                    if (lb > min * LB_RELAX_FACTOR) break;
                    exchangeBaseCase(pos, in, varMap);
                    auto s = size(in); total_min = std::min(total_min, s); total_max = std::max(total_max, s);
                    assert(is_locally_consistent_dd(in)); --pos;
                    if (s < min) { min = s; optimalPos = pos; }
                }
                // sifting up with relaxed LB pruning
                while (pos < n - 1) {
                    uint64_t lb = computeLowerBoundUp(varMap, pos);
                    if (lb > min * LB_RELAX_FACTOR) break;
                    exchangeBaseCase(pos + 1, in, varMap);
                    auto s = size(in); total_min = std::min(total_min, s); total_max = std::max(total_max, s);
                    assert(is_locally_consistent_dd(in)); ++pos;
                    if (s < min) { min = s; optimalPos = pos; }
                }
                // back to optimal
                while (pos > optimalPos) {
                    exchangeBaseCase(pos, in, varMap);
                    auto s = size(in); total_min = std::min(total_min, s); total_max = std::max(total_max, s);
                    assert(is_locally_consistent_dd(in)); --pos;
                }
                while (pos < optimalPos) {
                    exchangeBaseCase(pos + 1, in, varMap);
                    auto s = size(in); total_min = std::min(total_min, s); total_max = std::max(total_max, s);
                    assert(is_locally_consistent_dd(in)); ++pos;
                }
            } else {
                // sifting up with relaxed LB pruning
                while (pos < n - 1) {
                    uint64_t lb = computeLowerBoundUp(varMap, pos);
                    if (lb > min * LB_RELAX_FACTOR) break;
                    exchangeBaseCase(pos + 1, in, varMap);
                    auto s = size(in); total_min = std::min(total_min, s); total_max = std::max(total_max, s);
                    assert(is_locally_consistent_dd(in)); ++pos;
                    if (s < min) { min = s; optimalPos = pos; }
                }
                // sifting down with relaxed LB pruning
                while (pos > 0) {
                    uint64_t lb = computeLowerBoundDown(varMap, pos);
                    if (lb > min * LB_RELAX_FACTOR) break;
                    exchangeBaseCase(pos, in, varMap);
                    assert(is_locally_consistent_dd(in));
                    auto s = size(in); total_min = std::min(total_min, s); total_max = std::max(total_max, s);
                    --pos;
                    if (s < min) { min = s; optimalPos = pos; }
                }
                // back to optimal
                while (pos < optimalPos) {
                    exchangeBaseCase(pos + 1, in, varMap);
                    auto s = size(in); total_min = std::min(total_min, s); total_max = std::max(total_max, s);
                    assert(is_locally_consistent_dd(in)); ++pos;
                }
                while (pos > optimalPos) {
                    exchangeBaseCase(pos, in, varMap);
                    auto s = size(in); total_min = std::min(total_min, s); total_max = std::max(total_max, s);
                    assert(is_locally_consistent_dd(in)); --pos;
                }
            }

            initComputeTable();
            if (unnormalizedNodes > 0) {
                auto oldroot = root; root = renormalize(root);
                decRef(oldroot); incRef(root); in.p = root.p; in.w = root.w;
            }
            computeMatrixProperties = Enabled;
            markForMatrixPropertyRecomputation(root);
            recomputeMatrixProperties(root);

            // Rebuild invVarMap from varMap (3-arg exchangeBaseCase already updated varMap)
            invVarMap.clear();
            for (const auto& kv : varMap) invVarMap[kv.second] = kv.first;
        }
        return {in, total_min, total_max};
    }

} // namespace dd
