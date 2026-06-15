/*
 * Group Sifting for QMDD variable reordering.
 *
 * Based on the classic BDD group sifting algorithm (Panda & Somenzi, 1995):
 *   1. Detect symmetric groups via interaction graph profiles
 *   2. Gather group members to be adjacent (contiguous block)
 *   3. Sift the entire group as a single unit through all positions
 *   4. At the best position, optimize group-internal order via adjacent swaps
 *
 * Two variants:
 *   groupSifting:   symmetry detection + group sifting (no IG direction)
 *   igGroupSifting: symmetry detection + IG-guided direction + LB pruning
 */

#include "DDpackage.h"
#include <algorithm>
#include <numeric>
#include <cassert>
#include <cstdlib>
#include <vector>
#include <map>

namespace dd {

    // Move entire group one position down (toward level 0).
    // Group occupies [lo..hi]. After call: [lo-1..hi-1].
    static void moveGroupDown(Package& pkg, Edge& in,
                              std::map<unsigned short, unsigned short>& varMap,
                              short lo, short hi)
    {
        for (short p = lo; p <= hi; ++p) {
            pkg.exchangeBaseCase(static_cast<unsigned short>(p), in, varMap);
        }
    }

    // Move entire group one position up (toward higher levels).
    // Group occupies [lo..hi]. After call: [lo+1..hi+1].
    static void moveGroupUp(Package& pkg, Edge& in,
                            std::map<unsigned short, unsigned short>& varMap,
                            short lo, short hi)
    {
        for (short p = hi; p >= lo; --p) {
            pkg.exchangeBaseCase(static_cast<unsigned short>(p + 1), in, varMap);
        }
    }

    // Get sorted positions (ascending) of group members in current ordering.
    static std::vector<short> getGroupPositions(
        const std::vector<short>& groupVars,
        const std::map<unsigned short, unsigned short>& varMap)
    {
        std::vector<short> positions;
        positions.reserve(groupVars.size());
        for (short var : groupVars) {
            auto it = varMap.find(static_cast<unsigned short>(var));
            if (it != varMap.end()) {
                positions.push_back(static_cast<short>(it->second));
            }
        }
        std::sort(positions.begin(), positions.end());
        return positions;
    }

    // Gather group members into a contiguous block.
    static void gatherGroup(Package& pkg, Edge& in,
                            std::map<unsigned short, unsigned short>& varMap,
                            const std::vector<short>& groupVars)
    {
        if (groupVars.size() <= 1) return;

        std::vector<std::pair<short, short>> varPos;
        for (short var : groupVars) {
            auto it = varMap.find(static_cast<unsigned short>(var));
            if (it != varMap.end()) {
                varPos.push_back({static_cast<short>(it->second), var});
            }
        }
        if (varPos.size() <= 1) return;

        std::sort(varPos.begin(), varPos.end());

        for (size_t i = 1; i < varPos.size(); ++i) {
            short targetPos = varPos[i - 1].first + 1;
            short curVar = varPos[i].second;
            short curPos = static_cast<short>(varMap[static_cast<unsigned short>(curVar)]);

            while (curPos > targetPos) {
                pkg.exchangeBaseCase(static_cast<unsigned short>(curPos), in, varMap);
                --curPos;
            }
            while (curPos < targetPos) {
                pkg.exchangeBaseCase(static_cast<unsigned short>(curPos + 1), in, varMap);
                ++curPos;
            }
            varPos[i].first = curPos;
        }
    }

    // Optimize internal order of a contiguous group via bubble-sort style swaps.
    static void optimizeGroupInternalOrder(Package& pkg, Edge& in,
                                           std::map<unsigned short, unsigned short>& varMap,
                                           const std::vector<short>& groupVars)
    {
        if (groupVars.size() <= 1) return;

        auto positions = getGroupPositions(groupVars, varMap);
        if (positions.size() <= 1) return;

        short lo = positions.front();
        short hi = positions.back();
        if (hi - lo + 1 != static_cast<short>(positions.size())) return;

        bool improved = true;
        int maxPasses = static_cast<int>(groupVars.size());

        while (improved && maxPasses-- > 0) {
            improved = false;
            for (short p = lo; p < hi; ++p) {
                unsigned int beforeSwap = pkg.size(in);
                pkg.exchangeBaseCase(static_cast<unsigned short>(p + 1), in, varMap);
                unsigned int afterSwap = pkg.size(in);

                if (afterSwap < beforeSwap) {
                    improved = true;
                } else {
                    pkg.exchangeBaseCase(static_cast<unsigned short>(p + 1), in, varMap);
                }
            }
        }
    }

    // ================== groupSifting ==================

    std::tuple<Edge, unsigned int, unsigned int> Package::groupSifting(
        Edge in, std::map<unsigned short, unsigned short>& varMap,
        InteractionGraph& ig)
    {
        if (ig.n <= 0) return sifting(in, varMap);

        const auto n = static_cast<short>(in.p->v);
        if (n <= 1) return {in, size(in), size(in)};

        ig.detectSymmetry();

        computeMatrixProperties = Disabled;
        Edge root{in};
        unsigned int total_max = size(in);
        unsigned int total_min = total_max;

        // Build sift units: groups (size>=2) + singletons
        std::vector<std::vector<short>> siftUnits;
        std::vector<bool> inGroup(n, false);

        for (auto& group : ig.symmetricGroups) {
            std::vector<short> valid;
            for (short var : group) {
                if (var < n && varMap.count(static_cast<unsigned short>(var))) {
                    valid.push_back(var);
                    inGroup[var] = true;
                }
            }
            if (valid.size() >= 2) {
                siftUnits.push_back(std::move(valid));
            } else {
                for (short v : valid) inGroup[v] = false;
            }
        }
        for (short var = 0; var < n; ++var) {
            if (!inGroup[var] && varMap.count(static_cast<unsigned short>(var))) {
                siftUnits.push_back({var});
            }
        }

        // Process order: largest active node count first
        std::sort(siftUnits.begin(), siftUnits.end(),
            [this, &varMap](const std::vector<short>& a, const std::vector<short>& b) {
                int maxA = 0, maxB = 0;
                for (short v : a) {
                    auto it = varMap.find(static_cast<unsigned short>(v));
                    if (it != varMap.end()) maxA = std::max(maxA, active.at(it->second));
                }
                for (short v : b) {
                    auto it = varMap.find(static_cast<unsigned short>(v));
                    if (it != varMap.end()) maxB = std::max(maxB, active.at(it->second));
                }
                return maxA > maxB;
            });

        for (auto& unit : siftUnits) {
            if (unit.size() == 1) {
                // Standard singleton sift
                short var = unit[0];
                short pos = static_cast<short>(varMap[static_cast<unsigned short>(var)]);
                short optimalPos = pos;
                unsigned long min = size(in);

                if (pos < n / 2) {
                    while (pos > 0) {
                        exchangeBaseCase(pos, in, varMap);
                        auto s = size(in);
                        total_min = std::min(total_min, s); total_max = std::max(total_max, s);
                        --pos;
                        if (s < min) { min = s; optimalPos = pos; }
                    }
                    while (pos < n - 1) {
                        exchangeBaseCase(pos + 1, in, varMap);
                        auto s = size(in);
                        total_min = std::min(total_min, s); total_max = std::max(total_max, s);
                        ++pos;
                        if (s < min) { min = s; optimalPos = pos; }
                    }
                } else {
                    while (pos < n - 1) {
                        exchangeBaseCase(pos + 1, in, varMap);
                        auto s = size(in);
                        total_min = std::min(total_min, s); total_max = std::max(total_max, s);
                        ++pos;
                        if (s < min) { min = s; optimalPos = pos; }
                    }
                    while (pos > 0) {
                        exchangeBaseCase(pos, in, varMap);
                        auto s = size(in);
                        total_min = std::min(total_min, s); total_max = std::max(total_max, s);
                        --pos;
                        if (s < min) { min = s; optimalPos = pos; }
                    }
                }
                while (pos < optimalPos) { exchangeBaseCase(pos + 1, in, varMap); ++pos; }
                while (pos > optimalPos) { exchangeBaseCase(pos, in, varMap); --pos; }
            } else {
                // Group sift: gather then sift as block
                gatherGroup(*this, in, varMap, unit);

                auto positions = getGroupPositions(unit, varMap);
                if (positions.empty()) continue;

                short lo = positions.front();
                short hi = positions.back();
                short optimalLo = lo;
                unsigned long min = size(in);

                short groupSize = static_cast<short>(positions.size());
                if (lo < (n - groupSize) / 2) {
                    // Lower half: sift down first
                    while (lo > 0) {
                        moveGroupDown(*this, in, varMap, lo, hi);
                        --lo; --hi;
                        auto s = size(in);
                        total_min = std::min(total_min, s); total_max = std::max(total_max, s);
                        if (s < min) { min = s; optimalLo = lo; }
                    }
                    while (hi < n - 1) {
                        moveGroupUp(*this, in, varMap, lo, hi);
                        ++lo; ++hi;
                        auto s = size(in);
                        total_min = std::min(total_min, s); total_max = std::max(total_max, s);
                        if (s < min) { min = s; optimalLo = lo; }
                    }
                } else {
                    // Upper half: sift up first
                    while (hi < n - 1) {
                        moveGroupUp(*this, in, varMap, lo, hi);
                        ++lo; ++hi;
                        auto s = size(in);
                        total_min = std::min(total_min, s); total_max = std::max(total_max, s);
                        if (s < min) { min = s; optimalLo = lo; }
                    }
                    while (lo > 0) {
                        moveGroupDown(*this, in, varMap, lo, hi);
                        --lo; --hi;
                        auto s = size(in);
                        total_min = std::min(total_min, s); total_max = std::max(total_max, s);
                        if (s < min) { min = s; optimalLo = lo; }
                    }
                }

                while (lo < optimalLo) { moveGroupUp(*this, in, varMap, lo, hi); ++lo; ++hi; }
                while (lo > optimalLo) { moveGroupDown(*this, in, varMap, lo, hi); --lo; --hi; }

                optimizeGroupInternalOrder(*this, in, varMap, unit);
            }

            initComputeTable();
            if (unnormalizedNodes > 0) {
                auto oldroot = root; root = renormalize(root);
                decRef(oldroot); incRef(root); in.p = root.p; in.w = root.w;
            }
            computeMatrixProperties = Enabled;
            markForMatrixPropertyRecomputation(root);
            recomputeMatrixProperties(root);
        }

        total_min = std::min(total_min, size(in));
        total_max = std::max(total_max, size(in));
        return {in, total_min, total_max};
    }

    // ================== igGroupSifting ==================

    std::tuple<Edge, unsigned int, unsigned int> Package::igGroupSifting(
        Edge in, std::map<unsigned short, unsigned short>& varMap,
        InteractionGraph& ig)
    {
        if (ig.n <= 0) return sifting(in, varMap);

        const auto n = static_cast<short>(in.p->v);
        if (n <= 1) return {in, size(in), size(in)};

        ig.detectSymmetry();

        std::map<unsigned short, unsigned short> invVarMap{};
        for (const auto& i : varMap) invVarMap[i.second] = i.first;

        computeMatrixProperties = Disabled;
        Edge root{in};
        unsigned int total_max = size(in);
        unsigned int total_min = total_max;

        // Build sift units
        std::vector<std::vector<short>> siftUnits;
        std::vector<bool> inGroup(n, false);

        for (auto& group : ig.symmetricGroups) {
            std::vector<short> valid;
            for (short var : group) {
                if (var < n && varMap.count(static_cast<unsigned short>(var))) {
                    valid.push_back(var);
                    inGroup[var] = true;
                }
            }
            if (valid.size() >= 2) {
                siftUnits.push_back(std::move(valid));
            } else {
                for (short v : valid) inGroup[v] = false;
            }
        }
        for (short var = 0; var < n; ++var) {
            if (!inGroup[var] && varMap.count(static_cast<unsigned short>(var))) {
                siftUnits.push_back({var});
            }
        }

        // Sort by IG degree (highest first)
        std::sort(siftUnits.begin(), siftUnits.end(),
            [&ig](const std::vector<short>& a, const std::vector<short>& b) {
                int maxA = 0, maxB = 0;
                for (short v : a) { if (v < ig.n) maxA = std::max(maxA, ig.degree[v]); }
                for (short v : b) { if (v < ig.n) maxB = std::max(maxB, ig.degree[v]); }
                return maxA > maxB;
            });

        for (auto& unit : siftUnits) {
            if (unit.size() == 1) {
                // Singleton with IG direction + LB pruning
                short var = unit[0];
                short pos = static_cast<short>(varMap[static_cast<unsigned short>(var)]);
                short optimalPos = pos;
                unsigned long min = size(in);

                bool upFirst = ig.shouldSiftUpFirst(var, pos, n, invVarMap);

                if (!upFirst) {
                    while (pos > 0) {
                        uint64_t lb = computeLowerBoundDown(varMap, pos);
                        if (lb > min) break;
                        exchangeBaseCase(pos, in, varMap);
                        auto s = size(in);
                        total_min = std::min(total_min, s); total_max = std::max(total_max, s);
                        --pos;
                        if (s < min) { min = s; optimalPos = pos; }
                    }
                    while (pos < n - 1) {
                        uint64_t lb = computeLowerBoundUp(varMap, pos);
                        if (lb > min) break;
                        exchangeBaseCase(pos + 1, in, varMap);
                        auto s = size(in);
                        total_min = std::min(total_min, s); total_max = std::max(total_max, s);
                        ++pos;
                        if (s < min) { min = s; optimalPos = pos; }
                    }
                } else {
                    while (pos < n - 1) {
                        uint64_t lb = computeLowerBoundUp(varMap, pos);
                        if (lb > min) break;
                        exchangeBaseCase(pos + 1, in, varMap);
                        auto s = size(in);
                        total_min = std::min(total_min, s); total_max = std::max(total_max, s);
                        ++pos;
                        if (s < min) { min = s; optimalPos = pos; }
                    }
                    while (pos > 0) {
                        uint64_t lb = computeLowerBoundDown(varMap, pos);
                        if (lb > min) break;
                        exchangeBaseCase(pos, in, varMap);
                        auto s = size(in);
                        total_min = std::min(total_min, s); total_max = std::max(total_max, s);
                        --pos;
                        if (s < min) { min = s; optimalPos = pos; }
                    }
                }
                while (pos < optimalPos) { exchangeBaseCase(pos + 1, in, varMap); ++pos; }
                while (pos > optimalPos) { exchangeBaseCase(pos, in, varMap); --pos; }
            } else {
                // Group sift with IG direction
                gatherGroup(*this, in, varMap, unit);

                auto positions = getGroupPositions(unit, varMap);
                if (positions.empty()) continue;

                short lo = positions.front();
                short hi = positions.back();
                short optimalLo = lo;
                unsigned long min = size(in);

                // IG direction from representative (highest degree member)
                short rep = unit[0];
                for (short v : unit) {
                    if (v < ig.n && ig.degree[v] > ig.degree[rep]) rep = v;
                }
                short repPos = static_cast<short>(varMap[static_cast<unsigned short>(rep)]);
                bool upFirst = ig.shouldSiftUpFirst(rep, repPos, n, invVarMap);

                if (!upFirst) {
                    while (lo > 0) {
                        moveGroupDown(*this, in, varMap, lo, hi);
                        --lo; --hi;
                        auto s = size(in);
                        total_min = std::min(total_min, s); total_max = std::max(total_max, s);
                        if (s < min) { min = s; optimalLo = lo; }
                    }
                    while (hi < n - 1) {
                        moveGroupUp(*this, in, varMap, lo, hi);
                        ++lo; ++hi;
                        auto s = size(in);
                        total_min = std::min(total_min, s); total_max = std::max(total_max, s);
                        if (s < min) { min = s; optimalLo = lo; }
                    }
                } else {
                    while (hi < n - 1) {
                        moveGroupUp(*this, in, varMap, lo, hi);
                        ++lo; ++hi;
                        auto s = size(in);
                        total_min = std::min(total_min, s); total_max = std::max(total_max, s);
                        if (s < min) { min = s; optimalLo = lo; }
                    }
                    while (lo > 0) {
                        moveGroupDown(*this, in, varMap, lo, hi);
                        --lo; --hi;
                        auto s = size(in);
                        total_min = std::min(total_min, s); total_max = std::max(total_max, s);
                        if (s < min) { min = s; optimalLo = lo; }
                    }
                }

                while (lo < optimalLo) { moveGroupUp(*this, in, varMap, lo, hi); ++lo; ++hi; }
                while (lo > optimalLo) { moveGroupDown(*this, in, varMap, lo, hi); --lo; --hi; }

                optimizeGroupInternalOrder(*this, in, varMap, unit);
            }

            initComputeTable();
            if (unnormalizedNodes > 0) {
                auto oldroot = root; root = renormalize(root);
                decRef(oldroot); incRef(root); in.p = root.p; in.w = root.w;
            }
            computeMatrixProperties = Enabled;
            markForMatrixPropertyRecomputation(root);
            recomputeMatrixProperties(root);

            // Rebuild invVarMap
            invVarMap.clear();
            for (const auto& kv : varMap) invVarMap[kv.second] = kv.first;
        }

        total_min = std::min(total_min, size(in));
        total_max = std::max(total_max, size(in));
        return {in, total_min, total_max};
    }

} // namespace dd
