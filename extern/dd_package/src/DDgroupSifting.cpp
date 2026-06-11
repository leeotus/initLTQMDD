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
#include <cstdlib>

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

                constexpr int MAX_GROUP_PLACE = 8;
                int placedCount = 0;

                for (short member : members) {
                    if (placedCount >= MAX_GROUP_PLACE) break;
                    if (member >= n || placed[member]) continue;

                    auto memberIt = varMap.find((unsigned short)member);
                    if (memberIt == varMap.end()) continue;
                    short memberPos = static_cast<short>(memberIt->second);

                    // Refresh representative's current DD position
                    short repPos = static_cast<short>(varMap[pos]);
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

                    // Commit only if strictly improves
                    if (bestTarget != memberPos && bestSize < size(in)) {
                        short tmpPos = memberPos;
                        while (tmpPos < bestTarget) { exchangeBaseCase(tmpPos + 1, in, varMap); ++tmpPos; }
                        while (tmpPos > bestTarget) { exchangeBaseCase(tmpPos, in, varMap); --tmpPos; }

                        // Safety: bail out if DD blew up
                        unsigned long afterSize = size(in);
                        if (afterSize > sizeBeforeGroup * 2) break;
                    }

                    total_min = std::min(total_min, size(in));
                    total_max = std::max(total_max, size(in));
                    placed[member] = true;
                    free.at(varMap[member]) = false;
                    ++placedCount;
                }

                placed[pos] = true;
            }
        }

        return {in, total_min, total_max};
    }

    // ================== igGroupSifting ==================
    // Phase 1: Full igLbSifting for good base ordering
    // Phase 2: Conservative group adjustment for symmetric members

    std::tuple<Edge, unsigned int, unsigned int> Package::igGroupSifting(
        Edge in, std::map<unsigned short, unsigned short>& varMap,
        InteractionGraph& ig)
    {
        if (ig.n <= 0) return sifting(in, varMap);

        // Phase 1: Use igLbSifting to get a solid base ordering
        auto [result, tmin, tmax] = igLbSifting(in, varMap, ig);
        in = result;
        unsigned int total_min = tmin;
        unsigned int total_max = tmax;

        // Phase 2: Detect symmetry and try conservative group adjustments
        ig.detectSymmetry();
        if (ig.symmetricGroups.empty()) {
            return {in, total_min, total_max};
        }

        const auto n = static_cast<short>(in.p->v);
        constexpr short MAX_MOVE_DIST = 2;

        for (auto& group : ig.symmetricGroups) {
            if (group.size() < 2) continue;

            // Find the group representative (first member in the group)
            short rep = group[0];
            auto repIt = varMap.find((unsigned short)rep);
            if (repIt == varMap.end() || rep >= n) continue;

            for (size_t gi = 1; gi < group.size() && gi <= 4; ++gi) {
                short member = group[gi];
                if (member >= n) continue;
                auto memberIt = varMap.find((unsigned short)member);
                if (memberIt == varMap.end()) continue;

                short memberPos = static_cast<short>(memberIt->second);
                short repPos = static_cast<short>(varMap[(unsigned short)rep]);

                short dist = (memberPos > repPos) ? (memberPos - repPos) : (repPos - memberPos);
                if (dist <= 1 || dist > MAX_MOVE_DIST + 1) continue;

                unsigned long curSize = size(in);
                short bestTarget = memberPos;
                unsigned long bestSize = curSize;

                // Try repPos+1
                if (repPos + 1 <= n - 1 && std::abs(memberPos - (repPos + 1)) <= MAX_MOVE_DIST) {
                    short target = repPos + 1;
                    short tmpPos = memberPos;
                    while (tmpPos < target) { exchangeBaseCase(tmpPos + 1, in, varMap); ++tmpPos; }
                    while (tmpPos > target) { exchangeBaseCase(tmpPos, in, varMap); --tmpPos; }
                    auto s = size(in);
                    if (s < bestSize) { bestSize = s; bestTarget = target; }
                    while (tmpPos < memberPos) { exchangeBaseCase(tmpPos + 1, in, varMap); ++tmpPos; }
                    while (tmpPos > memberPos) { exchangeBaseCase(tmpPos, in, varMap); --tmpPos; }
                }

                // Try repPos-1
                if (repPos - 1 >= 0 && std::abs(memberPos - (repPos - 1)) <= MAX_MOVE_DIST) {
                    short target = repPos - 1;
                    short tmpPos = memberPos;
                    while (tmpPos < target) { exchangeBaseCase(tmpPos + 1, in, varMap); ++tmpPos; }
                    while (tmpPos > target) { exchangeBaseCase(tmpPos, in, varMap); --tmpPos; }
                    auto s = size(in);
                    if (s < bestSize) { bestSize = s; bestTarget = target; }
                    while (tmpPos < memberPos) { exchangeBaseCase(tmpPos + 1, in, varMap); ++tmpPos; }
                    while (tmpPos > memberPos) { exchangeBaseCase(tmpPos, in, varMap); --tmpPos; }
                }

                // Commit only if strictly improves
                if (bestTarget != memberPos && bestSize < curSize) {
                    if (memberPos < bestTarget) {
                        while (memberPos < bestTarget) { exchangeBaseCase(memberPos + 1, in, varMap); ++memberPos; }
                    } else {
                        while (memberPos > bestTarget) { exchangeBaseCase(memberPos, in, varMap); --memberPos; }
                    }
                    unsigned long afterSize = size(in);
                    total_min = std::min(total_min, (unsigned int)afterSize);
                    if (afterSize > curSize) break;
                }
            }
        }

        total_min = std::min(total_min, size(in));
        total_max = std::max(total_max, size(in));
        return {in, total_min, total_max};
    }

} // namespace dd
