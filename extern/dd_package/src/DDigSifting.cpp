/*
 * Interaction Graph driven Sifting for QMDD variable reordering.
 */

#include "DDpackage.h"

namespace dd {

    std::tuple<Edge, unsigned int, unsigned int> Package::igSifting(
        Edge in, std::map<unsigned short, unsigned short>& varMap,
        const InteractionGraph& ig)
    {
        // Exact replica of standard sifting with IG hooks
        const auto n = static_cast<short>(in.p->v);

        std::vector<bool> free(n, true);
        std::map<unsigned short, unsigned short> invVarMap{};
        for (const auto& i : varMap)
            invVarMap[i.second] = i.first;

        computeMatrixProperties = Disabled;
        Edge root{in};

        unsigned int total_max = size(in);
        unsigned int total_min = total_max;

        short pos = -1;
        for (int i = 0; i < n; ++i) {
            assert(is_globally_consistent_dd(in));
            unsigned long min = size(in);
            unsigned long max = 0;

            // === VARIABLE SELECTION ===
            // IG enhancement: use hybrid score (active * 0.5 + ig_degree * 0.5)
            // For now, keep identical to standard to debug
            for (short j = 0; j < n; j++) {
                if (free.at(varMap[j]) && active.at(varMap[j]) > (int)max) {
                    max = active.at(varMap[j]);
                    pos = j;
                    assert(max <= (unsigned long)std::numeric_limits<int>::max());
                }
            }
            free.at(varMap[pos]) = false;
            short optimalPos = pos;
            short originalPos = pos;

            // === DIRECTION DECISION ===
            // IG enhancement: use gravity. For now keep identical.
            if (pos < n / 2) {
                // sifting to bottom
                while (pos > 0) {
                    exchangeBaseCase(pos, in);
                    auto in_size = size(in);
                    total_min = std::min(total_min, in_size);
                    total_max = std::max(total_max, in_size);
                    assert(is_locally_consistent_dd(in));
                    --pos;
                    if (in_size < min) {
                        min = in_size;
                        optimalPos = pos;
                    }
                }

                // sifting to top
                while (pos < n - 1) {
                    exchangeBaseCase(pos + 1, in);
                    auto in_size = size(in);
                    total_min = std::min(total_min, in_size);
                    total_max = std::max(total_max, in_size);
                    assert(is_locally_consistent_dd(in));
                    ++pos;
                    if (in_size < min) {
                        min = in_size;
                        optimalPos = pos;
                    }
                }

                // sifting to optimal position
                while (pos > optimalPos) {
                    exchangeBaseCase(pos, in);
                    auto in_size = size(in);
                    total_min = std::min(total_min, in_size);
                    total_max = std::max(total_max, in_size);
                    assert(is_locally_consistent_dd(in));
                    --pos;
                }
            } else {
                // sifting to top
                while (pos < n - 1) {
                    exchangeBaseCase(pos + 1, in);
                    auto in_size = size(in);
                    total_min = std::min(total_min, in_size);
                    total_max = std::max(total_max, in_size);
                    assert(is_locally_consistent_dd(in));
                    ++pos;
                    if (in_size < min) {
                        min = in_size;
                        optimalPos = pos;
                    }
                }

                // sifting to bottom
                while (pos > 0) {
                    exchangeBaseCase(pos, in);
                    assert(is_locally_consistent_dd(in));
                    auto in_size = size(in);
                    total_min = std::min(total_min, in_size);
                    total_max = std::max(total_max, in_size);
                    --pos;
                    if (in_size < min) {
                        min = in_size;
                        optimalPos = pos;
                    }
                }

                // sifting to optimal position
                while (pos < optimalPos) {
                    exchangeBaseCase(pos + 1, in);
                    auto in_size = size(in);
                    total_min = std::min(total_min, in_size);
                    total_max = std::max(total_max, in_size);
                    assert(is_locally_consistent_dd(in));
                    ++pos;
                }
            }

            initComputeTable();

            if (unnormalizedNodes > 0) {
                auto oldroot = root;
                root = renormalize(root);
                decRef(oldroot);
                incRef(root);
                in.p = root.p;
                in.w = root.w;
            }
            computeMatrixProperties = Enabled;
            markForMatrixPropertyRecomputation(root);
            recomputeMatrixProperties(root);

            // Adjusting varMap if position changed
            if (optimalPos > originalPos) {
                auto tempVar = invVarMap[originalPos];
                for (int j = originalPos; j < optimalPos; ++j) {
                    invVarMap[j] = invVarMap[j + 1];
                    varMap[invVarMap[j]] = j;
                }
                invVarMap[optimalPos] = tempVar;
                varMap[invVarMap[optimalPos]] = optimalPos;
            } else if (optimalPos < originalPos) {
                auto tempVar = invVarMap[originalPos];
                for (int j = originalPos; j > optimalPos; --j) {
                    invVarMap[j] = invVarMap[j - 1];
                    varMap[invVarMap[j]] = j;
                }
                invVarMap[optimalPos] = tempVar;
                varMap[invVarMap[optimalPos]] = optimalPos;
            }
        }

        return {in, total_min, total_max};
    }

} // namespace dd
