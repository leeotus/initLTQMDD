/*
 * This file is part of IIC-JKU DD package which is released under the MIT license.
 * See file README.md or go to http://iic.jku.at/eda/research/quantum_dd/ for more information.
 */

#include "DDpackage.h"
#include <iomanip>
#include <random>

namespace dd {

	void Package::recomputeMatrixProperties(Edge in) {
		if (isTerminal(in))
			return;

		if (in.p->computeMatrixProperties == Recompute) {
			for (const auto & e : in.p->e)
				recomputeMatrixProperties(e);

			in.p->computeMatrixProperties = computeMatrixProperties;
			checkSpecialMatrices(in.p);
		}
	}

	void Package::markForMatrixPropertyRecomputation(Edge in) {
		if (isTerminal(in))
			return;

		if (in.p->computeMatrixProperties != Recompute) {
			for (const auto & e : in.p->e)
				markForMatrixPropertyRecomputation(e);

			in.p->computeMatrixProperties = Recompute;
		}
	}

	void Package::resetNormalizationFactor(Edge in, Complex defaultValue) {
		if (isTerminal(in))
			return;

		if (in.p->normalizationFactor == defaultValue)
			return;

		for (const auto & e : in.p->e)
			resetNormalizationFactor(e, defaultValue);

		if (defaultValue == CN::ZERO && in.p->normalizationFactor != CN::ONE) {
		    assert(!CN::equalsOne(in.p->normalizationFactor));
            unnormalizedNodes--;
		}


		in.p->normalizationFactor = defaultValue;
	}

	Edge Package::renormalize(Edge in) {
        assert(is_locally_consistent_dd(in));

		const auto before = cn.cacheCount;
		in = renormalize2(in);
		resetNormalizationFactor(in, CN::ZERO);
		resetNormalizationFactor(in, CN::ONE);
		const auto after = cn.cacheCount;

		assert(before == after);
		return in;
	}

	Edge Package::renormalize2(Edge in) {
		if (isTerminal(in))
			return in;

		nOps[CTkind::renormalize]++;

		Complex weight = in.w;
		in.w = CN::ONE;

		Edge r = CTlookup(in, in, CTkind::renormalize);

		if (r.p!= nullptr) {
			if (r.w != CN::ZERO) {
				auto c = cn.getTempCachedComplex();
				CN::mul(c, r.w, weight);
				r.w = cn.lookup(c);
				if (CN::equalsZero(r.w)) {
					return DDzero;
				}
			}
			return r;
		}

		std::array<Edge, NEDGE> e{};
		for (int i=0; i<NEDGE; ++i) {
			e[i] = renormalize2(in.p->e[i]);
		}

		if (in.p->normalizationFactor != CN::ONE) {
		    assert(!CN::equalsOne(in.p->normalizationFactor));
			const Complex factor = in.p->normalizationFactor;
			in.p->normalizationFactor = CN::ONE;
			r = makeNonterminal(in.p->v, e);
			in.p->normalizationFactor = factor;
			auto c = cn.getTempCachedComplex();
			CN::mul(c, r.w, factor);
			r.w = cn.lookup(c);
		} else {
			r = makeNonterminal(in.p->v, e);
		}

		CTinsert(in, in, r, CTkind::renormalize);

		auto c = cn.getTempCachedComplex();
		CN::mul(c, r.w, weight);
		r.w = cn.lookup(c);

		return r;
	}

	void Package::substitute_node(short index, Edge original, Edge substitute, Edge in) {
        /*std::clog << "BEGIN: substitute_node(" << index << ", "
                  << debugnode_line(original.p) << ", "
                  << debugnode_line(substitute.p)
                  << ")\n";*/
        node_substitutions++;
        substitute_node_ut(index, original, substitute, in);
        substitute_node_dd(index, in, original, substitute, in);

        /*std::clog << "END: substitute_node(" << index << ", "
                  << debugnode_line(original.p) << ", "
                  << debugnode_line(substitute.p)
                  << ")\n";*/
	}

	void Package::substitute_node_ut(short index, Edge original, Edge substitute, Edge in) {
	    assert(index == original.p->v);
	    assert(index == substitute.p->v);
        /*std::clog << "begin: substitute_node_ut(" << index << ", "
                  << debugnode_line(original.p) << ", "
                  << debugnode_line(substitute.p)
                  << ")\n";*/

        assert(original.p != substitute.p);

        while(original.p->ref > 0) {
            decRef(original);
            incRef(substitute);
        }

        for (auto &bucket : Unique[index+1]) {
            NodePtr parent = bucket;
            while (parent != nullptr) {
                bool substituted = false;
                const std::size_t parent_original_hash = UThash(parent);
                for (auto& edge : parent->e){
                    if(parent->v - 1 != edge.p->v && !isTerminal(edge)) {
                        debugnode(parent);
                        throw std::runtime_error("p->v - 1 != p->e[_].p->v in substitute_node_ut");
                    }
                    if (edge.p == original.p) {
                        assert(edge.p->v == substitute.p->v);
                        assert(edge.p->v == original.p->v);
                        assert(parent->v == index + 1);
                        edge.p = substitute.p;
                        substituted = true;
                    }
                }

                if (substituted) {
                    assert(parent->v == index+1);
                    if(UThash(parent) != parent_original_hash) {
                        // at this point parent is still the same pointer but the hash has changed
                        Edge newEdge = UT_update_node({parent, CN::ZERO}, parent_original_hash, in);
                        //std::clog << "After UT_update_node: " << debugnode_line(parent) << "\n";

                        if(newEdge.p != parent) {
                            assert(newEdge.p->v == parent->v);
                            substitute_node(index+1, {parent, CN::ZERO}, newEdge, in);
                        }
                    }
                }

                parent = parent->next;
            }
        }

        /*std::clog << "end: substitute_node_ut(" << index << ", "
	        << debugnode_line(original.p) << ", "
	        << debugnode_line(substitute.p)
	        << ")\n";*/
	}

    void Package::substitute_node_dd(short index, Edge parent, const Edge original, const Edge substitute, Edge in) {
        assert(original.p->v == substitute.p->v);
        assert(parent.p != original.p);

	    if(isTerminal(parent) || parent.p->v <= original.p->v) {
            return;
        }

        for (auto& child : parent.p->e) {
            if(child.p == original.p) {
                auto ut_check_parent = UTcheck(parent);
                auto parent_pre_hash = UThash(parent.p);

                assert(ut_check_parent[0] != '!' && ut_check_parent != "0");
                assert(UTcheck(child) == "not_found");

                //std::clog << "DD_substitute_node(" << debugnode_line(parent.p) << ", " << debugnode_line(original.p) << ", " << debugnode_line(substitute.p) << ") found a node\n";
                //debugnode(parent.p);
                child.p = substitute.p;

                if(UThash(parent.p) != parent_pre_hash) {
                    Edge newEdge = UT_update_node(parent, parent_pre_hash, in);
                    //std::clog << "After UT_update_node: " << debugnode_line(parent.p) << "\n";

                    if(newEdge.p != parent.p) {
                        assert(newEdge.p->v == parent.p->v);
                        substitute_node(parent.p->v, parent, newEdge, in);
                    }
                }
            } else {
                substitute_node_dd(index, child, original, substitute, in);
            }
        }
	}

	void Package::reuseNonterminal(short v, const Edge *edge, NodePtr p, Edge in) {
		Edge new_e{p, CN::ONE}; //新建一条指向一层结点的边
        new_e.p->computeMatrixProperties = computeMatrixProperties;
		assert(v == p->v); //确认正常
        new_e.p->v = v;
		std::memcpy(new_e.p->e, edge, NEDGE * sizeof(Edge));

        assert(v-1 == edge[0].p->v || isTerminal(edge[0]));
        assert(v-1 == edge[1].p->v || isTerminal(edge[1]));
        assert(v-1 == edge[2].p->v || isTerminal(edge[2]));
        assert(v-1 == edge[3].p->v || isTerminal(edge[3]));

        assert(is_locally_consistent_dd(new_e));

        {
            auto olde = new_e;
            new_e = normalize(new_e, false); //对new_e规范化

            if (olde.p != new_e.p) {
                throw std::runtime_error("Normalized edge changed to different node.");
            }
        }

        assert(is_locally_consistent_dd(new_e));

		if (new_e.w != CN::ONE) {
			if (new_e.p->normalizationFactor == CN::ONE) {
				unnormalizedNodes++;
				assert(!CN::equalsOne(new_e.w));
                new_e.p->normalizationFactor = new_e.w;
			} else {
				auto c = cn.getTempCachedComplex();
				CN::mul(c, new_e.p->normalizationFactor, new_e.w);
                new_e.p->normalizationFactor = cn.lookup(c);
			}
            new_e.w = CN::ONE;

			if (new_e.p->normalizationFactor == CN::ONE)
				unnormalizedNodes--;
		}
		// !problematic if lookup would change NodePtr
        assert(is_locally_consistent_dd(new_e));

        Edge lookedup_e = UTlookup(new_e, true);

        assert(v-1 == new_e.p->e[0].p->v || isTerminal(new_e.p->e[0]));
        assert(v-1 == new_e.p->e[1].p->v || isTerminal(new_e.p->e[1]));
        assert(v-1 == new_e.p->e[2].p->v || isTerminal(new_e.p->e[2]));
        assert(v-1 == new_e.p->e[3].p->v || isTerminal(new_e.p->e[3]));

        assert(v-1 == lookedup_e.p->e[0].p->v || isTerminal(lookedup_e.p->e[0]));
        assert(v-1 == lookedup_e.p->e[1].p->v || isTerminal(lookedup_e.p->e[1]));
        assert(v-1 == lookedup_e.p->e[2].p->v || isTerminal(lookedup_e.p->e[2]));
        assert(v-1 == lookedup_e.p->e[3].p->v || isTerminal(lookedup_e.p->e[3]));

        if(lookedup_e.p != new_e.p) {
            node_collapses++;
            substitute_node(v, new_e, lookedup_e, in);
            check_node_is_really_gone(new_e.p, in);
            assert(is_locally_consistent_dd(lookedup_e));
		} else {
            assert(is_locally_consistent_dd(new_e));
        }
	}

	void Package::exchange(unsigned short i, unsigned short j) {
	    //std::clog << "    exchange(" << i << ", " << j << ")\n";
		if (i == j) {
			return;
		} else if (i > j) {
			return exchange(j, i);
		}

		if ((i + 1) == j)
			return exchangeBaseCase(j, {});

		auto g = static_cast<short>(i + 1);

		// shuffeling the lower level i up until it is in its position
		while (g < j)
            exchangeBaseCase(g++, {});
        exchangeBaseCase(g, {});

		// shuffeling the upper level j down until it is in its position
		while (g > i+1)
            exchangeBaseCase(--g, {});
	}


    dd::Edge Package::exchange2(unsigned short i, unsigned short j, std::map<unsigned short, unsigned short>& varMap, Edge in) {
        //std::clog << "    exchange2(" << i << ", " << j << ")\n";

        std::map<unsigned short, unsigned short> invVarMap;
        for (const auto& entry : varMap)
            invVarMap[entry.second] = entry.first;

        if (i == j) {
            return in;
        } else if (i > j) {
            return exchange2(j, i, varMap, in);
        }

        auto g = static_cast<short>(i + 1);

        // shuffeling the lower level i up until it is in its position
        while (g < j)
            exchangeBaseCase(g++, in);
        exchangeBaseCase(g, in);

        // shuffeling the upper level j down until it is in its position
        while (g > i+1)
            exchangeBaseCase(--g, in);

        if (unnormalizedNodes > 0) {
            //std::clog << "{" << unnormalizedNodes << "} ";
            auto oldroot = in;
            in = renormalize(in);
            decRef(oldroot);
            incRef(in);
            assert(is_locally_consistent_dd(in));
            if (unnormalizedNodes > 0) {
                throw std::runtime_error("Renormalization failed. " + std::to_string(unnormalizedNodes) + " unnormalized nodes remaining.");
            }
        }

        auto tempVar = varMap[invVarMap[i]];
        varMap[invVarMap[i]] = varMap[invVarMap[j]];
        varMap[invVarMap[j]] = tempVar;
        assert(is_locally_consistent_dd(in));
        return in;
    }

	void Package::exchangeBaseCase(unsigned short i, Edge in) { //传入电路变量 和 dd的根指针
		// std::clog << "ex:" << i << ", ";
		
	    exchange_base_cases++;
		// copy unique table from higher variable and empty it
		std::array<NodePtr, NBUCKET> table{};
		for (unsigned short bucket=0; bucket < NBUCKET; ++bucket) {
			table.at(bucket) = Unique[i][bucket]; //把i变量的唯一表拿出来
			Unique[i][bucket] = nullptr; //把属于i变量的唯一表清空
		} //到此，已经把i变量的唯一表拿出来了，并且清空了它

		initComputeTable(); //初始化计算表

		// iterate over all obtained nodes
		for (unsigned short bucket=0; bucket < NBUCKET; ++bucket) { //对i变量的每个哈希桶进行遍历
			NodePtr p = table[bucket]; //取出当前桶里的冲突链
			while (p != nullptr) { //遍历冲突链的每个结点
				NodePtr pnext = p->next;
				assert(p->v == i); //确认当前指针指向的结点就是属于i变量的
                assert(i-1 == p->e[0].p->v || isTerminal(p->e[0])); //确认这个结点的四个孩子结点都是属于i-1变量，或是它的孩子是终端结点
                assert(i-1 == p->e[1].p->v || isTerminal(p->e[1]));
                assert(i-1 == p->e[2].p->v || isTerminal(p->e[2]));
                assert(i-1 == p->e[3].p->v || isTerminal(p->e[3]));
				if (p->ref != 0) { //确认这个结点引用计数不是0
                    exchangeBaseCase2(p, i, in); //传入该结点的指针，变量索引，DD指针
				}
                assert(p->v == i);
                assert(i-1 == p->e[0].p->v || isTerminal(p->e[0]));
                assert(i-1 == p->e[1].p->v || isTerminal(p->e[1]));
                assert(i-1 == p->e[2].p->v || isTerminal(p->e[2]));
                assert(i-1 == p->e[3].p->v || isTerminal(p->e[3]));
				p = pnext; //轮到冲突链的下一个结点
			}
		}
		if (node_substitutions > 0) {
            //std::clog << "♻" << node_substitutions;
        }
	}

	void Package::exchangeBaseCase2(NodePtr p, unsigned short index, Edge in) {
		Edge t[NEDGE][NEDGE]{ }; //创建一个矩阵
		assert(index > 0);
		assert(index == p->v); //确定传入结点是属于index变量

		// creating matrix T
		for (int i = 0; i < NEDGE; i++) { //遍历要处理结点的每条出边
			for (int j = 0; j < NEDGE; j++) { //遍历它的孩子结点的每一条出边
				if (p->e[i].p->v == index - 1) { //这里是图变量，如果要处理结点的孩子结点是属于它的下一个变量的（如本结点是i变量的，它的孩子结点是i-1变量的）
				    assert(!isTerminal(p->e[i]));

					t[j][i] = p->e[i].p->e[j]; //结点第i条边的孩子的第j条边放入矩阵中
					auto c = cn.getTempCachedComplex(); //拿到一个临时复数c
					CN::mul(c, p->e[i].p->e[j].w, p->e[i].w); //复数乘法，结点出边权值*它的孩子的出边权值=c
					if (p->e[i].p->normalizationFactor != CN::ONE) { //如果改结点的规范因子不是1，就把c乘上规范因子
						CN::mul(c, c, p->e[i].p->normalizationFactor);
					}
					t[j][i].w = cn.lookup(c); //给记录的边附上权值c
				} else if (isTerminal(p->e[i])) { //要处理结点的孩子结点是终端结点
					// edge pointing to a terminal
					t[j][i] = p->e[i]; //矩阵记录的就是结点的出边
                    assert(t[j][i].p->normalizationFactor == CN::ONE); //这样的话确定这条边的结点规范因子是1
				} else { //上面两种情况都不是的话，就出问题了，出错处理
				    debugnode(p);
				    std::stringstream hex_addr;
				    hex_addr << "0x" << std::hex << reinterpret_cast<std::uintptr_t>(p);
				    throw std::runtime_error("Edge " + std::to_string(i)
				        + " of " + hex_addr.str()
				        + " pointing to a skipped variable: "
				        + std::to_string(index) + " --> " + std::to_string(p->e[i].p->v));
				}
			}
		} //到这里，要处理的结点和它的孩子结点全都处理完，并把结点的孩子结点的出边和联合权值放入矩阵
        assert(is_locally_consistent_dd({p, CN::ZERO}));
		// creating new nodes and appending corresponding edges
		Edge newEdges[NEDGE]{ }; //创建要处理结点（P）的一组（4条）出边
		for (int x = 0; x < NEDGE; ++x) { //搞定p结点的四个孩子
			newEdges[x] = makeNonterminal(static_cast<short>(index - 1), t[x]); //新建结点，输入index-1和index-1（4个）结点里面的第x条边
            incRef(newEdges[x]); //给新结点加引用计数
			assert(is_locally_consistent_dd(newEdges[x]));
		}//到此，变量索引为index-1的四个结点完成规范化

		for (dd::Edge& x : p->e) //把结点原来只想它的孩子的边撤掉
            decRef(x);
		// reuse p to build new top node
		assert(p->ref > 0);
        reuseNonterminal(static_cast<short>(index), newEdges, p, in); //重用结点
		// p might be discarded at this point if nodes were substituted
	}

/// Dynamically reorder a given decision diagram with the current variable map using the specific strategy
/// \param in decision diagram to reorder
/// \param varMap stores the variable mapping. varMap[circuit qubit] = corresponding DD qubit, e.g.
///			given the varMap (reversed var. order):
/// 			0->2,
/// 			1->1,
/// 			2->0
/// 		the circuit operation "H q[0]" leads to the DD equivalent to "H q[varMap[0]]" = "H q[2]".
///			the qubits in the decision diagram are always ordered as n-1 > n-2 > ... > 1 > 0
/// \param strat strategy to apply
/// \return the resulting decision diagram (and the changed variable map and output permutation, which are returned as reference)
	Edge Package::dynamicReorder(Edge in, std::map<unsigned short, unsigned short>& varMap, DynamicReorderingStrategy strat) {
		switch (strat) {
			case None: return in;
			case Sifting: return std::get<0>(sifting(in, varMap));
			case Random:
                {
                    std::mt19937_64 mt;
                    return random(in, varMap, mt);
                }
			case Window3: return window3(in, varMap);
			case linearSifting: return linearAndSiftingAux(in, varMap, 1);
			case lbLinearSifting: return lbLinearAndSiftingAux(in, varMap);
			case mixlinearSifting: return mixLinearAndSiftingAux(in, varMap, 1);
		}

		return in;
	}

/// Apply sifting dynamic reordering to a decision diagram given the
/// current variable map \param in decision diagram to apply sifting to
/// \param varMap stores the variable mapping (cf. dynamicReorder(...))
/// \return the resulting decision diagram (and the changed variable map and output permutation, which are returned as reference)
    std::tuple<Edge, unsigned int, unsigned int> Package::sifting(Edge in, std::map<unsigned short, unsigned short>& varMap) {
		// const auto n = static_cast<short>(in.p->v + 1); //变量个数
		const auto n = static_cast<short>(in.p->v);

		std::vector<bool> free(n, true); //记录数组
		std::map<unsigned short, unsigned short> invVarMap{}; //DD qubit（变量） 到 电路 qubit（变量） 的映射
		for (const auto & i : varMap)
			invVarMap[i.second] = i.first; //DD qubit 到 电路 qubit 的映射

		computeMatrixProperties = Disabled;
		Edge root{in}; //声明root边

		//std::clog << "  Start Sifting. n=" << std::setw(2) << n << " -- ";
//		for (auto &entry: varMap) {
		    //std::clog << entry.second << " ";
//		}
		//std::clog << "\n";

		unsigned int total_max = size(in); //DD的大小
		unsigned int total_min = total_max;

        short pos = -1; //circuit variable (index)
        for (int i = 0; i < n; ++i) { //遍历各个变量
            assert(is_globally_consistent_dd(in));
            unsigned long min = size(in);
            unsigned long max = 0;

            // std::clog << "    " << i << "/" << n << " size=" << min << " | ";
            for (short j = 0; j < n; j++) {
                if (free.at(varMap[j]) && active.at(varMap[j]) > max) { //该变量没有被处理过并该变量存在结点
                    max = active.at(varMap[j]); //更新max
                    pos = j; //更新pos
                    assert(max <= std::numeric_limits<int>::max());
                }
            } //到此找到拥有最大结点数的变量，和该变量的索引（位置）
            free.at(varMap[pos]) = false; //设置选中的变量为处理状态
            short optimalPos = pos; 
            short originalPos = pos;

            if (pos < n / 2) {  // variable is in lower half -> sifting to bottom first
                // sifting to bottom
                while (pos > 0) {
                    exchangeBaseCase(pos, in);
                    auto in_size = size(in);
                    total_min = std::min(total_min, in_size);
                    total_max = std::max(total_max, in_size);

                    // std::clog << "↓" << in_size << " ";
                    assert(is_locally_consistent_dd(in));
                    --pos;
                    if (in_size < min) {
                        min = in_size;
                        optimalPos = pos;
                    }
                } //到这里被选中的变量走到了最下面

                // sifting to top
                while (pos < n - 1) {
                    exchangeBaseCase(pos + 1, in);
                    auto in_size = size(in);
                    total_min = std::min(total_min, in_size);
                    total_max = std::max(total_max, in_size);
                    // std::clog << "↑" << in_size << " ";
                    assert(is_locally_consistent_dd(in));
                    ++pos;
                    if (in_size < min) {
                        min = in_size;
                        optimalPos = pos;
                    }
                }

                //std::clog << "[" << min << "] ";

                // sifting to optimal position
                while (pos > optimalPos) {
                    exchangeBaseCase(pos, in);
                    auto in_size = size(in);
                    total_min = std::min(total_min, in_size);
                    total_max = std::max(total_max, in_size);
                    // std::clog << "↓" << in_size << " ";
                    assert(is_locally_consistent_dd(in));
                    --pos;
                }
            } else {  // variable is in upper half -> sifting to top first
                // sifting to top
                while (pos < n - 1) {
                    exchangeBaseCase(pos + 1, in);
                    auto in_size = size(in);
                    total_min = std::min(total_min, in_size);
                    total_max = std::max(total_max, in_size);
                    //std::clog << "↑" << in_size << " ";
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
                    //std::clog << "↓" << in_size << " ";
                    --pos;
                    if (in_size < min) {
                        min = in_size;
                        optimalPos = pos;
                    }
                }

                //std::clog << "[" << min << "] ";

                // sifting to optimal position
                while (pos < optimalPos) {
                    exchangeBaseCase(pos + 1, in);
                    auto in_size = size(in);
                    total_min = std::min(total_min, in_size);
                    total_max = std::max(total_max, in_size);
                    //std::clog << "↑" << size(in) << " ";
                    assert(is_locally_consistent_dd(in));
                    ++pos;
                }
            }

            initComputeTable();

                        // there are nodes which need to renormalized
            if (unnormalizedNodes > 0) {
                // std::clog << "{" << unnormalizedNodes << "} ";
                auto oldroot = root;
                root = renormalize(root);
                decRef(oldroot);
                incRef(root);
                in.p = root.p;
                in.w = root.w;
                if (unnormalizedNodes > 0) {
                    // throw std::runtime_error("Renormalization failed. " + std::to_string(unnormalizedNodes) + " unnormalized nodes remaining.");
                }
            }
            computeMatrixProperties = Enabled;
            markForMatrixPropertyRecomputation(root); //标记
            recomputeMatrixProperties(root);

            // Adjusting varMap if position changed
            if (optimalPos != originalPos) {
                //std::clog << "| " << originalPos << "-->" << optimalPos << " (min=" << min << "; real size=" << size(in) << ")\n";
            } else {
                //std::clog << "| ##### (min=" << min << "; real size=" << size(in) << ")\n";
            }

			//交换DD的层，即map的key
            if (optimalPos > originalPos) { //向上到最佳位置
                auto tempVar = invVarMap[originalPos]; //暂存最佳位置对应的电路变量
                for (int j = originalPos; j < optimalPos; ++j) {
                    invVarMap[j] = invVarMap[j + 1]; //调整 电路变量（qubit）
                    varMap[invVarMap[j]] = j; //更新DD变量（qubit）
                }
                invVarMap[optimalPos] = tempVar; //在合适位置放入最佳位置对应的电路变量
                varMap[invVarMap[optimalPos]] = optimalPos; //在合适位置放入最佳位置对应的DD变量
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
		return {in, total_min, total_max}; //返回DD指针和这次sifting过程中最大和最小的DD size
	}

// // 我的测试
//     std::tuple<Edge, unsigned int, unsigned int> Package::sifting(Edge in, std::map<unsigned short, unsigned short>& varMap) {
// 		const auto n = static_cast<short>(in.p->v + 1); //变量个数
// 		// const auto n = static_cast<short>(in.p->v);
// 		// std::clog<<varMap[in.p->v]<<active.at(varMap[in.p->v])<<std::endl;

// 		std::vector<bool> free(n, true); //记录数组
// 		std::map<unsigned short, unsigned short> invVarMap{}; //DD qubit（变量） 到 电路 qubit（变量） 的映射
// 		for (const auto & i : varMap)
// 			invVarMap[i.second] = i.first; //DD qubit 到 电路 qubit 的映射

// 		computeMatrixProperties = Disabled;
// 		Edge root{in}; //声明root边

// 		//std::clog << "  Start Sifting. n=" << std::setw(2) << n << " -- ";
// //		for (auto &entry: varMap) {
// 		    //std::clog << entry.second << " ";
// //		}
// 		//std::clog << "\n";

// 		unsigned int total_max = size(in); //DD的大小
// 		unsigned int total_min = total_max;

//         short pos = -1; //DD level （层）
//         for (int i = 0; i < n; ++i) { //遍历各个变量
//             assert(is_globally_consistent_dd(in));
//             unsigned long min = size(in);
//             unsigned long max = 0;

//             std::clog << "    " << i << "/" << n << " size=" << min << " | ";
//             for (short j = 0; j < n; j++) {
//                 if (free.at(varMap[j]) && active.at(varMap[j]) > max) { //该变量没有被处理过并该变量存在结点
//                     max = active.at(varMap[j]); //更新max
//                     pos = j; //更新pos
//                     assert(max <= std::numeric_limits<int>::max());
//                 }
//             } //到此找到拥有最大结点数的变量，和该变量的索引（位置）
//             free.at(varMap[pos]) = false; //设置选中的变量为处理状态
//             short optimalPos = pos; 
//             short originalPos = pos;

//             if(pos == 0) {
// 				// sifting to top
//                 while (pos < n - 1) {
//                     // exchangeBaseCase(pos + 1, in, varMap);
// 					linearInPlace(pos+1, in, varMap);
//                     auto in_size = size(in);
//                     total_min = std::min(total_min, in_size);
//                     total_max = std::max(total_max, in_size);
//                     std::clog << "↑" << in_size << " ";
//                     assert(is_locally_consistent_dd(in));
//                     ++pos;
//                     if (in_size < min) {
//                         min = in_size;
//                         optimalPos = pos;
//                     }
//                 }

//                 std::clog << "[" << min << "] ";

//                 // sifting to optimal position
//                 while (pos > optimalPos) {
                    
// 					linearInPlace(pos, in, varMap);
// 					// exchangeBaseCase(pos, in, varMap);
//                     auto in_size = size(in);
//                     total_min = std::min(total_min, in_size);
//                     total_max = std::max(total_max, in_size);
//                     std::clog << "↓" << in_size << " ";
//                     assert(is_locally_consistent_dd(in));
//                     --pos;
//                 }
// 			} 
// 			else if(pos == n-1) {
// 				// sifting to bottom
//                 while (pos > 0) {
//                     // exchangeBaseCase(pos, in, varMap);
// 					linearInPlace(pos, in, varMap);
//                     assert(is_locally_consistent_dd(in));
//                     auto in_size = size(in);
//                     total_min = std::min(total_min, in_size);
//                     total_max = std::max(total_max, in_size);
//                     std::clog << "↓" << in_size << " ";
//                     --pos;
//                     if (in_size < min) {
//                         min = in_size;
//                         optimalPos = pos;
//                     }
//                 }

//                 std::clog << "[" << min << "] ";

//                 // sifting to optimal position
//                 while (pos < optimalPos) {
// 					linearInPlace(pos+1, in, varMap);
//                     // exchangeBaseCase(pos + 1, in, varMap);
					
//                     auto in_size = size(in);
//                     total_min = std::min(total_min, in_size);
//                     total_max = std::max(total_max, in_size);
//                     std::clog << "↑" << size(in) << " ";
//                     assert(is_locally_consistent_dd(in));
//                     ++pos;
//                 }
// 			}
// 			else if (pos < n / 2) {  // variable is in lower half -> sifting to bottom first
//                 // sifting to bottom
//                 while (pos > 0) {
// 					linearInPlace(pos, in, varMap);
//                     // exchangeBaseCase(pos, in, varMap);
					
//                     auto in_size = size(in);
//                     total_min = std::min(total_min, in_size);
//                     total_max = std::max(total_max, in_size);

//                     std::clog << "↓" << in_size << " ";
//                     assert(is_locally_consistent_dd(in));
//                     --pos;
//                     if (in_size < min) {
//                         min = in_size;
//                         optimalPos = pos;
//                     }
//                 } //到这里被选中的变量走到了最下面

//                 // sifting to top
//                 while (pos < n - 1) {
					
//                     // exchangeBaseCase(pos + 1, in, varMap);
// 					linearInPlace(pos+1, in, varMap);
//                     auto in_size = size(in);
//                     total_min = std::min(total_min, in_size);
//                     total_max = std::max(total_max, in_size);
//                     std::clog << "↑" << in_size << " ";
//                     assert(is_locally_consistent_dd(in));
//                     ++pos;
//                     if (in_size < min) {
//                         min = in_size;
//                         optimalPos = pos;
//                     }
//                 }

//                 std::clog << "[" << min << "] ";

//                 // sifting to optimal position
//                 while (pos > optimalPos) {
                    
// 					linearInPlace(pos, in, varMap);
// 					// exchangeBaseCase(pos, in, varMap);
//                     auto in_size = size(in);
//                     total_min = std::min(total_min, in_size);
//                     total_max = std::max(total_max, in_size);
//                     std::clog << "↓" << in_size << " ";
//                     assert(is_locally_consistent_dd(in));
//                     --pos;
//                 }
//             } else {  // variable is in upper half -> sifting to top first
//                 // sifting to top
//                 while (pos < n - 1) {
// 					linearInPlace(pos+1, in, varMap);
//                     // exchangeBaseCase(pos + 1, in, varMap);
					
//                     auto in_size = size(in);
//                     total_min = std::min(total_min, in_size);
//                     total_max = std::max(total_max, in_size);
//                     std::clog << "↑" << in_size << " ";
//                     assert(is_locally_consistent_dd(in));
//                     ++pos;
//                     if (in_size < min) {
//                         min = in_size;
//                         optimalPos = pos;
//                     }
//                 }

//                 // sifting to bottom
//                 while (pos > 0) {
					
//                     // exchangeBaseCase(pos, in, varMap);
// 					linearInPlace(pos, in, varMap);
//                     assert(is_locally_consistent_dd(in));
//                     auto in_size = size(in);
//                     total_min = std::min(total_min, in_size);
//                     total_max = std::max(total_max, in_size);
//                     std::clog << "↓" << in_size << " ";
//                     --pos;
//                     if (in_size < min) {
//                         min = in_size;
//                         optimalPos = pos;
//                     }
//                 }

//                 std::clog << "[" << min << "] ";

//                 // sifting to optimal position
//                 while (pos < optimalPos) {
// 					linearInPlace(pos+1, in, varMap);
//                     // exchangeBaseCase(pos + 1, in, varMap);
					
//                     auto in_size = size(in);
//                     total_min = std::min(total_min, in_size);
//                     total_max = std::max(total_max, in_size);
//                     std::clog << "↑" << size(in) << " ";
//                     assert(is_locally_consistent_dd(in));
//                     ++pos;
//                 }
//             }

//             initComputeTable();

//                         // there are nodes which need to renormalized
//             if (unnormalizedNodes > 0) {
//                 std::clog << "{" << unnormalizedNodes << "} ";
//                 auto oldroot = root;
//                 root = renormalize(root);
//                 decRef(oldroot);
//                 incRef(root);
//                 in.p = root.p;
//                 in.w = root.w;
//                 if (unnormalizedNodes > 0) {
//                     throw std::runtime_error("Renormalization failed. " + std::to_string(unnormalizedNodes) + " unnormalized nodes remaining.");
//                 }
//             }
//             computeMatrixProperties = Enabled;
//             markForMatrixPropertyRecomputation(root); //标记
//             recomputeMatrixProperties(root);

//             // Adjusting varMap if position changed
//             if (optimalPos != originalPos) {
//                 //std::clog << "| " << originalPos << "-->" << optimalPos << " (min=" << min << "; real size=" << size(in) << ")\n";
//             } else {
//                 //std::clog << "| ##### (min=" << min << "; real size=" << size(in) << ")\n";
//             }

// 			//交换DD的层，即map的key
//             // if (optimalPos > originalPos) { //向上到最佳位置
//             //     auto tempVar = invVarMap[originalPos]; //暂存最佳位置对应的电路变量
//             //     for (int j = originalPos; j < optimalPos; ++j) {
//             //         invVarMap[j] = invVarMap[j + 1]; //调整 电路变量（qubit）
//             //         varMap[invVarMap[j]] = j; //更新DD变量（qubit）
//             //     }
//             //     invVarMap[optimalPos] = tempVar; //在合适位置放入最佳位置对应的电路变量
//             //     varMap[invVarMap[optimalPos]] = optimalPos; //在合适位置放入最佳位置对应的DD变量
//             // } else if (optimalPos < originalPos) {
//             //     auto tempVar = invVarMap[originalPos];
//             //     for (int j = originalPos; j > optimalPos; --j) {
//             //         invVarMap[j] = invVarMap[j - 1];
//             //         varMap[invVarMap[j]] = j;
//             //     }
//             //     invVarMap[optimalPos] = tempVar;
//             //     varMap[invVarMap[optimalPos]] = optimalPos;
//             // }



//         }
// 		return {in, total_min, total_max}; //返回DD指针和这次sifting过程中最大和最小的DD size
// 	}

	/// First counts the number of nodes in the given DD.
	/// Then a loop is executed nodeCount-many times and inside
	/// this loop two randomly selcted levels are swap.
    Edge Package::random(Edge in, std::map<unsigned short, unsigned short> &varMap, std::mt19937_64 &mt) {
		int n = (in.p->v + 1);
		unsigned long min = activeNodeCount;
		std::queue<Edge> q{ };
		int nodeCount = 0;
		static std::unordered_set<NodePtr> visited{NODECOUNT_BUCKETS}; // 2e6

		visited.clear();
		q.push(in);
		while (!q.empty()) {
			Edge e = q.front();
			if (visited.insert(e.p).second) ++nodeCount;
			q.pop();

			for (auto& x : e.p->e)
				if (x.p != nullptr && !isTerminal(x)) q.push(x);
		}
        std::uniform_int_distribution<int> dist(0, n-1);
		for (int x = 0; x < nodeCount; x++) {
			int i = dist(mt);
			int j = dist(mt);

			exchange(varMap[i], varMap[j]);

			if (min > activeNodeCount) {
				min = activeNodeCount;

				unsigned short temp = varMap[i];
				varMap[i] = varMap[j];
				varMap[j] = temp;
			} else {
				exchange(varMap[j], varMap[i]);
			}
		}

		return in;
	}

	Edge Package::window3(Edge in, std::map<unsigned short, unsigned short>& varMap) {
		std::map<unsigned short, unsigned short> invVarMap{ };
		int n = in.p->v;

		for (const auto& i : varMap) invVarMap[i.second] = i.first;

		for (int i = 0; i + 1 < n; i++) {
			int x = i;
			int y = x + 1;
			int z = y + 1;
			auto min = activeNodeCount;
			int best = 1;  // ABC

			exchange(x, y);  // BAC
			auto tempVar = varMap[invVarMap[x]];
			varMap[invVarMap[x]] = varMap[invVarMap[y]];
			varMap[invVarMap[y]] = tempVar;
			if (min > activeNodeCount) {
				best = 2;
				min = activeNodeCount;
			}

			exchange(y, z);  // BCA
			tempVar = varMap[invVarMap[z]];
			varMap[invVarMap[z]] = varMap[invVarMap[y]];
			varMap[invVarMap[y]] = tempVar;
			if (min > activeNodeCount) {
				best = 3;
				min = activeNodeCount;
			}

			exchange(x, y);  // CBA
			tempVar = varMap[invVarMap[x]];
			varMap[invVarMap[x]] = varMap[invVarMap[y]];
			varMap[invVarMap[y]] = tempVar;
			if (min > activeNodeCount) {
				best = 4;
				min = activeNodeCount;
			}

			exchange(y, z);  // CAB
			tempVar = varMap[invVarMap[z]];
			varMap[invVarMap[z]] = varMap[invVarMap[y]];
			varMap[invVarMap[y]] = tempVar;
			if (min > activeNodeCount) {
				best = 5;
				min = activeNodeCount;
			}

			exchange(x, y);  // ACB
			tempVar = varMap[invVarMap[x]];
			varMap[invVarMap[x]] = varMap[invVarMap[y]];
			varMap[invVarMap[y]] = tempVar;
			if (min > activeNodeCount) {
				best = 6;
				min = activeNodeCount;
			}

			switch (best) {
				case 3:  // BCA
					exchange(y, z);
					tempVar = varMap[invVarMap[z]];
					varMap[invVarMap[z]] = varMap[invVarMap[y]];
					varMap[invVarMap[y]] = tempVar;
					break;
				case 4:  // CBA
					exchange(x, y);
					tempVar = varMap[invVarMap[x]];
					varMap[invVarMap[x]] = varMap[invVarMap[y]];
					varMap[invVarMap[y]] = tempVar;
					break;
				case 1:  // ABC
					exchange(y, z);
					tempVar = varMap[invVarMap[z]];
					varMap[invVarMap[z]] = varMap[invVarMap[y]];
					varMap[invVarMap[y]] = tempVar;
					break;
				case 6:  // ACB
					break;
				case 2:  // BAC
					exchange(y, z);
					tempVar = varMap[invVarMap[z]];
					varMap[invVarMap[z]] = varMap[invVarMap[y]];
					varMap[invVarMap[y]] = tempVar;
					break;
				case 5:  // CAB
					exchange(x, y);
					tempVar = varMap[invVarMap[x]];
					varMap[invVarMap[x]] = varMap[invVarMap[y]];
					varMap[invVarMap[y]] = tempVar;
					break;
				default:
					break;
			}
		}
		return in;
	}

}
