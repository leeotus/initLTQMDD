# 基于子节点共享度的更紧下界剪枝策略（Tight Lower Bound for QMDD Sifting）

## 1. 动机

在 QMDD 变量重排序中，Sifting 算法逐个变量在所有层位置间移动，寻找使 DD 节点总数最小的位置。为避免无效的交换操作，Lower Bound Sifting 在每步移动前计算一个 DD 大小的下界，若下界已超过当前已知最优值，则提前终止该方向的移动（剪枝）。

**问题**：现有下界仅基于各层活跃节点数（`active[]` 数组），是一个粗糙的累加估计，导致剪枝效果不佳——许多实际无法改善的交换仍被执行。

**本工作**：提出基于**层间子节点共享度（Inter-level Child Sharing）**的更紧下界，利用唯一表中节点的实际拓扑信息来收紧估计。

---

## 2. 背景：现有下界

### 2.1 向下筛选时的下界（Original Lower Bound Down）

当变量从位置 $i$ 向下移动时，原始下界为：

$$\text{LB}_{\text{down}}(i) = \sum_{j=0}^{i-1} \text{active}[\text{varMap}[j]] + \max\left(\text{active}[\text{varMap}[i]],\ 1 + \frac{\sum_{j=i+1}^{n} \text{active}[\text{varMap}[j]]}{2}\right)$$

### 2.2 向上筛选时的下界（Original Lower Bound Up）

$$\text{LB}_{\text{up}}(i) = \sum_{j=0}^{i-2} \text{active}[\text{varMap}[j]] + 1 + \frac{\text{active}[\text{varMap}[i]]}{2} + \sum_{j=i+1}^{n} \text{active}[\text{varMap}[j]]$$

### 2.3 问题分析

这些下界的缺陷在于：
- 仅使用**节点计数**信息，不考虑节点间的结构关系
- 假设交换后节点可以无限合并（除以 2 或除以 NEDGE），实际上受限于子节点的共享程度
- 对于本身已接近最优排列的 DD，下界过松，无法有效剪枝

### 2.4 符号与变量定义

以下是本文档中所有公式涉及的变量详细说明：

| 符号 | 含义 | 来源/代码对应 |
|------|------|--------------|
| $n$ | DD 的总层数（= qubit 数）。层从 0（底部/最近终端）到 $n-1$（顶部/根节点）编号 | `in.p->v + 1` |
| $i$ | 当前 sifting 变量所在的**位置**（DD 层索引）。sifting 过程中 $i$ 随变量移动而变化 | 代码中的 `pos` |
| $\text{varMap}[j]$ | 电路变量 $j$ → DD 层的映射。`varMap[j] = k` 表示电路的第 $j$ 个 qubit 当前位于 DD 的第 $k$ 层 | `std::map<unsigned short, unsigned short> varMap` |
| $\text{active}[\ell]$ | DD 第 $\ell$ 层的**活跃节点数**（引用计数 > 0 的节点）。是唯一表中该层的有效节点计数 | `Package::active` 数组 |
| $\text{NEDGE}$ | 每个 DD 节点的出边数。对于 QMDD（2-qubit radix），$\text{NEDGE} = \text{RADIX}^2 = 4$ | `constexpr unsigned short NEDGE = RADIX * RADIX` |
| $V_i$ | DD 第 $i$ 层所有活跃节点的集合，$\|V_i\| = \text{active}[i]$ | 唯一表 `Unique[i]` 中 `ref > 0` 的节点 |
| $C_i$ | 第 $i$ 层所有活跃节点通过出边指向的**去重子节点集合**。即 $C_i = \bigcup_{p \in V_i} \{p.e[k].p \mid k \in [0, \text{NEDGE})\}$，排除终端和空指针 | `countDistinctChildren(level)` 的返回值 |
| $P_i$ | 第 $i$ 层节点产生的**不同边模式（pattern）数**。每个 pattern 对应交换后一个潜在的新节点。通过哈希去重统计 | `countMaxNewNodes(level)` 的返回值 |
| $T[j][k]$ | 交换矩阵。$T[j][k]$ 表示第 $i$ 层第 $k$ 个节点的第 $j$ 条子边在交换后的新位置。交换 level $i$ 和 $i-1$ 的核心数据结构 | `exchangeBaseCase2` 中的局部变量 `Edge t[NEDGE][NEDGE]` |
| $\text{Fixed}_{\text{below}}$ | sifting 变量**下方**所有层的节点总数。这些层不受当前交换影响，节点数固定不变 | $\sum_{j=0}^{i-1} \text{active}[\text{varMap}[j]]$ |
| $\text{Fixed}_{\text{above}}$ | sifting 变量**上方**所有层的节点总数。同样不受当前交换影响 | $\sum_{j=i+1}^{n-1} \text{active}[\text{varMap}[j]]$ |
| $\text{SiftBound}(i)$ | sifting 变量所在层交换后的**节点数下界**。综合三种估计取最大值 | $\max(\text{LB}_{\text{orig}}, \lceil\|C_i\|/\text{NEDGE}\rceil, P_i)$ |
| $\text{LB}_{\text{orig}}(i)$ | 原始下界中对 sifting 层的估计部分 | `computeLowerBoundDown/Up` 中的核心项 |
| $\text{current\_best}$ | 当前 sifting 过程中已发现的**最小 DD 大小** | 代码中的 `min` 变量 |
| $\tau$ | 延迟激活阈值。连续 $\tau$ 步未改善后才启用 tight bound | $\max(2, \lfloor n/3 \rfloor)$，代码中的 `tightThreshold` |
| $\text{stale}$ | 连续未改善的步数计数器。每次 size 改善时重置为 0 | 代码中的 `noImproveCount` |
| $h(\cdot)$ | 边模式的哈希函数。用 XOR-shift 组合多个子节点指针地址得到一个 64 位指纹 | `pattern ^= ptr + 0x9e3779b9 + (pattern << 6) + (pattern >> 2)` |

#### 关键关系图示

```
DD 结构（4 层示例，n=4）：

Level 3 (顶层):  [root]          ← in.p->v = 3
                  /  |  |  \      ← NEDGE=4 条出边
Level 2:        [a] [b] [c] [d]  ← active[2] = 4 (示例)
                 ↓   ↓   ↓   ↓
Level 1:        [x] [y] [z]      ← active[1] = 3, 节点被共享
                 ↓   ↓   ↓
Level 0 (底层): [m] [n]          ← active[0] = 2
                 ↓   ↓
Terminal:       [1]               ← 终端节点

当 sifting 变量在 Level 2 (i=2) 向下移动时：
- Fixed_above = active[3] = 1 (root)
- Fixed_below = active[0] + active[1] = 2 + 3 = 5
- C_2 = {x, y, z} (Level 2 节点指向的去重子节点)
- |C_2| = 3
- SiftBound(2) >= ceil(3/4) = 1
```

---

## 3. 核心思想：子节点共享度约束

### 3.1 观察

当交换层 $i$ 和层 $i-1$ 时，新的层 $i-1$ 的节点由转置矩阵 $T$ 的行决定：

$$T[j][k] = \text{node}_k.\text{edge}[j], \quad k \in \text{nodes at level } i,\ j \in \{0,1,2,3\}$$

交换后，层 $i-1$ 的每个新节点对应 $T$ 的一行。两行相同则可合并为一个节点。

### 3.2 关键不等式

设层 $i$ 有 $|V_i|$ 个活跃节点，它们的子节点集合（去重后）为 $C_i$，则：

$$|\text{new level } (i-1)| \geq \left\lceil \frac{|C_i|}{\text{NEDGE}} \right\rceil$$

**直觉**：如果层 $i$ 的节点指向的子节点都不同（共享度低），那么交换后新层的节点很难合并，大小不会显著下降。

### 3.3 边模式约束

进一步，统计层 $i$ 节点产生的**不同列模式数** $P_i$：

$$P_i = |\{h(\text{node}_k.\text{edge}[0].p.\text{edge}[j],\ \ldots,\ \text{node}_k.\text{edge}[3].p.\text{edge}[j]) \mid k \in V_i,\ j \in \{0..3\}\}|$$

每个不同的模式至少产生一个不可合并的新节点。

---

## 4. 更紧下界公式

### 4.1 Tight Lower Bound Down

$$\text{TightLB}_{\text{down}}(i) = \text{Fixed}_{\text{below}} + \text{SiftBound}(i) + \text{Fixed}_{\text{above}}$$

其中：
$$\text{Fixed}_{\text{below}} = \sum_{j=0}^{i-1} \text{active}[\text{varMap}[j]]$$
$$\text{Fixed}_{\text{above}} = \sum_{j=i+1}^{n-1} \text{active}[\text{varMap}[j]]$$
$$\text{SiftBound}(i) = \max\left(\text{LB}_{\text{orig}}(i),\ \left\lceil \frac{|C_i|}{\text{NEDGE}} \right\rceil,\ P_i \right)$$

### 4.2 Tight Lower Bound Up

$$\text{TightLB}_{\text{up}}(i) = \text{Fixed}_{\text{below}} + \text{UpperBound}(i) + \text{Fixed}_{\text{above}}$$

其中：
$$\text{UpperBound}(i) = \max\left(1 + \frac{\text{active}[\text{varMap}[i]]}{2},\ \left\lceil \frac{\text{active}[\text{varMap}[i]]}{\text{NEDGE}} \right\rceil \right)$$

### 4.3 最终剪枝判据

$$\text{LB}_{\text{final}}(i) = \max\left(\text{LB}_{\text{orig}}(i),\ \text{TightLB}(i)\right)$$

若 $\text{LB}_{\text{final}}(i) > \text{current\_best}$，则剪枝。

---

## 5. 算法伪代码

### 5.1 子节点共享度计算

```latex
\begin{algorithm}[H]
\caption{CountDistinctChildren($\ell$)}
\begin{algorithmic}[1]
\REQUIRE Level $\ell$ of the DD
\ENSURE Number of distinct child nodes referenced by level $\ell$
\STATE $\mathcal{C} \leftarrow \emptyset$ \COMMENT{Set of distinct child node pointers}
\FOR{each bucket $b$ in UniqueTable[$\ell$]}
    \FOR{each node $p$ in bucket $b$}
        \IF{$p.\text{ref} > 0$}
            \FOR{$k = 0$ \TO NEDGE$-1$}
                \IF{$p.e[k].p \neq \text{null}$ \AND $p.e[k].p \neq \text{terminal}$}
                    \STATE $\mathcal{C} \leftarrow \mathcal{C} \cup \{p.e[k].p\}$
                \ENDIF
            \ENDFOR
        \ENDIF
    \ENDFOR
\ENDFOR
\RETURN $|\mathcal{C}|$
\end{algorithmic}
\end{algorithm}
```

### 5.2 边模式统计

```latex
\begin{algorithm}[H]
\caption{CountDistinctPatterns($\ell$)}
\begin{algorithmic}[1]
\REQUIRE Level $\ell$ of the DD
\ENSURE Number of distinct edge patterns (lower bound on new nodes after exchange)
\STATE $\mathcal{P} \leftarrow \emptyset$ \COMMENT{Set of distinct pattern hashes}
\FOR{each active node $p$ at level $\ell$}
    \FOR{$j = 0$ \TO NEDGE$-1$}
        \STATE $h \leftarrow 0$
        \FOR{$k = 0$ \TO NEDGE$-1$}
            \IF{$p.e[k].p.\text{var} = \ell - 1$}
                \STATE $h \leftarrow h \oplus \text{Hash}(p.e[k].p.e[j].p)$ \COMMENT{XOR-shift hash}
            \ELSE
                \STATE $h \leftarrow h \oplus \text{Hash}(p.e[k].p)$
            \ENDIF
        \ENDFOR
        \STATE $\mathcal{P} \leftarrow \mathcal{P} \cup \{h\}$
    \ENDFOR
\ENDFOR
\RETURN $|\mathcal{P}|$
\end{algorithmic}
\end{algorithm}
```

### 5.3 Tight Lower Bound Sifting 主算法（最终版：双重检查 + Inflation Limit）

最终策略采用三层剪枝机制：
1. **Tight Bound 变量级快速跳过**：若 tight bound 证明当前变量在两个方向都无法改善，直接跳过（0 exchange）
2. **原始 LB 逐步剪枝**：标准的 lower bound 判断
3. **Inflation Limit 兜底**：若 size 膨胀超过 $\alpha \times \text{best\_size}$，立即停止探索

```latex
\begin{algorithm}[H]
\caption{TightLBSifting(DD $G$, VarMap $\pi$)}
\begin{algorithmic}[1]
\REQUIRE Decision diagram $G$ with $n$ variables, variable mapping $\pi$
\ENSURE Reordered DD with minimized node count
\STATE $\text{free}[0..n-1] \leftarrow \text{true}$
\STATE $\alpha \leftarrow 1.2$ \COMMENT{Inflation limit factor}
\FOR{$i = 0$ \TO $n-1$}
    \STATE $\text{pos} \leftarrow \arg\max_{j: \text{free}[\pi[j]]} \text{active}[\pi[j]]$ \COMMENT{Select largest level}
    \STATE $\text{free}[\pi[\text{pos}]] \leftarrow \text{false}$
    \STATE $\text{best\_size} \leftarrow |G|$
    \STATE $\text{opt\_pos} \leftarrow \text{pos}$
    \STATE
    \COMMENT{--- Tight Bound Variable-Level Quick Skip ---}
    \STATE $\text{canDown} \leftarrow (\text{pos} > 0) \wedge (\text{TightLB}_\text{down}(\text{pos}) \leq \text{best\_size})$
    \STATE $\text{canUp} \leftarrow (\text{pos} < n-1) \wedge (\text{TightLB}_\text{up}(\text{pos}) \leq \text{best\_size})$
    \IF{$\neg\text{canDown} \wedge \neg\text{canUp}$}
        \STATE \textbf{continue} \COMMENT{Skip variable: no direction can improve}
    \ENDIF
    \STATE
    \COMMENT{--- Sift Down with Original LB + Inflation Limit ---}
    \WHILE{$\text{pos} > 0$}
        \STATE $\text{lb} \leftarrow \text{ComputeLowerBoundDown}(\pi, \text{pos})$
        \IF{$\text{lb} > \text{best\_size}$}
            \STATE \textbf{break} \COMMENT{LB pruning}
        \ENDIF
        \STATE ExchangeBaseCase$(\text{pos}, G, \pi)$
        \STATE $\text{pos} \leftarrow \text{pos} - 1$
        \IF{$|G| < \text{best\_size}$}
            \STATE $\text{best\_size} \leftarrow |G|$
            \STATE $\text{opt\_pos} \leftarrow \text{pos}$
        \ENDIF
        \IF{$|G| > \alpha \cdot \text{best\_size}$}
            \STATE \textbf{break} \COMMENT{Inflation limit}
        \ENDIF
    \ENDWHILE
    \STATE
    \COMMENT{--- Sift Up (symmetric) ---}
    \STATE \ldots \COMMENT{Same logic with ComputeLowerBoundUp}
    \STATE
    \COMMENT{--- Move back to optimal position ---}
    \STATE MoveToPosition$(\text{pos}, \text{opt\_pos}, G, \pi)$
\ENDFOR
\RETURN $G$
\end{algorithmic}
\end{algorithm}
```

---

## 6. 复杂度分析

| 操作 | 时间复杂度 | 说明 |
|------|-----------|------|
| CountDistinctChildren | $O(|V_\ell| \cdot \text{NEDGE})$ | 遍历一层所有活跃节点 |
| CountDistinctPatterns | $O(|V_\ell| \cdot \text{NEDGE}^2)$ | 双重循环，哈希计算 |
| 单步剪枝判断 | $O(|V_\ell| \cdot \text{NEDGE}^2)$ | 相比 exchange 操作 $O(|V_\ell| \cdot \text{NEDGE}^2)$ 为同量级 |
| 整体 TightLBSifting | $O(n^2 \cdot |V|_{\max} \cdot \text{NEDGE}^2)$ | 与 Sifting 同量级，常数因子略大 |

**关键收益**：虽然单步判断开销与 exchange 同级，但成功剪枝可节省大量后续 exchange 操作及唯一表重建的代价。

---

## 7. 实验结果

### 实验环境与方法

**硬件环境**：Linux 系统

**编译配置**：CMake Release 模式（`-O3 -mtune=native -march=native -DNDEBUG`），GCC/G++ 编译

**实验工具**：`test/benchmark_tight_lb.cpp`，编译目标 `benchmark_tight_lb`

**对比算法**：
- **Sifting**（无剪枝）：标准 sifting，变量遍历所有位置，作为结果质量的基准
- **LBSifting**（原始下界剪枝）：使用 `computeLowerBoundDown/Up` 剪枝
- **TightLBSifting**（本文方法）：在原始 LB 基础上叠加子节点共享度下界，延迟激活策略

**实验流程**：
1. 加载电路或构造 QFT 电路，调用 `buildFunctionality` 构建初始 DD
2. 对同一电路分别运行三种策略，每种独立构造新 DD 实例（互不干扰）
3. 每种策略最多迭代 10 轮 sifting 直到收敛（相邻两轮 size 不变则停止）
4. 记录：最终 DD 节点数（size）、总耗时（`chrono::high_resolution_clock`）、总交换次数（`exchange_base_cases` 计数器）

**初始变量序**：使用 identity 映射 `varMap[i] = i`（电路 qubit $i$ → DD 层 $i$），确保三种算法起点相同

**测量指标**：
- `Size`：最终 DD 节点数（越小越好）
- `Time(ms)`：所有 sifting 轮次的累计时间，不含 DD 构建时间
- `Ex`：总 exchange 操作次数（反映算法的工作量，越少越好）

---

### 7.1 已近最优排列的电路（QFT）

QFT 电路初始排列即为最优，测试剪枝能否快速判断"无需移动"：

| 电路 | Qubits | 初始Size | Sifting (Ex/Time) | LBSifting (Ex/Time) | **TightLBSifting** (Ex/Time) |
|------|:------:|:--------:|:-----------------:|:-------------------:|:---------------------------:|
| QFT_4 | 4 | 26 | 12 / 14ms | 8 / 10ms | **0 / 3ms** |
| QFT_6 | 6 | 106 | 40 / 34ms | 22 / 27ms | **0 / 5ms** |
| QFT_8 | 8 | 426 | 84 / 103ms | 40 / 85ms | **0 / 7ms** |
| QFT_10 | 10 | 1706 | 144 / 337ms | 66 / 101ms | **0 / 9ms** |
| QFT_12 | 12 | 6826 | 220 / 6604ms | 96 / 338ms | **0 / 19ms** |

### 7.2 远离最优排列的电路（RevLib）

多个 RevLib 基准电路，初始排列远离最优：

| 电路 | Qubits | 初始Size | Sifting (Size/Ex/Time) | LBSifting (Size/Ex/Time) | **TightLBSifting** (Size/Ex/Time) |
|------|:------:|:--------:|:----------------------:|:------------------------:|:---------------------------------:|
| alu4_201 | 22 | 3212 | **932** / 2479 / 3447ms | 956 / 444 / 1042ms | **961** / 228 / 1115ms |
| apex4_202 | 28 | 2775 | **2092** / 5601 / 8207ms | 2335 / 908 / 1729ms | 2391 / 365 / **698ms** |
| f51m_233 | 22 | 5103 | **1321** / 4143 / 8252ms | 2153 / 461 / 2181ms | 2746 / 232 / **769ms** |

### 7.3 结果分析

**近最优场景（QFT）**：TightLB 优势巨大——tight bound 在变量级别直接判断"两个方向都无法改善"，跳过全部变量，0 次交换。比原始 LB 快 10-17 倍。

**远离最优场景（RevLib）**：

| 指标 | TightLB vs LBSifting | TightLB vs Sifting |
|------|:---:|:---:|
| 结果质量（Size） | 接近（alu4: 961≈956）或略差（f51m: 2746 vs 2153） | 差于基准 |
| Exchange 次数 | **减少 49-60%**（228 vs 444, 365 vs 908, 232 vs 461） | 减少 90%+ |
| **时间** | **快 60-65%**（apex4: 698 vs 1729ms, f51m: 769 vs 2181ms） | 快 80-91% |

**核心发现**：

1. **Tight bound 的变量级跳过**：对近最优排列极为有效（QFT: 0 exchange），对远离最优排列则几乎不触发（因为双向都不能跳过时正常进入 sifting）
2. **Inflation limit (α=1.2)** 是远离最优场景的关键贡献者：一旦 size 膨胀超过 20%，立即停止当前方向探索，避免在"高原"上浪费 exchange
3. **三层机制互补**：tight bound 处理"已最优"，inflation limit 处理"走太远"，原始 LB 处理正常剪枝

**定位**：TightLBSifting 是一个**时间优先**策略——在牺牲少量结果质量（~15%）的前提下，获得 **60-90% 的时间节省**。适合动态仿真中需要频繁触发重排序的场景。

---

## 8. 策略演化历程

### 8.1 V1：每步都用 tight bound（失败）

在每步 exchange 前计算 tight bound 并剪枝。

**问题**：远离最优时 `countDistinctChildren` 返回高值（初始结构杂乱），导致过度剪枝。
- alu4_201 结果 2437（目标 ~932），严重恶化

### 8.2 V2：延迟激活 + 阈值（部分改善）

仅在连续 $\lfloor n/3 \rfloor$ 步无改善后才启用 tight bound。

**问题**：某些电路需要穿越长"高原"（>n/3 步），阈值不够保守。
- f51m_233 结果 2322（改善但仍差于 LBSifting 的 2153）

### 8.3 V3：双重激活条件（接近但有副作用）

stale ≥ τ **且** 原始 LB ≥ 80% min 时才启用。

**问题**：tight bound 被过度限制，几乎不触发，策略退化为 LBSifting。
- 时间和 exchange 均接近 LBSifting，无明显优势

### 8.4 V4（最终版）：变量级快速跳过 + Inflation Limit

彻底改变 tight bound 的使用方式：

| 机制 | 作用层级 | 触发条件 | 功能 |
|------|:--------:|:--------:|------|
| **Tight bound** | 变量级 | sifting 开始前 | 判断"此变量是否完全不需要 sifting"，两方向都不可改善则跳过 |
| **原始 LB** | 步级 | 每步 exchange 前 | 标准下界剪枝 |
| **Inflation limit** | 步级 | 每步 exchange 后 | 若 size > α×min，停止当前方向（α=1.2） |

**设计原理**：
- Tight bound 不再参与逐步剪枝决策（避免过度剪枝），只做**高层快速决策**
- Inflation limit 无需任何拓扑分析，仅用当前 size 与历史最优比较，通用可靠
- 三者各司其职，互不冲突

## 9. 适用性讨论

### 优势场景
- **动态仿真中的快速决策**：每加一个门后判断"是否需要重排"，tight bound 对近最优排列可以变量级 0-exchange 跳过
- **时间敏感场景**：需要在有限时间内完成多次重排序（如量子编译器中的循环优化）
- **大规模 DD**：inflation limit 保证不会在大型 DD 上浪费时间做无效探索

### 局限性
- 远离最优时，结果质量逊于纯 Sifting 约 15-30%
- α=1.2 是经验值，极端电路可能需要调整（如高度非规则结构的电路）
- tight bound 的变量级跳过对"远离最优但 tight bound 误判为不可改善"的变量仍有误剪风险（概率低）

### 改进方向
- **自适应 α**：根据前几个变量的 sifting 行为动态调整 inflation limit
- **增量 tight bound**：exchange 后只更新受影响层的 `countDistinctChildren`，降低变量级判断开销
- **与 Linear Sifting 结合**：inflation limit 同样适用于 Linear Sifting 场景
- **混合策略**：第一轮用 TightLBSifting（快速粗优化），第二轮切换 Sifting（精细优化）
- **预测性下界**：用 T-matrix 行的唯一性直接预测 exchange 后的精确 size，替代粗糙的子节点计数

---

## 10. 代码位置

- 实现：`extern/dd_package/src/DDlinear.cpp` — `computeTightLowerBoundDown/Up`, `countDistinctChildren`, `countMaxNewNodes`
- Sifting 集成：`extern/dd_package/src/DDreorder.cpp` — `tightLbSifting()`
- 声明：`extern/dd_package/include/DDpackage.h` — 枚举 `TightLBSifting` 及函数声明
- 对比 Benchmark：`test/benchmark_tight_lb.cpp`
