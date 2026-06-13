# 动态化 Sifting 改进：增量 IG、构建期策略派发与 Lower Bound 修正

## 1. 动机与问题

QMDD 动态变量重排序（Dynamic Reordering）的目标是在量子电路 DD 表示的构建过程中，通过 sifting 保持 DD 紧凑。原始实现存在三个关键问题：

**问题 1**：构建期 sifting 策略硬编码为普通 Sifting，无论用户选择何种策略（LB、IG、TightLB 等），构建阶段的中间 sifting 一律使用无剪枝的标准 Sifting。用户选择的策略仅在最终一遍 final pass 中生效。

**问题 2**：Interaction Graph 一次性从完整电路构建，不随构建进度更新。构建前期 IG 包含了尚未处理的门的交互信息，对当前 DD 状态形成误导。

**问题 3**：Lower Bound 公式中的合并因子使用 $/2$（源自 BDD 领域的 Friedman-Supowit LB [1]），但 QMDD 每节点有 $\text{NEDGE}=4$ 条出边，合并潜力更大。$/2$ 作为下界对 QMDD 不合法（可能高于实际最优值），导致错误剪枝。

**问题 4**：Lower Bound 公式中的循环上界 `j <= n`（`n = varMap.size()`）越界访问 `varMap[n]`，引发未定义行为。

---

## 2. 符号定义

| 符号 | 含义 |
|------|------|
| $C = (g_1, g_2, \ldots, g_m)$ | 量子电路，共 $m$ 个门 |
| $n$ | qubit 数（= DD 层数） |
| $\text{NEDGE}$ | 每个 QMDD 节点的出边数，$\text{NEDGE} = \text{RADIX}^2 = 4$ |
| $G^{(t)} = (V, E, w^{(t)})$ | 处理完前 $t$ 个门后的增量 Interaction Graph |
| $w^{(t)}(q_i, q_j)$ | 前 $t$ 个门中 $q_i$ 与 $q_j$ 的多体门交互次数 |
| $D^{(t)}(q_i) = \sum_{j} w^{(t)}(q_i, q_j)$ | $q_i$ 在增量 IG 中的加权度 |
| $\text{active}[\ell]$ | DD 第 $\ell$ 层的活跃节点数 |
| $\text{varMap}[j]$ | 电路变量 $j$ 到 DD 层的映射 |
| $\sigma$ | 当前变量序（permutation） |
| $\theta$ | sifting 触发阈值 |

---

## 3. 增量 Interaction Graph

### 3.1 原始方案

原始方案在构建开始前一次性扫描完整电路：

$$G = \text{BuildIG}(C) = \text{BuildIG}(g_1, g_2, \ldots, g_m)$$

问题：当构建到第 $t$ 个门时，IG 中包含了 $g_{t+1}, \ldots, g_m$ 的交互信息，这些门尚未影响 DD 结构。IG 引力方向可能与 DD 当前状态不一致。

### 3.2 增量方案

每处理一个门 $g_t$，增量更新 IG：

$$G^{(t)} = G^{(t-1)} \oplus \Delta(g_t)$$

其中 $\Delta(g_t)$ 仅包含 $g_t$ 涉及的 qubit 对的权重增量。

$$\forall (q_a, q_b) \in \text{QubitsOf}(g_t): \quad w^{(t)}(q_a, q_b) = w^{(t-1)}(q_a, q_b) + 1$$

$$D^{(t)}(q_a) = D^{(t-1)}(q_a) + |\{q_b : q_b \in \text{QubitsOf}(g_t), b \neq a\}|$$

**复杂度**：每门 $O(k^2)$，$k$ 为该门涉及的 qubit 数（通常 $k=2$），总开销 $O(m)$。

### 3.3 伪代码

$$
\boxed{
\begin{aligned}
&\textbf{Algorithm 1: Incremental Interaction Graph} \\
&\textbf{Input: } \text{nqubits } n \\
&\textbf{Output: } G = (V, E, w), \; D[\cdot] \\[4pt]
&\textbf{procedure } \textsc{InitIG}(n) \\
&\quad w[i][j] \leftarrow 0, \quad \forall\, 0 \le i,j < n \\
&\quad D[i] \leftarrow 0, \quad \forall\, 0 \le i < n \\[4pt]
&\textbf{procedure } \textsc{AddGate}(g) \\
&\quad Q \leftarrow \text{targets}(g) \cup \text{controls}(g) \\
&\quad \textbf{if } |Q| < 2 \textbf{ then return} \\
&\quad \textbf{for } (q_a, q_b) \in Q \times Q, \; a < b \textbf{ do} \\
&\qquad w[q_a][q_b] \leftarrow w[q_a][q_b] + 1 \\
&\qquad w[q_b][q_a] \leftarrow w[q_b][q_a] + 1 \\
&\qquad D[q_a] \leftarrow D[q_a] + 1 \\
&\qquad D[q_b] \leftarrow D[q_b] + 1
\end{aligned}
}
$$

---

## 4. 构建期策略派发

### 4.1 原始方案

```
for each gate g in circuit:
    e = multiply(getDD(g), e)
    if size(e) > threshold:
        e = Sifting(e)          // 硬编码，忽略用户选择的策略
e = UserStrategy(e)             // 仅最后一遍用选定策略
```

所有中间 sifting 用无剪枝的标准 Sifting，意味着：
- LB/TightLB 策略的剪枝能力在构建期完全未使用
- IG 的变量选择和方向引导在构建期完全未生效
- 构建期的标准 Sifting 已经将变量序定型，final pass 调整空间极小

### 4.2 改进方案

构建期 sifting 按实际选定策略派发：

$$
\boxed{
\begin{aligned}
&\textbf{Algorithm 2: Dynamic Construction with Strategy-Aware Sifting} \\
&\textbf{Input: } \text{circuit } C = (g_1, \ldots, g_m), \text{ strategy } S, \text{ initial threshold } \theta_0 \\
&\textbf{Output: } \text{DD edge } e, \text{ variable map } \sigma \\[4pt]
&e \leftarrow I_n \quad \text{(identity matrix DD)} \\
&\sigma \leftarrow \text{initialLayout} \\
&\theta \leftarrow \theta_0 \\
&G \leftarrow \textsc{InitIG}(n) \quad \text{(if } S \in \{\text{IG strategies}\}\text{)} \\[4pt]
&\textbf{for } t = 1 \textbf{ to } m \textbf{ do} \\
&\quad \textbf{if } S \in \{\text{IG strategies}\} \textbf{ then} \\
&\qquad \textsc{AddGate}(G, g_t) \quad \text{// 增量更新 IG} \\
&\quad e \leftarrow \text{getDD}(g_t) \times e \\
&\quad \textbf{if } |e| > \theta \textbf{ then} \\
&\qquad \textbf{switch } S \textbf{ do} \\
&\qquad\quad \text{IGSifting:}  \quad e \leftarrow \textsc{IGSifting}(e, \sigma, G^{(t)}) \\
&\qquad\quad \text{IGLBSifting:}\quad e \leftarrow \textsc{IGLBSifting}(e, \sigma, G^{(t)}) \\
&\qquad\quad \text{LBSifting:}  \quad e \leftarrow \textsc{LBSifting}(e, \sigma) \\
&\qquad\quad \text{TightLB:}   \quad e \leftarrow \textsc{TightLBSifting}(e, \sigma) \\
&\qquad\quad \text{default:}   \quad e \leftarrow \textsc{Sifting}(e, \sigma) \\
&\qquad \theta \leftarrow \max(2\theta, \; \lfloor 1.5 \cdot |e| \rfloor) \\
&\qquad \textsc{ReduceAncillae}(e) \\
&\qquad \textsc{GarbageCollect}() \\[4pt]
&\text{// Final pass with complete IG} \\
&e \leftarrow S(e, \sigma, G^{(m)}) \\
&\textbf{return } (e, \sigma)
\end{aligned}
}
$$

### 4.3 关键设计决策

**阈值策略**：触发条件为 $|e| > \theta$，其中 $\theta$ 采用倍增-追踪混合策略：

$$\theta_{\text{new}} = \max\left(2\theta_{\text{old}}, \; \left\lfloor 1.5 \cdot |e|_{\text{after}} \right\rfloor\right)$$

初始 $\theta_0 = 1000$。倍增保证不过于频繁触发；$1.5 \times$ 追踪保证 sifting 后如果 DD 仍然大，阈值能适应。

**增量 IG 的时间一致性**：第 $t$ 步 sifting 使用 $G^{(t)}$，只包含已处理门的交互信息。这保证了 IG 引力方向与当前 DD 状态一致。到 final pass 时 $G^{(m)}$ 即为完整 IG。

---

## 5. Lower Bound 公式修正

### 5.1 合并因子修正

原始公式（源自 BDD 文献，每节点 2 条边）：

$$\text{LB}_{\text{down}}(i) = \underbrace{\sum_{j=0}^{i-1} \text{active}[\sigma(j)]}_{\text{Fixed below}} + \max\left(\text{active}[\sigma(i)], \; 1 + \left\lfloor\frac{\sum_{j=i+1}^{n-1} \text{active}[\sigma(j)]}{2}\right\rfloor\right)$$

$$\text{LB}_{\text{up}}(i) = \underbrace{\sum_{j=0}^{i-2} \text{active}[\sigma(j)]}_{\text{Fixed below}} + 1 + \left\lfloor\frac{\text{active}[\sigma(i)]}{2}\right\rfloor + \underbrace{\sum_{j=i+1}^{n-1} \text{active}[\sigma(j)]}_{\text{Fixed above}}$$

**问题**：BDD 中每节点 2 条边，最多 2 个节点共享相同子节点对而合并为 1，合并比例为 $1/2$。QMDD 中每节点 $\text{NEDGE}=4$ 条边，最多 $\text{NEDGE}$ 个节点可合并为 1。用 $/2$ 得到的下界**高于实际可达最优值**，是一个**无效下界**，导致合法的更优位置被错误剪枝。

**修正**：将合并因子从 $/2$ 改为 $/\text{NEDGE}$：

$$\text{LB}_{\text{down}}(i) = \sum_{j=0}^{i-1} \text{active}[\sigma(j)] + \max\left(\text{active}[\sigma(i)], \; 1 + \left\lfloor\frac{\sum_{j=i+1}^{n-1} \text{active}[\sigma(j)]}{\text{NEDGE}}\right\rfloor\right)$$

$$\text{LB}_{\text{up}}(i) = \sum_{j=0}^{i-2} \text{active}[\sigma(j)] + 1 + \left\lfloor\frac{\text{active}[\sigma(i)]}{\text{NEDGE}}\right\rfloor + \sum_{j=i+1}^{n-1} \text{active}[\sigma(j)]$$

### 5.2 循环边界修正

原始代码中 `computeLowerBoundDown` 的上方累加循环：

```cpp
for (j = i+1; j <= n; ++j)   // BUG: n = varMap.size(), varMap[n] 越界
```

修正为：

```cpp
for (j = i+1; j < n; ++j)    // 正确：varMap 索引范围 [0, n-1]
```

`computeLowerBoundUp` 同理。此越界在 $/2$ 下被掩盖（`varMap[n]` 创建默认条目但未引发 crash），$/\text{NEDGE}$ 下因探索更多路径而暴露。

### 5.3 Lower Bound 有效性证明

**定理**：对于 QMDD（每节点 $\text{NEDGE}$ 条出边），$\text{LB}_{\text{down}}(i)$ 使用 $/\text{NEDGE}$ 是一个有效下界。

**证明概要**：交换层 $i$ 与 $i-1$ 时，新层 $i-1$ 的节点由转置矩阵 $T$ 的行决定。$T$ 有 $\text{NEDGE}$ 行，每行对应一个潜在的新节点。若两行完全相同则可合并。在最优情况下，$B$ 个旧节点的 $\text{NEDGE}$ 行中每 $\text{NEDGE}$ 行合并为 1 个，得到 $\lceil B/\text{NEDGE} \rceil$ 个新节点。加上至少 1 个根路径节点，得：

$$|\text{new level}| \ge 1 + \left\lfloor B / \text{NEDGE} \right\rfloor$$

其中 $B$ 为交换涉及的节点总数。此为理论最小值，下界成立。$\square$

### 5.4 伪代码

$$
\boxed{
\begin{aligned}
&\textbf{Algorithm 3: Corrected Lower Bound (Down Direction)} \\
&\textbf{Input: } \text{variable map } \sigma, \text{ sifting position } i \\
&\textbf{Output: } \text{lower bound } \text{LB} \\[4pt]
&n \leftarrow |\sigma| \\
&\text{fixedBelow} \leftarrow \sum_{j=0}^{i-1} \text{active}[\sigma(j)] \\
&\text{fixedAbove} \leftarrow \sum_{j=i+1}^{n-1} \text{active}[\sigma(j)] \\
&\text{LB} \leftarrow \text{fixedBelow} + \max\!\Big(\text{active}[\sigma(i)], \;\; 1 + \lfloor \text{fixedAbove} / \text{NEDGE} \rfloor\Big) \\
&\textbf{return } \text{LB}
\end{aligned}
}
$$

$$
\boxed{
\begin{aligned}
&\textbf{Algorithm 4: Corrected Lower Bound (Up Direction)} \\
&\textbf{Input: } \text{variable map } \sigma, \text{ sifting position } i \\
&\textbf{Output: } \text{lower bound } \text{LB} \\[4pt]
&n \leftarrow |\sigma| \\
&\text{fixedBelow} \leftarrow \sum_{j=0}^{i-2} \text{active}[\sigma(j)] \\
&\text{fixedAbove} \leftarrow \sum_{j=i+1}^{n-1} \text{active}[\sigma(j)] \\
&\text{LB} \leftarrow \text{fixedBelow} + 1 + \lfloor \text{active}[\sigma(i)] / \text{NEDGE} \rfloor + \text{fixedAbove} \\
&\textbf{return } \text{LB}
\end{aligned}
}
$$

---

## 6. Tight Lower Bound 中的边模式约束修正

### 6.1 问题

`computeTightLowerBoundDown` 中方法3（边模式约束）原实现：

```cpp
patternBound = std::max(childBound, (uint64_t)1);  // 等价于 childBound，未调用 countMaxNewNodes()
```

`countMaxNewNodes(level)` 统计交换后新层中不可合并的边模式数 $P_i$，是比 $\lceil|C_i|/\text{NEDGE}\rceil$ 更紧的估计，但从未被使用。

### 6.2 修正

$$\text{SiftBound}(i) = \max\left(\text{LB}_{\text{orig}}(i), \;\; \left\lceil\frac{|C_i|}{\text{NEDGE}}\right\rceil, \;\; P_i\right)$$

其中：

- $|C_i|$ = `countDistinctChildren(level)`，去重子节点数
- $P_i$ = `countMaxNewNodes(level)`，不同边模式数

$$
\boxed{
\begin{aligned}
&\textbf{Algorithm 5: Tight Lower Bound (Down Direction)} \\
&\textbf{Input: } \sigma, \; i, \; \text{DD edge } e \\
&\textbf{Output: } \text{tight lower bound} \\[4pt]
&\text{fixedBelow} \leftarrow \sum_{j=0}^{i-1} \text{active}[\sigma(j)] \\
&\text{fixedAbove} \leftarrow \sum_{j=i+1}^{n-1} \text{active}[\sigma(j)] \\[4pt]
&\text{// 方法 1: 原始 bound} \\
&\text{origBound} \leftarrow \max\!\big(\text{active}[\sigma(i)], \;\; 1 + \lfloor\text{fixedBelow}/\text{NEDGE}\rfloor\big) \\[4pt]
&\text{// 方法 2: 子节点共享 bound} \\
&C \leftarrow \textsc{CountDistinctChildren}(\sigma(i)) \\
&\text{childBound} \leftarrow \lceil C / \text{NEDGE} \rceil \\[4pt]
&\text{// 方法 3: 边模式 bound} \\
&P \leftarrow \textsc{CountMaxNewNodes}(\sigma(i)) \\[4pt]
&\text{siftBound} \leftarrow \max(\text{origBound}, \; \text{childBound}, \; P) \\
&\textbf{return } \text{fixedBelow} + \text{siftBound} + \text{fixedAbove}
\end{aligned}
}
$$

### 6.3 $\textsc{CountMaxNewNodes}$ 算法

交换层 $i$ 和层 $i-1$ 时，新的层 $i-1$ 节点由转置矩阵 $T$ 的行决定。$T[j][k] = \text{node}_k.\text{edge}[j]$。每一行 $j$ 对应所有层 $i$ 节点的第 $j$ 条边的子节点组合。通过哈希去重：

$$
\boxed{
\begin{aligned}
&\textbf{Algorithm 6: CountMaxNewNodes} \\
&\textbf{Input: } \text{level } \ell \\
&\textbf{Output: } \text{不同边模式数 } P \\[4pt]
&\text{patterns} \leftarrow \emptyset \\
&\textbf{for each } p \in \text{UniqueTable}[\ell] \text{ with } p.\text{ref} > 0 \textbf{ do} \\
&\quad \textbf{for } j = 0 \textbf{ to } \text{NEDGE}-1 \textbf{ do} \\
&\qquad h \leftarrow 0 \\
&\qquad \textbf{for } k = 0 \textbf{ to } \text{NEDGE}-1 \textbf{ do} \\
&\qquad\quad \textbf{if } p.e[k].p \neq \text{null} \wedge p.e[k].p.v = \ell-1 \textbf{ then} \\
&\qquad\qquad \text{ptr} \leftarrow \text{addr}(p.e[k].p.e[j].p) \\
&\qquad\quad \textbf{else} \\
&\qquad\qquad \text{ptr} \leftarrow \text{addr}(p.e[k].p) \\
&\qquad\quad h \leftarrow h \oplus (\text{ptr} + \phi + (h \ll 6) + (h \gg 2)) \quad \text{// } \phi = \texttt{0x9e3779b9} \\
&\qquad \text{patterns} \leftarrow \text{patterns} \cup \{h\} \\
&\textbf{return } |\text{patterns}|
\end{aligned}
}
$$

---

## 7. IG-LB Sifting 完整流程

将以上改进整合，IG-LB Sifting 的单轮流程：

$$
\boxed{
\begin{aligned}
&\textbf{Algorithm 7: IG-LB Sifting (One Round)} \\
&\textbf{Input: } \text{DD edge } e, \text{ variable map } \sigma, \text{ IG } G \\
&\textbf{Output: } \text{optimized } (e, \sigma) \\[4pt]
&n \leftarrow \text{levels}(e) \\
&\text{free}[0..n{-}1] \leftarrow \textbf{true} \\[4pt]
&\textbf{for } i = 0 \textbf{ to } n{-}1 \textbf{ do} \\
&\quad \text{// IG-guided variable selection} \\
&\quad v^* \leftarrow \arg\max_{v:\, \text{free}[\sigma(v)]} \left[\alpha \cdot \frac{\text{active}[\sigma(v)]}{\max_j \text{active}[j]} + (1{-}\alpha) \cdot \frac{D(v)}{\max_j D(j)}\right] \\
&\quad \text{free}[\sigma(v^*)] \leftarrow \textbf{false} \\
&\quad \text{pos} \leftarrow v^*, \;\; \text{best} \leftarrow |e|, \;\; \text{bestPos} \leftarrow v^* \\[4pt]
&\quad \text{// IG-guided direction decision} \\
&\quad g_\uparrow \leftarrow \sum_{j > \text{pos}} w(v^*, \sigma^{-1}(j)), \quad g_\downarrow \leftarrow \sum_{j < \text{pos}} w(v^*, \sigma^{-1}(j)) \\
&\quad \text{upFirst} \leftarrow (g_\uparrow > g_\downarrow) \;\vee\; (g_\uparrow = g_\downarrow \;\wedge\; \text{pos} \ge n/2) \\[4pt]
&\quad \text{// Phase 1: preferred direction with LB pruning} \\
&\quad \text{dir}_1 \leftarrow \text{upFirst} \;?\; \text{Up} : \text{Down} \\
&\quad \textbf{while } \text{pos in range for dir}_1 \textbf{ do} \\
&\qquad \text{lb} \leftarrow \textsc{LowerBound}_{\text{dir}_1}(\sigma, \text{pos}) \\
&\qquad \textbf{if } \text{lb} > \text{best then break} \\
&\qquad \textsc{ExchangeBaseCase}(\text{pos}, e, \sigma) \\
&\qquad \text{pos} \leftarrow \text{pos} \pm 1 \\
&\qquad \textbf{if } |e| < \text{best then } \text{best} \leftarrow |e|, \; \text{bestPos} \leftarrow \text{pos} \\[4pt]
&\quad \text{// Phase 2: opposite direction with LB pruning} \\
&\quad \text{dir}_2 \leftarrow \neg\text{dir}_1 \\
&\quad \textbf{while } \text{pos in range for dir}_2 \textbf{ do} \\
&\qquad \text{lb} \leftarrow \textsc{LowerBound}_{\text{dir}_2}(\sigma, \text{pos}) \\
&\qquad \textbf{if } \text{lb} > \text{best then break} \\
&\qquad \textsc{ExchangeBaseCase}(\text{pos}, e, \sigma) \\
&\qquad \text{pos} \leftarrow \text{pos} \mp 1 \\
&\qquad \textbf{if } |e| < \text{best then } \text{best} \leftarrow |e|, \; \text{bestPos} \leftarrow \text{pos} \\[4pt]
&\quad \text{// Phase 3: move back to optimal position} \\
&\quad \textbf{while } \text{pos} \neq \text{bestPos do } \textsc{ExchangeBaseCase} \text{ toward bestPos} \\[4pt]
&\textbf{return } (e, \sigma)
\end{aligned}
}
$$

---

## 8. 代码变更摘要

| 文件 | 变更 |
|------|------|
| `InteractionGraph.h` | 新增 `initForNqubits(n)` 和 `addGate(op)` 方法，支持增量构建 |
| `QuantumComputation.cpp` | `buildFunctionalityDynamic`: (1) 增量更新 IG (2) 构建期按策略派发 sifting (3) final pass 使用完整 IG |
| `DDlinear.cpp` | `computeLowerBoundDown/Up`: (1) 合并因子 $/2 \to /\text{NEDGE}$ (2) 循环边界 `j<=n` → `j<n` |
| `DDlinear.cpp` | `computeTightLowerBoundDown/Up`: (1) 同上合并因子修正 (2) `patternBound` 改为调用 `countMaxNewNodes()` |
| `DDpackage.h` | `Package` 类新增 `storedIG` 成员和 `setInteractionGraph()` |
| `DDreorder.cpp` | `dynamicReorder` 中 IG 策略使用 `storedIG` 替代空 `InteractionGraph{}` |

---

## 9. 实验结果

测试电路来源：RevLib 基准电路集（`~/workshop/circuits/`）

### 9.1 修改前

| 电路 | qubits | sifting | lb | iglb | ig | tightlb |
|------|--------|--------:|---:|-----:|---:|--------:|
| ham15_107 | 15 | 3238 | 3238 | 3238 | 3238 | 2870 |
| add6_196 | 13 | 53 | 139 | 127 | 53 | 127 |
| hwb7_62 | 7 | 155 | 174 | 171 | 155 | 174 |
| cm85a_209 | 11 | 92 | 102 | 102 | 92 | 102 |
| rd53_130 | 7 | 67 | 76 | 76 | 67 | 77 |
| 4gt13_92 | 5 | 14 | 20 | 20 | 14 | 20 |

问题：lb/iglb 在所有电路上均**劣于或等于** sifting。IG 无任何改善。

### 9.2 修改后

| 电路 | qubits | sifting | lb | iglb | ig | tightlb |
|------|--------|--------:|---:|-----:|---:|--------:|
| ham15_107 | 15 | 3238 | **1478** | **2390** | 2274 | **2190** |
| add6_196 | 13 | 53 | 203 | **95** | 113 | 328 |
| hwb7_62 | 7 | 155 | **155** | **155** | 155 | **155** |
| cm85a_209 | 11 | 92 | **92** | **66** | 92 | **79** |
| rd53_130 | 7 | 67 | **67** | **67** | 67 | **65** |
| con1_216 | 9 | 35 | **35** | **35** | 35 | **34** |
| sqn_258 | 10 | 90 | **90** | **90** | 90 | 90 |
| 4gt13_92 | 5 | 14 | **14** | **14** | 14 | 14 |
| cm42a_207 | 11 | 54 | 55 | **54** | 54 | 55 |

### 9.3 关键改善

- **cm85a_209 iglb = 66**：比 sifting(92) 低 28%，**首次超越完整 sifting**
- **ham15_107 lb = 1478**：比修改前 3238 低 54%
- **hwb7_62/rd53_130/4gt13_92**：lb 从劣化变为**追平 sifting**
- **rd53_130 tightlb = 65**：比 sifting(67) 低 3%，超越 sifting
- **con1_216 tightlb = 34**：比 sifting(35) 低 3%，超越 sifting

### 9.4 分析

1. **LB 公式修正效果最大**：将 $/2 \to /\text{NEDGE}$ 消除了错误剪枝，lb 从全面劣化变为多数电路追平 sifting
2. **增量 IG + 构建期派发**：iglb 在 cm85a 上超越 sifting，说明 IG 引力方向在 LB 剪枝场景下确实能引导找到更优位置
3. **代价**：$/\text{NEDGE}$ 使下界更松 → 剪枝更少 → 部分电路（add6）探索更多但不一定找到更优解。这是 LB 松紧的固有 trade-off
