# QST 理论文档：概念、定义、用途与算法伪代码

## 1. 量子态层析（QST）概述

量子态层析（Quantum State Tomography, QST）是量子信息领域中用于**完整重建量子系统状态**的实验方法。其核心思想是：通过对大量相同量子态副本在不同测量基下进行测量，利用统计分析从测量数据中推断出系统的完整量子态描述——密度矩阵 $\rho$。

---

## 2. 核心数学定义

### 2.1 密度矩阵

对于 $N$ 量子比特系统，量子态由**密度矩阵** $\rho$ 描述：

$$\rho \in \mathbb{C}^{2^N \times 2^N}$$

满足三个物理约束：

$$\text{Tr}(\rho) = 1 \quad \text{（归一性）}$$

$$\rho = \rho^\dagger \quad \text{（厄米性）}$$

$$\rho \succeq 0 \quad \text{（半正定性，即所有特征值} \geq 0\text{）}$$

对于纯态 $|\psi\rangle$，密度矩阵退化为秩一矩阵：

$$\rho = |\psi\rangle\langle\psi|$$

密度矩阵的自由度为 $4^N - 1$（考虑厄米性后实自由度减半，再减去迹归一约束）。

### 2.2 Pauli 算符基

$N$ 量子比特系统的测量基由 **Pauli 算符的张量积**构成。单量子比特的三个 Pauli 矩阵为：

$$\sigma_X = \begin{pmatrix}0 & 1 \\ 1 & 0\end{pmatrix}, \quad
\sigma_Y = \begin{pmatrix}0 & -i \\ i & 0\end{pmatrix}, \quad
\sigma_Z = \begin{pmatrix}1 & 0 \\ 0 & -1\end{pmatrix}$$

$N$ 量子比特的 Pauli 基元素为：

$$P_{b_1 b_2 \cdots b_N} = \sigma_{b_1} \otimes \sigma_{b_2} \otimes \cdots \otimes \sigma_{b_N}, \quad b_i \in \{X, Y, Z\}$$

共 $3^N$ 个独立测量基。完整层析需要枚举全部 $3^N$ 个基。

### 2.3 投影测量与 Born 定理

对于 Pauli 基 $b = (b_1, \ldots, b_N)$ 和比特串结果 $s = (s_1, \ldots, s_N) \in \{0,1\}^N$，对应的**投影算子**为：

$$\Pi_{b,s} = \bigotimes_{i=1}^{N} \Pi_{b_i, s_i}^{(1)}$$

其中单量子比特投影算子为：

$$\Pi_{Z,0} = |0\rangle\langle 0| = \begin{pmatrix}1&0\\0&0\end{pmatrix}, \quad
\Pi_{Z,1} = |1\rangle\langle 1| = \begin{pmatrix}0&0\\0&1\end{pmatrix}$$

$$\Pi_{X,0} = |{+}\rangle\langle{+}| = \frac{1}{2}\begin{pmatrix}1&1\\1&1\end{pmatrix}, \quad
\Pi_{X,1} = |{-}\rangle\langle{-}| = \frac{1}{2}\begin{pmatrix}1&-1\\-1&1\end{pmatrix}$$

$$\Pi_{Y,0} = \frac{1}{2}\begin{pmatrix}1&-i\\i&1\end{pmatrix}, \quad
\Pi_{Y,1} = \frac{1}{2}\begin{pmatrix}1&i\\-i&1\end{pmatrix}$$

满足 $\Pi_{b,0} + \Pi_{b,1} = I$（完备性）以及 $\Pi_{b,s}^2 = \Pi_{b,s}$（幂等性）。

**Born 定理**给出测量结果 $(b,s)$ 的概率：

$$p(b,s) = \text{Tr}\!\left(\Pi_{b,s}\,\rho\right)$$

### 2.4 QST 的逆问题

QST 的核心是求解以下**线性逆问题**：

给定 $K$ 个测量 $\{(\Pi_k, f_k)\}_{k=1}^K$，其中 $f_k$ 是观测频率（由实验或模拟计算），求密度矩阵 $\rho$ 使得：

$$f_k \approx p_k(\rho) = \text{Tr}(\Pi_k\,\rho), \quad k = 1, \ldots, K$$

约束条件：$\rho \succeq 0$，$\text{Tr}(\rho) = 1$，$\rho = \rho^\dagger$。

当 $K \geq 4^N - 1$ 且测量集合完备时，方程有唯一解（完整层析）；当 $K < 4^N - 1$ 时系统欠定（局部层析）。

---

## 3. 最大似然估计（MLE）

### 3.1 目标函数

最大似然估计在所有物理可行密度矩阵中，最大化对数似然函数：

$$\hat{\rho} = \arg\max_{\rho \succeq 0,\, \text{Tr}(\rho)=1} \mathcal{L}(\rho),
\quad \mathcal{L}(\rho) = \sum_{k=1}^{K} f_k \log p_k(\rho)$$

这等价于最小化观测分布与理论分布之间的 KL 散度：

$$D_{\text{KL}}(f \| p(\rho)) = \sum_k f_k \log\frac{f_k}{p_k(\rho)}$$

### 3.2 $R\rho R$ 迭代算法

由于密度矩阵的半正定约束，直接梯度下降会破坏物理性质。本项目采用 **$R\rho R$ 迭代**（Hradil 1997），每步自动保持半正定性：

**迭代公式：**

$$R^{(t)} = \sum_{k=1}^{K} \frac{f_k}{p_k^{(t)}} \Pi_k, \qquad p_k^{(t)} = \text{Tr}\!\left(\Pi_k\,\rho^{(t)}\right)$$

$$\rho^{(t+1)} = \frac{R^{(t)}\,\rho^{(t)}\,R^{(t)}}{\text{Tr}\!\left(R^{(t)}\,\rho^{(t)}\,R^{(t)}\right)}$$

**收敛性质：**

- 若 $\rho^{(0)} \succ 0$，则每步 $\rho^{(t)} \succ 0$（严格半正定不变量）
- $\mathcal{L}(\rho^{(t)})$ 单调不减，收敛到局部极大值
- 当测量完备时收敛到全局最优解

**停止条件：**$\|R^{(t)} - I\|_F < \varepsilon$ 或达到最大迭代次数

### 3.3 计算复杂度

每轮迭代的主要开销：

| 操作 | 次数 | 单次代价 |
|------|------|---------|
| $\text{Tr}(\Pi_k\,\rho)$ | $K$ 次 | DD 矩阵乘法 + 迹 |
| $R = \sum_k w_k \Pi_k$ | $K$ 次加权累加 | DD 加法 |
| $R\rho R$ | 2 次 | DD 矩阵乘法 |
| 归一化 | 1 次 | 迹 + 标量乘法 |

总时间复杂度：$O(T \cdot K \cdot C_{\text{DD}})$，其中 $T$ 为迭代次数，$K = \text{ValidBases}$，$C_{\text{DD}}$ 为单次 DD 矩阵乘法代价。

---

## 4. 保真度与迹距离

### 4.1 保真度

保真度衡量重建态与真实态的相似程度。对纯态 $|\psi\rangle$ 和密度矩阵 $\rho_{\text{recon}}$：

$$F = \langle\psi|\,\rho_{\text{recon}}\,|\psi\rangle = \text{Tr}(|\psi\rangle\langle\psi|\,\rho_{\text{recon}})$$

$$F \in [0, 1], \quad F = 1 \Leftrightarrow \rho_{\text{recon}} = |\psi\rangle\langle\psi|$$

### 4.2 迹距离

本实现中迹距离采用近似公式：

$$D \approx \sqrt{1 - F}$$

（严格定义为 $D = \frac{1}{2}\text{Tr}|\rho_1 - \rho_2|$，此近似在纯态情形下对 Uhlmann 保真度精确成立。）

---

## 5. 项目算法伪代码

### 5.1 主流程

```
算法：QST-QMDD（量子态层析，QMDD 实现）
输入：量子电路 C，测量基数 nBases，迭代次数 T，策略 strat，
       同步重排间隔 syncInterval，同步重排阈值 syncThreshold
输出：CSV 结果（保真度、RhoSize、PeakDD、时间）

1.  构建电路功能矩阵：U ← buildFunctionality(C)
2.  计算真实输出态：|ψ⟩ ← U |0⟩^N
3.  生成测量基列表：
      若 nBases ≥ 3^N：枚举全部 3^N 个 Pauli 基（完整层析）
      否则：随机抽取 nBases 个 Pauli 基（局部层析）
4.  构建测量集合：
      对每个基 b，枚举所有 2^N 个比特串 s：
        计算概率 f_{b,s} ← ⟨ψ|Π_{b,s}|ψ⟩
        若 f_{b,s} > 10^{-10}：将 (Π_{b,s}, f_{b,s}) 加入 M
      归一化：f_k ← f_k / Σ_k f_k
5.  对每种 Sifting 策略 σ：
      rho ← MLE-RhoR(M, T, σ, syncInterval, syncThreshold)
      F ← ⟨ψ|rho|ψ⟩
      记录（F, RhoSize, PeakDD, 时间）到 CSV
6.  输出 CSV
```

### 5.2 MLE $R\rho R$ 迭代（含同步重排）

```
算法：MLE-RhoR
输入：测量集合 M = {(Π_k, f_k)}，迭代次数 T，Sifting 策略 σ，
       syncInterval（轮次触发），syncThreshold（节点数触发）
输出：重建密度矩阵 rho

1.  rho ← I / 2^N                   // 初始化为最大混合态（DD 表示）
2.  for t = 0, 1, ..., T-1 do:
3.    // ── 同步重排触发判断 ──
4.    byInterval  ← (syncInterval  > 0) AND (t > 0) AND (t mod syncInterval = 0)
5.    byThreshold ← (syncThreshold > 0) AND (|rho|_DD > syncThreshold)
6.    if (byInterval OR byThreshold) AND σ ≠ None:
7.      rho ← Sifting(rho, σ)       // 更新 Unique Table，projs 自动同步
8.    end if
9.    // ── 构建 R 算子 ──
10.   R ← 0
11.   for each (Π_k, f_k) in M do:
12.     p_k ← Tr(Π_k × rho)         // DD 矩阵乘法 + trace
13.     p_k ← max(p_k, 10^{-12})    // 数值稳定性下界
14.     w_k ← min(f_k / p_k, 10^6)  // 权重（截断防溢出）
15.     R ← R + w_k · Π_k           // DD 加权累加
16.   end for
17.   if R = 0: break               // 退化保护
18.   // ── 更新 rho ──
19.   rho' ← R × rho × R            // DD 矩阵乘法（两次）
20.   τ ← Tr(rho')
21.   if τ ≤ 0 OR τ > 10^{20}: break  // 数值发散保护
22.   rho ← rho' / τ                // 归一化
23.   记录 peakDD ← max(peakDD, pkg.maxActive)
24. end for
25. // ── MLE 完成后做一次 Sifting（用于测量 DD 压缩效果）──
26. if σ ≠ None: rho ← Sifting(rho, σ)
27. return rho
```

### 5.3 投影算子构建

```
算法：BuildProjector
输入：N 量子比特，Pauli 基向量 b = (b_1,...,b_N)，比特串 s = (s_1,...,s_N)
输出：N 变量矩阵 DD 表示的 Π_{b,s}

1.  P ← I_N                    // N 量子比特单位矩阵 DD
2.  for q = 0, 1, ..., N-1 do:
3.    根据 (b_q, s_q) 选择单量子比特投影矩阵 M_q：
        (Z,0): [[1,0],[0,0]]
        (Z,1): [[0,0],[0,1]]
        (X,0): [[1,1],[1,1]]/2
        (X,1): [[1,-1],[-1,1]]/2
        (Y,0): [[1,-i],[i,1]]/2
        (Y,1): [[1,i],[-i,1]]/2
4.    G_q ← makeGateDD(M_q, N, target=q)  // 嵌入 N 量子比特空间
5.    P ← P × G_q                          // DD 矩阵乘法
6.  end for
7.  return P
```

### 5.4 同步重排机制

```
算法：SynchronousReorder
输入：rho（DD），策略 σ，变量映射 vm，交互图 ig
输出：重排后的 rho（projs 等所有活跃 DD 自动同步）

注：DD package 的 Unique Table 是 package 级别共享的。
    exchange2(i, j, vm, rho) 对所有引用计数>0 的 DD 节点做变量交换，
    因此只需对 rho 调用 Sifting，projs 自动被同步重排。

1.  prev_size ← |rho|_DD
2.  for p = 1, 2, 3 do:           // 最多 3 遍以趋近收敛
3.    rho ← dynamicReorder(rho, vm, σ)
      // dynamicReorder 内部通过 exchange2 调用 exchangeBaseCase，
      // exchangeBaseCase 遍历 Unique[i]（变量 i 的唯一表），
      // 对所有 ref > 0 的节点做 i 与 i+1 的交换
      // → projs 中所有 Π_k 自动被同步更新
4.    new_size ← |rho|_DD
5.    if new_size = prev_size: break  // 已收敛
6.    prev_size ← new_size
7.  end for
8.  return rho
```

---

## 6. QST 的用途

### 6.1 量子算法验证

通过比较实验输出态 $\rho_{\text{exp}}$ 与理论预期态 $|\psi_{\text{ideal}}\rangle$ 的保真度，验证量子算法是否正确执行：

$$F = \langle\psi_{\text{ideal}}|\,\rho_{\text{exp}}\,|\psi_{\text{ideal}}\rangle \approx 1 \Rightarrow \text{算法正确}$$

典型用例：Bell 态制备验证、Grover 搜索算法、量子隐形传态。

### 6.2 量子硬件校准

QST 可以暴露量子门的误差模式。差的保真度意味着：
- 量子门操作误差（unitary error）
- 量子比特退相干（decoherence）
- 测量误差（readout error）

### 6.3 QMDD 上层实验基准

本项目的具体目标：**以 QST 作为 QMDD 的上层应用**，评估不同 Sifting 变量排序策略对密度矩阵 DD 表示的影响：

- 测量 $\rho_{\text{recon}}$ 的 DD 节点数（RhoSize）→ **空间效益**
- 测量 MLE 迭代时间（TimeMs）→ **时间效益**  
- 测量保真度（Fidelity）→ **数值精度影响**

---

## 7. 未来工作方向

### 7.1 同步重排优化（已部分实现）

**现状**：MLE 完成后做一次 Sifting（只影响最终 RhoSize）；可选参数 `--sync-interval` / `--sync-threshold` 在迭代中途触发。

**待探索**：
- 最优触发策略：何时做同步重排收益最大？（节点膨胀率 vs. 重排开销）
- 是否每次迭代都做（完全在线重排）能获得更好的数值稳定性？
- 不同 Sifting 算法（IGSift、LBSift 等）在同步重排场景下的比较

### 7.2 压缩感知 QST（Compressed Sensing）

**动机**：对于稀疏量子态（大多数矩阵元素接近 0），可以用远少于 $4^N - 1$ 个测量恢复密度矩阵。

**方向**：
- 将 DD 表示的稀疏性先验与 $L_1$ 正则化 MLE 结合
- 对于 QMDD 节点数远小于 $2^{2N}$ 的态，理论上可以大幅减少所需测量数量

### 7.3 大规模局部层析的 Sifting 效益分析

**目标**：系统测试 $N = 5 \sim 10$ 的电路，比较：
- 无 Sifting（NoReorder）：内存大、速度慢
- 带 Sifting：内存小、但每次重排有 overhead
- 找到 Sifting 开始值得的"临界点"（RhoSize 多大时 Sifting 有净收益）

### 7.4 基于 DD 结构的自适应测量基选取

**动机**：随机选取的 Pauli 基对某些量子态效率很低（ValidBases 比例低）。

**方向**：利用 QMDD 结构信息预判哪些 Pauli 基最"信息量大"，优先选取这些基，减少所需测量数量。

### 7.5 噪声模型集成

**动机**：真实量子硬件输出的是混合态（噪声态），而非纯态。

**方向**：
- 在 QST 流程中加入信道噪声模型（去极化噪声、振幅阻尼等）
- 重建混合态密度矩阵并计算与理想纯态的距离
- 用 $F = \text{Tr}(\rho_{\text{ideal}} \rho_{\text{noisy}})$ 评估硬件保真度

### 7.6 Mode B 的实用化

**现状**：Mode B（2N 变量密度矩阵）在 $N \geq 6$ 时因变量数翻倍而不可行。

**方向**：
- 研究 Mode B 的 DD 是否有特殊稀疏结构可利用
- 探索是否可以通过部分重排降低 Mode B 的内存峰值
- Mode B 在混合态验证场景中的理论优势是否可以通过优化实现

---

## 8. 噪声信道集成

### 8.1 量子噪声信道

实际量子硬件中的噪声通过**量子信道**（Quantum Channel）描述，数学上表示为一组 **Kraus 算子** $\{K_k\}$ 对密度矩阵的完全正迹保持映射（CPTP map）：

$$\mathcal{E}(\rho) = \sum_k K_k \rho K_k^\dagger, \quad \sum_k K_k^\dagger K_k = I$$

**静态态假设**：QST 中噪声信道只作用一次，生成固定的含噪混合态 $\rho_{\text{noisy}} = \mathcal{E}(|\psi\rangle\langle\psi|)$，MLE 迭代不改变噪声。

### 8.2 三种实现的噪声信道

**去极化信道（Depolarizing, D）**：

$$\mathcal{E}_D(\rho) = (1-p)\rho + \frac{p}{3}(X\rho X + Y\rho Y + Z\rho Z)$$

Kraus 算子：$\sqrt{1-p}\,I$，$\sqrt{p/3}\,X$，$\sqrt{p/3}\,Y$，$\sqrt{p/3}\,Z$

**振幅阻尼信道（Amplitude Damping, A）**：

$$K_0 = \begin{pmatrix}1&0\\0&\sqrt{1-\gamma}\end{pmatrix}, \quad K_1 = \begin{pmatrix}0&\sqrt{\gamma}\\0&0\end{pmatrix}$$

模拟能量耗散（$|1\rangle \to |0\rangle$ 的自发辐射），$\gamma$ 为跃迁概率。

**相位翻转信道（Phase Flip, P）**：

$$\mathcal{E}_P(\rho) = (1-p)\rho + p\,Z\rho Z$$

Kraus 算子：$\sqrt{1-p}\,I$，$\sqrt{p}\,Z$

### 8.3 含噪概率计算

利用迹的线性性，避免显式构建密度矩阵：

$$p(b,s) = \mathrm{Tr}\!\left(\Pi_{b,s}\,\mathcal{E}(|\psi\rangle\langle\psi|)\right) = \sum_{\mathbf{k}} \langle\psi|\,K_{\mathbf{k}}^\dagger\,\Pi_{b,s}\,K_{\mathbf{k}}\,|\psi\rangle$$

其中 $\mathbf{k} = (k_1, \ldots, k_N)$ 遍历所有 qubit Kraus 分支的组合。

### 8.4 与理想态的距离

| 度量 | 公式 | 说明 |
|------|------|------|
| 保真度 | $F = \langle\psi\|\rho_{\text{noisy}}\|\psi\rangle$ | 最常用，本项目实现 |
| 迹距离 | $D \approx \sqrt{1-F}$ | 本项目近似实现 |
| 冯·诺依曼熵 | $S(\rho) = -\mathrm{Tr}(\rho\log\rho)$ | 衡量混合程度，纯态时 $S=0$ |

---

*文档生成时间：2026-06-17*  
*项目路径：`~/workshop/initLTQMDD`*  
*核心文件：`src/algorithms/QST.cpp`，`apps/qst_app.cpp`，`test/run_qst.sh`*
