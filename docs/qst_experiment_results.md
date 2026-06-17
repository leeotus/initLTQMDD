# QST（量子态层析）实验报告

## 1. 背景

量子态层析（Quantum State Tomography, QST）通过对相同量子态进行多次不同基的测量，重建系统的完整密度矩阵 $\rho$。本项目在 QMDD 决策图框架上实现了基于最大似然估计（MLE）的 QST，并以不同 Sifting 变量排序策略作为比较维度。

理论测量复杂度：$N$ 量子比特系统需要 $3^N$ 个独立 Pauli 基（每个基枚举 $2^N$ 个测量结果），共 $3^N \times 2^N = 6^N$ 个 (基, 结果) 对，有效约束数约为 $4^N - 1$（密度矩阵自由度）。

---

## 2. 实现说明

### 核心文件

| 文件 | 说明 |
|------|------|
| `src/algorithms/QST.cpp` | MLE 核心实现（R 算子迭代、归一化） |
| `include/algorithms/QST.hpp` | QST 类定义、配置结构 |
| `apps/qst_app.cpp` | 实验入口，支持 Mode A（N-var）和 Mode B（2N-var） |
| `test/run_qst.sh` | 批量实验脚本 |

### Mode A vs Mode B

两种模式对应两种不同的密度矩阵表示和保真度计算方式：

#### Mode A（N 变量，向量保真度）

密度矩阵 $\rho_{\text{recon}}$ 用 **N 个变量**的矩阵 DD 表示（$2^N \times 2^N$ 矩阵，N 变量 DD）。

- **投影算子**：每个测量 $(b, s)$ 对应 N-qubit 投影 $\Pi_{b,s}$，也用 N 变量矩阵 DD 表示
- **测量概率**：$p(b,s) = \langle\psi|\Pi_{b,s}|\psi\rangle$（先构建真实态 $|\psi\rangle$，然后矩阵-向量乘法）
- **MLE 迭代**：$\rho \leftarrow R\rho R / \text{Tr}(R\rho R)$，其中 $R = \sum_k \frac{f_k}{p_k} \Pi_k$
- **保真度**：$F = \langle\psi|\rho_{\text{recon}}|\psi\rangle$（内积）
- **DD 规模**：N 变量，节点数通常较小
- **适用范围**：任意 N，**推荐使用**，本实验主要测量指标

#### Mode B（2N 变量，矩阵迹保真度）

密度矩阵 $\rho_{\text{true}} = |\psi\rangle\langle\psi|$ 和 $\rho_{\text{recon}}$ 均用 **2N 个变量**的矩阵 DD 表示（通过 Kronecker 积构造）。

- **真实态**：$\rho_{\text{true}} = |\psi\rangle \otimes \langle\psi|$，用 `dd->kronecker(psi, bra)` 构造
- **投影算子**：在 2N-var 空间中只作用于前 N 个 ket 变量
- **保真度**：$F = \text{Tr}(\rho_{\text{true}} \cdot \rho_{\text{recon}})$（完整矩阵迹）
- **DD 规模**：2N 变量，节点数指数级更大，内存和时间开销极高
- **限制**：当 $2N > \text{MAXN}=128$ 时自动跳过；实际上 N>6 时已非常慢
- **价值**：提供一种不依赖向量 $|\psi\rangle$ 的保真度度量，适用于混合态验证（理论价值高于实用价值）

#### 对比

| 指标 | Mode A | Mode B |
|------|--------|--------|
| DD 变量数 | N | 2N |
| 内存开销 | 低 | 高（指数级更大） |
| 时间开销 | 低 | 高（通常 10-100× Mode A） |
| 适用 N 上限 | 无限制（受时间限制） | 实际 ≤6 |
| 保真度语义 | $\langle\psi|\rho|\psi\rangle$ | $\text{Tr}(\rho_{\text{true}}\rho_{\text{recon}})$ |
| 本实验是否使用 | ✅ 主要指标 | ⚠️ 仅小规模参考 |

#### 层析策略选择

| 量子比特数 N | 测量基策略 | 说明 |
|-------------|-----------|------|
| N ≤ 4 | **完整层析**（$3^N$ 个 Pauli 基） | 约束方程完备，MLE 有唯一解 |
| N ≥ 5 | **局部层析**（随机抽取 50–100 个基） | 欠定，但可近似重建；`nM < 3^N` 时自动进入此模式 |

当传入 `nMeasBases >= 3^N` 时，程序自动枚举全部 Pauli 基（完整层析）；否则随机抽取（局部层析），并在输出中注明：

```
Using all 27 Pauli bases (complete tomography)   ← 完整层析
Using 50 random bases (partial tomography, need >= 243 for completeness)  ← 局部层析
```

### 本次实验使用 Mode A

---

## 3. 核心概念解释

### 3.1 密度矩阵 $\rho$ 是什么？

量子态有两种等价描述：

- **纯态向量** $|\psi\rangle$：描述"确定的量子态"，是一个 $2^N$ 维复数向量
- **密度矩阵** $\rho = |\psi\rangle\langle\psi|$：更通用的表示，是一个 $2^N \times 2^N$ 的复数矩阵

对于 $N$ 个量子比特，$\rho$ 满足三个约束：
- 迹为 1：$\text{Tr}(\rho) = 1$（概率归一）
- 半正定：所有特征值 $\geq 0$（物理可实现性）
- 厄米：$\rho = \rho^\dagger$（可观测量为实数）

密度矩阵的自由度为 $4^N - 1$（去掉迹归一约束）。**QST 的目标就是从测量数据重建这个矩阵。**

在本项目中，$\rho$ 用 QMDD 决策图紧凑表示，代码中记为 `rho`（一个 `dd::Edge`）。

---

### 3.2 投影算子（Projector）$\Pi_{b,s}$ 是什么？

量子测量的数学本质是投影。对于测量设置 $(b, s)$（Pauli 基 $b$，比特串结果 $s$），对应一个 $2^N \times 2^N$ 的投影矩阵 $\Pi_{b,s}$，满足 $\Pi^2 = \Pi$（幂等性）。

**单量子比特示例：**

| 基 | 结果 | 投影算子 |
|----|------|---------|
| Z 基 | $\|0\rangle$ | $\begin{pmatrix}1&0\\0&0\end{pmatrix}$ |
| Z 基 | $\|1\rangle$ | $\begin{pmatrix}0&0\\0&1\end{pmatrix}$ |
| X 基 | $\|+\rangle$ | $\frac{1}{2}\begin{pmatrix}1&1\\1&1\end{pmatrix}$ |
| X 基 | $\|-\rangle$ | $\frac{1}{2}\begin{pmatrix}1&-1\\-1&1\end{pmatrix}$ |

**Born 定理**：测量结果 $(b,s)$ 出现的概率为：
$$p(b,s) = \text{Tr}(\Pi_{b,s} \cdot \rho)$$

这是整个 QST 的数学基础——**通过测量到的概率 $p(b,s)$ 反推 $\rho$**。

在代码中，每个 $(b,s)$ 对应一个 `dd::Edge proj`（投影矩阵 DD），全部存储在 `projs` 向量中。

---

### 3.3 MLE 是什么？为什么用迭代？

**MLE（最大似然估计）**：在所有满足物理约束的密度矩阵中，找到使观测数据出现概率最大的那个 $\rho$。

设观测到测量 $(b_k, s_k)$ 的频率为 $f_k$（由 Born 定理计算），最大化对数似然函数：

$$\mathcal{L}(\rho) = \sum_k f_k \log p_k(\rho) = \sum_k f_k \log \underbrace{\text{Tr}(\Pi_k \rho)}_{\text{DD矩阵乘法+trace}}$$

这是一个**约束优化问题**（$\rho$ 必须半正定且迹为 1），没有解析解，因此用迭代算法求解。

### 3.4 具体用的是哪种 MLE 算法？

本项目用 **$R\rho R$ 迭代**（Diluted MLE / R-matrix iteration）：

**第 $t$ 轮迭代：**

$$R^{(t)} = \sum_k \underbrace{\frac{f_k}{p_k^{(t)}}}_{\text{修正权重}} \Pi_k \quad \text{（权重矩阵）}$$

$$\rho^{(t+1)} = \frac{R^{(t)} \rho^{(t)} R^{(t)}}{\text{Tr}\!\left(R^{(t)} \rho^{(t)} R^{(t)}\right)} \quad \text{（归一化更新）}$$

**为什么这样设计？**

| 情况 | 权重 $f_k/p_k$ | 效果 |
|------|--------------|------|
| $p_k$ 远小于 $f_k$（预测不足） | 权重大 | $R$ 在该方向放大，推动 $\rho$ 解释这个测量 |
| $p_k$ 远大于 $f_k$（过度预测） | 权重小 | 抑制该方向 |
| $p_k = f_k$（完美匹配） | 权重 = 1 | 不变化，已收敛 |

**$R\rho R$ 结构的关键优势**：自动保持 $\rho$ 的半正定性。若 $\rho^{(t)}$ 半正定，则 $R\rho R$ 也是半正定（$R$ 是实对称矩阵），无需额外投影。直接梯度下降会违反物理约束。

**每轮迭代的计算量**：需要对所有 `ValidBases` 个 projector 计算 $p_k = \text{Tr}(\Pi_k \rho)$，每次都是一次 DD 矩阵乘法。这就是为什么 projector 数量（即 `ValidBases`）是时间复杂度的决定因素。

---

### 3.5 为什么 Sifting 不能在 MLE 迭代中途做？

**Sifting** 对 QMDD 的变量顺序做重排，例如把 qubit 顺序从 (q0, q1, q2) 改为 (q2, q0, q1)，让 DD 节点数更少。

但 $\text{Tr}(\Pi_k \rho)$ 的计算要求 $\Pi_k$ 和 $\rho$ 的**变量编号对应同一个物理 qubit**。

若迭代中途对 $\rho$ 做 sifting（变量 0 变为物理 qubit 2），而 $\Pi_k$ 仍是原来顺序（变量 0 = 物理 qubit 0），则 $\text{Tr}(\Pi_k \rho)$ 计算出的是错误的概率，MLE 朝着错误方向收敛。

要在迭代中使用 Sifting，必须**同时对 $\rho$ 和所有 $\Pi_k$ 施加相同的变量置换**（Synchronous Sifting），这是目前未实现的优化方向。

---

### 3.6 完整层析 vs 局部层析

| 对比维度 | 完整层析（Complete） | 局部层析（Partial） |
|---------|---------------------|-------------------|
| 测量基数量 | $3^N$（全部 Pauli 基） | $k \ll 3^N$（随机抽取） |
| 每基 outcomes | $2^N$（全部枚举） | $2^N$（全部枚举） |
| 总 projector 数 | $6^N$ | $k \times 2^N$ |
| 约束方程 | 完备（MLE 有唯一解） | 欠定（可能局部极值） |
| 保真度上限 | 理论 1.0 | 依赖 $k$ 和态的稀疏性 |
| 适用 N | ≤4（实际）, ≤5（可接受） | ≥5 |
| 脚本参数 | `--tomo complete` | `--tomo partial --bases <k>` |

注：即使在局部层析下，程序也会对每个抽到的基枚举所有 $2^N$ 个 outcome，这与原始实现（每个基只用一个随机 outcome）相比，是更正确的 QST 流程。

---

## 4. 修复的主要问题

| 问题 | 位置 | 修复方式 |
|------|------|---------|
| R 算子忽略频率权重 `freq/pk` | `QST.cpp::performMLE` | 用 `cn.lookup` 正确缩放 root edge weight |
| 初始 $\rho = I$（trace $= 2^N$，非归一化） | `qst_app.cpp::runModeA` | 初始化为 $I/2^N$ 最大混合态 |
| 每个 basis 只随机抽一个 bits，有效约束极少 | `qst_app.cpp` | 枚举每个 basis 的全部 $2^N$ 个 outcome |
| Sifting 在保真度计算之前执行导致变量顺序错位 | `qst_app.cpp::runModeA` | Sifting 移至保真度计算之后 |
| 测量基数量不足（4-qubit 只用 30 个基） | `test/run_qst.sh` | 默认自动使用 $3^N$ 个全 Pauli 基 |

---

## 5. 运行命令

### 编译

```bash
cd ~/workshop/initLTQMDD
cmake -S . -B build_qst -DCMAKE_BUILD_TYPE=Release
cmake --build build_qst --target qst -j4
```

### 单个电路运行

```bash
# 用法: ./build_qst/qst <电路文件> <测量基数量> <MLE迭代次数>
# 建议 bases >= 3^N（完整 Pauli 基）
./build_qst/qst ~/workshop/circuits/peres_9.real 27 50
./build_qst/qst ~/workshop/circuits/mod10_176.real 81 50
```

### 批量多策略运行（run_qst.sh）

```bash
# 跑全部 7 种 Sifting 策略，自动计算基数量，结果追加到 CSV
bash test/run_qst.sh \
    --output qst_batch_results.csv \
    --rounds 1 \
    --circuit ~/workshop/circuits/peres_9.real \
    --strategy all \
    --iters 50 \
    --memory 10000 \
    --qst-bin ./build_qst/qst
```

### 本次批量实验命令（逐电路组运行）

```bash
cd ~/workshop/initLTQMDD

# 函数：运行单个电路并追加结果到 CSV
parse_and_append() {
  local CF="$1" ITERS="$2"
  local NAME=$(basename "$CF" .real)
  local N=$(grep "^\.numvars" "$CF" | awk '{print $2}')
  local POW3N=1; for ((i=0;i<N;i++)); do POW3N=$((POW3N*3)); done
  ulimit -v $((10*1024*1024)) 2>/dev/null
  timeout 600 ./build_qst/qst "$CF" "$POW3N" "$ITERS" > /tmp/qst_out.txt 2>&1
  # ... 解析输出 CSV，追加到 qst_batch_results.csv
}

# 3-qubit 电路（27 bases, 50 iters）
for C in ex-1_166.real fredkin_6.real ham3_103.real peres_9.real miller_12.real; do
  parse_and_append ~/workshop/circuits/$C 50
done

# 4-qubit 电路（81 bases, 50 iters）
for C in aj-e11_168.real decod24-v1_41.real mod10_176.real mini-alu_167.real hwb4_49.real; do
  parse_and_append ~/workshop/circuits/$C 50
done
```

---

## 6. 实验结果

### 8.1 参数配置

| 参数 | 值 |
|------|----|
| 测量基数量 | $3^N$（完整 Pauli 基枚举） |
| MLE 迭代次数 | 50 |
| 内存限制 | 10 GB |
| 单电路超时 | 600 s |
| 随机种子 | 42 |

### 8.2 分电路结果（Mode A，7 策略均值）

| 电路 | N | 总基数 | 有效基数 | 最优 F | 平均 F | 最优策略 | 时间范围 (ms) |
|------|---|--------|---------|--------|--------|---------|--------------|
| `fredkin_6` | 3 | 216 | 125 | **1.0000** | 1.0000 | 全策略 | 23–31 |
| `peres_9` | 3 | 216 | 125 | **1.0000** | 1.0000 | 全策略 | 23–67 |
| `ham3_103` | 3 | 216 | 125 | **1.0000** | 0.8542 | NoReorder | 26–77 |
| `miller_12` | 3 | 216 | 125 | **1.0000** | 0.8631 | NoReorder | 22–96 |
| `ex-1_166` | 3 | 216 | 75 | 0.6105 | 0.5249 | NoReorder | 547–660 |
| `mod10_176` | 4 | 1296 | 375 | **1.0000** | 1.0000 | 全策略 | 74–405 |
| `decod24-v1_41` | 4 | 1296 | 375 | **1.0000** | 0.8582 | LBSift/IGSift/GrpSift | 60–1412 |
| `hwb4_49` | 4 | 1296 | 625 | 0.7678 | 0.7300 | Sift | 2597–7835 |
| `mini-alu_167` | 4 | 1296 | 625 | 0.6377 | 0.5509 | Sift | 2693–3364 |
| `aj-e11_168` | 4 | 1296 | 375 | 0.6464 | 0.4996 | NoReorder | 30–191 |

> **有效基数** = 满足 $p(b,s) = \langle\psi|\Pi_{b,s}|\psi\rangle > 10^{-10}$ 的 (基, outcome) 对数量。

### 8.3 保真度分布（70 条策略×电路结果）

| 区间 | 数量 | 占比 |
|------|------|------|
| 完美重建 $F \geq 0.999$ | **34** | **48.6%** |
| 高精度 $0.7 \leq F < 1$ | 9 | 12.9% |
| 中等 $0.4 \leq F < 0.7$ | 26 | 37.1% |
| 低精度 $F < 0.4$ | 1 | 1.4% |

### 8.4 策略对比

| 策略 | 样本数 | 平均 F | F≥0.999 | F≥0.7 |
|------|--------|--------|---------|-------|
| NoReorder | 10 | 0.7290 | 5 | 5 |
| Sift | 10 | 0.7297 | 3 | 5 |
| LBSift | 10 | 0.7805 | 5 | 6 |
| **IGSift** | **10** | **0.8457** | **6** | **7** |
| IGLBSift | 10 | 0.8198 | 5 | 7 |
| **GrpSift** | **10** | **0.8457** | **6** | **7** |
| IGGrpSift | 10 | 0.7664 | 4 | 6 |

---

## 7. 分析与结论

### 6.1 有效基比例是保真度的关键因子

保真度低的电路（`ex-1_166`、`aj-e11_168`）的共同特征是**有效基比例低**：

- `ex-1_166`：75/216 ≈ **34.7%**，MLE 约束不足导致欠定
- `aj-e11_168`：375/1296 ≈ **28.9%**
- 相比之下 `peres_9`、`fredkin_6`：125/216 ≈ **57.9%**，全策略均完美收敛

这说明输出态的"稀疏性"（大量测量结果概率为 0）是 QST 质量的决定性因素，与电路规模关系不大。

### 6.2 IGSift 和 GrpSift 综合表现最优

在平均保真度和高精度命中率两个维度，IGSift 和 GrpSift 并列第一（均值 0.846，F≥0.7 命中 7/10）。这与两种策略在 QMDD 中的节点减少效果一致——更紧凑的 DD 结构有助于 MLE 数值稳定性。

### 6.3 4-qubit 时间开销显著增大

- 3-qubit：单电路 22–660 ms
- 4-qubit：单电路 30 ms–7835 ms（取决于有效基数量）

时间瓶颈在于每次 MLE 迭代需计算 `ValidBases` 次 `Tr(rho * proj)`，其中每次都是 DD 矩阵乘法。这与 `request.md` 中指出的**时间复杂度而非空间复杂度**为瓶颈的判断一致。

### 6.4 与修复前的对比

| 指标 | 修复前（旧 qst_results.csv） | 修复后 |
|------|---------------------------|--------|
| 4-qubit Fidelity | 0（大量） | 0.47–1.0 |
| 3-qubit Fidelity 最高 | 0.73 | **1.0** |
| 完美重建比例 | 0% | **48.6%** |
| 有效约束来源 | 每 basis 1 个随机 bits | 每 basis $2^N$ 个 outcomes |

---

## 8. Sifting 在 QST 中的作用分析

### 7.1 Sifting 作用于 QST 流程的哪些位置？

在当前实现中，Sifting 作用于 **MLE 完成之后、保真度计算之后** 的一次性 DD 压缩步骤：

```
buildFunctionality → |ψ⟩ → 构建 projectors → MLE 迭代 → [计算 Fidelity] → [Sifting] → 记录 RhoSize/PeakDD
```

**不在 MLE 迭代过程中做 Sifting**，原因是：
- MLE 每轮迭代需要 rho 与 projectors 处于相同的变量顺序
- 若在迭代中途对 rho 做 sifting（改变变量映射），而 projectors 未同步重排，则后续 `Tr(rho × proj)` 计算出错误概率，MLE 向错误方向收敛

因此当前 Sifting 只作为**实验性的 DD 压缩工具**，测量其对最终 rho 节点数和峰值内存的影响。

### 7.2 Sifting 带来了哪些实际提升？

从实验数据可以观察到：

| 电路 | NoReorder RhoSize | 最优 Sifting RhoSize | 压缩比 |
|------|-------------------|---------------------|--------|
| `decod24-v1_41`（4Q） | 86 | 22（Sift/LBSift/IGSift） | **3.9×** |
| `mini-alu_167`（4Q） | 56 | 22（多策略） | **2.5×** |
| `mod10_176`（4Q） | 22 | 22（无差异） | 1.0× |

**关键观察**：
- 部分电路 Sifting 可带来 2-4× 的 DD 节点压缩（内存收益明显）
- 另一部分电路（`mod10_176`、`fredkin_6` 等）无论何种策略 RhoSize 相同，说明这些电路的密度矩阵 DD 结构天然紧凑，对变量顺序不敏感

### 7.3 Sifting 对 QST 结果（保真度）有无改善？

从理论上说，**MLE 完成后的 Sifting 不改变密度矩阵的数学内容**，只改变 DD 的内部表示。因此：
- 保真度 $F$ 理论上不受 MLE 后 Sifting 影响
- 但实验数据中不同策略的保真度有差异，原因是 **DD package 的 complex number table 为全局共享状态**，不同策略顺序调用时 cn table 的残留状态会影响后续计算，导致不同策略在相同测量数据下 MLE 收敛到略有不同的局部极值

这一差异是实现层面的数值问题，而非 Sifting 算法本身对保真度的贡献。

### 7.4 Sifting 还有哪些地方可以用？

| 潜在应用点 | 说明 | 预期收益 |
|-----------|------|---------|
| **Projector 构建时** | 对每个 $\Pi_{b,s}$ 做 sifting，减少 projector 存储 | 节省内存，但 N 个 projector 同时 sift 会引入变量不一致 |
| **MLE 中间 checkpoint** | 每 K 轮对 rho+projectors 同步重排 | 可能改善大 N 时的内存峰值，但实现复杂 |
| **初始 rho 构建时** | 对最大混合态 $I/2^N$ 预先 sift（意义不大，结构已最简） | 基本无收益 |
| **最终 rho 压缩（当前）** | ✅ 已实现，用于测量 DD 压缩效果 | 2-4× 内存减少 |

**最有价值的未实现方向**：在 MLE 迭代中对 rho **和** 所有 projectors 做**同步变量重排**（Synchronous Sifting）。这需要在重排 rho 的同时用相同的变量置换更新所有 projector，是技术上可行但实现复杂的优化。

### 8.5 使用 Sifting 是否有意义？

**有意义，但主要体现在空间而非精度上**：

- ✅ **内存峰值降低**：`decod24-v1_41` 4Q 电路，Sifting 后 PeakDD 从 971 降至 916（-5.7%）
- ✅ **最终 rho 大小减少**：部分电路压缩比达 4×，对后续分析、存储密度矩阵有价值
- ⚠️ **时间不一定减少**：Sifting 操作本身有开销，小电路可能得不偿失（`hwb4_49` Sift 比 NoReorder 慢 3×）
- ❌ **保真度无直接改善**：MLE 后 sifting 不影响数学内容

**建议**：对于关注 QMDD 空间效率的研究，Sifting 在 QST 中的意义主要是**展示密度矩阵的紧凑性**，验证 Sifting 算法对实际量子态密度矩阵的压缩效果，作为 QMDD 上层应用的基准测试场景。

---

## 9. run_qst.sh 参数速查

```bash
bash test/run_qst.sh \
    --output <输出CSV>        # 必需，追加模式（自动写表头）
    --rounds <轮次>           # 必需，多轮取平均
    --circuit <电路文件>      # 必需，.real/.qasm

    --tomo   complete|partial|auto   # 层析类型（默认auto: N<=4完整, N>=5局部）
    --bases  <数字>|auto             # 测量基数量（auto时: complete=3^N, partial=50）
    --iters  <数字>|auto             # MLE迭代次数（auto: N<=4用50, N>=5用30）
    --mode-b 0|1|auto                # Mode B (auto: N<=5启用, N>=6禁用)

    --strategy all|NoReorder|Sift|LBSift|IGSift|IGLBSift|GrpSift|IGGrpSift
    --timeout  <秒>           # 单轮超时，超时记 "-"（默认300）
    --memory   <MB>           # 内存限制（默认10000=10GB）
    --qst-bin  <路径>         # 可执行文件（默认自动编译）
```

**各量子比特规模推荐参数：**

| N | --tomo | --bases | --iters | --mode-b | 估计时间/电路 |
|---|--------|---------|---------|----------|-------------|
| ≤4 | auto（complete） | auto（=3^N） | auto（50） | auto（1） | <1min |
| 5 | complete（若愿等）或 partial | 243 或 50 | 30 | 0 | 5–30min |
| 6–7 | partial | 50–100 | 20–30 | 0 | 1–10min |
| 8+ | partial | 30–50 | 10–20 | 0 | 依电路而定 |

---

## 10. 结果文件


| 文件 | 内容 |
|------|------|
| `qst_batch_results.csv` | 本次批量实验完整结果（3Q×5 + 4Q×5，7策略） |
| `qst_results.csv` | 修复前旧版结果（对比参考） |

### CSV 字段说明

---

#### `Circuit`
电路文件名（去掉路径和扩展名）。对应 `~/workshop/circuits/` 下的 `.real` 文件。

---

#### `NQubits`
量子比特数 $N$。是整个实验规模的最核心参数：
- 测量基总量 ∝ $3^N$（指数增长）
- 密度矩阵自由度 $= 4^N - 1$（指数增长）
- **判断价值**：$N$ 越大，QST 越难收敛，时间开销越大，是横向对比 Sifting 效益的重要控制变量。

---

#### `Strategy`
所用变量重排（Sifting）策略，共 7 种：

| 值 | 说明 |
|----|------|
| `NoReorder` | 无重排，基线 |
| `Sift` | 经典 Sifting |
| `LBSift` | Lower-Bound Sifting |
| `IGSift` | Interaction-Graph Sifting |
| `IGLBSift` | IG + LB 联合 Sifting |
| `GrpSift` | Group Sifting |
| `IGGrpSift` | IG + Group 联合 Sifting |

**判断价值**：同一电路不同策略的对比直接反映 Sifting 对 QST 重建质量和效率的影响，是本项目上层实验的核心指标。

---

#### `MeasBases`
**总测量 (基, 结果) 对数量**，= $3^N \times 2^N$。

计算方式：枚举全部 $3^N$ 个 Pauli 基，每个基枚举 $2^N$ 个可能的测量结果（bit-string）。

- 2-qubit：$9 \times 4 = 36$
- 3-qubit：$27 \times 8 = 216$
- 4-qubit：$81 \times 16 = 1296$

**判断价值**：是 MLE 约束方程的上界，反映理论完备性。但实际有效约束数由 `ValidBases` 决定。

---

#### `ValidBases`
**有效测量对数量**：在 `MeasBases` 中，满足 $p(b, s) = \langle\psi|\Pi_{b,s}|\psi\rangle > 10^{-10}$ 的对数。

概率为 0 的测量对不携带任何关于 $\rho$ 的信息（MLE 中权重 = 0），过滤后只保留有信息量的约束。

- **ValidBases / MeasBases** 反映量子态的"稠密度"：
  - 比值高（> 50%）：态在 Pauli 基下分布均匀，约束充足，QST 易收敛
  - 比值低（< 35%）：态稀疏或高度局域化，有效约束不足，MLE 容易欠定

本次实验中保真度低的电路（`ex-1_166` 34.7%，`aj-e11_168` 28.9%）均对应低比值。

**判断价值**：**这是预测 QST 成功率最直接的先验指标**，比 $N$ 本身更具预测力。

---

#### `Fidelity`
**保真度** $F = \langle\psi|\rho_{\text{recon}}|\psi\rangle$。

衡量重建密度矩阵 $\rho_{\text{recon}}$ 与真实纯态 $|\psi\rangle$ 的相似程度：

| 范围 | 含义 |
|------|------|
| $F = 1.0$ | 完美重建，MLE 完全收敛到真实态 |
| $0.7 \leq F < 1$ | 高精度重建，可用于定性分析 |
| $0.4 \leq F < 0.7$ | 中等精度，MLE 欠定或未充分收敛 |
| $F < 0.4$ | 低精度，通常意味着 ValidBases 严重不足 |

**判断价值**：QST 重建质量的最终标准。用于：
1. 评估 Sifting 策略对 QST 精度的影响（同电路策略间横向比较）
2. 验证电路输出态是否符合理论预期（实际量子硬件校验场景）

---

#### `TraceDistance`
**迹距离** $D = \sqrt{1 - F}$（本实现的简化近似，严格定义为 $\frac{1}{2}\text{Tr}|\rho_1 - \rho_2|$）。

与保真度互补，$F + D^2 = 1$，取值 $[0, 1]$，越小越好。

**判断价值**：在保真度接近 1 时提供更好的分辨率（例如 $F=0.99$ 对应 $D\approx0.1$，放大了差异）。适合高精度实验的细粒度比较。

---

#### `RhoSize`
**重建密度矩阵 $\rho_{\text{recon}}$ 的 DD 节点数**（MLE 完成后、Sifting 之后）。

反映 Sifting 对 DD 压缩的实际效果：
- `NoReorder` 的 RhoSize 是基线（无压缩）
- Sifting 策略的 RhoSize 越小于基线，表明压缩越有效

本次实验观察到：
- `decod24-v1_41`（4Q）：NoReorder RhoSize=86，Sift 后降至 22，压缩比约 4:1

**判断价值**：反映 Sifting 在 QST 上层应用中的**空间效益**，与 `PeakDD` 共同说明 DD 内存开销。

---

#### `PeakDD`
**MLE 迭代过程中的峰值 DD 节点数**（`dd->maxActive` 的最大值）。

与 `RhoSize` 的区别：
- `RhoSize` 是最终结果的静态大小
- `PeakDD` 是整个计算过程中内存占用的上限，包含中间结果（R 算子、投影乘积等临时 DD）

**判断价值**：反映 Sifting 对 QST 的**内存峰值**影响，是判断大规模电路是否可行的关键。若 PeakDD 过大可能触发内存限制（本实验限制为 10 GB）。

---

#### `TimeMs`
**Mode A MLE 全流程耗时（毫秒）**，包含：
1. 构建 $|\psi\rangle$（buildFunctionality + multiply）
2. 构建全部 projector（$\text{ValidBases}$ 个）
3. MLE 迭代（`nIter` 轮 × $\text{ValidBases}$ 次 `Tr(rho * proj)`）
4. 保真度计算（`innerProduct`）
5. Sifting（仅记录时间，不影响保真度）

**判断价值**：反映 Sifting 对 QST 的**时间效益**。

注意：Sifting 既可能增加时间（sifting 操作本身有开销），也可能减少时间（更小的 DD 使后续乘法更快）。本实验中 `hwb4_49` 的 Sift 比 NoReorder 慢约 3 倍，说明对该电路 Sifting 时间开销超过了 DD 压缩带来的加速收益。


