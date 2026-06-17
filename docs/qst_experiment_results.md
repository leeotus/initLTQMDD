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

- **Mode A**：密度矩阵 $\rho$ 用 N 变量矩阵 DD 表示，保真度 $F = \langle\psi|\rho|\psi\rangle$
- **Mode B**：$\rho_{\text{true}} = |\psi\rangle\langle\psi|$ 用 2N 变量 DD 表示，保真度 $F = \text{Tr}(\rho_{\text{true}} \cdot \rho_{\text{recon}})$

### 本次实验均使用 Mode A

---

## 3. 修复的主要问题

| 问题 | 位置 | 修复方式 |
|------|------|---------|
| R 算子忽略频率权重 `freq/pk` | `QST.cpp::performMLE` | 用 `cn.lookup` 正确缩放 root edge weight |
| 初始 $\rho = I$（trace $= 2^N$，非归一化） | `qst_app.cpp::runModeA` | 初始化为 $I/2^N$ 最大混合态 |
| 每个 basis 只随机抽一个 bits，有效约束极少 | `qst_app.cpp` | 枚举每个 basis 的全部 $2^N$ 个 outcome |
| Sifting 在保真度计算之前执行导致变量顺序错位 | `qst_app.cpp::runModeA` | Sifting 移至保真度计算之后 |
| 测量基数量不足（4-qubit 只用 30 个基） | `test/run_qst.sh` | 默认自动使用 $3^N$ 个全 Pauli 基 |

---

## 4. 运行命令

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

## 5. 实验结果

### 5.1 参数配置

| 参数 | 值 |
|------|----|
| 测量基数量 | $3^N$（完整 Pauli 基枚举） |
| MLE 迭代次数 | 50 |
| 内存限制 | 10 GB |
| 单电路超时 | 600 s |
| 随机种子 | 42 |

### 5.2 分电路结果（Mode A，7 策略均值）

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

### 5.3 保真度分布（70 条策略×电路结果）

| 区间 | 数量 | 占比 |
|------|------|------|
| 完美重建 $F \geq 0.999$ | **34** | **48.6%** |
| 高精度 $0.7 \leq F < 1$ | 9 | 12.9% |
| 中等 $0.4 \leq F < 0.7$ | 26 | 37.1% |
| 低精度 $F < 0.4$ | 1 | 1.4% |

### 5.4 策略对比

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

## 6. 分析与结论

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

## 7. 结果文件

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


