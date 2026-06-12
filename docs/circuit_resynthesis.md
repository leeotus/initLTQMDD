# 量子电路重综合 (Circuit Resynthesis via DD + Sifting + Givens Decomposition)

## 1. 概述

本工具实现了一种基于决策图 (QMDD) 的量子电路重综合分析流程：

1. 从输入电路构建 QMDD
2. 应用 IGGroup Sifting 等策略压缩 DD
3. 从压缩后的 DD 提取完整酉矩阵
4. 通过 Givens 旋转分解分析重综合所需的门数

核心价值：**紧凑的 DD 表示（通过 sifting 获得）可减少 Givens 旋转数量，从而指导产生更精简的电路分解**。

## 2. 原理

### 2.1 从 DD 到酉矩阵

QMDD 是量子酉矩阵的紧凑表示。一个 n-qubit 电路对应 2^n × 2^n 的酉矩阵 U。QMDD 通过共享子图实现指数级压缩。

矩阵元素提取：对于 U[row][col]，从 DD 根节点到终端节点的路径编码了行/列的二进制位。在每一层 (qubit) q：

```
edge_index = row_bit(q) * 2 + col_bit(q)
```

沿路径累乘边权重即得矩阵元素值。

### 2.2 Sifting 压缩

变量重排序 (sifting) 改变 DD 中 qubit 的层序，使得子图共享最大化。不同策略的效果：

| 策略 | 特点 |
|------|------|
| none | 不重排序，保持原始层序 |
| sifting | 标准 sifting，逐个变量尝试所有位置 |
| iggroup | IG 引导 + 对称组合探测，通常压缩率最优 |

**关键洞察**：DD 节点数与矩阵的结构复杂度正相关。压缩后的 DD 意味着矩阵具有更多可利用的结构冗余。

### 2.3 Givens 旋转分解

任何酉矩阵可通过一系列 Givens 旋转分解为对角矩阵：

```
G_N · ... · G_2 · G_1 · U = D (对角)
```

因此：

```
U = G_1† · G_2† · ... · G_N† · D
```

每个 G_k 是一个 "two-level unitary"——仅在两个基态 |i⟩, |j⟩ 之间非平凡的酉变换。

### 2.4 Two-Level Unitary → 量子门

当 |i⟩ 和 |j⟩ 仅在 1 bit 上不同时，该 two-level unitary 直接对应一个**多控制单 qubit 门**：

- 不同的 bit → 目标 qubit
- 其余 bit 的值 → 控制条件

控制数量决定了所需 CNOT 门数：

| 控制数 | 预估 CNOT 数 | 总门数 |
|--------|-------------|--------|
| 0      | 0           | 1      |
| 1      | 2           | 5      |
| 2      | 8           | 17     |
| k (k≥3)| ~4k         | ~8k    |

## 3. 算法伪代码

### 3.1 整体流程

```
输入: 量子电路 C, sifting 策略 S
输出: 重综合分析报告

1. dd = BuildDD(C)                    // 构建 QMDD
2. dd = ApplySifting(dd, S)           // 压缩 DD
3. U = ExtractMatrix(dd, n_qubits)    // 提取 2^n × 2^n 酉矩阵
4. stats = GivensAnalysis(U, n_qubits) // Givens 分解分析
5. OutputReport(stats)
```

### 3.2 矩阵提取

```
function ExtractMatrix(dd, n):
    dim = 2^n
    for row in 0..dim-1:
        for col in 0..dim-1:
            path = ""
            for q in 0..n-1:
                row_bit = (row >> q) & 1
                col_bit = (col >> q) & 1
                path[q] = row_bit * 2 + col_bit
            U[row][col] = dd.getValueByPath(path)
    return U
```

### 3.3 Givens 分解 (1-bit-adjacent 约束)

```
function GivensDecompose(U, n):
    dim = 2^n
    work = copy(U)
    ops = []

    for col in 0..dim-2:
        for row in dim-1 downto col+1:
            if |work[row][col]| < ε: continue

            // 找从 row 到 col 的 1-bit 翻转路径
            diff = row XOR col
            path = [row]
            cur = row
            for each bit b set in diff:
                cur = cur XOR (1 << b)
                path.append(cur)

            // 沿路径级联 Givens 消元
            for step in 0..len(path)-2:
                r = path[step]      // 待消除行
                p = path[step+1]    // 伙伴行 (与 r 相差 1 bit)
                
                x = work[p][col]
                y = work[r][col]
                mag = sqrt(|x|² + |y|²)
                
                cos = conj(x) / mag
                sin = conj(y) / mag
                
                // 应用 Givens: 消除 work[r][col]
                for j in 0..dim-1:
                    (work[p][j], work[r][j]) = 
                        (cos†·work[p][j] + sin†·work[r][j],
                         -sin·work[p][j] + cos·work[r][j])
                
                ops.append(TwoLevelOp(p, r, G†))
    
    // work 现为对角矩阵
    return ops, diag(work)
```

### 3.4 门数估算

```
function EstimateGateCount(ops, n):
    cx_count = 0
    for each op in ops:
        num_ctrls = n - 1  // 除目标 bit 外所有 bit 为控制
        cx_count += CX_cost(num_ctrls)
    return cx_count
```

## 4. 使用方法

### 4.1 编译

```bash
cd release
cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_POLICY_VERSION_MINIMUM=3.5
make -j$(nproc) circuit_resynth
```

### 4.2 运行

```bash
# 基本用法
./circuit_resynth <input_circuit> [strategy]

# 示例
./circuit_resynth circuits/rd53_130.real iggroup
./circuit_resynth circuits/hwb8_114.real sifting
./circuit_resynth my_circuit.qasm none
```

### 4.3 输出解读

```
=== 量子电路重综合分析 ===
电路: circuits/rd53_130.real (7 qubits, 32 gates)
策略: iggroup
---
DD nodes:       103 -> 67 (35.0% 压缩)      ← sifting 压缩效果
Givens rotations: 202                        ← 所需 Givens 旋转数
  0-ctrl: 0, 1-ctrl: 0, 2-ctrl: 0, 3+ctrl: 202
预估 CX:        4860                         ← 预估 CNOT 门数
预估总门数:     9732                         ← 预估总门数
原始门数:       32                           ← 原始电路门数
耗时:           0.084s
```

### 4.4 支持的输入格式

- `.real` — 可逆电路格式
- `.qasm` — OpenQASM 2.0
- `.tfc` — TFC 格式

### 4.5 支持的 Sifting 策略

| 参数 | 策略 |
|------|------|
| `none` | 不做重排序 |
| `sifting` | 标准 Sifting |
| `lb` | Lower Bound Sifting |
| `ig` | Interaction Graph Sifting |
| `iglb` | IG + LB Sifting |
| `group` | Group Sifting |
| `iggroup` | IG Group Sifting (推荐) |

## 5. 限制

1. **规模限制**: 最多支持 10 qubits。矩阵提取需要 O(4^n) 空间和时间。
   - 8 qubits: 256×256 矩阵, 1MB, ~0.1s
   - 10 qubits: 1024×1024 矩阵, 16MB, ~10s
   - 12+ qubits: 不可行

2. **门数膨胀**: Givens 分解是通用方法，不利用电路的特定结构（如可逆性、稀疏性）。对于已经高度优化的小电路，重综合后门数会显著增加。

3. **多控制门分解**: 当前仅提供门数估算。精确的 C^k U 门级分解（尤其 k≥2）实现复杂，涉及辅助比特和递归分解。

4. **全局相位**: 对角矩阵中的全局相位被忽略，可能导致重综合电路与原电路差一个全局相位（物理上等价）。

5. **变量序与量子比特映射**: Sifting 改变了 DD 的变量顺序，重综合后的电路隐含了不同的 qubit 排列。需要在输出时添加 SWAP 门恢复原始映射。

## 6. 与 DD 压缩的关系

DD 节点数与 Givens 旋转数的关系：

- DD 节点少 → 矩阵结构冗余多 → 零元素更多或重复模式更多
- 零元素跳过 → Givens 旋转数减少
- Givens 少 → 重综合门数少

实验验证 (rd53_130, 7 qubits):

| 策略 | DD 节点 | Givens 数 | 预估 CX |
|------|---------|-----------|---------|
| none | 103 | 214 | 5150 |
| iggroup | 67 | 202 | 4860 |
| 减少 | 35% | 5.6% | 5.6% |

**结论**: IGGroup sifting 压缩 DD 35%, 直接导致 Givens 旋转减少 ~6%, CNOT 门数等比例减少。

## 7. 优化方向

### 7.1 短期可实现

1. **精确多控制门分解**: 实现 Barenco et al. (1995) 的递归 C^n U 分解，使用辅助 qubit 降低 CNOT 复杂度：
   - C²U: 从 ~8 CX 降至 6 CX (标准 Toffoli 分解)
   - C^k U: 从 O(4k) 降至 O(k²) 使用线性辅助

2. **对角矩阵优化**: 对角相位矩阵有专门的高效分解方法 (Welch et al. 2014)，可将 2^n 个相位合并为 O(2^n) 个 CNOT。

3. **输出 QASM 文件**: 将 Givens 分解结果直接输出为可执行的 QASM 门序列（≤4 qubits 已验证正确）。

### 7.2 中期研究

4. **利用 DD 结构直接分解**: 不提取完整矩阵，直接从 DD 的图结构递归生成门序列。DD 每层对应一个 qubit 的 Shannon 展开，可映射为 multiplexed rotation：
   ```
   DD 层 k → UCR_y(angles) + UCR_z(angles) on qubit k
   ```
   这种方法可处理 >10 qubits（不需要完整矩阵）。

5. **近似重综合**: 允许一定误差 ε，通过截断 DD 中权重小于阈值的边来减少 Givens 数量。对 NISQ 设备场景（本身有噪声），近似重综合有实际价值。

6. **结合 ZX-calculus**: 将 Givens 分解结果转化为 ZX-diagram，利用 ZX-calculus 规则进一步化简。

### 7.3 长期目标

7. **可扩展架构**: 对于大规模电路 (>10 qubits)：
   - 分块策略: 将电路切分为子电路（每块 ≤10 qubits），分别重综合后拼接
   - 基于 DD 遍历的直接分解: 避免矩阵提取，复杂度从 O(4^n) 降至 O(DD_size × n)

8. **Sifting 指导的分解顺序**: 利用 IG (Interaction Graph) 信息决定 Givens 消元顺序——优先消除交互最强的 qubit 对，可减少中间态膨胀。

9. **与量子编译器集成**: 将本工具作为编译 pass 集成到 Qiskit/tket 等框架中：
   - 前端: 接收优化后电路
   - DD 分析: 评估是否值得重综合（DD 压缩率超过阈值才触发）
   - 后端: 输出平台适配的门集 (CX+U3 或 native gate set)
