# QEC Benchmark 使用指南

## 概述

`qec_benchmark` 是一个基于 QMDD (Quantum Multiple-valued Decision Diagram) 的**量子纠错码模拟与评估工具**。它利用 DD 的图结构压缩 + IG (Interaction Graph) 对称性检测 + ddsim 噪声引擎，高效模拟 Steane [[7,1,3]] 等 stabilizer 码在噪声下的行为。

## 编译

### 方式一：集成到根项目编译

```bash
cd /home/flareon/workshop/lbqmdd

# 配置：开启 QEC 模块
cmake -S . -B build_qec \
  -DBUILD_QEC=ON \
  -DBUILD_DDSIM=OFF \
  -DBUILD_QFR_TESTS=OFF \
  -DGIT_SUBMODULE=OFF \
  -DCMAKE_BUILD_TYPE=Release

# 编译
cmake --build build_qec --target qec_benchmark -j$(nproc)
```

### 方式二：直接编译（如果已有 debug 构建目录）

```bash
cd build_qec
cmake .. -DBUILD_QEC=ON
make qec_benchmark -j$(nproc)
```

### 编译产物

```
build_qec/qec/qec_benchmark   # 可执行文件
build_qec/qec/libqec.a        # QEC 静态库
```

## 运行

```bash
# 直接运行所有 5 个实验
./build_qec/qec/qec_benchmark
```

### 实验列表

| 编号 | 实验名称 | 内容 |
|------|---------|------|
| 1 | **编码正确性验证** | 构建 Steane 编码电路，1000 shots 验证 |0⟩_L 解码为 logical 0 |
| 2 | **DD 压缩效果** | 5 轮 QEC 电路，对比 sifting vs no-sifting 的 DD size |
| 3 | **逻辑错误率 vs 物理错误率** | 扫描 p ∈ [1e-5, 1e-1]，输出 LER ± 2σ 曲线 |
| 4 | **DD 可扩展性** | 1-20 轮 QEC 的 DD size 增长 vs State Vector 大小 |
| 5 | **IG 对称性分析** | 检测 QEC 电路的 qubit 对称组，输出 IG 权重矩阵 |

### 运行单个实验

编辑 `qec/apps/qec_benchmark.cpp`，在 `main()` 中注释掉不需要的实验：

```cpp
int main() {
    printHeader();
    experiment1_encoding_verification();   // 实验 1
    // experiment2_dd_compression();       // 跳过实验 2
    experiment3_noise_threshold();         // 实验 3
    experiment4_round_scalability();       // 实验 4
    experiment5_ig_symmetry();            // 实验 5
    return 0;
}
```

## 预期输出示例

### 实验 1：编码正确性验证

```
══════ Experiment 1: Encoding Correctness Verification ══════

Encoding circuit: 12 gates
Physical circuit: 42 gates, 13 qubits

Decoding results (1000 shots):
  Logical 0: 1000 (100.0%)
  Logical 1: 0 (0.0%)
  Status: ✓ PASSED
```

### 实验 2：DD 压缩效果

```
══════ Experiment 2: DD Compression with QEC Circuit ══════

Circuit: 350 gates, 13 qubits
         5 QEC rounds (each round = 6 stabilizer measurements)

DD Size (no sifting):  4521 nodes (peak: 5200)
DD Size (with sifting): 1234 nodes (peak: 1800)
Compression ratio:     3.66x
Reordering calls:      7
```

### 实验 3：噪声阈值扫描

```
══════ Experiment 3: Logical Error Rate vs Physical Error Rate ══════

QEC rounds: 5 | Physical qubits: 13 | Shots per point: 500

p_physical   LER             +/- 2σ          DD size         Time(s)
--------------------------------------------------------------------------
1.00e-05     0.0000          0.0000          850             0.12
3.16e-05     0.0020          0.0040          920             0.14
1.00e-04     0.0080          0.0080          1050            0.15
...
1.00e-02     0.1240          0.0294          2100            0.22
1.00e-01     0.4860          0.0446          2800            0.28
```

### 实验 4：DD 可扩展性

```
══════ Experiment 4: DD Size Growth vs QEC Rounds ══════

Rounds     Gates      DD Size         StateVec Size   Ratio
------------------------------------------------------------------
1          70         245             8192            0.03x
2          140        452             8192            0.06x
5          350        1103            8192            0.13x
10         700        1980            8192            0.24x
20         1400       3567            8192            0.44x
```

### 实验 5：IG 对称性分析

```
══════ Experiment 5: IG Symmetry Analysis for QEC Codes ══════

Qubits: 13 (data: 7, ancilla: 6)
Symmetry groups detected: 3

Group 0 (3 qubits): {data0, data2, data4}  degree_range=[12-12]
Group 1 (4 qubits): {data1, data3, data5, data6}  degree_range=[12-12]
Group 2 (6 qubits): {anc7, anc8, anc9, anc10, anc11, anc12}  degree_range=[6-6]

IG Weight Matrix (first 7 data qubits):
     0   1   2   3   4   5   6
  0  0   6   4   4   6   4   4   degree=12
  1  6   0   4   4   4   6   4   degree=12
  ...
```

## API 使用（嵌入到自己的项目）

### 基本用法

```cpp
#include "QECCode.hpp"
#include "SteaneCode.hpp"
#include "QECSimulator.hpp"

using namespace qec;

int main() {
    // 1. 创建 Steane [[7,1,3]] 编码
    auto code = std::make_unique<SteaneCode>();

    // 2. 创建逻辑电路
    auto logical = std::make_unique<qc::QuantumComputation>(1);
    // 添加 Hadamard 门
    logical->emplace_back<qc::StandardOperation>(1, 0, qc::H);

    // 3. 编码：逻辑电路 → 物理电路（5 轮 QEC）
    auto physical = code->encodeLogicalCircuit(*logical, 5, 0);

    // 4. 用 QMDD 模拟
    auto dd = std::make_unique<dd::Package>();
    auto root = dd->makeZeroState(code->nPhysical() + code->nAncilla());
    dd->incRef(root);

    std::array<short, qc::MAX_QUBITS> line;
    line.fill(qc::LINE_DEFAULT);
    std::map<unsigned short, unsigned short> varMap;
    for (unsigned short i = 0; i < 13; ++i) varMap[i] = i;

    for (auto& op : *physical) {
        if (op->isStandardOperation()) {
            auto dd_op = op->getDD(dd, line, varMap);
            auto tmp = dd->multiply(dd_op, root);
            dd->incRef(tmp);
            dd->decRef(root);
            root = tmp;
        }
    }

    // 5. 测量并解码
    // ... 执行 measurement loop ...
    int logical_result = code->decodeMeasurement(measurement_bitstring);

    return 0;
}
```

### 生成单个子电路

```cpp
auto code = std::make_unique<SteaneCode>();

// 编码电路
auto enc = code->generateEncodingCircuit();
std::cout << enc->getNops() << " gates\n";

// 一轮 syndrome extraction
auto synd = code->generateSyndromeExtraction();
std::cout << synd->getNops() << " gates\n";

// 根据 syndrome 生成纠错电路
auto corr = code->generateCorrection({0, 1, 0, 0, 0, 0});
std::cout << corr->getNops() << " gates\n";

// 解码电路
auto dec = code->decodeLogicalState(true); // Z-basis
std::cout << dec->getNops() << " gates\n";
```

### 噪声模拟

```cpp
QECSimulator sim(std::make_unique<SteaneCode>());

// 单次模拟
auto config = QECExperimentConfig{
    .base_noise = QECSimulator::makeDepolarizingNoise(0.001),
    .shots = 1000,
    .qec_rounds = {5},
};
config.logical_circuit_path = ""; // 使用 identity

auto results = sim.runThresholdSweep(*logical, config);
for (auto& r : results) {
    std::cout << "LER = " << r.logical_error_rate
              << ", Fidelity = " << r.output_fidelity << "\n";
}
```

## 添加新的 QEC 编码

继承 `QECCode` 基类并实现所有虚函数：

```cpp
class MyCode : public QECCode {
public:
    int nPhysical() const override { return n; }
    int kLogical() const override { return k; }
    int distance() const override { return d; }
    int nAncilla() const override { return m; }

    const std::vector<StabilizerGenerator>& stabilizers() const override {
        return stabs_;
    }

    // 实现所有纯虚函数...
    std::unique_ptr<qc::QuantumComputation> generateEncodingCircuit() const override;
    std::unique_ptr<qc::QuantumComputation> generateSyndromeExtraction() const override;
    std::unique_ptr<qc::QuantumComputation> generateCorrection(
        const std::vector<int>& syndrome) const override;
    // ... etc

private:
    std::vector<StabilizerGenerator> stabs_;
};
```

## 文件结构

```
qec/
├── include/
│   ├── QECCode.hpp           # 抽象基类
│   ├── SteaneCode.hpp        # Steane [[7,1,3]]
│   └── QECSimulator.hpp      # 模拟器接口
├── src/
│   ├── QECCode.cpp           # 基类实现（encodeLogicalCircuit 等）
│   └── SteaneCode.cpp        # Steane 完整实现
├── apps/
│   └── qec_benchmark.cpp     # 基准测试主程序
├── test/
│   └── (预留测试目录)
└── CMakeLists.txt