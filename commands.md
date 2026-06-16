# 密度矩阵仿真实验使用指南

## 快速开始

```bash
cd /home/lijianxian/workshop/initLTQMDD
bash test/run_dm_experiments.sh
```

默认参数：n=4~12，p=0.05，10轮噪声，QFT+Clifford 合成电路。

---

## 参数说明

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--nmin N` | 最小 qubit 数 | 4 |
| `--nmax N` | 最大 qubit 数 | 12 |
| `--noise-p P` | depolarizing 噪声强度 | 0.05 |
| `--rounds K` | rounds 实验噪声轮数 | 10 |
| `--seeds "s1,s2,..."` | Clifford 随机电路种子 | "42,123" |
| `--circuits PATH` | 外部电路文件或目录 (.real/.qasm) | 无 |
| `--experiments "e1 e2"` | 运行哪些实验 | "scale rounds compress expval" |
| `--output DIR` | 结果输出目录 | results/density_matrix |
| `--build-type TYPE` | Release 或 Debug | Release |

---

## 常用命令

### 1. 仅合成电路，扩大 qubit 范围
```bash
bash test/run_dm_experiments.sh --nmin 4 --nmax 16 --noise-p 0.05
```

### 2. 加入外部电路：目录 或 单个文件
```bash
# 整个目录
bash test/run_dm_experiments.sh --nmin 4 --nmax 12 --circuits ~/workshop/circuits

# 单个文件（qubit 数会自动读取，超出 nmin/nmax 范围时仍运行并给出警告）
bash test/run_dm_experiments.sh --circuits ~/workshop/circuits/C17_204.real
```

### 3. 只运行某个实验
```bash
# 仅压缩对比
bash test/run_dm_experiments.sh --experiments "compress" --nmin 5 --nmax 10

# 仅规模对比
bash test/run_dm_experiments.sh --experiments "scale" --nmax 15
```

### 4. 大规模扫描（多种子 + 真实电路）
```bash
bash test/run_dm_experiments.sh \
  --nmin 4 --nmax 12 \
  --noise-p 0.1 \
  --rounds 15 \
  --seeds "42,123,777,999" \
  --circuits ~/workshop/circuits \
  --experiments "scale rounds compress expval"
```

### 5. 直接调用 benchmark 可执行文件（单次实验）
```bash
# 可执行文件位置
BIN=./build_dm_release/test/benchmark_density_matrix

# scale 实验（n=6~10，整个电路目录）
$BIN scale --nmin 6 --nmax 10 --noise-p 0.05 --circuits ~/workshop/circuits

# scale 实验（指定单个文件）
$BIN scale --circuits ~/workshop/circuits/C17_204.real

# compress 实验（单个文件，忽略 nmin/nmax 限制）
$BIN compress --circuits ~/workshop/circuits/ham15_107.real

# rounds 实验（15轮噪声，仅合成电路）
$BIN rounds --nmin 5 --nmax 9 --rounds 15

# compress 实验（仅合成电路）
$BIN compress --nmin 5 --nmax 10 --seeds "42,123"

# expval 实验
$BIN expval --nmin 4 --nmax 8
```

---

## 4 个实验含义

| 实验 | 文件名 | 关键列 | 展示的对比 |
|------|--------|--------|------------|
| **scale** | exp1_scale_latest.csv | `rho_noisy_nodes` vs `dense_nodes` | DD 节点数线性增长，Dense 指数爆炸 |
| **rounds** | exp2_rounds_latest.csv | `dd_size` vs `dense_nodes` per round | 多轮噪声后 DD 不爆炸，trace 精确=1 |
| **compress** | exp3_compress_rho_latest.csv | `ig_ratio`, `sift_ratio` | igGroupSifting 对含噪密度矩阵的压缩效果 |
| **expval** | exp4_expval_latest.csv | `expval_Z0` vs `analytical` | 期望值随噪声衰减，单 qubit 解析解完全吻合 |

---

## 查看与对比结果

结果保存在 `results/density_matrix/`：
- `*_latest.csv`：最新一次运行结果
- `*_YYYYMMDD_HHMMSS.csv`：带时间戳的历史记录

### 基础查看
```bash
# 列表查看所有最新文件
ls -lh results/density_matrix/*_latest.csv

# 表格形式查看
column -t -s',' results/density_matrix/exp1_scale_latest.csv
column -t -s',' results/density_matrix/exp3_compress_rho_latest.csv
```

### 关键对比指令

```bash
# --- Exp1: DD节点 vs Dense节点，计算压缩比 ---
awk -F',' 'NR>1 {printf "n=%-3s  DD=%-6s  Dense=%-12s  ratio=%s\n", $3,$5,$8,$9}' \
  results/density_matrix/exp1_scale_latest.csv

# --- Exp2: 某 n 的多轮噪声演化（n=7）---
awk -F',' 'NR==1 || $3==7' results/density_matrix/exp2_rounds_latest.csv | column -t -s','

# --- Exp2: DD size 全程不超过初始值的验证 ---
awk -F',' 'NR>1' results/density_matrix/exp2_rounds_latest.csv | \
  awk -F',' '{key=$1; if(NR==1||$4>max[key]) max[key]=$4; if(NR==1||$4<min[key]) min[key]=$4}
             END{for(k in max) print k, "max_dd="max[k], "min_dd="min[k]}' | sort

# --- Exp3: igGroupSifting 优于 Sifting 的所有行 ---
awk -F',' 'NR==1 || ($6!="" && $7!="" && $6+0 < $7+0)' \
  results/density_matrix/exp3_compress_rho_latest.csv | column -t -s','

# --- Exp3: 按 ig_ratio 排序，找压缩最显著的 Top10 ---
tail -n +2 results/density_matrix/exp3_compress_rho_latest.csv | \
  sort -t',' -k7 -n | head -10 | column -t -s','

# --- Exp3: 各 n 的平均压缩率 ---
awk -F',' 'NR>1 {sum[$3]+=$7; cnt[$3]++}
           END{for(n in sum) printf "n=%s  avg_ig_ratio=%.4f\n", n, sum[n]/cnt[n]}' \
  results/density_matrix/exp3_compress_rho_latest.csv | sort -t= -k2 -n

# --- Exp4: 解析解误差（单 qubit 行）---
awk -F',' 'NR>1 && $2=="single_qubit" {diff=$5-$6; if(diff<0)diff=-diff;
           printf "p=%-5s  ev=%-10s  analytical=%-10s  error=%s\n", $4,$5,$6,diff}' \
  results/density_matrix/exp4_expval_latest.csv

# --- 两次运行结果对比（时间戳不同）---
diff <(tail -n +2 results/density_matrix/exp1_scale_20260616_160654.csv) \
     <(tail -n +2 results/density_matrix/exp1_scale_latest.csv)
```

---

## 与开源项目的对比实验

### 对比对象

| 项目 | 类型 | 噪声仿真 | 内存模型 |
|------|------|----------|---------|
| **本项目** | DD 密度矩阵 | Kraus 算子 | DD 节点（线性） |
| **Qiskit Aer** `DensityMatrixSimulator` | Dense 密度矩阵 | Kraus/Lindblad | 4^n 复数（指数） |
| **QuTiP** `mesolve` | Dense 密度矩阵 | Lindblad | 4^n 复数（指数） |
| **MQT DDSIM** | DD 纯态仿真 | 不支持混态 | DD 节点（无噪声） |

---

### 对比1：内存占用 vs Qiskit Aer（无需安装，理论计算）

密度矩阵 Dense 内存 = `4^n × 16 bytes`（complex128）：

```bash
# 用我们的 scale 实验直接得出对比数字
column -t -s',' results/density_matrix/exp1_scale_latest.csv | \
  awk 'NR==1 || /QFT/'
# dense_MB 列即为 Qiskit Aer 所需内存，dd_size 列为本项目节点数
```

| n | 本项目 DD 节点 | Qiskit Aer 内存 |
|---|--------------|----------------|
| 8 | ~9 | 1 MB |
| 10 | ~11 | 16 MB |
| 12 | ~13 | 256 MB |
| 14 | ~15 | 4 GB → OOM |
| 16 | ~17 | 64 GB → OOM |

---

### 对比2：与 Qiskit Aer 实测（需要 Python 环境）

**安装依赖（脚本会自动尝试安装，也可手动）：**
```bash
pip install numpy qiskit qiskit-aer psutil
```

**运行 Qiskit 对比脚本：**
```bash
# scale 实验（输出 dense 内存、仿真时间、peak 内存）
python3 test/compare_qiskit.py scale \
  --nmin 4 --nmax 10 --noise-p 0.05 \
  --output results/density_matrix/qiskit_scale.csv

# rounds 实验（多轮噪声 purity 演化）
python3 test/compare_qiskit.py rounds \
  --nmin 4 --nmax 8 --noise-p 0.05 --rounds 10 \
  --output results/density_matrix/qiskit_rounds.csv

# expval 实验（<Z> 解析解验证）
python3 test/compare_qiskit.py expval \
  --nmin 1 --nmax 7 --noise-p 0.05 \
  --output results/density_matrix/qiskit_expval.csv
```

**Qiskit CSV 输出列含义：**

| 实验 | 关键列 |
|------|--------|
| scale  | `dense_MB`=理论内存需求, `sim_time_ms`=仿真时间, `peak_mem_MB`=实测峰值内存 |
| rounds | `purity`=每轮后纯度（应与DD结果一致）, `trace`=每轮后迹 |
| expval | `expval_Z0`=测量值, `analytical`=解析解（单qubit可验证） |

**合并 DD 与 Qiskit 结果为对比表：**
```bash
# 先生成 DD 的 scale 结果（单 seed，便于合并）
./build_dm_release/test/benchmark_density_matrix scale \
  --nmin 4 --nmax 10 --noise-p 0.05 --seeds 42 \
  > /tmp/dd_scale.csv

# 生成 Qiskit 的 scale 结果
python3 test/compare_qiskit.py scale \
  --nmin 4 --nmax 10 --noise-p 0.05 --seeds 42 \
  --output /tmp/qiskit_scale.csv

# 合并对比（输出 dd_nodes, qiskit_dense_nodes, memory_ratio 等）
python3 test/merge_comparison.py \
  --dd /tmp/dd_scale.csv \
  --qiskit /tmp/qiskit_scale.csv \
  --output results/density_matrix/comparison_dd_vs_qiskit.csv

# 查看最终对比表
column -t -s',' results/density_matrix/comparison_dd_vs_qiskit.csv
```


---

### 对比3：与 MQT DDSIM（纯态 vs 含噪）

MQT DDSIM 只支持纯态仿真（不支持 Kraus 噪声通道）。对比维度：

```bash
# 本项目：含噪声的密度矩阵节点数
awk -F',' 'NR>1 && $2=="QFT"' results/density_matrix/exp2_rounds_latest.csv | \
  awk -F',' 'BEGIN{print "n,round,dd_noisy_nodes"} {print $3","$4","$5}' | \
  column -t -s','

# MQT DDSIM 纯态参考：纯态 DD 节点数 = rho_pure_nodes 列
awk -F',' 'NR>1 && $2=="QFT" {print $3","$4}' \
  results/density_matrix/exp1_scale_latest.csv | \
  awk -F',' 'BEGIN{print "n,pure_dd_nodes"} {print}' | column -t -s','

# 结论：纯态 DD (DDSIM 能做) vs 含噪 DD (本项目独有)
# rounds 实验显示：即使经历 10 轮噪声，本项目 DD 节点数仍接近纯态
```

---

### 对比4：igGroupSifting vs 标准 Sifting（本项目内部对比）

```bash
# ig 压缩率 vs sift 压缩率，找 ig 有优势的场景
awk -F',' 'NR>1 {
  ig=$7+0; sift=$8+0;
  if(ig < sift) status="ig_wins";
  else if(ig == sift) status="tie";
  else status="sift_wins";
  print $1","$3","$4","ig","sift","status
}' results/density_matrix/exp3_compress_rho_latest.csv | \
  awk -F',' 'BEGIN{print "source,n,p,ig_ratio,sift_ratio,winner"} {print}' | \
  column -t -s','

# 统计胜负
awk -F',' 'NR>1 {ig=$7+0; sift=$8+0;
  if(ig<sift) ig_win++; else if(ig==sift) tie++; else sift_win++}
  END{print "ig_wins="ig_win" tie="tie" sift_wins="sift_win}' \
  results/density_matrix/exp3_compress_rho_latest.csv
```

---

## 编译说明

脚本会自动检测并编译，也可手动编译：

```bash
mkdir -p build_dm_release && cd build_dm_release
cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_POLICY_VERSION_MINIMUM=3.5 -Wno-dev
make benchmark_density_matrix -j$(nproc)
```

编译日志：`build_dm_release/cmake_configure.log`，`build_dm_release/cmake_build.log`

---

## gtest 单元测试（14 个密度矩阵测试）

```bash
cd build_dm_release
make qfr_test -j$(nproc)
./test/qfr_test --gtest_filter="DensityMatrix*"
```
