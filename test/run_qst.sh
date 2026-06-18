#!/bin/bash
#==============================================================================
# run_qst.sh — QST (Quantum State Tomography) 实验脚本
#
# 用法:
#   bash test/run_qst.sh --help
#   bash test/run_qst.sh --output <csv> --rounds <N> --circuit <file> [选项...]
#   bash test/run_qst.sh --output <csv> --rounds <N> --circuit <dir>  [选项...]
#
# ────────────────────────────────────────────────────────────────
# 必需参数
# ────────────────────────────────────────────────────────────────
#   --output <path>
#       输出 CSV 文件路径。若文件不存在自动创建并写表头；
#       若文件已存在则追加（可多次调用累积多个电路的结果）。
#
#   --rounds <N>
#       每个策略运行 N 轮，对所有成功轮取平均值。
#       失败轮（超时/OOM/错误）不计入平均，对应结果记 "-"。
#
#   --circuit <file|dir>
#       量子电路文件路径，支持 .real（RevLib 格式）和 .qasm（OpenQASM）。
#       若传入目录，则自动处理目录内所有 .real 文件。
#
# ────────────────────────────────────────────────────────────────
# 层析策略参数
# ────────────────────────────────────────────────────────────────
#   --strategy <name>  （默认: all）
#       指定要测试的变量重排（Sifting）策略。
#       "all" 会依次运行全部 7 种策略：
#         NoReorder  : 不做任何变量重排（基线）
#         Sift       : 经典 Sifting
#         LBSift     : Lower-Bound Sifting（依下界剪枝）
#         IGSift     : Interaction-Graph Sifting（按 qubit 交互图排序）
#         IGLBSift   : IG + LB 联合 Sifting
#         GrpSift    : Group Sifting（按对称群分组）
#         IGGrpSift  : IG + Group 联合 Sifting
#
#   --tomo <complete|partial|auto>  （默认: auto）
#       控制测量基的选取方式：
#         complete : 强制完整层析，枚举全部 3^N 个 Pauli 基
#                    每个基再枚举所有 2^N 个 outcome，总 projector 数 = 6^N
#                    约束方程完备，MLE 有唯一解，保真度可达 1.0
#                    适合 N≤4；N=5 可用但每轮约需数分钟
#         partial  : 强制局部层析，随机抽取 --bases 个 Pauli 基
#                    约束欠定，保真度是近似值；适合 N≥5 快速实验
#         auto     : N≤4 自动用 complete，N≥5 自动用 partial（默认）
#
#   --bases <N|auto>  （默认: auto）
#       测量 Pauli 基数量：
#         auto     : complete 时取 3^N，partial 时取 50
#         数字     : 直接指定；若 >= 3^N 则自动识别为完整层析
#       注意：程序对每个基都会枚举所有 2^N 个 outcome，
#             实际 projector 总数 = bases × 2^N
#
#   --iters <N|auto>  （默认: auto）
#       MLE（最大似然估计）迭代次数：
#         auto     : N≤4 取 50，N≥5 取 30
#         数字     : 直接指定；次数越多收敛越好，但时间线性增长
#
#   --mode-b <0|1|auto>  （默认: auto）
#       是否运行 Mode B（2N 变量密度矩阵层析）：
#         0        : 禁用 Mode B（推荐大规模实验使用）
#         1        : 启用 Mode B
#         auto     : N≤5 时启用，N≥6 时自动禁用（极慢）
#       Mode B 在 2N 变量空间重建密度矩阵，理论更严格但开销远高于 Mode A，
#       实际使用中 N>5 时几乎不可行。
#
# ────────────────────────────────────────────────────────────────
# 同步重排参数（MLE 迭代中途对 rho 做 Sifting）
# ────────────────────────────────────────────────────────────────
#   --sync-interval <N>  （默认: 0）
#       每隔 N 轮 MLE 迭代触发一次同步重排。
#       0 = 不使用此触发条件。
#       推荐值：5-10。
#
#   --sync-threshold <N>  （默认: 0）
#       当前 rho 的 DD 节点数超过 N 时触发同步重排。
#       0 = 不使用此触发条件。
#       推荐值：当 NoReorder 的 RhoSize 为 X 时，可设置 2X-3X。
#
#   两个参数满足任一条件即触发，均为 0 则不做同步重排（默认）。
#   触发后对 rho 做 --strategy 指定的 Sifting，
#   由于 DD package 的 Unique Table 是 package 级别共享的，
#   所有 projector 会自动随 rho 同步重排，无需额外处理。
#   CSV 的 Tomography 字段会附加 +syncI<N> 或 +syncT<N> 标记。
#   注意：NoReorder 策略不受同步重排参数影响（无法重排）。
#
# ────────────────────────────────────────────────────────────────
# 噪声信道参数
# ────────────────────────────────────────────────────────────────
#   --noise <type>        （默认: "" = 无噪声）
#       噪声信道类型，字符串组合：
#         D : 去极化信道（Depolarizing）
#             Kraus: sqrt(1-p)*I, sqrt(p/3)*X, sqrt(p/3)*Y, sqrt(p/3)*Z
#         A : 振幅阻尼（Amplitude Damping）
#             Kraus: [[1,0],[0,sqrt(1-p)]], [[0,sqrt(p)],[0,0]]
#         P : 相位翻转（Phase Flip）
#             Kraus: sqrt(1-p)*I, sqrt(p)*Z
#       噪声作用于电路输出纯态，生成含噪混合态，作为 QST 重建的目标。
#       保真度 F = <psi|rho_recon|psi> 衡量重建的含噪态与理想纯态的距离。
#
#   --noise-prob <p>      （默认: 0.01）
#       噪声概率（每个 qubit 独立同分布，范围 (0,1)）。
#       CSV 的 Tomography 字段会附加 +noise<type><p> 标记。
#
# ────────────────────────────────────────────────────────────────
# 资源限制参数
# ────────────────────────────────────────────────────────────────
#   --timeout <sec>   （默认: 300）
#       单轮超时秒数。超时的轮次记录为 "-" 而非报错退出。
#       建议大规模实验设置较大值（600-1800）。
#
#   --memory <MB>     （默认: 10000 = 10GB）
#       进程虚拟内存限制（MB）。超出时进程被杀，该轮记 "-"。
#       实际物理内存使用通常低于虚拟内存限制。
#
#   --qst-bin <path>  （默认: 自动编译到 ./build_qst/qst）
#       指定已编译好的 qst 可执行文件路径，跳过自动编译。
#
# ────────────────────────────────────────────────────────────────
# 所有参数速查表
# ────────────────────────────────────────────────────────────────
#
#   参数               默认值          说明
#   ──────────────────────────────────────────────────────────────
#   --output <path>    （必需）        输出CSV文件路径（追加模式）
#   --rounds <N>       （必需）        运行轮次，取平均值
#   --circuit <file|dir> （必需）      量子电路文件(.real/.qasm)或目录
#   ──────────────────────────────────────────────────────────────
#   --strategy <name>  all             Sifting策略
#                                      all/NoReorder/Sift/LBSift/
#                                      IGSift/IGLBSift/GrpSift/IGGrpSift
#   --tomo <type>      auto            层析类型
#                                      complete/partial/auto
#   --bases <N|auto>   auto            测量Pauli基数量
#   --iters <N|auto>   auto            MLE迭代次数
#   --mode-b <0|1|auto> auto           Mode B (2N变量层析)
#                                      auto: N<=5启用，N>=6禁用
#   ──────────────────────────────────────────────────────────────
#   --sync-interval  <N>  0            每N轮MLE触发同步重排（0=关）
#   --sync-threshold <N>  0            rho节点数>N时触发重排（0=关）
#                                      两条件OR关系，均0则不重排
#   ──────────────────────────────────────────────────────────────
#   --noise <type>     ""              噪声信道类型（空=无噪声）
#                                      D=去极化 A=振幅阻尼 P=相位翻转
#   --noise-prob <p>   0.01            噪声概率（每qubit，范围(0,1)）
#   ──────────────────────────────────────────────────────────────
#   --timeout <sec>    300             单轮超时（超时记"-"）
#   --memory <MB>      10000           内存限制（OOM记"-"）
#   --qst-bin <path>   auto-build      qst可执行文件路径
#   --help                             显示本帮助
#
# ────────────────────────────────────────────────────────────────
# 输出 CSV 格式说明
# ────────────────────────────────────────────────────────────────
#   Circuit      : 电路名称（从文件名提取）
#   NQubits      : 量子比特数 N
#   Tomography   : 层析标签，格式示例：
#                    complete
#                    partial
#                    complete+syncI5            （轮次触发重排）
#                    complete+syncT30           （节点数触发重排）
#                    complete+syncI5+syncT30    （双触发）
#                    complete+noiseD0.02        （去极化噪声）
#                    partial+noiseA0.05+syncI5  （组合）
#   Strategy     : Sifting 策略名
#   MeasBases    : 实际投影算子数 = bases × 2^N
#   ValidBases   : 有效投影数（概率 > 1e-10），是MLE约束数量
#   Fidelity     : 保真度 F = <ψ|ρ_recon|ψ>，越接近1越好
#                  含噪实验中：F衡量重建含噪态与理想纯态的距离
#   TraceDistance: 迹距离 ≈ sqrt(1-F)，越小越好
#   RhoSize      : 最终 ρ 的 DD 节点数（Sifting压缩后）
#   PeakDD       : MLE过程峰值 DD 节点数（内存占用指标）
#   TimeMs       : 总耗时（毫秒，含建图+MLE+Sifting+保真度计算）
#   MaxRSS_MB    : 进程最大物理内存（MB），多轮取平均
#   TotalTimeSec : 脚本处理该电路的总耗时（秒）
#
# ────────────────────────────────────────────────────────────────
# 使用示例
# ────────────────────────────────────────────────────────────────
#   # 1. 3Q 完整层析，全部策略（auto 自动配置）
#   bash test/run_qst.sh \
#       --output results/3q.csv --rounds 1 \
#       --circuit ~/workshop/circuits/peres_9.real
#
#   # 2. 4Q 完整层析，同步重排（每5轮触发）
#   bash test/run_qst.sh \
#       --output results/4q_sync.csv --rounds 1 \
#       --circuit ~/workshop/circuits/mini-alu_167.real \
#       --tomo complete --sync-interval 5 --mode-b 0
#
#   # 3. 4Q 完整层析，节点数超30时触发重排
#   bash test/run_qst.sh \
#       --output results/4q_thresh.csv --rounds 1 \
#       --circuit ~/workshop/circuits/mini-alu_167.real \
#       --tomo complete --sync-threshold 30 --mode-b 0
#
#   # 4. 双触发条件（满足任一即重排）
#   bash test/run_qst.sh \
#       --output results/4q_both.csv --rounds 1 \
#       --circuit ~/workshop/circuits/mini-alu_167.real \
#       --tomo complete --sync-interval 10 --sync-threshold 50 --mode-b 0
#
#   # 5. 5Q 局部层析，禁用 Mode B，增大超时
#   bash test/run_qst.sh \
#       --output results/5q.csv --rounds 1 \
#       --circuit ~/workshop/circuits/4gt4-v0_73.real \
#       --tomo partial --bases 80 --iters 30 --mode-b 0 --timeout 600
#
#   # 6. 含噪声实验（去极化，p=2%）
#   bash test/run_qst.sh \
#       --output results/noise_D.csv --rounds 3 \
#       --circuit ~/workshop/circuits/peres_9.real \
#       --tomo complete --mode-b 0 \
#       --noise D --noise-prob 0.02
#
#   # 7. 含噪声 + 同步重排组合实验
#   bash test/run_qst.sh \
#       --output results/noise_sync.csv --rounds 1 \
#       --circuit ~/workshop/circuits/mini-alu_167.real \
#       --tomo complete --mode-b 0 \
#       --noise A --noise-prob 0.05 \
#       --sync-interval 5 --sync-threshold 30
#
#   # 8. 批量多电路追加到同一 CSV
#   for C in ~/workshop/circuits/3_17_13.real \
#             ~/workshop/circuits/4_49_16.real \
#             ~/workshop/circuits/4gt4-v0_73.real; do
#     bash test/run_qst.sh \
#         --output results/batch.csv --rounds 1 \
#         --circuit "$C" --mode-b 0 --timeout 300
#   done
#
#   # 9. 仅测单个策略（快速验证）
#   bash test/run_qst.sh \
#       --output results/quick.csv --rounds 1 \
#       --circuit ~/workshop/circuits/peres_9.real \
#       --strategy IGSift
#
#   # 10. 比较有无噪声（追加到同一CSV便于对比）
#   for NOISE in "" "D" "A"; do
#     bash test/run_qst.sh \
#         --output results/noise_compare.csv --rounds 1 \
#         --circuit ~/workshop/circuits/peres_9.real \
#         --tomo complete --mode-b 0 \
#         --noise "$NOISE" --noise-prob 0.03
#   done
#
#   # 11. 目录模式：处理文件夹内所有 .real 文件
#   bash test/run_qst.sh \
#       --output results/q3_batch.csv --rounds 5 \
#       --circuit ./assets/q3 --strategy IGSift --tomo complete --mode-b 0
#==============================================================================

set -o pipefail

# ============================================================
# 自动编译
# ============================================================
auto_build_qst() {
    local build_dir="${QST_BUILD_DIR:-build_qst}"
    local qst_bin="$build_dir/qst"
    if [ -x "$qst_bin" ]; then echo "$qst_bin"; return 0; fi
    echo "  [auto] 编译 QST..." >&2
    cmake -S . -B "$build_dir" -DCMAKE_BUILD_TYPE=Release > /dev/null 2>&1
    cmake --build "$build_dir" --target qst > /dev/null 2>&1
    if [ -x "$qst_bin" ]; then echo "$qst_bin"; return 0; fi
    echo "  错误: 编译失败" >&2; return 1
}

show_help() {
    sed -n '/^#======/,/^#======/p' "$0" | grep -v '^#!/' | sed 's/^#//; s/^ //'
    exit 0
}

pow3() {
    local n=$1 r=1
    for ((i=0;i<n;i++)); do r=$((r*3)); done
    echo $r
}

# ============================================================
# 参数解析
# ============================================================
OUTPUT_CSV="" ROUNDS="" CIRCUIT="" STRATEGY="all"
TIMEOUT_SEC=300 MEM_LIMIT_MB=10000
BASES_ARG="auto" ITERS_ARG="auto" MODE_B_ARG="auto"
TOMO_ARG="auto" SYNC_INTERVAL=0 SYNC_THRESHOLD=0
NOISE_TYPE="" NOISE_PROB=0.01
QST_BIN=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --help)     show_help ;;
        --output)   OUTPUT_CSV="$2";   shift 2 ;;
        --rounds)   ROUNDS="$2";       shift 2 ;;
        --circuit)  CIRCUIT="$2";      shift 2 ;;
        --strategy) STRATEGY="$2";     shift 2 ;;
        --timeout)  TIMEOUT_SEC="$2";  shift 2 ;;
        --memory)   MEM_LIMIT_MB="$2"; shift 2 ;;
        --bases)    BASES_ARG="$2";    shift 2 ;;
        --iters)    ITERS_ARG="$2";    shift 2 ;;
        --mode-b)   MODE_B_ARG="$2";   shift 2 ;;
        --tomo)           TOMO_ARG="$2";       shift 2 ;;
        --sync-interval)  SYNC_INTERVAL="$2";  shift 2 ;;
        --sync-threshold) SYNC_THRESHOLD="$2"; shift 2 ;;
        --noise)          NOISE_TYPE="$2";     shift 2 ;;
        --noise-prob)     NOISE_PROB="$2";     shift 2 ;;
        --qst-bin)        QST_BIN="$2";        shift 2 ;;
        *) echo "未知参数: $1 (用 --help 查看帮助)"; exit 1 ;;
    esac
done

# 必需参数检查
if [ -z "$OUTPUT_CSV" ] || [ -z "$ROUNDS" ] || [ -z "$CIRCUIT" ]; then
    echo "错误: 缺少必需参数 (--output, --rounds, --circuit)"; exit 1
fi
if [ ! -f "$CIRCUIT" ] && [ ! -d "$CIRCUIT" ]; then
    echo "错误: 电路文件或目录不存在: $CIRCUIT"; exit 1
fi
if ! [[ "$TOMO_ARG" =~ ^(complete|partial|auto)$ ]]; then
    echo "错误: --tomo 只能是 complete / partial / auto"; exit 1
fi

# 自动编译
if [ -z "$QST_BIN" ]; then QST_BIN=$(auto_build_qst) || exit 1; fi

# 策略校验
VALID_STRATS="all None NoReorder Sift LBSift IGSift IGLBSift GrpSift IGGrpSift"
if ! echo "$VALID_STRATS" | grep -qw "$STRATEGY"; then
    echo "错误: 无效 --strategy: $STRATEGY"; exit 1
fi

ALL_STRATS=("NoReorder" "Sift" "LBSift" "IGSift" "IGLBSift" "GrpSift" "IGGrpSift")
if [ "$STRATEGY" = "all" ]; then
    RUN_STRATS=("${ALL_STRATS[@]}")
else
    [ "$STRATEGY" = "None" ] && STRATEGY="NoReorder"
    RUN_STRATS=("$STRATEGY")
fi

MEM_LIMIT_KB=$((MEM_LIMIT_MB * 1024))

# ============================================================
# CSV 表头（首次创建）
# ============================================================
OUTDIR=$(dirname "$OUTPUT_CSV")
mkdir -p "$OUTDIR"
if [ ! -f "$OUTPUT_CSV" ]; then
    echo "Circuit,NQubits,Tomography,Strategy,MeasBases,ValidBases,Fidelity,TraceDistance,RhoSize,PeakDD,TimeMs,MaxRSS_MB,TotalTimeSec" \
        > "$OUTPUT_CSV"
fi

# ============================================================
# run_one_circuit — 处理单个电路文件的核心函数
# ============================================================
run_one_circuit() {
    local CIRCUIT_FILE="$1"
    if [ ! -f "$CIRCUIT_FILE" ]; then
        echo "  跳过: 文件不存在 $CIRCUIT_FILE" >&2
        return 1
    fi

    local CIRCUIT_NAME
    CIRCUIT_NAME=$(basename "$CIRCUIT_FILE" | sed 's/\.[^.]*$//')

    # ── 解析 qubit 数 ──
    local NQUBITS
    NQUBITS=$(grep "^\.numvars" "$CIRCUIT_FILE" 2>/dev/null | head -1 | awk '{print $2}')
    if [ -z "$NQUBITS" ]; then
        NQUBITS=$(grep -E "^qreg" "$CIRCUIT_FILE" 2>/dev/null | grep -oP '\[\K[0-9]+' \
                  | awk '{s+=$1} END{print s}')
    fi
    if [ -z "$NQUBITS" ] || ! [[ "$NQUBITS" =~ ^[0-9]+$ ]]; then
        NQUBITS=$("$QST_BIN" "$CIRCUIT_FILE" 1 1 2>/dev/null | grep -oP 'N=\K[0-9]+' | head -1)
        [ -z "$NQUBITS" ] && NQUBITS="?"
    fi

    # ── 确定层析类型 / bases / iters / mode-b ──
    local POW3N="?"
    local TOMO_TYPE="partial"

    if [[ "$NQUBITS" =~ ^[0-9]+$ ]]; then
        POW3N=$(pow3 "$NQUBITS")

        if [ "$TOMO_ARG" = "complete" ]; then
            TOMO_TYPE="complete"
        elif [ "$TOMO_ARG" = "partial" ]; then
            TOMO_TYPE="partial"
        else
            [ "$NQUBITS" -le 4 ] && TOMO_TYPE="complete" || TOMO_TYPE="partial"
        fi

        local BASES_ARG ITERS_ARG MODE_B_ARG
        BASES_ARG="${SCRIPT_BASES_ARG:-auto}"
        if [ "$BASES_ARG" = "auto" ]; then
            [ "$TOMO_TYPE" = "complete" ] && BASES_ARG=$POW3N || BASES_ARG=50
        elif [[ "$BASES_ARG" =~ ^[0-9]+$ ]] && [ "$BASES_ARG" -ge "$POW3N" ]; then
            TOMO_TYPE="complete"
        fi
        ITERS_ARG="${SCRIPT_ITERS_ARG:-auto}"
        [ "$ITERS_ARG" = "auto" ] && { [ "$NQUBITS" -le 4 ] && ITERS_ARG=50 || ITERS_ARG=30; }
        MODE_B_ARG="${SCRIPT_MODE_B_ARG:-auto}"
        [ "$MODE_B_ARG" = "auto" ] && { [ "$NQUBITS" -le 5 ] && MODE_B_ARG=1 || MODE_B_ARG=0; }
    else
        BASES_ARG="${SCRIPT_BASES_ARG:-30}"
        ITERS_ARG="${SCRIPT_ITERS_ARG:-50}"
        MODE_B_ARG="${SCRIPT_MODE_B_ARG:-1}"
    fi

    # ── 打印配置 ──
    echo "=============================================="
    echo "  QST 实验脚本"
    echo "=============================================="
    echo "  电路:         $CIRCUIT_NAME (${NQUBITS}Q)"
    echo "  层析类型:     $TOMO_TYPE  (3^N=${POW3N}, 实际基数=${BASES_ARG})"
    echo "  MLE 迭代:     $ITERS_ARG"
    echo "  同步重排:     interval=${SYNC_INTERVAL} threshold=${SYNC_THRESHOLD}  (均为0则不重排)"
    if [ -n "$NOISE_TYPE" ]; then
        echo "  噪声信道:     type=${NOISE_TYPE} prob=${NOISE_PROB}"
    else
        echo "  噪声信道:     无（理想纯态）"
    fi
    echo "  Mode B:       $MODE_B_ARG  (0=禁用, 1=启用)"
    echo "  策略:         $STRATEGY"
    echo "  轮次:         $ROUNDS"
    echo "  超时:         ${TIMEOUT_SEC}s/轮  (超时记 \"-\")"
    echo "  内存限制:     ${MEM_LIMIT_MB}MB"
    echo "  输出CSV:      $OUTPUT_CSV"
    echo "=============================================="
    echo ""

    # ── TOMO_LABEL ──
    local TOMO_LABEL="$TOMO_TYPE"
    [ "${SYNC_INTERVAL:-0}" -gt 0 ] 2>/dev/null && TOMO_LABEL="${TOMO_LABEL}+syncI${SYNC_INTERVAL}"
    [ "${SYNC_THRESHOLD:-0}" -gt 0 ] 2>/dev/null && TOMO_LABEL="${TOMO_LABEL}+syncT${SYNC_THRESHOLD}"
    [ -n "$NOISE_TYPE" ] && TOMO_LABEL="${TOMO_LABEL}+noise${NOISE_TYPE}${NOISE_PROB}"

    # ── 关联数组 ──
    declare -A S_FID S_TD S_RHO S_PEAK S_TIME S_VALID S_MEAS S_CNT S_FAIL S_RSS
    for strat in "${RUN_STRATS[@]}"; do
        S_FID[$strat]=0; S_TD[$strat]=0; S_RHO[$strat]=0
        S_PEAK[$strat]=0; S_TIME[$strat]=0
        S_VALID[$strat]=0; S_MEAS[$strat]=0
        S_CNT[$strat]=0; S_FAIL[$strat]=true
        S_RSS[$strat]=0
    done
    local GLOBAL_VB="-"; local GLOBAL_MEAS="-"
    local SCRIPT_START_SEC=$SECONDS

    # ── 运行 ROUNDS 轮 ──
    for ((r=1; r<=ROUNDS; r++)); do
        echo "  [轮 $r/$ROUNDS] 运行 qst..."
        local FULL_LOG="/tmp/qst_main_${CIRCUIT_NAME}_r${r}_$$.log"
        local TIME_LOG="/tmp/qst_time_${CIRCUIT_NAME}_r${r}_$$.log"

        (
            ulimit -v "$MEM_LIMIT_KB" 2>/dev/null
            /usr/bin/time -v -o "$TIME_LOG" timeout "${TIMEOUT_SEC}" "$QST_BIN" \
                "$CIRCUIT_FILE" "$BASES_ARG" "$ITERS_ARG" "$MODE_B_ARG" \
                "$SYNC_INTERVAL" "$SYNC_THRESHOLD" \
                "$NOISE_TYPE" "$NOISE_PROB" 2>/dev/null
        ) > "$FULL_LOG" 2>&1
        local EXIT_CODE=$?
        local RSS_MB=0
        if [ -f "$TIME_LOG" ]; then
            RSS_MB=$(grep "Maximum resident set size" "$TIME_LOG" 2>/dev/null | awk '{printf "%.1f", $NF/1024}')
            [ -z "$RSS_MB" ] && RSS_MB=0
            rm -f "$TIME_LOG"
        fi
        echo "    RSS: ${RSS_MB} MB"

        if [ $EXIT_CODE -ne 0 ]; then
            case $EXIT_CODE in
                124) echo "    超时 (${TIMEOUT_SEC}s)" ;;
                137|139) echo "    内存超出/段错误" ;;
                *) echo "    运行失败 (exit=$EXIT_CODE)" ;;
            esac
            rm -f "$FULL_LOG"; continue
        fi

        # 定位 qst_app 生成的 CSV
        local SAFE_NAME
        SAFE_NAME=$(echo "$CIRCUIT_FILE" | sed 's|/|_|g; s|\.[^.]*$||')
        local QST_CSV=""
        [ -f "qst_${SAFE_NAME}.csv"        ] && QST_CSV="qst_${SAFE_NAME}.csv"
        [ -z "$QST_CSV" ] && [ -f "qst_${CIRCUIT_NAME}.csv" ] && QST_CSV="qst_${CIRCUIT_NAME}.csv"
        [ -z "$QST_CSV" ] && QST_CSV=$(ls qst_*${CIRCUIT_NAME}*.csv 2>/dev/null | head -1)

        if [ -z "$QST_CSV" ] || [ ! -f "$QST_CSV" ]; then
            echo "    未找到 QST 输出CSV"; rm -f "$FULL_LOG"; continue
        fi

        # 解析 ValidBases
        local VB MB
        VB=$(grep "Valid measurement bases:" "$FULL_LOG" | grep -oP '\d+(?= /)' | head -1)
        MB=$(grep "Valid measurement bases:" "$FULL_LOG" | grep -oP '(?<= / )\d+' | head -1)
        [ -n "$VB" ] && GLOBAL_VB="$VB"
        [ -n "$MB" ] && GLOBAL_MEAS="$MB"

        # 逐策略读取
        for strat in "${RUN_STRATS[@]}"; do
            local LINE
            LINE=$(grep "^${strat}," "$QST_CSV" 2>/dev/null | head -1)
            if [ -z "$LINE" ]; then
                echo "    [$strat] 未在 CSV 中找到"
                continue
            fi
            local FID TD RHO PEAK TM
            FID=$(echo "$LINE" | cut -d, -f2)
            TD=$(echo  "$LINE" | cut -d, -f3)
            RHO=$(echo "$LINE" | cut -d, -f4)
            PEAK=$(echo "$LINE"| cut -d, -f5)
            TM=$(echo  "$LINE" | cut -d, -f6)
            echo "    [$strat] Fid=$FID Rho=$RHO Time=${TM}ms"

            S_FID[$strat]=$(awk "BEGIN{printf \"%.6f\", ${S_FID[$strat]} + $FID}")
            S_TD[$strat]=$(awk "BEGIN{printf \"%.6f\", ${S_TD[$strat]} + $TD}")
            S_RHO[$strat]=$((${S_RHO[$strat]}   + ${RHO:-0}))
            S_PEAK[$strat]=$((${S_PEAK[$strat]} + ${PEAK:-0}))
            S_TIME[$strat]=$(awk "BEGIN{printf \"%.2f\", ${S_TIME[$strat]} + $TM}")
            S_VALID[$strat]=$((${S_VALID[$strat]} + ${VB:-0}))
            S_MEAS[$strat]=$((${S_MEAS[$strat]}  + ${MB:-0}))
            S_CNT[$strat]=$((${S_CNT[$strat]} + 1))
            S_RSS[$strat]=$(awk "BEGIN{printf \"%.1f\", ${S_RSS[$strat]} + ${RSS_MB:-0}}")
            S_FAIL[$strat]=false
        done

        rm -f "$QST_CSV" "$FULL_LOG"
        echo ""
    done

    # ── 计算本轮总耗时 ──
    local SCRIPT_TOTAL_SEC=$(( SECONDS - SCRIPT_START_SEC ))

    # ── 写各策略结果行 ──
    for strat in "${RUN_STRATS[@]}"; do
        local n=${S_CNT[$strat]}
        if ${S_FAIL[$strat]} || [ "$n" -eq 0 ]; then
            echo "  [$strat] 所有轮次均失败"
            echo "${CIRCUIT_NAME},${NQUBITS},${TOMO_LABEL},${strat},${BASES_ARG},-,-,-,-,-,-,-,${SCRIPT_TOTAL_SEC}" \
                >> "$OUTPUT_CSV"
        else
            local AVG_FID AVG_TD AVG_RHO AVG_PEAK AVG_TIME AVG_VB AVG_MB AVG_RSS
            AVG_FID=$(awk "BEGIN{printf \"%.6f\", ${S_FID[$strat]} / $n}")
            AVG_TD=$(awk "BEGIN{printf \"%.6f\", ${S_TD[$strat]} / $n}")
            AVG_RHO=$((${S_RHO[$strat]}  / n))
            AVG_PEAK=$((${S_PEAK[$strat]} / n))
            AVG_TIME=$(awk "BEGIN{printf \"%.2f\", ${S_TIME[$strat]} / $n}")
            AVG_VB=$((${S_VALID[$strat]} / n))
            AVG_MB=$((${S_MEAS[$strat]}  / n))
            AVG_RSS=$(awk "BEGIN{printf \"%.1f\", ${S_RSS[$strat]} / $n}")
            echo "${CIRCUIT_NAME},${NQUBITS},${TOMO_LABEL},${strat},${AVG_MB},${AVG_VB},${AVG_FID},${AVG_TD},${AVG_RHO},${AVG_PEAK},${AVG_TIME},${AVG_RSS},${SCRIPT_TOTAL_SEC}" \
                >> "$OUTPUT_CSV"
            echo "  [$strat] 平均保真度=$AVG_FID  有效轮次=$n/$ROUNDS  RSS=${AVG_RSS}MB  总耗时=${SCRIPT_TOTAL_SEC}s"
        fi
    done

    echo ""
    echo "  本电路耗时: ${SCRIPT_TOTAL_SEC}s"
    echo "  结果已写入: $OUTPUT_CSV"
}

# ============================================================
# 入口：判断 --circuit 是目录还是文件
# ============================================================
if [ -d "$CIRCUIT" ]; then
    CIRCUIT_FILES=()
    while IFS= read -r -d '' f; do
        CIRCUIT_FILES+=("$f")
    done < <(find "$CIRCUIT" -maxdepth 1 -name '*.real' -print0 | sort -z)
    if [ ${#CIRCUIT_FILES[@]} -eq 0 ]; then
        echo "错误: 目录中没有 .real 文件: $CIRCUIT"
        exit 1
    fi
    echo "=============================================="
    echo "  [目录模式] 找到 ${#CIRCUIT_FILES[@]} 个 .real 文件"
    echo "  输出CSV: $OUTPUT_CSV"
    echo "=============================================="
    echo ""
    OVERALL_START=$SECONDS
    for cf in "${CIRCUIT_FILES[@]}"; do
        echo "==== 处理: $(basename "$cf") ===="
        run_one_circuit "$cf"
        echo ""
    done
    OVERALL_SEC=$(( SECONDS - OVERALL_START ))
    echo "=============================================="
    echo "  全部完成! 总共 ${#CIRCUIT_FILES[@]} 个电路, 总耗时: ${OVERALL_SEC}s"
    echo "=============================================="
else
    run_one_circuit "$CIRCUIT"
fi