#!/bin/bash
#==============================================================================
# run_qst.sh — QST (Quantum State Tomography) 实验脚本
#
# 用法:
#   bash test/run_qst.sh --help
#   bash test/run_qst.sh --output <csv> --rounds <N> --circuit <file> [选项...]
#
# 必需参数:
#   --output    输出CSV路径
#   --rounds    运行轮次 (取平均值；超时/失败轮次记 "-")
#   --circuit   电路文件路径 (.real / .qasm 等)
#
# 可选参数:
#   --strategy  Sifting策略 (默认 all)
#               可选: all None Sift LBSift IGSift IGLBSift GrpSift IGGrpSift
#
#   --tomo      层析类型 (默认 auto)
#               complete : 强制完整层析，使用全部 3^N 个 Pauli 基
#                          适合 N<=5，时间长但约束完备
#               partial  : 强制局部层析，使用 --bases 指定数量的随机基
#                          适合 N>=5，速度快但结果近似
#               auto     : N<=4 完整层析，N>=5 局部层析（默认）
#
#   --bases     测量基数量 (默认 auto)
#               auto     : tomo=complete 时为 3^N，tomo=partial 时为 50
#               数字     : 直接指定基数量（若 >= 3^N 自动视为完整层析）
#
#   --iters     MLE 迭代次数 (默认 auto: N<=4 用 50，N>=5 用 30)
#
#   --mode-b    是否运行 Mode B (2N-var 层析) (默认 auto)
#               0=禁用, 1=启用
#               auto: N<=5 时启用，N>=6 时禁用（Mode B 在大 N 时极慢）
#
#   --timeout   单轮超时秒数 (默认 300)；超时记录为 "-"
#   --memory    内存限制 MB (默认 10000)
#   --qst-bin   QST 可执行文件路径 (默认自动编译到 ./build_qst/qst)
#   --help      显示帮助
#
# 输出CSV格式:
#   Circuit,NQubits,Tomography,Strategy,MeasBases,ValidBases,
#   Fidelity,TraceDistance,RhoSize,PeakDD,TimeMs
#
# 示例:
#   # 3Q 完整层析（auto 自动选择）
#   bash test/run_qst.sh --output results/3q.csv --rounds 1 \
#       --circuit ~/workshop/circuits/peres_9.real
#
#   # 5Q 完整层析（强制，时间较长，timeout 保护）
#   bash test/run_qst.sh --output results/5q_full.csv --rounds 1 \
#       --circuit ~/workshop/circuits/4gt4-v0_73.real \
#       --tomo complete --mode-b 0 --timeout 600 --memory 10000
#
#   # 6Q 局部层析，100 个随机基
#   bash test/run_qst.sh --output results/6q.csv --rounds 1 \
#       --circuit ~/workshop/circuits/ex3_229.real \
#       --tomo partial --bases 100 --mode-b 0 --timeout 300
#
#   # 大规模批量（多电路追加到同一 CSV）
#   for C in ~/workshop/circuits/3_17_13.real \
#             ~/workshop/circuits/4_49_16.real \
#             ~/workshop/circuits/4gt4-v0_73.real; do
#     bash test/run_qst.sh --output results/batch.csv --rounds 1 --circuit "$C" \
#         --mode-b 0 --timeout 120
#   done
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
TOMO_ARG="auto"
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
        --tomo)     TOMO_ARG="$2";     shift 2 ;;
        --qst-bin)  QST_BIN="$2";      shift 2 ;;
        *) echo "未知参数: $1 (用 --help 查看帮助)"; exit 1 ;;
    esac
done

# 必需参数检查
if [ -z "$OUTPUT_CSV" ] || [ -z "$ROUNDS" ] || [ -z "$CIRCUIT" ]; then
    echo "错误: 缺少必需参数 (--output, --rounds, --circuit)"; exit 1
fi
if [ ! -f "$CIRCUIT" ]; then
    echo "错误: 电路文件不存在: $CIRCUIT"; exit 1
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

# ============================================================
# 解析量子比特数
# ============================================================
CIRCUIT_NAME=$(basename "$CIRCUIT" | sed 's/\.[^.]*$//')

NQUBITS=$(grep "^\.numvars" "$CIRCUIT" 2>/dev/null | head -1 | awk '{print $2}')
if [ -z "$NQUBITS" ]; then
    NQUBITS=$(grep -E "^qreg" "$CIRCUIT" 2>/dev/null | grep -oP '\[\K[0-9]+' \
              | awk '{s+=$1} END{print s}')
fi
if [ -z "$NQUBITS" ] || ! [[ "$NQUBITS" =~ ^[0-9]+$ ]]; then
    NQUBITS=$("$QST_BIN" "$CIRCUIT" 1 1 2>/dev/null | grep -oP 'N=\K[0-9]+' | head -1)
    [ -z "$NQUBITS" ] && NQUBITS="?"
fi

# ============================================================
# 决定层析类型 / bases / iters / mode-b
# ============================================================
POW3N="?"
TOMO_TYPE="partial"

if [[ "$NQUBITS" =~ ^[0-9]+$ ]]; then
    POW3N=$(pow3 "$NQUBITS")

    # --- 确定层析类型 ---
    if [ "$TOMO_ARG" = "complete" ]; then
        TOMO_TYPE="complete"
    elif [ "$TOMO_ARG" = "partial" ]; then
        TOMO_TYPE="partial"
    else
        # auto: N<=4 完整，N>=5 局部
        [ "$NQUBITS" -le 4 ] && TOMO_TYPE="complete" || TOMO_TYPE="partial"
    fi

    # --- bases ---
    if [ "$BASES_ARG" = "auto" ]; then
        if [ "$TOMO_TYPE" = "complete" ]; then
            BASES_ARG=$POW3N
        else
            BASES_ARG=50    # 局部层析默认 50 个随机基
        fi
    else
        # 用户明确指定数字：若 >= 3^N 则实际是完整层析
        if [[ "$BASES_ARG" =~ ^[0-9]+$ ]] && [ "$BASES_ARG" -ge "$POW3N" ]; then
            TOMO_TYPE="complete"
        fi
    fi

    # --- iters ---
    if [ "$ITERS_ARG" = "auto" ]; then
        [ "$NQUBITS" -le 4 ] && ITERS_ARG=50 || ITERS_ARG=30
    fi

    # --- mode-b ---
    if [ "$MODE_B_ARG" = "auto" ]; then
        [ "$NQUBITS" -le 5 ] && MODE_B_ARG=1 || MODE_B_ARG=0
    fi
else
    [ "$BASES_ARG" = "auto" ] && BASES_ARG=30
    [ "$ITERS_ARG" = "auto" ] && ITERS_ARG=50
    [ "$MODE_B_ARG" = "auto" ] && MODE_B_ARG=1
fi

MEM_LIMIT_KB=$((MEM_LIMIT_MB * 1024))

# ============================================================
# 打印配置
# ============================================================
echo "=============================================="
echo "  QST 实验脚本"
echo "=============================================="
echo "  电路:         $CIRCUIT_NAME (${NQUBITS}Q)"
echo "  层析类型:     $TOMO_TYPE  (3^N=${POW3N}, 实际基数=${BASES_ARG})"
echo "  MLE 迭代:     $ITERS_ARG"
echo "  Mode B:       $MODE_B_ARG  (0=禁用, 1=启用)"
echo "  策略:         $STRATEGY"
echo "  轮次:         $ROUNDS"
echo "  超时:         ${TIMEOUT_SEC}s/轮  (超时记 \"-\")"
echo "  内存限制:     ${MEM_LIMIT_MB}MB"
echo "  输出CSV:      $OUTPUT_CSV"
echo "=============================================="
echo ""

# ============================================================
# CSV 表头（首次创建）
# ============================================================
OUTDIR=$(dirname "$OUTPUT_CSV")
mkdir -p "$OUTDIR"
if [ ! -f "$OUTPUT_CSV" ]; then
    echo "Circuit,NQubits,Tomography,Strategy,MeasBases,ValidBases,Fidelity,TraceDistance,RhoSize,PeakDD,TimeMs" \
        > "$OUTPUT_CSV"
fi

# ============================================================
# 单策略运行函数
# ============================================================
run_strategy() {
    local strat_csv="$1"
    echo "  --- 策略: $strat_csv ---"

    local SUM_FID=0 SUM_TD=0 SUM_RHO=0 SUM_PEAK=0 SUM_TIME=0
    local SUM_MEAS=0 SUM_VALID=0 VALID_ROUNDS=0 ALL_FAILED=true

    for ((r=1; r<=ROUNDS; r++)); do
        echo -n "    [轮 $r/$ROUNDS] 运行中... "

        local FULL_LOG="/tmp/qst_run_${CIRCUIT_NAME}_${strat_csv}_r${r}_$$.log"

        (
            ulimit -v "$MEM_LIMIT_KB" 2>/dev/null
            timeout "${TIMEOUT_SEC}" "$QST_BIN" \
                "$CIRCUIT" "$BASES_ARG" "$ITERS_ARG" "$MODE_B_ARG" 2>/dev/null
        ) > "$FULL_LOG" 2>&1
        local EXIT_CODE=$?

        if [ $EXIT_CODE -ne 0 ]; then
            case $EXIT_CODE in
                124) echo "超时 (${TIMEOUT_SEC}s)" ;;
                137|139) echo "内存超出/段错误" ;;
                *) echo "运行失败 (exit=$EXIT_CODE)" ;;
            esac
            rm -f "$FULL_LOG"; continue
        fi

        # 定位 qst_app 生成的 CSV
        # qst_app 命名规则: qst_<路径中/替换为_，去扩展名>.csv
        local SAFE_NAME QST_CSV
        SAFE_NAME=$(echo "$CIRCUIT" | sed 's|/|_|g; s|\.[^.]*$||')
        QST_CSV=""
        [ -f "qst_${SAFE_NAME}.csv"    ] && QST_CSV="qst_${SAFE_NAME}.csv"
        [ -z "$QST_CSV" ] && [ -f "qst_${CIRCUIT_NAME}.csv" ] \
            && QST_CSV="qst_${CIRCUIT_NAME}.csv"
        if [ -z "$QST_CSV" ]; then
            QST_CSV=$(ls qst_*${CIRCUIT_NAME}*.csv 2>/dev/null | head -1)
        fi

        if [ -z "$QST_CSV" ] || [ ! -f "$QST_CSV" ]; then
            echo "未找到 QST 输出CSV"; rm -f "$FULL_LOG"; continue
        fi

        # 解析策略行 (Strategy,Fidelity,TraceDistance,RhoSize,PeakDD,TimeMs)
        local LINE
        LINE=$(grep "^${strat_csv}," "$QST_CSV" 2>/dev/null | head -1)
        if [ -z "$LINE" ]; then
            echo "策略 $strat_csv 未找到"
            rm -f "$QST_CSV" "$FULL_LOG"; continue
        fi

        local FID TD RHO PEAK TM
        FID=$(echo  "$LINE" | cut -d, -f2)
        TD=$(echo   "$LINE" | cut -d, -f3)
        RHO=$(echo  "$LINE" | cut -d, -f4)
        PEAK=$(echo "$LINE" | cut -d, -f5)
        TM=$(echo   "$LINE" | cut -d, -f6)

        # 解析有效基数量
        local VALID_B MEAS_B
        VALID_B=$(grep "Valid measurement bases:" "$FULL_LOG" \
                  | grep -oP '\d+(?= /)' | head -1)
        MEAS_B=$(grep  "Valid measurement bases:" "$FULL_LOG" \
                  | grep -oP '(?<= / )\d+' | head -1)
        [ -z "$VALID_B" ] && VALID_B="$BASES_ARG"
        [ -z "$MEAS_B"  ] && MEAS_B="$BASES_ARG"

        echo "Fid=$FID VB=$VALID_B/$MEAS_B Rho=$RHO Time=${TM}ms"

        SUM_FID=$(echo  "$SUM_FID + $FID" | bc -l 2>/dev/null || echo "0")
        SUM_TD=$(echo   "$SUM_TD  + $TD"  | bc -l 2>/dev/null || echo "0")
        SUM_RHO=$((SUM_RHO   + ${RHO:-0}))
        SUM_PEAK=$((SUM_PEAK + ${PEAK:-0}))
        SUM_TIME=$(echo "$SUM_TIME + $TM" | bc -l 2>/dev/null || echo "0")
        SUM_VALID=$((SUM_VALID + ${VALID_B:-0}))
        SUM_MEAS=$((SUM_MEAS   + ${MEAS_B:-0}))
        VALID_ROUNDS=$((VALID_ROUNDS + 1))
        ALL_FAILED=false

        rm -f "$QST_CSV" "$FULL_LOG"
    done

    # 写结果行
    if $ALL_FAILED; then
        echo "  所有轮次均失败（超时/OOM/错误）"
        echo "${CIRCUIT_NAME},${NQUBITS},${TOMO_TYPE},${strat_csv},${BASES_ARG},-,-,-,-,-,-" \
            >> "$OUTPUT_CSV"
    else
        local AVG_FID AVG_TD AVG_RHO AVG_PEAK AVG_TIME AVG_MEAS AVG_VALID
        AVG_FID=$(echo  "scale=6; $SUM_FID / $VALID_ROUNDS" | bc -l 2>/dev/null || echo "0")
        AVG_TD=$(echo   "scale=6; $SUM_TD  / $VALID_ROUNDS" | bc -l 2>/dev/null || echo "0")
        AVG_RHO=$((SUM_RHO   / VALID_ROUNDS))
        AVG_PEAK=$((SUM_PEAK / VALID_ROUNDS))
        AVG_TIME=$(echo "scale=2; $SUM_TIME / $VALID_ROUNDS" | bc -l 2>/dev/null || echo "0")
        AVG_MEAS=$((SUM_MEAS  / VALID_ROUNDS))
        AVG_VALID=$((SUM_VALID / VALID_ROUNDS))
        echo "${CIRCUIT_NAME},${NQUBITS},${TOMO_TYPE},${strat_csv},${AVG_MEAS},${AVG_VALID},${AVG_FID},${AVG_TD},${AVG_RHO},${AVG_PEAK},${AVG_TIME}" \
            >> "$OUTPUT_CSV"
        echo "    -> 平均保真度: $AVG_FID  有效轮次: $VALID_ROUNDS/$ROUNDS"
    fi
    echo ""
}

# ============================================================
# 主循环
# ============================================================
for strat in "${RUN_STRATS[@]}"; do
    run_strategy "$strat"
done

echo "  结果已写入: $OUTPUT_CSV"
