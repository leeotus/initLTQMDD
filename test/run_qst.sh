#!/bin/bash
#==============================================================================
# run_qst.sh — QST (Quantum State Tomography) 实验脚本
#
# 用法:
#   bash test/run_qst.sh --help
#   bash test/run_qst.sh --output <csv> --rounds <N> --circuit <file> \
#                        [--strategy <name>] [--timeout <sec>] [--memory <MB>]
#
# 参数 (--xxx 格式):
#   --output    输出CSV路径         (必需)
#   --rounds    运行轮次            (必需, 取平均值)
#   --circuit   电路文件路径         (必需, .real / .qasm / 等)
#   --strategy  Sifting策略         (可选, 默认 IGGrpSift)
#               可选: None Sift LBSift IGSift IGLBSift GrpSift IGGrpSift
#   --timeout   单轮超时(秒)        (可选, 默认 300)
#   --memory    内存限制(MB)        (可选, 默认 80000)
#   --qst-bin   QST可执行文件路径   (可选, 默认自动编译到 ./build_qst/qst)
#   --help      显示此帮助信息
#
# 输出CSV格式:
#   Circuit,NQubits,Strategy,Rounds,Fidelity,TraceDistance,RhoSize,PeakDD,TimeMs
#   超时/内存超出/失败时对应列写 "-"
#
# 示例:
#   bash test/run_qst.sh --output results/alu4.csv --rounds 3 \
#       --circuit ~/workshop/ltqmdd/circuits/revLib/alu4_201.real \
#       --strategy IGGrpSift --timeout 300 --memory 80000
#
#   bash test/run_qst.sh --output results/qft6.csv --rounds 1 \
#       --circuit test/circuits/bell.qasm --strategy Sift
#==============================================================================

set -o pipefail

# ============================================================
# 自动编译 QST 二进制文件
# ============================================================
auto_build_qst() {
    local build_dir="${QST_BUILD_DIR:-build_qst}"
    local qst_bin="$build_dir/qst"

    if [ -x "$qst_bin" ]; then
        echo "$qst_bin"
        return 0
    fi

    echo "  [auto] QST 二进制不存在，正在编译 (Release)..."
    if ! command -v cmake &>/dev/null; then
        echo "  错误: cmake 未安装"
        return 1
    fi

    cmake -S . -B "$build_dir" -DCMAKE_BUILD_TYPE=Release > /dev/null 2>&1
    cmake --build "$build_dir" --config Release --target qst > /dev/null 2>&1

    if [ -x "$qst_bin" ]; then
        echo "  [auto] 编译完成: $qst_bin"
        echo "$qst_bin"
        return 0
    else
        echo "  错误: 编译失败"
        return 1
    fi
}

# ============================================================
# 帮助信息
# ============================================================
show_help() {
    sed -n '/^#======/,/^#======/p' "$0" | grep -v '^#!/' | sed 's/^#//; s/^ //'
    exit 0
}

# ============================================================
# 参数解析
# ============================================================
OUTPUT_CSV=""
ROUNDS=""
CIRCUIT=""
STRATEGY="IGGrpSift"
TIMEOUT_SEC=300
MEM_LIMIT_MB=80000
QST_BIN=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --help)    show_help ;;
        --output)  OUTPUT_CSV="$2";  shift 2 ;;
        --rounds)  ROUNDS="$2";      shift 2 ;;
        --circuit) CIRCUIT="$2";     shift 2 ;;
        --strategy) STRATEGY="$2";   shift 2 ;;
        --timeout) TIMEOUT_SEC="$2"; shift 2 ;;
        --memory)  MEM_LIMIT_MB="$2";shift 2 ;;
        --qst-bin) QST_BIN="$2";    shift 2 ;;
        *) echo "未知参数: $1 (用 --help 查看帮助)"; exit 1 ;;
    esac
done

# 检查必需参数
if [ -z "$OUTPUT_CSV" ] || [ -z "$ROUNDS" ] || [ -z "$CIRCUIT" ]; then
    echo "错误: 缺少必需参数 (--output, --rounds, --circuit)"
    echo "使用 --help 查看帮助"
    exit 1
fi

# 自动编译 QST 二进制
if [ -z "$QST_BIN" ]; then
    QST_BIN=$(auto_build_qst) || exit 1
fi

# ============================================================
# 参数校验
# ============================================================
if [ ! -f "$CIRCUIT" ]; then
    echo "错误: 电路文件不存在: $CIRCUIT"
    exit 1
fi

if ! [[ "$ROUNDS" =~ ^[0-9]+$ ]] || [ "$ROUNDS" -lt 1 ]; then
    echo "错误: --rounds 必须是正整数: $ROUNDS"
    exit 1
fi

if ! [[ "$TIMEOUT_SEC" =~ ^[0-9]+$ ]]; then
    echo "错误: --timeout 必须是正整数: $TIMEOUT_SEC"
    exit 1
fi

if ! [[ "$MEM_LIMIT_MB" =~ ^[0-9]+$ ]]; then
    echo "错误: --memory 必须是正整数: $MEM_LIMIT_MB"
    exit 1
fi

VALID_STRATS="None Sift LBSift IGSift IGLBSift GrpSift IGGrpSift"
if ! echo "$VALID_STRATS" | grep -qw "$STRATEGY"; then
    echo "错误: 无效的 --strategy: $STRATEGY"
    echo "  有效选项: $VALID_STRATS"
    exit 1
fi

# ============================================================
# 电路信息
# ============================================================
CIRCUIT_NAME=$(basename "$CIRCUIT" | sed 's/\.[^.]*$//')
NQUBITS=$(grep "^\.numvars" "$CIRCUIT" 2>/dev/null | head -1 | awk '{print $2}')
if [ -z "$NQUBITS" ]; then
    NQUBITS=$("$QST_BIN" "$CIRCUIT" 1 1 2>&1 | grep -oP 'Qubits:\s*\K\d+' | head -1)
    [ -z "$NQUBITS" ] && NQUBITS="?"
fi

MEM_LIMIT_KB=$((MEM_LIMIT_MB * 1024))

# ============================================================
# 打印配置
# ============================================================
echo "=============================================="
echo "  QST 实验脚本"
echo "=============================================="
echo "  电路:       $CIRCUIT_NAME  (${NQUBITS}Q)"
echo "  策略:       $STRATEGY"
echo "  轮次:       $ROUNDS"
echo "  超时:       ${TIMEOUT_SEC}s"
echo "  内存限制:   ${MEM_LIMIT_MB}MB"
echo "  输出CSV:    $OUTPUT_CSV"
echo "  QST程序:    $QST_BIN"
echo "=============================================="
echo ""

# ============================================================
# 写 CSV 表头
# ============================================================
OUTDIR=$(dirname "$OUTPUT_CSV")
mkdir -p "$OUTDIR"

if [ ! -f "$OUTPUT_CSV" ]; then
    echo "Circuit,NQubits,Strategy,Rounds,Fidelity,TraceDistance,RhoSize,PeakDD,TimeMs" > "$OUTPUT_CSV"
fi

# ============================================================
# 运行实验
# ============================================================
SUM_FID=0; SUM_TD=0; SUM_RHO=0; SUM_PEAK=0; SUM_TIME=0; VALID_ROUNDS=0
ALL_FAILED=true

for ((r=1; r<=ROUNDS; r++)); do
    echo -n "  [轮 $r/$ROUNDS] 运行中... "

    FULL_LOG="/tmp/qst_run_${CIRCUIT_NAME}_r${r}.log"
    START_TS=$(date +%s%N)

    (
        ulimit -v "$MEM_LIMIT_KB" 2>/dev/null
        timeout "${TIMEOUT_SEC}" "$QST_BIN" "$CIRCUIT" 10 10 2>/dev/null
    ) > "$FULL_LOG" 2>&1
    EXIT_CODE=$?

    END_TS=$(date +%s%N)

    if [ $EXIT_CODE -ne 0 ]; then
        case $EXIT_CODE in
            124) echo "超时 (${TIMEOUT_SEC}s)" ;;
            137|139) echo "内存超出/段错误" ;;
            *) echo "运行失败 (exit=$EXIT_CODE)" ;;
        esac
        rm -f "$FULL_LOG"
        continue
    fi

    # 查找 QST 程序生成的 CSV (格式: qst_<safe_circuit_name>.csv)
    SAFE_NAME=$(echo "$CIRCUIT" | sed 's|/|_|g; s|\.[^.]*$||')
    QST_CSV=$(ls "qst_${CIRCUIT_NAME}.csv" 2>/dev/null | head -1)
    [ -z "$QST_CSV" ] && QST_CSV=$(ls "qst_${SAFE_NAME}.csv" 2>/dev/null | head -1)

    if [ -z "$QST_CSV" ] || [ ! -f "$QST_CSV" ]; then
        echo "未找到 QST 输出CSV (stdout见 $FULL_LOG)"
        continue
    fi

    LINE=$(grep "^${STRATEGY}," "$QST_CSV" 2>/dev/null | head -1)
    if [ -z "$LINE" ]; then
        echo "策略 $STRATEGY 在输出中未找到"
        rm -f "$QST_CSV" "$FULL_LOG"
        continue
    fi

    FID=$(echo "$LINE" | cut -d, -f2)
    TD=$(echo  "$LINE" | cut -d, -f3)
    RHO=$(echo "$LINE" | cut -d, -f4)
    PEAK=$(echo "$LINE"| cut -d, -f5)
    TM=$(echo  "$LINE" | cut -d, -f6)

    echo "Fid=$FID Rho=$RHO Time=${TM}ms"

    SUM_FID=$(echo "$SUM_FID + $FID" | bc -l 2>/dev/null || echo "0")
    SUM_TD=$(echo  "$SUM_TD  + $TD"  | bc -l 2>/dev/null || echo "0")
    SUM_RHO=$((SUM_RHO + RHO))
    SUM_PEAK=$((SUM_PEAK + PEAK))
    SUM_TIME=$(echo "$SUM_TIME + $TM" | bc -l 2>/dev/null || echo "0")
    VALID_ROUNDS=$((VALID_ROUNDS + 1))
    ALL_FAILED=false

    rm -f "$QST_CSV" "$FULL_LOG"
done

# ============================================================
# 写结果行
# ============================================================
if $ALL_FAILED; then
    echo "  所有轮次均失败, 写入占位行"
    echo "${CIRCUIT_NAME},${NQUBITS},${STRATEGY},0/${ROUNDS},-,-,-,-,-" >> "$OUTPUT_CSV"
else
    AVG_FID=$(echo "scale=6; $SUM_FID / $VALID_ROUNDS" | bc -l 2>/dev/null || echo "0")
    AVG_TD=$(echo  "scale=6; $SUM_TD  / $VALID_ROUNDS" | bc -l 2>/dev/null || echo "0")
    AVG_RHO=$((SUM_RHO / VALID_ROUNDS))
    AVG_PEAK=$((SUM_PEAK / VALID_ROUNDS))
    AVG_TIME=$(echo "scale=2; $SUM_TIME / $VALID_ROUNDS" | bc -l 2>/dev/null || echo "0")
    echo "${CIRCUIT_NAME},${NQUBITS},${STRATEGY},${VALID_ROUNDS}/${ROUNDS},${AVG_FID},${AVG_TD},${AVG_RHO},${AVG_PEAK},${AVG_TIME}" >> "$OUTPUT_CSV"
fi

echo ""
echo "  结果已写入: $OUTPUT_CSV"
echo "  有效轮次: $VALID_ROUNDS / $ROUNDS"