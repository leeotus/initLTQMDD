#!/bin/bash
# IG Sifting vs Sifting 对比测试
#
# 用法:
#   ./run_benchmark_ig_sifting.sh <电路文件或文件夹> [输出csv文件]
#
# 示例:
#   ./run_benchmark_ig_sifting.sh ~/circuits/alu4_201.real                # 单文件
#   ./run_benchmark_ig_sifting.sh ~/circuits/                             # 整个目录
#   ./run_benchmark_ig_sifting.sh ~/circuits/ ig_results.csv              # 指定输出文件

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

# 查找可执行文件
BENCHMARK=""
for candidate in \
    "$SCRIPT_DIR/../build_test/benchmark_ig_sifting" \
    "$SCRIPT_DIR/../build_test/apps/benchmark_ig_sifting" \
    "$SCRIPT_DIR/../release/benchmark_ig_sifting" \
    "$SCRIPT_DIR/../build/benchmark_ig_sifting"; do
    if [ -x "$candidate" ]; then
        BENCHMARK="$candidate"
        break
    fi
done

if [ -z "$BENCHMARK" ]; then
    echo "错误: 找不到 benchmark_ig_sifting 可执行文件，请先编译:" >&2
    echo "  cd build_test && make benchmark_ig_sifting -j\$(nproc)" >&2
    exit 1
fi

INPUT="${1:?请提供电路文件或文件夹路径}"
OUTPUT="${2:-ig_sifting_results.csv}"

# CSV Header
echo "circuit,qubits,init_size,sift_size,sift_time_ms,sift_ex,lb_size,lb_time_ms,lb_ex,iglb_size,iglb_time_ms,iglb_ex" > "$OUTPUT"

# 收集文件列表
FILES=()
if [ -d "$INPUT" ]; then
    while IFS= read -r -d '' f; do
        FILES+=("$f")
    done < <(find "$INPUT" -type f \( -name "*.real" -o -name "*.qasm" -o -name "*.tfc" -o -name "*.qc" \) -print0 | sort -z)
elif [ -f "$INPUT" ]; then
    FILES=("$INPUT")
else
    echo "错误: '$INPUT' 不是有效的文件或文件夹" >&2
    exit 1
fi

TOTAL=${#FILES[@]}
if [ "$TOTAL" -eq 0 ]; then
    echo "错误: 在 '$INPUT' 中未找到电路文件 (.real/.qasm/.tfc/.qc)" >&2
    exit 1
fi

echo "可执行文件: $BENCHMARK"
echo "输入: $INPUT ($TOTAL 个电路文件)"
echo "输出: $OUTPUT"
echo ""

COUNT=0
FAILED=0
for circuit in "${FILES[@]}"; do
    COUNT=$((COUNT + 1))
    basename=$(basename "$circuit")
    printf "[%d/%d] %-40s " "$COUNT" "$TOTAL" "$basename"

    if result=$("$BENCHMARK" "$circuit" 2>/dev/null); then
        echo "$result" >> "$OUTPUT"
        sift_sz=$(echo "$result" | cut -d',' -f4)
        lb_sz=$(echo "$result" | cut -d',' -f7)
        iglb_sz=$(echo "$result" | cut -d',' -f10)
        printf "sift=%s  lb=%s  iglb=%s\n" "$sift_sz" "$lb_sz" "$iglb_sz"
    else
        printf "失败\n"
        FAILED=$((FAILED + 1))
    fi
done

echo ""
echo "========================================="
echo "完成! 成功: $((COUNT - FAILED))/$TOTAL, 失败: $FAILED"
echo "结果: $OUTPUT"
