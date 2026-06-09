#!/bin/bash
# 用法:
#   ./run_benchmark.sh <电路文件或文件夹> [输出csv文件]
#
# 示例:
#   ./run_benchmark.sh ~/workshop/circuits/                    # 跑整个文件夹
#   ./run_benchmark.sh ~/workshop/circuits/alu4_201.real       # 跑单个文件
#   ./run_benchmark.sh ~/workshop/circuits/ results.csv        # 指定输出文件

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
# 查找 benchmark 可执行文件
BENCHMARK=""
for candidate in \
    "$SCRIPT_DIR/../build_lb/benchmark" \
    "$SCRIPT_DIR/../build/benchmark" \
    "$SCRIPT_DIR/../release/benchmark" \
    "$SCRIPT_DIR/../debugDynamic/benchmark"; do
    if [ -x "$candidate" ]; then
        BENCHMARK="$candidate"
        break
    fi
done

if [ -z "$BENCHMARK" ]; then
    echo "错误: 找不到 benchmark 可执行文件，请先编译项目。" >&2
    echo "  cd build_lb && make benchmark" >&2
    exit 1
fi

INPUT="${1:?请提供电路文件或文件夹路径}"
OUTPUT="${2:-benchmark_results.csv}"

# CSV Header
echo "circuit,initial_size,build_time(s),sifting_size,sifting_time(s),lb_sifting_size,lb_sifting_time(s),upper_linear_size,upper_linear_time(s),lb_upper_linear_size,lb_upper_linear_time(s),lower_linear_size,lower_linear_time(s),lb_lower_linear_size,lb_lower_linear_time(s),mix_linear_size,mix_linear_time(s),lb_mix_linear_size,lb_mix_linear_time(s)" > "$OUTPUT"

# 收集文件列表
FILES=()
if [ -d "$INPUT" ]; then
    while IFS= read -r -d '' f; do
        FILES+=("$f")
    done < <(find "$INPUT" -type f -print0 | sort -z)
elif [ -f "$INPUT" ]; then
    FILES=("$INPUT")
else
    echo "错误: '$INPUT' 不是有效的文件或文件夹" >&2
    exit 1
fi

TOTAL=${#FILES[@]}
if [ "$TOTAL" -eq 0 ]; then
    echo "错误: 在 '$INPUT' 中未找到文件" >&2
    exit 1
fi

echo "共 $TOTAL 个电路文件，结果输出到: $OUTPUT"
echo ""

COUNT=0
FAILED=0
for circuit in "${FILES[@]}"; do
    COUNT=$((COUNT + 1))
    basename=$(basename "$circuit")
    printf "[%d/%d] %s ... " "$COUNT" "$TOTAL" "$basename"

    if result=$("$BENCHMARK" "$circuit" 2>/dev/null); then
        echo "$result" >> "$OUTPUT"
        printf "完成\n"
    else
        printf "失败\n"
        FAILED=$((FAILED + 1))
    fi
done

echo ""
echo "========================================="
echo "完成! 成功: $((TOTAL - FAILED)), 失败: $FAILED"
echo "结果文件: $OUTPUT"
