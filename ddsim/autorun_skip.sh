#!/bin/bash

# ================= 配置区 =================
EXECUTABLE="./build/ddsim_noise_aware"
TARGET_DIR="../ltqmdd/circuits/revLib"
OUTPUT_CSV="benchmark_results.csv"
TIMEOUT_SECONDS=3600 # 1小时 = 3600秒
# ==========================================

# 如果CSV文件不存在，则写入CSV表头；若已存在则保留原数据
if [ ! -f "$OUTPUT_CSV" ]; then
    echo "Filename,Static_Duration,DynamicSifting_Duration" > "$OUTPUT_CSV"
fi

# 遍历目标目录下的所有文件
for file in "$TARGET_DIR"/*; do
    # 提取不带路径的文件名
    filename=$(basename "$file")
    
    # 检查该文件名是否已经在CSV中存在记录，避免重复运行
    if grep -q "^${filename}," "$OUTPUT_CSV"; then
        echo "跳过已完成的文件: $filename"
        continue
    fi

    echo "正在处理: $filename ..."

    # 初始化两个变量的值，默认为"-"（代表超时或失败）
    duration_0="-"
    duration_1="-"

    # ---------------- 运行命令 1 (dynamic_reorder 0) ----------------
    output_0=$(timeout "$TIMEOUT_SECONDS" "$EXECUTABLE" --noise_effects APD --stoch_runs 20000 --noise_prob 0.005 --dynamic_reorder 0 --post_reorder 1 --simulate_file "$file" 2>&1)
    exit_status_0=$?

    if [ $exit_status_0 -eq 0 ]; then
        # 提取时间，并去除所有的回车符(\r)和换行符(\n)
        extracted_0=$(echo "$output_0" | grep -oP 'Duration simulation: \K.*' | tr -d '\r\n')
        # xargs 可以去除字符串首尾的空格和制表符，确保字符串干净
        if [ -n "$extracted_0" ]; then
            duration_0=$(echo "$extracted_0" | xargs)
        fi
    else
        echo "  [!] 命令1 (reorder 0) 超时或执行失败，退出码: $exit_status_0"
    fi

    # ---------------- 运行命令 2 (dynamic_reorder 1) ----------------
    output_1=$(timeout "$TIMEOUT_SECONDS" "$EXECUTABLE" --noise_effects APD --stoch_runs 20000 --noise_prob 0.005 --dynamic_reorder 1 --post_reorder 1 --simulate_file "$file" 2>&1)
    exit_status_1=$?

    if [ $exit_status_1 -eq 0 ]; then
        # 提取时间，并去除所有的回车符(\r)和换行符(\n)
        extracted_1=$(echo "$output_1" | grep -oP 'Duration simulation: \K.*' | tr -d '\r\n')
        if [ -n "$extracted_1" ]; then
            duration_1=$(echo "$extracted_1" | xargs)
        fi
    else
        echo "  [!] 命令2 (reorder 1) 超时或执行失败，退出码: $exit_status_1"
    fi

    # ---------------- 写入 CSV ----------------
    # 将时间字符串中可能存在的逗号替换为空格，防止CSV列错位
    duration_0_csv=$(echo "$duration_0" | tr ',' ' ')
    duration_1_csv=$(echo "$duration_1" | tr ',' ' ')
    
    # 使用 printf 替代 echo 写入，避免 echo 自身的换行行为干扰
    printf "%s,%s,%s\n" "$filename" "$duration_0_csv" "$duration_1_csv" >> "$OUTPUT_CSV"

done

echo "================================================"
echo "所有测试已完成！结果已保存至: $OUTPUT_CSV"