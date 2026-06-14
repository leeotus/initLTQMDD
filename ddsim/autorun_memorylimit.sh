#!/bin/bash

# ================= 配置区 =================
EXECUTABLE="./build/ddsim_noise_aware"
TARGET_DIR="../ltqmdd/circuits/revLib"
OUTPUT_CSV="benchmark_results.csv"
TIMEOUT_SECONDS=3600      # 1小时 = 3600秒
MEMORY_LIMIT_KB=31457280  # 30GB = 30 * 1024 * 1024 KB = 31457280 KB
# ==========================================

# 写入CSV表头
echo "Filename,Static_Duration,DynamicSifting_Duration" > "$OUTPUT_CSV"

# 遍历目标目录下的所有文件
for file in "$TARGET_DIR"/*; do
    filename=$(basename "$file")
    echo "正在处理: $filename ..."

    duration_0="-"
    duration_1="-"

    # ---------------- 运行命令 1 (dynamic_reorder 0) ----------------
    # 使用子进程 () 执行，ulimit 只影响该子进程，防止主脚本被限制
    # 2>/dev/null 用于屏蔽由于内存超限导致的 Bash 内存分配错误提示
    output_0=$(
        ulimit -v "$MEMORY_LIMIT_KB"
        timeout "$TIMEOUT_SECONDS" "$EXECUTABLE" --noise_effects APD --stoch_runs 20000 --noise_prob 0.005 --dynamic_reorder 0 --simulate_file "$file" 2>&1
    )
    exit_status_0=$?

    if [ $exit_status_0 -eq 0 ]; then
        extracted_0=$(echo "$output_0" | grep -oP 'Duration simulation: \K.*' | tr -d '\r\n')
        if [ -n "$extracted_0" ]; then
            duration_0=$(echo "$extracted_0" | xargs)
        fi
    else
        # 137 = 被SIGKILL杀掉 (可能是内存超限被OOM Killer杀掉，或timeout超时被杀)
        # 134 = 被SIGABRT中止 (C++的std::bad_alloc通常触发此信号)
        echo "  [!] 命令1 (reorder 0) 超时或内存超限，退出码: $exit_status_0"
    fi

    # ---------------- 运行命令 2 (dynamic_reorder 1) ----------------
    output_1=$(
        ulimit -v "$MEMORY_LIMIT_KB"
        timeout "$TIMEOUT_SECONDS" "$EXECUTABLE" --noise_effects APD --stoch_runs 20000 --noise_prob 0.005 --dynamic_reorder 1 --simulate_file "$file" 2>&1
    )
    exit_status_1=$?

    if [ $exit_status_1 -eq 0 ]; then
        extracted_1=$(echo "$output_1" | grep -oP 'Duration simulation: \K.*' | tr -d '\r\n')
        if [ -n "$extracted_1" ]; then
            duration_1=$(echo "$extracted_1" | xargs)
        fi
    else
        echo "  [!] 命令2 (reorder 1) 超时或内存超限，退出码: $exit_status_1"
    fi

    # ---------------- 写入 CSV ----------------
    duration_0_csv=$(echo "$duration_0" | tr ',' ' ')
    duration_1_csv=$(echo "$duration_1" | tr ',' ' ')
    printf "%s,%s,%s\n" "$filename" "$duration_0_csv" "$duration_1_csv" >> "$OUTPUT_CSV"

done

echo "================================================"
echo "所有测试已完成！结果已保存至: $OUTPUT_CSV"