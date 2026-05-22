#!/bin/bash

# ==========================================
# 配置区域
# ==========================================
# 可执行文件路径
EXEC="./releaseStatic/apps/ltqmdd"
# 测试用例目录
CIRCUIT_DIR="../ltqmdd/circuits/revLib"
# 输出的CSV文件名（带时间戳防覆盖）
OUTPUT_CSV="static_sifting_results_$(date +%Y%m%d_%H%M%S).csv"
# 超时时间设置（1小时 = 3600秒）
TIMEOUT_LIMIT=3600

# ==========================================
# 初始化CSV文件并写入表头
# ==========================================
echo "File,Initial_DD_Size,Initial_Time,Sifting_DD_Size,Sifting_Time" > "$OUTPUT_CSV"

# ==========================================
# 遍历目录下的所有文件
# ==========================================
for file in "$CIRCUIT_DIR"/*; do
    # 提取文件名用于显示和记录
    filename=$(basename "$file")
    echo "正在处理: $filename ..."

    # 运行程序，使用 timeout 限制最大执行时间为 TIMEOUT_LIMIT 秒
    # 2>&1 将标准错误合并到标准输出
    output=$(timeout "$TIMEOUT_LIMIT" $EXEC "$file" 2>&1)
    exit_code=$?

    # 判断退出状态码
    # 124 是 timeout 命令专用的超时退出码
    if [ $exit_code -eq 124 ]; then
        echo "  -> 运行超时（超过 $((TIMEOUT_LIMIT/3600)) 小时），记录为 '-'"
        echo "$filename,-,-,-,-" >> "$OUTPUT_CSV"
        continue
    fi

    # 非0退出码表示运行出错（如段错误、断言失败等）
    if [ $exit_code -ne 0 ]; then
        echo "  -> 运行出错，记录为 '-'"
        echo "$filename,-,-,-,-" >> "$OUTPUT_CSV"
        continue
    fi

    # 检查输出中是否包含预期的特征字符串，防止程序未报错退出但输出异常
    if ! echo "$output" | grep -q "初始DD大小"; then
        echo "  -> 输出异常，记录为 '-'"
        echo "$filename,-,-,-,-" >> "$OUTPUT_CSV"
        continue
    fi

    # ==========================================
    # 解析输出数据
    # ==========================================
    init_info=$(echo "$output" | grep "初始DD大小")
    init_size=$(echo "$init_info" | grep -oP '初始DD大小:\s*\K[0-9]+')
    init_time=$(echo "$init_info" | grep -oP '构造时间:\s*\K[0-9.]+')

    sift_info=$(echo "$output" | grep "sifting直到收敛")
    sift_size=$(echo "$sift_info" | grep -oP 'sifting直到收敛:\s*\K[0-9]+')
    sift_time=$(echo "$sift_info" | grep -oP '花费时间:\s*\K[0-9.]+')

    # 将解析到的数据写入CSV
    echo "$filename,$init_size,${init_time}s,$sift_size,${sift_time}s" >> "$OUTPUT_CSV"
    echo "  -> 成功: Init=$init_size, Sift=$sift_size"

done

echo "=========================================="
echo "所有文件处理完毕！结果已保存至: $OUTPUT_CSV"