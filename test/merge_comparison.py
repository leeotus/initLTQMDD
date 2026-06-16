#!/usr/bin/env python3
"""
merge_comparison.py — 合并 DD benchmark 与 Qiskit Aer baseline 结果

用法:
  python3 merge_comparison.py --dd <DD的scale csv> --qiskit <Qiskit的scale csv> [--output FILE]

输出列:
  source, type, n, noise_p,
  dd_nodes, qiskit_dense_nodes,
  qiskit_dense_MB, qiskit_sim_ms, qiskit_peak_MB,
  purity_dd, purity_qiskit, trace_dd, trace_qiskit,
  memory_ratio   (dense_nodes / dd_nodes)
"""

import argparse, csv, sys

def load_csv(path):
    with open(path) as f:
        reader = csv.DictReader(f)
        return list(reader)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dd",     required=True, help="DD benchmark CSV (scale)")
    parser.add_argument("--qiskit", required=True, help="Qiskit baseline CSV (scale)")
    parser.add_argument("--output", default="-")
    args = parser.parse_args()

    dd_rows = {(r["source"], r["n"]): r for r in load_csv(args.dd)}
    qk_rows = {(r["source"], r["n"]): r for r in load_csv(args.qiskit)}

    out = open(args.output, "w") if args.output != "-" else sys.stdout
    writer = csv.writer(out)
    writer.writerow([
        "source", "type", "n", "noise_p",
        "dd_nodes", "qiskit_dense_nodes", "qiskit_dense_MB",
        "dd_sim_ms", "qiskit_sim_ms", "qiskit_peak_MB",
        "purity_dd", "purity_qiskit",
        "trace_dd",  "trace_qiskit",
        "memory_ratio"
    ])

    for key, qkr in sorted(qk_rows.items(), key=lambda x: int(x[0][1])):
        ddr = dd_rows.get(key)
        if ddr is None:
            continue
        n = int(qkr["n"])
        dd_nodes = int(ddr.get("rho_noisy_nodes", ddr.get("rho_pure_nodes", 1)))
        dense_nodes = int(qkr["dense_nodes"])
        ratio = dense_nodes / max(1, dd_nodes)

        writer.writerow([
            key[0],
            qkr.get("type", ""),
            n,
            qkr["noise_p"],
            dd_nodes,
            dense_nodes,
            qkr["dense_MB"],
            ddr.get("sim_time_ms", "N/A"),
            qkr.get("sim_time_ms", "N/A"),
            qkr.get("peak_mem_MB", "N/A"),
            ddr.get("purity", "N/A"),
            qkr.get("purity", "N/A"),
            ddr.get("trace", "N/A"),
            qkr.get("trace", "N/A"),
            f"{ratio:.1f}x"
        ])

    if args.output != "-":
        out.close()
        print(f"[INFO] 对比结果: {args.output}", file=sys.stderr)

if __name__ == "__main__":
    main()
