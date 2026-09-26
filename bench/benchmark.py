#!/usr/bin/env python3
"""Simulate several regimes, run flumi on each (filter and mark mode) and
score the output against the truth.

    python3 bench/benchmark.py --work /tmp/flumi-bench [--only typical_d3]
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
import evaluate  # noqa: E402

SCENARIOS = {
    # name: simulate.py arguments
    "typical_d3": ["--profile", "typical", "--depth", "3"],
    "noisy_d3": ["--profile", "noisy", "--depth", "3"],
    "typical_d8": ["--profile", "typical", "--depth", "8"],
    "good_d1.5": ["--profile", "good", "--depth", "1.5"],
}


def run(cmd):
    t = time.perf_counter()
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return time.perf_counter() - t


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--work", required=True, type=Path)
    ap.add_argument("--genes", type=int, default=800)
    ap.add_argument("--hot-molecules", type=int, default=20000)
    ap.add_argument("--seed", type=int, default=11)
    ap.add_argument("--only", action="append", help="run only these scenarios")
    ap.add_argument("--flumi", type=Path, default=ROOT / "target" / "release" / "flumi")
    ap.add_argument("--extra", default="", help="extra arguments for flumi")
    args = ap.parse_args(argv)

    rows = []
    for name, sim_args in SCENARIOS.items():
        if args.only and name not in args.only:
            continue
        d = args.work / name
        if not (d / "reads.bam").exists():
            subprocess.run(
                [sys.executable, str(HERE / "simulate.py"), "--out", str(d), "--genes", str(args.genes),
                 "--hot-genes", "1", "--hot-molecules", str(args.hot_molecules), "--seed", str(args.seed)] + sim_args,
                check=True,
            )
        bam, gtf = d / "reads.bam", d / "genes.gtf"
        cmd = [str(args.flumi), "--bam", str(bam), "--gtf", str(gtf), "--threads", "2"] + args.extra.split()
        t = run(cmd + ["--out", str(d / "filter.bam"), "--stats", str(d / "stats.tsv")])
        run(cmd + ["--out", str(d / "mark.bam"), "--mode", "mark"])

        truth = evaluate.load_truth(d / "truth.tsv")
        obs = evaluate.observable(bam, 1)
        gene_n = evaluate.Counter(truth[m][0] for m in obs)
        hot = [g for g, _ in gene_n.most_common(1)]
        r = evaluate.score(d / "filter.bam", truth, obs, hot)
        r["overmerge"] = evaluate.score(d / "mark.bam", truth, obs, hot).get("overmerge_mi", "")
        rows.append((name, len(obs), r, t))

    keys = ["molecules", "split", "lost", "overmerge", "chimeric", "count_err_pct", "gene_mape", "gene_ape_p90",
            "hot_gene_err", "rep_exact"]
    print("| scenario | true molecules | " + " | ".join(keys) + " | runtime s |")
    print("|---|---|" + "---|" * (len(keys) + 1))
    for name, n, r, t in rows:
        print(f"| {name} | {n} | " + " | ".join(str(r.get(k, "")) for k in keys) + f" | {t:.1f} |")
    return 0


if __name__ == "__main__":
    sys.exit(main())
