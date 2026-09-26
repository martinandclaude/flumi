#!/usr/bin/env python3
"""Score deduplication outputs against simulate.py's truth.

A molecule is *observable* when at least one of its primary, non-chimeric
reads (MAPQ >= --min-mapq) carries an RX tag: a perfect tool emits exactly
one representative for each observable molecule and nothing else.

Per run (``--run NAME=out.bam``; representatives are records without the
duplicate flag that carry uS/MI):

* molecules      representatives written
* split          extra representatives of a molecule already represented
                 (under-merging: the molecule is counted more than once)
* lost           observable molecules with no representative (over-merged
                 into another molecule, or all their reads dropped)
* chimeric       representatives that are PCR chimeras (spurious molecules)
* count_err_pct  (molecules - observable) / observable
* gene_mape      median absolute % error of per-gene molecule counts
* gene_ape_p90   90th percentile of that error over genes
* hot_gene_err   % error on the most highly expressed genes
* rep_exact      representatives whose junction chain equals the true
                 molecule's junctions within the aligned span
"""
from __future__ import annotations

import argparse
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path

import pysam


def mol_of(name: str):
    parts = name.split(":")
    return int(parts[0][1:]), len(parts) > 2 and parts[2] == "x"


def junctions(rec):
    out = []
    pos = rec.reference_start
    for op, n in rec.cigartuples or ():
        if op == 3:
            out.append((pos, pos + n))
            pos += n
        elif op in (0, 2, 7, 8):
            pos += n
    return out


def load_truth(path: Path):
    truth = {}
    with open(path) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            j = [] if f[8] == "." else [tuple(int(x) for x in p.split("-")) for p in f[8].split(",")]
            truth[int(f[0])] = (f[1], j)
    return truth


def observable(path: Path, min_mapq: int):
    obs = set()
    with pysam.AlignmentFile(str(path)) as bf:
        for r in bf:
            if r.is_unmapped or r.is_secondary or r.is_supplementary or r.mapping_quality < min_mapq:
                continue
            m, chim = mol_of(r.query_name)
            if not chim and r.has_tag("RX"):
                obs.add(m)
    return obs


def score(path: Path, truth, obs, hot_genes):
    reps = 0
    chim = 0
    seen = Counter()
    exact = 0
    groups = defaultdict(set)
    mols_mi = defaultdict(set)
    with pysam.AlignmentFile(str(path)) as bf:
        for r in bf:
            if r.is_secondary or r.is_supplementary:
                continue
            m, is_chim = mol_of(r.query_name)
            if r.has_tag("MI"):
                mi = r.get_tag("MI")
                groups[mi].add(m)
                mols_mi[m].add(mi)
            if r.is_duplicate or not (r.has_tag("uS") or r.has_tag("MI")):
                continue
            reps += 1
            if is_chim:
                chim += 1
                continue
            seen[m] += 1
            s, e = r.reference_start, r.reference_end
            true_in_span = [j for j in truth[m][1] if j[0] >= s and j[1] <= e]
            exact += junctions(r) == true_in_span
    distinct = len(seen)
    split = sum(seen.values()) - distinct
    lost = len(obs - set(seen))
    true_gene = Counter(truth[m][0] for m in obs)
    pred_gene = Counter()
    for m, k in seen.items():
        pred_gene[truth[m][0]] += k
    errs = [abs(pred_gene[g] - n) / n for g, n in true_gene.items() if n >= 5]
    hot = [100.0 * (pred_gene[g] - true_gene[g]) / true_gene[g] for g in hot_genes if true_gene[g]]
    row = {
        "molecules": reps,
        "split": split,
        "lost": lost,
        "chimeric": chim,
        "count_err_pct": round(100.0 * (reps - len(obs)) / len(obs), 2),
        "gene_mape": round(100.0 * statistics.median(errs), 2) if errs else 0.0,
        "gene_ape_p90": round(100.0 * sorted(errs)[int(0.9 * (len(errs) - 1))], 2) if errs else 0.0,
        "hot_gene_err": ",".join(f"{h:+.1f}" for h in hot) or ".",
        "rep_exact": round(100.0 * exact / max(1, sum(seen.values())), 2),
    }
    if groups:
        row["overmerge_mi"] = sum(len(v) - 1 for v in groups.values())
        row["split_mi"] = sum(len(mols_mi[m]) - 1 for m in obs if mols_mi.get(m))
    return row


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--truth", required=True, type=Path)
    ap.add_argument("--input", required=True, type=Path, help="simulated input BAM")
    ap.add_argument("--run", action="append", default=[], help="NAME=output.bam")
    ap.add_argument("--min-mapq", type=int, default=1)
    ap.add_argument("--hot", type=int, default=3, help="report the N most expressed genes")
    args = ap.parse_args(argv)

    truth = load_truth(args.truth)
    obs = observable(args.input, args.min_mapq)
    gene_n = Counter(truth[m][0] for m in obs)
    hot_genes = [g for g, _ in gene_n.most_common(args.hot)]
    print(f"observable molecules: {len(obs)}  (hot genes: {', '.join(f'{g}={gene_n[g]}' for g in hot_genes)})")
    rows = []
    for spec in args.run:
        name, path = spec.split("=", 1)
        rows.append((name, score(Path(path), truth, obs, hot_genes)))
    keys = []
    for _, r in rows:
        keys += [k for k in r if k not in keys]
    print("| run | " + " | ".join(keys) + " |")
    print("|---|" + "---|" * len(keys))
    for name, r in rows:
        print(f"| {name} | " + " | ".join(str(r.get(k, "")) for k in keys) + " |")
    return 0


if __name__ == "__main__":
    sys.exit(main())
