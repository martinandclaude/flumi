#!/usr/bin/env python3
"""Summarise flumi runs on real data, one directory per sample.

Each directory is the output of run_sample.sh. Per sample: tag rates, flumi statistics, an independent check of gene labels
and full-length calls against the representatives' CIGAR blocks, and a UMI
neighbour test in the densest 5' neighbourhood (singletons 1 or 2
substitutions from a molecule with >= 3 reads, against UMIs drawn from other
loci). Across samples: agreement of per-gene molecule counts.

    summarize.py --gtf genes.gtf.gz [--chrom chr22] [--json results.json] fc1_bc01 ...
"""
import argparse
import collections
import csv
import gzip
import itertools
import json
import math
import os
import random
import re
import subprocess
import sys


def stats(path):
    with open(path) as f:
        return {r[0]: r[1] for r in csv.reader(f, delimiter="\t") if r and r[0] != "metric"}


def blocks_of(cigar, pos):
    p, b = pos - 1, []
    for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        n = int(n)
        if op in "M=X":
            b.append((p, p + n))
            p += n
        elif op == "D":
            if b:
                b[-1] = (b[-1][0], p + n)
            p += n
        elif op == "N":
            p += n
    return b


def ov(a, b):
    return sum(max(0, min(x1, y1) - max(x0, y0)) for x0, x1 in a for y0, y1 in b)


def load_gtf(path, chrom):
    tx, gene = collections.defaultdict(list), collections.defaultdict(list)
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("#") or (chrom and not line.startswith(chrom + "\t")):
                continue
            c = line.split("\t")
            if c[2] != "exon":
                continue
            t = re.search(r'transcript_id "([^"]+)"', c[8]).group(1)
            g = re.search(r'gene_id "([^"]+)"', c[8]).group(1)
            e = (int(c[3]) - 1, int(c[4]))
            tx[t].append(e)
            gene[g].append(e)
    return tx, gene


def neighbours(u, alphabet="ACG"):
    for i in range(len(u)):
        for a in alphabet:
            if a == u[i]:
                continue
            v = u[:i] + a + u[i + 1:]
            yield v, 1
            for j in range(i + 1, len(u)):
                for b in alphabet:
                    if b != v[j]:
                        yield v[:j] + b + v[j + 1:], 2


def near_parents(group, umis):
    parents = {u for u, r in zip(umis, group) if int(r["reads"]) >= 3}
    singles = [u for u, r in zip(umis, group) if r["reads"] == "1"]
    c = collections.Counter()
    for u in singles:
        ds = {d for v, d in neighbours(u) if v in parents}
        if 1 in ds:
            c[1] += 1
        elif 2 in ds:
            c[2] += 1
    return len(singles), c


def pct(a, b):
    return 100.0 * a / b if b else float("nan")


def first(d, *names):
    for n in names:
        if os.path.exists(f"{d}/{n}"):
            return f"{d}/{n}"
    raise FileNotFoundError(f"{d}: none of {names}")


def sample(d, tx, gene):
    s = stats(f"{d}/filter.stats.tsv")
    i = lambda k: int(float(s[k]))
    total = i("total_records") - i("dropped_secondary_supplementary")  # primary records
    rows = list(csv.DictReader(open(f"{d}/molecules.tsv"), delimiter="\t"))
    view = subprocess.run(["samtools", "view", first(d, "molecules.bam", "chr22.molecules.bam")], capture_output=True, text=True, check=True)
    blocks = {}
    for line in view.stdout.splitlines():
        f = line.split("\t", 6)
        blocks[f[0]] = blocks_of(f[5], int(f[3]))
    fl = [r for r in rows if r["full_length"] == "Y"]
    cover = [ov(blocks[r["rep_read"]], tx[r["transcript"]]) / sum(e - s0 for s0, e in tx[r["transcript"]]) for r in fl]
    bad_gene = sum(
        1 for r in rows if r["genes"] != "." for g in r["genes"].split(",") if ov(blocks[r["rep_read"]], gene[g]) == 0
    )
    # densest 5' neighbourhood (20 bp buckets), undamaged UMIs only
    und = [r for r in rows if r["umi_damaged"] == "N" and len(r["umi"]) == 16]
    grp = collections.defaultdict(list)
    for r in und:
        grp[(r["strand"], int(r["anchor"]) // 20)].append(r)
    key, big = max(grp.items(), key=lambda kv: len(kv[1]))
    ns, obs = near_parents(big, [r["umi"] for r in big])
    pool = [r["umi"] for k, g in grp.items() if k != key for r in g]
    random.seed(1)
    _, null = near_parents(big, random.sample(pool, len(big))) if len(pool) >= len(big) else (0, collections.Counter())
    log = open(f"{d}/flumi.filter.log").read()
    version = open(f"{d}/flumi.version").read().split()[-1] if os.path.exists(f"{d}/flumi.version") else "?"
    wall = re.search(r"([\d.]+) real", log)
    mem = re.search(r"(\d+)\s+peak memory footprint", log)
    per_gene, per_gene_reads = collections.Counter(), collections.Counter()
    for r in rows:
        if r["genes"] != "." and "," not in r["genes"]:
            per_gene[r["genes"]] += 1
            per_gene_reads[r["genes"]] += int(r["reads"])
    return {
        "sample": os.path.basename(os.path.normpath(d)),
        "flumi": version,
        "primary reads": f"{total:,}",
        "RX present": f"{pct(total - i('dropped_no_RX'), total):.1f} %",
        "UMI block-damaged": f"{float(s['umi_block_damaged_pct']):.1f} %",
        "chimeras dropped": f"{pct(i('dropped_chimeric_umi_end'), total):.2f} %",
        "reads used": f"{i('reads_used'):,}",
        "molecules": f"{i('molecules_out'):,}",
        "PCR duplicates": f"{float(s['pcr_duplicate_pct']):.1f} %",
        "singleton molecules": f"{pct(i('singleton_molecules'), i('molecules_out')):.1f} %",
        "reads/molecule p90 / max": f"{s['cluster_size_p90']} / {s['cluster_size_max']}",
        "full length (GTF)": f"{float(s['full_length_pct']):.1f} %",
        "poly(A) anchor": f"{float(s['polya_anchor_pct']):.1f} %",
        "check: FL calls <80 % covered": f"{sum(c <= 0.8 for c in cover):,}",
        "check: listed genes w/o exon base": f"{bad_gene:,}",
        "densest locus molecules": f"{len(big):,}",
        "  singletons 1 sub from parent": f"{pct(obs[1], ns):.2f} % (null {pct(null[1], ns):.2f} %)",
        "  singletons 2 subs from parent": f"{pct(obs[2], ns):.2f} % (null {pct(null[2], ns):.2f} %)",
        "flumi wall / peak memory": f"{wall.group(1) if wall else '?'} s / {int(mem.group(1)) / 1e6:.0f} MB" if mem else "?",
        "_genes": per_gene,
        "_gene_reads": per_gene_reads,
    }


def spearman(x, y):
    def rank(v):
        order = sorted(range(len(v)), key=lambda k: v[k])
        r = [0.0] * len(v)
        i = 0
        while i < len(v):
            j = i
            while j + 1 < len(v) and v[order[j + 1]] == v[order[i]]:
                j += 1
            for k in range(i, j + 1):
                r[order[k]] = (i + j) / 2
            i = j + 1
        return r

    rx, ry = rank(x), rank(y)
    mx, my = sum(rx) / len(rx), sum(ry) / len(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = math.sqrt(sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry))
    return num / den


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--chrom", default="", help="restrict the GTF to one chromosome (faster)")
    ap.add_argument("--json", help="also write the table and comparisons as JSON (for the docs page)")
    ap.add_argument("--min-molecules", type=int, default=20, help="genes compared across samples")
    ap.add_argument("samples", nargs="+")
    a = ap.parse_args()
    tx, gene = load_gtf(a.gtf, a.chrom)
    res = [sample(d, tx, gene) for d in a.samples]
    keys = [k for k in res[0] if not k.startswith("_")]
    w = max(len(k) for k in keys)
    print("| metric | " + " | ".join(r["sample"] for r in res) + " |")
    print("|---|" + "---|" * len(res))
    for k in keys[1:]:
        print(f"| {k} | " + " | ".join(r[k] for r in res) + " |")
    comparisons = []
    if len(res) > 1:
        print("\nPer-gene molecule counts between samples (single-gene molecules; genes with >= "
              f"{a.min_molecules} in both): Spearman rho of molecule and of read counts, median |log2 ratio| of "
              "molecule counts after scaling to equal totals")
        for x, y in itertools.combinations(res, 2):
            gx, gy = x["_genes"], y["_genes"]
            genes = [g for g in gx if gx[g] >= a.min_molecules and gy.get(g, 0) >= a.min_molecules]
            if len(genes) < 3:
                continue
            sx, sy = sum(gx[g] for g in genes), sum(gy[g] for g in genes)
            lr = sorted(abs(math.log2((gx[g] / sx) / (gy[g] / sy))) for g in genes)
            rho = spearman([gx[g] for g in genes], [gy[g] for g in genes])
            rx, ry = x["_gene_reads"], y["_gene_reads"]
            rho_reads = spearman([rx[g] for g in genes], [ry[g] for g in genes])
            comparisons.append({"a": x["sample"], "b": y["sample"], "genes": len(genes), "rho": round(rho, 3),
                                "rho_reads": round(rho_reads, 3), "median_abs_log2": round(lr[len(lr) // 2], 3)})
            print(f"  {x['sample']} vs {y['sample']}: {len(genes)} genes, rho molecules {rho:.3f}, reads {rho_reads:.3f}, "
                  f"median |log2| {lr[len(lr) // 2]:.3f}")


    if a.json:
        with open(a.json, "w") as f:
            json.dump({"region": a.chrom or "genome", "metrics": keys[1:],
                       "samples": [{k: r[k] for k in keys} for r in res],
                       "comparisons": comparisons}, f, indent=1, ensure_ascii=False)


if __name__ == "__main__":
    sys.exit(main())
