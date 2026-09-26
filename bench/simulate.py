#!/usr/bin/env python3
"""Simulate aligned ONT PCR-cDNA (SQK-PCB114.24) reads with known molecules.

Every read's name is ``m<molecule>:<copy>`` (``:x`` appended for PCR chimeras),
so any deduplication output can be scored against the truth.

What is modelled, and why it matters for UMI deduplication:

* The UMI tag exactly as dorado builds it: the read front (adapter/barcode
  stand-in, SSP primer, UMI, GGG, insert) carries nanopore errors
  (substitutions, indels, homopolymer length errors); the SSP primer is found
  with edlib HW alignment (score >= 0.8), a 40-nt window starting 6 nt before
  the primer end is searched for ``TTTVVVVTTVVVVTTVVVVTTVVVVTTT`` with V
  wildcards, and the hit (score >= 0.8, i.e. <= 5 edits) is reported as RX.
* TS:A read orientation (dorado), ts:A on spliced reads (minimap2), pt:i
  poly(A) evidence, qs:f.
* PCR: a copy-number distribution per molecule and early-cycle polymerase
  errors that give a fraction of copies a variant UMI.
* Alignment geometry: TSS heterogeneity between molecules (sharp and broad
  promoters), few-bp 5' jitter between copies, short first exons that the
  aligner soft-clips in some copies, junction wobble, fabricated junctions,
  3' truncation of forward reads, reverse reads that stop before the UMI.
* PCR chimeras: the UMI end belongs to one molecule, the primary alignment
  to another gene (5' soft clip + SA), with a supplementary record.
* Antisense-overlapping genes and a few very highly expressed genes, where
  UMI collisions between distinct molecules actually happen.
"""
from __future__ import annotations

import argparse
import math
import random
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional, Tuple

import edlib
import pysam

SSP = "TTTCTGTTGGTGCTGATATTGCTTT"
UMI_PATTERN = "TTTVVVVTTVVVVTTVVVVTTVVVVTTT"
V_EQ = [("V", "A"), ("V", "C"), ("V", "G")]
OTHER = {b: [c for c in "ACGT" if c != b] for b in "ACGT"}

PROFILES = {
    # Per-base substitution, insertion, deletion, and homopolymer length error
    # per extra base of a run. "typical" reproduces the ~12 % of dorado RX tags
    # that are not 28 nt long reported for real SQK-PCB114.24 data; "good" and
    # "noisy" give ~7 % and ~27 %.
    "good": (0.002, 0.0005, 0.001, 0.004),
    "typical": (0.004, 0.001, 0.002, 0.0065),
    "noisy": (0.010, 0.003, 0.006, 0.015),
}


def mutate(seq: str, rng: random.Random, prof) -> str:
    sub, ins, dele, hp = prof
    out = []
    i, n = 0, len(seq)
    while i < n:
        j = i
        while j < n and seq[j] == seq[i]:
            j += 1
        base, run = seq[i], j - i
        if run >= 2 and rng.random() < hp * (run - 1):
            run += -1 if rng.random() < 0.6 else 1
        for _ in range(run):
            r = rng.random()
            if r < dele:
                pass
            elif r < dele + sub:
                out.append(rng.choice(OTHER[base]))
            else:
                out.append(base)
            if rng.random() < ins:
                out.append(rng.choice("ACGT"))
        i = j
    return "".join(out)


def dorado_rx(front: str) -> Optional[str]:
    """dorado's UMI extraction on the read front (forward orientation)."""
    region = front[:150]
    res = edlib.align(SSP, region, mode="HW", task="locations")
    if res["editDistance"] < 0 or 1.0 - res["editDistance"] / len(SSP) < 0.8:
        return None
    first = res["locations"][0][1] + 1
    a = first - 6
    if a < 0 or a + 40 >= len(front):
        return None
    window = front[a : a + 40]
    res = edlib.align(UMI_PATTERN, window, mode="HW", task="locations", additionalEqualities=V_EQ)
    if res["editDistance"] < 0 or 1.0 - res["editDistance"] / len(UMI_PATTERN) < 0.8:
        return None
    s, e = res["locations"][0]
    return window[s : e + 1]


def random_seq(rng: random.Random, n: int) -> str:
    return "".join(rng.choice("ACGT") for _ in range(n))


def umi_tag(blocks: List[str]) -> str:
    return "TTT" + "TT".join(blocks) + "TTT"


@dataclass
class Gene:
    gid: str
    tid: int
    strand: str
    exons: List[Tuple[int, int]]  # genomic, ascending, canonical transcript
    isoforms: List[List[int]]  # exon index subsets (ascending)
    broad_tss: bool
    expr: int = 0


@dataclass
class Molecule:
    mid: int
    gene: Gene
    iso: int
    blocks: List[str]
    exons: List[Tuple[int, int]]  # genomic, in transcript order, after TSS/polyA shifts
    tail: int


def transcript_exons(g: Gene, iso: int) -> List[Tuple[int, int]]:
    ex = [g.exons[i] for i in g.isoforms[iso]]
    return ex if g.strand == "+" else ex[::-1]


def make_genes(rng: random.Random, n_genes: int, contig_len: int, n_contigs: int, n_hot: int, hot_mols: int, mean_mols: float):
    genes: List[Gene] = []
    per = math.ceil(n_genes / n_contigs)
    for c in range(n_contigs):
        pos = 5000
        for k in range(per):
            if len(genes) >= n_genes:
                break
            strand = rng.choice("+-")
            n_ex = 1 if rng.random() < 0.1 else min(14, 2 + int(rng.expovariate(1 / 4.0)))
            lens = []
            for e in range(n_ex):
                if e == 0 and n_ex > 1 and rng.random() < 0.15:
                    lens.append(rng.randint(10, 30))  # short first exon (transcript order)
                elif e == n_ex - 1:
                    lens.append(rng.randint(300, 1500))
                else:
                    lens.append(rng.randint(60, 250))
            introns = [rng.randint(150, 4000) for _ in range(n_ex - 1)]
            if strand == "-":
                lens = lens[::-1]
                introns = introns[::-1]
            overlap_prev = genes and rng.random() < 0.08 and genes[-1].tid == c
            start = pos
            if overlap_prev:
                prev = genes[-1]
                start = max(prev.exons[0][0] + 200, prev.exons[-1][1] - 1500)
                strand_prev = prev.strand
                if strand == strand_prev:
                    strand = "-" if strand_prev == "+" else "+"
                    lens = lens[::-1]
                    introns = introns[::-1]
            exons = []
            p = start
            for i, L in enumerate(lens):
                exons.append((p, p + L))
                p += L
                if i < len(introns):
                    p += introns[i]
            if p + 5000 >= contig_len:
                break
            isoforms = [list(range(n_ex))]
            if n_ex >= 4:
                for _ in range(rng.randint(0, 2)):
                    skip = rng.randint(1, n_ex - 2)
                    iso = [i for i in range(n_ex) if i != skip]
                    if iso not in isoforms:
                        isoforms.append(iso)
            g = Gene(f"G{len(genes):05d}", c, strand, exons, isoforms, rng.random() < 0.2)
            g.expr = max(1, int(rng.lognormvariate(math.log(mean_mols), 1.2)))
            genes.append(g)
            pos = max(pos, p + rng.randint(2000, 15000))
    hot = rng.sample(range(len(genes)), min(n_hot, len(genes)))
    for h in hot:
        genes[h].expr = hot_mols
    return genes


def make_molecule(rng: random.Random, mid: int, g: Gene) -> Molecule:
    iso = 0 if len(g.isoforms) == 1 or rng.random() < 0.6 else rng.randrange(len(g.isoforms))
    ex = transcript_exons(g, iso)
    # TSS offset in transcript direction (positive = downstream).
    sd = 25.0 if g.broad_tss else 1.5
    d = int(round(rng.gauss(0, sd)))
    L0 = ex[0][1] - ex[0][0]
    d = max(-40, min(d, L0 - 6))
    s, e = ex[0]
    ex[0] = (s + d, e) if g.strand == "+" else (s, e - d)
    # poly(A) site jitter / internal priming at the 3' end.
    s, e = ex[-1]
    Ll = e - s
    cut = int(Ll * rng.uniform(0.3, 0.7)) if rng.random() < 0.03 else -int(round(rng.gauss(0, 2)))
    cut = max(-20, min(cut, Ll - 30))
    ex[-1] = (s, e - cut) if g.strand == "+" else (s + cut, e)
    blocks = ["".join(rng.choice("ACG") for _ in range(4)) for _ in range(4)]
    tail = max(10, int(rng.gauss(90, 30)))
    return Molecule(mid, g, iso, blocks, ex, tail)


def tx_len(ex) -> int:
    return sum(e - s for s, e in ex)


def blocks_for(m: Molecule, t_from: int, t_to: int) -> List[Tuple[int, int]]:
    """Genomic intervals (ascending) covering transcript coordinates [t_from, t_to)."""
    out = []
    t = 0
    strand = m.gene.strand
    for s, e in m.exons:
        L = e - s
        a, b = max(t_from, t), min(t_to, t + L)
        if a < b:
            if strand == "+":
                out.append((s + (a - t), s + (b - t)))
            else:
                out.append((e - (b - t), e - (a - t)))
        t += L
    out.sort()
    return out


def cigar_from(blocks, clip_left: int, clip_right: int) -> List[Tuple[int, int]]:
    cig = []
    if clip_left > 0:
        cig.append((4, clip_left))
    for i, (s, e) in enumerate(blocks):
        if i:
            cig.append((3, s - blocks[i - 1][1]))
        cig.append((0, e - s))
    if clip_right > 0:
        cig.append((4, clip_right))
    return cig


@dataclass
class Rec:
    tid: int
    pos: int
    name: str
    flag: int
    cigar: list
    tags: list
    mapq: int = 60


def clip_prob(first_exon_len: int) -> float:
    if first_exon_len < 12:
        return 0.9
    if first_exon_len < 20:
        return 0.5
    if first_exon_len < 30:
        return 0.2
    return 0.0


def perturb_blocks(rng: random.Random, blocks):
    """Junction wobble and the occasional fabricated junction."""
    blocks = [list(b) for b in blocks]
    for i in range(1, len(blocks)):
        if rng.random() < 0.02:
            d = rng.choice([-3, -2, -1, 1, 2, 3])
            if blocks[i - 1][1] + d > blocks[i - 1][0] + 1 and blocks[i][0] + d < blocks[i][1] - 1:
                blocks[i - 1][1] += d
                blocks[i][0] += d
    if rng.random() < 0.01:
        cands = [i for i, (s, e) in enumerate(blocks) if e - s >= 200]
        if cands:
            i = rng.choice(cands)
            s, e = blocks[i]
            cut = rng.randint(s + 60, e - 120)
            gap = rng.randint(20, 60)
            blocks[i] = [s, cut]
            blocks.insert(i + 1, [cut + gap, e])
    return [tuple(b) for b in blocks]


def read_records(rng, prof, m: Molecule, copy: int, pcr_blocks, other_mols, p_chim, seq_mode) -> List[Rec]:
    g = m.gene
    strand = g.strand
    L = tx_len(m.exons)
    forward = rng.random() < 0.55
    name = f"m{m.mid}:{copy}"
    classified = rng.random() > 0.03

    # Transcript span covered by the read.
    if forward:
        t_from = 0
        full = rng.random() < 0.7
        t_to = L if full else rng.randint(min(L, 150), L)
        has_umi_end = True
    else:
        t_to = L
        full = rng.random() < 0.7
        t_from = 0 if full else rng.randint(0, max(0, L - 150))
        has_umi_end = t_from == 0
    has_3p = t_to == L

    rx = None
    if classified and has_umi_end:
        front = random_seq(rng, 60) + SSP + umi_tag(pcr_blocks)[3:] + "GGG" + random_seq(rng, 40)
        # Only the primer/UMI stretch is mutated; the random flanks already
        # stand in for arbitrary sequence.
        cut = 60 - 5
        rx = dorado_rx(front[:cut] + mutate(front[cut:], rng, prof))

    # 5' geometry.
    clip5 = 3
    first_len = m.exons[0][1] - m.exons[0][0]
    if t_from == 0 and len(m.exons) > 1 and rng.random() < clip_prob(first_len):
        t_from = first_len
        clip5 += first_len
    elif t_from == 0 and rng.random() < 0.05 and first_len > 20:
        k = rng.randint(1, 8)
        t_from = k
        clip5 += k
    blocks = blocks_for(m, t_from, t_to)
    if t_from == 0 and rng.random() < 0.2:
        g_ext = rng.randint(1, 2)
        clip5 -= g_ext
        if strand == "+":
            blocks[0] = (blocks[0][0] - g_ext, blocks[0][1])
        else:
            blocks[-1] = (blocks[-1][0], blocks[-1][1] + g_ext)
    blocks = perturb_blocks(rng, blocks)
    clip3 = min(m.tail, 200) if has_3p else rng.randint(0, 5)

    chimera = rng.random() < p_chim and forward and rx is not None and other_mols
    recs = []
    if chimera:
        b = rng.choice(other_mols)
        bl = tx_len(b.exons)
        seg = min(bl, rng.randint(400, 1200))
        a_len = min(L, rng.randint(150, 400))
        b_blocks = blocks_for(b, bl - seg, bl)
        a_blocks = blocks_for(m, 0, a_len)
        q_len = a_len + seg + 3
        # Primary: the B part; the A part (with the UMI end) is soft-clipped.
        bs = b.gene.strand
        rev = bs == "-"
        c5, c3 = a_len + 3, min(b.tail, 200)
        cig = cigar_from(b_blocks, c5 if bs == "+" else c3, c3 if bs == "+" else c5)
        sa_strand = "+" if strand == "+" else "-"
        sa = f"{'chr%d' % (g.tid + 1)},{a_blocks[0][0] + 1},{sa_strand},{a_len}M{seg}S,60,0;"
        tags = [("RX", rx, "Z"), ("TS", "+", "A"), ("qs", round(rng.gauss(18, 3), 2), "f"), ("SA", sa, "Z")]
        if len(b_blocks) > 1:
            tags.append(("ts", "+", "A"))
        recs.append(Rec(b.gene.tid, b_blocks[0][0], name + ":x", 16 if rev else 0, cig, tags))
        sup_cig = cigar_from(a_blocks, 3 if strand == "+" else seg, seg if strand == "+" else 3)
        recs.append(Rec(g.tid, a_blocks[0][0], name + ":x", 2048 | (16 if strand == "-" else 0), sup_cig, [("qs", 18.0, "f")]))
        return recs

    # Clip placement in genomic order.
    left, right = (clip5, clip3) if strand == "+" else (clip3, clip5)
    cig = cigar_from(blocks, max(0, left), max(0, right))
    rev = (strand == "+") != forward
    tags = []
    if rx is not None:
        tags.append(("RX", rx, "Z"))
    if classified:
        tags.append(("TS", "+" if forward else "-", "A"))
    tags.append(("qs", round(min(40.0, max(5.0, rng.gauss(18, 3))), 2), "f"))
    if has_3p:
        r = rng.random()
        pt = max(1, m.tail + int(round(rng.gauss(0, 5)))) if r < 0.9 else (0 if r < 0.93 else -1)
    else:
        pt = -1
    tags.append(("pt", pt, "i"))
    if len(blocks) > 1:
        tags.append(("ts", "+" if forward else "-", "A"))
    mapq = 60 if rng.random() < 0.97 else rng.randint(0, 5)
    recs.append(Rec(g.tid, blocks[0][0], name, 16 if rev else 0, cig, tags, mapq))
    return recs


def pcr_variant(rng: random.Random, blocks: List[str]) -> List[str]:
    b = [list(x) for x in blocks]
    k, i = rng.randrange(4), rng.randrange(4)
    b[k][i] = rng.choice(OTHER[b[k][i]])
    return ["".join(x) for x in b]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True, type=Path, help="output directory")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--genes", type=int, default=2000)
    ap.add_argument("--mean-molecules", type=float, default=25.0, help="median molecules per gene")
    ap.add_argument("--hot-genes", type=int, default=2)
    ap.add_argument("--hot-molecules", type=int, default=20000)
    ap.add_argument("--depth", type=float, default=3.0, help="mean reads per molecule")
    ap.add_argument("--profile", choices=sorted(PROFILES), default="typical")
    ap.add_argument("--pcr-error", type=float, default=0.02, help="molecules with an early-cycle UMI error")
    ap.add_argument("--chimera", type=float, default=0.02, help="forward reads that are PCR chimeras")
    ap.add_argument("--contigs", type=int, default=4)
    ap.add_argument("--with-seq", action="store_true", help="write random SEQ/QUAL (realistic record size)")
    args = ap.parse_args(argv)

    rng = random.Random(args.seed)
    prof = PROFILES[args.profile]
    contig_len = 250_000_000
    genes = make_genes(rng, args.genes, contig_len, args.contigs, args.hot_genes, args.hot_molecules, args.mean_molecules)
    args.out.mkdir(parents=True, exist_ok=True)

    # GTF with the canonical (unshifted) transcripts.
    with open(args.out / "genes.gtf", "w") as fh:
        for g in genes:
            chrom = f"chr{g.tid + 1}"
            fh.write(f'{chrom}\tsim\tgene\t{g.exons[0][0] + 1}\t{g.exons[-1][1]}\t.\t{g.strand}\t.\tgene_id "{g.gid}";\n')
            for k, iso in enumerate(g.isoforms):
                tid = f"{g.gid}.{k}"
                for i in iso:
                    s, e = g.exons[i]
                    fh.write(
                        f'{chrom}\tsim\texon\t{s + 1}\t{e}\t.\t{g.strand}\t.\tgene_id "{g.gid}"; transcript_id "{tid}";\n'
                    )

    header = pysam.AlignmentHeader.from_dict(
        {
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [{"SN": f"chr{i + 1}", "LN": contig_len} for i in range(args.contigs)],
            "PG": [{"ID": "simulate", "PN": "simulate.py"}],
        }
    )
    records: List[Rec] = []
    truth = open(args.out / "truth.tsv", "w")
    truth.write("molecule\tgene\tisoform\tchrom\tstrand\tanchor\tumi\treads\tjunctions\n")
    mid = 0
    per_gene: List[List[Molecule]] = []
    for g in genes:
        mols = []
        for _ in range(g.expr):
            mols.append(make_molecule(rng, mid, g))
            mid += 1
        per_gene.append(mols)
    all_mols = [m for mols in per_gene for m in mols]
    pool = rng.sample(all_mols, min(5000, len(all_mols)))
    lam = max(0.0, args.depth - 1.0)
    for g, mols in zip(genes, per_gene):
        others = [x for x in pool if x.gene is not g]
        for m in mols:
            # Overdispersed copy number with mean `depth` (gamma-mixed).
            n = 1 + (int(rng.gammavariate(2.0, lam / 2.0) + 0.5) if lam > 0 else 0)
            variant = None
            frac = 0.0
            if n >= 2 and rng.random() < args.pcr_error:
                variant = pcr_variant(rng, m.blocks)
                frac = rng.choice([0.5, 0.25, 0.125])
            for c in range(n):
                blocks = variant if variant is not None and rng.random() < frac else m.blocks
                records.extend(read_records(rng, prof, m, c, blocks, others, args.chimera, args.with_seq))
            anchor = m.exons[0][0] + 1 if g.strand == "+" else m.exons[0][1]
            ex_sorted = sorted(m.exons)
            juncs = ",".join(f"{ex_sorted[i][1]}-{ex_sorted[i + 1][0]}" for i in range(len(ex_sorted) - 1))
            truth.write(
                f"{m.mid}\t{g.gid}\t{m.iso}\tchr{g.tid + 1}\t{g.strand}\t{anchor}\t{''.join(m.blocks)}\t{n}\t{juncs or '.'}\n"
            )
    truth.close()

    records.sort(key=lambda r: (r.tid, r.pos, r.name))
    path = args.out / "reads.bam"
    pool = random_seq(rng, 1 << 20) if args.with_seq else ""
    qual_pool = pysam.qualitystring_to_array("".join(chr(33 + rng.randint(8, 30)) for _ in range(1 << 16))) if args.with_seq else None
    with pysam.AlignmentFile(str(path), "wb", header=header, threads=4) as bf:
        for r in records:
            a = pysam.AlignedSegment(header)
            a.query_name = r.name
            a.flag = r.flag
            a.reference_id = r.tid
            a.reference_start = r.pos
            a.mapping_quality = r.mapq
            a.cigartuples = r.cigar
            if args.with_seq:
                qlen = sum(n for op, n in r.cigar if op in (0, 1, 4, 7, 8))
                off = rng.randrange(len(pool) - qlen)
                a.query_sequence = pool[off : off + qlen]
                qo = rng.randrange(len(qual_pool) - qlen)
                a.query_qualities = qual_pool[qo : qo + qlen]
            a.set_tags(r.tags)
            bf.write(a)
    pysam.index(str(path))
    print(f"{len(genes)} genes, {mid} molecules, {len(records)} records -> {path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
