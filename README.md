# flumi

**F**ul**l**-length **UMI** deduplication for Oxford Nanopore PCR-cDNA
(**SQK-PCB114.24**): collapse reads to **one representative read per original
RNA molecule**, using the UMI that strand switching attaches to every molecule
*before* PCR, and flag full-length molecules.

flumi reads a coordinate-sorted, splice-aligned BAM and writes a
coordinate-sorted BAM, so it drops into a pipeline between the aligner and the
isoform quantifier.

```bash
flumi --bam sample.sorted.bam --out sample.molecules.bam \
      --gtf gencode.gtf.gz --stats sample.stats.tsv --molecules sample.molecules.tsv \
      --threads 8
samtools index sample.molecules.bam
```

On simulated reads with known molecules, flumi's molecule count is within
0.5 % of the truth at typical nanopore error rates (1.5 % with noisy UMI
reads). It deduplicates 2.1 M reads in 19 s using 74 MB of memory
([benchmarks](#benchmarks)).

## Input

* Basecalled with **dorado ≥ 0.9.5**, `--kit-name SQK-PCB114-24`, trimming on.
  dorado moves the UMI into `RX:Z` and records the read orientation it inferred
  from the primers in `TS:A`. With `--estimate-poly-a` it also writes `pt:i`
  (optional, used for the poly(A) evidence).
* Aligned spliced and **coordinate-sorted**, with those tags carried over:
  `dorado aligner`, or
  `samtools fastq -T RX,TS,pt,qs calls.bam | minimap2 -ax splice -y ref.fa -`.
  The input is read twice, so it must be a file, not a pipe.
* A GTF is optional. Deduplication never uses it; it only labels genes and
  calls full-length transcripts.

## How it works

**1. Decoding `RX` the way dorado built it.** dorado finds the UMI by aligning
`TTTVVVVTTVVVVTTVVVVTTVVVVTTT` (V = A/C/G) to a 40-nt window behind the SSP
primer with edlib in infix mode. It accepts the hit with up to 5 edits and
stores the aligned read substring in `RX`. So `RX` carries substitutions and
indels anywhere in the tag. And because edlib keeps the co-optimal hit that
ends first, the trailing anchor can lose a `T` with no sequencing error at all.

flumi therefore re-aligns `RX` globally to the pattern and reads the four blocks
off the alignment, instead of cutting at fixed offsets. Costs follow nanopore
error rates. A T-run gaining or losing a T is cheapest, then a block base
substituted, gained or lost, then a substitution inside an anchor, which
silently discards a base. Outer anchor Ts missing from either end of `RX` are
free. Anchor errors of any kind therefore cost nothing.

A tag with a `T` or a wrong length inside a block cannot be a real UMI. It is
kept, marked **damaged**, and later attached to its molecule rather than
founding one or being thrown away. UMIs are compared by **blockwise edit
distance**: the anchors resynchronise the tag, so an indel never shifts bases
across blocks.

**2. Molecule identity = strand + 5′ end + UMI, no annotation.**
- Strand is dorado's `TS:A` (read orientation) combined with the alignment
  strand, falling back to minimap2's `ts:A`.
- dorado only finds the UMI if the read reaches the SSP primer, which sits on
  the molecule's 5′ end. So every UMI-bearing copy contains that end, and copies
  agree on it within `--pos-tol` (20 bp).
- One exception is handled explicitly: a short first exon (≤ `--bridge-exon-max`)
  that the aligner soft-clipped in some copies. A copy starting on the exon-2
  acceptor of another copy, with a soft clip long enough to hold the exon, is
  bridged to it.
- Reads are grouped into strand-specific loci of overlapping alignments. No gene
  annotation is involved, so reads outside annotated genes count, sense and
  antisense molecules never share a UMI namespace, and one molecule cannot be
  counted once per overlapping gene.
- A read whose UMI-side soft clip exceeds `--chimera-clip` and which has a
  supplementary alignment is a PCR chimera: its UMI belongs to another locus,
  so it is dropped.

**3. Clustering within a locus.**
- *Nodes* are reads with the same decoded UMI whose 5′ ends agree. The same UMI
  at a distant 5′ end is a different molecule.
- *Damaged nodes* are credited to their nearest compatible undamaged node, so
  the next step sees every read of a molecule.
- *Directional merging* (UMI-tools' rule, `n_hub ≥ 2·n_sub − 1`) joins
  positionally compatible undamaged nodes.
- Merges **beyond one edit are density-gated.** A UMI two or three edits from a
  hub is either a read with that many errors or a different molecule that
  happens to lie close. Which one is likelier depends on how many distinct UMIs
  share the 5′ neighbourhood. Among `m` random UMIs, about `m·N_d/3¹⁶` lie
  within distance `d` of any given UMI (`N_d` = 33, 513, 4993 for d = 1, 2, 3).
  A merge at d ≥ 2 is allowed only while that expectation stays within a
  budget:
  - `--collision-budget` 0.01 for undamaged UMIs, most of which are real
    molecules. This allows 2 edits in neighbourhoods of up to ~840 molecules.
  - `--damaged-collision-budget` 0.5 for damaged ones, which are certainly
    errors of *some* molecule. This allows 2 edits up to ~42k molecules and 3
    edits up to ~4.3k.

  In an ordinary gene this recovers reads with two or three UMI errors, which
  are common at nanopore error rates and would otherwise each count as an extra
  molecule. In a gene with tens of thousands of molecules at one TSS it falls
  back to single edits, where larger distances would mostly be coincidences.
- Neighbour search is exact and near-linear: UMIs within blockwise distance
  `d` share ≥ `4 − d` blocks verbatim, so candidates come from sorted keys over
  every `(4 − d)`-block subset.

**4. Representative read.** PCR copies share one true structure, so
disagreements between them are technical. Each junction is judged by majority
among the copies whose alignment spans it. A candidate is penalised for:
- each contradicted junction it carries,
- each majority junction it spans but lacks,
- each junction no other copy confirms, unless it is an annotated intron
  (turn this off with `--lenient-junctions`).

Ties go to the read that reached the poly(A) anchor (`pt ≥ 0`), then one
without a supplementary alignment, then the largest exonic footprint (not
genomic span), then basecall quality. `--rep longest` ranks by footprint alone.

**5. Full-length calls** (unless `--no-full-length`).
- `uF`/`uT`: against the annotation, with UMImap's criterion: more than
  `--fl-cov` of a transcript's exonic length covered, and at least
  `--fl-terminal` bases into both terminal exons.
- `uA`/`uP`: from dorado's poly(A) anchor, pooled over the molecule's copies.
  One copy that found the 3′ anchor proves the molecule reached it. `uP` is the
  median tail estimate over copies.

**6. Two streaming passes.**
- Pass 1 keeps no BAM records, only ~56 bytes per read of the loci still open.
  It decides each locus's molecules as soon as the sorted stream moves past it.
- Pass 2 re-reads the input and writes records in input order, so the output is
  coordinate-sorted without a reorder buffer.
- Memory is the largest open locus, plus ~70 bytes per molecule found, plus
  4 bytes per input record in `--mode mark`. That is 74 MB for the 2.1 M-read
  benchmark below; a run with 50 M molecules needs about 3.5 GB.

## Output

**`--mode filter`** (default) writes one record per molecule: the
representative. **`--mode mark`** writes every input record. Members of a
molecule carry `MI`, `uN` and `uS`, and every member except the representative
gets the duplicate flag `0x400`. Reads without a usable UMI, and unmapped,
secondary and supplementary records, pass through unchanged, apart from losing
tags and a duplicate flag left by an earlier flumi run.
`samtools view -F 0x400 -d MI` then selects one read per molecule.

| Tag | Written on | Meaning |
|---|---|---|
| `MI:Z` | molecule members | molecule id, numbered 1… in the order representatives appear |
| `uN:Z` | molecule members | the molecule's UMI (the 16 informative bases when undamaged) |
| `uS:i` | molecule members | reads in the molecule |
| `uJ:i` | representative | copies sharing the representative's exact junction chain |
| `uD:i` | representative | the representative's disagreements with the copies' junction consensus |
| `uA:A` | representative | poly(A) anchor pooled over copies: `Y`, `N`, or `?` (no `pt` tags) |
| `uP:i` | representative | median poly(A) tail length over copies |
| `uF:A`, `uT:Z` | representative | full-length against the annotation, and the transcript covered (GTF with exons) |
| `uG:Z` | representative | genes overlapping the representative on the molecule's strand (GTF) |

`--molecules` writes one row per molecule: id, chromosome, strand, 5′ end
(median over copies, 1-based), UMI, whether it is damaged, reads, `uJ`, `uD`,
the representative's name and span, poly(A) call and length, full-length call,
transcript, genes.

### Stats TSV

- **Drops:** `dropped_*` counts reads that cannot be assigned to a molecule, by
  reason: unmapped, secondary/supplementary, MAPQ, qs, no `RX`, undecodable `RX`,
  chimeric UMI end.
- **UMI decoding:** `umi_clean`, `umi_anchor_repaired` and `umi_block_damaged`
  split the reads that were used.
- **Strand source:** `strand_from_TS`, `strand_from_ts` and `strand_unknown`.
  Reads without a strand are grouped by UMI within their locus only (no 5′-end
  check); flumi warns when they are the majority.
- **Clustering:** `umi_nodes`, `umi_position_splits` (same UMI, different 5′
  end), `short_exon_bridges`, `directional_merges`, `density_gated_edges`,
  `damaged_to_valid`, `damaged_attached_beyond_1_edit`, `damaged_to_damaged`,
  `damaged_seeds`.
- **Result:** `molecules_out`, `pcr_duplicate_reads`, `pcr_duplicate_pct`,
  `singleton_molecules`, cluster-size percentiles.
- **Representatives:** `rep_switched_by_consensus`, `rep_disagreements`.
- **Full length:** `full_length_molecules`, `polya_anchor_molecules`,
  `polya_missing_tag_molecules`, `full_length_both_signals`.

## Benchmarks

There is no real dataset with known molecules, so accuracy is measured on
simulated reads (`bench/simulate.py`). The simulator:
- runs dorado's own UMI extraction (SSP located with edlib infix alignment, the
  40-nt window, the V-wildcard pattern, the 0.8 score thresholds, edlib's own
  tie rules) on read fronts carrying nanopore errors. The `typical` profile
  reproduces the ~12 % of `RX` tags that are not 28 nt long seen in real
  SQK-PCB114.24 data; `good` and `noisy` give ~7 % and ~27 %;
- models overdispersed PCR copy numbers with early-cycle UMI errors;
- varies TSS between molecules (sharp and broad promoters), with 5′ jitter
  between copies and short first exons clipped in some copies;
- adds junction wobble and fabricated junctions, 3′-truncated forward reads,
  reverse reads that stop before the UMI, and PCR chimeras with supplementary
  alignments;
- includes antisense-overlapping genes and one gene with 20,000 molecules at
  one TSS, where genuine UMI collisions occur.

Each read's name encodes its true molecule. `bench/evaluate.py` scores an
output against that truth:
- *split*: extra representatives of an already-represented molecule;
- *lost*: molecules with no representative;
- *chimeric*: representatives that are PCR chimeras;
- *gene err*: absolute per-gene count error, median and 90th percentile;
- *rep exact*: representatives whose junctions equal the true molecule's within
  the aligned span.

Four regimes, each with ~58 k molecules over 800 genes plus one gene at 20 k
molecules. "typical/noisy/good" is the UMI error profile, `dN` is the mean
number of reads per molecule:

| regime | true molecules | molecules out | split | lost | chimeric | count err | gene err median / p90 | rep exact |
|---|---|---|---|---|---|---|---|---|
| typical, d3 | 57,603 | 57,583 | 164 | 184 | 0 | −0.03 % | 0.0 / 0.0 % | 98.4 % |
| noisy, d3 | 57,652 | 58,517 | 1,204 | 339 | 0 | +1.50 % | 0.0 / 2.8 % | 98.4 % |
| typical, d8 | 59,065 | 59,312 | 384 | 137 | 0 | +0.42 % | 0.0 / 0.0 % | 99.7 % |
| good, d1.5 | 52,281 | 52,131 | 14 | 164 | 0 | −0.29 % | 0.0 / 0.0 % | 95.5 % |

Most *lost* molecules are genuine collisions in the 20 k-molecule gene. The
noisy regime is the hardest: about 27 % of its tags have indels, and several
percent carry two or three block errors.

Where a representative is wrong, it is almost always unavoidable. In
`typical_d3`, 862 of the 919 inexact representatives belong to a single-read
molecule or to a molecule with no correctly aligned copy.

**Large-scale run.** 2.1 M records (2.5 GB BAM with sequences and qualities),
673 k true molecules, three genes with ~29 k molecules each, `typical` profile,
4 threads:

| mode | molecules | split | lost | chimeric | count err | gene err median | rep exact | wall time | peak RSS |
|---|---|---|---|---|---|---|---|---|---|
| `--mode filter` | 672,599 | 1,109 | 1,627 | 0 | −0.08 % | 0.0 % | 98.0 % | 19.0 s | 74 MB |
| `--mode mark` | same | | | | | | | 31.1 s | 70 MB |

In the three 29 k-molecule genes flumi undercounts by 0.5–0.6 %. At that
density UMIs one edit apart are sometimes different molecules, and neither the
UMI nor the shared 5′ end can tell them apart. That is the capacity limit of a
16-base UMI.

Reproduce with:

```bash
pip install pysam edlib
cargo build --release
python3 bench/benchmark.py --work /tmp/flumi-bench          # the four regimes

python3 bench/simulate.py --out /tmp/flumi-scale --genes 12000 --mean-molecules 25 \
    --hot-genes 3 --hot-molecules 30000 --depth 3 --with-seq --seed 5
python3 bench/timed.py filter -- target/release/flumi --bam /tmp/flumi-scale/reads.bam \
    --gtf /tmp/flumi-scale/genes.gtf --out /tmp/flumi-scale/filter.bam --threads 4
python3 bench/evaluate.py --truth /tmp/flumi-scale/truth.tsv \
    --input /tmp/flumi-scale/reads.bam --run filter=/tmp/flumi-scale/filter.bam
```

## Install

Every [release](https://github.com/martinandclaude/flumi/releases) has prebuilt
binaries, each with a SHA-256 checksum:

| Archive | Runs on |
|---|---|
| `flumi-<tag>-x86_64-unknown-linux-gnu.tar.gz` | Linux x86-64 with glibc ≥ 2.17 (RHEL/CentOS 7 and later) |
| `flumi-<tag>-aarch64-unknown-linux-gnu.tar.gz` | Linux ARM64 with glibc ≥ 2.17 |
| `flumi-<tag>-aarch64-apple-darwin.tar.gz` | macOS 11+ on Apple silicon |
| `flumi-<tag>-x86_64-apple-darwin.tar.gz` | macOS 11+ on Intel |

```bash
tar xzf flumi-v1.0.0-x86_64-unknown-linux-gnu.tar.gz
./flumi-v1.0.0-x86_64-unknown-linux-gnu/flumi --help
```

## Build and test

```bash
cargo build --release          # binary at target/release/flumi
cargo install --path .         # or install it
cargo test                     # unit tests and end-to-end runs over generated BAMs
```

Building compiles htslib and libdeflate from source (through rust-htslib), so a
C compiler and CMake are needed. The minimum supported Rust is 1.88. Release
notes are in [CHANGELOG.md](CHANGELOG.md).

Pushing a `v*` tag builds the release binaries and publishes them
(`.github/workflows/release.yml`); for a tag that already exists, run the
Release workflow by hand with the tag as input.

## Limitations

- The accuracy figures come from simulation. The simulator follows dorado's
  code path and reproduces the reported `RX` length drift, but real libraries
  can differ in error profile and UMI synthesis bias. Where possible, check
  against spike-ins with known molecule counts.
- The density gate assumes UMIs are drawn uniformly. Biased UMI synthesis makes
  collisions likelier than it predicts; lower the budgets if you suspect that.
- Reads without `RX` (dorado did not reach or recognise the UMI) cannot be
  assigned. They are dropped by `--mode filter` and passed through unmarked by
  `--mode mark`.
- In `--mode mark`, secondary and supplementary records of duplicate reads are
  not flagged.

## License

MIT — see [LICENSE](LICENSE).
