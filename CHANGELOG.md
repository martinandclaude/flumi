# Changelog

## 1.0.1 (unreleased)

- **Gene labels and full-length calls use the representative's aligned
  blocks, not its genomic span.** Before, a gene lying inside one of the read's
  introns on the same strand was listed in `uG`, and the molecule could be
  called full-length (`uF`/`uT`) against it without a single base on it. On
  real HG002 PCB114-24 data (chr22, one barcode) 27 % of full-length calls were
  spurious in this way, and 30 % of molecules listed more than five genes; the
  full-length rate falls from 96.1 % to 88.6 %. Deduplication never used the
  annotation, so molecules, representatives and every deduplication statistic
  are unchanged.

## 1.0.0 (2026-09-26)

First release. flumi collapses ONT PCR-cDNA (SQK-PCB114.24) reads to one
representative read per original RNA molecule, using the pre-PCR UMI that
dorado reports in `RX:Z`.

- **UMI decoding follows dorado's `RX` extraction.** `RX` is aligned to the UMI
  pattern the way dorado built it, so substitutions and indels anywhere in the
  tag are handled. Block-damaged UMIs are kept and attached to their molecule
  instead of being discarded.
- **A molecule is identified by strand, 5′ end and UMI**, within
  strand-specific loci. Strand comes from dorado's `TS:A`, with minimap2's
  `ts:A` as the fallback. No gene annotation is needed.
- **UMI merging.** Damaged UMIs are credited to their source before UMI-tools'
  directional merge, and UMIs two or three edits apart are merged where the
  local molecule density makes a coincidence unlikely.
- **PCR chimeras** are dropped: reads whose UMI end is soft-clipped with a
  supplementary alignment.
- **Representatives** are chosen by majority-vote junction consensus among the
  copies, then poly(A) completeness, then exonic footprint.
- **Output**: `--mode filter` (one read per molecule) or `--mode mark` (every
  record kept, `MI` tags, duplicate flag on non-representatives), an optional
  `--molecules` table, and output verified to be a complete BGZF file. An
  optional GTF adds gene labels and full-length calls.
- Two streaming passes over the input; 2.1 M reads in 19 s and 74 MB on the
  simulated benchmark in the README.
- Prebuilt binaries for Linux (x86-64, ARM64; glibc ≥ 2.17) and macOS (Apple
  silicon, Intel). The minimum supported Rust is 1.88.
