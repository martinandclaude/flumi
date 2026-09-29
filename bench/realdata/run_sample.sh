#!/usr/bin/env bash
# flumi on the public ONT dataset UHRR_HG002_2026.06, HG002 SQK-PCB114-24.
#
# Those reads were basecalled with dorado 2.0.0 --kit-name SQK-PCB114-24
# --estimate-poly-a --no-trim, so primers and UMI are still in the sequence
# and RX/TS are missing. The wf-transcriptomes alignment of each barcode is
# public and indexed; this script takes its primary reads (one chromosome, or
# all of them), adds RX/TS with dorado trim, realigns and runs flumi in filter
# and mark mode.
#
#   run_sample.sh <flowcell 1|2> <barcode 01..04> [chromosome]
#
# With a chromosome, reads are streamed for that chromosome only and realigned
# to it alone (supplementary alignments elsewhere are then not seen). Without
# one, the whole barcode (~36-44 GB) is downloaded and aligned genome-wide.
#
# Environment: REF (GRCh38 no-alt analysis set FASTA), GTF (Ensembl GTF with
# chr names), DORADO, FLUMI, MINIMAP2 (default: on PATH), THREADS (default 4),
# OUT (default ./fc<flowcell>_bc<barcode>[_<chromosome>]).
set -euo pipefail
FC=$1 BC=$2 CHROM=${3:-}
: "${REF:?set REF to the GRCh38 FASTA}" "${GTF:?set GTF to an Ensembl GTF with chr names}"
DORADO=${DORADO:-dorado} FLUMI=${FLUMI:-flumi} MINIMAP2=${MINIMAP2:-minimap2} THREADS=${THREADS:-4}
OUT=${OUT:-fc${FC}_bc${BC}${CHROM:+_$CHROM}}
URL=https://ont-open-data.s3.amazonaws.com/UHRR_HG002_2026.06/analysis/cDNA/HAC/HG002/PCB114_24/HG002_PCB114_24_PolyA_${FC}/samples/barcode${BC}/alignment/reads.bam
T=(/usr/bin/time -v)
[[ $(uname) == Darwin ]] && T=(/usr/bin/time -l)
mkdir -p "$OUT"
cd "$OUT"

if [[ ! -s primary.bam ]]; then
  # samtools caches the remote index in the working directory.
  samtools view -b -F 0x900 -@ 2 -o primary.bam.part "$URL" ${CHROM:+"$CHROM"}
  rm -f reads.bam.bai
  mv primary.bam.part primary.bam
fi
if [[ -n $CHROM ]]; then
  samtools faidx "$REF" "$CHROM" > target.fa
  TARGET=target.fa
else
  TARGET=$REF
fi

"${T[@]}" "$DORADO" trim --sequencing-kit SQK-PCB114-24 primary.bam > trim.bam 2> trim.log
samtools fastq -T RX,TS,pt,qs,BC trim.bam 2> /dev/null > trim.fq
"${T[@]}" "$MINIMAP2" -ax splice -y -t "$THREADS" "$TARGET" trim.fq 2> mm2.log \
  | samtools sort -@ 2 -m 1G -o sorted.bam -
samtools index sorted.bam
rm trim.fq

"${T[@]}" "$FLUMI" --bam sorted.bam --out molecules.bam --gtf "$GTF" \
    --stats filter.stats.tsv --molecules molecules.tsv --threads "$THREADS" 2> flumi.filter.log
"${T[@]}" "$FLUMI" --bam sorted.bam --out marked.bam --mode mark --gtf "$GTF" \
    --stats mark.stats.tsv --threads "$THREADS" 2> flumi.mark.log
samtools quickcheck molecules.bam marked.bam
"$FLUMI" --version > flumi.version
echo "$OUT done"
