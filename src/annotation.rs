//! Optional GTF annotation: gene labels (`uG`), full-length calls against
//! annotated transcripts (`uF`/`uT`) and the set of annotated introns used to
//! vouch for junctions seen in a single copy. Deduplication itself never
//! depends on it.

use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

use anyhow::{Context, Result};
use flate2::read::MultiGzDecoder;
use rustc_hash::{FxHashMap, FxHashSet};

use crate::features::Strand;

#[derive(Clone, Debug)]
struct Gene {
    start: u32,
    end: u32,
    strand: Strand,
    idx: u32,
}

#[derive(Default)]
struct ChromGenes {
    genes: Vec<Gene>,
    /// Running maximum of `end` over `genes[..=i]`.
    max_end: Vec<u32>,
}

#[derive(Clone, Debug)]
pub struct Transcript {
    pub id: String,
    pub start: u32,
    pub end: u32,
    /// Sorted half-open exons.
    pub exons: Box<[(u32, u32)]>,
    pub exonic_len: u32,
}

pub struct Annotation {
    by_tid: Vec<ChromGenes>,
    introns: Vec<FxHashSet<(u32, u32)>>,
    gene_ids: Vec<String>,
    /// Union of each gene's transcript exons; empty for a gene known only
    /// from a `gene` record, which then overlaps by its span.
    gene_exons: Vec<Box<[(u32, u32)]>>,
    transcripts: Vec<Transcript>,
    /// Transcripts of each gene, as a range into `tx_of_gene`.
    gene_tx: Vec<(u32, u32)>,
    tx_of_gene: Vec<u32>,
}

fn attr<'a>(attrs: &'a str, key: &str) -> Option<&'a str> {
    let mut rest = attrs;
    while let Some(i) = rest.find(key) {
        let after = &rest[i + key.len()..];
        let at_word_start = i == 0 || matches!(rest.as_bytes()[i - 1], b' ' | b';' | b'\t');
        if at_word_start {
            if let Some(v) = after.strip_prefix(" \"") {
                return v.find('"').map(|j| &v[..j]);
            }
        }
        rest = after;
    }
    None
}

/// Bases shared by two sorted lists of half-open intervals, each free of
/// internal overlaps.
fn overlap(a: &[(u32, u32)], b: &[(u32, u32)]) -> u64 {
    let (mut i, mut j, mut n) = (0, 0, 0u64);
    while i < a.len() && j < b.len() {
        let (s, e) = (a[i].0.max(b[j].0), a[i].1.min(b[j].1));
        if s < e {
            n += u64::from(e - s);
        }
        if a[i].1 < b[j].1 {
            i += 1;
        } else {
            j += 1;
        }
    }
    n
}

fn parse_strand(s: &str) -> Strand {
    match s {
        "+" => Strand::Plus,
        "-" => Strand::Minus,
        _ => Strand::Unknown,
    }
}

struct TxAcc {
    tid: u32,
    gene: u32,
    exons: Vec<(u32, u32)>,
}

impl Annotation {
    /// `names` maps BAM reference names to tids; other chromosomes are skipped.
    pub fn load(path: &Path, names: &FxHashMap<Vec<u8>, u32>, n_tids: usize) -> Result<Annotation> {
        let file = File::open(path).with_context(|| format!("opening GTF {}", path.display()))?;
        let raw: Box<dyn Read> = if path.extension().is_some_and(|e| e == "gz") {
            Box::new(MultiGzDecoder::new(file))
        } else {
            Box::new(file)
        };
        let mut reader = BufReader::with_capacity(1 << 20, raw);

        // Per chromosome: gene_id -> gene index. Spans from `gene` records win
        // over spans derived from transcripts/exons.
        let mut gene_idx: Vec<FxHashMap<String, u32>> = (0..n_tids).map(|_| FxHashMap::default()).collect();
        let mut gene_ids: Vec<String> = Vec::new();
        let mut spans: Vec<(u32, u32, u32, Strand, bool)> = Vec::new(); // tid, start, end, strand, from_gene_record
                                                                        // Per chromosome too: some GTFs reuse transcript ids across the X/Y
                                                                        // pseudoautosomal regions.
        let mut txs: Vec<FxHashMap<String, TxAcc>> = (0..n_tids).map(|_| FxHashMap::default()).collect();

        let mut line = String::new();
        loop {
            line.clear();
            if reader
                .read_line(&mut line)
                .with_context(|| format!("reading {}", path.display()))?
                == 0
            {
                break;
            }
            let text = line.trim_end_matches(['\n', '\r']);
            if text.is_empty() || text.starts_with('#') {
                continue;
            }
            let mut f = text.splitn(9, '\t');
            let (Some(chrom), Some(_), Some(feature), Some(s1), Some(e1), Some(_), Some(strand), Some(_), Some(attrs)) = (
                f.next(),
                f.next(),
                f.next(),
                f.next(),
                f.next(),
                f.next(),
                f.next(),
                f.next(),
                f.next(),
            ) else {
                continue;
            };
            if feature != "gene" && feature != "transcript" && feature != "exon" {
                continue;
            }
            let Some(&tid) = names.get(chrom.as_bytes()) else {
                continue;
            };
            let (Ok(s1), Ok(e1)) = (s1.parse::<u32>(), e1.parse::<u32>()) else {
                continue;
            };
            if s1 == 0 || e1 < s1 {
                continue;
            }
            let (start, end) = (s1 - 1, e1);
            let Some(gene_id) = attr(attrs, "gene_id") else {
                continue;
            };
            let strand = parse_strand(strand);

            let genes_here = &mut gene_idx[tid as usize];
            let g = match genes_here.get(gene_id) {
                Some(&g) => g,
                None => {
                    let g = gene_ids.len() as u32;
                    gene_ids.push(gene_id.to_string());
                    spans.push((tid, start, end, strand, false));
                    genes_here.insert(gene_id.to_string(), g);
                    g
                }
            };
            let span = &mut spans[g as usize];
            if feature == "gene" {
                if !span.4 {
                    *span = (tid, start, end, strand, true);
                } else {
                    span.1 = span.1.min(start);
                    span.2 = span.2.max(end);
                }
            } else if !span.4 {
                span.1 = span.1.min(start);
                span.2 = span.2.max(end);
            }
            if feature == "exon" {
                if let Some(tx_id) = attr(attrs, "transcript_id") {
                    let here = &mut txs[tid as usize];
                    let acc = match here.get_mut(tx_id) {
                        Some(acc) => acc,
                        None => here.entry(tx_id.to_string()).or_insert(TxAcc {
                            tid,
                            gene: g,
                            exons: Vec::new(),
                        }),
                    };
                    acc.exons.push((start, end));
                }
            }
        }

        let mut by_tid: Vec<ChromGenes> = (0..n_tids).map(|_| ChromGenes::default()).collect();
        for (idx, &(tid, start, end, strand, _)) in spans.iter().enumerate() {
            by_tid[tid as usize].genes.push(Gene {
                start,
                end,
                strand,
                idx: idx as u32,
            });
        }
        for cg in &mut by_tid {
            cg.genes.sort_unstable_by_key(|g| (g.start, g.end, g.idx));
            let mut m = 0;
            cg.max_end = cg
                .genes
                .iter()
                .map(|g| {
                    m = m.max(g.end);
                    m
                })
                .collect();
        }

        let mut introns: Vec<FxHashSet<(u32, u32)>> = (0..n_tids).map(|_| FxHashSet::default()).collect();
        let mut all_txs: Vec<(u32, String, TxAcc)> = txs
            .into_iter()
            .enumerate()
            .flat_map(|(tid, m)| m.into_iter().map(move |(id, acc)| (tid as u32, id, acc)))
            .collect();
        all_txs.sort_unstable_by(|a, b| (a.0, &a.1).cmp(&(b.0, &b.1)));
        let mut transcripts = Vec::with_capacity(all_txs.len());
        let mut per_gene: Vec<Vec<u32>> = vec![Vec::new(); gene_ids.len()];
        for (_, id, mut acc) in all_txs {
            acc.exons.sort_unstable();
            acc.exons.dedup();
            let exonic_len = acc.exons.iter().map(|&(s, e)| e - s).sum();
            for w in acc.exons.windows(2) {
                if w[0].1 < w[1].0 {
                    introns[acc.tid as usize].insert((w[0].1, w[1].0));
                }
            }
            per_gene[acc.gene as usize].push(transcripts.len() as u32);
            transcripts.push(Transcript {
                id,
                start: acc.exons[0].0,
                end: acc.exons.iter().map(|e| e.1).max().expect("non-empty"),
                exons: acc.exons.into_boxed_slice(),
                exonic_len,
            });
        }
        let mut gene_tx = Vec::with_capacity(per_gene.len());
        let mut tx_of_gene = Vec::new();
        let mut gene_exons = Vec::with_capacity(per_gene.len());
        for list in per_gene {
            let mut exons: Vec<(u32, u32)> = list
                .iter()
                .flat_map(|&t| transcripts[t as usize].exons.iter().copied())
                .collect();
            exons.sort_unstable();
            let mut merged: Vec<(u32, u32)> = Vec::with_capacity(exons.len());
            for (s, e) in exons {
                match merged.last_mut() {
                    Some(last) if s <= last.1 => last.1 = last.1.max(e),
                    _ => merged.push((s, e)),
                }
            }
            gene_exons.push(merged.into_boxed_slice());
            gene_tx.push((tx_of_gene.len() as u32, list.len() as u32));
            tx_of_gene.extend(list);
        }
        Ok(Annotation {
            by_tid,
            introns,
            gene_ids,
            gene_exons,
            transcripts,
            gene_tx,
            tx_of_gene,
        })
    }

    pub fn n_genes(&self) -> usize {
        self.gene_ids.len()
    }

    pub fn n_transcripts(&self) -> usize {
        self.transcripts.len()
    }

    pub fn gene_id(&self, idx: u32) -> &str {
        &self.gene_ids[idx as usize]
    }

    pub fn transcript(&self, idx: u32) -> &Transcript {
        &self.transcripts[idx as usize]
    }

    /// Genes on `strand` (any strand if unknown) whose exons share a base
    /// with the aligned `blocks`, sorted by index. A gene inside one of the
    /// read's introns does not count.
    pub fn overlapping_genes(&self, tid: u32, blocks: &[(u32, u32)], strand: Strand, out: &mut Vec<u32>) {
        out.clear();
        let (Some(cg), Some(&(start, _)), Some(&(_, end))) =
            (self.by_tid.get(tid as usize), blocks.first(), blocks.last())
        else {
            return;
        };
        let hi = cg.genes.partition_point(|g| g.start < end);
        for i in (0..hi).rev() {
            if cg.max_end[i] <= start {
                break;
            }
            let g = &cg.genes[i];
            let strand_ok = strand == Strand::Unknown || g.strand == Strand::Unknown || g.strand == strand;
            if g.end <= start || !strand_ok {
                continue;
            }
            let exons = &self.gene_exons[g.idx as usize];
            let hit = if exons.is_empty() {
                overlap(blocks, &[(g.start, g.end)]) > 0
            } else {
                overlap(blocks, exons) > 0
            };
            if hit {
                out.push(g.idx);
            }
        }
        out.sort_unstable();
    }

    pub fn is_intron(&self, tid: u32, junction: (u32, u32)) -> bool {
        self.introns.get(tid as usize).is_some_and(|s| s.contains(&junction))
    }

    /// Full-length criterion: the read's aligned `blocks` cover more than
    /// `cov` of a transcript's exonic bases and at least `term` bases of each
    /// of its terminal exons. Returns the longest qualifying transcript.
    pub fn full_length(&self, genes: &[u32], blocks: &[(u32, u32)], cov: f64, term: u32) -> Option<u32> {
        let term = u64::from(term);
        let mut best: Option<u32> = None;
        for &g in genes {
            let (first, n) = self.gene_tx[g as usize];
            for &t in &self.tx_of_gene[first as usize..(first + n) as usize] {
                let tx = &self.transcripts[t as usize];
                if overlap(blocks, &tx.exons) as f64 <= cov * f64::from(tx.exonic_len) {
                    continue;
                }
                let (first_exon, last_exon) = (tx.exons[0], tx.exons[tx.exons.len() - 1]);
                if overlap(blocks, &[first_exon]) < term || overlap(blocks, &[last_exon]) < term {
                    continue;
                }
                if best.is_none_or(|b| tx.exonic_len > self.transcripts[b as usize].exonic_len) {
                    best = Some(t);
                }
            }
        }
        best
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn load(text: &str) -> Annotation {
        let mut f = tempfile::Builder::new().suffix(".gtf").tempfile().unwrap();
        f.write_all(text.as_bytes()).unwrap();
        let mut names = FxHashMap::default();
        names.insert(b"chr1".to_vec(), 0);
        names.insert(b"chrY".to_vec(), 1);
        Annotation::load(f.path(), &names, 2).unwrap()
    }

    #[test]
    fn transcript_ids_reused_across_chromosomes_stay_apart() {
        let gtf = "\
chr1\tt\texon\t101\t200\t.\t+\t.\tgene_id \"P\"; transcript_id \"TP\";
chr1\tt\texon\t401\t500\t.\t+\t.\tgene_id \"P\"; transcript_id \"TP\";
chrY\tt\texon\t9001\t9100\t.\t+\t.\tgene_id \"P\"; transcript_id \"TP\";
";
        let a = load(gtf);
        assert_eq!(a.n_transcripts(), 2);
        assert_eq!(a.n_genes(), 2);
        let mut out = Vec::new();
        a.overlapping_genes(0, &[(150, 450)], Strand::Plus, &mut out);
        let t = a.full_length(&out, &[(100, 200), (400, 500)], 0.8, 25).unwrap();
        assert_eq!((a.transcript(t).start, a.transcript(t).end), (100, 500));
    }

    #[test]
    fn overlap_of_interval_lists() {
        assert_eq!(overlap(&[(0, 10), (20, 30)], &[(5, 25)]), 10);
        assert_eq!(overlap(&[(0, 10)], &[(10, 20)]), 0);
        assert_eq!(overlap(&[(0, 100)], &[(10, 20), (30, 40), (90, 120)]), 30);
        assert_eq!(overlap(&[], &[(0, 1)]), 0);
    }

    #[test]
    fn gene_inside_an_intron_is_not_the_reads_gene() {
        // Host H spliced 100-200 / 900-1000; N is a single-exon gene inside
        // H's intron on the same strand.
        let gtf = "\
chr1\tt\texon\t101\t200\t.\t+\t.\tgene_id \"H\"; transcript_id \"TH\";
chr1\tt\texon\t901\t1000\t.\t+\t.\tgene_id \"H\"; transcript_id \"TH\";
chr1\tt\texon\t401\t480\t.\t+\t.\tgene_id \"N\"; transcript_id \"TN\";
";
        let a = load(gtf);
        let read = [(100, 200), (900, 1000)];
        let mut out = Vec::new();
        a.overlapping_genes(0, &read, Strand::Plus, &mut out);
        assert_eq!(out.iter().map(|&g| a.gene_id(g)).collect::<Vec<_>>(), vec!["H"]);
        // Even offered N, the read is not full-length for it: no base on it.
        let all: Vec<u32> = (0..a.n_genes() as u32).collect();
        let t = a.full_length(&all, &read, 0.8, 25).unwrap();
        assert_eq!(a.transcript(t).id, "TH");
        assert_eq!(a.full_length(&all, &[(100, 200), (950, 1000)], 0.8, 25), None);
    }

    #[test]
    fn full_length_counts_bases_on_the_transcript() {
        let gtf = "\
chr1\tt\texon\t101\t200\t.\t+\t.\tgene_id \"G\"; transcript_id \"T\";
chr1\tt\texon\t401\t500\t.\t+\t.\tgene_id \"G\"; transcript_id \"T\";
";
        let a = load(gtf);
        // 170 of 200 exonic bases, 70 into the last exon.
        assert!(a.full_length(&[0], &[(100, 200), (430, 500)], 0.8, 25).is_some());
        // Enough coverage overall would need the last exon: 20 < 25 bases.
        assert_eq!(a.full_length(&[0], &[(100, 200), (480, 500)], 0.8, 25), None);
        // A retained intron still covers every exonic base.
        assert!(a.full_length(&[0], &[(100, 500)], 0.8, 25).is_some());
        // A long footprint elsewhere does not count.
        assert_eq!(a.full_length(&[0], &[(100, 130), (600, 2000)], 0.8, 25), None);
    }

    #[test]
    fn attr_matches_whole_keys_only() {
        let a = r#"gene_id "G1"; transcript_id "T1"; havana_gene_id "X";"#;
        assert_eq!(attr(a, "gene_id"), Some("G1"));
        assert_eq!(attr(a, "transcript_id"), Some("T1"));
        let b = r#"havana_gene_id "X"; gene_id "G2";"#;
        assert_eq!(attr(b, "gene_id"), Some("G2"));
    }

    #[test]
    fn genes_transcripts_and_introns() {
        let gtf = "\
chr1\tt\texon\t101\t200\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";
chr1\tt\texon\t401\t500\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";
chr1\tt\tgene\t450\t900\t.\t-\t.\tgene_id \"G2\";
chr2\tt\tgene\t1\t900\t.\t-\t.\tgene_id \"G3\";
";
        let a = load(gtf);
        assert_eq!(a.n_genes(), 2);
        assert_eq!(a.n_transcripts(), 1);
        let mut out = Vec::new();
        a.overlapping_genes(0, &[(150, 460)], Strand::Plus, &mut out);
        assert_eq!(out.iter().map(|&g| a.gene_id(g)).collect::<Vec<_>>(), vec!["G1"]);
        a.overlapping_genes(0, &[(150, 460)], Strand::Unknown, &mut out);
        assert_eq!(out.len(), 2);
        // G2 has no exons, so its span counts; G1's intron alone does not.
        a.overlapping_genes(0, &[(250, 300)], Strand::Unknown, &mut out);
        assert!(out.is_empty());
        a.overlapping_genes(0, &[(600, 700)], Strand::Unknown, &mut out);
        assert_eq!(out.iter().map(|&g| a.gene_id(g)).collect::<Vec<_>>(), vec!["G2"]);
        assert!(a.is_intron(0, (200, 400)));
        assert!(!a.is_intron(0, (200, 401)));
        let t = a.full_length(&[0], &[(100, 200), (400, 500)], 0.8, 25).unwrap();
        assert_eq!(a.transcript(t).id, "T1");
        assert_eq!(a.full_length(&[0], &[(290, 500)], 0.8, 25), None);
    }
}
