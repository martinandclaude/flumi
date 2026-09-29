//! The two passes over a coordinate-sorted BAM.
//!
//! Pass 1 keeps no BAM records: it streams the input, holds a compact entry
//! per read of each open bundle, and when a bundle closes decides its
//! molecules and representatives. Pass 2 streams the input again and writes
//! either the representatives (`filter`) or every record, tagged and with the
//! duplicate flag on non-representatives (`mark`). Output order is input
//! order, so it is coordinate-sorted without a reorder buffer. Memory is the
//! largest open bundle plus a compact record per molecule (and, in `mark`
//! mode, 4 bytes per input record); no BAM record is ever held.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use std::time::Instant;

use anyhow::{bail, Context, Result};
use log::{info, warn};
use rust_htslib::bam::{self, record::Aux, Read};
use rust_htslib::tpool::ThreadPool;
use rustc_hash::FxHashMap;

use crate::annotation::Annotation;
use crate::bundle::{Bundle, BundleRead};
use crate::cluster::{self, Params};
use crate::consensus::{self, RepParams};
use crate::features::{self, aligned_blocks, Geometry, Strand, StrandSource, NO_POS, PT_NO_ANCHOR};
use crate::stats::Stats;
use crate::umi;

const NONE: u32 = u32::MAX;
const FLAG_UNMAPPED: u16 = 0x4;
const FLAG_REVERSE: u16 = 0x10;
const FLAG_SECONDARY: u16 = 0x100;
const FLAG_DUPLICATE: u16 = 0x400;
const FLAG_SUPPLEMENTARY: u16 = 0x800;

/// Tags this tool writes; stale copies are removed first.
const OUR_TAGS: [&[u8; 2]; 10] = [b"MI", b"uN", b"uS", b"uJ", b"uA", b"uP", b"uF", b"uT", b"uG", b"uD"];

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Mode {
    /// Write one representative record per molecule.
    Filter,
    /// Write every record; tag molecule members and flag non-representatives.
    Mark,
}

#[derive(Clone, Debug)]
pub struct Config {
    pub input: PathBuf,
    pub output: PathBuf,
    pub gtf: Option<PathBuf>,
    pub stats: Option<PathBuf>,
    pub molecules: Option<PathBuf>,
    pub mode: Mode,
    pub umi_dist: u32,
    pub damaged_umi_dist: u32,
    pub valid_budget: f64,
    pub damaged_budget: f64,
    pub max_umi_edits: u8,
    pub pos_tol: u32,
    pub use_position: bool,
    pub use_strand: bool,
    pub bridge_exon_max: u32,
    pub chimera_clip: u32,
    pub min_mapq: u8,
    pub min_qs: f32,
    pub rep: RepParams,
    pub full_length: bool,
    pub fl_cov: f64,
    pub fl_terminal: u32,
    pub threads: usize,
    pub command_line: String,
}

struct MoleculeRecord {
    rep_ordinal: u64,
    umi: u64,
    n_reads: u32,
    same_chain: u32,
    disagreements: u32,
    /// b'Y', b'N' or b'?' for pooled poly(A)-anchor evidence.
    polya: u8,
    /// Median tail length over copies that estimated one; 0 if none.
    polya_len: u32,
    /// b'Y'/b'N' against the annotation, 0 when not assessed.
    fl: u8,
    fl_tx: u32,
    genes: u32,
    tid: u32,
    strand: Strand,
    anchor: u32,
}

struct Plan {
    molecules: Vec<MoleculeRecord>,
    mi: Vec<u32>,
    /// Mark mode: molecule of each input record, `NONE` if unassigned.
    assign: Vec<u32>,
    /// Filter mode: (representative ordinal, molecule), sorted.
    reps: Vec<(u64, u32)>,
    gene_lists: Vec<String>,
    total_records: u64,
}

struct Pass1<'a> {
    cfg: &'a Config,
    annotation: Option<&'a Annotation>,
    params: Params,
    params_unstranded: Params,
    bundles: [Option<Bundle>; 3],
    cur_tid: i64,
    last_pos: u32,
    geom: Geometry,
    gene_buf: Vec<u32>,
    block_buf: Vec<(u32, u32)>,
    gene_list_index: FxHashMap<Vec<u32>, u32>,
    plan: Plan,
    stats: Stats,
}

fn median(v: &mut [u32]) -> u32 {
    if v.is_empty() {
        return 0;
    }
    v.sort_unstable();
    let n = v.len();
    if n % 2 == 1 {
        v[n / 2]
    } else {
        ((u64::from(v[n / 2 - 1]) + u64::from(v[n / 2])) / 2) as u32
    }
}

impl<'a> Pass1<'a> {
    fn new(cfg: &'a Config, annotation: Option<&'a Annotation>) -> Self {
        let params = Params {
            umi_dist: cfg.umi_dist,
            damaged_dist: cfg.damaged_umi_dist,
            pos_tol: cfg.pos_tol,
            use_position: cfg.use_position,
            valid_budget: cfg.valid_budget,
            damaged_budget: cfg.damaged_budget,
        };
        let params_unstranded = Params {
            use_position: false,
            ..params.clone()
        };
        Pass1 {
            cfg,
            annotation,
            params,
            params_unstranded,
            bundles: [None, None, None],
            cur_tid: -1,
            last_pos: 0,
            geom: Geometry::default(),
            gene_buf: Vec::new(),
            block_buf: Vec::new(),
            gene_list_index: FxHashMap::default(),
            plan: Plan {
                molecules: Vec::new(),
                mi: Vec::new(),
                assign: Vec::new(),
                reps: Vec::new(),
                gene_lists: Vec::new(),
                total_records: 0,
            },
            stats: Stats::default(),
        }
    }

    fn add(&mut self, rec: &bam::Record, ordinal: u64) -> Result<()> {
        let cfg = self.cfg;
        let st = &mut self.stats;
        st.total_records += 1;
        if cfg.mode == Mode::Mark {
            self.plan.assign.push(NONE);
        }
        let flags = rec.flags();
        if flags & FLAG_UNMAPPED != 0 || rec.tid() < 0 {
            st.unmapped += 1;
            return Ok(());
        }
        if flags & (FLAG_SECONDARY | FLAG_SUPPLEMENTARY) != 0 {
            st.secondary_supplementary += 1;
            return Ok(());
        }
        if rec.mapq() < cfg.min_mapq {
            st.low_mapq += 1;
            return Ok(());
        }
        let tags = features::scan_tags(rec, cfg.max_umi_edits);
        if tags.qs.unwrap_or(0.0) < cfg.min_qs {
            st.low_qs += 1;
            return Ok(());
        }
        let umi = match tags.umi {
            None => {
                st.no_umi += 1;
                return Ok(());
            }
            Some(Err(_)) => {
                st.umi_undecodable += 1;
                return Ok(());
            }
            Some(Ok(u)) => u,
        };

        let reverse = flags & FLAG_REVERSE != 0;
        let resolve = |c: Option<u8>| c.map_or(Strand::Unknown, |c| features::genomic_strand(c, reverse));
        let (strand, source) = if !cfg.use_strand {
            (Strand::Unknown, StrandSource::None)
        } else if resolve(tags.dorado_ts) != Strand::Unknown {
            (resolve(tags.dorado_ts), StrandSource::Dorado)
        } else if resolve(tags.minimap_ts) != Strand::Unknown {
            (resolve(tags.minimap_ts), StrandSource::Minimap2)
        } else {
            (Strand::Unknown, StrandSource::None)
        };
        let pos = rec.pos();
        if !(0..i64::from(u32::MAX / 2)).contains(&pos) {
            bail!("alignment position {pos} out of range");
        }
        self.geom.parse(pos as u32, rec.raw_cigar());
        let clip5 = self.geom.clip5(strand);
        if tags.has_sa && strand != Strand::Unknown && clip5 > cfg.chimera_clip {
            st.chimeric_umi_end += 1;
            return Ok(());
        }

        if umi.edits == 0 {
            st.umi_clean += 1;
        } else if umi.is_damaged() {
            st.umi_damaged += 1;
        } else {
            st.umi_anchor_repaired += 1;
        }
        match (strand, source) {
            (Strand::Unknown, _) => st.strand_unknown += 1,
            (_, StrandSource::Dorado) => st.strand_dorado += 1,
            (_, StrandSource::Minimap2) => st.strand_minimap2 += 1,
            (_, StrandSource::None) => st.strand_unknown += 1,
        }
        st.reads_used += 1;

        let tid = rec.tid() as u32;
        let start = self.geom.start;
        if i64::from(tid) != self.cur_tid {
            if i64::from(tid) < self.cur_tid {
                bail!(
                    "input BAM is not coordinate-sorted (reference {tid} after {})",
                    self.cur_tid
                );
            }
            self.close_all();
            self.cur_tid = i64::from(tid);
        } else if start < self.last_pos {
            bail!(
                "input BAM is not coordinate-sorted (position {start} after {})",
                self.last_pos
            );
        }
        self.last_pos = start;
        for s in 0..3 {
            if self.bundles[s]
                .as_ref()
                .is_some_and(|b| start >= b.max_end.saturating_add(cfg.pos_tol))
            {
                self.close(s);
            }
        }

        let geom = &self.geom;
        let (alt_anchor, alt_exon) = geom
            .alt_anchor(strand, cfg.bridge_exon_max)
            .map_or((NO_POS, 0), |(p, e)| (p, e.min(u32::from(u16::MAX)) as u16));
        let bundle = self.bundles[strand as usize].get_or_insert_with(|| Bundle::new(tid, strand));
        let chain = bundle.chains.intern(&geom.junctions);
        bundle.reads.push(BundleRead {
            ordinal,
            umi: umi.key,
            start: geom.start,
            end: geom.end,
            anchor: geom.anchor(strand),
            alt_anchor,
            alt_exon,
            clip5: clip5.min(u32::from(u16::MAX)) as u16,
            chain,
            exonic_len: geom.exonic_len,
            qs: tags.qs.unwrap_or(0.0),
            pt: tags.pt,
            has_sa: tags.has_sa,
        });
        bundle.max_end = bundle.max_end.max(geom.end);
        Ok(())
    }

    fn close_all(&mut self) {
        for s in 0..3 {
            self.close(s);
        }
    }

    fn gene_list(&mut self, genes: &[u32]) -> u32 {
        let Some(ann) = self.annotation else { return NONE };
        if genes.is_empty() {
            return NONE;
        }
        if let Some(&i) = self.gene_list_index.get(genes) {
            return i;
        }
        let i = self.plan.gene_lists.len() as u32;
        let label = genes.iter().map(|&g| ann.gene_id(g)).collect::<Vec<_>>().join(",");
        self.plan.gene_lists.push(label);
        self.gene_list_index.insert(genes.to_vec(), i);
        i
    }

    fn close(&mut self, s: usize) {
        let Some(bundle) = self.bundles[s].take() else { return };
        self.stats.bundles += 1;
        self.stats.max_bundle_reads = self.stats.max_bundle_reads.max(bundle.reads.len() as u64);
        let params = if bundle.strand == Strand::Unknown {
            &self.params_unstranded
        } else {
            &self.params
        };
        let molecules = cluster::cluster(&bundle.reads, params, &mut self.stats.cluster);
        let ann = self.annotation;
        let tid = bundle.tid;
        let annotated = |j: (u32, u32)| ann.is_some_and(|a| a.is_intron(tid, j));
        let mut lens = Vec::new();
        let mut anchors = Vec::new();
        for m in molecules {
            let choice = consensus::choose(&bundle.reads, &m.reads, &bundle.chains, &annotated, &self.cfg.rep);
            let rep = &bundle.reads[choice.rep as usize];

            lens.clear();
            let (mut anchor_found, mut anchor_missed) = (false, false);
            for &i in &m.reads {
                let pt = bundle.reads[i as usize].pt;
                if pt >= 0 {
                    anchor_found = true;
                    if pt > 0 {
                        lens.push(pt as u32);
                    }
                } else if pt == PT_NO_ANCHOR {
                    anchor_missed = true;
                }
            }
            // The molecule's 5' end: the median over its copies.
            anchors.clear();
            anchors.extend(m.reads.iter().map(|&i| bundle.reads[i as usize].anchor));
            let anchor = median(&mut anchors);

            let polya = if anchor_found {
                b'Y'
            } else if anchor_missed {
                b'N'
            } else {
                b'?'
            };

            let (mut fl, mut fl_tx, mut genes) = (0u8, NONE, NONE);
            if let Some(a) = ann {
                let mut buf = std::mem::take(&mut self.gene_buf);
                aligned_blocks(rep.start, rep.end, bundle.chains.get(rep.chain), &mut self.block_buf);
                a.overlapping_genes(tid, &self.block_buf, bundle.strand, &mut buf);
                genes = self.gene_list(&buf);
                if self.cfg.full_length && a.n_transcripts() > 0 {
                    let t = a.full_length(&buf, &self.block_buf, self.cfg.fl_cov, self.cfg.fl_terminal);
                    fl = if t.is_some() { b'Y' } else { b'N' };
                    fl_tx = t.unwrap_or(NONE);
                }
                self.gene_buf = buf;
            }

            let st = &mut self.stats;
            st.molecules += 1;
            *st.size_hist.entry(m.reads.len() as u32).or_default() += 1;
            st.rep_switched += u64::from(choice.switched);
            st.rep_disagreements += u64::from(choice.disagreements);
            if fl == b'Y' {
                st.full_length_annot += 1;
            }
            if self.cfg.full_length {
                match polya {
                    b'Y' => {
                        st.polya_anchor += 1;
                        if fl == b'Y' {
                            st.full_length_both += 1;
                        }
                    }
                    b'?' => st.polya_missing += 1,
                    _ => {}
                }
            }

            let idx = self.plan.molecules.len() as u32;
            self.plan.molecules.push(MoleculeRecord {
                rep_ordinal: rep.ordinal,
                umi: m.umi,
                n_reads: m.reads.len() as u32,
                same_chain: choice.same_chain,
                disagreements: choice.disagreements,
                polya,
                polya_len: median(&mut lens),
                fl,
                fl_tx,
                genes,
                tid,
                strand: bundle.strand,
                anchor,
            });
            match self.cfg.mode {
                Mode::Mark => {
                    for &i in &m.reads {
                        self.plan.assign[bundle.reads[i as usize].ordinal as usize] = idx;
                    }
                }
                Mode::Filter => self.plan.reps.push((rep.ordinal, idx)),
            }
        }
    }

    fn finish(mut self) -> (Plan, Stats) {
        self.close_all();
        let mut plan = self.plan;
        plan.total_records = self.stats.total_records;
        plan.reps.sort_unstable();
        let mut order: Vec<u32> = (0..plan.molecules.len() as u32).collect();
        order.sort_unstable_by_key(|&m| plan.molecules[m as usize].rep_ordinal);
        plan.mi = vec![0; plan.molecules.len()];
        for (rank, &m) in order.iter().enumerate() {
            plan.mi[m as usize] = rank as u32 + 1;
        }
        (plan, self.stats)
    }
}

const BGZF_EOF: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00, 0x1b, 0x00, 0x03,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];

fn check_bgzf_eof(path: &std::path::Path) -> Result<()> {
    use std::io::{Read as _, Seek, SeekFrom};
    let mut f = File::open(path).with_context(|| format!("reopening {}", path.display()))?;
    let mut tail = [0u8; 28];
    let complete = f.seek(SeekFrom::End(-28)).is_ok() && f.read_exact(&mut tail).is_ok() && tail == BGZF_EOF;
    if !complete {
        bail!(
            "{} is truncated (no BGZF end-of-file block): the final write failed",
            path.display()
        );
    }
    Ok(())
}

fn check_sorted(header: &bam::HeaderView) -> Result<()> {
    let text = String::from_utf8_lossy(header.as_bytes());
    let so = text
        .lines()
        .find(|l| l.starts_with("@HD"))
        .and_then(|l| l.split('\t').find_map(|f| f.strip_prefix("SO:")))
        .unwrap_or("unknown");
    if so != "coordinate" {
        bail!("input BAM must be coordinate-sorted (@HD SO:coordinate), found SO:{so}; run `samtools sort` first");
    }
    Ok(())
}

fn pg_record(header: &bam::HeaderView) -> (String, Option<String>) {
    let text = String::from_utf8_lossy(header.as_bytes());
    let ids: Vec<String> = text
        .lines()
        .filter(|l| l.starts_with("@PG"))
        .filter_map(|l| l.split('\t').find_map(|f| f.strip_prefix("ID:")).map(str::to_string))
        .collect();
    let mut id = "flumi".to_string();
    let mut k = 1;
    while ids.contains(&id) {
        id = format!("flumi.{k}");
        k += 1;
    }
    (id, ids.last().cloned())
}

fn set_tag(rec: &mut bam::Record, tag: &[u8], value: Aux<'_>) -> Result<()> {
    rec.push_aux_unchecked(tag, value)
        .with_context(|| format!("writing tag {}", String::from_utf8_lossy(tag)))
}

/// Remove tags a previous flumi run wrote; true if there were any. One scan
/// in the common case where there are none.
fn strip_stale_tags(rec: &mut bam::Record) -> bool {
    if !rec
        .aux_iter()
        .flatten()
        .any(|(tag, _)| OUR_TAGS.iter().any(|t| t.as_slice() == tag))
    {
        return false;
    }
    for tag in OUR_TAGS {
        while rec.remove_aux(tag).is_ok() {}
    }
    true
}

fn tag_record(
    rec: &mut bam::Record,
    plan: &Plan,
    ann: Option<&Annotation>,
    m: u32,
    is_rep: bool,
    polya_on: bool,
) -> Result<()> {
    strip_stale_tags(rec);
    let mol = &plan.molecules[m as usize];
    set_tag(rec, b"MI", Aux::String(&plan.mi[m as usize].to_string()))?;
    set_tag(rec, b"uN", Aux::String(&umi::to_string(mol.umi)))?;
    set_tag(rec, b"uS", Aux::I32(mol.n_reads as i32))?;
    if !is_rep {
        return Ok(());
    }
    set_tag(rec, b"uJ", Aux::I32(mol.same_chain as i32))?;
    set_tag(rec, b"uD", Aux::I32(mol.disagreements as i32))?;
    if polya_on {
        set_tag(rec, b"uA", Aux::Char(mol.polya))?;
        if mol.polya_len > 0 {
            set_tag(rec, b"uP", Aux::I32(mol.polya_len as i32))?;
        }
    }
    if mol.fl != 0 {
        set_tag(rec, b"uF", Aux::Char(mol.fl))?;
        if let (Some(a), true) = (ann, mol.fl_tx != NONE) {
            set_tag(rec, b"uT", Aux::String(&a.transcript(mol.fl_tx).id))?;
        }
    }
    if mol.genes != NONE {
        set_tag(rec, b"uG", Aux::String(&plan.gene_lists[mol.genes as usize]))?;
    }
    Ok(())
}

fn write_molecule_row(
    w: &mut impl Write,
    rec: &bam::Record,
    header: &bam::HeaderView,
    plan: &Plan,
    ann: Option<&Annotation>,
    m: u32,
) -> Result<()> {
    let mol = &plan.molecules[m as usize];
    let fl = if mol.fl == 0 { '.' } else { mol.fl as char };
    let tx = match (ann, mol.fl_tx) {
        (Some(a), t) if t != NONE => a.transcript(t).id.as_str(),
        _ => ".",
    };
    let genes = if mol.genes == NONE {
        "."
    } else {
        plan.gene_lists[mol.genes as usize].as_str()
    };
    writeln!(
        w,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        plan.mi[m as usize],
        String::from_utf8_lossy(header.tid2name(mol.tid)),
        mol.strand.as_char(),
        u64::from(mol.anchor) + 1,
        umi::to_string(mol.umi),
        if umi::is_damaged(mol.umi) { 'Y' } else { 'N' },
        mol.n_reads,
        mol.same_chain,
        mol.disagreements,
        String::from_utf8_lossy(rec.qname()),
        rec.pos() + 1,
        rec.cigar().end_pos(),
        mol.polya as char,
        mol.polya_len,
        fl,
        tx,
        genes,
    )?;
    Ok(())
}

pub fn run(cfg: &Config) -> Result<Stats> {
    let t0 = Instant::now();
    let threads = cfg.threads.max(1);
    let pool = ThreadPool::new(threads as u32).context("creating thread pool")?;

    let mut reader = bam::Reader::from_path(&cfg.input).with_context(|| format!("opening {}", cfg.input.display()))?;
    reader.set_thread_pool(&pool)?;
    check_sorted(reader.header())?;

    let annotation = match &cfg.gtf {
        Some(path) => {
            let header = reader.header();
            let names: FxHashMap<Vec<u8>, u32> = (0..header.target_count())
                .map(|t| (header.tid2name(t).to_vec(), t))
                .collect();
            let a = Annotation::load(path, &names, header.target_count() as usize)?;
            info!(
                "annotation: {} genes, {} transcripts on the BAM's references",
                a.n_genes(),
                a.n_transcripts()
            );
            if cfg.full_length && a.n_transcripts() == 0 {
                warn!("--full-length: the GTF has no exon records, so uF/uT are not assessed");
            }
            Some(a)
        }
        None => None,
    };

    let mut pass1 = Pass1::new(cfg, annotation.as_ref());
    let mut rec = bam::Record::new();
    let mut ordinal = 0u64;
    while let Some(r) = reader.read(&mut rec) {
        r.with_context(|| format!("reading record {ordinal}"))?;
        pass1.add(&rec, ordinal)?;
        ordinal += 1;
    }
    drop(reader);
    let (plan, mut stats) = pass1.finish();
    if stats.reads_used > 0 && stats.strand_unknown * 2 > stats.reads_used && cfg.use_strand {
        warn!(
            "{} of {} reads have no strand (no TS:A from dorado, no ts:A from minimap2); \
             they are grouped by locus and UMI only",
            stats.strand_unknown, stats.reads_used
        );
    }
    info!(
        "pass 1: {} records, {} reads with a UMI, {} molecules in {:.1}s",
        stats.total_records,
        stats.reads_used,
        stats.molecules,
        t0.elapsed().as_secs_f64()
    );

    // Pass 2.
    let t1 = Instant::now();
    let mut reader = bam::Reader::from_path(&cfg.input)?;
    reader.set_thread_pool(&pool)?;
    let mut header = bam::Header::from_template(reader.header());
    let (pg_id, prev) = pg_record(reader.header());
    let cl = cfg.command_line.replace(['\t', '\n'], " ");
    let mut pg = bam::header::HeaderRecord::new(b"PG");
    pg.push_tag(b"ID", &pg_id).push_tag(b"PN", "flumi");
    if let Some(p) = &prev {
        pg.push_tag(b"PP", p);
    }
    pg.push_tag(b"VN", env!("CARGO_PKG_VERSION")).push_tag(b"CL", &cl);
    header.push_record(&pg);
    let mut writer = bam::Writer::from_path(&cfg.output, &header, bam::Format::Bam)
        .with_context(|| format!("creating {}", cfg.output.display()))?;
    writer.set_thread_pool(&pool)?;
    let header_view = reader.header().clone();
    let mut tsv = match &cfg.molecules {
        Some(p) => {
            let mut w = BufWriter::new(File::create(p).with_context(|| format!("creating {}", p.display()))?);
            writeln!(
                w,
                "MI\tchrom\tstrand\tanchor\tumi\tumi_damaged\treads\tsame_chain\tdisagreements\t\
                 rep_read\trep_start\trep_end\tpolyA\tpolyA_len\tfull_length\ttranscript\tgenes"
            )?;
            Some(w)
        }
        None => None,
    };

    let ann = annotation.as_ref();
    let polya_on = cfg.full_length;
    let mut next_rep = 0usize;
    let mut ordinal = 0u64;
    while let Some(r) = reader.read(&mut rec) {
        r.with_context(|| format!("reading record {ordinal} (pass 2)"))?;
        match cfg.mode {
            Mode::Filter => {
                if plan.reps.get(next_rep).is_some_and(|&(o, _)| o == ordinal) {
                    let m = plan.reps[next_rep].1;
                    next_rep += 1;
                    tag_record(&mut rec, &plan, ann, m, true, polya_on)?;
                    if let Some(w) = tsv.as_mut() {
                        write_molecule_row(w, &rec, &header_view, &plan, ann, m)?;
                    }
                    writer.write(&rec)?;
                    stats.records_written += 1;
                }
            }
            Mode::Mark => {
                let m = *plan
                    .assign
                    .get(ordinal as usize)
                    .context("input changed between passes")?;
                if m != NONE {
                    let is_rep = plan.molecules[m as usize].rep_ordinal == ordinal;
                    tag_record(&mut rec, &plan, ann, m, is_rep, polya_on)?;
                    let flags = rec.flags();
                    rec.set_flags(if is_rep {
                        flags & !FLAG_DUPLICATE
                    } else {
                        flags | FLAG_DUPLICATE
                    });
                    if is_rep {
                        if let Some(w) = tsv.as_mut() {
                            write_molecule_row(w, &rec, &header_view, &plan, ann, m)?;
                        }
                    }
                } else if strip_stale_tags(&mut rec) {
                    // Assigned by an earlier run, not by this one: its
                    // duplicate flag was ours too.
                    rec.set_flags(rec.flags() & !FLAG_DUPLICATE);
                }
                writer.write(&rec)?;
                stats.records_written += 1;
            }
        }
        ordinal += 1;
    }
    if ordinal != plan.total_records || next_rep != plan.reps.len() {
        bail!(
            "input changed between passes ({} records, expected {})",
            ordinal,
            plan.total_records
        );
    }
    if let Some(mut w) = tsv {
        w.flush()?;
    }
    // rust-htslib's Drop ignores hts_close errors, so check that the final
    // flush reached the file: a complete BGZF stream ends in the EOF block.
    drop(writer);
    check_bgzf_eof(&cfg.output)?;
    info!(
        "pass 2: wrote {} records in {:.1}s",
        stats.records_written,
        t1.elapsed().as_secs_f64()
    );
    if let Some(path) = &cfg.stats {
        stats
            .write_tsv(path)
            .with_context(|| format!("writing {}", path.display()))?;
    }
    Ok(stats)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn truncated_output_is_detected() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("x.bam");
        let mut header = bam::Header::new();
        let mut sq = bam::header::HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", "chr1").push_tag(b"LN", 1000);
        header.push_record(&sq);
        drop(bam::Writer::from_path(&path, &header, bam::Format::Bam).unwrap());
        check_bgzf_eof(&path).unwrap();
        let bytes = std::fs::read(&path).unwrap();
        std::fs::write(&path, &bytes[..bytes.len() - 28]).unwrap();
        assert!(check_bgzf_eof(&path).is_err());
        std::fs::write(&path, b"").unwrap();
        assert!(check_bgzf_eof(&path).is_err());
    }
}
