use std::path::PathBuf;
use std::process::ExitCode;

use clap::{Parser, ValueEnum};
use log::error;

use flumi::consensus::{RepMode, RepParams};
use flumi::dedup::{self, Config, Mode};

#[derive(Clone, Copy, Debug, ValueEnum)]
enum ModeArg {
    /// One representative record per molecule.
    Filter,
    /// Every record; molecule members tagged, non-representatives flagged 0x400.
    Mark,
}

#[derive(Clone, Copy, Debug, ValueEnum)]
enum RepArg {
    /// Majority-vote junction consensus, then poly(A) completeness, then length.
    JunctionConsensus,
    /// Largest exonic footprint.
    Longest,
}

/// Collapse ONT PCR-cDNA (SQK-PCB114.24) reads to one read per original RNA
/// molecule using the pre-PCR UMI dorado reports in RX:Z.
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Cli {
    /// Coordinate-sorted BAM, spliced alignments (minimap2 -ax splice), with dorado's RX:Z and TS:A.
    #[arg(long, short = 'i')]
    bam: PathBuf,
    /// Output BAM (coordinate-sorted, same order as the input).
    #[arg(long, short = 'o')]
    out: PathBuf,
    /// Optional GTF (.gz OK): gene labels (uG), full-length calls (uF/uT), annotated introns.
    #[arg(long)]
    gtf: Option<PathBuf>,
    /// Run statistics TSV.
    #[arg(long)]
    stats: Option<PathBuf>,
    /// Per-molecule TSV (one row per representative).
    #[arg(long)]
    molecules: Option<PathBuf>,
    /// Output mode.
    #[arg(long, value_enum, default_value_t = ModeArg::Filter)]
    mode: ModeArg,
    /// Largest blockwise edit distance between two undamaged UMIs of one molecule
    /// (directional rule; beyond 1 edit only where the density gate allows).
    #[arg(long, default_value_t = 2)]
    umi_dist: u32,
    /// Largest blockwise edit distance between a damaged UMI and its molecule
    /// (beyond 1 edit only where the density gate allows).
    #[arg(long, default_value_t = 3)]
    damaged_umi_dist: u32,
    /// Density gate for undamaged merges beyond one edit: the most unrelated
    /// UMIs expected within that distance by chance among the molecules sharing
    /// the 5' neighbourhood.
    #[arg(long, default_value_t = 0.01)]
    collision_budget: f64,
    /// Density gate for damaged UMIs beyond one edit (looser: a damaged UMI is
    /// certainly an error of some molecule).
    #[arg(long, default_value_t = 0.5)]
    damaged_collision_budget: f64,
    /// Maximum edits between RX and the UMI pattern (dorado accepts at most 5).
    #[arg(long, default_value_t = 5)]
    max_umi_edits: u8,
    /// Largest 5'-end difference (bp) between PCR copies of one molecule.
    #[arg(long, default_value_t = 20)]
    pos_tol: u32,
    /// Ignore 5'-end positions: merge by UMI within a strand-specific locus only.
    #[arg(long)]
    no_position: bool,
    /// Ignore read strand (TS:A / ts:A); implies --no-position.
    #[arg(long)]
    ignore_strand: bool,
    /// Longest first exon for which a copy that soft-clipped it is bridged to one that aligned it.
    #[arg(long, default_value_t = 50)]
    bridge_exon_max: u32,
    /// Drop reads with a supplementary alignment whose UMI-side soft clip exceeds this (PCR chimeras).
    #[arg(long, default_value_t = 60)]
    chimera_clip: u32,
    /// Minimum mapping quality.
    #[arg(long, default_value_t = 1)]
    min_mapq: u8,
    /// Minimum mean basecall quality (qs tag).
    #[arg(long, default_value_t = 0.0)]
    min_qs: f32,
    /// How to choose each molecule's representative read.
    #[arg(long, value_enum, default_value_t = RepArg::JunctionConsensus)]
    rep: RepArg,
    /// Do not penalise junctions only one copy shows (default: penalise unless annotated).
    #[arg(long)]
    lenient_junctions: bool,
    /// Skip full-length calling (uF/uT from the GTF, uA/uP from dorado's pt:i).
    #[arg(long)]
    no_full_length: bool,
    /// Fraction of a transcript's exonic length a full-length read must exceed.
    #[arg(long, default_value_t = 0.8)]
    fl_cov: f64,
    /// Bases a full-length read must reach into the transcript's first and last exon.
    #[arg(long, default_value_t = 25)]
    fl_terminal: u32,
    /// Threads for BAM compression and decompression.
    #[arg(long, short = 't', default_value_t = 4)]
    threads: usize,
    /// Debug logging.
    #[arg(long, short = 'v')]
    verbose: bool,
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    env_logger::Builder::new()
        .filter_level(if cli.verbose {
            log::LevelFilter::Debug
        } else {
            log::LevelFilter::Info
        })
        .format_timestamp_secs()
        .init();
    let cfg = Config {
        input: cli.bam,
        output: cli.out,
        gtf: cli.gtf,
        stats: cli.stats,
        molecules: cli.molecules,
        mode: match cli.mode {
            ModeArg::Filter => Mode::Filter,
            ModeArg::Mark => Mode::Mark,
        },
        umi_dist: cli.umi_dist,
        damaged_umi_dist: cli.damaged_umi_dist,
        valid_budget: cli.collision_budget,
        damaged_budget: cli.damaged_collision_budget,
        max_umi_edits: cli.max_umi_edits,
        pos_tol: cli.pos_tol,
        use_position: !cli.no_position && !cli.ignore_strand,
        use_strand: !cli.ignore_strand,
        bridge_exon_max: cli.bridge_exon_max,
        chimera_clip: cli.chimera_clip,
        min_mapq: cli.min_mapq,
        min_qs: cli.min_qs,
        rep: RepParams {
            mode: match cli.rep {
                RepArg::JunctionConsensus => RepMode::Consensus,
                RepArg::Longest => RepMode::Longest,
            },
            strict: !cli.lenient_junctions,
        },
        full_length: !cli.no_full_length,
        fl_cov: cli.fl_cov,
        fl_terminal: cli.fl_terminal,
        threads: cli.threads,
        command_line: std::env::args().collect::<Vec<_>>().join(" "),
    };
    match dedup::run(&cfg) {
        Ok(_) => ExitCode::SUCCESS,
        Err(e) => {
            error!("{e:#}");
            ExitCode::FAILURE
        }
    }
}
