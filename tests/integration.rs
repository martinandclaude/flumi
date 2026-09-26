//! End-to-end runs over small BAMs written with rust-htslib.

use std::path::{Path, PathBuf};

use rust_htslib::bam::{self, record::Aux, record::CigarString, Read};

use flumi::consensus::{RepMode, RepParams};
use flumi::dedup::{self, Config, Mode};

const B: [&str; 4] = ["ACGA", "CGAC", "GACG", "ACGC"];

fn rx(anchors: [&str; 5], blocks: [&str; 4]) -> String {
    let mut s = String::new();
    for k in 0..4 {
        s.push_str(anchors[k]);
        s.push_str(blocks[k]);
    }
    s + anchors[4]
}

fn clean(blocks: [&str; 4]) -> String {
    rx(["TTT", "TT", "TT", "TT", "TTT"], blocks)
}

#[derive(Clone)]
struct R {
    name: String,
    tid: i32,
    pos: i64,
    cigar: String,
    flags: u16,
    rx: Option<String>,
    ts: Option<u8>,
    minimap_ts: Option<u8>,
    pt: Option<i32>,
    sa: bool,
    mapq: u8,
}

fn r(name: &str, pos: i64, cigar: &str, rx: Option<String>) -> R {
    R {
        name: name.into(),
        tid: 0,
        pos,
        cigar: cigar.into(),
        flags: 0,
        rx,
        ts: Some(b'+'),
        minimap_ts: None,
        pt: None,
        sa: false,
        mapq: 60,
    }
}

fn write_bam(dir: &Path, reads: &[R], sorted: bool) -> PathBuf {
    let mut header = bam::Header::new();
    let mut hd = bam::header::HeaderRecord::new(b"HD");
    hd.push_tag(b"VN", "1.6")
        .push_tag(b"SO", if sorted { "coordinate" } else { "unsorted" });
    header.push_record(&hd);
    for name in ["chr1", "chr2"] {
        let mut sq = bam::header::HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", name).push_tag(b"LN", 1_000_000);
        header.push_record(&sq);
    }
    let path = dir.join("in.bam");
    let mut w = bam::Writer::from_path(&path, &header, bam::Format::Bam).unwrap();
    let mut reads = reads.to_vec();
    if sorted {
        reads.sort_by_key(|x| (x.tid, x.pos));
    }
    for x in &reads {
        let cigar = CigarString::try_from(x.cigar.as_str()).unwrap();
        let qlen = cigar
            .iter()
            .map(|c| match c {
                bam::record::Cigar::Match(n) | bam::record::Cigar::Ins(n) | bam::record::Cigar::SoftClip(n) => *n,
                _ => 0,
            })
            .sum::<u32>() as usize;
        let mut rec = bam::Record::new();
        rec.set(x.name.as_bytes(), Some(&cigar), &vec![b'A'; qlen], &vec![30u8; qlen]);
        rec.set_tid(x.tid);
        rec.set_pos(x.pos);
        rec.set_mapq(x.mapq);
        rec.set_flags(x.flags);
        rec.set_mtid(-1);
        rec.set_mpos(-1);
        if let Some(s) = &x.rx {
            rec.push_aux(b"RX", Aux::String(s)).unwrap();
        }
        if let Some(c) = x.ts {
            rec.push_aux(b"TS", Aux::Char(c)).unwrap();
        }
        if let Some(c) = x.minimap_ts {
            rec.push_aux(b"ts", Aux::Char(c)).unwrap();
        }
        rec.push_aux(b"qs", Aux::Float(20.0)).unwrap();
        if let Some(pt) = x.pt {
            rec.push_aux(b"pt", Aux::I32(pt)).unwrap();
        }
        if x.sa {
            rec.push_aux(b"SA", Aux::String("chr2,5000,+,300M,60,0;")).unwrap();
        }
        w.write(&rec).unwrap();
    }
    path
}

fn config(dir: &Path, input: PathBuf, mode: Mode) -> Config {
    Config {
        input,
        output: dir.join("out.bam"),
        gtf: None,
        stats: Some(dir.join("stats.tsv")),
        molecules: Some(dir.join("molecules.tsv")),
        mode,
        umi_dist: 2,
        damaged_umi_dist: 3,
        valid_budget: 0.01,
        damaged_budget: 0.5,
        max_umi_edits: 5,
        pos_tol: 20,
        use_position: true,
        use_strand: true,
        bridge_exon_max: 50,
        chimera_clip: 60,
        min_mapq: 1,
        min_qs: 0.0,
        rep: RepParams {
            mode: RepMode::Consensus,
            strict: true,
        },
        full_length: true,
        fl_cov: 0.8,
        fl_terminal: 25,
        threads: 1,
        command_line: "flumi test".into(),
    }
}

struct Out {
    name: String,
    pos: i64,
    dup: bool,
    tags: Vec<(String, String)>,
}

impl Out {
    fn tag(&self, t: &str) -> Option<&str> {
        self.tags.iter().find(|(k, _)| k == t).map(|(_, v)| v.as_str())
    }
}

fn read_out(path: &Path) -> Vec<Out> {
    let mut reader = bam::Reader::from_path(path).unwrap();
    reader
        .records()
        .map(|rec| {
            let rec = rec.unwrap();
            let tags = rec
                .aux_iter()
                .map(|t| {
                    let (k, v) = t.unwrap();
                    let v = match v {
                        Aux::String(s) => s.to_string(),
                        Aux::Char(c) => (c as char).to_string(),
                        Aux::I8(n) => n.to_string(),
                        Aux::U8(n) => n.to_string(),
                        Aux::I16(n) => n.to_string(),
                        Aux::U16(n) => n.to_string(),
                        Aux::I32(n) => n.to_string(),
                        Aux::U32(n) => n.to_string(),
                        Aux::Float(f) => f.to_string(),
                        _ => String::new(),
                    };
                    (String::from_utf8_lossy(k).into_owned(), v)
                })
                .collect();
            Out {
                name: String::from_utf8_lossy(rec.qname()).into_owned(),
                pos: rec.pos(),
                dup: rec.is_duplicate(),
                tags,
            }
        })
        .collect()
}

fn stat(dir: &Path, key: &str) -> String {
    let text = std::fs::read_to_string(dir.join("stats.tsv")).unwrap();
    text.lines()
        .find_map(|l| l.strip_prefix(&format!("{key}\t")))
        .unwrap_or_else(|| panic!("{key}"))
        .to_string()
}

fn pcr_family() -> Vec<R> {
    vec![
        r("short", 100, "3S200M", Some(clean(B))),
        r("long", 101, "3S400M", Some(clean(B))),
        // Anchor indel: RX is 27 nt but the blocks are intact.
        r("drift", 102, "3S300M", Some(rx(["TT", "TT", "TT", "TT", "TTT"], B))),
        // Block deletion: damaged UMI one edit from the family's.
        r("damaged", 100, "3S250M", Some(clean(["ACG", "CGAC", "GACG", "ACGC"]))),
        r("no_rx", 100, "3S300M", None),
    ]
}

#[test]
fn filter_mode_writes_one_tagged_representative() {
    let dir = tempfile::tempdir().unwrap();
    let input = write_bam(dir.path(), &pcr_family(), true);
    dedup::run(&config(dir.path(), input, Mode::Filter)).unwrap();
    let out = read_out(&dir.path().join("out.bam"));
    assert_eq!(out.len(), 1);
    let rep = &out[0];
    assert_eq!(rep.name, "long");
    assert_eq!(rep.tag("uS"), Some("4"));
    assert_eq!(rep.tag("uN"), Some("ACGACGACGACGACGC"));
    assert_eq!(rep.tag("MI"), Some("1"));
    assert_eq!(rep.tag("uJ"), Some("4"));
    assert_eq!(stat(dir.path(), "molecules_out"), "1");
    assert_eq!(stat(dir.path(), "dropped_no_RX"), "1");
    assert_eq!(stat(dir.path(), "umi_block_damaged"), "1");
    assert_eq!(stat(dir.path(), "damaged_to_valid"), "1");
    let tsv = std::fs::read_to_string(dir.path().join("molecules.tsv")).unwrap();
    let row: Vec<&str> = tsv.lines().nth(1).unwrap().split('\t').collect();
    assert_eq!(&row[..7], &["1", "chr1", "+", "101", "ACGACGACGACGACGC", "N", "4"]);
    assert_eq!(row[9], "long");
}

#[test]
fn mark_mode_keeps_every_record_and_flags_duplicates() {
    let dir = tempfile::tempdir().unwrap();
    let input = write_bam(dir.path(), &pcr_family(), true);
    dedup::run(&config(dir.path(), input, Mode::Mark)).unwrap();
    let out = read_out(&dir.path().join("out.bam"));
    assert_eq!(out.len(), 5);
    let positions: Vec<i64> = out.iter().map(|o| o.pos).collect();
    assert!(positions.windows(2).all(|w| w[0] <= w[1]));
    for o in &out {
        match o.name.as_str() {
            "no_rx" => assert!(o.tag("MI").is_none() && !o.dup),
            "long" => assert!(!o.dup && o.tag("uJ").is_some()),
            _ => assert!(o.dup && o.tag("MI") == Some("1") && o.tag("uJ").is_none()),
        }
    }
}

#[test]
fn same_umi_at_distant_five_prime_ends_is_two_molecules() {
    let dir = tempfile::tempdir().unwrap();
    let reads = vec![
        r("a1", 100, "3S300M", Some(clean(B))),
        r("a2", 104, "3S300M", Some(clean(B))),
        r("b1", 5000, "3S300M", Some(clean(B))),
    ];
    let input = write_bam(dir.path(), &reads, true);
    dedup::run(&config(dir.path(), input, Mode::Filter)).unwrap();
    assert_eq!(read_out(&dir.path().join("out.bam")).len(), 2);
    assert_eq!(stat(dir.path(), "umi_position_splits"), "0");

    // Overlapping alignments, 5' ends 300 bp apart: one bundle, split by anchor.
    let reads = vec![
        r("a1", 100, "3S600M", Some(clean(B))),
        r("b1", 400, "3S600M", Some(clean(B))),
    ];
    let input = write_bam(dir.path(), &reads, true);
    dedup::run(&config(dir.path(), input, Mode::Filter)).unwrap();
    assert_eq!(read_out(&dir.path().join("out.bam")).len(), 2);
    assert_eq!(stat(dir.path(), "umi_position_splits"), "1");
}

#[test]
fn minus_strand_molecules_are_anchored_at_the_alignment_end() {
    let dir = tempfile::tempdir().unwrap();
    // Forward reads (TS:+) aligned in reverse: transcript on the minus strand,
    // 5' end at the right. Copies truncated at the 3' end start at different
    // positions but end together.
    let mut reads = vec![
        r("c1", 1000, "400M3S", Some(clean(B))),
        r("c2", 1200, "200M3S", Some(clean(B))),
        r("c3", 1300, "100M3S", Some(clean(B))),
        // Same UMI, different 5' end on the minus strand.
        r("other", 1250, "100M3S", Some(clean(B))),
    ];
    for x in &mut reads {
        x.flags = 0x10;
    }
    let input = write_bam(dir.path(), &reads, true);
    dedup::run(&config(dir.path(), input, Mode::Filter)).unwrap();
    let out = read_out(&dir.path().join("out.bam"));
    let names: Vec<&str> = out.iter().map(|o| o.name.as_str()).collect();
    assert_eq!(names, vec!["c1", "other"]);
    assert_eq!(out[0].tag("uS"), Some("3"));
}

#[test]
fn antisense_molecules_with_the_same_umi_stay_apart() {
    let dir = tempfile::tempdir().unwrap();
    let mut minus = r("minus", 100, "300M3S", Some(clean(B)));
    minus.ts = Some(b'-');
    minus.flags = 0;
    let mut minus_rev = r("plus_rev", 100, "3S300M", Some(clean(B)));
    minus_rev.ts = Some(b'-');
    minus_rev.flags = 0x10;
    let reads = vec![r("plus", 100, "3S300M", Some(clean(B))), minus, minus_rev];
    let input = write_bam(dir.path(), &reads, true);
    dedup::run(&config(dir.path(), input, Mode::Filter)).unwrap();
    // TS:- unreversed = minus strand; TS:- reversed = plus strand.
    let out = read_out(&dir.path().join("out.bam"));
    assert_eq!(out.len(), 2);
    assert_eq!(stat(dir.path(), "strand_from_TS"), "3");
}

#[test]
fn chimeric_umi_end_is_dropped() {
    let dir = tempfile::tempdir().unwrap();
    let mut chim = r(
        "chimera",
        100,
        "250S300M",
        Some(clean(["CCGA", "GGAC", "GACG", "ACGC"])),
    );
    chim.sa = true;
    let mut clipped_no_sa = r(
        "clipped",
        100,
        "250S300M",
        Some(clean(["GAGA", "GGAC", "GACG", "ACGC"])),
    );
    clipped_no_sa.sa = false;
    let reads = vec![r("real", 100, "3S300M", Some(clean(B))), chim, clipped_no_sa];
    let input = write_bam(dir.path(), &reads, true);
    dedup::run(&config(dir.path(), input, Mode::Filter)).unwrap();
    let names: Vec<String> = read_out(&dir.path().join("out.bam"))
        .into_iter()
        .map(|o| o.name)
        .collect();
    assert_eq!(names, vec!["real", "clipped"]);
    assert_eq!(stat(dir.path(), "dropped_chimeric_umi_end"), "1");
}

#[test]
fn unsorted_input_is_rejected() {
    let dir = tempfile::tempdir().unwrap();
    let input = write_bam(dir.path(), &pcr_family(), false);
    let err = dedup::run(&config(dir.path(), input, Mode::Filter)).unwrap_err();
    assert!(format!("{err:#}").contains("coordinate-sorted"));

    // Header claims sorted, records are not.
    let reads = vec![
        r("b", 500, "3S100M", Some(clean(B))),
        r("a", 100, "3S100M", Some(clean(B))),
    ];
    let dir2 = tempfile::tempdir().unwrap();
    let path = write_bam(dir2.path(), &reads, true);
    let mut reader = bam::Reader::from_path(&path).unwrap();
    let header = bam::Header::from_template(reader.header());
    let unsorted = dir2.path().join("lying.bam");
    let mut w = bam::Writer::from_path(&unsorted, &header, bam::Format::Bam).unwrap();
    let mut recs: Vec<bam::Record> = reader.records().map(|x| x.unwrap()).collect();
    recs.reverse();
    for rec in &recs {
        w.write(rec).unwrap();
    }
    drop(w);
    let err = dedup::run(&config(dir2.path(), unsorted, Mode::Filter)).unwrap_err();
    assert!(format!("{err:#}").contains("not coordinate-sorted"));
}

#[test]
fn reads_without_strand_are_grouped_by_umi_within_locus() {
    let dir = tempfile::tempdir().unwrap();
    let mut reads = pcr_family();
    for x in &mut reads {
        x.ts = None;
    }
    reads.push(r("far_start", 150, "3S300M", Some(clean(B))));
    reads.last_mut().unwrap().ts = None;
    let input = write_bam(dir.path(), &reads, true);
    dedup::run(&config(dir.path(), input, Mode::Filter)).unwrap();
    let out = read_out(&dir.path().join("out.bam"));
    // Without a strand the 5' end is unknown, so position cannot split them.
    assert_eq!(out.len(), 1);
    assert_eq!(out[0].tag("uS"), Some("5"));
    assert_eq!(stat(dir.path(), "strand_unknown"), "5");
}

#[test]
fn annotation_and_polya_tags() {
    let dir = tempfile::tempdir().unwrap();
    let gtf = dir.path().join("genes.gtf");
    std::fs::write(
        &gtf,
        "chr1\tt\texon\t101\t300\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";\n\
         chr1\tt\texon\t401\t600\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";\n\
         chr1\tt\tgene\t50\t700\t.\t-\t.\tgene_id \"ANTI\";\n",
    )
    .unwrap();
    let mut a = r("full", 100, "3S200M100N200M40S", Some(clean(B)));
    a.pt = Some(62);
    let mut b = r("copy", 100, "3S200M100N150M", Some(clean(B)));
    b.pt = Some(-1);
    let mut c = r("copy2", 100, "3S200M100N150M", Some(clean(B)));
    c.pt = Some(58);
    let input = write_bam(dir.path(), &[a, b, c], true);
    let mut cfg = config(dir.path(), input, Mode::Filter);
    cfg.gtf = Some(gtf);
    dedup::run(&cfg).unwrap();
    let out = read_out(&dir.path().join("out.bam"));
    assert_eq!(out.len(), 1);
    let rep = &out[0];
    assert_eq!(rep.name, "full");
    assert_eq!(rep.tag("uG"), Some("G1"));
    assert_eq!(rep.tag("uF"), Some("Y"));
    assert_eq!(rep.tag("uT"), Some("T1"));
    // Poly(A) evidence pooled over copies; length is their median.
    assert_eq!(rep.tag("uA"), Some("Y"));
    assert_eq!(rep.tag("uP"), Some("60"));
    assert_eq!(stat(dir.path(), "full_length_both_signals"), "1");
}

#[test]
fn polya_anchor_outranks_a_longer_read_without_it() {
    let dir = tempfile::tempdir().unwrap();
    let mut long = r("long_no_anchor", 100, "3S450M", Some(clean(B)));
    long.pt = Some(-1);
    let mut fl = r("reached_polya", 100, "3S400M60S", Some(clean(B)));
    fl.pt = Some(55);
    let input = write_bam(dir.path(), &[long, fl], true);
    dedup::run(&config(dir.path(), input, Mode::Filter)).unwrap();
    let out = read_out(&dir.path().join("out.bam"));
    assert_eq!(out[0].name, "reached_polya");
    assert_eq!(stat(dir.path(), "rep_switched_by_consensus"), "1");
}

#[test]
fn rerunning_on_own_output_is_stable() {
    let dir = tempfile::tempdir().unwrap();
    let input = write_bam(dir.path(), &pcr_family(), true);
    dedup::run(&config(dir.path(), input, Mode::Mark)).unwrap();
    let first = dir.path().join("first.bam");
    std::fs::rename(dir.path().join("out.bam"), &first).unwrap();
    dedup::run(&config(dir.path(), first.clone(), Mode::Mark)).unwrap();
    let a = read_out(&first);
    let b = read_out(&dir.path().join("out.bam"));
    assert_eq!(a.len(), b.len());
    for (x, y) in a.iter().zip(&b) {
        assert_eq!(
            (&x.name, x.dup, x.tag("MI"), x.tag("uS")),
            (&y.name, y.dup, y.tag("MI"), y.tag("uS"))
        );
        // No duplicated tags after the second pass.
        let mut keys: Vec<&String> = y.tags.iter().map(|(k, _)| k).collect();
        keys.sort();
        keys.dedup();
        assert_eq!(keys.len(), y.tags.len());
    }
    let header = bam::Reader::from_path(dir.path().join("out.bam"))
        .unwrap()
        .header()
        .clone();
    let text = String::from_utf8_lossy(header.as_bytes()).into_owned();
    assert!(text.contains("ID:flumi\t") && text.contains("ID:flumi.1\tPN:flumi\tPP:flumi"));
}

#[test]
fn unresolved_dorado_orientation_falls_back_to_minimap2() {
    let dir = tempfile::tempdir().unwrap();
    let mut a = r("a", 100, "3S200M100N200M", Some(clean(B)));
    a.ts = Some(b'?');
    a.minimap_ts = Some(b'+');
    let mut b = r("b", 101, "3S200M100N100M", Some(clean(B)));
    b.minimap_ts = Some(b'+');
    let input = write_bam(dir.path(), &[a, b], true);
    dedup::run(&config(dir.path(), input, Mode::Filter)).unwrap();
    assert_eq!(stat(dir.path(), "strand_from_ts"), "1");
    assert_eq!(stat(dir.path(), "strand_from_TS"), "1");
    assert_eq!(read_out(&dir.path().join("out.bam")).len(), 1);
}

#[test]
fn remarking_with_stricter_filters_clears_stale_marks() {
    let dir = tempfile::tempdir().unwrap();
    let input = write_bam(dir.path(), &pcr_family(), true);
    dedup::run(&config(dir.path(), input, Mode::Mark)).unwrap();
    let first = dir.path().join("first.bam");
    std::fs::rename(dir.path().join("out.bam"), &first).unwrap();
    assert!(read_out(&first).iter().any(|o| o.dup));
    let mut cfg = config(dir.path(), first, Mode::Mark);
    cfg.min_mapq = 61; // nothing passes
    dedup::run(&cfg).unwrap();
    for o in read_out(&dir.path().join("out.bam")) {
        assert!(!o.dup, "{} keeps a stale duplicate flag", o.name);
        assert!(
            o.tag("MI").is_none() && o.tag("uS").is_none(),
            "{} keeps stale tags",
            o.name
        );
    }
}
