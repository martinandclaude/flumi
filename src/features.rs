//! Per-read features: the tags and alignment facts the deduplication needs.

use rust_htslib::bam::record::Aux;
use rust_htslib::bam::Record;

use crate::umi::{self, DecodeError, Umi};

/// `pt:i` was not written (poly(A) estimation off).
pub const PT_ABSENT: i32 = -2;
/// `pt:i == -1`: dorado did not find the poly(A) primer anchor.
pub const PT_NO_ANCHOR: i32 = -1;

pub const NO_POS: u32 = u32::MAX;

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Strand {
    Plus = 0,
    Minus = 1,
    Unknown = 2,
}

impl Strand {
    pub fn as_char(self) -> char {
        match self {
            Strand::Plus => '+',
            Strand::Minus => '-',
            Strand::Unknown => '.',
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum StrandSource {
    /// dorado's `TS:A`: read orientation from primer classification.
    Dorado,
    /// minimap2's `ts:A`: transcript strand from splice motifs.
    Minimap2,
    None,
}

/// Transcript strand on the genome from a read-relative orientation.
///
/// Both dorado's `TS:A` and minimap2's `ts:A` describe the read as
/// basecalled; the alignment may have reverse-complemented it.
pub fn genomic_strand(read_orientation: u8, is_reverse: bool) -> Strand {
    match (read_orientation, is_reverse) {
        (b'+', false) | (b'-', true) => Strand::Plus,
        (b'+', true) | (b'-', false) => Strand::Minus,
        _ => Strand::Unknown,
    }
}

#[derive(Debug, Default)]
pub struct Tags {
    pub umi: Option<Result<Umi, DecodeError>>,
    pub dorado_ts: Option<u8>,
    pub minimap_ts: Option<u8>,
    pub qs: Option<f32>,
    pub pt: i32,
    pub has_sa: bool,
}

fn aux_int(a: &Aux<'_>) -> Option<i64> {
    Some(match *a {
        Aux::I8(v) => v.into(),
        Aux::U8(v) => v.into(),
        Aux::I16(v) => v.into(),
        Aux::U16(v) => v.into(),
        Aux::I32(v) => v.into(),
        Aux::U32(v) => v.into(),
        _ => return None,
    })
}

fn aux_float(a: &Aux<'_>) -> Option<f32> {
    match *a {
        Aux::Float(v) => Some(v),
        Aux::Double(v) => Some(v as f32),
        _ => aux_int(a).map(|v| v as f32),
    }
}

/// One pass over the aux fields. Iterating (rather than looking tags up by
/// name) matters for `ts`: dorado's `ts:i` (signal trim) and minimap2's
/// `ts:A` (transcript strand) can both be present, and a lookup returns
/// whichever comes first.
pub fn scan_tags(rec: &Record, max_umi_edits: u8) -> Tags {
    let mut t = Tags {
        pt: PT_ABSENT,
        ..Tags::default()
    };
    for item in rec.aux_iter() {
        let Ok((tag, aux)) = item else { break };
        match tag {
            b"RX" => {
                if let Aux::String(s) = aux {
                    t.umi = Some(umi::decode(s.as_bytes(), max_umi_edits));
                }
            }
            b"TS" => {
                if let Aux::Char(c) = aux {
                    t.dorado_ts = Some(c);
                }
            }
            b"ts" => {
                if let Aux::Char(c) = aux {
                    t.minimap_ts = Some(c);
                }
            }
            b"qs" => t.qs = aux_float(&aux),
            b"pt" => {
                if let Some(v) = aux_int(&aux) {
                    t.pt = v.clamp(i64::from(PT_NO_ANCHOR), i64::from(i32::MAX)) as i32;
                }
            }
            b"SA" => t.has_sa = true,
            _ => {}
        }
    }
    t
}

/// Alignment geometry from the CIGAR.
#[derive(Debug, Default, Clone)]
pub struct Geometry {
    pub start: u32,
    pub end: u32,
    pub clip_left: u32,
    pub clip_right: u32,
    /// Reference bases covered by M/=/X/D: the exonic footprint.
    pub exonic_len: u32,
    /// Introns (N operations) as half-open reference intervals.
    pub junctions: Vec<(u32, u32)>,
}

const OP_M: u32 = 0;
const OP_D: u32 = 2;
const OP_N: u32 = 3;
const OP_S: u32 = 4;
const OP_H: u32 = 5;
const OP_EQ: u32 = 7;
const OP_X: u32 = 8;

impl Geometry {
    /// Parse `cigar` for an alignment starting at `pos`, reusing `self`.
    pub fn parse(&mut self, pos: u32, cigar: &[u32]) {
        self.junctions.clear();
        self.clip_left = 0;
        self.clip_right = 0;
        self.exonic_len = 0;
        self.start = pos;
        let mut p = pos;
        let mut seen_ref = false;
        for &c in cigar {
            let (op, len) = (c & 0xf, c >> 4);
            match op {
                OP_M | OP_EQ | OP_X | OP_D => {
                    p += len;
                    self.exonic_len += len;
                    seen_ref = true;
                    self.clip_right = 0;
                }
                OP_N => {
                    self.junctions.push((p, p + len));
                    p += len;
                    seen_ref = true;
                    self.clip_right = 0;
                }
                OP_S | OP_H => {
                    if seen_ref {
                        self.clip_right += len;
                    } else {
                        self.clip_left += len;
                    }
                }
                _ => {}
            }
        }
        self.end = p.max(pos + 1);
    }

    /// Unaligned bases on the transcript's 5' (UMI) side.
    pub fn clip5(&self, strand: Strand) -> u32 {
        match strand {
            Strand::Minus => self.clip_right,
            _ => self.clip_left,
        }
    }

    /// The molecule's 5' end: where strand switching attached the UMI.
    pub fn anchor(&self, strand: Strand) -> u32 {
        match strand {
            Strand::Minus => self.end - 1,
            _ => self.start,
        }
    }

    /// If the first aligned exon is short, where the alignment would begin
    /// had the aligner soft-clipped that exon instead: the acceptor of the
    /// first junction. Returns (position, first exon length).
    pub fn alt_anchor(&self, strand: Strand, max_exon: u32) -> Option<(u32, u32)> {
        let (pos, exon) = match strand {
            Strand::Plus => {
                let &(lo, hi) = self.junctions.first()?;
                (hi, lo - self.start)
            }
            Strand::Minus => {
                // Exon 1 lies right of the last intron; a read that clipped it
                // ends on the base before that intron.
                let &(lo, hi) = self.junctions.last()?;
                (lo.saturating_sub(1), self.end - hi)
            }
            Strand::Unknown => return None,
        };
        (exon <= max_exon).then_some((pos, exon))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cig(ops: &[(u32, u32)]) -> Vec<u32> {
        ops.iter().map(|&(op, len)| (len << 4) | op).collect()
    }

    #[test]
    fn geometry_of_spliced_alignment() {
        let mut g = Geometry::default();
        g.parse(
            100,
            &cig(&[
                (OP_S, 3),
                (OP_M, 20),
                (OP_N, 500),
                (OP_M, 50),
                (OP_D, 2),
                (OP_M, 10),
                (OP_S, 40),
            ]),
        );
        assert_eq!(g.start, 100);
        assert_eq!(g.end, 100 + 20 + 500 + 50 + 2 + 10);
        assert_eq!(g.exonic_len, 82);
        assert_eq!(g.junctions, vec![(120, 620)]);
        assert_eq!((g.clip_left, g.clip_right), (3, 40));
        assert_eq!(g.anchor(Strand::Plus), 100);
        assert_eq!(g.anchor(Strand::Minus), g.end - 1);
        assert_eq!(g.clip5(Strand::Plus), 3);
        assert_eq!(g.clip5(Strand::Minus), 40);
        assert_eq!(g.alt_anchor(Strand::Plus, 30), Some((620, 20)));
        assert_eq!(g.alt_anchor(Strand::Plus, 10), None);
        // Minus strand: exon 1 is the last 62 aligned reference bases.
        assert_eq!(g.alt_anchor(Strand::Minus, 70), Some((119, 62)));
    }

    #[test]
    fn minus_strand_alt_anchor_matches_clipped_read_anchor() {
        // Minus-strand transcript: exon 2 = [1000, 1100), exon 1 = [1600, 1615).
        let mut full = Geometry::default();
        full.parse(1000, &cig(&[(OP_M, 100), (OP_N, 500), (OP_M, 15), (OP_S, 3)]));
        let mut clipped = Geometry::default();
        clipped.parse(1000, &cig(&[(OP_M, 100), (OP_S, 18)]));
        let (alt, exon) = full.alt_anchor(Strand::Minus, 40).unwrap();
        assert_eq!(exon, 15);
        assert_eq!(alt, clipped.anchor(Strand::Minus));
    }

    #[test]
    fn strand_from_read_orientation() {
        assert_eq!(genomic_strand(b'+', false), Strand::Plus);
        assert_eq!(genomic_strand(b'+', true), Strand::Minus);
        assert_eq!(genomic_strand(b'-', false), Strand::Minus);
        assert_eq!(genomic_strand(b'-', true), Strand::Plus);
        assert_eq!(genomic_strand(b'?', true), Strand::Unknown);
    }
}
