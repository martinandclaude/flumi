//! Decoding and comparison of the SQK-PCB114.24 UMI carried in dorado's `RX:Z`.
//!
//! Dorado finds the UMI by aligning the pattern `TTTVVVVTTVVVVTTVVVVTTVVVVTTT`
//! (V = A/C/G) against a 40-nt window behind the SSP primer with edlib in infix
//! mode, accepts the hit with up to 5 edits, and stores the aligned read
//! substring in `RX` (reverse-complemented for reverse reads). `RX` therefore
//! carries substitutions *and* indels, anywhere in the tag.
//!
//! Decoding re-aligns `RX` globally to the pattern and reads the four V blocks
//! off that alignment. Anchor errors of any kind (substitution, insertion,
//! deletion, several at once) cost nothing: the blocks come out intact. Block
//! errors are kept rather than guessed away: a block holding a `T` or with a
//! length other than 4 cannot occur in a real UMI, so the UMI is marked
//! *damaged* and the clustering attaches it to the molecule it came from
//! instead of letting it found one of its own.
//!
//! A decoded UMI is packed into a `u64`: four 15-bit blocks (3-bit length plus
//! up to six 2-bit bases) and a damaged flag.

use std::fmt;

pub const PATTERN: &[u8; PLEN] = b"TTTVVVVTTVVVVTTVVVVTTVVVVTTT";
const PLEN: usize = 28;
pub const N_BLOCKS: usize = 4;
pub const BLOCK_LEN: usize = 4;
pub const MAX_BLOCK_LEN: usize = 6;
const MAX_RX_LEN: usize = 48;

const BLOCK_BITS: u32 = 15;
const BLOCK_MASK: u64 = (1 << BLOCK_BITS) - 1;
const DAMAGED: u64 = 1 << 60;

/// Cost of one edit. On top of it each operation carries a small tie-break
/// ranking explanations with the same number of edits by how likely they are
/// on nanopore reads, and by how much they risk: a T-run gaining or losing a
/// T (+0) is the commonest error; a block base gained, lost or substituted
/// (+1) keeps every informative base; a substitution inside an anchor (+3)
/// silently discards one. Keeping bases matters when an interpretation is
/// ambiguous: the kept base leaves a damaged block one edit from the truth,
/// whereas discarding it can shift a block into a clean-looking UMI two or
/// three edits away.
const EDIT: u16 = 32;

/// dorado keeps the co-optimal edlib hit that ends first, so a T of the
/// trailing anchor can be missing from `RX` without any sequencing error
/// (edlib then extends the start as far as the score allows, so leading slack
/// shows up as extra bases instead). Deleting outer anchor Ts costs only a
/// tie-break, at either end.
const EDGE_DEL: u16 = 1;
const LEAD: usize = 3;
const TRAIL: usize = PLEN - 3;

#[derive(Clone, Copy, PartialEq, Eq)]
enum Gap {
    /// Between two V positions of block k.
    Block(u8),
    /// Between two anchor positions, or outside the pattern.
    Anchor,
    /// Between an anchor and block k.
    Boundary(u8),
}

const fn is_v(j: usize) -> bool {
    PATTERN[j] == b'V'
}

const fn block_of(j: usize) -> u8 {
    // Blocks occupy 3..7, 9..13, 15..19, 21..25.
    ((j - 3) / 6) as u8
}

const fn gap(g: usize) -> Gap {
    let left_v = g > 0 && is_v(g - 1);
    let right_v = g < PLEN && is_v(g);
    match (left_v, right_v) {
        (true, true) => Gap::Block(block_of(g)),
        (false, false) => Gap::Anchor,
        (true, false) => Gap::Boundary(block_of(g - 1)),
        (false, true) => Gap::Boundary(block_of(g)),
    }
}

#[inline]
fn sub_cost(b: u8, j: usize) -> u16 {
    match (b, is_v(j)) {
        (b'T', false) | (b'A' | b'C' | b'G', true) => 0,
        (_, true) => EDIT + 1,
        (_, false) => EDIT + 3,
    }
}

/// Inserting base `b` into gap `g`. A `T` belongs in an anchor and anything
/// else in a block, so equal-edit alignments keep informative bases and shed
/// anchor slack.
#[inline]
fn ins_cost(b: u8, g: usize) -> u16 {
    EDIT + match (gap(g), b == b'T') {
        (Gap::Anchor | Gap::Boundary(_), true) => 0,
        (Gap::Block(_) | Gap::Boundary(_), false) => 1,
        (Gap::Block(_), true) | (Gap::Anchor, false) => 2,
    }
}

/// Deleting pattern position `j` after `i` of `n` RX bases: losing a T from
/// an anchor run is cheaper than losing a block base, and an outer anchor T
/// beyond either end of `RX` was never an error at all.
#[inline]
fn del_cost(j: usize, i: usize, n: usize) -> u16 {
    if (i == 0 && j < LEAD) || (i == n && j >= TRAIL) {
        EDGE_DEL
    } else {
        EDIT + u16::from(is_v(j))
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DecodeError {
    Empty,
    /// More edits against the pattern than allowed.
    TooManyEdits,
    /// A block came out empty or longer than `MAX_BLOCK_LEN`.
    BlockLength,
    /// A block holds a base outside ACGT.
    BadBase,
}

/// A decoded UMI.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Umi {
    pub key: u64,
    /// Edits between `RX` and the pattern (anchor and block errors).
    pub edits: u8,
}

impl Umi {
    pub fn is_damaged(&self) -> bool {
        is_damaged(self.key)
    }
}

#[inline]
fn base_code(b: u8) -> Option<u16> {
    match b {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

fn pack(blocks: &[[u8; MAX_BLOCK_LEN]; N_BLOCKS], lens: &[usize; N_BLOCKS]) -> Result<u64, DecodeError> {
    let mut key = 0u64;
    let mut damaged = false;
    for k in 0..N_BLOCKS {
        let len = lens[k];
        if len == 0 || len > MAX_BLOCK_LEN {
            return Err(DecodeError::BlockLength);
        }
        let mut code = (len as u16) << 12;
        for (i, &b) in blocks[k][..len].iter().enumerate() {
            let c = base_code(b).ok_or(DecodeError::BadBase)?;
            damaged |= c == 3;
            code |= c << (2 * i);
        }
        damaged |= len != BLOCK_LEN;
        key |= u64::from(code) << (BLOCK_BITS * k as u32);
    }
    if damaged {
        key |= DAMAGED;
    }
    Ok(key)
}

/// Decode an `RX` value into its four informative blocks.
///
/// `max_edits` bounds the edits between `RX` and the pattern; dorado itself
/// accepts at most 5.
pub fn decode(rx: &[u8], max_edits: u8) -> Result<Umi, DecodeError> {
    let n = rx.len();
    if n == 0 {
        return Err(DecodeError::Empty);
    }
    if n > MAX_RX_LEN.min(PLEN + usize::from(max_edits)) || n + LEAD + 3 + usize::from(max_edits) < PLEN {
        return Err(DecodeError::TooManyEdits);
    }
    let mut seq = [0u8; MAX_RX_LEN];
    for (d, s) in seq.iter_mut().zip(rx) {
        *d = s.to_ascii_uppercase();
    }
    let seq = &seq[..n];

    let mut blocks = [[0u8; MAX_BLOCK_LEN]; N_BLOCKS];
    let mut lens = [0usize; N_BLOCKS];

    // Fast path: a clean 28-nt tag.
    if n == PLEN && (0..PLEN).all(|j| sub_cost(seq[j], j) == 0) {
        for j in (0..PLEN).filter(|&j| is_v(j)) {
            let k = block_of(j) as usize;
            blocks[k][lens[k]] = seq[j];
            lens[k] += 1;
        }
        return Ok(Umi {
            key: pack(&blocks, &lens)?,
            edits: 0,
        });
    }

    let mut cost = [[0u16; PLEN + 1]; MAX_RX_LEN + 1];
    for j in 1..=PLEN {
        cost[0][j] = cost[0][j - 1] + del_cost(j - 1, 0, n);
    }
    for i in 1..=n {
        let b = seq[i - 1];
        cost[i][0] = cost[i - 1][0] + ins_cost(b, 0);
        for j in 1..=PLEN {
            let diag = cost[i - 1][j - 1] + sub_cost(b, j - 1);
            let ins = cost[i - 1][j] + ins_cost(b, j);
            let del = cost[i][j - 1] + del_cost(j - 1, i, n);
            cost[i][j] = diag.min(ins).min(del);
        }
    }

    // Traceback, filling blocks back to front.
    let (mut i, mut j) = (n, PLEN);
    let mut edits = 0u32;
    let mut push = |k: u8, b: u8, blocks: &mut [[u8; MAX_BLOCK_LEN]; N_BLOCKS]| -> Result<(), DecodeError> {
        let k = k as usize;
        if lens[k] == MAX_BLOCK_LEN {
            return Err(DecodeError::BlockLength);
        }
        blocks[k][lens[k]] = b;
        lens[k] += 1;
        Ok(())
    };
    while i > 0 || j > 0 {
        let c = cost[i][j];
        if i > 0 && j > 0 && c == cost[i - 1][j - 1] + sub_cost(seq[i - 1], j - 1) {
            if sub_cost(seq[i - 1], j - 1) != 0 {
                edits += 1;
            }
            if is_v(j - 1) {
                push(block_of(j - 1), seq[i - 1], &mut blocks)?;
            }
            i -= 1;
            j -= 1;
        } else if i > 0 && c == cost[i - 1][j] + ins_cost(seq[i - 1], j) {
            edits += 1;
            let b = seq[i - 1];
            match gap(j) {
                Gap::Block(k) => push(k, b, &mut blocks)?,
                Gap::Boundary(k) if b != b'T' => push(k, b, &mut blocks)?,
                _ => {}
            }
            i -= 1;
        } else {
            if del_cost(j - 1, i, n) != EDGE_DEL {
                edits += 1;
            }
            j -= 1;
        }
        if edits > u32::from(max_edits) {
            return Err(DecodeError::TooManyEdits);
        }
    }
    for k in 0..N_BLOCKS {
        blocks[k][..lens[k]].reverse();
    }
    Ok(Umi {
        key: pack(&blocks, &lens)?,
        edits: edits as u8,
    })
}

#[inline]
pub fn is_damaged(key: u64) -> bool {
    key & DAMAGED != 0
}

#[inline]
pub fn block(key: u64, k: usize) -> u16 {
    ((key >> (BLOCK_BITS * k as u32)) & BLOCK_MASK) as u16
}

#[inline]
fn block_len(code: u16) -> usize {
    usize::from(code >> 12)
}

/// Edit distance between two blocks, exact up to `budget`; above it only a
/// lower bound is guaranteed.
fn block_distance(x: u16, y: u16, budget: u32) -> u32 {
    let (lx, ly) = (block_len(x), block_len(y));
    if lx == ly {
        let z = (x ^ y) & 0x0FFF;
        let hamming = ((z | (z >> 1)) & 0x0555).count_ones();
        // Equal lengths: one edit can only be a substitution, and an
        // insertion plus a deletion can undercut a Hamming distance >= 2.
        if hamming <= 1 || budget < 2 {
            return hamming.min(2);
        }
    } else if lx.abs_diff(ly) as u32 > budget {
        return lx.abs_diff(ly) as u32;
    }
    let base = |code: u16, i: usize| (code >> (2 * i)) & 3;
    let mut prev = [0u32; MAX_BLOCK_LEN + 1];
    let mut cur = [0u32; MAX_BLOCK_LEN + 1];
    for (j, p) in prev.iter_mut().enumerate().take(ly + 1) {
        *p = j as u32;
    }
    for i in 1..=lx {
        cur[0] = i as u32;
        for j in 1..=ly {
            let s = prev[j - 1] + u32::from(base(x, i - 1) != base(y, j - 1));
            cur[j] = s.min(prev[j] + 1).min(cur[j - 1] + 1);
        }
        std::mem::swap(&mut prev, &mut cur);
    }
    prev[ly]
}

/// Blockwise edit distance: the sum of per-block edit distances. Anchors
/// resynchronise the tag, so an indel never shifts bases across blocks.
/// Returns early once the running sum exceeds `limit`.
pub fn distance(a: u64, b: u64, limit: u32) -> u32 {
    let mut d = 0;
    for k in 0..N_BLOCKS {
        let (x, y) = (block(a, k), block(b, k));
        if x != y {
            d += block_distance(x, y, limit - d);
            if d > limit {
                break;
            }
        }
    }
    d
}

/// Bitmask of the blocks two UMIs share exactly.
#[inline]
pub fn shared_blocks(a: u64, b: u64) -> u8 {
    (0..N_BLOCKS).fold(0u8, |m, k| if block(a, k) == block(b, k) { m | (1 << k) } else { m })
}

/// Key made of the blocks selected by `keep`; UMIs within blockwise distance
/// `N_BLOCKS - popcount(keep)` share at least one such key exactly.
#[inline]
pub fn masked_key(key: u64, keep: u8) -> u64 {
    let mut m = u64::from(keep) << 60;
    for k in 0..N_BLOCKS {
        if keep & (1 << k) != 0 {
            m |= u64::from(block(key, k)) << (BLOCK_BITS * k as u32);
        }
    }
    m
}

/// All block masks keeping exactly `n_keep` blocks.
pub fn keep_masks(n_keep: usize) -> Vec<u8> {
    (0u8..16).filter(|m| m.count_ones() as usize == n_keep).collect()
}

/// The lowest `n_keep` blocks of `shared`: the one mask under which a pair is
/// examined, so a pair sharing several masks is not compared twice.
#[inline]
pub fn canonical_keep(shared: u8, n_keep: usize) -> u8 {
    let mut out = 0u8;
    let mut left = n_keep;
    for k in 0..N_BLOCKS {
        if left == 0 {
            break;
        }
        if shared & (1 << k) != 0 {
            out |= 1 << k;
            left -= 1;
        }
    }
    out
}

/// The informative bases, blocks concatenated: 16 nt for an undamaged UMI.
pub struct Display(pub u64);

impl fmt::Display for Display {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        const BASES: &[u8; 4] = b"ACGT";
        let mut buf = [0u8; N_BLOCKS * MAX_BLOCK_LEN];
        let mut n = 0;
        for k in 0..N_BLOCKS {
            let code = block(self.0, k);
            for i in 0..block_len(code) {
                buf[n] = BASES[usize::from((code >> (2 * i)) & 3)];
                n += 1;
            }
        }
        f.write_str(std::str::from_utf8(&buf[..n]).expect("ASCII"))
    }
}

pub fn to_string(key: u64) -> String {
    Display(key).to_string()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Visible damage: extra or missing bases plus `T`s summed over the blocks.
    fn damage(key: u64) -> u32 {
        (0..N_BLOCKS)
            .map(|k| {
                let code = block(key, k);
                let len = block_len(code);
                let ts = (0..len).filter(|i| (code >> (2 * i)) & 3 == 3).count();
                (len.abs_diff(BLOCK_LEN) + ts) as u32
            })
            .sum()
    }

    fn rx(anchors: [&str; 5], blocks: [&str; 4]) -> String {
        let mut s = String::new();
        for k in 0..4 {
            s.push_str(anchors[k]);
            s.push_str(blocks[k]);
        }
        s.push_str(anchors[4]);
        s
    }

    const A: [&str; 5] = ["TTT", "TT", "TT", "TT", "TTT"];
    const B: [&str; 4] = ["ACGA", "CGAC", "GACG", "ACGC"];

    fn dec(s: &str) -> Umi {
        decode(s.as_bytes(), 5).unwrap_or_else(|e| panic!("{s}: {e:?}"))
    }

    #[test]
    fn clean_tag_takes_fast_path() {
        let u = dec(&rx(A, B));
        assert_eq!(to_string(u.key), "ACGACGACGACGACGC");
        assert_eq!(u.edits, 0);
        assert!(!u.is_damaged());
    }

    #[test]
    fn anchor_errors_cost_nothing() {
        let clean = dec(&rx(A, B)).key;
        for anchors in [
            ["TT", "TT", "TT", "TT", "TTT"],   // deletion
            ["TTT", "TTT", "TT", "TT", "TTT"], // insertion
            ["TTT", "TA", "TT", "TT", "TTT"],  // substitution, 28 nt
            ["TTT", "TAT", "TT", "TT", "TT"],  // substitution + deletion
            ["TT", "T", "TTT", "TT", "TTTT"],  // several at once
            ["TTTT", "TTT", "TTT", "TTT", "TTTT"],
        ] {
            let u = dec(&rx(anchors, B));
            assert_eq!(u.key, clean, "{anchors:?}");
            assert!(!u.is_damaged());
        }
    }

    #[test]
    fn compensating_indels_do_not_shift_blocks() {
        // One extra anchor T and one missing anchor T keep the length at 28 but
        // shift the middle blocks: a fixed-offset parse would read garbage.
        let s = rx(["TTT", "TTT", "TT", "T", "TTT"], B);
        assert_eq!(s.len(), 28);
        assert_eq!(dec(&s).key, dec(&rx(A, B)).key);
    }

    #[test]
    fn block_substitution_to_v_stays_valid() {
        let u = dec(&rx(A, ["ACGC", "CGAC", "GACG", "ACGC"]));
        assert!(!u.is_damaged());
        assert_eq!(distance(u.key, dec(&rx(A, B)).key, 4), 1);
    }

    #[test]
    fn block_damage_is_kept_and_flagged() {
        let clean = dec(&rx(A, B)).key;
        for blocks in [
            ["ACTA", "CGAC", "GACG", "ACGC"],  // V -> T
            ["ACG", "CGAC", "GACG", "ACGC"],   // deletion
            ["ACGAA", "CGAC", "GACG", "ACGC"], // insertion
            ["AGA", "CGAC", "GACG", "ACGC"],   // deletion inside block
            ["ACGA", "CGAC", "GACG", "ACGCG"], // insertion before final anchor
        ] {
            let u = dec(&rx(A, blocks));
            assert!(u.is_damaged(), "{blocks:?}");
            assert_eq!(distance(u.key, clean, 4), 1, "{blocks:?}");
            assert_eq!(damage(u.key), 1, "{blocks:?}");
        }
    }

    #[test]
    fn inserted_t_at_block_edge_goes_to_anchor() {
        // ACGA + T + TT: the extra T is anchor slack, not block content.
        let u = dec(&rx(["TTT", "TTT", "TT", "TT", "TTT"], B));
        assert!(!u.is_damaged());
    }

    #[test]
    fn garbage_is_rejected() {
        assert_eq!(decode(b"", 5), Err(DecodeError::Empty));
        assert!(decode(b"ACGTACGTACGTACGTACGTACGTACGT", 5).is_err());
        assert!(decode(&[b'T'; 28], 5).is_err());
        assert!(decode(rx(A, B).as_bytes().get(..20).unwrap(), 5).is_err());
        let mut long = rx(A, B);
        long.push_str("ACGAACGA");
        assert!(decode(long.as_bytes(), 5).is_err());
        // N inside a block cannot be packed.
        assert_eq!(
            decode(rx(A, ["ACNA", "CGAC", "GACG", "ACGC"]).as_bytes(), 5),
            Err(DecodeError::BadBase)
        );
    }

    #[test]
    fn max_edits_is_enforced() {
        let s = rx(["TTT", "T", "T", "T", "TTT"], B); // 3 internal anchor deletions
        assert_eq!(decode(s.as_bytes(), 3).unwrap().edits, 3);
        assert_eq!(decode(s.as_bytes(), 2), Err(DecodeError::TooManyEdits));
    }

    #[test]
    fn outer_anchor_ts_missing_from_rx_are_not_errors() {
        // dorado keeps the co-optimal hit that ends first, so trailing Ts can
        // be cut; leading ones only go missing through real deletions.
        let u = dec(&rx(["T", "TT", "TT", "TT", "T"], B));
        assert_eq!(u.key, dec(&rx(A, B)).key);
        assert_eq!(u.edits, 0);
    }

    #[test]
    fn block_insertion_with_trimmed_trailing_t_stays_one_edit_away() {
        // Observed in simulation: true last block AGCA gained an A (AGACA) and
        // dorado's hit stopped before the final anchor T. Reading AGAC as the
        // block would be a clean UMI two edits from the truth.
        let truth = dec(&rx(A, ["ACGG", "CAAG", "CGGG", "AGCA"])).key;
        let u = dec("TTTACGGTTCAAGTTCGGGTTAGACATT");
        assert!(u.is_damaged());
        assert_eq!(distance(u.key, truth, 4), 1);
        // Same RX if the truth was AGAC with an anchor substitution: still 1.
        let alt = dec(&rx(A, ["ACGG", "CAAG", "CGGG", "AGAC"])).key;
        assert_eq!(distance(u.key, alt, 4), 1);
    }

    #[test]
    fn lowercase_is_accepted() {
        assert_eq!(dec(&rx(A, B).to_lowercase()).key, dec(&rx(A, B)).key);
    }

    #[test]
    fn distance_is_symmetric_and_blockwise() {
        let a = dec(&rx(A, B)).key;
        let b = dec(&rx(A, ["CCGA", "CGAC", "GACG", "ACGA"])).key;
        let c = dec(&rx(A, ["ACG", "CGAC", "GACG", "ACGC"])).key;
        assert_eq!(distance(a, b, 8), 2);
        assert_eq!(distance(b, a, 8), 2);
        assert_eq!(distance(a, c, 8), 1);
        assert_eq!(distance(b, c, 8), 3);
        assert_eq!(distance(a, a, 8), 0);
        // Levenshtein inside a block: ACGA vs CGAC is 2, not Hamming 4.
        let d = dec(&rx(A, ["CGAC", "CGAC", "GACG", "ACGC"])).key;
        assert_eq!(distance(a, d, 8), 2);
    }

    #[test]
    fn masked_keys_find_all_neighbours() {
        let a = dec(&rx(A, B)).key;
        let b = dec(&rx(A, ["ACGA", "CGAC", "GAGG", "ACGC"])).key;
        let shared = shared_blocks(a, b);
        assert_eq!(shared, 0b1011);
        let keep = canonical_keep(shared, 3);
        assert_eq!(masked_key(a, keep), masked_key(b, keep));
        assert_eq!(keep_masks(3).len(), 4);
        assert_eq!(keep_masks(2).len(), 6);
        assert_eq!(canonical_keep(0b1111, 2), 0b0011);
    }
}
