//! Choosing the read that represents a molecule.
//!
//! PCR copies of one molecule share one true structure, so any disagreement
//! between them is technical. Each junction is judged by majority among the
//! copies whose alignment spans it: *supported* when at least two copies and
//! a strict majority carry it, *contradicted* when a strict majority of the
//! spanning copies do not. A candidate is penalised for every contradicted
//! junction it carries, every supported junction it spans but lacks, and (in
//! strict mode) every junction no other copy confirms unless the annotation
//! knows it. Among the candidates with fewest disagreements the read that
//! reached the poly(A) anchor wins, then one without a supplementary
//! alignment, then the largest exonic footprint.

use std::cmp::Ordering;

use rustc_hash::FxHashMap;

use crate::bundle::{BundleRead, ChainStore};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RepMode {
    /// Junction consensus, then completeness, then length.
    Consensus,
    /// Exonic footprint alone.
    Longest,
}

#[derive(Clone, Debug)]
pub struct RepParams {
    pub mode: RepMode,
    /// Penalise junctions seen in one copy only (unless annotated).
    pub strict: bool,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RepChoice {
    /// Index into the bundle's reads.
    pub rep: u32,
    /// Copies sharing the representative's exact junction chain.
    pub same_chain: u32,
    /// Disagreements between the representative and the copies' consensus.
    pub disagreements: u32,
    /// Consensus picked a different read than length alone would have.
    pub switched: bool,
}

fn longer(a: &BundleRead, b: &BundleRead) -> Ordering {
    a.exonic_len
        .cmp(&b.exonic_len)
        .then(a.qs.total_cmp(&b.qs))
        .then(b.ordinal.cmp(&a.ordinal))
}

pub fn choose(
    reads: &[BundleRead],
    members: &[u32],
    chains: &ChainStore,
    annotated: &dyn Fn((u32, u32)) -> bool,
    p: &RepParams,
) -> RepChoice {
    let read = |i: u32| &reads[i as usize];
    let longest = *members
        .iter()
        .max_by(|&&a, &&b| longer(read(a), read(b)))
        .expect("non-empty molecule");
    let same_chain = |rep: u32| members.iter().filter(|&&m| read(m).chain == read(rep).chain).count() as u32;

    if members.len() == 1 || p.mode == RepMode::Longest {
        return RepChoice {
            rep: longest,
            same_chain: same_chain(longest),
            disagreements: 0,
            switched: false,
        };
    }

    let mut chain_count: FxHashMap<u32, u32> = FxHashMap::default();
    for &m in members {
        *chain_count.entry(read(m).chain).or_default() += 1;
    }
    let mut support: FxHashMap<(u32, u32), u32> = FxHashMap::default();
    for (&c, &k) in &chain_count {
        for &j in chains.get(c) {
            *support.entry(j).or_default() += k;
        }
    }

    // Junction -> (supported, bad-to-carry).
    let mut verdict: FxHashMap<(u32, u32), (bool, bool)> = FxHashMap::default();
    let mut supported: Vec<(u32, u32)> = Vec::new();
    for (&j, &s) in &support {
        let spanning = members
            .iter()
            .filter(|&&m| read(m).start <= j.0 && read(m).end >= j.1)
            .count() as u32;
        let is_supported = s >= 2 && 2 * s > spanning;
        let contradicted = 2 * s < spanning;
        let unconfirmed = s == 1 && p.strict && !annotated(j);
        verdict.insert(j, (is_supported, contradicted || unconfirmed));
        if is_supported {
            supported.push(j);
        }
    }

    let mut chain_bad: FxHashMap<u32, u32> = FxHashMap::default();
    for &c in chain_count.keys() {
        let bad = chains.get(c).iter().filter(|j| verdict[j].1).count() as u32;
        chain_bad.insert(c, bad);
    }
    let disagreements = |r: &BundleRead| {
        let own = chains.get(r.chain);
        let missing = supported
            .iter()
            .filter(|j| r.start <= j.0 && r.end >= j.1 && own.binary_search(j).is_err())
            .count() as u32;
        chain_bad[&r.chain] + missing
    };

    let mut best = members[0];
    let mut best_bad = disagreements(read(best));
    for &m in &members[1..] {
        let (r, b) = (read(m), read(best));
        let bad = disagreements(r);
        let ord = best_bad
            .cmp(&bad)
            .then((r.pt >= 0).cmp(&(b.pt >= 0)))
            .then(b.has_sa.cmp(&r.has_sa))
            .then(longer(r, b));
        if ord == Ordering::Greater {
            best = m;
            best_bad = bad;
        }
    }
    RepChoice {
        rep: best,
        same_chain: same_chain(best),
        disagreements: best_bad,
        switched: best != longest,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::features::NO_POS;

    fn read(ordinal: u64, start: u32, end: u32, chain: u32, exonic_len: u32) -> BundleRead {
        BundleRead {
            ordinal,
            umi: 0,
            start,
            end,
            anchor: start,
            alt_anchor: NO_POS,
            alt_exon: 0,
            clip5: 0,
            chain,
            exonic_len,
            qs: 20.0,
            pt: -2,
            has_sa: false,
        }
    }

    const CONS: RepParams = RepParams {
        mode: RepMode::Consensus,
        strict: true,
    };

    fn none(_: (u32, u32)) -> bool {
        false
    }

    #[test]
    fn minority_junction_loses_even_when_longest() {
        let mut chains = ChainStore::default();
        let good = chains.intern(&[(200, 300)]);
        let bad = chains.intern(&[(150, 200)]);
        let reads = vec![
            read(0, 100, 400, good, 200),
            read(1, 100, 400, good, 200),
            read(2, 100, 400, good, 200),
            read(3, 100, 500, bad, 350),
        ];
        let c = choose(&reads, &[0, 1, 2, 3], &chains, &none, &CONS);
        assert_ne!(c.rep, 3);
        assert_eq!(c.same_chain, 3);
        assert_eq!(c.disagreements, 0);
        assert!(c.switched);
        let l = choose(
            &reads,
            &[0, 1, 2, 3],
            &chains,
            &none,
            &RepParams {
                mode: RepMode::Longest,
                strict: true,
            },
        );
        assert_eq!(l.rep, 3);
    }

    #[test]
    fn read_missing_a_majority_junction_is_penalised() {
        let mut chains = ChainStore::default();
        let spliced = chains.intern(&[(200, 300)]);
        // Read 2 claims the intron is retained: longest footprint, but wrong.
        let reads = vec![
            read(0, 100, 400, spliced, 200),
            read(1, 100, 400, spliced, 200),
            read(2, 100, 400, 0, 300),
        ];
        let c = choose(&reads, &[0, 1, 2], &chains, &none, &CONS);
        assert_ne!(c.rep, 2);
    }

    #[test]
    fn longest_wins_when_structures_agree() {
        let mut chains = ChainStore::default();
        let ch = chains.intern(&[(200, 300)]);
        let reads = vec![
            read(0, 100, 400, ch, 200),
            read(1, 100, 400, ch, 200),
            read(2, 100, 600, ch, 400),
        ];
        assert_eq!(choose(&reads, &[0, 1, 2], &chains, &none, &CONS).rep, 2);
    }

    #[test]
    fn full_length_read_preferred_over_longer_truncated_one() {
        let mut chains = ChainStore::default();
        let ch = chains.intern(&[(200, 300)]);
        let mut fl = read(0, 100, 400, ch, 200);
        fl.pt = 60;
        let trunc = read(1, 100, 400, ch, 210);
        assert_eq!(choose(&[fl, trunc], &[0, 1], &chains, &none, &CONS).rep, 0);
    }

    #[test]
    fn annotated_singleton_junction_is_not_penalised() {
        let mut chains = ChainStore::default();
        let a = chains.intern(&[(200, 300), (350, 450)]);
        let b = chains.intern(&[(200, 300)]);
        // Read 1 does not span (350, 450), so nothing contradicts it.
        let reads = vec![read(0, 100, 600, a, 300), read(1, 100, 320, b, 150)];
        assert_eq!(choose(&reads, &[0, 1], &chains, &none, &CONS).rep, 1);
        let known = |j: (u32, u32)| j == (350, 450);
        assert_eq!(choose(&reads, &[0, 1], &chains, &known, &CONS).rep, 0);
        let lenient = RepParams {
            mode: RepMode::Consensus,
            strict: false,
        };
        assert_eq!(choose(&reads, &[0, 1], &chains, &none, &lenient).rep, 0);
    }
}
