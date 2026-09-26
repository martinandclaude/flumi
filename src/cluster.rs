//! Grouping the reads of one bundle into molecules.
//!
//! 1. **Nodes.** Reads with the same decoded UMI are split by 5' anchor:
//!    single linkage within `pos_tol`, plus a bridge between a read that
//!    aligned a short first exon and one that soft-clipped it (the second
//!    read's anchor sits on the first read's exon-2 acceptor). Two molecules
//!    that drew the same UMI but start in different places stay apart.
//! 2. **Damaged UMIs** (a `T` or a wrong length inside a block) cannot be real
//!    UMIs. Each is credited to its nearest compatible undamaged node, so the
//!    directional comparison below sees every read of a molecule, and follows
//!    that node into its molecule.
//! 3. **Directional merging** of undamaged UMIs (UMI-tools' rule,
//!    `n_hub >= 2 n_sub - 1`) between positionally compatible nodes.
//!
//! Merges beyond one edit are **density-gated**. A UMI two or three edits from
//! a hub is either a read with two or three errors or a different molecule
//! that happens to lie close; which is likelier depends on how many distinct
//! UMIs share the 5' neighbourhood. With `m` of them drawn from the 3^16 UMI
//! space, about `m * N_d / 3^16` lie within distance `d` of any UMI by chance
//! (`N_d` = 33, 513, 4993 for d = 1, 2, 3). A merge at distance `d >= 2` is
//! allowed only while that expectation stays within a budget: a strict one
//! for undamaged UMIs (most of which are real molecules), a looser one for
//! damaged UMIs (which are certainly errors of *some* molecule). In an
//! ordinary gene this recovers reads with two errors; in a gene with tens of
//! thousands of molecules at one TSS it falls back to single edits.
//!
//! Neighbour search is exact: UMIs within blockwise distance `d` share at
//! least `4 - d` blocks verbatim, so candidates come from sorted keys built
//! from every `(4 - d)`-block subset.

use std::cmp::Reverse;

use rustc_hash::FxHashMap;

use crate::bundle::BundleRead;
use crate::features::NO_POS;
use crate::umi;

/// Soft clip a read needs before it can stand for a clipped first exon.
const MIN_BRIDGE_CLIP: u16 = 5;

#[derive(Clone, Debug)]
pub struct Params {
    /// Blockwise edits between undamaged UMIs of one molecule.
    pub umi_dist: u32,
    /// Blockwise edits between a damaged UMI and its molecule.
    pub damaged_dist: u32,
    /// Largest 5'-anchor difference between reads of one molecule.
    pub pos_tol: u32,
    /// Compare 5' anchors at all (off for unknown-strand bundles).
    pub use_position: bool,
    /// Chance-neighbour budget for undamaged merges beyond one edit.
    pub valid_budget: f64,
    /// Chance-neighbour budget for damaged attachments beyond one edit.
    pub damaged_budget: f64,
}

/// Undamaged UMIs within blockwise substitution distance d of a UMI, itself
/// included: 1, 1 + 16*2, + C(16,2)*2^2, + C(16,3)*2^3.
const NEIGHBOURHOOD: [f64; 4] = [1.0, 33.0, 513.0, 4993.0];
/// 3^16 possible undamaged UMIs.
const UMI_SPACE: f64 = 43_046_721.0;

/// Expected number of UMIs within distance `d` of a given one by chance among
/// `m` distinct UMIs.
fn chance(d: u32, m: u32) -> f64 {
    f64::from(m) * NEIGHBOURHOOD[d.min(3) as usize] / UMI_SPACE
}

/// Largest merge distance, at most `max`, whose chance stays within `budget`.
fn reach(max: u32, m: u32, budget: f64) -> u32 {
    let mut d = max.min(1);
    while d < max.min(3) && chance(d + 1, m) <= budget {
        d += 1;
    }
    d
}

#[derive(Default, Debug, Clone)]
pub struct ClusterStats {
    pub nodes: u64,
    pub damaged_nodes: u64,
    /// Extra nodes from splitting one exact UMI across distant 5' anchors.
    pub position_splits: u64,
    /// Read pairs joined through a clipped short first exon.
    pub bridged: u64,
    pub valid_merged: u64,
    /// Directional edges beyond one edit that the density gate admitted.
    pub far_edges: u64,
    pub damaged_to_valid: u64,
    /// Damaged nodes attached beyond one edit.
    pub damaged_far: u64,
    pub damaged_to_damaged: u64,
    pub damaged_seeds: u64,
}

impl ClusterStats {
    pub fn add(&mut self, o: &ClusterStats) {
        self.nodes += o.nodes;
        self.damaged_nodes += o.damaged_nodes;
        self.position_splits += o.position_splits;
        self.bridged += o.bridged;
        self.valid_merged += o.valid_merged;
        self.far_edges += o.far_edges;
        self.damaged_to_valid += o.damaged_to_valid;
        self.damaged_far += o.damaged_far;
        self.damaged_to_damaged += o.damaged_to_damaged;
        self.damaged_seeds += o.damaged_seeds;
    }
}

#[derive(Debug)]
pub struct Molecule {
    /// UMI of the node that seeded the molecule (its most abundant).
    pub umi: u64,
    /// Indices into the bundle's reads.
    pub reads: Vec<u32>,
}

struct Node {
    umi: u64,
    count: u32,
    amin: u32,
    amax: u32,
    max_clip5: u16,
    first: u32,
    len: u32,
    alt_first: u32,
    alt_len: u32,
}

#[derive(Default)]
struct UnionFind {
    parent: Vec<u32>,
}

impl UnionFind {
    fn reset(&mut self, n: usize) {
        self.parent.clear();
        self.parent.extend(0..n as u32);
    }

    fn find(&mut self, mut x: u32) -> u32 {
        while self.parent[x as usize] != x {
            let p = self.parent[x as usize];
            self.parent[x as usize] = self.parent[p as usize];
            x = p;
        }
        x
    }

    fn union(&mut self, a: u32, b: u32) -> bool {
        let (ra, rb) = (self.find(a), self.find(b));
        if ra == rb {
            return false;
        }
        let (lo, hi) = if ra < rb { (ra, rb) } else { (rb, ra) };
        self.parent[hi as usize] = lo;
        true
    }
}

#[inline]
fn range_gap(a0: u32, a1: u32, b0: u32, b1: u32) -> u32 {
    b0.saturating_sub(a1).max(a0.saturating_sub(b1))
}

#[inline]
fn min_bridge_clip(exon: u16) -> u16 {
    (exon / 2).max(MIN_BRIDGE_CLIP)
}

struct Graph<'a> {
    nodes: &'a [Node],
    alts: &'a [(u32, u16)],
    p: &'a Params,
}

impl Graph<'_> {
    fn bridged(&self, a: &Node, b: &Node) -> bool {
        let alts = &self.alts[a.alt_first as usize..(a.alt_first + a.alt_len) as usize];
        alts.iter().any(|&(pos, exon)| {
            range_gap(pos, pos, b.amin, b.amax) <= self.p.pos_tol && b.max_clip5 >= min_bridge_clip(exon)
        })
    }

    fn compatible(&self, u: u32, v: u32) -> bool {
        if !self.p.use_position {
            return true;
        }
        let (a, b) = (&self.nodes[u as usize], &self.nodes[v as usize]);
        range_gap(a.amin, a.amax, b.amin, b.amax) <= self.p.pos_tol || self.bridged(a, b) || self.bridged(b, a)
    }
}

fn push_node(
    reads: &[BundleRead],
    run: impl Iterator<Item = u32>,
    nodes: &mut Vec<Node>,
    members: &mut Vec<u32>,
    alts: &mut Vec<(u32, u16)>,
) {
    let first = members.len();
    let alt_first = alts.len();
    let mut node = Node {
        umi: 0,
        count: 0,
        amin: u32::MAX,
        amax: 0,
        max_clip5: 0,
        first: first as u32,
        len: 0,
        alt_first: alt_first as u32,
        alt_len: 0,
    };
    for i in run {
        let r = &reads[i as usize];
        node.umi = r.umi;
        node.amin = node.amin.min(r.anchor);
        node.amax = node.amax.max(r.anchor);
        node.max_clip5 = node.max_clip5.max(r.clip5);
        if r.alt_anchor != NO_POS {
            alts.push((r.alt_anchor, r.alt_exon));
        }
        members.push(i);
    }
    node.len = (members.len() - first) as u32;
    node.count = node.len;
    let node_alts = &mut alts[alt_first..];
    node_alts.sort_unstable();
    let mut w = 0;
    for r in 0..node_alts.len() {
        if w == 0 || node_alts[r].0 != node_alts[w - 1].0 {
            node_alts[w] = node_alts[r];
            w += 1;
        }
    }
    alts.truncate(alt_first + w);
    node.alt_len = w as u32;
    nodes.push(node);
}

/// Build nodes: exact UMI, split by 5' anchor when positions are in use.
fn build_nodes(reads: &[BundleRead], p: &Params, st: &mut ClusterStats) -> (Vec<Node>, Vec<u32>, Vec<(u32, u16)>) {
    let n = reads.len();
    let mut order: Vec<u32> = (0..n as u32).collect();
    order.sort_unstable_by_key(|&i| {
        let r = &reads[i as usize];
        (r.umi, r.anchor, r.ordinal)
    });

    let mut nodes = Vec::new();
    let mut members = Vec::with_capacity(n);
    let mut alts = Vec::new();
    let mut uf = UnionFind::default();
    let mut comp_of: Vec<u32> = Vec::new();
    let mut i = 0;
    while i < n {
        let umi = reads[order[i] as usize].umi;
        let mut j = i + 1;
        while j < n && reads[order[j] as usize].umi == umi {
            j += 1;
        }
        let run = &order[i..j];
        i = j;
        if !p.use_position || run.len() == 1 {
            push_node(reads, run.iter().copied(), &mut nodes, &mut members, &mut alts);
            continue;
        }
        let m = run.len();
        uf.reset(m);
        let anchor = |k: usize| reads[run[k] as usize].anchor;
        for k in 1..m {
            if anchor(k) - anchor(k - 1) <= p.pos_tol {
                uf.union(k as u32 - 1, k as u32);
            }
        }
        for k in 0..m {
            let r = &reads[run[k] as usize];
            if r.alt_anchor == NO_POS {
                continue;
            }
            let (lo, hi) = (
                r.alt_anchor.saturating_sub(p.pos_tol),
                r.alt_anchor.saturating_add(p.pos_tol),
            );
            let from = run.partition_point(|&x| reads[x as usize].anchor < lo);
            for (q, &x) in run.iter().enumerate().skip(from) {
                let s = &reads[x as usize];
                if s.anchor > hi {
                    break;
                }
                if s.clip5 >= min_bridge_clip(r.alt_exon) && uf.union(k as u32, q as u32) {
                    st.bridged += 1;
                }
            }
        }
        // Number components by first appearance, then bucket members in one
        // pass (a common UMI can sit at thousands of 5' ends).
        comp_of.clear();
        comp_of.resize(m, u32::MAX);
        let mut n_comp = 0u32;
        let mut comp: Vec<u32> = Vec::with_capacity(m);
        for k in 0..m {
            let root = uf.find(k as u32) as usize;
            if comp_of[root] == u32::MAX {
                comp_of[root] = n_comp;
                n_comp += 1;
            }
            comp.push(comp_of[root]);
        }
        st.position_splits += u64::from(n_comp - 1);
        let mut by_comp: Vec<(u32, u32)> = comp.iter().zip(run).map(|(&c, &r)| (c, r)).collect();
        by_comp.sort_unstable();
        for group in by_comp.chunk_by(|a, b| a.0 == b.0) {
            push_node(
                reads,
                group.iter().map(|&(_, r)| r),
                &mut nodes,
                &mut members,
                &mut alts,
            );
        }
    }
    (nodes, members, alts)
}

/// Sorted `(masked key, node)` entries over every `n_keep`-block subset.
fn block_index(nodes: &[Node], ids: &[u32], masks: &[u8]) -> Vec<(u64, u32)> {
    let mut entries = Vec::with_capacity(ids.len() * masks.len());
    for &v in ids {
        for &m in masks {
            entries.push((umi::masked_key(nodes[v as usize].umi, m), v));
        }
    }
    entries.sort_unstable();
    entries
}

struct ClusterInfo {
    seed: u32,
    total: u64,
}

/// Distinct undamaged UMIs whose 5' anchors fall within `pos_tol` of each
/// node's: the population a chance neighbour would come from.
fn local_density(nodes: &[Node], valid: &[u32], p: &Params) -> Vec<u32> {
    if !p.use_position {
        return vec![valid.len() as u32; nodes.len()];
    }
    let mut amins: Vec<u32> = valid.iter().map(|&v| nodes[v as usize].amin).collect();
    amins.sort_unstable();
    nodes
        .iter()
        .map(|n| {
            let (lo, hi) = (n.amin.saturating_sub(p.pos_tol), n.amax.saturating_add(p.pos_tol));
            (amins.partition_point(|&a| a <= hi) - amins.partition_point(|&a| a < lo)) as u32
        })
        .collect()
}

/// Candidates sharing `n_keep` blocks with `key`, each pair visited once.
fn candidates<'a>(
    index: &'a [(u64, u32)],
    key: u64,
    n_keep: usize,
    nodes: &'a [Node],
) -> impl Iterator<Item = u32> + 'a {
    umi::keep_masks(n_keep).into_iter().flat_map(move |m| {
        let mk = umi::masked_key(key, m);
        let lo = index.partition_point(|e| e.0 < mk);
        index[lo..]
            .iter()
            .take_while(move |e| e.0 == mk)
            .filter_map(move |&(_, v)| {
                let vu = nodes[v as usize].umi;
                (umi::canonical_keep(umi::shared_blocks(key, vu), n_keep) == m).then_some(v)
            })
    })
}

pub fn cluster(reads: &[BundleRead], p: &Params, st: &mut ClusterStats) -> Vec<Molecule> {
    if reads.is_empty() {
        return Vec::new();
    }
    let (nodes, members, alts) = build_nodes(reads, p, st);
    let graph = Graph {
        nodes: &nodes,
        alts: &alts,
        p,
    };
    st.nodes += nodes.len() as u64;

    let (mut valid, mut damaged): (Vec<u32>, Vec<u32>) =
        (0..nodes.len() as u32).partition(|&v| !umi::is_damaged(nodes[v as usize].umi));
    st.damaged_nodes += damaged.len() as u64;
    damaged.sort_unstable_by_key(|&v| {
        let n = &nodes[v as usize];
        (Reverse(n.count), n.umi, n.amin)
    });
    let density = local_density(&nodes, &valid, p);
    let damaged_reach = |d: u32| reach(p.damaged_dist, density[d as usize], p.damaged_budget);

    // Credit each damaged node to its nearest compatible undamaged node.
    const NONE: u32 = u32::MAX;
    let mut weight: Vec<u64> = nodes.iter().map(|n| u64::from(n.count)).collect();
    let mut parent = vec![NONE; nodes.len()];
    let widest = damaged.iter().map(|&d| damaged_reach(d)).max().unwrap_or(0);
    if widest > 0 && !valid.is_empty() {
        let n_keep = umi::N_BLOCKS - widest as usize;
        let vindex = block_index(&nodes, &valid, &umi::keep_masks(n_keep));
        for &d in &damaged {
            let (du, limit) = (nodes[d as usize].umi, damaged_reach(d));
            let mut best: Option<(u32, Reverse<u32>, u32)> = None;
            for v in candidates(&vindex, du, n_keep, &nodes) {
                let dist = umi::distance(du, nodes[v as usize].umi, limit);
                if dist > limit || !graph.compatible(d, v) {
                    continue;
                }
                let cand = (dist, Reverse(nodes[v as usize].count), v);
                if best.is_none_or(|b| cand < b) {
                    best = Some(cand);
                }
            }
            if let Some((dist, _, v)) = best {
                parent[d as usize] = v;
                weight[v as usize] += u64::from(nodes[d as usize].count);
                st.damaged_far += u64::from(dist >= 2);
            }
        }
    }
    valid.sort_unstable_by_key(|&v| {
        let n = &nodes[v as usize];
        (Reverse(weight[v as usize]), n.umi, n.amin)
    });

    // Directional edges between undamaged nodes. The index only needs keys
    // as wide as the largest distance the density gate allows anywhere here.
    let k = valid
        .iter()
        .map(|&v| reach(p.umi_dist, density[v as usize], p.valid_budget))
        .max()
        .unwrap_or(0)
        .min(3);
    let n_keep = umi::N_BLOCKS - k as usize;
    let mut edges: Vec<(u32, u32)> = Vec::new();
    if valid.len() >= 2 && k > 0 {
        let entries = block_index(&nodes, &valid, &umi::keep_masks(n_keep));
        for group in entries.chunk_by(|a, b| a.0 == b.0) {
            if group.len() < 2 {
                continue;
            }
            let keep = (group[0].0 >> 60) as u8;
            for x in 0..group.len() {
                for y in x + 1..group.len() {
                    let (u, v) = (group[x].1, group[y].1);
                    let (a, b) = (&nodes[u as usize], &nodes[v as usize]);
                    if umi::canonical_keep(umi::shared_blocks(a.umi, b.umi), n_keep) != keep {
                        continue;
                    }
                    let dist = umi::distance(a.umi, b.umi, k);
                    if dist > k || !graph.compatible(u, v) {
                        continue;
                    }
                    let (wa, wb) = (weight[u as usize], weight[v as usize]);
                    if dist >= 2 {
                        // A read with two errors is one read: the far node
                        // must be small, and the neighbourhood sparse.
                        let m = density[u as usize].max(density[v as usize]);
                        if dist > reach(p.umi_dist, m, p.valid_budget) || wa.min(wb) > 2 {
                            continue;
                        }
                    }
                    let before = edges.len();
                    if wa >= 2 * wb - 1 {
                        edges.push((u, v));
                    }
                    if wb >= 2 * wa - 1 {
                        edges.push((v, u));
                    }
                    if dist >= 2 && edges.len() > before {
                        st.far_edges += 1;
                    }
                }
            }
        }
    }
    edges.sort_unstable();
    edges.dedup();
    let mut offsets = vec![0u32; nodes.len() + 1];
    for &(u, _) in &edges {
        offsets[u as usize + 1] += 1;
    }
    for i in 0..nodes.len() {
        offsets[i + 1] += offsets[i];
    }
    let neighbours = |u: u32| &edges[offsets[u as usize] as usize..offsets[u as usize + 1] as usize];

    let mut cluster_of = vec![NONE; nodes.len()];
    let mut clusters: Vec<ClusterInfo> = Vec::new();
    let mut stack = Vec::new();
    for &s in &valid {
        if cluster_of[s as usize] != NONE {
            continue;
        }
        let cid = clusters.len() as u32;
        clusters.push(ClusterInfo { seed: s, total: 0 });
        cluster_of[s as usize] = cid;
        stack.push(s);
        while let Some(x) = stack.pop() {
            clusters[cid as usize].total += weight[x as usize];
            if x != s {
                st.valid_merged += 1;
            }
            for &(_, y) in neighbours(x) {
                if cluster_of[y as usize] == NONE {
                    cluster_of[y as usize] = cid;
                    stack.push(y);
                }
            }
        }
    }

    // Damaged nodes follow their parent; the rest join the nearest damaged
    // node within reach or found a molecule of their own.
    let mut dindex: FxHashMap<u64, Vec<u32>> = FxHashMap::default();
    let single = umi::keep_masks(1);
    for &d in &damaged {
        let v = parent[d as usize];
        if v != NONE {
            cluster_of[d as usize] = cluster_of[v as usize];
            st.damaged_to_valid += 1;
            for &m in &single {
                dindex
                    .entry(umi::masked_key(nodes[d as usize].umi, m))
                    .or_default()
                    .push(d);
            }
        }
    }
    for &d in damaged.iter().filter(|&&d| parent[d as usize] == NONE) {
        let (du, limit) = (nodes[d as usize].umi, damaged_reach(d));
        let mut best: Option<(u32, Reverse<u64>, Reverse<u32>, u32)> = None;
        for &m in &single {
            let Some(list) = dindex.get(&umi::masked_key(du, m)) else {
                continue;
            };
            for &v in list {
                let vu = nodes[v as usize].umi;
                if umi::canonical_keep(umi::shared_blocks(du, vu), 1) != m {
                    continue;
                }
                let dist = umi::distance(du, vu, limit);
                if dist > limit || !graph.compatible(d, v) {
                    continue;
                }
                let c = &clusters[cluster_of[v as usize] as usize];
                let cand = (dist, Reverse(c.total), Reverse(nodes[v as usize].count), v);
                if best.is_none_or(|b| cand < b) {
                    best = Some(cand);
                }
            }
        }
        match best {
            Some((_, _, _, v)) => {
                let cid = cluster_of[v as usize];
                cluster_of[d as usize] = cid;
                clusters[cid as usize].total += u64::from(nodes[d as usize].count);
                st.damaged_to_damaged += 1;
            }
            None => {
                cluster_of[d as usize] = clusters.len() as u32;
                clusters.push(ClusterInfo {
                    seed: d,
                    total: u64::from(nodes[d as usize].count),
                });
                st.damaged_seeds += 1;
            }
        }
        for &m in &single {
            dindex.entry(umi::masked_key(du, m)).or_default().push(d);
        }
    }

    let mut molecules: Vec<Molecule> = clusters
        .iter()
        .map(|c| Molecule {
            umi: nodes[c.seed as usize].umi,
            reads: Vec::new(),
        })
        .collect();
    for (v, node) in nodes.iter().enumerate() {
        let span = node.first as usize..(node.first + node.len) as usize;
        molecules[cluster_of[v] as usize]
            .reads
            .extend_from_slice(&members[span]);
    }
    molecules
}

#[cfg(test)]
mod tests {
    use super::*;

    fn umi_of(blocks: [&str; 4]) -> u64 {
        let rx = format!("TTT{}TT{}TT{}TT{}TTT", blocks[0], blocks[1], blocks[2], blocks[3]);
        umi::decode(rx.as_bytes(), 5).unwrap().key
    }

    fn read(ordinal: u64, umi: u64, anchor: u32) -> BundleRead {
        BundleRead {
            ordinal,
            umi,
            start: anchor,
            end: anchor + 500,
            anchor,
            alt_anchor: NO_POS,
            alt_exon: 0,
            clip5: 3,
            chain: 0,
            exonic_len: 500,
            qs: 20.0,
            pt: -2,
            has_sa: false,
        }
    }

    fn params() -> Params {
        Params {
            umi_dist: 1,
            damaged_dist: 2,
            pos_tol: 20,
            use_position: true,
            valid_budget: 0.01,
            damaged_budget: 0.5,
        }
    }

    fn sizes(mols: &[Molecule]) -> Vec<usize> {
        let mut s: Vec<usize> = mols.iter().map(|m| m.reads.len()).collect();
        s.sort_unstable();
        s
    }

    const B: [&str; 4] = ["ACGA", "CGAC", "GACG", "ACGC"];

    #[test]
    fn directional_merges_skewed_but_not_equal_counts() {
        let a = umi_of(B);
        let b = umi_of(["ACGC", "CGAC", "GACG", "ACGC"]);
        let mut reads: Vec<BundleRead> = (0..3).map(|i| read(i, a, 100)).collect();
        reads.push(read(3, b, 101));
        let mut st = ClusterStats::default();
        assert_eq!(sizes(&cluster(&reads, &params(), &mut st)), vec![4]);

        let mut reads: Vec<BundleRead> = (0..2).map(|i| read(i, a, 100)).collect();
        reads.extend((2..4).map(|i| read(i, b, 100)));
        assert_eq!(sizes(&cluster(&reads, &params(), &mut st)), vec![2, 2]);
    }

    #[test]
    fn same_umi_far_apart_is_two_molecules() {
        let a = umi_of(B);
        let reads = vec![read(0, a, 100), read(1, a, 105), read(2, a, 5000)];
        let mut st = ClusterStats::default();
        assert_eq!(sizes(&cluster(&reads, &params(), &mut st)), vec![1, 2]);
        assert_eq!(st.position_splits, 1);
        let p = Params {
            use_position: false,
            ..params()
        };
        assert_eq!(sizes(&cluster(&reads, &p, &mut st)), vec![3]);
    }

    #[test]
    fn neighbour_umi_far_apart_is_not_merged() {
        let a = umi_of(B);
        let b = umi_of(["ACGC", "CGAC", "GACG", "ACGC"]);
        let reads = vec![read(0, a, 100), read(1, a, 100), read(2, a, 100), read(3, b, 900)];
        let mut st = ClusterStats::default();
        assert_eq!(sizes(&cluster(&reads, &params(), &mut st)), vec![1, 3]);
    }

    #[test]
    fn damaged_umis_join_their_source() {
        let a = umi_of(B);
        let far = umi_of(["CCGC", "GGAC", "GACG", "ACGC"]);
        let del = umi_of(["ACG", "CGAC", "GACG", "ACGC"]);
        let two = umi_of(["ACTA", "CGAC", "GACG", "ACGCC"]);
        let reads = vec![read(0, a, 100), read(1, far, 100), read(2, del, 102), read(3, two, 99)];
        let mut st = ClusterStats::default();
        let mols = cluster(&reads, &params(), &mut st);
        assert_eq!(sizes(&mols), vec![1, 3]);
        let big = mols.iter().find(|m| m.reads.len() == 3).unwrap();
        assert_eq!(big.umi, a);
        assert_eq!(st.damaged_to_valid, 2);
    }

    #[test]
    fn damaged_umis_without_source_group_together() {
        let d1 = umi_of(["ACTA", "CGAC", "GACG", "ACGC"]);
        let d2 = umi_of(["ACGA", "CGC", "GACG", "ACGC"]);
        let reads = vec![read(0, d1, 100), read(1, d2, 100)];
        let mut st = ClusterStats::default();
        assert_eq!(sizes(&cluster(&reads, &params(), &mut st)), vec![2]);
        assert_eq!((st.damaged_seeds, st.damaged_to_damaged), (1, 1));
    }

    #[test]
    fn directional_does_not_bridge_two_hubs() {
        let hub1 = umi_of(["CAAA", "AAAA", "AAAA", "AAAA"]);
        let hub2 = umi_of(["ACAA", "AAAA", "AAAA", "AAAA"]);
        let low = umi_of(["AAAA", "AAAA", "AAAA", "AAAA"]);
        let mut reads = Vec::new();
        let mut o = 0;
        for (u, n) in [(hub1, 10), (hub2, 10), (low, 3)] {
            for _ in 0..n {
                reads.push(read(o, u, 100));
                o += 1;
            }
        }
        let mut st = ClusterStats::default();
        assert_eq!(sizes(&cluster(&reads, &params(), &mut st)), vec![10, 13]);
    }

    #[test]
    fn clipped_short_first_exon_is_bridged() {
        let a = umi_of(B);
        // Read 0 aligned a 15-nt first exon at 100 and reaches exon 2 at 900.
        let mut full = read(0, a, 100);
        full.alt_anchor = 900;
        full.alt_exon = 15;
        // Read 1 soft-clipped that exon: it starts on exon 2.
        let mut clipped = read(1, a, 900);
        clipped.clip5 = 18;
        // Read 2 genuinely starts at exon 2 (no clip): a different molecule.
        let other = read(2, umi_of(["ACGC", "CGAC", "GACG", "ACGC"]), 900);
        let mut st = ClusterStats::default();
        let mols = cluster(&[full.clone(), clipped.clone()], &params(), &mut st);
        assert_eq!(sizes(&mols), vec![2]);
        assert_eq!(st.bridged, 1);
        let mols = cluster(&[full, other], &params(), &mut st);
        assert_eq!(sizes(&mols), vec![1, 1]);
    }
}
