//! Strand-specific loci and the compact per-read record kept while one is open.
//!
//! Every PCR copy of a molecule carries the UMI *and* the molecule's 5' end
//! (dorado only finds the UMI if the read reaches the SSP primer, which sits
//! on the 5' end), so all copies of a molecule overlap on the same transcript
//! strand. A bundle is a maximal run of such overlapping alignments; it is
//! closed once the coordinate-sorted stream moves past its end, and never
//! needs a gene annotation.

use rustc_hash::FxHashMap;

use crate::features::Strand;

#[derive(Clone, Debug)]
pub struct BundleRead {
    /// Record index in the input file, the key for the second pass.
    pub ordinal: u64,
    pub umi: u64,
    pub start: u32,
    pub end: u32,
    /// Genomic position of the molecule's 5' end as seen by this read.
    pub anchor: u32,
    /// `anchor` if the aligner had clipped a short first exon, else `NO_POS`.
    pub alt_anchor: u32,
    pub alt_exon: u16,
    pub clip5: u16,
    /// Interned junction chain.
    pub chain: u32,
    pub exonic_len: u32,
    pub qs: f32,
    pub pt: i32,
    pub has_sa: bool,
}

/// Junction chains interned per bundle: PCR copies share them, so each read
/// only stores an id. Id 0 is the unspliced chain.
pub struct ChainStore {
    index: FxHashMap<Box<[(u32, u32)]>, u32>,
    chains: Vec<Box<[(u32, u32)]>>,
}

impl Default for ChainStore {
    fn default() -> Self {
        Self {
            index: FxHashMap::default(),
            chains: vec![Box::new([])],
        }
    }
}

impl ChainStore {
    pub fn intern(&mut self, junctions: &[(u32, u32)]) -> u32 {
        if junctions.is_empty() {
            return 0;
        }
        if let Some(&id) = self.index.get(junctions) {
            return id;
        }
        let id = self.chains.len() as u32;
        let chain: Box<[(u32, u32)]> = junctions.into();
        self.chains.push(chain.clone());
        self.index.insert(chain, id);
        id
    }

    pub fn get(&self, id: u32) -> &[(u32, u32)] {
        &self.chains[id as usize]
    }
}

pub struct Bundle {
    pub tid: u32,
    pub strand: Strand,
    pub reads: Vec<BundleRead>,
    pub chains: ChainStore,
    pub max_end: u32,
}

impl Bundle {
    pub fn new(tid: u32, strand: Strand) -> Self {
        Self {
            tid,
            strand,
            reads: Vec::new(),
            chains: ChainStore::default(),
            max_end: 0,
        }
    }
}
