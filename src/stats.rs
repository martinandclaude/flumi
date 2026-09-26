//! Run statistics and their TSV rendering.

use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::Result;

use crate::cluster::ClusterStats;

#[derive(Default, Debug, Clone)]
pub struct Stats {
    pub total_records: u64,
    pub unmapped: u64,
    pub secondary_supplementary: u64,
    pub low_mapq: u64,
    pub low_qs: u64,
    pub no_umi: u64,
    pub umi_undecodable: u64,
    pub chimeric_umi_end: u64,
    pub umi_clean: u64,
    pub umi_anchor_repaired: u64,
    pub umi_damaged: u64,
    pub strand_dorado: u64,
    pub strand_minimap2: u64,
    pub strand_unknown: u64,
    pub reads_used: u64,
    pub bundles: u64,
    pub max_bundle_reads: u64,
    pub cluster: ClusterStats,
    pub molecules: u64,
    pub size_hist: BTreeMap<u32, u64>,
    pub rep_switched: u64,
    pub rep_disagreements: u64,
    pub full_length_annot: u64,
    pub polya_anchor: u64,
    pub polya_missing: u64,
    pub full_length_both: u64,
    pub records_written: u64,
}

fn pct(num: u64, den: u64) -> f64 {
    if den == 0 {
        0.0
    } else {
        (10000.0 * num as f64 / den as f64).round() / 100.0
    }
}

fn percentile(hist: &BTreeMap<u32, u64>, q: f64) -> u32 {
    let n: u64 = hist.values().sum();
    if n == 0 {
        return 0;
    }
    let rank = ((n - 1) as f64 * q).round() as u64;
    let mut seen = 0;
    for (&size, &k) in hist {
        seen += k;
        if rank < seen {
            return size;
        }
    }
    0
}

impl Stats {
    pub fn rows(&self) -> Vec<(&'static str, String)> {
        let umi_reads = self.umi_clean + self.umi_anchor_repaired + self.umi_damaged;
        let dup = self.reads_used - self.molecules.min(self.reads_used);
        let singletons = self.size_hist.get(&1).copied().unwrap_or(0);
        let polya_called = self.molecules - self.polya_missing;
        let c = &self.cluster;
        vec![
            ("total_records", self.total_records.to_string()),
            ("dropped_unmapped", self.unmapped.to_string()),
            (
                "dropped_secondary_supplementary",
                self.secondary_supplementary.to_string(),
            ),
            ("dropped_low_mapq", self.low_mapq.to_string()),
            ("dropped_low_qs", self.low_qs.to_string()),
            ("dropped_no_RX", self.no_umi.to_string()),
            ("dropped_umi_undecodable", self.umi_undecodable.to_string()),
            ("dropped_chimeric_umi_end", self.chimeric_umi_end.to_string()),
            ("umi_clean", self.umi_clean.to_string()),
            ("umi_anchor_repaired", self.umi_anchor_repaired.to_string()),
            ("umi_block_damaged", self.umi_damaged.to_string()),
            ("umi_block_damaged_pct", pct(self.umi_damaged, umi_reads).to_string()),
            ("strand_from_TS", self.strand_dorado.to_string()),
            ("strand_from_ts", self.strand_minimap2.to_string()),
            ("strand_unknown", self.strand_unknown.to_string()),
            ("reads_used", self.reads_used.to_string()),
            ("bundles", self.bundles.to_string()),
            ("max_bundle_reads", self.max_bundle_reads.to_string()),
            ("umi_nodes", c.nodes.to_string()),
            ("umi_nodes_damaged", c.damaged_nodes.to_string()),
            ("umi_position_splits", c.position_splits.to_string()),
            ("short_exon_bridges", c.bridged.to_string()),
            ("directional_merges", c.valid_merged.to_string()),
            ("density_gated_edges", c.far_edges.to_string()),
            ("damaged_to_valid", c.damaged_to_valid.to_string()),
            ("damaged_attached_beyond_1_edit", c.damaged_far.to_string()),
            ("damaged_to_damaged", c.damaged_to_damaged.to_string()),
            ("damaged_seeds", c.damaged_seeds.to_string()),
            ("molecules_out", self.molecules.to_string()),
            ("pcr_duplicate_reads", dup.to_string()),
            ("pcr_duplicate_pct", pct(dup, self.reads_used).to_string()),
            ("singleton_molecules", singletons.to_string()),
            ("cluster_size_p50", percentile(&self.size_hist, 0.5).to_string()),
            ("cluster_size_p90", percentile(&self.size_hist, 0.9).to_string()),
            (
                "cluster_size_max",
                self.size_hist.keys().next_back().copied().unwrap_or(0).to_string(),
            ),
            ("rep_switched_by_consensus", self.rep_switched.to_string()),
            ("rep_disagreements", self.rep_disagreements.to_string()),
            ("full_length_molecules", self.full_length_annot.to_string()),
            (
                "full_length_pct",
                pct(self.full_length_annot, self.molecules).to_string(),
            ),
            ("polya_anchor_molecules", self.polya_anchor.to_string()),
            ("polya_anchor_pct", pct(self.polya_anchor, polya_called).to_string()),
            ("polya_missing_tag_molecules", self.polya_missing.to_string()),
            ("full_length_both_signals", self.full_length_both.to_string()),
            ("records_written", self.records_written.to_string()),
        ]
    }

    pub fn write_tsv(&self, path: &Path) -> Result<()> {
        let mut w = BufWriter::new(File::create(path)?);
        writeln!(w, "metric\tvalue")?;
        for (k, v) in self.rows() {
            writeln!(w, "{k}\t{v}")?;
        }
        w.flush()?;
        Ok(())
    }
}
