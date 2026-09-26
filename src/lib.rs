//! flumi: collapse ONT PCR-cDNA reads to one read per original molecule using
//! the pre-PCR UMI that dorado reports in `RX:Z`.

pub mod annotation;
pub mod bundle;
pub mod cluster;
pub mod consensus;
pub mod dedup;
pub mod features;
pub mod stats;
pub mod umi;
