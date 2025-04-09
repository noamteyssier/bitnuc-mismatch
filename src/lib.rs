//! # bitnuc-mismatch
//!
//! This is a library for generating unambiguous one-off mismatches for bitnuc scalars.
//! The library provides a function to generate all possible one-off mismatches for a given bitnuc scalar.
//! The library also provides a function to build a mismatch table, which maps mismatches to their parent sequences.
//! The library is designed to be used in the context of generating mismatches for a set of parent sequences while avoiding ambiguous mismatches.
//! Ambiguous mismatches are mismatches that are within the one-off distance of multiple parent sequences.
//!
//! Note that parent sequences will be members of the mismatch table.
//!
//! ## Example
//! ```rust
//!
//! use bitnuc_mismatch::build_mismatch_table;
//!
//! let parent_sequences = vec![
//!     b"ACTG",
//!     b"ACCG",
//! ];
//! let parent_scalars: Vec<u64> = parent_sequences
//!     .into_iter()
//!     .map(|seq| bitnuc::as_2bit(seq).unwrap())
//!     .collect();
//!
//! let mismatch_table = build_mismatch_table(&parent_scalars, 4).unwrap();
//!
//! // Test some expected mismatches
//! let gcta = bitnuc::as_2bit(b"GCTG").unwrap();
//! assert_eq!(mismatch_table.get(&gcta), Some(&parent_scalars[0]));
//!
//! // Validate that unexpected mismatches are not present
//! let acgg = bitnuc::as_2bit(b"ACGG").unwrap();
//! assert!(mismatch_table.get(&acgg).is_none());
//!
//! // Validate that parent sequences are members of the mismatch table
//! assert!(mismatch_table.contains_key(&parent_scalars[0]));
//! assert!(mismatch_table.contains_key(&parent_scalars[1]));
//! ```

mod error;
mod utils;

pub use error::MismatchError;
pub use utils::{build_mismatch_table, generate_mismatches};

/// Type alias for a mismatch table mapping mismatches to their parent sequences
pub type MismatchTable = hashbrown::HashMap<BitSeqMut, BitSeq>;

/// Type alias for a bitnuc scalar sequence (child or parent sequence)
pub type BitSeqMut = u64;

/// Type alias for a bitnuc scalar sequence (parent sequence)
pub type BitSeq = u64;
