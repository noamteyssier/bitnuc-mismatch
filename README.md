# bitnuc-mismatch

Create unambiguous one-off mismatch hash tables from bitnuc scalars.

This library adapts my work in [`disambiseq`](https://crates.io/crates/disambiseq) to operating in 2-bit space.

Note that this is for sequences that are represented as [`bitnuc`](https://crates.io/crates/bitnuc) scalars.
By definition this limits the sequences to a maximum length of 32 nucleotides.

Future work will include support for longer sequences.

## Usage

This is a library for generating unambiguous one-off mismatches for bitnuc scalars.
The library provides a function to generate all possible one-off mismatches for a given bitnuc scalar.
It also provides a function to build a mismatch table, which maps mismatches to their parent sequences.
The library is designed to be used in the context of generating mismatches for a set of parent sequences while avoiding ambiguous mismatches.
Ambiguous mismatches are mismatches that are within the one-off distance of multiple parent sequences.

This builds on the [`bitnuc`](https://crates.io/crates/bitnuc) library, which provides functions for converting nucleotide sequences to bitnuc scalars.

## Example

```rust

use bitnuc_mismatch::build_mismatch_table;

// Define a set of parent sequences
let parent_sequences = vec![
    b"ACTG",
    b"ACCG",
];

// Convert the parent sequences to bitnuc scalars
let parent_scalars: Vec<u64> = parent_sequences
    .into_iter()
    .map(|seq| bitnuc::as_2bit(seq).unwrap())
    .collect();

// Build a mismatch table
let mismatch_table = build_mismatch_table(&parent_scalars, 4).unwrap();

// Test some expected mismatches
let gcta = bitnuc::as_2bit(b"GCTG").unwrap();
assert_eq!(mismatch_table.get(&gcta), Some(&parent_scalars[0]));

// Validate that unexpected mismatches are not present
let acgg = bitnuc::as_2bit(b"ACGG").unwrap();
assert!(mismatch_table.get(&acgg).is_none());
```
