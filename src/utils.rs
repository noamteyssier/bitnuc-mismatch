use crate::AmbiguousMismatchTable;

use super::{MismatchError, MismatchTable};
use hashbrown::{HashMap, HashSet};

/// Generates all one-off mismatches for a given bitnuc scalar sequence
///
/// This function will generate 3 mismatches for each position in the sequence,
/// as each base can be mutated to any of the other 3 possible nucleotides.
///
/// # Arguments
/// * sequence: The bitnuc scalar (2-bit encoded nucleotide sequence)
/// * slen: Length of the sequence in bases
/// * mm_buffer: Buffer to store generated mismatches
///
/// # Returns
/// * Ok(()) if successful
/// * Err(MismatchError) if sequence length > 32
pub fn generate_mismatches(
    sequence: u64,
    slen: usize,
    mm_buffer: &mut Vec<u64>,
) -> Result<(), MismatchError> {
    if slen > 32 {
        return Err(MismatchError::InvalidSequenceLength(slen));
    }

    // Clear the buffer
    mm_buffer.clear();

    // For each position in the sequence
    for pos in 0..slen {
        // Get the bit position
        let bit_pos = pos * 2;

        // Get the current base at this position (2 bits)
        let base = (sequence >> bit_pos) & 0b11;

        // Generate the other 3 possible bases
        for other_base in 0..4u64 {
            if other_base != base {
                // Create a mask to clear the current base
                let clear_mask = !(0b11u64 << bit_pos);

                // Clear the current base and set the new base
                let mutated = (sequence & clear_mask) | (other_base << bit_pos);

                mm_buffer.push(mutated);
            }
        }
    }

    Ok(())
}

/// Builds a table mapping all possible one-off mismatches to their parent sequences
///
/// For each input sequence, this function:
/// 1. Adds the sequence itself as its own parent
/// 2. Generates all one-off mismatches and maps them to the parent sequence
/// 3. Handles collisions by removing ambiguous mappings
///
/// # Arguments
/// * sequences: Slice of bitnuc scalar sequences
/// * slen: Length of each sequence in bases
///
/// # Returns
/// * Ok(HashMap) mapping mismatches to their unambiguous parent
/// * Err(MismatchError) if sequence length > 32
pub fn build_mismatch_table(
    sequences: &[u64],
    slen: usize,
) -> Result<MismatchTable, MismatchError> {
    if slen > 32 {
        return Err(MismatchError::InvalidSequenceLength(slen));
    }

    let mut mm_table = HashMap::new();
    let mut ambiguous = HashSet::new();
    let mut mm_buffer = Vec::with_capacity(slen * 3);

    // Initialize the MM table with all the sequences mapping to themselves.
    // This is the case where there are no mismatches.
    for &seq in sequences {
        mm_table.insert(seq, seq);
        ambiguous.insert(seq);
    }

    // Process each sequence
    for &parent in sequences {
        // generate mismatches for the parent sequence
        generate_mismatches(parent, slen, &mut mm_buffer)?;

        // process each mismatch
        for &mm_seq in &mm_buffer {
            // Skip if this mismatch is marked ambiguous
            if ambiguous.contains(&mm_seq) {
                continue;
            }

            // If we've seen this mismatch before with a different parent
            // mark it as ambiguous and remove it from the table
            if let Some(_existing_parent) = mm_table.get(&mm_seq) {
                ambiguous.insert(mm_seq);
                mm_table.remove(&mm_seq);
            } else {
                // First time seeing this mismatch - add it to the table
                mm_table.insert(mm_seq, parent);
            }
        }
    }

    Ok(mm_table)
}

/// Builds a table mapping all possible one-off mismatches to their parent sequences
/// and an ambiguous mismatch table which maps ambiguous sequences to all respective parent sequences.
///
/// For each input sequence, this function:
/// 1. Adds the sequence itself as its own parent
/// 2. Generates all one-off mismatches and maps them to the parent sequence
/// 3. Handles collisions by removing ambiguous mappings
///
/// # Arguments
/// * sequences: Slice of bitnuc scalar sequences
/// * slen: Length of each sequence in bases
///
/// # Returns
/// * Ok(HashMap) mapping mismatches to their unambiguous parent
/// * Err(MismatchError) if sequence length > 32
pub fn build_mismatch_table_with_ambiguous(
    sequences: &[u64],
    slen: usize,
) -> Result<(MismatchTable, AmbiguousMismatchTable), MismatchError> {
    if slen > 32 {
        return Err(MismatchError::InvalidSequenceLength(slen));
    }

    let mut mm_table = HashMap::new();
    let mut ambiguous = HashMap::new();
    let mut mm_buffer = Vec::with_capacity(slen * 3);

    // Initialize the MM table with all the sequences mapping to themselves.
    // This is the case where there are no mismatches.
    for &seq in sequences {
        mm_table.insert(seq, seq);

        // Initialize the ambiguous table with an empty vector for parent sequences
        ambiguous.insert(seq, Vec::new());
    }

    // Process each sequence
    for &parent in sequences {
        // generate mismatches for the parent sequence
        generate_mismatches(parent, slen, &mut mm_buffer)?;

        // process each mismatch
        for &mm_seq in &mm_buffer {
            // Skip if this mismatch is marked ambiguous
            if ambiguous.contains_key(&mm_seq) {
                ambiguous.get_mut(&mm_seq).unwrap().push(parent); // insert parent into ambiguous list
                continue;
            }

            // If we've seen this mismatch before with a different parent
            // mark it as ambiguous and remove it from the table
            if let Some(existing_parent) = mm_table.get(&mm_seq) {
                // Insert both parents into the ambiguous list
                let list = ambiguous.entry(mm_seq).or_default();
                list.push(*existing_parent);
                list.push(parent);

                // Remove the mismatch from the table
                mm_table.remove(&mm_seq);
            } else {
                // First time seeing this mismatch - add it to the table
                mm_table.insert(mm_seq, parent);
            }
        }
    }

    Ok((mm_table, ambiguous))
}

#[cfg(test)]
mod tests {
    use super::*;
    use bitnuc::as_2bit;
    use nucgen::Sequence;
    use rand::thread_rng;

    fn generate_random_bitnuc_sequence(len: usize) -> u64 {
        let mut seq = Sequence::new();
        let mut rng = thread_rng();
        seq.fill_buffer(&mut rng, len);
        bitnuc::as_2bit(seq.bytes()).unwrap()
    }

    #[test]
    fn test_generate_mismatches() {
        // Test sequence: "ACGT"
        let seq = as_2bit(b"ACGT").unwrap();
        let mut mm_buffer = Vec::new();

        generate_mismatches(seq, 4, &mut mm_buffer).unwrap();

        // Should generate 12 mismatches (3 per position)
        assert_eq!(mm_buffer.len(), 12);

        // Verify all expected mismatches are present
        let expected_mismatches = vec![
            b"CCGT", b"GCGT", b"TCGT", b"AAGT", b"AGGT", b"ATGT", b"ACAT", b"ACCT", b"ACTT",
            b"ACGA", b"ACGC", b"ACGG",
        ];

        for expected in expected_mismatches {
            let expected_seq = as_2bit(expected).unwrap();
            assert!(mm_buffer.contains(&expected_seq));
        }

        // Verify that the original sequence is not present
        let original_seq = as_2bit(b"ACGT").unwrap();
        assert!(!mm_buffer.contains(&original_seq));
    }

    #[test]
    fn test_generate_mismatches_sizing() {
        let mut mm_buffer = Vec::new();
        for slen in 1..=32 {
            // clear the buffer
            mm_buffer.clear();

            // Generate a random sequence
            let seq = generate_random_bitnuc_sequence(slen);

            // Generate mismatches
            generate_mismatches(seq, slen, &mut mm_buffer).unwrap();

            // Validate the correct size of the mismatch buffer
            assert_eq!(mm_buffer.len(), slen * 3);

            // Validate that the original sequence is not present
            assert!(mm_buffer.iter().all(|&x| x != seq));
        }
    }

    #[test]
    fn test_build_mismatch_table() {
        // Test with two similar sequences: "ACGT" and "AGGT"
        let seq1 = as_2bit(b"ACGT").unwrap();
        let seq2 = as_2bit(b"AGGT").unwrap();
        let sequences = vec![seq1, seq2];

        let table = build_mismatch_table(&sequences, 4).unwrap();

        // Should contain the original sequences
        assert_eq!(table.get(&seq1), Some(&seq1));
        assert_eq!(table.get(&seq2), Some(&seq2));

        // Test some unambiguous mismatches
        let ccgt = as_2bit(b"CCGT").unwrap();
        let tggt = as_2bit(b"TGGT").unwrap();

        assert_eq!(table.get(&ccgt), Some(&seq1));
        assert_eq!(table.get(&tggt), Some(&seq2));

        // Test an ambiguous mismatch (should not be present)
        let atgt = as_2bit(b"ATGT").unwrap();
        assert!(table.get(&atgt).is_none());
    }

    #[test]
    fn test_sequence_length_validation() {
        let mut mm_buffer = Vec::new();
        assert!(generate_mismatches(0, 33, &mut mm_buffer).is_err());
        assert!(build_mismatch_table(&[0], 33).is_err());
    }

    #[test]
    fn test_edge_cases() {
        // Test single-base sequences
        let seq_a = as_2bit(b"A").unwrap();
        let seq_t = as_2bit(b"T").unwrap();
        let sequences = vec![seq_a, seq_t];

        let table = build_mismatch_table(&sequences, 1).unwrap();

        // Should have 2 entries total (2 originals + 0 unambiguous mutations)
        assert_eq!(table.len(), 2);

        // Original sequences map to themselves
        assert_eq!(table.get(&seq_a), Some(&seq_a));
        assert_eq!(table.get(&seq_t), Some(&seq_t));

        // No mismatches should be present since they are all ambiguous
        let seq_c = as_2bit(b"C").unwrap();
        let seq_g = as_2bit(b"G").unwrap();
        assert!(table.get(&seq_c).is_none());
        assert!(table.get(&seq_g).is_none());
    }

    #[test]
    fn test_edge_cases_with_ambiguous() {
        // Test single-base sequences
        let seq_a = as_2bit(b"A").unwrap();
        let seq_t = as_2bit(b"T").unwrap();
        let sequences = vec![seq_a, seq_t];

        let (table, atable) = build_mismatch_table_with_ambiguous(&sequences, 1).unwrap();

        // Should have 2 entries total (2 originals + 0 unambiguous mutations)
        assert_eq!(table.len(), 2);

        // Original sequences map to themselves
        assert_eq!(table.get(&seq_a), Some(&seq_a));
        assert_eq!(table.get(&seq_t), Some(&seq_t));

        // No mismatches should be present since they are all ambiguous
        let seq_c = as_2bit(b"C").unwrap();
        let seq_g = as_2bit(b"G").unwrap();
        assert!(table.get(&seq_c).is_none());
        assert!(table.get(&seq_g).is_none());

        assert_eq!(atable.get(&seq_c), Some(&sequences));
        assert_eq!(atable.get(&seq_g), Some(&sequences));
    }

    #[test]
    fn test_max_length_sequence() {
        // Test with maximum length sequence (32 bases)
        let max_seq = generate_random_bitnuc_sequence(32);
        let mut mm_buffer = Vec::new();

        // Should succeed for length 32
        assert!(generate_mismatches(max_seq, 32, &mut mm_buffer).is_ok());
        assert_eq!(mm_buffer.len(), 32 * 3); // 3 mutations per position

        // Should succeed for table building too
        let table = build_mismatch_table(&[max_seq], 32).unwrap();
        assert!(table.contains_key(&max_seq));
    }

    #[test]
    fn test_ambiguity_corner_cases() {
        // Test three sequences that form an ambiguity triangle
        let seq1 = as_2bit(b"AAAA").unwrap();
        let seq2 = as_2bit(b"AAAT").unwrap();
        let seq3 = as_2bit(b"AATA").unwrap();

        let sequences = vec![seq1, seq2, seq3];
        let table = build_mismatch_table(&sequences, 4).unwrap();

        // Original sequences should map to themselves
        assert_eq!(table.get(&seq1), Some(&seq1));
        assert_eq!(table.get(&seq2), Some(&seq2));
        assert_eq!(table.get(&seq3), Some(&seq3));

        // AATT should be ambiguous (1 mutation from both seq2 and seq3)
        let ambiguous = as_2bit(b"AATT").unwrap();
        assert!(table.get(&ambiguous).is_none());
    }

    #[test]
    fn test_ambiguity_corner_cases_with_ambiguous() {
        // Test three sequences that form an ambiguity triangle
        let seq1 = as_2bit(b"AAAA").unwrap();
        let seq2 = as_2bit(b"AAAT").unwrap();
        let seq3 = as_2bit(b"AATA").unwrap();

        let sequences = vec![seq1, seq2, seq3];
        let (table, atable) = build_mismatch_table_with_ambiguous(&sequences, 4).unwrap();

        // Original sequences should map to themselves
        assert_eq!(table.get(&seq1), Some(&seq1));
        assert_eq!(table.get(&seq2), Some(&seq2));
        assert_eq!(table.get(&seq3), Some(&seq3));

        // AATT should be ambiguous (1 mutation from both seq2 and seq3)
        let expected_ambiguous = vec![seq2, seq3];
        let ambiguous = as_2bit(b"AATT").unwrap();
        assert!(table.get(&ambiguous).is_none());
        assert!(atable.get(&ambiguous).is_some());
        assert_eq!(atable.get(&ambiguous), Some(&expected_ambiguous));
    }

    #[test]
    fn test_zero_length() {
        // While this isn't a valid use case, we should test how we handle it
        let table = build_mismatch_table(&[], 0).unwrap();
        assert!(table.is_empty());

        let mut buffer = Vec::new();
        assert!(generate_mismatches(0, 0, &mut buffer).is_ok());
        assert!(buffer.is_empty());
    }

    #[test]
    fn test_zero_length_with_ambiguous() {
        let (table, atable) = build_mismatch_table_with_ambiguous(&[], 0).unwrap();
        assert!(table.is_empty());
        assert!(atable.is_empty());
    }

    #[test]
    fn test_parents() {
        let sequences = vec![
            as_2bit(b"AAAA").unwrap(),
            as_2bit(b"AAAT").unwrap(),
            as_2bit(b"AATA").unwrap(),
        ];
        let table = build_mismatch_table(&sequences, 4).unwrap();

        // Original sequences should map to themselves
        assert_eq!(table.get(&sequences[0]), Some(&sequences[0]));
        assert_eq!(table.get(&sequences[1]), Some(&sequences[1]));
        assert_eq!(table.get(&sequences[2]), Some(&sequences[2]));

        // AATT should be ambiguous (1 mutation from both seq2 and seq3)
        let ambiguous = as_2bit(b"AATT").unwrap();
        assert!(table.get(&ambiguous).is_none());
    }

    #[test]
    fn test_parents_with_ambiguous() {
        let sequences = vec![
            as_2bit(b"AAAA").unwrap(),
            as_2bit(b"AAAT").unwrap(),
            as_2bit(b"AATA").unwrap(),
        ];
        let (table, atable) = build_mismatch_table_with_ambiguous(&sequences, 4).unwrap();

        // Original sequences should map to themselves
        assert_eq!(table.get(&sequences[0]), Some(&sequences[0]));
        assert_eq!(table.get(&sequences[1]), Some(&sequences[1]));
        assert_eq!(table.get(&sequences[2]), Some(&sequences[2]));

        // AATT should be ambiguous (1 mutation from both seq2 and seq3)
        let ambiguous = as_2bit(b"AATT").unwrap();
        assert!(table.get(&ambiguous).is_none());

        // AATT should map to expected parents
        let expected_parents = vec![sequences[1], sequences[2]];
        assert_eq!(atable.get(&ambiguous), Some(&expected_parents));
    }
}
