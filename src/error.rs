/// Error type for mismatch generation
///
/// This error is returned when the input sequence length is greater than 32.
#[derive(thiserror::Error, Debug)]
pub enum MismatchError {
    #[error("Invalid sequence length: {0}")]
    InvalidSequenceLength(usize),
}
