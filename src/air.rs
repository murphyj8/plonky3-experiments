//! Fibonacci AIR — the canonical "hello world" for Plonky3 STARKs.
//!
//! Two main columns `(left, right)` implementing the recurrence
//! ```text
//!     left'  = right
//!     right' = left + right
//! ```
//! plus three public inputs `[a, b, x]` constraining the boundary:
//! ```text
//!     row 0:    left == a, right == b
//!     row N-1:  right == x
//! ```
//!
//! The AIR struct is non-parameterised (per the slice's "do not generalise
//! over fields" constraint), but the trait impls remain generic over the
//! field/builder so `p3-uni-stark::prove` can instantiate them with both
//! `SymbolicAirBuilder<Goldilocks>` (constraint counting) and
//! `ProverConstraintFolder<'_, StarkConfig>` (real evaluation) without
//! requiring duplicate code.
//!
//! Sources (Plonky3 v0.5.2):
//! - https://github.com/Plonky3/Plonky3/blob/v0.5.2/air/src/air.rs
//! - https://github.com/Plonky3/Plonky3/blob/v0.5.2/uni-stark/tests/fib_air.rs

use p3_air::{Air, AirBuilder, BaseAir, WindowAccess};
use p3_goldilocks::Goldilocks;
use p3_matrix::dense::RowMajorMatrix;

use crate::field_arith::from_u64;

/// Number of main-trace columns used by [`FibonacciAir`].
pub const NUM_FIBONACCI_COLS: usize = 2;

/// Two-column Fibonacci AIR with three public inputs `[a, b, x]`.
pub struct FibonacciAir;

impl<F> BaseAir<F> for FibonacciAir {
    fn width(&self) -> usize {
        NUM_FIBONACCI_COLS
    }

    fn num_public_values(&self) -> usize {
        3
    }

    fn max_constraint_degree(&self) -> Option<usize> {
        // is_first_row / is_transition / is_last_row are degree-1 selectors
        // applied to degree-1 (Var - PublicVar) expressions, giving degree 2.
        Some(2)
    }
}

impl<AB: AirBuilder> Air<AB> for FibonacciAir {
    fn eval(&self, builder: &mut AB) {
        let main = builder.main();
        let pis = builder.public_values();
        let a = pis[0];
        let b = pis[1];
        let x = pis[2];

        let local = main.current_slice();
        let next = main.next_slice();
        let local_left = local[0];
        let local_right = local[1];
        let next_left = next[0];
        let next_right = next[1];

        let mut when_first_row = builder.when_first_row();
        when_first_row.assert_eq(local_left, a);
        when_first_row.assert_eq(local_right, b);

        let mut when_transition = builder.when_transition();
        when_transition.assert_eq(local_right, next_left);
        when_transition.assert_eq(local_left + local_right, next_right);

        builder.when_last_row().assert_eq(local_right, x);
    }
}

/// Build the first `n` rows of the Fibonacci sequence starting from `(a, b)`
/// as a 2-column trace matrix over Goldilocks. `n` must be a power of two.
pub fn generate_fibonacci_trace(a: u64, b: u64, n: usize) -> RowMajorMatrix<Goldilocks> {
    assert!(n.is_power_of_two(), "trace length must be a power of two (got {n})");
    let mut values = Vec::with_capacity(n * NUM_FIBONACCI_COLS);
    let (mut left, mut right) = (from_u64(a), from_u64(b));
    for _ in 0..n {
        values.push(left);
        values.push(right);
        let new_right = left + right;
        left = right;
        right = new_right;
    }
    RowMajorMatrix::new(values, NUM_FIBONACCI_COLS)
}

#[cfg(test)]
mod tests {
    use super::*;
    use p3_field::PrimeCharacteristicRing;

    #[test]
    fn trace_starts_at_initial_pair() {
        let trace = generate_fibonacci_trace(0, 1, 8);
        assert_eq!(trace.values[0], Goldilocks::ZERO);
        assert_eq!(trace.values[1], Goldilocks::ONE);
    }

    #[test]
    fn last_row_right_is_eighth_fibonacci() {
        // F_8 = 21 with F_0 = 0, F_1 = 1.
        let trace = generate_fibonacci_trace(0, 1, 8);
        let last_right = trace.values[trace.values.len() - 1];
        assert_eq!(last_right, from_u64(21));
    }

    #[test]
    #[should_panic(expected = "power of two")]
    fn non_power_of_two_length_rejected() {
        let _ = generate_fibonacci_trace(0, 1, 5);
    }
}
