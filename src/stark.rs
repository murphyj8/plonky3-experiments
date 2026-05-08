//! Concrete `StarkConfig` and prove/verify glue for Goldilocks AIRs.
//!
//! Reuses the FRI-PCS scaffolding from `crate::fri` (and therefore the
//! `crate::merkle` MMCS plumbing underneath) so a single seed pins every
//! cryptographic primitive in the proof: input MMCS, FRI commit-phase MMCS,
//! and Fiat-Shamir challenger.
//!
//! `prove_stark` / `verify_stark` are direct re-exports of the upstream
//! functions; the wrapper API surface here is the type aliases plus
//! [`make_stark_config`].
//!
//! Sources (Plonky3 v0.5.2):
//! - https://github.com/Plonky3/Plonky3/blob/v0.5.2/uni-stark/src/lib.rs
//! - https://github.com/Plonky3/Plonky3/blob/v0.5.2/uni-stark/src/config.rs
//! - https://github.com/Plonky3/Plonky3/blob/v0.5.2/uni-stark/src/prover.rs
//! - https://github.com/Plonky3/Plonky3/blob/v0.5.2/uni-stark/src/verifier.rs
//! - https://github.com/Plonky3/Plonky3/blob/v0.5.2/uni-stark/tests/fib_air.rs (template)

use p3_uni_stark::{Proof as InnerProof, VerificationError};

use crate::fri::{Challenge, Challenger, Pcs, VerifyError, make_challenger, make_pcs};

pub use p3_uni_stark::{prove as prove_stark, verify as verify_stark};

/// Concrete `StarkConfig` over Goldilocks + FRI.
pub type StarkConfig = p3_uni_stark::StarkConfig<Pcs, Challenge, Challenger>;

/// STARK proof produced by [`prove_stark`].
pub type StarkProof = InnerProof<StarkConfig>;

/// Error returned by [`verify_stark`] on rejection.
pub type StarkError = VerificationError<VerifyError>;

/// Build a deterministic [`StarkConfig`] keyed off `seed`. The seed pins the
/// Poseidon2 round constants for the input MMCS, the FRI commit-phase MMCS,
/// and the Fiat-Shamir challenger — same seed → same proof for the same trace
/// and public inputs.
pub fn make_stark_config(seed: u64, log_final_poly_len: usize) -> StarkConfig {
    let pcs = make_pcs(seed, log_final_poly_len);
    let challenger = make_challenger(seed);
    StarkConfig::new(pcs, challenger)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::air::{FibonacciAir, SumAir, generate_fibonacci_trace, generate_sum_trace};
    use crate::field_arith::from_u64;
    use p3_field::PrimeCharacteristicRing;
    use p3_goldilocks::Goldilocks;
    use p3_matrix::dense::RowMajorMatrix;

    fn fib_setup(n: usize, x: u64) -> (StarkConfig, Vec<Goldilocks>, RowMajorMatrix<Goldilocks>) {
        let trace = generate_fibonacci_trace(0, 1, n);
        let config = make_stark_config(1, 0);
        let pis = vec![Goldilocks::ZERO, Goldilocks::ONE, from_u64(x)];
        (config, pis, trace)
    }

    #[test]
    fn fibonacci_round_trip() {
        // Fibonacci over Goldilocks starting at (0, 1): row 7 has right = F_8 = 21.
        let (config, pis, trace) = fib_setup(8, 21);
        let proof = prove_stark(&config, &FibonacciAir, trace, &pis);
        verify_stark(&config, &FibonacciAir, &proof, &pis)
            .expect("honest STARK proof should verify");
    }

    /// Tampering with the public inputs trips the prover's own debug-build
    /// constraint check, mirroring the BabyBear `test_incorrect_public_value`.
    /// In release mode the prover would succeed and the verifier would reject;
    /// the panic path is the canonical Plonky3 negative test.
    #[cfg(debug_assertions)]
    #[test]
    #[should_panic(expected = "constraints not satisfied on row")]
    fn wrong_public_input_is_rejected() {
        let (config, pis, trace) = fib_setup(8, 999);
        let _ = prove_stark(&config, &FibonacciAir, trace, &pis);
    }

    fn sum_setup(
        n: usize,
        init: u64,
        claimed_final: u64,
    ) -> (StarkConfig, Vec<Goldilocks>, RowMajorMatrix<Goldilocks>) {
        let trace = generate_sum_trace(init, n);
        let config = make_stark_config(1, 0);
        let pis = vec![from_u64(init), from_u64(claimed_final)];
        (config, pis, trace)
    }

    #[test]
    fn running_sum_round_trip() {
        // init = 0, n = 8 → row 7 holds acc = 7. Honest claim: acc_final = 7.
        let (config, pis, trace) = sum_setup(8, 0, 7);
        let proof = prove_stark(&config, &SumAir, trace, &pis);
        verify_stark(&config, &SumAir, &proof, &pis)
            .expect("honest running-sum STARK should verify");
    }

    #[cfg(debug_assertions)]
    #[test]
    #[should_panic(expected = "constraints not satisfied on row")]
    fn running_sum_wrong_final_is_rejected() {
        // Honest trace ends at acc=7, but we claim acc_final=999. The prover's
        // debug-build constraint check fails on the last row.
        let (config, pis, trace) = sum_setup(8, 0, 999);
        let _ = prove_stark(&config, &SumAir, trace, &pis);
    }

    /// Release-mode negative test: produce an honest proof, then mutate one
    /// element of `opened_values.trace_local`. The verifier observes the
    /// tampered claim, but the FRI opening proof was generated against the
    /// original value and the commitment binds it — `pcs.verify` must reject
    /// before the constraint check ever runs, yielding
    /// `InvalidOpeningArgument(_)`. The `OodEvaluationMismatch` arm is
    /// matched too so the test stays correct if upstream ever reorders the
    /// verifier's checks.
    ///
    /// Source: https://github.com/Plonky3/Plonky3/blob/v0.5.2/uni-stark/src/verifier.rs
    #[test]
    fn tampered_trace_local_is_rejected() {
        use crate::fri::Challenge;
        let (config, pis, trace) = fib_setup(8, 21);
        let mut proof = prove_stark(&config, &FibonacciAir, trace, &pis);
        proof.opened_values.trace_local[0] += Challenge::ONE;

        let result = verify_stark(&config, &FibonacciAir, &proof, &pis);
        assert!(
            matches!(
                result,
                Err(StarkError::InvalidOpeningArgument(_))
                    | Err(StarkError::OodEvaluationMismatch { .. })
            ),
            "expected InvalidOpeningArgument or OodEvaluationMismatch, got {result:?}"
        );
    }
}
