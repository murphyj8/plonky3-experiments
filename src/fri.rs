//! FRI low-degree-test scaffolding for Goldilocks evaluations.
//!
//! This module wires the existing [`crate::merkle::GoldilocksMmcs`] into the
//! `p3-fri` v0.5.2 polynomial commitment scheme `TwoAdicFriPcs` — the canonical
//! user-facing entry point that internally drives `prove_fri` / `verify_fri`.
//!
//! Increment scope (slice 1):
//! - Type aliases pinning the full FRI stack for Goldilocks.
//! - Deterministic builders for the FRI parameters, the PCS, and the
//!   matching Fiat-Shamir challenger, all parameterised by a single `seed`.
//! - Setup-invariant tests confirming the parameter math, deterministic
//!   challenger sampling, and that the PCS itself is constructible.
//!
//! The end-to-end commit / open / verify round-trip lands in slice 2 — see
//! `agent_prompts/` for the next prompt.
//!
//! Sources (Plonky3 v0.5.2):
//! - https://github.com/Plonky3/Plonky3/blob/v0.5.2/fri/src/lib.rs
//! - https://github.com/Plonky3/Plonky3/blob/v0.5.2/fri/src/config.rs
//! - https://github.com/Plonky3/Plonky3/blob/v0.5.2/fri/src/two_adic_pcs.rs
//! - https://github.com/Plonky3/Plonky3/blob/v0.5.2/fri/tests/fri.rs (BabyBear template)
//! - https://github.com/Plonky3/Plonky3/blob/v0.5.2/goldilocks/src/extension.rs

use p3_challenger::DuplexChallenger;
use p3_commit::ExtensionMmcs;
use p3_dft::Radix2Dit;
use p3_field::extension::BinomialExtensionField;
use p3_fri::{FriParameters, TwoAdicFriPcs, create_test_fri_params};
use p3_goldilocks::Goldilocks;
use rand::SeedableRng;
use rand::rngs::SmallRng;

use crate::merkle::{GoldilocksMmcs, Perm, RATE, WIDTH, make_mmcs};

/// Base field for trace polynomials.
pub type Val = Goldilocks;

/// Soundness-amplification field for FRI challenges. Goldilocks has a degree-2
/// binomial extension (W = 7, see `goldilocks/src/extension.rs`); 2 × 64 ≈ 128
/// bits of soundness is the standard target for STARK provers over Goldilocks.
pub type Challenge = BinomialExtensionField<Val, 2>;

/// MMCS for FRI commit-phase commitments — same Merkle tree shape as the
/// input MMCS, but lifted over the extension field via [`ExtensionMmcs`].
pub type ChallengeMmcs = ExtensionMmcs<Val, Challenge, GoldilocksMmcs>;

/// Fiat-Shamir transcript driver, sponging over the same width-8 Poseidon2
/// used for the Merkle hashes so the prover and verifier deterministically
/// agree on every sampled challenge.
pub type Challenger = DuplexChallenger<Val, Perm, WIDTH, RATE>;

/// Concrete FRI-based polynomial commitment scheme over Goldilocks.
pub type Pcs = TwoAdicFriPcs<Val, Radix2Dit<Val>, GoldilocksMmcs, ChallengeMmcs>;

/// Build FRI parameters using the upstream `create_test_fri_params` defaults
/// (log_blowup = 2, num_queries = 2, max_log_arity = 1, PoW bits = 1).
/// The seed deterministically picks the Poseidon2 round constants used by the
/// commit-phase MMCS.
pub fn make_fri_params(seed: u64, log_final_poly_len: usize) -> FriParameters<ChallengeMmcs> {
    let challenge_mmcs = ChallengeMmcs::new(make_mmcs(seed));
    create_test_fri_params(challenge_mmcs, log_final_poly_len)
}

/// Build a fully-configured FRI-PCS. The same `seed` parameterises both the
/// input MMCS and the FRI commit-phase MMCS, so a single seed pins the entire
/// hash configuration of the protocol.
pub fn make_pcs(seed: u64, log_final_poly_len: usize) -> Pcs {
    let input_mmcs = make_mmcs(seed);
    let fri_params = make_fri_params(seed, log_final_poly_len);
    Pcs::new(Radix2Dit::default(), input_mmcs, fri_params)
}

/// Build a fresh Fiat-Shamir challenger seeded from the same Poseidon2
/// permutation as [`make_mmcs`] / [`make_pcs`] when called with the same
/// `seed`.
pub fn make_challenger(seed: u64) -> Challenger {
    let mut rng = SmallRng::seed_from_u64(seed);
    let perm = Perm::new_from_rng_128(&mut rng);
    Challenger::new(perm)
}

#[cfg(test)]
mod tests {
    use super::*;
    use p3_challenger::CanSampleBits;

    #[test]
    fn fri_params_have_expected_blowup() {
        let params = make_fri_params(0, 0);
        assert_eq!(params.log_blowup, 2);
        assert_eq!(params.blowup(), 4);
    }

    #[test]
    fn fri_params_have_expected_query_count() {
        let params = make_fri_params(0, 0);
        assert_eq!(params.num_queries, 2);
    }

    #[test]
    fn conjectured_soundness_bits_match_formula() {
        let params = make_fri_params(0, 0);
        // log_blowup * num_queries + query_proof_of_work_bits = 2*2 + 1 = 5
        assert_eq!(params.conjectured_soundness_bits(), 5);
    }

    #[test]
    fn final_poly_len_tracks_parameter() {
        for log_final in 0..4 {
            let params = make_fri_params(0, log_final);
            assert_eq!(params.final_poly_len(), 1 << log_final);
            assert_eq!(params.log_final_poly_len, log_final);
        }
    }

    #[test]
    fn pcs_is_constructible_for_multiple_seeds() {
        for seed in [0u64, 1, 7, 42, 1 << 32] {
            let _pcs = make_pcs(seed, 0);
        }
    }

    #[test]
    fn challenger_is_deterministic_per_seed() {
        let mut a = make_challenger(11);
        let mut b = make_challenger(11);
        for _ in 0..4 {
            assert_eq!(a.sample_bits(16), b.sample_bits(16));
        }
    }

    #[test]
    fn challenger_diverges_on_different_seeds() {
        let mut a = make_challenger(1);
        let mut b = make_challenger(2);
        // Sample several rounds — the probability of all-equal across 4 16-bit
        // samples by chance is 2^-64, well below any reasonable test flake rate.
        let collisions = (0..4).filter(|_| a.sample_bits(16) == b.sample_bits(16)).count();
        assert!(collisions < 4, "challengers from different seeds should not all agree");
    }
}
