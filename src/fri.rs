//! FRI low-degree-test scaffolding for Goldilocks evaluations.
//!
//! Wires the existing [`crate::merkle::GoldilocksMmcs`] into the `p3-fri`
//! v0.5.2 polynomial commitment scheme `TwoAdicFriPcs` — the canonical
//! user-facing entry point that internally drives `prove_fri` / `verify_fri`.
//!
//! Module surface:
//! - **Slice 1 — setup**: type aliases pinning the full FRI stack for
//!   Goldilocks, deterministic builders for the FRI parameters, PCS, and
//!   matching Fiat-Shamir challenger, all parameterised by a single `seed`.
//! - **Slice 2 — round-trip**: [`prove_low_degree`] and [`verify_low_degree`]
//!   implement a single-polynomial commit → open → verify path on top of
//!   the slice-1 scaffolding, mirroring the BabyBear template at
//!   `fri/tests/fri.rs` adapted to a single matrix / single column / single
//!   opening point.
//!
//! Sources (Plonky3 v0.5.2):
//! - https://github.com/Plonky3/Plonky3/blob/v0.5.2/fri/src/lib.rs
//! - https://github.com/Plonky3/Plonky3/blob/v0.5.2/fri/src/config.rs
//! - https://github.com/Plonky3/Plonky3/blob/v0.5.2/fri/src/two_adic_pcs.rs
//! - https://github.com/Plonky3/Plonky3/blob/v0.5.2/fri/tests/fri.rs (BabyBear template)
//! - https://github.com/Plonky3/Plonky3/blob/v0.5.2/commit/src/pcs.rs (Pcs trait)
//! - https://github.com/Plonky3/Plonky3/blob/v0.5.2/goldilocks/src/extension.rs

use p3_challenger::{CanObserve, DuplexChallenger, FieldChallenger};
use p3_commit::{ExtensionMmcs, Pcs as PcsTrait};
use p3_dft::Radix2Dit;
use p3_field::extension::BinomialExtensionField;
use p3_fri::{FriParameters, TwoAdicFriPcs, create_test_fri_params};
use p3_goldilocks::Goldilocks;
use p3_matrix::dense::RowMajorMatrix;
use rand::SeedableRng;
use rand::rngs::SmallRng;

use crate::merkle::{GoldilocksMmcs, Perm, RATE, WIDTH, make_mmcs};
use crate::poly::DensePoly;

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

/// PCS-side associated types, surfaced as crate-level aliases for ergonomic
/// signatures in [`LowDegreeProof`] / [`prove_low_degree`] / [`verify_low_degree`].
pub type Domain = <Pcs as PcsTrait<Challenge, Challenger>>::Domain;
pub type Commitment = <Pcs as PcsTrait<Challenge, Challenger>>::Commitment;
pub type Proof = <Pcs as PcsTrait<Challenge, Challenger>>::Proof;
pub type VerifyError = <Pcs as PcsTrait<Challenge, Challenger>>::Error;

/// Output of [`prove_low_degree`]. Bundles everything the verifier needs
/// (commitment, claimed evaluation, FRI proof, trace size) plus the
/// prover-side `zeta` for diagnostics.
pub struct LowDegreeProof {
    /// Sent to the verifier. The Merkle cap over the LDE evaluations.
    pub commitment: Commitment,
    /// Claimed value `f(zeta)` at the FRI opening point, sent to the verifier.
    pub opened_value: Challenge,
    /// FRI low-degree-test proof produced by `TwoAdicFriPcs::open`.
    pub proof: Proof,
    /// `log2` of the trace domain size — public, both sides need it to
    /// reconstruct the natural domain.
    pub log_height: usize,
    /// Opening point sampled by the prover. The verifier re-derives this from
    /// its own challenger and ignores the field; included for diagnostics.
    pub zeta: Challenge,
}

/// Commit to a [`DensePoly`] via its NTT evaluations on the natural 2-adic
/// subgroup, observe the commitment in the transcript, sample an opening
/// point ζ via Fiat-Shamir, open at ζ, and produce a FRI proof.
///
/// The matrix is single-column (one polynomial), the batch is single-element
/// (one commitment), and there is one opening point (`zeta`). For multi-poly
/// batching, see the upstream BabyBear template at `fri/tests/fri.rs`.
pub fn prove_low_degree(
    pcs: &Pcs,
    poly: &DensePoly,
    challenger: &mut Challenger,
) -> LowDegreeProof {
    let evals = poly.ntt();
    assert!(
        evals.len() >= 2,
        "polynomial must yield at least 2 NTT evaluations (got {})",
        evals.len()
    );
    let height = evals.len();
    let log_height = height.trailing_zeros() as usize;

    let domain =
        <Pcs as PcsTrait<Challenge, Challenger>>::natural_domain_for_degree(pcs, height);
    let mat = RowMajorMatrix::new(evals, 1);

    let (commitment, prover_data) =
        <Pcs as PcsTrait<Challenge, Challenger>>::commit(pcs, vec![(domain, mat)]);

    challenger.observe(commitment.clone());
    let zeta: Challenge = challenger.sample_algebra_element();

    let open_data = vec![(&prover_data, vec![vec![zeta]])];
    let (opened_values, proof) = pcs.open(open_data, challenger);
    let opened_value = opened_values[0][0][0][0];

    LowDegreeProof { commitment, opened_value, proof, log_height, zeta }
}

/// Verify a [`LowDegreeProof`]. The verifier re-derives `zeta` from its own
/// challenger after observing the commitment; `proof.zeta` is not consulted.
pub fn verify_low_degree(
    pcs: &Pcs,
    proof: &LowDegreeProof,
    challenger: &mut Challenger,
) -> Result<(), VerifyError> {
    challenger.observe(proof.commitment.clone());
    let zeta: Challenge = challenger.sample_algebra_element();

    let domain = <Pcs as PcsTrait<Challenge, Challenger>>::natural_domain_for_degree(
        pcs,
        1 << proof.log_height,
    );

    let commitments_with_opening_points = vec![(
        proof.commitment.clone(),
        vec![(domain, vec![(zeta, vec![proof.opened_value])])],
    )];

    pcs.verify(commitments_with_opening_points, &proof.proof, challenger)
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

    use crate::field_arith::from_u64;
    use p3_field::PrimeCharacteristicRing;

    fn small_poly() -> DensePoly {
        // 8 coefficients → log_height = 3 trace, log_lde = 5 with log_blowup = 2.
        DensePoly::new((1u64..=8).map(from_u64).collect())
    }

    #[test]
    fn round_trip_proves_and_verifies() {
        let seed = 7;
        let pcs = make_pcs(seed, 0);
        let mut p_chal = make_challenger(seed);
        let proof = prove_low_degree(&pcs, &small_poly(), &mut p_chal);

        let mut v_chal = make_challenger(seed);
        verify_low_degree(&pcs, &proof, &mut v_chal)
            .expect("honest proof should verify");
    }

    #[test]
    fn prover_and_verifier_derive_the_same_zeta() {
        let seed = 13;
        let pcs = make_pcs(seed, 0);
        let mut p_chal = make_challenger(seed);
        let proof = prove_low_degree(&pcs, &small_poly(), &mut p_chal);

        // Replay the verifier's transcript prefix to confirm zeta agreement.
        let mut v_chal = make_challenger(seed);
        v_chal.observe(proof.commitment.clone());
        let zeta_v: Challenge = v_chal.sample_algebra_element();
        assert_eq!(proof.zeta, zeta_v, "Fiat-Shamir transcripts must agree on zeta");
    }

    #[test]
    fn tampered_opened_value_is_rejected() {
        let seed = 7;
        let pcs = make_pcs(seed, 0);
        let mut p_chal = make_challenger(seed);
        let mut proof = prove_low_degree(&pcs, &small_poly(), &mut p_chal);

        proof.opened_value += Challenge::ONE;

        let mut v_chal = make_challenger(seed);
        verify_low_degree(&pcs, &proof, &mut v_chal)
            .expect_err("tampering with opened_value must fail verification");
    }

    #[test]
    fn tampered_proof_final_poly_is_rejected() {
        let seed = 7;
        let pcs = make_pcs(seed, 0);
        let mut p_chal = make_challenger(seed);
        let mut proof = prove_low_degree(&pcs, &small_poly(), &mut p_chal);

        proof.proof.final_poly[0] += Challenge::ONE;

        let mut v_chal = make_challenger(seed);
        verify_low_degree(&pcs, &proof, &mut v_chal)
            .expect_err("tampering with proof.final_poly must fail verification");
    }
}
