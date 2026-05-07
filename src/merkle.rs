//! Minimal Merkle commitment over polynomial evaluations using `p3-merkle-tree` v0.5.2.
//!
//! Configuration (matches the canonical Plonky3 v0.5.2 test pattern in
//! `merkle-tree/src/mmcs.rs` adapted for Goldilocks):
//! - Permutation: width-8 Poseidon2 over Goldilocks (`Poseidon2Goldilocks<8>`),
//!   constructed via `new_from_rng_128` with a seeded `SmallRng` so commitments
//!   are deterministic for a given (seed, evaluations) pair.
//! - Leaf hash: `PaddingFreeSponge<perm, WIDTH=8, RATE=4, OUT=4>`.
//! - Internal compression: `TruncatedPermutation<perm, N=2, CHUNK=4, WIDTH=8>`.
//! - MMCS: arity 2, 4-element digest, `cap_height = 0` so the commitment
//!   contains the single root.
//!
//! Sources (Plonky3 v0.5.2):
//! - https://github.com/Plonky3/Plonky3/blob/v0.5.2/merkle-tree/src/lib.rs
//! - https://github.com/Plonky3/Plonky3/blob/v0.5.2/merkle-tree/src/mmcs.rs
//! - https://docs.rs/p3-poseidon2/0.5.2/p3_poseidon2/struct.Poseidon2.html
//! - https://docs.rs/p3-goldilocks/0.5.2/p3_goldilocks/type.Poseidon2Goldilocks.html

use p3_commit::Mmcs;
use p3_field::Field;
use p3_goldilocks::{Goldilocks, Poseidon2Goldilocks};
use p3_matrix::dense::RowMajorMatrix;
use p3_merkle_tree::{MerkleCap, MerkleTree, MerkleTreeMmcs};
use p3_symmetric::{PaddingFreeSponge, TruncatedPermutation};
use rand::SeedableRng;
use rand::rngs::SmallRng;

/// Width of the Poseidon2 permutation backing every hash/compression/challenger
/// over Goldilocks in this crate. Shared with `crate::fri` so a single seed
/// produces the same `Perm` instance for the MMCS and the Fiat-Shamir transcript.
pub const WIDTH: usize = 8;
/// Sponge / challenger absorption rate (capacity = `WIDTH - RATE` = 4).
pub const RATE: usize = 4;
const DIGEST_ELEMS: usize = 4;

/// Width-8 Poseidon2 over Goldilocks. Shared between this module and `crate::fri`.
pub type Perm = Poseidon2Goldilocks<WIDTH>;
type Hasher = PaddingFreeSponge<Perm, WIDTH, RATE, DIGEST_ELEMS>;
type Compress = TruncatedPermutation<Perm, 2, DIGEST_ELEMS, WIDTH>;

/// MMCS specialised for Goldilocks evaluations with the configuration documented
/// at the module level. Use [`make_mmcs`] to construct an instance.
pub type GoldilocksMmcs = MerkleTreeMmcs<
    <Goldilocks as Field>::Packing,
    <Goldilocks as Field>::Packing,
    Hasher,
    Compress,
    2,
    DIGEST_ELEMS,
>;

/// Verifier-side commitment (cap-height-0 Merkle root, 4 Goldilocks elements).
pub type Commitment = MerkleCap<Goldilocks, [Goldilocks; DIGEST_ELEMS]>;

/// Prover-side state needed to later open the commitment.
pub type ProverData =
    MerkleTree<Goldilocks, Goldilocks, RowMajorMatrix<Goldilocks>, 2, DIGEST_ELEMS>;

/// Bundle of prover/verifier outputs from a single commitment.
pub struct CommitResult {
    /// Sent to the verifier.
    pub commitment: Commitment,
    /// Retained by the prover for later opening proofs.
    pub prover_data: ProverData,
}

/// Build a deterministic [`GoldilocksMmcs`]. The seed parameterises the
/// Poseidon2 round constants — same seed → same commitment for the same input.
pub fn make_mmcs(seed: u64) -> GoldilocksMmcs {
    let mut rng = SmallRng::seed_from_u64(seed);
    let perm = Perm::new_from_rng_128(&mut rng);
    let hasher = Hasher::new(perm.clone());
    let compress = Compress::new(perm);
    GoldilocksMmcs::new(hasher, compress, 0)
}

/// Commit to a single evaluation vector (treated as one column).
pub fn commit_evaluations(mmcs: &GoldilocksMmcs, evaluations: Vec<Goldilocks>) -> CommitResult {
    let (commitment, prover_data) = mmcs.commit_vec(evaluations);
    CommitResult { commitment, prover_data }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::field_arith::from_u64;
    use crate::poly::DensePoly;

    fn evals(vs: &[u64]) -> Vec<Goldilocks> {
        vs.iter().copied().map(from_u64).collect()
    }

    #[test]
    fn commitment_is_deterministic() {
        let mmcs = make_mmcs(0);
        let v = evals(&[1, 2, 3, 4, 5, 6, 7, 8]);
        let c1 = commit_evaluations(&mmcs, v.clone());
        let c2 = commit_evaluations(&mmcs, v);
        assert_eq!(c1.commitment, c2.commitment);
    }

    #[test]
    fn different_inputs_yield_different_commitments() {
        let mmcs = make_mmcs(0);
        let a = commit_evaluations(&mmcs, evals(&[1, 2, 3, 4, 5, 6, 7, 8]));
        let b = commit_evaluations(&mmcs, evals(&[1, 2, 3, 4, 5, 6, 7, 9]));
        assert_ne!(a.commitment, b.commitment);
    }

    #[test]
    fn different_seeds_yield_different_commitments() {
        let v = evals(&[1, 2, 3, 4, 5, 6, 7, 8]);
        let a = commit_evaluations(&make_mmcs(1), v.clone());
        let b = commit_evaluations(&make_mmcs(2), v);
        assert_ne!(a.commitment, b.commitment);
    }

    #[test]
    fn commits_polynomial_ntt_evaluations() {
        let mmcs = make_mmcs(42);
        let p1 = DensePoly::new(evals(&[7, 11, 13, 17]));
        let p2 = DensePoly::new(evals(&[7, 11, 13, 18]));
        assert_ne!(
            commit_evaluations(&mmcs, p1.ntt()).commitment,
            commit_evaluations(&mmcs, p2.ntt()).commitment,
        );
    }
}
