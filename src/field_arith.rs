//! Basic field arithmetic over the Goldilocks prime field
//! `p = 2^64 - 2^32 + 1`, backed by Plonky3 v0.5.2.
//!
//! Sources (Plonky3 v0.5.2):
//! - https://github.com/Plonky3/Plonky3/blob/v0.5.2/goldilocks/src/goldilocks.rs
//! - https://docs.rs/p3-field/0.5.2/p3_field/trait.PrimeCharacteristicRing.html
//! - https://docs.rs/p3-field/0.5.2/p3_field/trait.Field.html

use p3_field::{Field, PrimeField64, TwoAdicField};
use p3_goldilocks::Goldilocks;

/// The Goldilocks prime, `p = 2^64 - 2^32 + 1` = `0xFFFF_FFFF_0000_0001`.
/// Source: `impl PrimeField64 for Goldilocks` in goldilocks.rs.
pub const GOLDILOCKS_PRIME: u64 = <Goldilocks as PrimeField64>::ORDER_U64;

/// Construct a field element from any `u64` (no reduction needed:
/// Goldilocks uses a non-canonical internal representation).
#[inline]
pub fn from_u64(v: u64) -> Goldilocks {
    Goldilocks::new(v)
}

/// Canonical `u64` representative in `[0, p)`.
#[inline]
pub fn to_canonical_u64(x: Goldilocks) -> u64 {
    x.as_canonical_u64()
}

/// Multiplicative inverse, or `None` if `x == 0`. Source: `Field::try_inverse`.
#[inline]
pub fn try_inverse(x: Goldilocks) -> Option<Goldilocks> {
    x.try_inverse()
}

/// Generator of the order-`2^bits` subgroup of `F_p^*` (`TWO_ADICITY = 32`).
#[inline]
pub fn two_adic_generator(bits: usize) -> Goldilocks {
    Goldilocks::two_adic_generator(bits)
}

#[cfg(test)]
mod tests {
    use super::*;
    use p3_field::PrimeCharacteristicRing;

    #[test]
    fn prime_modulus_value() {
        assert_eq!(GOLDILOCKS_PRIME as u128, (1u128 << 64) - (1u128 << 32) + 1);
        assert_eq!(GOLDILOCKS_PRIME, 0xFFFF_FFFF_0000_0001);
    }

    #[test]
    fn add_wraps_around_prime() {
        let a = from_u64(GOLDILOCKS_PRIME - 1);
        let b = from_u64(2);
        assert_eq!(to_canonical_u64(a + b), 1);
    }

    #[test]
    fn neg_and_sub_agree() {
        let a = from_u64(123_456_789);
        let b = from_u64(987_654_321);
        assert_eq!(a - b, a + (-b));
        assert_eq!(to_canonical_u64(a - a), 0);
    }

    #[test]
    fn mul_and_inverse_round_trip() {
        let a = from_u64(0xDEAD_BEEF_DEAD_BEEF);
        let inv = try_inverse(a).expect("nonzero elements are invertible");
        assert_eq!(to_canonical_u64(a * inv), 1);
        assert!(try_inverse(Goldilocks::ZERO).is_none());
    }

    #[test]
    fn two_adic_generator_has_correct_order() {
        let bits = 5usize;
        let g = two_adic_generator(bits);
        assert_eq!(g.exp_power_of_2(bits), Goldilocks::ONE);
        assert_eq!(g.exp_power_of_2(bits - 1), Goldilocks::NEG_ONE);
    }

    #[test]
    fn sum_iter_matches_fold() {
        let xs: Vec<Goldilocks> = (1u64..=10).map(from_u64).collect();
        let folded = xs.iter().copied().fold(Goldilocks::ZERO, |a, b| a + b);
        let summed: Goldilocks = xs.iter().copied().sum();
        assert_eq!(summed, folded);
        assert_eq!(to_canonical_u64(summed), 55);
    }
}
