//! Dense univariate polynomials over Goldilocks.
//!
//! Coefficient ordering: `coeffs[i]` = coefficient of `x^i` (low-to-high),
//! matching the Plonky3 v0.5.2 `p3-dft` convention `y_i = sum_j c_j (s g^i)^j`
//! so a future `RowMajorMatrix::new_col(coeffs)` handoff is a no-op re-index.
//!
//! Sources (Plonky3 v0.5.2):
//! - https://github.com/Plonky3/Plonky3/blob/v0.5.2/dft/src/traits.rs
//! - https://github.com/Plonky3/Plonky3/blob/v0.5.2/goldilocks/src/goldilocks.rs

use core::ops::{Add, Mul, Neg, Sub};

use p3_field::{Field, PrimeCharacteristicRing};
use p3_goldilocks::Goldilocks;

/// Dense univariate polynomial over Goldilocks; canonical (no trailing zeros).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DensePoly {
    coeffs: Vec<Goldilocks>,
}

impl DensePoly {
    /// Build from coefficients (low-to-high). Trailing zeros are stripped so
    /// that `degree` and `Eq` operate on a canonical representation.
    pub fn new(mut coeffs: Vec<Goldilocks>) -> Self {
        while coeffs.last().is_some_and(|c| *c == Goldilocks::ZERO) {
            coeffs.pop();
        }
        Self { coeffs }
    }

    /// The zero polynomial.
    pub fn zero() -> Self {
        Self { coeffs: Vec::new() }
    }

    /// True iff this is the zero polynomial.
    pub fn is_zero(&self) -> bool {
        self.coeffs.is_empty()
    }

    /// Degree, or `None` for the zero polynomial (mathematically `-∞`).
    pub fn degree(&self) -> Option<usize> {
        self.coeffs.len().checked_sub(1)
    }

    /// Evaluate via Horner's rule.
    pub fn evaluate(&self, x: Goldilocks) -> Goldilocks {
        self.coeffs
            .iter()
            .rev()
            .copied()
            .fold(Goldilocks::ZERO, |acc, c| acc * x + c)
    }

    /// Polynomial long division: returns `(quotient, remainder)` such that
    /// `self = quotient * divisor + remainder` and `degree(remainder) <
    /// degree(divisor)`. Panics if `divisor` is the zero polynomial.
    pub fn divmod(&self, divisor: &DensePoly) -> (DensePoly, DensePoly) {
        assert!(!divisor.is_zero(), "division by zero polynomial");
        if self.coeffs.len() < divisor.coeffs.len() {
            return (DensePoly::zero(), self.clone());
        }
        let div_deg = divisor.coeffs.len() - 1;
        let div_lead_inv = divisor.coeffs[div_deg].inverse();
        let mut rem = self.coeffs.clone();
        let mut quot = vec![Goldilocks::ZERO; rem.len() - div_deg];
        while rem.len() > div_deg {
            let factor = *rem.last().unwrap() * div_lead_inv;
            let q_idx = rem.len() - 1 - div_deg;
            quot[q_idx] = factor;
            for (i, c) in divisor.coeffs.iter().enumerate() {
                rem[q_idx + i] -= factor * *c;
            }
            while rem.last().is_some_and(|c| *c == Goldilocks::ZERO) {
                rem.pop();
            }
        }
        (DensePoly::new(quot), DensePoly::new(rem))
    }
}

fn pad_get(xs: &[Goldilocks], i: usize) -> Goldilocks {
    xs.get(i).copied().unwrap_or(Goldilocks::ZERO)
}

impl Add for &DensePoly {
    type Output = DensePoly;
    fn add(self, rhs: Self) -> DensePoly {
        let n = self.coeffs.len().max(rhs.coeffs.len());
        DensePoly::new((0..n).map(|i| pad_get(&self.coeffs, i) + pad_get(&rhs.coeffs, i)).collect())
    }
}

impl Sub for &DensePoly {
    type Output = DensePoly;
    fn sub(self, rhs: Self) -> DensePoly {
        let n = self.coeffs.len().max(rhs.coeffs.len());
        DensePoly::new((0..n).map(|i| pad_get(&self.coeffs, i) - pad_get(&rhs.coeffs, i)).collect())
    }
}

impl Neg for &DensePoly {
    type Output = DensePoly;
    fn neg(self) -> DensePoly {
        // Negation preserves canonical form (nonzero stays nonzero), so we can
        // skip re-canonicalisation through `new`.
        DensePoly { coeffs: self.coeffs.iter().map(|c| -*c).collect() }
    }
}

impl Mul for &DensePoly {
    type Output = DensePoly;
    fn mul(self, rhs: Self) -> DensePoly {
        if self.is_zero() || rhs.is_zero() {
            return DensePoly::zero();
        }
        let mut out = vec![Goldilocks::ZERO; self.coeffs.len() + rhs.coeffs.len() - 1];
        for (i, a) in self.coeffs.iter().enumerate() {
            for (j, b) in rhs.coeffs.iter().enumerate() {
                out[i + j] += *a * *b;
            }
        }
        DensePoly::new(out)
    }
}

impl Mul<Goldilocks> for &DensePoly {
    type Output = DensePoly;
    fn mul(self, rhs: Goldilocks) -> DensePoly {
        DensePoly::new(self.coeffs.iter().map(|c| *c * rhs).collect())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::field_arith::from_u64;

    fn p(cs: &[u64]) -> DensePoly {
        DensePoly::new(cs.iter().copied().map(from_u64).collect())
    }

    #[test]
    fn zero_polynomial() {
        let z = DensePoly::zero();
        assert!(z.is_zero());
        assert_eq!(z.degree(), None);
        assert_eq!(z.evaluate(from_u64(42)), Goldilocks::ZERO);
    }

    #[test]
    fn new_strips_trailing_zeros() {
        let q = p(&[3, 0, 0]);
        assert_eq!(q.degree(), Some(0));
        assert_eq!(q, p(&[3]));
    }

    #[test]
    fn evaluate_via_horner() {
        // f(x) = 1 + 2x + 3x^2 ; f(5) = 1 + 10 + 75 = 86
        let f = p(&[1, 2, 3]);
        assert_eq!(f.evaluate(from_u64(5)), from_u64(86));
        assert_eq!(f.degree(), Some(2));
    }

    #[test]
    fn addition_pads_and_cancels() {
        // (1 + 2x + 3x^2) + (4 + 5x) = 5 + 7x + 3x^2
        assert_eq!(&p(&[1, 2, 3]) + &p(&[4, 5]), p(&[5, 7, 3]));
        // Cancellation strips leading zeros: x^2 + (-x^2) = 0
        let x2 = p(&[0, 0, 1]);
        let neg_x2 = DensePoly::new(vec![Goldilocks::ZERO, Goldilocks::ZERO, -Goldilocks::ONE]);
        assert!((&x2 + &neg_x2).is_zero());
    }

    #[test]
    fn neg_inverts_addition() {
        let f = p(&[1, 2, 3]);
        assert!((&f + &(-&f)).is_zero());
    }

    #[test]
    fn sub_matches_add_negation() {
        let f = p(&[1, 2, 3]);
        let g = p(&[4, 5]);
        assert_eq!(&f - &g, &f + &(-&g));
        assert!((&f - &f).is_zero());
    }

    #[test]
    fn mul_polynomial_schoolbook() {
        // (1 + x)(1 - x) = 1 - x^2
        let one_plus_x = p(&[1, 1]);
        let one_minus_x = DensePoly::new(vec![Goldilocks::ONE, -Goldilocks::ONE]);
        let expected = DensePoly::new(vec![Goldilocks::ONE, Goldilocks::ZERO, -Goldilocks::ONE]);
        assert_eq!(&one_plus_x * &one_minus_x, expected);
    }

    #[test]
    fn mul_distributes_over_evaluation() {
        // (f * g)(x) = f(x) * g(x) for any x
        let f = p(&[1, 2, 3]);
        let g = p(&[4, 5]);
        let x = from_u64(7);
        assert_eq!((&f * &g).evaluate(x), f.evaluate(x) * g.evaluate(x));
    }

    #[test]
    fn mul_with_zero_is_zero() {
        assert!((&p(&[1, 2, 3]) * &DensePoly::zero()).is_zero());
        assert!((&DensePoly::zero() * &p(&[1, 2, 3])).is_zero());
    }

    #[test]
    fn mul_scalar_scales_each_coefficient() {
        assert_eq!(&p(&[1, 2, 3]) * from_u64(3), p(&[3, 6, 9]));
        assert!((&p(&[1, 2, 3]) * Goldilocks::ZERO).is_zero());
    }

    #[test]
    fn divmod_exact_division() {
        // x^3 - 1 = (x - 1)(x^2 + x + 1)
        let x3_minus_1 = DensePoly::new(vec![
            -Goldilocks::ONE,
            Goldilocks::ZERO,
            Goldilocks::ZERO,
            Goldilocks::ONE,
        ]);
        let x_minus_1 = DensePoly::new(vec![-Goldilocks::ONE, Goldilocks::ONE]);
        let (q, r) = x3_minus_1.divmod(&x_minus_1);
        assert_eq!(q, p(&[1, 1, 1]));
        assert!(r.is_zero());
    }

    #[test]
    fn divmod_with_remainder_round_trips() {
        // From the Add test: 1 + 2x + 3x^2 = q*(4 + 5x) + r should round-trip
        let f = p(&[1, 2, 3]);
        let g = p(&[4, 5]);
        let (q, r) = f.divmod(&g);
        assert!(r.degree() < g.degree(), "deg(r) must be < deg(divisor)");
        assert_eq!(f, &(&q * &g) + &r);
    }

    #[test]
    fn divmod_smaller_dividend_returns_self() {
        let f = p(&[1, 2]);
        let g = p(&[1, 2, 3]);
        let (q, r) = f.divmod(&g);
        assert!(q.is_zero());
        assert_eq!(r, f);
    }

    #[test]
    #[should_panic(expected = "division by zero polynomial")]
    fn divmod_panics_on_zero_divisor() {
        let _ = p(&[1, 2, 3]).divmod(&DensePoly::zero());
    }
}
