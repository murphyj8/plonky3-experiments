use plonky3_experiments::field_arith::{GOLDILOCKS_PRIME, from_u64, to_canonical_u64};

fn main() {
    let a = from_u64(GOLDILOCKS_PRIME - 1);
    let b = from_u64(3);
    println!("Goldilocks prime p = {GOLDILOCKS_PRIME}");
    println!("(p - 1) + 3 mod p = {}", to_canonical_u64(a + b));
}
