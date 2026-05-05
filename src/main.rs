use plonky3_experiments::field_arith::{GOLDILOCKS_PRIME, from_u64};
use plonky3_experiments::poly::DensePoly;

fn main() {
    println!("Goldilocks prime p = {GOLDILOCKS_PRIME}");

    // f(x) = 1 + 2x + 3x^2 ; g(x) = 4 + 5x
    let f = DensePoly::new(vec![from_u64(1), from_u64(2), from_u64(3)]);
    let g = DensePoly::new(vec![from_u64(4), from_u64(5)]);
    let x = from_u64(7);

    println!("f({x})           = {}", f.evaluate(x));
    println!("(f + g)({x})     = {}", (&f + &g).evaluate(x));
    println!("(f - g)({x})     = {}", (&f - &g).evaluate(x));
    println!("(-g)({x})        = {}", (-&g).evaluate(x));
    println!("(f * g)({x})     = {}", (&f * &g).evaluate(x));
    println!("(10 * f)({x})    = {}", (&f * from_u64(10)).evaluate(x));

    let (q, r) = f.divmod(&g);
    println!("f / g => deg(q) = {:?}, deg(r) = {:?}", q.degree(), r.degree());

    let reconstructed = &(&q * &g) + &r;
    assert_eq!(f, reconstructed, "divmod must satisfy f == q*g + r");
    println!("verified: f == q*g + r");
}
