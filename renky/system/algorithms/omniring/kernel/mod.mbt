// ---------------------------------------------------------------------------
// A01  modpos     canonical representative of a residue class
// A02  addmod     addition in Z/p
// A03  submod     subtraction in Z/p
// A04  mulmod     multiplication in Z/p
// A05  negmod     additive inverse in Z/p
//
// All five carry machine-checked contracts (`moon prove`, Why3 + Z3).
// The runtime precondition `1 <= p` is the *same* hypothesis the prover sees,
// so the executable ring and the proved ring are the same ring.
// ---------------------------------------------------------------------------

///| Canonical representative: `modpos a p ∈ [0,p)` and `modpos a p ≡ a (mod p)`.
///|
///| Definition.  `modpos(a,p) = a - p*⌊a/p⌋` for `p ≥ 1` (Euclidean remainder,
///| corrected to be non-negative).
pub fn modpos(a : Int, p : Int) -> Int where {
  proof_require: pos(p),
  proof_ensure: result => in_range(result, p),
  proof_ensure: result => cong(a, result, p),
} {
  let r = a % p
  // Euclidean division identity — this is the axiom `Div_mod` of
  // `int.ComputerDivision`; stating it explicitly lets the solver close the
  // congruence goal.
  proof_assert a == p * (a / p) + r
  proof_assert (a - r) % p == 0
  if r < 0 {
    r + p
  } else {
    r
  }
}

///| Addition in `Z/p`.  Theorem (proved): `addmod a b p ∈ [0,p)` and
///| `addmod a b p ≡ a + b (mod p)`.
pub fn addmod(a : Int, b : Int, p : Int) -> Int where {
  proof_require: pos(p),
  proof_ensure: result => in_range(result, p),
  proof_ensure: result => cong(a + b, result, p),
} {
  modpos(a + b, p)
}

///| Subtraction in `Z/p`.
pub fn submod(a : Int, b : Int, p : Int) -> Int where {
  proof_require: pos(p),
  proof_ensure: result => in_range(result, p),
  proof_ensure: result => cong(a - b, result, p),
} {
  modpos(a - b, p)
}

///| Multiplication in `Z/p`.
pub fn mulmod(a : Int, b : Int, p : Int) -> Int where {
  proof_require: pos(p),
  proof_ensure: result => in_range(result, p),
  proof_ensure: result => cong(a * b, result, p),
} {
  modpos(a * b, p)
}

///| Additive inverse in `Z/p`.
pub fn negmod(a : Int, p : Int) -> Int where {
  proof_require: pos(p),
  proof_ensure: result => in_range(result, p),
  proof_ensure: result => cong(0 - a, result, p),
} {
  modpos(0 - a, p)
}
