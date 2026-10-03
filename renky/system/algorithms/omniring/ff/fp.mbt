// ===========================================================================
//  有限域骨架 · part 2 — the prime field F_p.
//
//  Definition.  For a prime p, F_p = Z/pZ with the operations
//      a + b := modpos(a+b, p),   a * b := modpos(a*b, p).
//  Theorem (machine-checked in `kernel`).  These operations make F_p a field:
//  every value has a canonical representative in [0,p), addition and
//  multiplication are congruent to the integer ones, and every nonzero element
//  has an inverse (proved exhaustively in `ff_wbtest`).
// ===========================================================================

///| An element of the prime field F_p.
pub struct Fp {
  p : Int
  v : Int
}

///| Canonical injection Z -> F_p.
pub fn fp_mk(p : Int, v : Int) -> Fp {
  { p, v: @kernel.modpos(v, p) }
}

pub fn Fp::zero(p : Int) -> Fp {
  { p, v: 0 }
}

pub fn Fp::one(p : Int) -> Fp {
  { p, v: @kernel.modpos(1, p) }
}

pub fn Fp::char(self : Fp) -> Int {
  self.p
}

pub fn Fp::repr(self : Fp) -> Int {
  self.v
}

pub fn Fp::is_zero(self : Fp) -> Bool {
  self.v == 0
}

pub fn Fp::equal(self : Fp, other : Fp) -> Bool {
  self.p == other.p && self.v == other.v
}

pub fn Fp::to_string(self : Fp) -> String {
  self.v.to_string() + "/" + self.p.to_string()
}

pub fn Fp::add(self : Fp, other : Fp) -> Fp {
  { p: self.p, v: @kernel.addmod(self.v, other.v, self.p) }
}

pub fn Fp::sub(self : Fp, other : Fp) -> Fp {
  { p: self.p, v: @kernel.submod(self.v, other.v, self.p) }
}

pub fn Fp::mul(self : Fp, other : Fp) -> Fp {
  { p: self.p, v: @kernel.mulmod(self.v, other.v, self.p) }
}

pub fn Fp::neg(self : Fp) -> Fp {
  { p: self.p, v: @kernel.negmod(self.v, self.p) }
}

///| Multiplicative inverse; aborts on zero (a field has no inverse for 0).
pub fn Fp::inv(self : Fp) -> Fp {
  if self.v == 0 {
    abort("Fp::inv of zero")
  }
  { p: self.p, v: @arith.invmod(self.v, self.p) }
}

pub fn Fp::div(self : Fp, other : Fp) -> Fp {
  self.mul(other.inv())
}

///| `self^e` for `e >= 0` (Fermat reduction is *not* applied, so this is the
///| literal power map).
pub fn Fp::pow(self : Fp, e : Int) -> Fp {
  { p: self.p, v: @arith.powmod(self.v, e, self.p) }
}

///| The Frobenius endomorphism `x |-> x^p`.  On F_p it is the identity
///| (Fermat's little theorem) — checked exhaustively in `ff_wbtest`.
pub fn Fp::frobenius(self : Fp) -> Fp {
  self.pow(self.p)
}

///| All elements of F_p, in canonical order.  A finite field is *enumerable*,
///| which is what makes exhaustive self-proof possible.
pub fn fp_elements(p : Int) -> Array[Fp] {
  let out : Array[Fp] = []
  for v = 0; v < p; v = v + 1 {
    out.push({ p, v })
  }
  out
}
