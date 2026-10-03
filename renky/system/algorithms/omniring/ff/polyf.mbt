// ===========================================================================
//  有限域骨架 · part 1 — univariate polynomials over the prime field F_p.
//
//  Definition.  F_p[x] is the polynomial ring in one indeterminate over the
//  field with p elements.  A value of type `FPoly` is a coefficient vector
//  (low order first) reduced modulo p and *normalised*: trailing zero
//  coefficients are stripped, so the zero polynomial is the empty vector.
//  `degree = length - 1`, and `degree(0) = -1`.
//
//  Every algorithm here is self-contained: no third-party package, only
//  `kernel` (machine-checked modular arithmetic) and `arith` (gcd/primality).
// ===========================================================================

pub struct FPoly {
  p : Int
  c : Array[Int]
}

///| Smart constructor: reduces coefficients mod `p` and normalises.
pub fn fpoly(p : Int, cs : Array[Int]) -> FPoly {
  let out : Array[Int] = []
  for i = 0; i < cs.length(); i = i + 1 {
    out.push(@kernel.modpos(cs[i], p))
  }
  fp_strip(out)
  { p, c: out }
}

fn fp_strip(c : Array[Int]) -> Unit {
  while c.length() > 0 && c[c.length() - 1] == 0 {
    ignore(c.pop())
  }
}

pub fn fp_zero(p : Int) -> FPoly {
  { p, c: [] }
}

pub fn fp_one(p : Int) -> FPoly {
  { p, c: [1] }
}

///| `a * x^d`.
pub fn fp_monom(p : Int, d : Int, a : Int) -> FPoly {
  if a % p == 0 {
    return fp_zero(p)
  }
  let c : Array[Int] = Array::make(d + 1, 0)
  c[d] = @kernel.modpos(a, p)
  { p, c }
}

///| The indeterminate `x`.
pub fn fp_x(p : Int) -> FPoly {
  fp_monom(p, 1, 1)
}

pub fn FPoly::degree(self : FPoly) -> Int {
  self.c.length() - 1
}

pub fn FPoly::is_zero(self : FPoly) -> Bool {
  self.c.is_empty()
}

pub fn FPoly::is_monic(self : FPoly) -> Bool {
  !(self.is_zero()) && self.c[self.degree()] == 1
}

pub fn FPoly::coeff(self : FPoly, i : Int) -> Int {
  if i < 0 || i >= self.c.length() {
    0
  } else {
    self.c[i]
  }
}

pub fn FPoly::lead(self : FPoly) -> Int {
  if self.is_zero() {
    0
  } else {
    self.c[self.degree()]
  }
}

pub fn FPoly::equal(self : FPoly, other : FPoly) -> Bool {
  self.p == other.p && self.c.equal(other.c)
}

pub fn FPoly::to_string(self : FPoly) -> String {
  if self.is_zero() {
    return "0"
  }
  let parts : Array[String] = []
  let mut i = self.degree()
  while i >= 0 {
    let a = self.c[i]
    if a != 0 {
      let coef = if a == 1 && i > 0 {
        ""
      } else {
        a.to_string()
      }
      let ind = if i == 0 {
        ""
      } else if i == 1 {
        "x"
      } else {
        "x^" + i.to_string()
      }
      parts.push(coef + ind)
    }
    i = i - 1
  }
  parts.join(" + ")
}

///| **A18 polynomial addition / subtraction / multiplication over F_p.**
pub fn FPoly::add(self : FPoly, other : FPoly) -> FPoly {
  let n = Int::max(self.c.length(), other.c.length())
  let out : Array[Int] = Array::make(n, 0)
  for i = 0; i < n; i = i + 1 {
    out[i] = @kernel.addmod(self.coeff(i), other.coeff(i), self.p)
  }
  fp_strip(out)
  { p: self.p, c: out }
}

pub fn FPoly::sub(self : FPoly, other : FPoly) -> FPoly {
  let n = Int::max(self.c.length(), other.c.length())
  let out : Array[Int] = Array::make(n, 0)
  for i = 0; i < n; i = i + 1 {
    out[i] = @kernel.submod(self.coeff(i), other.coeff(i), self.p)
  }
  fp_strip(out)
  { p: self.p, c: out }
}

pub fn FPoly::neg(self : FPoly) -> FPoly {
  let out : Array[Int] = Array::make(self.c.length(), 0)
  for i = 0; i < out.length(); i = i + 1 {
    out[i] = @kernel.negmod(self.c[i], self.p)
  }
  { p: self.p, c: out }
}

pub fn FPoly::scale(self : FPoly, a : Int) -> FPoly {
  let out : Array[Int] = Array::make(self.c.length(), 0)
  for i = 0; i < out.length(); i = i + 1 {
    out[i] = @kernel.mulmod(self.c[i], a, self.p)
  }
  fp_strip(out)
  { p: self.p, c: out }
}

pub fn FPoly::mul(self : FPoly, other : FPoly) -> FPoly {
  let p = self.p
  if self.is_zero() || other.is_zero() {
    return fp_zero(p)
  }
  let n = self.degree() + other.degree() + 1
  let out : Array[Int] = Array::make(n, 0)
  for i = 0; i <= self.degree(); i = i + 1 {
    for j = 0; j <= other.degree(); j = j + 1 {
      out[i + j] = @kernel.addmod(
        out[i + j],
        @kernel.mulmod(self.c[i], other.c[j], p),
        p,
      )
    }
  }
  fp_strip(out)
  { p, c: out }
}

///| Division with remainder: `self = q*other + r`, `degree(r) < degree(other)`.
///| Requires `other != 0`; `p` prime makes the leading coefficient invertible.
pub fn FPoly::divmod(self : FPoly, other : FPoly) -> (FPoly, FPoly) {
  let p = self.p
  if other.is_zero() {
    abort("FPoly::divmod by zero")
  }
  let rem = self.c.copy()
  let lio = @arith.invmod(other.lead(), p)
  let d = other.degree()
  let q : Array[Int] = Array::make(
    Int::max(0, self.degree() - d + 1),
    0,
  )
  let mut k = rem.length() - 1
  while k >= d {
    let a = if k < rem.length() {
      rem[k]
    } else {
      0
    }
    if a == 0 {
      k = k - 1
      continue
    }
    let t = @kernel.mulmod(a, lio, p)
    if k - d < q.length() {
      q[k - d] = t
    } else {
      q.push(t)
    }
    for j = 0; j <= d; j = j + 1 {
      rem[k - d + j] = @kernel.submod(
        rem[k - d + j],
        @kernel.mulmod(t, other.c[j], p),
        p,
      )
    }
    k = k - 1
  }
  fp_strip(rem)
  fp_strip(q)
  ({ p, c: q }, { p, c: rem })
}

pub fn FPoly::div(self : FPoly, other : FPoly) -> FPoly {
  self.divmod(other).0
}

pub fn FPoly::rem(self : FPoly, other : FPoly) -> FPoly {
  self.divmod(other).1
}

///| Make a polynomial monic (divide by its leading coefficient).
pub fn FPoly::monic(self : FPoly) -> FPoly {
  if self.is_zero() {
    return self
  }
  self.scale(@arith.invmod(self.lead(), self.p))
}

///| Evaluation `f(a)`.
pub fn FPoly::eval(self : FPoly, a : Int) -> Int {
  let mut acc = 0
  let mut i = self.c.length() - 1
  while i >= 0 {
    acc = @kernel.addmod(@kernel.mulmod(acc, a, self.p), self.c[i], self.p)
    i = i - 1
  }
  acc
}

///| Formal derivative `f'`.
pub fn FPoly::deriv(self : FPoly) -> FPoly {
  if self.c.length() <= 1 {
    return fp_zero(self.p)
  }
  let out : Array[Int] = Array::make(self.c.length() - 1, 0)
  for i = 1; i < self.c.length(); i = i + 1 {
    out[i - 1] = @kernel.mulmod(self.c[i], i, self.p)
  }
  fp_strip(out)
  { p: self.p, c: out }
}

///| **A19 polynomial gcd** (Euclid in F_p[x]).
pub fn FPoly::gcd(self : FPoly, other : FPoly) -> FPoly {
  let mut a = self
  let mut b = other
  while !(b.is_zero()) {
    let r = a.rem(b)
    a = b
    b = r
  }
  a.monic()
}

///| Extended gcd: returns `(g, s, t)` with `s*self + t*other = g`, `g` monic.
pub fn fp_egcd(f : FPoly, g : FPoly) -> (FPoly, FPoly, FPoly) {
  let p = f.p
  let mut r0 = f
  let mut r1 = g
  let mut s0 = fp_one(p)
  let mut s1 = fp_zero(p)
  let mut t0 = fp_zero(p)
  let mut t1 = fp_one(p)
  while !(r1.is_zero()) {
    let (q, r) = r0.divmod(r1)
    let s2 = s0.sub(q.mul(s1))
    let t2 = t0.sub(q.mul(t1))
    r0 = r1
    r1 = r
    s0 = s1
    s1 = s2
    t0 = t1
    t1 = t2
  }
  if r0.is_zero() {
    return (r0, s0, t0)
  }
  let inv = @arith.invmod(r0.lead(), p)
  (r0.scale(inv), s0.scale(inv), t0.scale(inv))
}

///| `f^e mod m` — binary exponentiation in the quotient ring F_p[x]/(m).
///| **A20 poly powmod**.
pub fn fp_powmod(f : FPoly, e : Int, m : FPoly) -> FPoly {
  let p = f.p
  let mut acc = fp_one(p).rem(m)
  let mut base = f.rem(m)
  let mut rest = e
  while rest > 0 {
    if rest % 2 == 1 {
      acc = acc.mul(base).rem(m)
    }
    base = base.mul(base).rem(m)
    rest = rest / 2
  }
  acc
}

///| Product of `(x - r)` for every root in `rs`.
pub fn fp_from_roots(p : Int, rs : Array[Int]) -> FPoly {
  let mut acc = fp_one(p)
  for i = 0; i < rs.length(); i = i + 1 {
    acc = acc.mul(fpoly(p, [@kernel.negmod(rs[i], p), 1]))
  }
  acc
}
