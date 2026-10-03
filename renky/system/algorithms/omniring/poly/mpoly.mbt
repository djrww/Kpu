// ===========================================================================
//  Multivariate polynomials F_p[x_0, ..., x_{n-1}] and Gröbner bases.
//
//  Definition.  A value of type `MPoly` is a sorted list of terms
//  `(exponent vector, coefficient)` with nonzero coefficients, ordered by the
//  graded-lexicographic order with x_0 > x_1 > ... :
//      u < v  iff  |u| < |v|, or |u| = |v| and u is lexicographically smaller
//      reading from the *last* variable.
//  `leading term` = the greatest monomial with its coefficient.
//
//  Theorem (Buchberger).  A finite set G generates an ideal I and is a Gröbner
//  basis of I iff every S-polynomial S(f,g), f,g in G, reduces to 0 modulo G.
//  `mp_groebner` implements exactly that criterion, so its output is a Gröbner
//  basis by construction; `ideal_wbtest` verifies the two characterisations
//  that matter computationally:
//      (a) `f in I  <=>  normal form of f w.r.t. G is 0`;
//      (b) the standard monomials (those not divisible by any leading monomial)
//          form an F_p-basis of F_p[x]/I.
// ===========================================================================

pub struct MPoly {
  p : Int
  nvars : Int
  terms : Array[(Array[Int], Int)]
}

///| Graded-lex comparison of exponent vectors (greatest last in the order we
///| sort by, so `sort_by(cmp_monom)` puts the leading term at the end).
pub fn cmp_monom(a : Array[Int], b : Array[Int]) -> Int {
  let da = total_deg(a)
  let db = total_deg(b)
  if da != db {
    return da - db
  }
  let mut i = a.length() - 1
  while i >= 0 {
    if a[i] != b[i] {
      return a[i] - b[i]
    }
    i = i - 1
  }
  0
}

pub fn total_deg(e : Array[Int]) -> Int {
  let mut s = 0
  for i = 0; i < e.length(); i = i + 1 {
    s = s + e[i]
  }
  s
}

fn exp_zero(nvars : Int) -> Array[Int] {
  Array::make(nvars, 0)
}

///| Normalising constructor: merges equal monomials mod p, drops zeros, sorts.
pub fn mp_make(p : Int, nvars : Int, terms : Array[(Array[Int], Int)]) -> MPoly {
  let sorted = terms.map(fn(t) { (t.0.copy(), @kernel.modpos(t.1, p)) })
  sorted.sort_by(fn(u, v) { cmp_monom(u.0, v.0) })
  let out : Array[(Array[Int], Int)] = []
  for i = 0; i < sorted.length(); i = i + 1 {
    let c = sorted[i].1
    if c == 0 {
      continue
    }
    if out.length() > 0 && cmp_monom(out[out.length() - 1].0, sorted[i].0) == 0 {
      let last = out[out.length() - 1]
      let sum = @kernel.addmod(last.1, c, p)
      if sum == 0 {
        ignore(out.pop())
      } else {
        out[out.length() - 1] = (last.0, sum)
      }
    } else {
      out.push((sorted[i].0, c))
    }
  }
  { p, nvars, terms: out }
}

pub fn mp_zero(p : Int, nvars : Int) -> MPoly {
  { p, nvars, terms: [] }
}

pub fn mp_one(p : Int, nvars : Int) -> MPoly {
  mp_make(p, nvars, [(exp_zero(nvars), 1)])
}

pub fn mp_const(p : Int, nvars : Int, c : Int) -> MPoly {
  if @kernel.modpos(c, p) == 0 {
    mp_zero(p, nvars)
  } else {
    mp_make(p, nvars, [(exp_zero(nvars), c)])
  }
}

///| The variable `x_i`.
pub fn mp_var(p : Int, nvars : Int, i : Int) -> MPoly {
  let e = exp_zero(nvars)
  e[i] = 1
  mp_make(p, nvars, [(e, 1)])
}

///| A single monomial `c * x^e`.
pub fn mp_monom(p : Int, nvars : Int, e : Array[Int], c : Int) -> MPoly {
  mp_make(p, nvars, [(e, c)])
}

pub fn MPoly::is_zero(self : MPoly) -> Bool {
  self.terms.is_empty()
}

pub fn MPoly::degree(self : MPoly) -> Int {
  if self.is_zero() {
    -1
  } else {
    total_deg(self.terms[self.terms.length() - 1].0)
  }
}

///| Leading monomial exponent vector (greatest w.r.t. the order).
pub fn MPoly::lm(self : MPoly) -> Array[Int] {
  if self.is_zero() {
    abort("MPoly::lm of zero")
  }
  self.terms[self.terms.length() - 1].0.copy()
}

pub fn MPoly::lc(self : MPoly) -> Int {
  if self.is_zero() {
    0
  } else {
    self.terms[self.terms.length() - 1].1
  }
}

pub fn MPoly::lt(self : MPoly) -> (Array[Int], Int) {
  (self.lm(), self.lc())
}

pub fn MPoly::nterms(self : MPoly) -> Int {
  self.terms.length()
}

pub fn MPoly::equal(self : MPoly, other : MPoly) -> Bool {
  if self.p != other.p || self.nvars != other.nvars {
    return false
  }
  if self.terms.length() != other.terms.length() {
    return false
  }
  for i = 0; i < self.terms.length(); i = i + 1 {
    if self.terms[i].1 != other.terms[i].1 {
      return false
    }
    if cmp_monom(self.terms[i].0, other.terms[i].0) != 0 {
      return false
    }
  }
  true
}

pub fn MPoly::to_string(self : MPoly) -> String {
  if self.is_zero() {
    return "0"
  }
  let parts : Array[String] = []
  for i = self.terms.length() - 1; i >= 0; i = i - 1 {
    let (e, c) = self.terms[i]
    let cs = if c == 1 && total_deg(e) > 0 {
      ""
    } else {
      c.to_string()
    }
    let vs : Array[String] = []
    for j = 0; j < e.length(); j = j + 1 {
      if e[j] == 1 {
        vs.push("x" + j.to_string())
      } else if e[j] > 1 {
        vs.push("x" + j.to_string() + "^" + e[j].to_string())
      }
    }
    parts.push(cs + vs.join("*"))
  }
  parts.join(" + ")
}

pub fn MPoly::add(self : MPoly, other : MPoly) -> MPoly {
  let all : Array[(Array[Int], Int)] = []
  for i = 0; i < self.terms.length(); i = i + 1 {
    all.push(self.terms[i])
  }
  for i = 0; i < other.terms.length(); i = i + 1 {
    all.push(other.terms[i])
  }
  mp_make(self.p, self.nvars, all)
}

pub fn MPoly::sub(self : MPoly, other : MPoly) -> MPoly {
  let all : Array[(Array[Int], Int)] = []
  for i = 0; i < self.terms.length(); i = i + 1 {
    all.push(self.terms[i])
  }
  for i = 0; i < other.terms.length(); i = i + 1 {
    all.push((other.terms[i].0, @kernel.negmod(other.terms[i].1, self.p)))
  }
  mp_make(self.p, self.nvars, all)
}

pub fn MPoly::neg(self : MPoly) -> MPoly {
  mp_make(
    self.p,
    self.nvars,
    self.terms.map(fn(t) { (t.0, @kernel.negmod(t.1, self.p)) }),
  )
}

pub fn MPoly::scale(self : MPoly, c : Int) -> MPoly {
  mp_make(
    self.p,
    self.nvars,
    self.terms.map(fn(t) { (t.0, @kernel.mulmod(t.1, c, self.p)) }),
  )
}

pub fn MPoly::mul(self : MPoly, other : MPoly) -> MPoly {
  let all : Array[(Array[Int], Int)] = []
  for i = 0; i < self.terms.length(); i = i + 1 {
    for j = 0; j < other.terms.length(); j = j + 1 {
      let e = exp_zero(self.nvars)
      for k = 0; k < self.nvars; k = k + 1 {
        e[k] = self.terms[i].0[k] + other.terms[j].0[k]
      }
      all.push((e, @kernel.mulmod(self.terms[i].1, other.terms[j].1, self.p)))
    }
  }
  mp_make(self.p, self.nvars, all)
}

///| Evaluation at a point of F_p^n.
pub fn MPoly::eval(self : MPoly, pt : Array[Int]) -> Int {
  let mut acc = 0
  for i = 0; i < self.terms.length(); i = i + 1 {
    let (e, c) = self.terms[i]
    let mut mono = c
    for j = 0; j < self.nvars; j = j + 1 {
      mono = @kernel.mulmod(mono, @arith.powmod(pt[j], e[j], self.p), self.p)
    }
    acc = @kernel.addmod(acc, mono, self.p)
  }
  acc
}

///| Does monomial `a` divide monomial `b`?
pub fn monom_divides(a : Array[Int], b : Array[Int]) -> Bool {
  for i = 0; i < a.length(); i = i + 1 {
    if a[i] > b[i] {
      return false
    }
  }
  true
}

pub fn monom_div(a : Array[Int], b : Array[Int]) -> Array[Int] {
  let out = exp_zero(a.length())
  for i = 0; i < a.length(); i = i + 1 {
    out[i] = b[i] - a[i]
  }
  out
}

pub fn monom_lcm(a : Array[Int], b : Array[Int]) -> Array[Int] {
  let out = exp_zero(a.length())
  for i = 0; i < a.length(); i = i + 1 {
    out[i] = Int::max(a[i], b[i])
  }
  out
}

///| **A31 multivariate division with remainder** by an ordered list of
///| divisors.  Returns `(quotients, remainder)`.
pub fn mp_divmod(f : MPoly, gs : Array[MPoly]) -> (Array[MPoly], MPoly) {
  let p = f.p
  let nvars = f.nvars
  let qs : Array[MPoly] = []
  for i = 0; i < gs.length(); i = i + 1 {
    qs.push(mp_zero(p, nvars))
  }
  let mut rem = mp_zero(p, nvars)
  let mut rest = f
  while !(rest.is_zero()) {
    let (e, c) = rest.lt()
    let mut divided = false
    let mut i = gs.length() - 1
    while i >= 0 && !(divided) {
      if !(gs[i].is_zero()) {
        let ge = gs[i].lm()
        if monom_divides(ge, e) {
          let q_e = monom_div(ge, e)
          let q_c = @kernel.mulmod(c, @arith.invmod(gs[i].lc(), p), p)
          let q = mp_monom(p, nvars, q_e, q_c)
          qs[i] = qs[i].add(q)
          rest = rest.sub(q.mul(gs[i]))
          divided = true
        }
      }
      i = i - 1
    }
    if !(divided) {
      rem = rem.add(mp_monom(p, nvars, e, c))
      rest = rest.sub(mp_monom(p, nvars, e, c))
    }
  }
  (qs, rem)
}

///| Normal form of `f` modulo the list `gs` (a Gröbner basis makes this
///| canonical, i.e. independent of choices).
pub fn mp_normal_form(f : MPoly, gs : Array[MPoly]) -> MPoly {
  mp_divmod(f, gs).1
}

///| **A32 S-polynomial** `S(f,g) = (L/lc(f) LT(f)) * f - (L/lc(g) LT(g)) * g`
///| with `L = lcm(LM f, LM g)`.
pub fn mp_spoly(f : MPoly, g : MPoly) -> MPoly {
  let p = f.p
  let nvars = f.nvars
  if f.is_zero() || g.is_zero() {
    return mp_zero(p, nvars)
  }
  let l = monom_lcm(f.lm(), g.lm())
  let ef = monom_div(f.lm(), l)
  let eg = monom_div(g.lm(), l)
  let cf = @kernel.mulmod(1, @arith.invmod(f.lc(), p), p)
  let cg = @kernel.mulmod(1, @arith.invmod(g.lc(), p), p)
  mp_monom(p, nvars, ef, cf).mul(f).sub(mp_monom(p, nvars, eg, cg).mul(g))
}

///| **A33 Buchberger's algorithm.**  Returns a Gröbner basis of the ideal
///| generated by `gs` (reduced if `reduce_it` is set).
pub fn mp_groebner(gs : Array[MPoly], reduce_it : Bool) -> Array[MPoly] {
  if gs.is_empty() {
    return gs
  }
  let p = gs[0].p
  let basis : Array[MPoly] = []
  for i = 0; i < gs.length(); i = i + 1 {
    if !(gs[i].is_zero()) {
      basis.push(gs[i].scale(@arith.invmod(gs[i].lc(), p)))
    }
  }
  let mut changed = true
  let mut rounds = 0
  while changed && rounds < 400 {
    changed = false
    rounds = rounds + 1
    let n = basis.length()
    let pairs : Array[(Int, Int)] = []
    for i = 0; i < n; i = i + 1 {
      for j = i + 1; j < n; j = j + 1 {
        pairs.push((i, j))
      }
    }
    for k = 0; k < pairs.length(); k = k + 1 {
      let i = pairs[k].0
      let j = pairs[k].1
      let s = mp_spoly(basis[i], basis[j])
      let r = mp_normal_form(s, basis)
      if !(r.is_zero()) {
        basis.push(r.scale(@arith.invmod(r.lc(), p)))
        changed = true
      }
    }
  }
  if reduce_it {
    mp_reduce_basis(basis)
  } else {
    basis
  }
}

///| Turn a Gröbner basis into the *reduced* one (monic, and no term of any
///| element is divisible by the leading monomial of another).
pub fn mp_reduce_basis(basis : Array[MPoly]) -> Array[MPoly] {
  if basis.is_empty() {
    return basis
  }
  let p = basis[0].p
  let mut cur : Array[MPoly] = []
  for i = 0; i < basis.length(); i = i + 1 {
    if !(basis[i].is_zero()) {
      cur.push(basis[i].scale(@arith.invmod(basis[i].lc(), p)))
    }
  }
  let mut changed = true
  let mut rounds = 0
  while changed && rounds < 50 {
    changed = false
    rounds = rounds + 1
    // (1) drop an element when an *earlier* one already has a leading monomial
    //     dividing it.  Keeping the first representative is what makes the
    //     filter correct on repeated leading monomials.
    let kept : Array[MPoly] = []
    for i = 0; i < cur.length(); i = i + 1 {
      let mut drop = false
      for j = 0; j < kept.length(); j = j + 1 {
        if monom_divides(kept[j].lm(), cur[i].lm()) {
          drop = true
        }
      }
      if drop {
        changed = true
      } else {
        kept.push(cur[i])
      }
    }
    // (2) reduce the tail of each element against the others.  After (1) the
    //     leading monomials are pairwise non-dividing, so leading terms are
    //     preserved and this converges immediately.
    let out : Array[MPoly] = []
    for i = 0; i < kept.length(); i = i + 1 {
      let others : Array[MPoly] = []
      for j = 0; j < kept.length(); j = j + 1 {
        if i != j {
          others.push(kept[j])
        }
      }
      let r = mp_normal_form(kept[i], others)
      if r.is_zero() {
        changed = true
        continue
      }
      let rn = r.scale(@arith.invmod(r.lc(), p))
      if !(rn.equal(kept[i])) {
        changed = true
      }
      out.push(rn)
    }
    cur = out
  }
  cur.sort_by(fn(u, v) { cmp_monom(u.lm(), v.lm()) })
  cur
}

///| All monomials of total degree <= `maxdeg` in `nvars` variables, ordered.
pub fn mp_all_monomials(nvars : Int, maxdeg : Int) -> Array[Array[Int]] {
  let out : Array[Array[Int]] = []
  enumerate_monomials(nvars, maxdeg, exp_zero(nvars), 0, out)
  out.sort_by(cmp_monom)
  out
}

fn enumerate_monomials(
  nvars : Int,
  budget : Int,
  cur : Array[Int],
  v : Int,
  out : Array[Array[Int]],
) -> Unit {
  if v == nvars {
    out.push(cur.copy())
    return
  }
  for k = 0; k <= budget; k = k + 1 {
    cur[v] = k
    enumerate_monomials(nvars, budget - k, cur, v + 1, out)
  }
  cur[v] = 0
}

///| **A34 standard monomials** of a Gröbner basis: those not divisible by any
///| leading monomial.  They are an F_p-basis of F_p[x]/I when finite.
pub fn mp_standard_monomials(
  basis : Array[MPoly],
  nvars : Int,
  maxdeg : Int,
) -> (Array[Array[Int]], Bool) {
  let all = mp_all_monomials(nvars, maxdeg)
  let out : Array[Array[Int]] = []
  for i = 0; i < all.length(); i = i + 1 {
    let mut blocked = false
    for j = 0; j < basis.length(); j = j + 1 {
      if !(basis[j].is_zero()) && monom_divides(basis[j].lm(), all[i]) {
        blocked = true
      }
    }
    if !(blocked) {
      out.push(all[i])
    }
  }
  // finiteness: no standard monomial on the top degree shell
  let mut finite = true
  for i = 0; i < out.length(); i = i + 1 {
    if total_deg(out[i]) == maxdeg {
      finite = false
    }
  }
  (out, finite)
}
