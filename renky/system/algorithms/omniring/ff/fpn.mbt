// ===========================================================================
//  有限域骨架 · part 3 — the extension field F_{p^n} = F_p[x]/(f),
//  with `f` a monic irreducible of degree n produced by Ben-Or's test.
//
//  Definition.  Elements are coefficient vectors of length <= n-1; addition is
//  coefficientwise mod p, multiplication is polynomial multiplication reduced
//  modulo f.  Because f is irreducible, F_p[x]/(f) is a field of size p^n.
//
//  Theorem (Wedderburn).  The multiplicative group of a finite division ring
//  is cyclic.  `ff_wbtest` exhibits an explicit generator for several
//  F_{p^n} and verifies that its order is p^n - 1.
// ===========================================================================

pub struct Fpn {
  p : Int
  n : Int
  irr : FPoly
  c : Array[Int]
}

fn fpn_wrap(p : Int, n : Int, irr : FPoly, c : Array[Int]) -> Fpn {
  let out : Array[Int] = Array::make(n, 0)
  for i = 0; i < n; i = i + 1 {
    out[i] = @kernel.modpos(if i < c.length() {
      c[i]
    } else {
      0
    }, p)
  }
  { p, n, irr, c: out }
}

///| Build F_{p^n}: searches the first monic irreducible of degree n over F_p.
pub fn fpn_create(p : Int, n : Int) -> Fpn {
  let f = fp_find_irreducible(p, n)
  fpn_wrap(p, n, f, if n == 1 {
    [1]
  } else {
    []
  })
}

///| Build F_{p^n} from a given monic irreducible.
pub fn fpn_of(p : Int, n : Int, irr : FPoly) -> Fpn {
  fpn_wrap(p, n, irr, [])
}

pub fn Fpn::zero(self : Fpn) -> Fpn {
  fpn_wrap(self.p, self.n, self.irr, [])
}

pub fn Fpn::one(self : Fpn) -> Fpn {
  fpn_wrap(self.p, self.n, self.irr, [1])
}

///| The image of `x` in F_p[x]/(f): a root of `f` inside the field.
pub fn Fpn::gen(self : Fpn) -> Fpn {
  fpn_wrap(self.p, self.n, self.irr, if self.n == 1 {
    [0]
  } else {
    [0, 1]
  })
}

///| Embed F_p into F_{p^n}.
pub fn Fpn::lift(self : Fpn, a : Int) -> Fpn {
  fpn_wrap(self.p, self.n, self.irr, [@kernel.modpos(a, self.p)])
}

pub fn Fpn::coeffs(self : Fpn) -> Array[Int] {
  self.c.copy()
}

pub fn Fpn::equal(self : Fpn, other : Fpn) -> Bool {
  self.p == other.p && self.n == other.n && self.c.equal(other.c)
}

pub fn Fpn::is_zero(self : Fpn) -> Bool {
  self.c.all(fn(v) { v == 0 })
}

pub fn Fpn::to_string(self : Fpn) -> String {
  let parts : Array[String] = []
  for i = 0; i < self.n; i = i + 1 {
    if self.c[i] != 0 {
      let ind = if i == 0 {
        ""
      } else if i == 1 {
        "w"
      } else {
        "w^" + i.to_string()
      }
      parts.push(self.c[i].to_string() + ind)
    }
  }
  if parts.is_empty() {
    "0"
  } else {
    parts.join(" + ")
  }
}

pub fn Fpn::add(self : Fpn, other : Fpn) -> Fpn {
  let out : Array[Int] = Array::make(self.n, 0)
  for i = 0; i < self.n; i = i + 1 {
    out[i] = @kernel.addmod(self.c[i], other.c[i], self.p)
  }
  fpn_wrap(self.p, self.n, self.irr, out)
}

pub fn Fpn::sub(self : Fpn, other : Fpn) -> Fpn {
  let out : Array[Int] = Array::make(self.n, 0)
  for i = 0; i < self.n; i = i + 1 {
    out[i] = @kernel.submod(self.c[i], other.c[i], self.p)
  }
  fpn_wrap(self.p, self.n, self.irr, out)
}

pub fn Fpn::neg(self : Fpn) -> Fpn {
  let out : Array[Int] = Array::make(self.n, 0)
  for i = 0; i < self.n; i = i + 1 {
    out[i] = @kernel.negmod(self.c[i], self.p)
  }
  fpn_wrap(self.p, self.n, self.irr, out)
}

///| Multiplication: polynomial product reduced modulo the irreducible.
pub fn Fpn::mul(self : Fpn, other : Fpn) -> Fpn {
  let a = fpoly(self.p, self.c)
  let b = fpoly(self.p, other.c)
  let r = a.mul(b).rem(self.irr)
  fpn_wrap(self.p, self.n, self.irr, r.c)
}

///| Inverse via the extended Euclidean algorithm in F_p[x]:
///| since `f` is irreducible and `a != 0`, `gcd(a,f)=1`, so `s*a + t*f = 1`
///| and `s` is the inverse of `a` modulo `f`.
pub fn Fpn::inv(self : Fpn) -> Fpn {
  if self.is_zero() {
    abort("Fpn::inv of zero")
  }
  let a = fpoly(self.p, self.c)
  let (_, s, _) = fp_egcd(a, self.irr)
  fpn_wrap(self.p, self.n, self.irr, s.rem(self.irr).c)
}

pub fn Fpn::div(self : Fpn, other : Fpn) -> Fpn {
  self.mul(other.inv())
}

pub fn Fpn::pow(self : Fpn, e : Int) -> Fpn {
  let mut acc = self.one()
  let mut base = self
  let mut rest = e
  while rest > 0 {
    if rest % 2 == 1 {
      acc = acc.mul(base)
    }
    base = base.mul(base)
    rest = rest / 2
  }
  acc
}

///| Frobenius automorphism `F(x) = x^p`.
pub fn Fpn::frobenius(self : Fpn) -> Fpn {
  self.pow(self.p)
}

///| `F^k`, the k-th iterate of Frobenius.
pub fn Fpn::frobenius_k(self : Fpn, k : Int) -> Fpn {
  let mut acc = self
  for _i = 0; _i < k; _i = _i + 1 {
    acc = acc.frobenius()
  }
  acc
}

///| Trace to F_p: `Tr(x) = x + x^p + ... + x^{p^{n-1}}`.
pub fn Fpn::trace(self : Fpn) -> Int {
  let mut acc = self.zero()
  let mut cur = self
  for _i = 0; _i < self.n; _i = _i + 1 {
    acc = acc.add(cur)
    cur = cur.frobenius()
  }
  acc.c[0]
}

///| Norm to F_p: `N(x) = x * x^p * ... * x^{p^{n-1}} = x^{(p^n-1)/(p-1)}`.
pub fn Fpn::norm(self : Fpn) -> Int {
  let mut acc = self.one()
  let mut cur = self
  for _i = 0; _i < self.n; _i = _i + 1 {
    acc = acc.mul(cur)
    cur = cur.frobenius()
  }
  acc.c[0]
}

///| The Galois conjugates `x, x^p, ..., x^{p^{n-1}}`.
pub fn Fpn::conjugates(self : Fpn) -> Array[Fpn] {
  let out : Array[Fpn] = []
  let mut cur = self
  for _i = 0; _i < self.n; _i = _i + 1 {
    out.push(cur)
    cur = cur.frobenius()
  }
  out
}

///| Minimal polynomial of `x` over F_p, computed as
///| `prod_{i=0}^{orbit-1} (X - F^i(x))` **inside** F_{p^n}[X].
///| Galois invariance forces every coefficient back into F_p; that is checked
///| on the fly, so the function is simultaneously the computation and a test of
///| the theorem.
pub fn Fpn::minpoly(self : Fpn) -> FPoly {
  let p = self.p
  let mut acc : Array[Fpn] = [self.one()]
  let mut cur = self
  let mut steps = 0
  while steps < self.n {
    let nxt : Array[Fpn] = []
    for i = 0; i <= acc.length(); i = i + 1 {
      nxt.push(self.zero())
    }
    for j = 0; j < acc.length(); j = j + 1 {
      nxt[j + 1] = nxt[j + 1].add(acc[j])
      nxt[j] = nxt[j].sub(acc[j].mul(cur))
    }
    acc = nxt
    cur = cur.frobenius()
    steps = steps + 1
    if cur.equal(self) {
      steps = self.n
    }
  }
  let cs : Array[Int] = Array::make(acc.length(), 0)
  for j = 0; j < acc.length(); j = j + 1 {
    let co = acc[j].coeffs()
    for t = 1; t < co.length(); t = t + 1 {
      if co[t] != 0 {
        abort("minpoly: coefficient escaped F_p (Galois invariance failed)")
      }
    }
    cs[j] = co[0]
  }
  fpoly(p, cs)
}

///| Multiplicative order of a nonzero element (divides p^n - 1).
pub fn Fpn::order(self : Fpn) -> Int {
  if self.is_zero() {
    return 0
  }
  let q = @arith.pow_int(self.p, self.n) - 1
  let ds = @arith.divisors(q)
  for i = 0; i < ds.length(); i = i + 1 {
    if self.pow(ds[i]).equal(self.one()) {
      return ds[i]
    }
  }
  q
}

///| A generator of the cyclic group F_{p^n}^* (Wedderburn's theorem, made
///| constructive by search).
pub fn Fpn::primitive_element(self : Fpn) -> Fpn {
  let q = @arith.pow_int(self.p, self.n) - 1
  let g = self.gen()
  if !(g.is_zero()) && g.order() == q {
    return g
  }
  // brute force over all elements in lexicographic order of coefficients
  let total = q + 1
  for k = 1; k < total; k = k + 1 {
    let el = fpn_from_index(self, k)
    if el.order() == q {
      return el
    }
  }
  abort("Fpn::primitive_element: none found (Wedderburn contradicted)")
}

///| The k-th element of F_{p^n} in mixed-radix (base p) order.
pub fn fpn_from_index(fld : Fpn, k : Int) -> Fpn {
  let ds = @arith.digits(k, fld.p)
  fpn_wrap(fld.p, fld.n, fld.irr, ds)
}

///| Number of elements, p^n.
pub fn Fpn::size(self : Fpn) -> Int {
  @arith.pow_int(self.p, self.n)
}

///| All elements (only for small fields; used by exhaustive self-proofs).
pub fn Fpn::elements(self : Fpn) -> Array[Fpn] {
  let out : Array[Fpn] = []
  for k = 0; k < self.size(); k = k + 1 {
    out.push(fpn_from_index(self, k))
  }
  out
}

///| Discrete logarithm base `g` by baby-step giant-step, using the
///| *machine-checked* binary search from package `select`/`kernel`.
///| **A21 BSGS discrete log**.  Returns -1 when `h` is not in `<g>`.
pub fn Fpn::dlog(self : Fpn, g : Fpn, h : Fpn) -> Int {
  let ord = g.order()
  if ord == 0 {
    return -1
  }
  let m = isqrt_ceil(ord)
  // baby steps: table[j] = g^j, j < m ; kept sorted by encoding elements
  let keys : FixedArray[Int] = FixedArray::make(m, 0)
  let idxs : FixedArray[Int] = FixedArray::make(m, 0)
  let mut cur = self.one()
  for j = 0; j < m; j = j + 1 {
    keys[j] = fpn_encode(cur)
    idxs[j] = j
    cur = cur.mul(g)
  }
  // sort (keys, idxs) by keys
  sort_pairs(keys, idxs)
  // giant steps: h * (g^{-m})^i
  let gm = g.pow(m).inv()
  let mut gamma = h
  for i = 0; i <= m; i = i + 1 {
    let target = fpn_encode(gamma)
    match @select.bsearch(keys, target) {
      Some(pos) => {
        let cand = i * m + idxs[pos]
        if cand < ord && g.pow(cand).equal(h) {
          return cand
        }
      }
      None => ()
    }
    gamma = gamma.mul(gm)
  }
  -1
}

///| Injective encoding of a field element into an Int (mixed radix, base p).
pub fn fpn_encode(el : Fpn) -> Int {
  @arith.from_digits(el.c, el.p)
}

fn sort_pairs(keys : FixedArray[Int], idxs : FixedArray[Int]) -> Unit {
  // insertion sort: m is small (<= a few thousand) and this keeps the pair
  // association explicit without allocating tuples.
  for i = 1; i < keys.length(); i = i + 1 {
    let k = keys[i]
    let v = idxs[i]
    let mut j = i
    while j > 0 && keys[j - 1] > k {
      keys[j] = keys[j - 1]
      idxs[j] = idxs[j - 1]
      j = j - 1
    }
    keys[j] = k
    idxs[j] = v
  }
}

fn isqrt_ceil(n : Int) -> Int {
  let mut r = 1
  while r * r < n {
    r = r + 1
  }
  r
}
