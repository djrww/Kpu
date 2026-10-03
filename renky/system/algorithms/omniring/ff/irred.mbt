// ===========================================================================
//  有限域骨架 · part 4 — irreducibility and factorisation in F_p[x].
//
//  Definition.  f ∈ F_p[x] monic of degree n is *irreducible* iff it has no
//  factor of degree in [1, n/2].
//
//  Theorem (Ben-Or).  f is irreducible iff gcd(f, x^{p^i} - x) = 1 for every
//  1 <= i <= n/2.  (Because x^{p^i} - x is the product of all monic
//  irreducibles whose degree divides i.)
//
//  Theorem (Gauss).  The number of monic irreducibles of degree d over F_p is
//      N_p(d) = (1/d) * sum_{k | d} mu(d/k) p^k.
//  `ff_wbtest` verifies this formula against a complete enumeration.
// ===========================================================================

///| `x^(p^k) mod f`, computed by iterating the Frobenius on polynomials.
pub fn fp_frobenius_pow(x : FPoly, k : Int, f : FPoly) -> FPoly {
  let mut acc = x
  for _i = 0; _i < k; _i = _i + 1 {
    acc = fp_powmod(acc, f.p, f)
  }
  acc
}

///| **A22 Ben-Or irreducibility test.**
pub fn fp_is_irreducible(f : FPoly) -> Bool {
  let n = f.degree()
  if n <= 0 {
    return false
  }
  if !(f.is_monic()) {
    return false
  }
  if n == 1 {
    return true
  }
  let x = fp_x(f.p)
  let mut h = fp_powmod(x, f.p, f)
  for i = 1; i <= n / 2; i = i + 1 {
    let g = f.gcd(h.sub(x))
    if g.degree() > 0 {
      return false
    }
    if i < n / 2 {
      h = fp_powmod(h, f.p, f)
    }
  }
  true
}

///| First monic irreducible of degree `n` over F_p in lexicographic order.
///| **A23 irreducible search.**
pub fn fp_find_irreducible(p : Int, n : Int) -> FPoly {
  if n <= 0 {
    abort("fp_find_irreducible: n must be positive")
  }
  let total = @arith.pow_int(p, n)
  for k = 0; k < total; k = k + 1 {
    let cs = @arith.digits(k, p)
    let c : Array[Int] = Array::make(n + 1, 0)
    for i = 0; i < n; i = i + 1 {
      c[i] = if i < cs.length() {
        cs[i]
      } else {
        0
      }
    }
    c[n] = 1
    let f = fpoly(p, c)
    if fp_is_irreducible(f) {
      return f
    }
  }
  abort("fp_find_irreducible: existence theorem contradicted")
}

///| Gauss's formula for the number of monic irreducibles of degree d.
pub fn fp_count_monic_irred(p : Int, d : Int) -> Int {
  let ds = @arith.divisors(d)
  let mut s = 0
  for i = 0; i < ds.length(); i = i + 1 {
    let k = ds[i]
    s = s + @arith.mobius(d / k) * @arith.pow_int(p, k)
  }
  s / d
}

///| Every monic irreducible of degree `d` over F_p, by complete enumeration.
pub fn fp_all_irred(p : Int, d : Int) -> Array[FPoly] {
  let out : Array[FPoly] = []
  let total = @arith.pow_int(p, d)
  for k = 0; k < total; k = k + 1 {
    let cs = @arith.digits(k, p)
    let c : Array[Int] = Array::make(d + 1, 0)
    for i = 0; i < d; i = i + 1 {
      c[i] = if i < cs.length() {
        cs[i]
      } else {
        0
      }
    }
    c[d] = 1
    let f = fpoly(p, c)
    if fp_is_irreducible(f) {
      out.push(f)
    }
  }
  out
}

///| **A24 distinct-degree factorisation.**  Input: monic squarefree `f`.
///| Output: pairs `(d, g_d)` where `g_d` is the product of all monic
///| irreducible factors of `f` of degree `d`.
pub fn fp_ddf(f : FPoly) -> Array[(Int, FPoly)] {
  let out : Array[(Int, FPoly)] = []
  let p = f.p
  let mut rest = f.monic()
  let x = fp_x(p)
  let mut g = x.rem(rest)
  let mut i = 1
  while rest.degree() > 0 && 2 * i <= rest.degree() {
    g = fp_powmod(g, p, rest)
    let h = rest.gcd(g.sub(x).rem(rest))
    if h.degree() > 0 {
      out.push((i, h))
      rest = rest.div(h)
      g = g.rem(rest)
    }
    i = i + 1
  }
  if rest.degree() > 0 {
    out.push((rest.degree(), rest))
  }
  out
}

///| **A25 Berlekamp factorisation** (deterministic; replaces the randomised
///| Cantor–Zassenhaus so that the whole library stays deterministic).
///|
///| Theorem (Berlekamp).  For `f` monic squarefree of degree `n`, the F_p-algebra
///| `A = F_p[x]/(f)` satisfies `dim_{F_p} {h in A : h^p = h} = r`, where `r` is
///| the number of distinct irreducible factors of `f`.  Each such `h` and each
///| `a in F_p` gives a factor `gcd(f, h - a)`, and together they split `f`.
///|
///| The kernel of `Q - I` (with `Q[i][j] = [x^i] x^{jp} mod f`) is computed by
///| the machine-tested `linalg.nullspace`, so this algorithm reuses verified
///| linear algebra instead of trusting a probabilistic split.
pub fn fp_berlekamp_kernel(f : FPoly) -> Array[Array[Int]] {
  let p = f.p
  let n = f.degree()
  let q = @linalg.mat_make(n, n)
  for j = 0; j < n; j = j + 1 {
    let h = fp_powmod(fp_monom(p, j, 1), p, f)
    for i = 0; i < n; i = i + 1 {
      q[i][j] = h.coeff(i)
    }
  }
  for i = 0; i < n; i = i + 1 {
    q[i][i] = @kernel.submod(q[i][i], 1, p)
  }
  @linalg.nullspace(q, p)
}

///| Split a monic squarefree polynomial into monic irreducibles.
pub fn fp_split(f : FPoly) -> Array[FPoly] {
  let p = f.p
  let n = f.degree()
  if n <= 1 {
    return [f.monic()]
  }
  if fp_is_irreducible(f) {
    return [f.monic()]
  }
  let ker = fp_berlekamp_kernel(f)
  // Berlekamp's theorem: the kernel has dimension r >= 2 for reducible f.
  if ker.length() < 1 {
    abort("fp_split: Berlekamp kernel is empty")
  }
  for k = 0; k < ker.length(); k = k + 1 {
    let h = fpoly(p, ker[k])
    if h.degree() <= 0 {
      continue
    }
    for a = 0; a < p; a = a + 1 {
      let g = f.gcd(h.sub(fpoly(p, [a]))).monic()
      if g.degree() > 0 && g.degree() < n {
        let left = fp_split(g)
        let right = fp_split(f.div(g).monic())
        return left.add(right)
      }
    }
  }
  abort("fp_split: no proper factor found for " + f.to_string())
}

///| Equal-degree factorisation: the irreducible factors of `f`, all of which
///| have degree `d` by hypothesis (checked).
pub fn fp_edf(f : FPoly, d : Int) -> Array[FPoly] {
  let parts = fp_split(f.monic())
  for i = 0; i < parts.length(); i = i + 1 {
    if parts[i].degree() != d {
      abort("fp_edf: factor of unexpected degree")
    }
  }
  parts
}

///| p-th root in characteristic p: if `f' == 0` then `f = g^p` with
///| `g_i = f_{i*p}`.
fn fp_pth_root(f : FPoly) -> FPoly {
  let p = f.p
  let n = f.degree() / p
  let c : Array[Int] = Array::make(n + 1, 0)
  for i = 0; i <= n; i = i + 1 {
    c[i] = f.coeff(i * p)
  }
  fpoly(p, c)
}

///| **A26 full factorisation into monic irreducibles with multiplicities.**
///|
///| Characteristic p breaks Yun's algorithm (the derivative of `h^p` vanishes),
///| so this routine combines three facts instead:
///|   1. `f' = 0`  <=>  `f = g^p` with `g_i = f_{i p}`  (freshman's dream);
///|   2. `f / gcd(f, f')` is the product of those irreducible factors whose
///|      multiplicity is *not* divisible by p;
///|   3. a squarefree polynomial splits completely by ddf + Berlekamp.
///| Multiplicities are then measured by repeated division, which is exact.
pub fn fp_factor(f : FPoly) -> Array[(FPoly, Int)] {
  let out : Array[(FPoly, Int)] = []
  fp_factor_acc(f.monic(), 1, out)
  sort_factors(out)
  out
}

fn fp_factor_acc(g : FPoly, mult : Int, out : Array[(FPoly, Int)]) -> Unit {
  if g.degree() <= 0 {
    return
  }
  let p = g.p
  if g.deriv().is_zero() {
    fp_factor_acc(fp_pth_root(g), mult * p, out)
    return
  }
  let d = g.gcd(g.deriv())
  if d.degree() <= 0 {
    let parts = fp_split(g)
    for i = 0; i < parts.length(); i = i + 1 {
      out.push((parts[i], mult))
    }
    return
  }
  let s = g.div(d).monic()
  let mut rest = g
  let parts = fp_split(s)
  for i = 0; i < parts.length(); i = i + 1 {
    let h = parts[i]
    let mut m = 0
    while rest.rem(h).is_zero() {
      rest = rest.div(h).monic()
      m = m + 1
    }
    if m > 0 {
      out.push((h, m * mult))
    }
  }
  if rest.degree() > 0 {
    fp_factor_acc(rest, mult, out)
  }
}

fn sort_factors(out : Array[(FPoly, Int)]) -> Unit {
  out.sort_by(
    fn(u, v) {
      if u.0.degree() != v.0.degree() {
        u.0.degree() - v.0.degree()
      } else {
        compare_coeffs(u.0, v.0)
      }
    },
  )
}

fn compare_coeffs(a : FPoly, b : FPoly) -> Int {
  let n = Int::max(a.degree(), b.degree())
  let mut i = n
  while i >= 0 {
    if a.coeff(i) != b.coeff(i) {
      return a.coeff(i) - b.coeff(i)
    }
    i = i - 1
  }
  0
}

///| All roots of `f` in F_p (complete enumeration — exact because F_p is finite).
pub fn fp_roots(f : FPoly) -> Array[Int] {
  let out : Array[Int] = []
  for a = 0; a < f.p; a = a + 1 {
    if f.eval(a) == 0 {
      out.push(a)
    }
  }
  out
}

///| Number of roots in F_p, counted with multiplicity via the squarefree part.
pub fn fp_root_count(f : FPoly) -> Int {
  fp_roots(f).length()
}
