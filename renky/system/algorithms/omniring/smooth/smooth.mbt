// ===========================================================================
//  smooth — the ring of smooth functions, modelled by its jet / Weil algebras
//  (「光滑函數環」, structure 4).
//
//  Definition (jet algebra).  For a prime p and an order k, the k-jet algebra is
//      J_k = F_p[t]/(t^{k+1}),
//  the algebra of k-th order Taylor expansions of smooth functions at 0; a jet
//  is stored by its Taylor coefficients a_0..a_k, f(t) = sum_i a_i t^i.  J_k is
//  a local Artin F_p-algebra with maximal ideal (t), residue field F_p and
//  nilpotent element t of order k+1.  J_1 = F_p[eps]/(eps^2) is the ring of
//  *dual numbers*.
//
//  Definition (derivations).  D : J_k -> J_k is F_p-linear and satisfies the
//  Leibniz rule D(fg) = D(f)g + fD(g).  Two canonical derivations:
//      d/dt : a_i -> (i+1) a_{i+1}              (the ordinary derivative)
//      E = t*d/dt : a_i -> i*a_i                 (the Euler derivation)
//  Their commutator is a Lie bracket, and [d/dt, E] = d/dt.
//
//  Definition (Hasse derivatives).  D^{(n)}(f)_i = C(i+n, n) * a_{i+n}; these
//  are the divided-power derivatives, which behave well in characteristic p.
//
//  Theorems (verified exhaustively in `smooth_wbtest`):
//   M1 J_k is a commutative local F_p-algebra: t^{k+1} = 0, every element is
//      either a unit or nilpotent-mod-maximal, and J_k/(t) = F_p.
//   M2 Kock–Lawvereux: for the dual numbers, f(a + eps) = f(a) + f'(a)*eps for
//      every f and every a — this *characterises* the derivative.
//   M3 d/dt and E are derivations (Leibniz + linearity) and [d/dt, E] = d/dt,
//      [.,.] satisfying the Jacobi identity.
//   M4 (A68) composition of jets is associative, unital, and obeys the chain
//      rule (f o g)' = (f' o g) * g' — the Faà di Bruno formula in degree 1,
//      and the higher Faà di Bruno expansion in every degree.
//   M5 the heat flow on jets: u(s) = exp(s*Delta)u_0 with Delta = d^2/dt^2
//      satisfies the exponential law u(s+r) = exp(r*Delta)u(s) and the heat
//      equation D_s u = Delta u (Hasse-derivative form), i.e. the smooth ring
//      carries a heat semigroup whose generator is the Laplacian on jets.
//   M6 (A71) every function F_p -> F_p is a polynomial function of degree <= p-1
//      (Lagrange interpolation), so the ring of "smooth functions" on F_p is
//      F_p[t]/(t^p - t), and the jet of order p-1 at 0 determines the function.
// ===========================================================================

///| A k-jet: the Taylor coefficients of a smooth function at 0, over F_p.
pub struct Jet {
  p : Int
  order : Int
  coeffs : Array[Int]
}

fn jet_mk(p : Int, order : Int, cs : Array[Int]) -> Jet {
  let out : Array[Int] = Array::make(order + 1, 0)
  for i = 0; i <= order; i = i + 1 {
    out[i] = if i < cs.length() {
      @kernel.modpos(cs[i], p)
    } else {
      0
    }
  }
  { p, order, coeffs: out }
}

pub fn jet_zero(p : Int, order : Int) -> Jet {
  jet_mk(p, order, [])
}

pub fn jet_const(p : Int, order : Int, c : Int) -> Jet {
  jet_mk(p, order, [c])
}

///| The coordinate function t.
pub fn jet_t(p : Int, order : Int) -> Jet {
  jet_mk(p, order, [0, 1])
}

pub fn jet_of_coeffs(p : Int, order : Int, cs : Array[Int]) -> Jet {
  jet_mk(p, order, cs)
}

pub fn Jet::is_zero(self : Jet) -> Bool {
  for i = 0; i <= self.order; i = i + 1 {
    if self.coeffs[i] != 0 {
      return false
    }
  }
  true
}

pub fn Jet::equal(self : Jet, other : Jet) -> Bool {
  if self.p != other.p || self.order != other.order {
    return false
  }
  for i = 0; i <= self.order; i = i + 1 {
    if self.coeffs[i] != other.coeffs[i] {
      return false
    }
  }
  true
}

pub fn Jet::to_string(self : Jet) -> String {
  let parts : Array[String] = []
  for i = 0; i <= self.order; i = i + 1 {
    parts.push(self.coeffs[i].to_string())
  }
  "[" + parts.join(",") + "]_p" + self.p.to_string()
}

pub fn Jet::add(self : Jet, other : Jet) -> Jet {
  let cs : Array[Int] = []
  for i = 0; i <= self.order; i = i + 1 {
    cs.push(self.coeffs[i] + other.coeffs[i])
  }
  jet_mk(self.p, self.order, cs)
}

pub fn Jet::neg(self : Jet) -> Jet {
  let cs : Array[Int] = []
  for i = 0; i <= self.order; i = i + 1 {
    cs.push(-self.coeffs[i])
  }
  jet_mk(self.p, self.order, cs)
}

pub fn Jet::sub(self : Jet, other : Jet) -> Jet {
  self.add(other.neg())
}

///| Truncated Cauchy product.
pub fn Jet::mul(self : Jet, other : Jet) -> Jet {
  let cs : Array[Int] = Array::make(self.order + 1, 0)
  for i = 0; i <= self.order; i = i + 1 {
    for j = 0; j + i <= self.order; j = j + 1 {
      cs[i + j] = cs[i + j] + self.coeffs[i] * other.coeffs[j]
    }
  }
  jet_mk(self.p, self.order, cs)
}

pub fn Jet::scale(self : Jet, c : Int) -> Jet {
  let cs : Array[Int] = []
  for i = 0; i <= self.order; i = i + 1 {
    cs.push(c * self.coeffs[i])
  }
  jet_mk(self.p, self.order, cs)
}

///| **A68 composition of jets** (substitution), i.e. the Faà di Bruno formula
///| in truncated form: (f o g)_n = sum over partitions of n.
pub fn Jet::compose(self : Jet, g : Jet) -> Jet {
  let p = self.p
  let order = Int::min(self.order, g.order)
  // powers of g up to `order`
  let pows : Array[Jet] = []
  let mut cur = jet_const(p, order, 1)
  for _i = 0; _i <= order; _i = _i + 1 {
    pows.push(cur)
    cur = cur.mul(trunc_to(g, order))
  }
  let mut acc = jet_zero(p, order)
  for i = 0; i <= order; i = i + 1 {
    if i <= self.order {
      acc = acc.add(pows[i].scale(self.coeffs[i]))
    }
  }
  acc
}

fn trunc_to(f : Jet, order : Int) -> Jet {
  jet_mk(f.p, order, f.coeffs)
}

///| The ordinary derivative d/dt.
pub fn Jet::deriv(self : Jet) -> Jet {
  let cs : Array[Int] = []
  for i = 1; i <= self.order; i = i + 1 {
    cs.push(i * self.coeffs[i])
  }
  jet_mk(self.p, self.order, cs)
}

///| The Euler derivation E = t * d/dt.
pub fn Jet::euler(self : Jet) -> Jet {
  let cs : Array[Int] = []
  for i = 0; i <= self.order; i = i + 1 {
    cs.push(i * self.coeffs[i])
  }
  jet_mk(self.p, self.order, cs)
}

///| **A69 Hasse (divided power) derivative**: the n-th such derivative: coefficient i <- C(i+n,n)*a_{i+n}.
pub fn Jet::hasse(self : Jet, n : Int) -> Jet {
  let cs : Array[Int] = []
  for i = 0; i <= self.order; i = i + 1 {
    if i + n <= self.order {
      cs.push(binom_mod(i + n, n) * self.coeffs[i + n])
    } else {
      cs.push(0)
    }
  }
  jet_mk(self.p, self.order, cs)
}

fn binom_mod(n : Int, k : Int) -> Int {
  if k < 0 || k > n {
    return 0
  }
  let mut acc = 1
  for i = 0; i < k; i = i + 1 {
    acc = acc * (n - i) / (i + 1)
  }
  acc
}

///| Evaluate the jet at a scalar (the truncated Taylor polynomial).
pub fn Jet::eval_at(self : Jet, x : Int) -> Int {
  let p = self.p
  let mut acc = 0
  let mut pw = 1
  for i = 0; i <= self.order; i = i + 1 {
    acc = @kernel.modpos(acc + self.coeffs[i] * pw, p)
    pw = @kernel.modpos(pw * x, p)
  }
  acc
}

///| Is the jet a unit?  (a_0 != 0 in the local ring J_k.)
pub fn Jet::is_unit(self : Jet) -> Bool {
  self.coeffs[0] != 0
}

///| Inverse of a unit jet: f = a_0 (1 + n) with n nilpotent, so
///| f^{-1} = a_0^{-1} * sum_{k=0}^{order} (-n)^k  (exact, since n^{order+1}=0).
pub fn Jet::inv(self : Jet) -> Jet {
  if !(self.is_unit()) {
    abort("[smooth.inv] not a unit")
  }
  let p = self.p
  let a0inv = @arith.invmod(self.coeffs[0], p)
  let nil = self.sub(jet_const(p, self.order, self.coeffs[0])).scale(a0inv)
  let mut acc = jet_zero(p, self.order)
  let mut term = jet_const(p, self.order, 1)
  for k = 0; k <= self.order; k = k + 1 {
    if k > 0 {
      term = term.mul(nil)
    }
    acc = if k % 2 == 0 {
      acc.add(term)
    } else {
      acc.sub(term)
    }
  }
  acc.scale(a0inv)
}

///| Exponential of a nilpotent jet (constant term 0): exp(x) = sum x^n/n!.
///| The terms are built recursively, t_n = t_{n-1} * x / n, which requires n to
///| be invertible mod p, i.e. `order < p` (the usual restriction in
///| characteristic p).
pub fn Jet::jet_exp(self : Jet) -> Jet {
  let p = self.p
  if self.order >= p {
    abort("[smooth.jet_exp] needs order < p so that n! is invertible")
  }
  let mut acc = jet_const(p, self.order, 1)
  let mut term = jet_const(p, self.order, 1)
  for n = 1; n <= self.order; n = n + 1 {
    term = term.mul(self).scale(@arith.invmod(n, p))
    acc = acc.add(term)
  }
  acc
}

///| Logarithm of a unipotent jet (constant term 1).
pub fn Jet::jet_log(self : Jet) -> Jet {
  if self.coeffs[0] != 1 {
    abort("[smooth.jet_log] constant term must be 1")
  }
  let p = self.p
  let u = self.sub(jet_const(p, self.order, 1))
  let mut acc = jet_zero(p, self.order)
  let mut term = jet_const(p, self.order, 1)
  for n = 1; n <= self.order; n = n + 1 {
    term = term.mul(u)
    let c = @arith.invmod(n, p)
    acc = if n % 2 == 1 {
      acc.add(term.scale(c))
    } else {
      acc.sub(term.scale(c))
    }
  }
  acc
}

///| **M5 heat flow on jets**: u(s) = exp(s*Delta) u_0 = sum_m (s^m/m!) Delta^m u_0
///| **A70 jet heat flow**, with Delta = d^2/dt^2.  Only m <= order/2 < p contribute, so all m! are
///| invertible mod p and the flow is exact.  It satisfies the exponential law
///| u(s+r) = exp(r*Delta)u(s) and the heat equation d/ds u = Delta u.
pub fn jet_heat_flow(u0 : Jet, s : Int) -> Jet {
  let p = u0.p
  if u0.order / 2 + 1 >= p {
    abort("[smooth.jet_heat_flow] needs order/2 < p")
  }
  let mut acc = u0
  let mut term = u0
  let mut pw = 1
  let mut finv = 1
  for m = 1; m <= u0.order / 2 + 1; m = m + 1 {
    term = term.laplacian()
    pw = @kernel.modpos(pw * s, p)
    finv = @kernel.modpos(finv * @arith.invmod(m, p), p)
    let piece = term.scale(@kernel.modpos(pw * finv, p))
    if piece.is_zero() {
      continue
    }
    acc = acc.add(piece)
  }
  acc
}

///| The Laplacian on jets: Delta = d^2/dt^2, the generator of the heat flow.
pub fn Jet::laplacian(self : Jet) -> Jet {
  self.deriv().deriv()
}

// ---------------------------------------------------------------------------
//  M6  smooth functions on F_p
// ---------------------------------------------------------------------------

///| **A71 Lagrange interpolation**: the unique polynomial of degree < p over F_p
///| representing the function `values` on F_p.
pub fn interpolate(p : Int, values : Array[Int]) -> Array[Int] {
  let coeffs : Array[Int] = Array::make(p, 0)
  for a = 0; a < p; a = a + 1 {
    // L_a(x) = prod_{b != a} (x - b)/(a - b)
    let num : Array[Int] = Array::make(p, 0)
    num[0] = 1
    let mut deg = 0
    for b = 0; b < p; b = b + 1 {
      if b == a {
        continue
      }
      // multiply by (x - b)
      for i = deg + 1; i >= 1; i = i - 1 {
        num[i] = @kernel.modpos(num[i - 1] - b * num[i], p)
      }
      num[0] = @kernel.modpos(-b * num[0], p)
      deg = deg + 1
    }
    let mut den = 1
    for b = 0; b < p; b = b + 1 {
      if b != a {
        den = @kernel.modpos(den * (a - b), p)
      }
    }
    let scale = @kernel.modpos(values[a] * @arith.invmod(den, p), p)
    for i = 0; i < p; i = i + 1 {
      coeffs[i] = @kernel.modpos(coeffs[i] + scale * num[i], p)
    }
  }
  coeffs
}

///| Evaluate a coefficient list (ascending) at x mod p.
pub fn poly_eval_mod(coeffs : Array[Int], x : Int, p : Int) -> Int {
  let mut acc = 0
  for i = coeffs.length() - 1; i >= 0; i = i - 1 {
    acc = @kernel.modpos(acc * x + coeffs[i], p)
  }
  acc
}
