// ===========================================================================
//  heat — the Maslov dequantization family (「轉熱」).
//
//  Definition (heat semifield).  For h > 0 set
//      a (+)_h b  =  -h * log( exp(-a/h) + exp(-b/h) )      (soft minimum)
//      a (*)_h b  =  a + b
//  Then (R, (+)_h, (*)_h) is a semifield, and
//      phi_h : R -> R_{>0},   phi_h(x) = exp(-x/h)
//  is an isomorphism of semifields onto (R_{>0}, +, *):
//      phi_h(a (+)_h b) = phi_h(a) + phi_h(b),      phi_h(a (*)_h b) = phi_h(a) * phi_h(b).
//
//  Theorem (Maslov dequantization / zero temperature limit).  As h -> 0+
//      a (+)_h b  ->  min(a, b),        a (*)_h b  =  a + b,
//  and the convergence is uniform with the explicit bound
//      -h * log 2 <= a (+)_h b - min(a,b) <= 0,
//  equality on the left exactly when a = b (so (+)_h is idempotent only in the
//  limit, as befits a semifield degenerating to an idempotent semiring).
//  Hence the tropical (min-plus) semiring is the h -> 0 limit of the heat
//  semifields: this is the morphism "turn to heat" of the axiom.
//
//  Theorem (free energy / entropy).  For energies E_0..E_{n-1},
//      F(h) = -h log sum_i exp(-E_i/h)      (free energy)
//      <E>(h) = sum_i E_i p_i,  p_i = exp(-E_i/h)/Z      (Gibbs state)
//      S(h) = -sum_i p_i log p_i                         (entropy)
//  satisfy the Legendre relation F = <E> - h*S, F is increasing in h, and
//      lim_{h->0+} F(h) = lim_{h->0+} <E>(h) = min_i E_i,   lim_{h->0+} S(h) = 0,
//  i.e. the zero-temperature limit of statistical mechanics is the tropical
//  minimum, and the ground state is the tropical argmin.
//
//  No third-party dependencies: exp and log are implemented here from their
//  defining series with range reduction.
// ===========================================================================

///| ln 2, to double precision (used for range reduction).
pub let ln2 : Double = 0.6931471805599453

///| **A44 exp** — Taylor series with range reduction x = k*ln2 + r, |r| <= ln2/2.
pub fn expd(x : Double) -> Double {
  if x > 700.0 {
    return 1.0e308
  }
  if x < -700.0 {
    return 0.0
  }
  let k = Double::round(x / ln2).to_int()
  let r = x - Double::from_int(k) * ln2
  // sum_{n>=0} r^n / n!
  let mut term = 1.0
  let mut acc = 1.0
  for n = 1; n <= 40; n = n + 1 {
    term = term * r / Double::from_int(n)
    acc = acc + term
    if Double::abs(term) < 1.0e-18 * Double::abs(acc) {
      break
    }
  }
  acc * pow2(k)
}

///| 2^k by binary exponentiation (self-contained: no `core/math`).
fn pow2(k : Int) -> Double {
  let mut n = if k < 0 {
    -k
  } else {
    k
  }
  let mut acc = 1.0
  let mut base = 2.0
  while n > 0 {
    if n % 2 == 1 {
      acc = acc * base
    }
    base = base * base
    n = n / 2
  }
  if k < 0 {
    1.0 / acc
  } else {
    acc
  }
}

///| **A45 log** — log x = e*ln2 + 2*atanh((m-1)/(m+1)) with x = m*2^e, m in [1/2,2).
pub fn logd(x : Double) -> Double {
  if x <= 0.0 {
    abort("[heat.logd] log of a non-positive number")
  }
  let mut e = 0
  let mut m = x
  while m >= 2.0 {
    m = m / 2.0
    e = e + 1
  }
  while m < 0.5 {
    m = m * 2.0
    e = e - 1
  }
  let s = (m - 1.0) / (m + 1.0)
  let s2 = s * s
  let mut term = s
  let mut acc = s
  for n = 1; n <= 60; n = n + 1 {
    term = term * s2
    acc = acc + term / Double::from_int(2 * n + 1)
    if Double::abs(term) < 1.0e-18 {
      break
    }
  }
  Double::from_int(e) * ln2 + 2.0 * acc
}

///| **A46 Maslov soft minimum** (heat addition at inverse temperature 1/h).
///| `softmin_h(a, b, h) = -h log(exp(-a/h) + exp(-b/h))`.
pub fn softmin_h(a : Double, b : Double, h : Double) -> Double {
  if h <= 0.0 {
    return Double::min(a, b)
  }
  let m = Double::min(a, b)
  // shift by the minimum for numerical stability (exact identity:
  // softmin(a,b) = m - h log(1 + exp(-|a-b|/h)))
  m - h * logd(1.0 + expd(-Double::abs(a - b) / h))
}

///| Soft maximum (heat addition in the max-plus convention).
pub fn softmax_h(a : Double, b : Double, h : Double) -> Double {
  if h <= 0.0 {
    return Double::max(a, b)
  }
  -(softmin_h(-a, -b, h))
}

///| n-ary soft minimum: -h log sum_i exp(-a_i/h).
pub fn softmin_all(xs : Array[Double], h : Double) -> Double {
  if xs.is_empty() {
    abort("[heat.softmin_all] empty family")
  }
  if h <= 0.0 {
    let mut m = xs[0]
    for i = 1; i < xs.length(); i = i + 1 {
      m = Double::min(m, xs[i])
    }
    return m
  }
  let mut m = xs[0]
  for i = 1; i < xs.length(); i = i + 1 {
    m = Double::min(m, xs[i])
  }
  let mut acc = 0.0
  for i = 0; i < xs.length(); i = i + 1 {
    acc = acc + expd(-(xs[i] - m) / h)
  }
  m - h * logd(acc)
}

///| Heat addition of the Maslov family, as a binary operation on R.
pub fn heat_add(a : Double, b : Double, h : Double) -> Double {
  softmin_h(a, b, h)
}

///| Heat multiplication of the Maslov family: ordinary addition.
pub fn heat_mul(a : Double, b : Double) -> Double {
  a + b
}

///| The dequantization isomorphism phi_h(x) = exp(-x/h).
pub fn phi(x : Double, h : Double) -> Double {
  if h <= 0.0 {
    abort("[heat.phi] h must be positive")
  }
  expd(-x / h)
}

///| Its inverse phi_h^{-1}(y) = -h log y.
pub fn phi_inv(y : Double, h : Double) -> Double {
  if h <= 0.0 {
    abort("[heat.phi_inv] h must be positive")
  }
  -(h * logd(y))
}

///| **A47 Gibbs state**: partition function Z(h) = sum_i exp(-E_i/h).
pub fn partition(energies : Array[Double], h : Double) -> Double {
  let mut acc = 0.0
  for i = 0; i < energies.length(); i = i + 1 {
    acc = acc + expd(-energies[i] / h)
  }
  acc
}

///| **Free energy** F(h) = -h log Z(h); identical to the n-ary soft minimum.
pub fn free_energy(energies : Array[Double], h : Double) -> Double {
  softmin_all(energies, h)
}

///| Gibbs probabilities p_i = exp(-E_i/h) / Z(h).
pub fn gibbs(energies : Array[Double], h : Double) -> Array[Double] {
  let z = partition(energies, h)
  let out : Array[Double] = []
  for i = 0; i < energies.length(); i = i + 1 {
    out.push(expd(-energies[i] / h) / z)
  }
  out
}

///| Mean energy <E>(h).
pub fn mean_energy(energies : Array[Double], h : Double) -> Double {
  let ps = gibbs(energies, h)
  let mut acc = 0.0
  for i = 0; i < energies.length(); i = i + 1 {
    acc = acc + energies[i] * ps[i]
  }
  acc
}

///| Gibbs entropy S(h) = -sum_i p_i log p_i.
pub fn entropy(energies : Array[Double], h : Double) -> Double {
  let ps = gibbs(energies, h)
  let mut acc = 0.0
  for i = 0; i < ps.length(); i = i + 1 {
    if ps[i] > 1.0e-300 {
      acc = acc - ps[i] * logd(ps[i])
    }
  }
  acc
}

///| Legendre transform of the entropy: F(beta) = min_E (E + S(E)/beta) with
///| beta = 1/h, evaluated through the Gibbs state.
pub fn legendre_check(energies : Array[Double], h : Double) -> (Double, Double) {
  let e = mean_energy(energies, h)
  let s = entropy(energies, h)
  (e - h * s, softmin_all(energies, h))
}

///| Zero-temperature (tropical) limit of a family of energies: the minimum,
///| together with the ground-state index (tropical argmin).
pub fn to_tropical(energies : Array[Double]) -> (Double, Int) {
  let mut m = energies[0]
  let mut k = 0
  for i = 1; i < energies.length(); i = i + 1 {
    if energies[i] < m {
      m = energies[i]
      k = i
    }
  }
  (m, k)
}
