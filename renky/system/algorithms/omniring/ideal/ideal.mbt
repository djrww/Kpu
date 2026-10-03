using @poly {
  type MPoly,
  mp_make,
  mp_one,
  mp_var,
  mp_monom,
  mp_normal_form,
  mp_divmod,
  mp_spoly,
  mp_reduce_basis,
  mp_standard_monomials,
  cmp_monom,
  total_deg,
  monom_divides,
}

// ===========================================================================
//  即時理想狀態 — the *live ideal state* engine.
//
//  An `IdealState` is a mutable record of an ideal I of R = F_p[x_0..x_{n-1}]
//  that keeps itself canonical at all times:
//
//      gens   : every generator ever fed in, in arrival order (the history)
//      basis  : the reduced Gröbner basis of the ideal they generate *now*
//      finite : whether R/I is finite dimensional over F_p
//      dim    : dim_{F_p} R/I  (-1 when infinite)
//      kdim   : Krull dimension of R/I
//
//  `add` performs an *incremental* Gröbner update (only the new pairs are
//  considered), so the state is always consistent — this is the "real-time
//  ideal state" of the axiom: at every moment the observer can read off
//  membership, quotient basis, Hilbert function and Krull dimension without
//  recomputing from scratch.
//
//  Theorems verified exhaustively in `ideal_wbtest`:
//    T1  f in I  <=>  nf_G(f) = 0                        (Gröbner characterisation)
//    T2  G is a Gröbner basis <=> every S-poly reduces to 0   (Buchberger)
//    T3  the standard monomials form an F_p-basis of R/I  (Macaulay)
//    T4  dim R/I = sum_d H(d)  and  H(d) = 0 for d >> 0 when kdim = 0
//    T5  kdim(R/I) = max |S| over variable sets S with supp(LM g) not contained in S
//    T6  I(V(I)) = rad(I + (x_i^p - x_i))                 (finite-field Nullstellensatz)
//    T7  f in rad(I) <=> 1 in (I, 1 - y f)                (Rabinowitsch)
//    T8  I maximal <=> R/I is a field <=> R/I has no zero divisors
// ===========================================================================

pub struct IdealState {
  p : Int
  nvars : Int
  gens : Array[MPoly]
  raw : Array[MPoly]
  mut basis : Array[MPoly]
  mut updates : Int
  mut pairs_done : Int
  mut finite : Bool
  mut dim : Int
  mut kdim : Int
  mut var_bounds : Array[Int]
  mut is_unit : Bool
}

pub fn ideal_new(p : Int, nvars : Int) -> IdealState {
  let st : IdealState = {
    p,
    nvars,
    gens: [],
    raw: [],
    basis: [],
    updates: 0,
    pairs_done: 0,
    finite: false,
    dim: -1,
    kdim: nvars,
    var_bounds: Array::make(nvars, 0),
    is_unit: false,
  }
  st.refresh()
  st
}

///| The zero ideal is proper; `add` of a unit collapses everything.
pub fn IdealState::add(self : IdealState, g : MPoly) -> Unit {
  self.gens.push(g)
  self.updates = self.updates + 1
  // incremental Buchberger: extend the current basis with the new element and
  // only process pairs that involve it (older pairs already reduced to zero).
  if g.is_zero() {
    return
  }
  let p = self.p
  let g1 = g.scale(@arith.invmod(g.lc(), p))
  // `raw` only ever grows, so the pair bookkeeping stays valid across updates;
  // the published `basis` is the reduction of `raw`, which is again a Gröbner
  // basis (reduction preserves the leading-monomial ideal).
  let old = self.raw.length()
  self.raw.push(g1)
  let mut changed = true
  let mut rounds = 0
  while changed && rounds < 200 {
    changed = false
    rounds = rounds + 1
    let n = self.raw.length()
    let pairs : Array[(Int, Int)] = []
    for i = 0; i < n; i = i + 1 {
      for j = i + 1; j < n; j = j + 1 {
        if i >= old || j >= old {
          pairs.push((i, j))
        }
      }
    }
    for k = 0; k < pairs.length(); k = k + 1 {
      self.pairs_done = self.pairs_done + 1
      let s = mp_spoly(self.raw[pairs[k].0], self.raw[pairs[k].1])
      let r = mp_normal_form(s, self.raw)
      if !(r.is_zero()) {
        self.raw.push(r.scale(@arith.invmod(r.lc(), p)))
        changed = true
      }
    }
  }
  self.basis = mp_reduce_basis(self.raw)
  self.refresh()
}

///| Recompute the derived invariants *exactly*.
///|
///| Zero-dimensionality criterion (exact, not a degree heuristic): R/I is a
///| finite dimensional F_p-vector space iff for every variable x_v some power
///| x_v^k is a leading monomial of the Gröbner basis, i.e. iff
///| I ∩ F_p[x_v] != 0 for all v.  `var_bounds[v]` records the least such k, and
///| then every standard monomial has exponent vector inside the box
///| `0 <= e_v < var_bounds[v]`, which is enumerated exhaustively.
fn IdealState::refresh(self : IdealState) -> Unit {
  // the unit ideal: a nonzero constant in the basis collapses R/I to zero
  for i = 0; i < self.basis.length(); i = i + 1 {
    let f = self.basis[i]
    if f.is_zero() {
      continue
    }
    let lm = f.lm()
    let mut allzero = true
    for v = 0; v < self.nvars; v = v + 1 {
      if lm[v] != 0 {
        allzero = false
      }
    }
    if allzero {
      self.is_unit = true
      self.finite = true
      self.dim = 0
      self.kdim = -1
      self.var_bounds = Array::make(self.nvars, 1)
      return
    }
  }
  self.is_unit = false
  let bounds : Array[Int] = Array::make(self.nvars, 0)
  for i = 0; i < self.basis.length(); i = i + 1 {
    if self.basis[i].is_zero() {
      continue
    }
    let lm = self.basis[i].lm()
    let mut v = -1
    let mut pure = true
    for j = 0; j < self.nvars; j = j + 1 {
      if lm[j] > 0 {
        if v >= 0 {
          pure = false
        }
        v = j
      }
    }
    if pure && v >= 0 {
      if bounds[v] == 0 || lm[v] < bounds[v] {
        bounds[v] = lm[v]
      }
    }
  }
  self.var_bounds = bounds
  let mut finite = true
  for v = 0; v < self.nvars; v = v + 1 {
    if bounds[v] == 0 {
      finite = false
    }
  }
  self.finite = finite
  self.kdim = mp_krull_dim(self.basis, self.nvars)
  if finite {
    self.dim = boxed_standard_monomials(self.basis, bounds).length()
  } else {
    self.dim = -1
  }
}

///| Standard monomials inside the box `0 <= e_v < bounds[v]` (exact basis of
///| R/I when the quotient is finite dimensional).
pub fn boxed_standard_monomials(
  basis : Array[MPoly],
  bounds : Array[Int],
) -> Array[Array[Int]] {
  let nvars = bounds.length()
  let out : Array[Array[Int]] = []
  let mut total = 1
  for v = 0; v < nvars; v = v + 1 {
    total = total * Int::max(1, bounds[v])
  }
  for k = 0; k < total; k = k + 1 {
    let e : Array[Int] = Array::make(nvars, 0)
    let mut rest = k
    let mut blocked = false
    for v = 0; v < nvars; v = v + 1 {
      let b = Int::max(1, bounds[v])
      e[v] = rest % b
      rest = rest / b
    }
    for i = 0; i < basis.length(); i = i + 1 {
      if !(basis[i].is_zero()) && monom_divides(basis[i].lm(), e) {
        blocked = true
      }
    }
    if !(blocked) {
      out.push(e)
    }
  }
  out.sort_by(cmp_monom)
  out
}

///| **A36 membership test and certificate**: `f in I` iff its normal form w.r.t.
///| the Gröbner basis is 0.
pub fn IdealState::in_ideal(self : IdealState, f : MPoly) -> Bool {
  mp_normal_form(f, self.basis).is_zero()
}

pub fn IdealState::normal_form(self : IdealState, f : MPoly) -> MPoly {
  mp_normal_form(f, self.basis)
}

///| Certificate of membership: `f = sum_i q_i g_i`.
pub fn IdealState::membership_witness(
  self : IdealState,
  f : MPoly,
) -> (Array[MPoly], MPoly) {
  mp_divmod(f, self.basis)
}

///| **A35 Krull dimension** of R/I from the initial ideal (T5).
pub fn mp_krull_dim(basis : Array[MPoly], nvars : Int) -> Int {
  let mut best = 0
  let total = 1 << nvars
  for mask = 0; mask < total; mask = mask + 1 {
    let mut ok = true
    for i = 0; i < basis.length(); i = i + 1 {
      if basis[i].is_zero() {
        continue
      }
      let lm = basis[i].lm()
      let mut supported = true
      let mut any = false
      for j = 0; j < nvars; j = j + 1 {
        if lm[j] > 0 {
          any = true
          if ((mask >> j) & 1) == 0 {
            supported = false
          }
        }
      }
      if any && supported {
        ok = false
      }
    }
    if ok {
      let cnt = popcount(mask)
      if cnt > best {
        best = cnt
      }
    }
  }
  best
}

fn popcount(x : Int) -> Int {
  let mut c = 0
  let mut y = x
  while y != 0 {
    c = c + (y & 1)
    y = y >> 1
  }
  c
}

///| F_p-basis of R/I as exponent vectors (exact; empty when R/I is infinite).
pub fn IdealState::quotient_basis(self : IdealState) -> Array[Array[Int]] {
  if !(self.finite) || self.is_unit {
    return []
  }
  boxed_standard_monomials(self.basis, self.var_bounds)
}

///| **A37 Hilbert function** H(d) = dim_{F_p} (R/I)_d for d = 0..maxdeg.
///| For a finite quotient this is exact for every d, since all standard
///| monomials have been enumerated.
pub fn IdealState::hilbert_fn(self : IdealState, maxdeg : Int) -> Array[Int] {
  let out : Array[Int] = Array::make(maxdeg + 1, 0)
  if self.is_unit {
    return out
  }
  let std = if self.finite {
    boxed_standard_monomials(self.basis, self.var_bounds)
  } else {
    mp_standard_monomials(self.basis, self.nvars, maxdeg).0
  }
  for i = 0; i < std.length(); i = i + 1 {
    let d = total_deg(std[i])
    if d <= maxdeg {
      out[d] = out[d] + 1
    }
  }
  out
}

///| **A38 variety**: all F_p-rational points of V(I) (complete enumeration:
///| F_p^n is finite).
pub fn IdealState::variety(self : IdealState) -> Array[Array[Int]] {
  let out : Array[Array[Int]] = []
  let total = @arith.pow_int(self.p, self.nvars)
  for k = 0; k < total; k = k + 1 {
    let pt = @arith.digits(k, self.p)
    let x : Array[Int] = Array::make(self.nvars, 0)
    for i = 0; i < self.nvars; i = i + 1 {
      x[i] = if i < pt.length() {
        pt[i]
      } else {
        0
      }
    }
    let mut vanishes = true
    for i = 0; i < self.gens.length(); i = i + 1 {
      if self.gens[i].eval(x) != 0 {
        vanishes = false
      }
    }
    if vanishes {
      out.push(x)
    }
  }
  out
}

///| The ideal `(x_0^p - x_0, ..., x_{n-1}^p - x_{n-1})`: cutting down to F_p-points.
pub fn field_equations(p : Int, nvars : Int) -> Array[MPoly] {
  let out : Array[MPoly] = []
  for i = 0; i < nvars; i = i + 1 {
    let e = Array::make(nvars, 0)
    e[i] = p
    out.push(mp_monom(p, nvars, e, 1).sub(mp_var(p, nvars, i)))
  }
  out
}

///| `I + (x_i^p - x_i)`, whose variety is exactly the set of F_p-points of I.
pub fn IdealState::with_field_equations(self : IdealState) -> IdealState {
  let st = ideal_new(self.p, self.nvars)
  for i = 0; i < self.gens.length(); i = i + 1 {
    st.add(self.gens[i])
  }
  let fe = field_equations(self.p, self.nvars)
  for i = 0; i < fe.length(); i = i + 1 {
    st.add(fe[i])
  }
  st
}

///| **A39 Seidenberg radical** of a zero-dimensional ideal: adjoin the squarefree part of the
///| minimal polynomial of each variable (Seidenberg's lemma).
pub fn IdealState::radical(self : IdealState) -> IdealState {
  let st = ideal_new(self.p, self.nvars)
  for i = 0; i < self.gens.length(); i = i + 1 {
    st.add(self.gens[i])
  }
  for v = 0; v < self.nvars; v = v + 1 {
    match self.univariate_minpoly(v) {
      Some(f) => {
        let sf = fp_squarefree(f)
        if sf.degree() < f.degree() {
          st.add(univ_to_mpoly(sf, self.nvars, v))
        }
      }
      None => ()
    }
  }
  st
}

///| **A40 univariate minimal polynomial**: the monic generator of `I ∩ F_p[x_v]`, if it exists.
pub fn IdealState::univariate_minpoly(self : IdealState, v : Int) -> @ff.FPoly? {
  let mut best : @ff.FPoly? = None
  for i = 0; i < self.basis.length(); i = i + 1 {
    let lm = self.basis[i].lm()
    let mut pure = true
    for j = 0; j < self.nvars; j = j + 1 {
      if j != v && lm[j] != 0 {
        pure = false
      }
    }
    if pure && lm[v] > 0 {
      match best {
        None => best = Some(mpoly_to_univ(self.basis[i], v))
        Some(b) =>
          if self.basis[i].degree() < b.degree() {
            best = Some(mpoly_to_univ(self.basis[i], v))
          }
      }
    }
  }
  best
}

fn mpoly_to_univ(f : MPoly, v : Int) -> @ff.FPoly {
  let cs : Array[Int] = Array::make(f.degree() + 1, 0)
  for i = 0; i < f.terms.length(); i = i + 1 {
    let (e, c) = f.terms[i]
    cs[e[v]] = c
  }
  @ff.fpoly(f.p, cs)
}

fn univ_to_mpoly(f : @ff.FPoly, nvars : Int, v : Int) -> MPoly {
  let terms : Array[(Array[Int], Int)] = []
  for i = 0; i <= f.degree(); i = i + 1 {
    if f.coeff(i) != 0 {
      let e = Array::make(nvars, 0)
      e[v] = i
      terms.push((e, f.coeff(i)))
    }
  }
  mp_make(f.p, nvars, terms)
}

///| **A40 (cont.) squarefree part** of a univariate polynomial over F_p.
pub fn fp_squarefree(f : @ff.FPoly) -> @ff.FPoly {
  if f.degree() <= 0 {
    return f
  }
  let d = f.deriv()
  if d.is_zero() {
    // char p: f = g^p
    return fp_squarefree(fp_pth_root_pub(f))
  }
  f.div(f.gcd(d)).monic()
}

fn fp_pth_root_pub(f : @ff.FPoly) -> @ff.FPoly {
  let p = f.p
  let n = f.degree() / p
  let c : Array[Int] = Array::make(n + 1, 0)
  for i = 0; i <= n; i = i + 1 {
    c[i] = f.coeff(i * p)
  }
  @ff.fpoly(p, c)
}

///| **A41 maximality test**.  Is I maximal?  Equivalent (T8) to R/I being a field, which for a finite
///| quotient is equivalent to having no zero divisors.
pub fn IdealState::is_maximal(self : IdealState) -> Bool {
  if !(self.finite) || self.dim <= 0 {
    return false
  }
  let bas = self.quotient_basis()
  if bas.length() != self.dim {
    return false
  }
  // A = R/I is an Artinian F_p-algebra, so every element is either a unit or a
  // zero divisor; hence A is a field iff no nonzero element is a zero divisor
  // iff for every nonzero a the F_p-linear map m_a : A -> A is injective
  // iff det(m_a) != 0 in F_p.  All p^dim elements are enumerated.
  let p = self.p
  let total = @arith.pow_int(p, self.dim)
  for c = 1; c < total; c = c + 1 {
    let ds = @arith.digits(c, p)
    let terms : Array[(Array[Int], Int)] = []
    for i = 0; i < bas.length(); i = i + 1 {
      terms.push((bas[i], if i < ds.length() {
        ds[i]
      } else {
        0
      }))
    }
    let a = mp_make(p, self.nvars, terms)
    if @linalg.det_mod(self.mul_matrix(a), p) == 0 {
      return false
    }
  }
  true
}

///| **A42 multiplication matrix** of `f` on the F_p-basis of R/I.
pub fn IdealState::mul_matrix(self : IdealState, f : MPoly) -> @linalg.Mat {
  let bas = self.quotient_basis()
  let n = bas.length()
  let m = @linalg.mat_make(n, n)
  for j = 0; j < n; j = j + 1 {
    let bj = mp_monom(self.p, self.nvars, bas[j], 1)
    let prod = self.normal_form(f.mul(bj))
    for i = 0; i < prod.terms.length(); i = i + 1 {
      let (e, c) = prod.terms[i]
      for t = 0; t < n; t = t + 1 {
        if cmp_monom(bas[t], e) == 0 {
          m[t][j] = c
        }
      }
    }
  }
  m
}

///| **A43 Rabinowitsch trick** (T7): `f in rad(I)` iff `1 in (I, 1 - y*f)` in one more
///| variable.  Returns true when f lies in the radical.
pub fn IdealState::in_radical(self : IdealState, f : MPoly) -> Bool {
  let p = self.p
  let n = self.nvars
  let st = ideal_new(p, n + 1)
  for i = 0; i < self.gens.length(); i = i + 1 {
    st.add(lift_vars(self.gens[i], n + 1))
  }
  let y = mp_var(p, n + 1, n)
  st.add(mp_one(p, n + 1).sub(y.mul(lift_vars(f, n + 1))))
  st.in_ideal(mp_one(p, n + 1))
}

fn lift_vars(f : MPoly, nvars : Int) -> MPoly {
  mp_make(
    f.p,
    nvars,
    f.terms.map(
      fn(t) {
        let e = Array::make(nvars, 0)
        for i = 0; i < t.0.length(); i = i + 1 {
          e[i] = t.0[i]
        }
        (e, t.1)
      },
    ),
  )
}

///| A printable snapshot of the live state (used by the closed-loop report).
pub fn IdealState::describe(self : IdealState) -> String {
  let h = self.hilbert_fn(6)
  let hs : Array[String] = []
  for i = 0; i < h.length(); i = i + 1 {
    hs.push(h[i].to_string())
  }
  "I = (" +
  self.gens.map(fn(g) { g.to_string() }).join(", ") +
  ") in F_" +
  self.p.to_string() +
  "[x_0.." +
  (self.nvars - 1).to_string() +
  "]\n" +
  "  updates        : " +
  self.updates.to_string() +
  "\n" +
  "  S-pairs done   : " +
  self.pairs_done.to_string() +
  "\n" +
  "  Groebner basis : " +
  self.basis.map(fn(g) { g.to_string() }).join(" | ") +
  "\n" +
  "  finite quotient: " +
  self.finite.to_string() +
  "\n" +
  "  dim_Fp R/I     : " +
  self.dim.to_string() +
  "\n" +
  "  Krull dim R/I  : " +
  self.kdim.to_string() +
  "\n" +
  "  Hilbert fn     : [" +
  hs.join(",") +
  "]" +
  "\n" +
  "  maximal        : " +
  self.is_maximal().to_string()
}
