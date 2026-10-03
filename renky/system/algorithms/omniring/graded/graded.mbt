// ===========================================================================
//  graded — graded rings, Hilbert series, Veronese subrings, associated graded
//  rings and the Koszul (exterior) algebra (「分次環」, structure 6).
//
//  Definition (graded algebra).  A commutative F_p-algebra A is *graded* when
//  A = (+)_d A_d with A_i * A_j subset A_{i+j}.  Here a finite dimensional
//  graded algebra is given by a homogeneous F_p-basis `basis` with
//  `degrees[i] = deg(basis_i)` and structure constants
//      basis_i * basis_j = sum_k mult[i][j][k] * basis_k.
//
//  Theorems (all verified exhaustively in `graded_wbtest`):
//   G1 grading:      mult[i][j][k] != 0  =>  degrees[k] = degrees[i] + degrees[j]
//   G2 associativity, commutativity, unitality and bilinearity of the
//      structure constants
//   G3 Hilbert function: H(d) = #\{i : degrees[i] = d\} = dim_{F_p} A_d, and
//      dim A = sum_d H(d)
//   G4 Hilbert–Serre: for a standard graded algebra the Hilbert function
//      agrees, for d >> 0, with a polynomial whose degree is dim(A) - 1
//      (Krull dimension minus one); detected here by finite differences.
//   G5 Veronese: A^{(k)} = (+)_d A_{kd} is a graded subalgebra and
//      H_{A^{(k)}}(d) = H_A(kd).
//   G6 associated graded: for an ideal I of a finite dimensional algebra,
//      gr_I(A) = (+)_k I^k/I^{k+1} has sum_k dim(I^k/I^{k+1}) = dim A.
//   G7 Koszul/exterior: Lambda(V) is graded-commutative with
//      a*b = (-1)^{deg a deg b} b*a, a*a = 0 for odd a, and H(d) = C(n,d).
// ===========================================================================

///| A finite dimensional graded commutative F_p-algebra by structure constants.
pub struct GradedAlg {
  p : Int
  n : Int
  degrees : Array[Int]
  mult : Array[Array[Array[Int]]]
  one : Int
}

///| Build the graded algebra R/I from a live ideal state, where variable `v`
///| has degree `w[v]`.
pub fn graded_from_ideal(
  st : @ideal.IdealState,
  w : Array[Int],
) -> GradedAlg {
  let p = st.p
  let bas = st.quotient_basis()
  let n = bas.length()
  let degrees : Array[Int] = []
  for i = 0; i < n; i = i + 1 {
    let mut d = 0
    for v = 0; v < bas[i].length(); v = v + 1 {
      d = d + w[v] * bas[i][v]
    }
    degrees.push(d)
  }
  let mult : Array[Array[Array[Int]]] = []
  for i = 0; i < n; i = i + 1 {
    let row : Array[Array[Int]] = []
    for j = 0; j < n; j = j + 1 {
      let bi = @poly.mp_monom(p, st.nvars, bas[i], 1)
      let bj = @poly.mp_monom(p, st.nvars, bas[j], 1)
      let prod = st.normal_form(bi.mul(bj))
      row.push(coeffs_in_basis(prod, bas, p))
    }
    mult.push(row)
  }
  // index of the unit: the basis element of degree 0 with all exponents zero
  let mut unit = -1
  for i = 0; i < n; i = i + 1 {
    let mut z = true
    for v = 0; v < bas[i].length(); v = v + 1 {
      if bas[i][v] != 0 {
        z = false
      }
    }
    if z {
      unit = i
    }
  }
  { p, n, degrees, mult, one: unit }
}

///| Coefficient vector of `f` (already in normal form) on the monomial basis.
fn coeffs_in_basis(
  f : @poly.MPoly,
  bas : Array[Array[Int]],
  p : Int,
) -> Array[Int] {
  let out : Array[Int] = Array::make(bas.length(), 0)
  for t = 0; t < f.terms.length(); t = t + 1 {
    let (e, c) = f.terms[t]
    for i = 0; i < bas.length(); i = i + 1 {
      if @poly.cmp_monom(bas[i], e) == 0 {
        out[i] = c % p
      }
    }
  }
  out
}

///| Hilbert function H(d) = dim A_d for d = 0..maxdeg.
pub fn GradedAlg::hilbert(self : GradedAlg, maxdeg : Int) -> Array[Int] {
  let out : Array[Int] = Array::make(maxdeg + 1, 0)
  for i = 0; i < self.n; i = i + 1 {
    if self.degrees[i] <= maxdeg {
      out[self.degrees[i]] = out[self.degrees[i]] + 1
    }
  }
  out
}

///| **A52 Hilbert series** as the coefficient array of sum_d H(d) t^d, plus the
///| cumulative dimension.
pub fn GradedAlg::hilbert_series(self : GradedAlg, maxdeg : Int) -> (Array[Int], Int) {
  let h = self.hilbert(maxdeg)
  let mut total = 0
  for d = 0; d <= maxdeg; d = d + 1 {
    total = total + h[d]
  }
  (h, total)
}

///| **A52 (cont.) Hilbert polynomial detection** by finite differences: the largest k
///| such that the k-th difference of H is eventually nonzero-free.  For a
///| standard graded algebra of Krull dimension d the answer is d - 1.
pub fn hilbert_poly_degree(h : Array[Int], start : Int) -> Int {
  let mut cur : Array[Int] = []
  for i = start; i < h.length(); i = i + 1 {
    cur.push(h[i])
  }
  let mut deg = -1
  while cur.length() > 1 {
    let mut constant = true
    for i = 1; i < cur.length(); i = i + 1 {
      if cur[i] != cur[0] {
        constant = false
      }
    }
    if constant {
      deg = deg + 1
      // constant sequence: the polynomial has degree `deg`
      return deg
    }
    deg = deg + 1
    let nxt : Array[Int] = []
    for i = 1; i < cur.length(); i = i + 1 {
      nxt.push(cur[i] - cur[i - 1])
    }
    cur = nxt
  }
  deg
}

///| **A53 Veronese subalgebra** A^{(k)} = (+)_d A_{kd}.
pub fn GradedAlg::veronese(self : GradedAlg, k : Int) -> GradedAlg {
  let idx : Array[Int] = []
  for i = 0; i < self.n; i = i + 1 {
    if self.degrees[i] % k == 0 {
      idx.push(i)
    }
  }
  let m = idx.length()
  let degrees : Array[Int] = Array::make(m, 0)
  for a = 0; a < m; a = a + 1 {
    degrees[a] = self.degrees[idx[a]] / k
  }
  let mult : Array[Array[Array[Int]]] = []
  for a = 0; a < m; a = a + 1 {
    let row : Array[Array[Int]] = []
    for b = 0; b < m; b = b + 1 {
      let src = self.mult[idx[a]][idx[b]]
      let dst : Array[Int] = Array::make(m, 0)
      for t = 0; t < m; t = t + 1 {
        dst[t] = src[idx[t]]
      }
      row.push(dst)
    }
    mult.push(row)
  }
  let mut unit = 0
  for a = 0; a < m; a = a + 1 {
    if idx[a] == self.one {
      unit = a
    }
  }
  { p: self.p, n: m, degrees, mult, one: unit }
}

///| Truncation A_{>= k}: the structure constants restricted to basis elements of
///| degree >= k (a non-unital graded algebra).
pub fn GradedAlg::truncate(self : GradedAlg, k : Int) -> GradedAlg {
  let idx : Array[Int] = []
  for i = 0; i < self.n; i = i + 1 {
    if self.degrees[i] >= k {
      idx.push(i)
    }
  }
  let m = idx.length()
  let degrees : Array[Int] = Array::make(m, 0)
  for a = 0; a < m; a = a + 1 {
    degrees[a] = self.degrees[idx[a]]
  }
  let mult : Array[Array[Array[Int]]] = []
  for a = 0; a < m; a = a + 1 {
    let row : Array[Array[Int]] = []
    for b = 0; b < m; b = b + 1 {
      let src = self.mult[idx[a]][idx[b]]
      let dst : Array[Int] = Array::make(m, 0)
      for t = 0; t < m; t = t + 1 {
        dst[t] = src[idx[t]]
      }
      row.push(dst)
    }
    mult.push(row)
  }
  { p: self.p, n: m, degrees, mult, one: -1 }
}

///| **A54 power of an ideal** inside the live state: I^k is generated by all
///| products of k generators.
pub fn ideal_power(st : @ideal.IdealState, k : Int) -> @ideal.IdealState {
  let out = @ideal.ideal_new(st.p, st.nvars)
  if k == 0 {
    out.add(@poly.mp_one(st.p, st.nvars))
    return out
  }
  let gens = st.gens
  if gens.length() == 0 {
    return out
  }
  let mut cur : Array[@poly.MPoly] = []
  for i = 0; i < gens.length(); i = i + 1 {
    cur.push(gens[i])
  }
  for step = 1; step < k; step = step + 1 {
    let nxt : Array[@poly.MPoly] = []
    for a = 0; a < cur.length(); a = a + 1 {
      for b = 0; b < gens.length(); b = b + 1 {
        nxt.push(cur[a].mul(gens[b]))
      }
    }
    cur = nxt
  }
  for i = 0; i < cur.length(); i = i + 1 {
    out.add(cur[i])
  }
  out
}

///| **A55/G6 associated graded (adic filtration)** dims: with D_k = dim_{F_p} R/J^k for J = (f) inside
///| the finite dimensional quotient R = F_p[x]/I, return
///|     dim(J^k/J^{k+1}) = D_{k+1} - D_k,   k = 0..maxk-1.
pub fn adic_grading(
  st : @ideal.IdealState,
  f : @poly.MPoly,
  maxk : Int,
) -> Array[Int] {
  let dims : Array[Int] = []
  let mut fk = @poly.mp_one(st.p, st.nvars)
  for k = 0; k <= maxk; k = k + 1 {
    let cur = @ideal.ideal_new(st.p, st.nvars)
    for i = 0; i < st.gens.length(); i = i + 1 {
      cur.add(st.gens[i])
    }
    cur.add(fk)
    dims.push(cur.dim)
    fk = fk.mul(f)
    ignore(k)
  }
  let out : Array[Int] = []
  for k = 0; k < maxk; k = k + 1 {
    out.push(dims[k + 1] - dims[k])
  }
  out
}

///| **A56/G7 Koszul (exterior) algebra** Lambda(V) on n generators over F_p,
///| graded by wedge degree.  Basis elements are subsets, encoded by bitmasks;
///| multiplication is the signed wedge.
pub fn exterior_alg(p : Int, n : Int) -> GradedAlg {
  let size = 1 << n
  let degrees : Array[Int] = Array::make(size, 0)
  for m = 0; m < size; m = m + 1 {
    degrees[m] = popcount(m)
  }
  let mult : Array[Array[Array[Int]]] = []
  for a = 0; a < size; a = a + 1 {
    let row : Array[Array[Int]] = []
    for b = 0; b < size; b = b + 1 {
      let v : Array[Int] = Array::make(size, 0)
      if (a & b) == 0 {
        let sgn = wedge_sign(a, b)
        v[a | b] = if sgn {
          1
        } else {
          ((-1) % p + p) % p
        }
      }
      row.push(v)
    }
    mult.push(row)
  }
  { p, n: size, degrees, mult, one: 0 }
}

fn popcount(x : Int) -> Int {
  let mut c = 0
  let mut y = x
  while y > 0 {
    c = c + (y & 1)
    y = y >> 1
  }
  c
}

///| Sign of the shuffle needed to sort the concatenated generators: the number
///| of inversions (i in a, j in b, i > j), taken mod 2.
fn wedge_sign(a : Int, b : Int) -> Bool {
  let mut inv = 0
  let mut j = 0
  while (1 << j) <= a || (1 << j) <= b {
    if ((b >> j) & 1) == 1 {
      let mut i = 0
      while (1 << i) <= a {
        if ((a >> i) & 1) == 1 && i > j {
          inv = inv + 1
        }
        i = i + 1
      }
    }
    j = j + 1
    if j > 30 {
      break
    }
  }
  inv % 2 == 0
}

///| Binomial coefficient C(n,k) by the multiplicative formula.
pub fn binom(n : Int, k : Int) -> Int {
  if k < 0 || k > n {
    return 0
  }
  let mut acc = 1
  for i = 0; i < k; i = i + 1 {
    acc = acc * (n - i) / (i + 1)
  }
  acc
}
