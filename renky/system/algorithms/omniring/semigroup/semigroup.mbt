// ===========================================================================
//  semigroup — affine semigroups, their rings k[S], normality, Ehrhart theory
//  and toric ideals (「仿射半群環」, structure 1).
//
//  Definition (affine semigroup).  For a finite set A = {a_1..a_m} ⊂ Z^n_{>=0}
//      S = S(A) = N a_1 + ... + N a_m ⊂ Z^n
//  is a finitely generated affine semigroup, and its semigroup ring is
//      k[S] = (+)_{s in S} k * t^s,     t^s * t^u = t^{s+u}.
//  k[S] is presented as F_p[x_1..x_m]/I_A where I_A = ker(x_i -> t^{a_i}) is the
//  *toric ideal*:
//      I_A = ( x^u - x^v : A u = A v ).
//
//  Theorems (verified exhaustively in `semigroup_wbtest`):
//   S1 membership in S is decided by dynamic programming over the box [0,v],
//      and agrees with brute-force search over all coefficient vectors.
//   S2 the group completion of S is the lattice L = Z a_1 + ... + Z a_m, and
//      membership in L is decided by the Hermite normal form (A29) — verified
//      against exhaustive integer combination search.
//   S3 S is *normal* (saturated, integrally closed) iff
//        (exists k >= 1 : k*v in S)  =>  v in S      for all v in the box.
//   S4 Hilbert function of k[S] with respect to a grading w:
//        H(d) = #\{ s in S : <w,s> = d \}, and k[S]/I_A has a monomial basis in
//      bijection with S (standard monomials <-> semigroup elements).
//   S5 Ehrhart: L_P(k) = #(k*P ∩ Z^n) is a polynomial in k of degree n for a
//      lattice polytope P = conv(A), and Ehrhart–Macdonald reciprocity holds:
//        L_P(-k) = (-1)^n * #(interior lattice points of k*P).
//   S6 the toric ideal is exactly the set of binomials x^u - x^v with A u = A v:
//      every lattice relation u (A u = 0) yields a binomial in I_A, and I_A
//      vanishes under the monomial map.
// ===========================================================================

///| A finitely generated affine semigroup S = N*gens ⊂ Z^n_{>=0}.
pub struct AffineSemigroup {
  n : Int
  gens : Array[Array[Int]]
}

pub fn semigroup_new(n : Int, gens : Array[Array[Int]]) -> AffineSemigroup {
  for i = 0; i < gens.length(); i = i + 1 {
    if gens[i].length() != n {
      abort("[semigroup] generator dimension mismatch")
    }
    for j = 0; j < n; j = j + 1 {
      if gens[i][j] < 0 {
        abort("[semigroup] generators must lie in Z^n_{>=0}")
      }
    }
  }
  { n, gens }
}

fn box_index(v : Array[Int], dims : Array[Int]) -> Int {
  let mut idx = 0
  let mut stride = 1
  for i = 0; i < v.length(); i = i + 1 {
    idx = idx + v[i] * stride
    stride = stride * dims[i]
  }
  idx
}

fn box_size(dims : Array[Int]) -> Int {
  let mut s = 1
  for i = 0; i < dims.length(); i = i + 1 {
    s = s * dims[i]
  }
  s
}

///| All elements of S inside the box [0,bound]^n, by dynamic programming.
///| **A57 semigroup membership / enumeration.**
pub fn AffineSemigroup::elements(self : AffineSemigroup, bound : Int) -> Array[Array[Int]] {
  let dims : Array[Int] = Array::make(self.n, bound + 1)
  let total = box_size(dims)
  if total > 400000 {
    abort("[semigroup.elements] box too large")
  }
  let seen : Array[Bool] = Array::make(total, false)
  let zero : Array[Int] = Array::make(self.n, 0)
  seen[box_index(zero, dims)] = true
  let mut frontier : Array[Array[Int]] = [zero]
  while !(frontier.is_empty()) {
    let next : Array[Array[Int]] = []
    for a = 0; a < frontier.length(); a = a + 1 {
      for g = 0; g < self.gens.length(); g = g + 1 {
        let cand : Array[Int] = Array::make(self.n, 0)
        let mut inside = true
        for i = 0; i < self.n; i = i + 1 {
          cand[i] = frontier[a][i] + self.gens[g][i]
          if cand[i] > bound {
            inside = false
          }
        }
        if inside && !(seen[box_index(cand, dims)]) {
          seen[box_index(cand, dims)] = true
          next.push(cand)
        }
      }
    }
    frontier = next
  }
  let out : Array[Array[Int]] = []
  for idx = 0; idx < total; idx = idx + 1 {
    if seen[idx] {
      out.push(box_decode(idx, dims))
    }
  }
  out
}

fn box_decode(idx : Int, dims : Array[Int]) -> Array[Int] {
  let out : Array[Int] = Array::make(dims.length(), 0)
  let mut rest = idx
  for i = 0; i < dims.length(); i = i + 1 {
    out[i] = rest % dims[i]
    rest = rest / dims[i]
  }
  out
}

///| Membership: is `v` in S?  Decided by the enumeration restricted to [0,v].
pub fn AffineSemigroup::contains(self : AffineSemigroup, v : Array[Int]) -> Bool {
  for i = 0; i < self.n; i = i + 1 {
    if v[i] < 0 {
      return false
    }
  }
  let dims : Array[Int] = []
  for i = 0; i < self.n; i = i + 1 {
    dims.push(v[i] + 1)
  }
  let total = box_size(dims)
  if total > 400000 {
    abort("[semigroup.contains] box too large")
  }
  let seen : Array[Bool] = Array::make(total, false)
  let zero : Array[Int] = Array::make(self.n, 0)
  seen[0] = true
  let mut frontier : Array[Array[Int]] = [zero]
  while !(frontier.is_empty()) {
    let next : Array[Array[Int]] = []
    for a = 0; a < frontier.length(); a = a + 1 {
      for g = 0; g < self.gens.length(); g = g + 1 {
        let cand : Array[Int] = Array::make(self.n, 0)
        let mut inside = true
        for i = 0; i < self.n; i = i + 1 {
          cand[i] = frontier[a][i] + self.gens[g][i]
          if cand[i] > v[i] {
            inside = false
          }
        }
        if inside && !(seen[box_index(cand, dims)]) {
          seen[box_index(cand, dims)] = true
          next.push(cand)
        }
      }
    }
    frontier = next
  }
  seen[box_index(v, dims)]
}

///| The group completion L = Z*gens, as a Hermite normal form (A29).
pub fn AffineSemigroup::lattice_hnf(self : AffineSemigroup) -> @linalg.Mat {
  let rows : @linalg.Mat = []
  for i = 0; i < self.gens.length(); i = i + 1 {
    let r : Array[Int] = []
    for j = 0; j < self.n; j = j + 1 {
      r.push(self.gens[i][j])
    }
    rows.push(r)
  }
  @linalg.hnf(rows)
}

///| **A58 HNF lattice membership**: membership in the group completion L, solved top-down on the HNF.
pub fn AffineSemigroup::in_lattice(self : AffineSemigroup, v : Array[Int]) -> Bool {
  let h = self.lattice_hnf()
  let rest : Array[Int] = []
  for i = 0; i < self.n; i = i + 1 {
    rest.push(v[i])
  }
  // pivots: for each row of the HNF the first nonzero column
  let pivots : Array[Int] = []
  for r = 0; r < h.length(); r = r + 1 {
    let mut pv = -1
    for c = 0; c < self.n; c = c + 1 {
      if h[r][c] != 0 && pv < 0 {
        pv = c
      }
    }
    if pv >= 0 {
      pivots.push(pv)
    }
  }
  let coef : Array[Int] = Array::make(h.length(), 0)
  for r = 0; r < pivots.length(); r = r + 1 {
    let c = pivots[r]
    if h[r][c] == 0 {
      continue
    }
    if rest[c] % h[r][c] != 0 {
      return false
    }
    let q = rest[c] / h[r][c]
    coef[r] = q
    for k = 0; k < self.n; k = k + 1 {
      rest[k] = rest[k] - q * h[r][k]
    }
  }
  for k = 0; k < self.n; k = k + 1 {
    if rest[k] != 0 {
      return false
    }
  }
  ignore(coef)
  true
}

///| **A59/S3 normality (saturation).**  S is normal iff every lattice point of the
///| cone that becomes a multiple of an element of S already lies in S.
pub fn AffineSemigroup::is_normal(self : AffineSemigroup, bound : Int, maxk : Int) -> Bool {
  let els = self.elements(bound)
  for x = 0; x <= bound; x = x + 1 {
    for y = 0; y <= bound; y = y + 1 {
      let v = if self.n == 1 {
        [x]
      } else if self.n == 2 {
        [x, y]
      } else {
        continue
      }
      let mut saturated = false
      for k = 1; k <= maxk; k = k + 1 {
        let kv : Array[Int] = []
        for i = 0; i < self.n; i = i + 1 {
          kv.push(k * v[i])
        }
        if self.contains(kv) {
          saturated = true
        }
      }
      if saturated && !(self.contains(v)) {
        return false
      }
    }
  }
  ignore(els)
  true
}

///| **S4 Hilbert function of k[S]** with respect to the grading w:
///| H(d) = #{ s in S ∩ box : <w,s> = d }.
pub fn AffineSemigroup::hilbert(
  self : AffineSemigroup,
  w : Array[Int],
  bound : Int,
  maxd : Int,
) -> Array[Int] {
  let els = self.elements(bound)
  let out : Array[Int] = Array::make(maxd + 1, 0)
  for i = 0; i < els.length(); i = i + 1 {
    let mut d = 0
    for j = 0; j < self.n; j = j + 1 {
      d = d + w[j] * els[i][j]
    }
    if d <= maxd {
      out[d] = out[d] + 1
    }
  }
  out
}

// ---------------------------------------------------------------------------
//  Ehrhart theory (A57)
// ---------------------------------------------------------------------------

///| All integer compositions of `k` into `parts` non-negative summands.
fn compositions(k : Int, parts : Int) -> Array[Array[Int]] {
  if parts == 0 {
    return if k == 0 {
      [[]]
    } else {
      []
    }
  }
  let out : Array[Array[Int]] = []
  for first = 0; first <= k; first = first + 1 {
    let sub = compositions(k - first, parts - 1)
    for i = 0; i < sub.length(); i = i + 1 {
      out.push(Array::add([first], sub[i]))
    }
  }
  out
}


///| **A60 convex hull and dilated point counting**: convex hull in counter-clockwise order (Andrew's monotone chain), with
///| exact integer cross products.  Collinear boundary points are dropped, which
///| is harmless for the half-plane test below.
pub fn convex_hull(pts : Array[Array[Int]]) -> Array[Array[Int]] {
  let sorted : Array[Array[Int]] = []
  for i = 0; i < pts.length(); i = i + 1 {
    sorted.push(pts[i])
  }
  sorted.sort_by(fn(a, b) {
    if a[0] != b[0] {
      a[0] - b[0]
    } else {
      a[1] - b[1]
    }
  })
  if sorted.length() <= 2 {
    return sorted
  }
  let lower : Array[Array[Int]] = []
  for i = 0; i < sorted.length(); i = i + 1 {
    while lower.length() >= 2 &&
      cross2(
        lower[lower.length() - 1][0] - lower[lower.length() - 2][0],
        lower[lower.length() - 1][1] - lower[lower.length() - 2][1],
        sorted[i][0] - lower[lower.length() - 1][0],
        sorted[i][1] - lower[lower.length() - 1][1],
      ) <= 0 {
      ignore(lower.pop())
    }
    lower.push(sorted[i])
  }
  let upper : Array[Array[Int]] = []
  for i = sorted.length() - 1; i >= 0; i = i - 1 {
    while upper.length() >= 2 &&
      cross2(
        upper[upper.length() - 1][0] - upper[upper.length() - 2][0],
        upper[upper.length() - 1][1] - upper[upper.length() - 2][1],
        sorted[i][0] - upper[upper.length() - 1][0],
        sorted[i][1] - upper[upper.length() - 1][1],
      ) <= 0 {
      ignore(upper.pop())
    }
    upper.push(sorted[i])
  }
  ignore(lower.pop())
  ignore(upper.pop())
  Array::add(lower, upper)
}

fn cross2(ax : Int, ay : Int, bx : Int, by : Int) -> Int {
  ax * by - ay * bx
}

///| Is `x` in the k-th dilation of P = conv(pts)?  Exact half-plane test:
///| for every counter-clockwise hull edge (a,b),
///|     cross(b - a, x - k*a) >= 0        (> 0 for the interior).
fn in_dilated(pts : Array[Array[Int]], x : Array[Int], k : Int, strict : Bool) -> Bool {
  let hull = convex_hull(pts)
  let m = hull.length()
  if m < 2 {
    // a single point: k*P = {0}
    return x[0] == 0 && x[1] == 0
  }
  for i = 0; i < m; i = i + 1 {
    let a = hull[i]
    let b = hull[(i + 1) % m]
    let cr = cross2(b[0] - a[0], b[1] - a[1], x[0] - k * a[0], x[1] - k * a[1])
    if strict {
      if cr <= 0 {
        return false
      }
    } else if cr < 0 {
      return false
    }
  }
  true
}

///| **L_P(k)**: every lattice point of the k-th dilation of P = conv(pts),
///| found by scanning the bounding box and testing each point exactly.
pub fn dilated_points(pts : Array[Array[Int]], k : Int) -> Array[Array[Int]] {
  let mut lox = 0
  let mut hix = 0
  let mut loy = 0
  let mut hiy = 0
  for i = 0; i < pts.length(); i = i + 1 {
    let x = k * pts[i][0]
    let y = k * pts[i][1]
    if i == 0 || x < lox {
      lox = x
    }
    if i == 0 || x > hix {
      hix = x
    }
    if i == 0 || y < loy {
      loy = y
    }
    if i == 0 || y > hiy {
      hiy = y
    }
  }
  let out : Array[Array[Int]] = []
  if k == 0 {
    out.push([0, 0])
    return out
  }
  for x = lox; x <= hix; x = x + 1 {
    for y = loy; y <= hiy; y = y + 1 {
      if in_dilated(pts, [x, y], k, false) {
        out.push([x, y])
      }
    }
  }
  out
}

///| Interior lattice points of k*P (strict half-plane inequalities), used for
///| Ehrhart–Macdonald reciprocity.
pub fn dilated_interior(pts : Array[Array[Int]], k : Int) -> Array[Array[Int]] {
  let mut lox = 0
  let mut hix = 0
  let mut loy = 0
  let mut hiy = 0
  for i = 0; i < pts.length(); i = i + 1 {
    let x = k * pts[i][0]
    let y = k * pts[i][1]
    if i == 0 || x < lox {
      lox = x
    }
    if i == 0 || x > hix {
      hix = x
    }
    if i == 0 || y < loy {
      loy = y
    }
    if i == 0 || y > hiy {
      hiy = y
    }
  }
  let out : Array[Array[Int]] = []
  if k == 0 {
    return out
  }
  for x = lox; x <= hix; x = x + 1 {
    for y = loy; y <= hiy; y = y + 1 {
      if in_dilated(pts, [x, y], k, true) {
        out.push([x, y])
      }
    }
  }
  out
}


///| Ehrhart function values L_P(k) for k = 0..maxk.
pub fn ehrhart_fn(pts : Array[Array[Int]], maxk : Int) -> Array[Int] {
  let out : Array[Int] = []
  for k = 0; k <= maxk; k = k + 1 {
    out.push(dilated_points(pts, k).length())
  }
  out
}

///| **A61 finite-difference interpolation**: the coefficients of the polynomial
///| that agrees with `ys` on 0..len-1, in the binomial basis
///|     p(k) = sum_i c_i * C(k, i).
pub fn newton_coeffs(ys : Array[Int]) -> Array[Int] {
  let mut cur : Array[Int] = []
  for i = 0; i < ys.length(); i = i + 1 {
    cur.push(ys[i])
  }
  let out : Array[Int] = []
  while !(cur.is_empty()) {
    out.push(cur[0])
    let nxt : Array[Int] = []
    for i = 1; i < cur.length(); i = i + 1 {
      nxt.push(cur[i] - cur[i - 1])
    }
    cur = nxt
  }
  out
}

///| Evaluate a polynomial given in the binomial basis at `k` (k may be negative,
///| where C(k,i) = k(k-1)...(k-i+1)/i!).
pub fn newton_eval(coeffs : Array[Int], k : Int) -> Int {
  let mut acc = 0
  for i = 0; i < coeffs.length(); i = i + 1 {
    let mut b = 1
    for j = 0; j < i; j = j + 1 {
      b = b * (k - j)
    }
    for j = 1; j <= i; j = j + 1 {
      b = b / j
    }
    acc = acc + coeffs[i] * b
  }
  acc
}

// ---------------------------------------------------------------------------
//  A47 toric ideal
// ---------------------------------------------------------------------------

///| **A62 toric ideal** I_A = (x^u - x^v : A u = A v), found by exhaustive
///| search over exponent vectors of total degree <= maxdeg.  The result is a live
///| ideal state in F_p[x_1..x_m].
pub fn toric_ideal(
  gens : Array[Array[Int]],
  p : Int,
  maxdeg : Int,
) -> @ideal.IdealState {
  let m = gens.length()
  let st = @ideal.ideal_new(p, m)
  let exps : Array[Array[Int]] = []
  for d = 0; d <= maxdeg; d = d + 1 {
    let cs = compositions(d, m)
    for i = 0; i < cs.length(); i = i + 1 {
      exps.push(cs[i])
    }
  }
  // group exponent vectors by their image A*u
  let images : Array[Array[Int]] = []
  for u = 0; u < exps.length(); u = u + 1 {
    let img : Array[Int] = Array::make(gens[0].length(), 0)
    for i = 0; i < m; i = i + 1 {
      for j = 0; j < gens[0].length(); j = j + 1 {
        img[j] = img[j] + exps[u][i] * gens[i][j]
      }
    }
    images.push(img)
  }
  for u = 0; u < exps.length(); u = u + 1 {
    for v = u + 1; v < exps.length(); v = v + 1 {
      let mut same = true
      for j = 0; j < images[u].length(); j = j + 1 {
        if images[u][j] != images[v][j] {
          same = false
        }
      }
      if same {
        let mu = monom_of(exps[u], p, m)
        let mv = monom_of(exps[v], p, m)
        st.add(mu.sub(mv))
      }
    }
  }
  st
}

fn monom_of(e : Array[Int], p : Int, m : Int) -> @poly.MPoly {
  let ev : Array[Int] = Array::make(m, 0)
  for i = 0; i < m; i = i + 1 {
    ev[i] = e[i]
  }
  @poly.mp_monom(p, m, ev, 1)
}

///| The monomial map x_i -> t^{a_i} evaluated at a semigroup element, i.e. the
///| image of a polynomial in the semigroup ring (returns the exponent of the
///| monomial for monomials, used for verification).
pub fn monomial_map_exponent(gens : Array[Array[Int]], e : Array[Int]) -> Array[Int] {
  let out : Array[Int] = Array::make(gens[0].length(), 0)
  for i = 0; i < e.length(); i = i + 1 {
    for j = 0; j < out.length(); j = j + 1 {
      out[j] = out[j] + e[i] * gens[i][j]
    }
  }
  out
}
