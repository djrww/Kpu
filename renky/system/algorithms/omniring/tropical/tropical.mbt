// ===========================================================================
//  tropical — the tropical (min-plus) semiring, tropical linear algebra, and
//  Newton polygon subdivision with balancing (「熱帶半環」, structure 3).
//
//  Definition (tropical semiring).  T = Z ∪ {+∞} with
//      a ⊕ b = min(a, b),      a ⊗ b = a + b,      0_⊕ = +∞,      1_⊗ = 0.
//  Axioms (all verified exhaustively in `tropical_wbtest`):
//      (T, ⊕) is an idempotent commutative monoid,
//      (T, ⊗) is a commutative monoid,
//      ⊗ distributes over ⊕,  0_⊕ ⊗ a = 0_⊕,  a ⊕ a = a.
//
//  Theorem (Maslov).  The heat semifields (R, (+)_h, (*)_h) of package `heat`
//  converge to (T, ⊕, ⊗) as h -> 0+; `trop_limit` is that limit and
//  `trop_dequantize` is the finite-h approximation.  (Numerically verified in
//  `heat` and cross-checked in `loop`.)
//
//  Definition (tropical polynomial / hypersurface).  For
//      f = ⊕_k  c_k ⊗ x^{e_k}   (c_k = valuations, e_k = exponent vectors)
//  the value at w is  f(w) = min_k (c_k + <e_k, w>) and w lies on the tropical
//  hypersurface V_trop(f) iff the minimum is attained at least twice.
//
//  Theorem (tropical determinant = optimal assignment).
//      det_trop(A) = min_{sigma in S_n} sum_i A[i][sigma(i)],
//  computed by the Hungarian algorithm (A48) and verified against complete
//  enumeration of S_n for n <= 5 and against random matrices for n = 6.
//
//  Theorem (Newton subdivision, A50).  The coherent subdivision of the Newton
//  polygon induced by the valuations is the set of argmin sets
//      C(w) = { k : c_k + <e_k,w> = min_j (c_j + <e_j,w>) },  w in R^2,
//  each of which is decided *exactly* by rational linear programming.
//
//  Theorem (balancing).  At every maximal cell C of the subdivision,
//      sum_{boundary edges {i,j} of C}  weight(i,j) * outward_normal(i,j,C) = 0,
//  where weight(i,j) = gcd(|de_x|, |de_y|) is the lattice length of the edge.
//  This is the balance condition making V_trop(f) a tropical curve.
// ===========================================================================

///| The additive zero of the min-plus semiring: +infinity.
pub let trop_zero : Int = 1000000000

///| The multiplicative unit of the min-plus semiring: 0.
pub let trop_one : Int = 0

///| ⊕ = min, absorbing +∞.
pub fn trop_add(a : Int, b : Int) -> Int {
  Int::min(a, b)
}

///| ⊗ = +, with ∞ absorbing.
pub fn trop_mul(a : Int, b : Int) -> Int {
  if a >= trop_zero || b >= trop_zero {
    trop_zero
  } else {
    a + b
  }
}

pub fn trop_is_zero(a : Int) -> Bool {
  a >= trop_zero
}

///| n-ary ⊕.
pub fn trop_sum(xs : Array[Int]) -> Int {
  let mut acc = trop_zero
  for i = 0; i < xs.length(); i = i + 1 {
    acc = trop_add(acc, xs[i])
  }
  acc
}

///| n-ary ⊗.
pub fn trop_prod(xs : Array[Int]) -> Int {
  let mut acc = trop_one
  for i = 0; i < xs.length(); i = i + 1 {
    acc = trop_mul(acc, xs[i])
  }
  acc
}

// ---------------------------------------------------------------------------
//  Tropical polynomials
// ---------------------------------------------------------------------------

///| A tropical polynomial: terms `(exponent vector, coefficient valuation)`,
///| kept sorted by graded-lex so that equality is syntactic.
pub struct TPoly {
  nvars : Int
  terms : Array[(Array[Int], Int)]
}

pub fn tp_make(
  nvars : Int,
  terms : Array[(Array[Int], Int)],
) -> TPoly {
  let kept : Array[(Array[Int], Int)] = []
  for i = 0; i < terms.length(); i = i + 1 {
    let (e, c) = terms[i]
    if trop_is_zero(c) {
      continue
    }
    let mut merged = false
    for j = 0; j < kept.length(); j = j + 1 {
      if exp_equal(kept[j].0, e) {
        kept[j] = (e, trop_add(kept[j].1, c))
        merged = true
      }
    }
    if !(merged) {
      kept.push((e, c))
    }
  }
  let fin : Array[(Array[Int], Int)] = []
  for i = 0; i < kept.length(); i = i + 1 {
    if !(trop_is_zero(kept[i].1)) {
      fin.push(kept[i])
    }
  }
  fin.sort_by(fn(a, b) { cmp_exp(a.0, b.0) })
  { nvars, terms: fin }
}

fn exp_equal(a : Array[Int], b : Array[Int]) -> Bool {
  if a.length() != b.length() {
    return false
  }
  for i = 0; i < a.length(); i = i + 1 {
    if a[i] != b[i] {
      return false
    }
  }
  true
}

fn cmp_exp(a : Array[Int], b : Array[Int]) -> Int {
  let da = deg_of(a)
  let db = deg_of(b)
  if da != db {
    return da - db
  }
  for i = 0; i < a.length(); i = i + 1 {
    if a[i] != b[i] {
      return b[i] - a[i]
    }
  }
  0
}

fn deg_of(e : Array[Int]) -> Int {
  let mut d = 0
  for i = 0; i < e.length(); i = i + 1 {
    d = d + e[i]
  }
  d
}

///| Tropical ⊕ of two tropical polynomials.
pub fn TPoly::add(self : TPoly, other : TPoly) -> TPoly {
  let all : Array[(Array[Int], Int)] = []
  for i = 0; i < self.terms.length(); i = i + 1 {
    all.push(self.terms[i])
  }
  for i = 0; i < other.terms.length(); i = i + 1 {
    all.push(other.terms[i])
  }
  tp_make(self.nvars, all)
}

///| Tropical ⊗ of two tropical polynomials: exponents add, valuations add,
///| and equal exponents combine with ⊕ = min.
pub fn TPoly::mul(self : TPoly, other : TPoly) -> TPoly {
  let all : Array[(Array[Int], Int)] = []
  for i = 0; i < self.terms.length(); i = i + 1 {
    for j = 0; j < other.terms.length(); j = j + 1 {
      let e : Array[Int] = Array::make(self.nvars, 0)
      for k = 0; k < self.nvars; k = k + 1 {
        e[k] = self.terms[i].0[k] + other.terms[j].0[k]
      }
      all.push((e, trop_mul(self.terms[i].1, other.terms[j].1)))
    }
  }
  tp_make(self.nvars, all)
}

///| f(w) = min_k (c_k + <e_k, w>).
pub fn TPoly::eval(self : TPoly, w : Array[Int]) -> Int {
  let mut acc = trop_zero
  for i = 0; i < self.terms.length(); i = i + 1 {
    let (e, c) = self.terms[i]
    let mut d = c
    for k = 0; k < e.length(); k = k + 1 {
      d = d + e[k] * w[k]
    }
    acc = trop_add(acc, d)
  }
  acc
}

///| Multiplicity of f at w: how many terms attain the minimum.
pub fn TPoly::multiplicity(self : TPoly, w : Array[Int]) -> Int {
  let m = self.eval(w)
  let mut cnt = 0
  for i = 0; i < self.terms.length(); i = i + 1 {
    let (e, c) = self.terms[i]
    let mut d = c
    for k = 0; k < e.length(); k = k + 1 {
      d = d + e[k] * w[k]
    }
    if d == m {
      cnt = cnt + 1
    }
  }
  cnt
}

///| w lies on the tropical hypersurface V_trop(f) iff the minimum is attained
///| at least twice.
pub fn TPoly::on_hypersurface(self : TPoly, w : Array[Int]) -> Bool {
  self.multiplicity(w) >= 2
}

pub fn TPoly::to_string(self : TPoly) -> String {
  let parts : Array[String] = []
  for i = 0; i < self.terms.length(); i = i + 1 {
    let (e, c) = self.terms[i]
    let mut mon = c.to_string()
    for k = 0; k < e.length(); k = k + 1 {
      if e[k] > 0 {
        mon = mon + "*x" + k.to_string() + "^" + e[k].to_string()
      }
    }
    parts.push(mon)
  }
  "min(" + parts.join(", ") + ")"
}

// ---------------------------------------------------------------------------
//  A36  tropical determinant = optimal assignment (Hungarian algorithm)
// ===========================================================================

///| **A48 Hungarian algorithm.**  Returns the optimal (minimum) assignment cost
///| and the assignment itself, for a square integer cost matrix.  O(n^3).
pub fn hungarian(cost : Array[Array[Int]]) -> (Int, Array[Int]) {
  let n = cost.length()
  if n == 0 {
    return (0, [])
  }
  let inf = 1 << 28
  // 1-based potentials as in the classical formulation
  let u : Array[Int] = Array::make(n + 1, 0)
  let v : Array[Int] = Array::make(n + 1, 0)
  let p : Array[Int] = Array::make(n + 1, 0)
  let way : Array[Int] = Array::make(n + 1, 0)
  for i = 1; i <= n; i = i + 1 {
    p[0] = i
    let mut j0 = 0
    let minv : Array[Int] = Array::make(n + 1, inf)
    let used : Array[Bool] = Array::make(n + 1, false)
    while true {
      used[j0] = true
      let i0 = p[j0]
      let mut delta = inf
      let mut j1 = -1
      for j = 1; j <= n; j = j + 1 {
        if !(used[j]) {
          let cur = cost[i0 - 1][j - 1] - u[i0] - v[j]
          if cur < minv[j] {
            minv[j] = cur
            way[j] = j0
          }
          if minv[j] < delta {
            delta = minv[j]
            j1 = j
          }
        }
      }
      for j = 0; j <= n; j = j + 1 {
        if used[j] {
          u[p[j]] = u[p[j]] + delta
          v[j] = v[j] - delta
        } else {
          minv[j] = minv[j] - delta
        }
      }
      j0 = j1
      if p[j0] == 0 {
        break
      }
    }
    while true {
      let j1 = way[j0]
      p[j0] = p[j1]
      j0 = j1
      if j0 == 0 {
        break
      }
    }
  }
  let ans : Array[Int] = Array::make(n, 0)
  for j = 1; j <= n; j = j + 1 {
    if p[j] > 0 {
      ans[p[j] - 1] = j - 1
    }
  }
  let mut total = 0
  for i = 0; i < n; i = i + 1 {
    total = total + cost[i][ans[i]]
  }
  (total, ans)
}

///| Brute-force tropical determinant: minimum over all n! permutations.
pub fn trop_det_bruteforce(m : Array[Array[Int]]) -> Int {
  let n = m.length()
  let idx : Array[Int] = []
  for i = 0; i < n; i = i + 1 {
    idx.push(i)
  }
  let mut best = trop_zero
  let perms = permutations(idx)
  for k = 0; k < perms.length(); k = k + 1 {
    let mut s = 0
    for i = 0; i < n; i = i + 1 {
      s = s + m[i][perms[k][i]]
    }
    if s < best {
      best = s
    }
  }
  best
}

fn permutations(xs : Array[Int]) -> Array[Array[Int]] {
  if xs.length() <= 1 {
    return [xs]
  }
  let out : Array[Array[Int]] = []
  for i = 0; i < xs.length(); i = i + 1 {
    let rest : Array[Int] = []
    for j = 0; j < xs.length(); j = j + 1 {
      if i != j {
        rest.push(xs[j])
      }
    }
    let sub = permutations(rest)
    for k = 0; k < sub.length(); k = k + 1 {
      let one : Array[Int] = [xs[i]]
      out.push(Array::add(one, sub[k]))
    }
  }
  out
}

///| Tropical determinant via the Hungarian algorithm.
pub fn trop_det(m : Array[Array[Int]]) -> Int {
  hungarian(m).0
}

///| **A49 tropical (min-plus) matrix product**.
pub fn trop_mat_mul(a : Array[Array[Int]], b : Array[Array[Int]]) -> Array[Array[Int]] {
  let n = a.length()
  let m = b[0].length()
  let k = b.length()
  let out : Array[Array[Int]] = []
  for i = 0; i < n; i = i + 1 {
    let row : Array[Int] = Array::make(m, trop_zero)
    for j = 0; j < m; j = j + 1 {
      let mut acc = trop_zero
      for t = 0; t < k; t = t + 1 {
        acc = trop_add(acc, trop_mul(a[i][t], b[t][j]))
      }
      row[j] = acc
    }
    out.push(row)
  }
  out
}

///| **A49 (cont.) all-pairs tropical shortest paths** (Floyd–Warshall in the min-plus
///| semiring): the tropical closure of a weight matrix.
pub fn trop_closure(w : Array[Array[Int]]) -> Array[Array[Int]] {
  let n = w.length()
  let d : Array[Array[Int]] = []
  for i = 0; i < n; i = i + 1 {
    let row : Array[Int] = Array::make(n, trop_zero)
    for j = 0; j < n; j = j + 1 {
      row[j] = w[i][j]
    }
    d.push(row)
  }
  for k = 0; k < n; k = k + 1 {
    for i = 0; i < n; i = i + 1 {
      for j = 0; j < n; j = j + 1 {
        let cand = trop_mul(d[i][k], d[k][j])
        if cand < d[i][j] {
          d[i][j] = cand
        }
      }
    }
  }
  d
}

// ---------------------------------------------------------------------------
//  A50/A51  Newton polygon subdivision (exact rational arithmetic) and balancing
// ---------------------------------------------------------------------------

///| **A50 Newton subdivision**: cells of the Newton subdivision — each cell is the list of term indices that
///| attain the minimum for some weight vector w.  Decided exactly, with
///| rational arithmetic realised as integer cross products.
pub fn newton_cells(terms : Array[(Array[Int], Int)]) -> Array[Array[Int]] {
  let n = terms.length()
  let cells : Array[Array[Int]] = []
  if n == 0 {
    return cells
  }
  if n == 1 {
    cells.push([0])
    return cells
  }
  for i = 0; i < n; i = i + 1 {
    for j = i + 1; j < n; j = j + 1 {
      let found = pair_cells(terms, i, j)
      for t = 0; t < found.length(); t = t + 1 {
        if !(contains_cell(cells, found[t])) {
          cells.push(found[t])
        }
      }
    }
  }
  // singletons that no pair line reaches are not lower-hull vertices; keep them
  // so that the cell collection covers every point
  for i = 0; i < n; i = i + 1 {
    let mut seen = false
    for k = 0; k < cells.length(); k = k + 1 {
      for t = 0; t < cells[k].length(); t = t + 1 {
        if cells[k][t] == i {
          seen = true
        }
      }
    }
    if !(seen) {
      cells.push([i])
    }
  }
  cells.sort_by(fn(a, b) {
    if a.length() != b.length() {
      a.length() - b.length()
    } else if a.is_empty() {
      0
    } else {
      a[0] - b[0]
    }
  })
  cells
}

fn contains_cell(cells : Array[Array[Int]], c : Array[Int]) -> Bool {
  for k = 0; k < cells.length(); k = k + 1 {
    if cells[k].length() == c.length() {
      let mut same = true
      for t = 0; t < c.length(); t = t + 1 {
        if cells[k][t] != c[t] {
          same = false
        }
      }
      if same {
        return true
      }
    }
  }
  false
}

///| All argmin cells along the line of weights where terms `i` and `j` tie.
///|
///| Writing w(s) = ((c*nx - s*ny)/d, (c*ny + s*nx)/d) the condition
///| "i and j are both minimal" is an interval in s, cut into chambers by the
///| finitely many roots where some other term becomes minimal.  Every chamber
///| endpoint, chamber midpoint and both unbounded rays is sampled, which
///| enumerates *all* cells that the line meets.
fn pair_cells(
  terms : Array[(Array[Int], Int)],
  i : Int,
  j : Int,
) -> Array[Array[Int]] {
  let out : Array[Array[Int]] = []
  let ei = terms[i].0
  let ej = terms[j].0
  let vi = terms[i].1
  let vj = terms[j].1
  if exp_equal(ei, ej) {
    if vi == vj {
      // identical lifted points: they always tie, sample a few weights
      for a = -3; a <= 3; a = a + 1 {
        for b = -3; b <= 3; b = b + 1 {
          let c = argmin_cell(terms, a, b, 1)
          if !(contains_cell(out, c)) {
            out.push(c)
          }
        }
      }
    }
    return out
  }
  let nx0 = ei[0] - ej[0]
  let ny0 = ei[1] - ej[1]
  let d = nx0 * nx0 + ny0 * ny0
  let c = vj - vi
  let big = 1 << 20
  let mut lo_n = -big
  let mut lo_d = 1
  let mut hi_n = big
  let mut hi_d = 1
  let roots : Array[(Int, Int)] = []
  for k = 0; k < terms.length(); k = k + 1 {
    let ek = terms[k].0
    let vk = terms[k].1
    let a = d * (vk - vi) + c * (nx0 * ek[0] + ny0 * ek[1] - nx0 * ei[0] - ny0 * ei[1])
    let b = nx0 * (ek[1] - ei[1]) - ny0 * (ek[0] - ei[0])
    if b == 0 {
      if a < 0 {
        return out
      }
    } else if b > 0 {
      let cn = -a
      let cd = b
      roots.push((cn, cd))
      if cn * lo_d > lo_n * cd {
        lo_n = cn
        lo_d = cd
      }
    } else {
      let cn = a
      let cd = -b
      roots.push((cn, cd))
      if cn * hi_d < hi_n * cd {
        hi_n = cn
        hi_d = cd
      }
    }
  }
  if lo_n * hi_d > hi_n * lo_d {
    return out
  }
  // candidate values of s: feasible roots, midpoints between consecutive ones,
  // the two endpoints, and one point beyond each endpoint
  let cands : Array[(Int, Int)] = []
  let inside : Array[(Int, Int)] = []
  for t = 0; t < roots.length(); t = t + 1 {
    let (rn, rd) = roots[t]
    if rn * lo_d >= lo_n * rd && rn * hi_d <= hi_n * rd {
      inside.push((rn, rd))
    }
  }
  inside.sort_by(fn(x, y) { x.0 * y.1 - y.0 * x.1 })
  cands.push((lo_n, lo_d))
  cands.push((hi_n, hi_d))
  cands.push((lo_n - lo_d, lo_d))
  cands.push((hi_n + hi_d, hi_d))
  for t = 0; t < inside.length(); t = t + 1 {
    cands.push(inside[t])
    if t + 1 < inside.length() {
      let (an, ad) = inside[t]
      let (bn, bd) = inside[t + 1]
      cands.push((an * bd + bn * ad, 2 * ad * bd))
    }
  }
  if inside.length() > 0 {
    let (an, ad) = inside[0]
    cands.push((lo_n * ad + an * lo_d, 2 * lo_d * ad))
    let (bn, bd) = inside[inside.length() - 1]
    cands.push((bn * hi_d + hi_n * bd, 2 * bd * hi_d))
  } else {
    cands.push((lo_n * 1 + hi_n * 1, 2 * lo_d * hi_d))
  }
  for t = 0; t < cands.length(); t = t + 1 {
    let (sn, sd) = cands[t]
    if sd <= 0 {
      continue
    }
    let wx = c * nx0 * sd - sn * ny0
    let wy = c * ny0 * sd + sn * nx0
    let den = d * sd
    let cell = argmin_cell(terms, wx, wy, den)
    if !(contains_cell(out, cell)) {
      out.push(cell)
    }
  }
  out
}

///| The argmin cell at the exact rational weight (nx/den, ny/den).
fn argmin_cell(
  terms : Array[(Array[Int], Int)],
  nx : Int,
  ny : Int,
  den : Int,
) -> Array[Int] {
  let mut best = 0
  let heights : Array[Int] = []
  for k = 0; k < terms.length(); k = k + 1 {
    let h = terms[k].1 * den + nx * terms[k].0[0] + ny * terms[k].0[1]
    heights.push(h)
    if k == 0 || h < best {
      best = h
    }
  }
  let out : Array[Int] = []
  for k = 0; k < terms.length(); k = k + 1 {
    if heights[k] == best {
      out.push(k)
    }
  }
  out
}

///| Dimension of the convex hull of a cell's exponent vectors: 0 for a point,
///| 1 for a collinear set, 2 for a genuine polygon.  A cell may contain lattice
///| points in the interior of an edge, so the number of indices is *not* the
///| dimension.
fn hull_dim(terms : Array[(Array[Int], Int)], c : Array[Int]) -> Int {
  if c.length() < 2 {
    return 0
  }
  let p0 = terms[c[0]].0
  let mut dx = 0
  let mut dy = 0
  let mut have = false
  for t = 1; t < c.length(); t = t + 1 {
    let p = terms[c[t]].0
    if p[0] != p0[0] || p[1] != p0[1] {
      if !(have) {
        dx = p[0] - p0[0]
        dy = p[1] - p0[1]
        have = true
      } else if cross(dx, dy, p[0] - p0[0], p[1] - p0[1]) != 0 {
        return 2
      }
    }
  }
  if have {
    1
  } else {
    0
  }
}

///| Maximal cells (hull dimension 2): the polygons of the Newton subdivision.
pub fn newton_faces(terms : Array[(Array[Int], Int)]) -> Array[Array[Int]] {
  let cells = newton_cells(terms)
  let out : Array[Array[Int]] = []
  for k = 0; k < cells.length(); k = k + 1 {
    if hull_dim(terms, cells[k]) == 2 {
      out.push(cells[k])
    }
  }
  out
}

///| Edge cells (hull dimension 1) of the Newton subdivision.  Interior lattice
///| points of an edge belong to the same cell, as they must.
pub fn newton_edges(terms : Array[(Array[Int], Int)]) -> Array[Array[Int]] {
  let cells = newton_cells(terms)
  let out : Array[Array[Int]] = []
  for k = 0; k < cells.length(); k = k + 1 {
    if hull_dim(terms, cells[k]) == 1 {
      out.push(cells[k])
    }
  }
  out
}

fn cross(ax : Int, ay : Int, bx : Int, by : Int) -> Int {
  ax * by - ay * bx
}

///| The two extreme points of a collinear cell, as a displacement vector.
fn extreme_displacement(
  terms : Array[(Array[Int], Int)],
  c : Array[Int],
) -> (Int, Int) {
  let mut bi = c[0]
  let mut bj = c[0]
  let mut best = -1
  for a = 0; a < c.length(); a = a + 1 {
    for b = a + 1; b < c.length(); b = b + 1 {
      let dx = terms[c[b]].0[0] - terms[c[a]].0[0]
      let dy = terms[c[b]].0[1] - terms[c[a]].0[1]
      let q = dx * dx + dy * dy
      if q > best {
        best = q
        bi = c[a]
        bj = c[b]
      }
    }
  }
  (terms[bj].0[0] - terms[bi].0[0], terms[bj].0[1] - terms[bi].0[1])
}

///| **Balancing condition.**  For every polygon C of the subdivision and every
///| **A51 Newton balancing**: for each
///| boundary edge cell E of C with displacement d, accumulate the rotated
///| displacement ±(d_y, -d_x) — this equals (lattice length of E) times the
///| primitive outward normal.  The cell is balanced when the sum vanishes.
pub fn newton_balancing(
  terms : Array[(Array[Int], Int)],
) -> Array[(Int, Int)] {
  let faces = newton_faces(terms)
  let edges = newton_edges(terms)
  let sums : Array[(Int, Int)] = []
  for f = 0; f < faces.length(); f = f + 1 {
    let c = faces[f]
    let mut sx = 0
    let mut sy = 0
    for e = 0; e < edges.length(); e = e + 1 {
      let ec = edges[e]
      // E must be a subset of C
      let mut subset = true
      for t = 0; t < ec.length(); t = t + 1 {
        let mut found = false
        for u = 0; u < c.length(); u = u + 1 {
          if c[u] == ec[t] {
            found = true
          }
        }
        if !(found) {
          subset = false
        }
      }
      if !(subset) || ec.length() == c.length() {
        continue
      }
      let (dx, dy) = extreme_displacement(terms, ec)
      // all points of C outside E must lie strictly on one side of the line E
      let mut sign = 0
      let mut boundary = true
      for t = 0; t < c.length(); t = t + 1 {
        let mut in_e = false
        for u = 0; u < ec.length(); u = u + 1 {
          if ec[u] == c[t] {
            in_e = true
          }
        }
        if in_e {
          continue
        }
        let o = terms[c[t]].0
        let base = terms[ec[0]].0
        let cr = cross(dx, dy, o[0] - base[0], o[1] - base[1])
        if cr == 0 {
          boundary = false
        } else if sign == 0 {
          sign = if cr > 0 {
            1
          } else {
            -1
          }
        } else if (cr > 0 && sign < 0) || (cr < 0 && sign > 0) {
          boundary = false
        }
      }
      if !(boundary) || sign == 0 {
        continue
      }
      // outward = away from the interior points of C
      let ox = if sign > 0 {
        dy
      } else {
        -dy
      }
      let oy = if sign > 0 {
        -dx
      } else {
        dx
      }
      sx = sx + ox
      sy = sy + oy
    }
    sums.push((sx, sy))
  }
  sums
}

///| Is the Newton subdivision of `terms` balanced (i.e. a tropical curve)?
pub fn is_balanced(terms : Array[(Array[Int], Int)]) -> Bool {
  let sums = newton_balancing(terms)
  for i = 0; i < sums.length(); i = i + 1 {
    if sums[i].0 != 0 || sums[i].1 != 0 {
      return false
    }
  }
  sums.length() > 0
}

// ---------------------------------------------------------------------------
//  The heat -> tropical bridge (Maslov dequantization)
// ---------------------------------------------------------------------------

///| Finite-h approximation of ⊕ on a family of real energies: the soft minimum.
pub fn trop_dequantize(xs : Array[Double], h : Double) -> Double {
  @heat.softmin_all(xs, h)
}

///| The h -> 0 limit of `trop_dequantize`: the tropical minimum, with argmin.
pub fn trop_limit(xs : Array[Double]) -> (Double, Int) {
  @heat.to_tropical(xs)
}
