// ===========================================================================
//  cohom — cohomology rings (「上同調環」, structure 2): graded-commutative
//  algebras with Poincaré duality, plus the algebraic de Rham DG algebra of a
//  jet algebra, whose cohomology is a graded ring.
//
//  Definition (cohomology ring).  A finite dimensional graded-commutative
//  F_p-algebra H^* = (+)_k H^k with a cup product
//      cup : H^i x H^j -> H^{i+j},   a cup b = (-1)^{ij} b cup a,
//  given by structure constants cup[i][j][k] on the basis of each degree, and a
//  distinguished top degree n with dim H^n = 1 (the fundamental class).
//
//  Definition (Poincaré duality).  H^* satisfies Poincaré duality in degree n if
//  for every k the pairing  H^k x H^{n-k} -> H^n ~= F_p  is non-degenerate,
//  i.e. its matrix has rank dim H^k.
//
//  Definition (DG algebra).  A graded algebra A^* with a differential
//  d : A^k -> A^{k+1}, d^2 = 0, satisfying the graded Leibniz rule
//      d(a * b) = d(a) * b + (-1)^{|a|} a * d(b).
//  Its cohomology H^*(A) = ker d / im d inherits a graded-commutative product.
//
//  The concrete DG algebra used here is the algebraic de Rham complex of the
//  jet algebra J = F_p[t]/(t^n):
//      Omega^0 = J,  Omega^1 = J dt,  d(f + g dt) = f' dt,
//      (f + g dt) ^ (f' + g' dt) = f f' + (f g' + g f') dt   (dt ^ dt = 0).
//
//  Theorems (verified exhaustively in `cohom_wbtest`):
//   C1 the models H^*(S^n), H^*(T^2), H^*(CP^2), H^*(RP^2;F_2),
//      H^*(S^2 x S^2) are associative, unital, graded-commutative algebras.
//   C2 Poincaré duality holds for all of them: the pairing matrix
//      H^k x H^{n-k} -> H^n has rank dim H^k.
//   C3 chi = sum_k (-1)^k dim H^k agrees with the topological values
//      (2 for S^n, 0 for T^2, 3 for CP^2, 1 for RP^2, 4 for S^2 x S^2).
//   C4 the de Rham complex of J: d^2 = 0 and the graded Leibniz rule hold for
//      every pair of forms.
//   C5 the product of cocycles is a cocycle and the product of a coboundary with
//      anything is a coboundary, so H^*(Omega) is a graded ring.
//   C6 dim H^0_dR = #{i < n : p | i} and dim H^1_dR = 1 + #{1 <= i < n : p | i},
//      both computed independently by linear algebra (rank of the derivative
//      map), and chi(Omega) = chi(H).
// ===========================================================================

pub type Mat = @linalg.Mat

///| **A77 cup product**: a finite dimensional graded-commutative F_p-algebra with a top degree.
pub struct CohomRing {
  name : String
  p : Int
  dims : Array[Int]
  cup : Array[Array[Array[Int]]]
  top : Int
}

///| A fresh n x n x n array of zeros.  (`Array::make(n, x)` would alias the same
///| inner array n times, so nested tables must be built with explicit loops.)
fn zeros3(n : Int) -> Array[Array[Array[Int]]] {
  let out : Array[Array[Array[Int]]] = []
  for _i = 0; _i < n; _i = _i + 1 {
    let row : Array[Array[Int]] = []
    for _j = 0; _j < n; _j = _j + 1 {
      row.push(Array::make(n, 0))
    }
    out.push(row)
  }
  out
}

fn ring_mk(
  name : String,
  p : Int,
  dims : Array[Int],
  cup : Array[Array[Array[Int]]],
  top : Int,
) -> CohomRing {
  { name, p, dims, cup, top }
}

///| Total number of basis elements.
pub fn CohomRing::size(self : CohomRing) -> Int {
  let mut acc = 0
  for i = 0; i < self.dims.length(); i = i + 1 {
    acc = acc + self.dims[i]
  }
  acc
}

///| Degree of the basis element with index i.
pub fn CohomRing::deg(self : CohomRing, i : Int) -> Int {
  let mut acc = 0
  for d = 0; d < self.dims.length(); d = d + 1 {
    if i < acc + self.dims[d] {
      return d
    }
    acc = acc + self.dims[d]
  }
  -1
}

///| The index range of the basis elements of degree d.
pub fn CohomRing::basis_range(self : CohomRing, d : Int) -> (Int, Int) {
  let mut acc = 0
  for k = 0; k < d; k = k + 1 {
    acc = acc + self.dims[k]
  }
  (acc, acc + self.dims[d])
}

pub fn CohomRing::euler(self : CohomRing) -> Int {
  let mut acc = 0
  for k = 0; k < self.dims.length(); k = k + 1 {
    acc = if k % 2 == 0 {
      acc + self.dims[k]
    } else {
      acc - self.dims[k]
    }
  }
  acc
}

pub fn CohomRing::to_string(self : CohomRing) -> String {
  self.name + " over F_" + self.p.to_string() + " dims " + self.dims.to_string()
}

///| H^*(S^n) = F_p[x]/(x^2), |x| = n.
pub fn sphere(p : Int, n : Int) -> CohomRing {
  let dims = Array::make(n + 1, 0)
  dims[0] = 1
  dims[n] = 1
  let size = 2
  let cup = zeros3(size)
  cup[0][0][0] = 1
  cup[0][1][1] = 1
  cup[1][0][1] = 1
  ring_mk("S^" + n.to_string(), p, dims, cup, n)
}

///| H^*(T^2) = Lambda(a, b) with |a| = |b| = 1: basis 1, a, b, ab.
pub fn torus2(p : Int) -> CohomRing {
  let dims = [1, 2, 1]
  let cup = zeros3(4)
  // indices: 0 = 1, 1 = a, 2 = b, 3 = ab
  for i = 0; i < 4; i = i + 1 {
    cup[0][i][i] = 1
    cup[i][0][i] = 1
  }
  cup[1][2][3] = 1
  cup[2][1][3] = @kernel.modpos(-1, p)
  ring_mk("T^2", p, dims, cup, 2)
}

///| H^*(CP^2) = F_p[h]/(h^3), |h| = 2: basis 1, h, h^2 in degrees 0, 2, 4.
pub fn cp2(p : Int) -> CohomRing {
  let dims = [1, 0, 1, 0, 1]
  let cup = zeros3(3)
  for i = 0; i < 3; i = i + 1 {
    for j = 0; j + i < 3; j = j + 1 {
      cup[i][j][i + j] = 1
    }
  }
  ring_mk("CP^2", p, dims, cup, 4)
}

///| H^*(RP^2; F_2) = F_2[x]/(x^3), |x| = 1: basis 1, x, x^2.
pub fn rp2(p : Int) -> CohomRing {
  let dims = [1, 1, 1]
  let cup = zeros3(3)
  for i = 0; i < 3; i = i + 1 {
    for j = 0; j + i < 3; j = j + 1 {
      cup[i][j][i + j] = 1
    }
  }
  ring_mk("RP^2", p, dims, cup, 2)
}

///| H^*(S^2 x S^2) = F_p[a,b]/(a^2, b^2), |a| = |b| = 2: basis 1, a, b, ab.
pub fn s2xs2(p : Int) -> CohomRing {
  let dims = [1, 0, 2, 0, 1]
  let cup = zeros3(4)
  for i = 0; i < 4; i = i + 1 {
    cup[0][i][i] = 1
    cup[i][0][i] = 1
  }
  cup[1][2][3] = 1
  cup[2][1][3] = 1
  ring_mk("S^2 x S^2", p, dims, cup, 4)
}

///| Associativity of the cup product on all basis triples.
pub fn CohomRing::associativity_ok(self : CohomRing) -> Bool {
  let n = self.cup.length()
  for i = 0; i < n; i = i + 1 {
    for j = 0; j < n; j = j + 1 {
      for t = 0; t < n; t = t + 1 {
        for m = 0; m < n; m = m + 1 {
          let mut lhs = 0
          let mut rhs = 0
          for k = 0; k < n; k = k + 1 {
            lhs = @kernel.addmod(lhs, self.cup[i][j][k] * self.cup[k][t][m], self.p)
            rhs = @kernel.addmod(rhs, self.cup[j][t][k] * self.cup[i][k][m], self.p)
          }
          if lhs != rhs {
            return false
          }
        }
      }
    }
  }
  true
}

///| Graded commutativity: cup[i][j][k] = (-1)^{|i||j|} cup[j][i][k].
pub fn CohomRing::graded_commutative_ok(self : CohomRing) -> Bool {
  let n = self.cup.length()
  for i = 0; i < n; i = i + 1 {
    let di = self.deg(i)
    for j = 0; j < n; j = j + 1 {
      let dj = self.deg(j)
      let sign = if (di * dj) % 2 == 0 {
        1
      } else {
        self.p - 1
      }
      for k = 0; k < n; k = k + 1 {
        if self.cup[i][j][k] != @kernel.modpos(sign * self.cup[j][i][k], self.p) {
          return false
        }
      }
    }
  }
  true
}

///| Degree preservation: cup[i][j][k] != 0 implies |k| = |i| + |j|.
pub fn CohomRing::degree_ok(self : CohomRing) -> Bool {
  let n = self.cup.length()
  for i = 0; i < n; i = i + 1 {
    for j = 0; j < n; j = j + 1 {
      for k = 0; k < n; k = k + 1 {
        if self.cup[i][j][k] != 0 && self.deg(k) != self.deg(i) + self.deg(j) {
          return false
        }
      }
    }
  }
  true
}

///| **A77 Poincaré pairing**: the pairing matrix H^k x H^{n-k} -> H^n ~= F_p.
pub fn CohomRing::pairing_matrix(self : CohomRing, k : Int) -> Mat {
  let (a0, a1) = self.basis_range(k)
  let (b0, b1) = self.basis_range(self.top - k)
  let (t0, _) = self.basis_range(self.top)
  let m = @linalg.mat_make(a1 - a0, b1 - b0)
  for i = a0; i < a1; i = i + 1 {
    for j = b0; j < b1; j = j + 1 {
      m[i - a0][j - b0] = self.cup[i][j][t0]
    }
  }
  m
}

///| **C2 Poincaré duality**: every pairing H^k x H^{n-k} -> H^n is perfect.
pub fn CohomRing::poincare_duality_ok(self : CohomRing) -> Bool {
  if self.dims[self.top] != 1 {
    return false
  }
  for k = 0; k <= self.top; k = k + 1 {
    let m = self.pairing_matrix(k)
    if @linalg.rank(m, self.p) != self.dims[k] {
      return false
    }
    if self.dims[k] != self.dims[self.top - k] {
      return false
    }
  }
  true
}

// ---------------------------------------------------------------------------
//  C4–C6  the algebraic de Rham DG algebra of the jet algebra F_p[t]/(t^n)
// ---------------------------------------------------------------------------

///| **A78 de Rham DG algebra of a jet algebra**: a differential form f + g*dt on
///| Spec F_p[t]/(t^n), with wedge product, exterior derivative and cohomology.
pub struct DRForm {
  p : Int
  n : Int
  f : Array[Int]
  g : Array[Int]
}

fn dr_mk(p : Int, n : Int, f : Array[Int], g : Array[Int]) -> DRForm {
  let ff : Array[Int] = Array::make(n, 0)
  let gg : Array[Int] = Array::make(n, 0)
  for i = 0; i < n; i = i + 1 {
    ff[i] = if i < f.length() {
      @kernel.modpos(f[i], p)
    } else {
      0
    }
    gg[i] = if i < g.length() {
      @kernel.modpos(g[i], p)
    } else {
      0
    }
  }
  { p, n, f: ff, g: gg }
}

pub fn dr_form(p : Int, n : Int, f : Array[Int], g : Array[Int]) -> DRForm {
  dr_mk(p, n, f, g)
}

pub fn DRForm::add(self : DRForm, other : DRForm) -> DRForm {
  let f : Array[Int] = []
  let g : Array[Int] = []
  for i = 0; i < self.n; i = i + 1 {
    f.push(self.f[i] + other.f[i])
    g.push(self.g[i] + other.g[i])
  }
  dr_mk(self.p, self.n, f, g)
}

pub fn DRForm::neg(self : DRForm) -> DRForm {
  let f : Array[Int] = []
  let g : Array[Int] = []
  for i = 0; i < self.n; i = i + 1 {
    f.push(-self.f[i])
    g.push(-self.g[i])
  }
  dr_mk(self.p, self.n, f, g)
}

pub fn DRForm::sub(self : DRForm, other : DRForm) -> DRForm {
  self.add(other.neg())
}

pub fn DRForm::is_zero(self : DRForm) -> Bool {
  for i = 0; i < self.n; i = i + 1 {
    if self.f[i] != 0 || self.g[i] != 0 {
      return false
    }
  }
  true
}

pub fn DRForm::to_string(self : DRForm) -> String {
  "(" + self.f.to_string() + ") + (" + self.g.to_string() + ")dt"
}

///| Truncated polynomial product in J = F_p[t]/(t^n).
fn jmul(a : Array[Int], b : Array[Int], p : Int, n : Int) -> Array[Int] {
  let out : Array[Int] = Array::make(n, 0)
  for i = 0; i < n; i = i + 1 {
    for j = 0; j + i < n; j = j + 1 {
      out[i + j] = @kernel.addmod(out[i + j], a[i] * b[j], p)
    }
  }
  out
}

///| Wedge product: (f + g dt) ^ (f' + g' dt) = f f' + (f g' + g f') dt.
pub fn DRForm::wedge(self : DRForm, other : DRForm) -> DRForm {
  let f0 = jmul(self.f, other.f, self.p, self.n)
  let g0 = jmul(self.f, other.g, self.p, self.n)
  let g1 = jmul(self.g, other.f, self.p, self.n)
  let gg : Array[Int] = []
  for i = 0; i < self.n; i = i + 1 {
    gg.push(g0[i] + g1[i])
  }
  dr_mk(self.p, self.n, f0, gg)
}

///| d(f + g dt) = f' dt.
pub fn DRForm::deriv(self : DRForm) -> DRForm {
  let gg : Array[Int] = []
  for i = 1; i < self.n; i = i + 1 {
    gg.push(@kernel.modpos(i * self.f[i], self.p))
  }
  dr_mk(self.p, self.n, [], gg)
}

///| The degree of a nonzero form: 0 if g = 0, else 1.
pub fn DRForm::deg(self : DRForm) -> Int {
  for i = 0; i < self.n; i = i + 1 {
    if self.g[i] != 0 {
      return 1
    }
  }
  0
}

///| All 2n-dimensional space of forms: p^(2n) of them (use small p, n).
pub fn all_forms(p : Int, n : Int) -> Array[DRForm] {
  let total = @arith.pow_int(p, 2 * n)
  let out : Array[DRForm] = []
  for k = 0; k < total; k = k + 1 {
    let ds = @arith.digits(k, p)
    let f : Array[Int] = []
    let g : Array[Int] = []
    for i = 0; i < n; i = i + 1 {
      f.push(if i < ds.length() {
        ds[i]
      } else {
        0
      })
      g.push(if n + i < ds.length() {
        ds[n + i]
      } else {
        0
      })
    }
    out.push(dr_mk(p, n, f, g))
  }
  out
}

///| **A78 (cont.)/C6 de Rham cohomology of the jet algebra**, by linear algebra:
///| H^0 = ker(d) and H^1 = coker(d) for d : J -> J, f |-> f'.
pub fn dr_cohomology_dims(p : Int, n : Int) -> (Int, Int) {
  let d = @linalg.mat_make(n, n)
  for i = 1; i < n; i = i + 1 {
    d[i - 1][i] = @kernel.modpos(i, p)
  }
  let r = @linalg.rank(d, p)
  (n - r, n - r)
}

///| The same dimensions by counting: dim H^0 = #{i < n : p | i} and
///| dim H^1 = 1 + #{1 <= i < n : p | i}.  Both count the p-th-power classes,
///| which are exactly the classes killed by the derivative in characteristic p.
pub fn dr_cohomology_count(p : Int, n : Int) -> (Int, Int) {
  let mut h0 = 0
  let mut h1 = 1
  for i = 0; i < n; i = i + 1 {
    if i % p == 0 {
      h0 = h0 + 1
    }
    if i >= 1 && i % p == 0 {
      h1 = h1 + 1
    }
  }
  (h0, h1)
}

///| Does d/dt descend to the jet algebra J = F_p[t]/(t^n)?  It does exactly
///| when d(t^n) = n t^{n-1} vanishes in J, i.e. when p | n.  This is the
///| criterion for the de Rham complex of J to be well defined.
pub fn dr_well_defined(p : Int, n : Int) -> Bool {
  // n * t^{n-1} = 0 in F_p[t]/(t^n)  <=>  p | n  (t^{n-1} != 0 in J)
  @kernel.modpos(n, p) == 0
}
