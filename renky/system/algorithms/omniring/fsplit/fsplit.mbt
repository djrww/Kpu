// ===========================================================================
//  fsplit — perfect F-split rings (「完美的 F-分裂環」, structure 5).
//
//  Setting.  R = F_p[x_1..x_n]/I with dim_{F_p} R < infinity.  Let
//  B = (e_0 = 1, e_1, ...) be the basis of standard monomials (A35) and let
//  mu[a][j][m] be the structure constants, e_a * e_j = sum_m mu[a][j][m] e_m.
//
//  Definition (Frobenius).  F : R -> R, F(r) = r^p, is a ring homomorphism; its
//  matrix in the basis B has as columns the coordinates of e_i^p.
//
//  Definition (Frobenius splitting).  A *splitting* is an R-linear map
//  phi : F_* R -> R with phi o F = id_R.  Writing the R-action on F_* R as
//  r . s = r^p s, R-linearity is
//      phi(r^p s) = r phi(s)      for all r, s in R,
//  and the splitting condition is phi(r^p) = r for all r (in particular
//  phi(1) = 1).  R is *F-split* if such a phi exists.
//
//  Definition (perfect).  R is perfect if F is bijective; for a finite F_p
//  algebra this is equivalent to F being surjective, and to R being reduced.
//
//  Definition (compatibly split ideal).  J <= R is compatibly split by phi if
//  phi(J) <= J.
//
//  Theorems (verified exhaustively in `fsplit_wbtest`):
//   X1 the splitting conditions are an affine linear system in the entries of
//      the matrix of phi; it is solved here by Gaussian elimination (A27), so
//      "R is F-split" is decidable by linear algebra.
//   X2 F-split => reduced (F injective), and for a finite dimensional F_p
//      algebra the converse holds: R is F-split if and only if R is reduced,
//      if and only if R is perfect, if and only if rank(F) = dim R.
//   X3 the set of splittings is a torsor: if phi is one splitting then every
//      splitting is phi + psi with psi R-linear and psi o F = 0, and the space
//      of such psi has dimension = dim nullspace of the homogeneous system.
//   X4 for R = F_p^n (a product of n copies of F_p) the ideals compatibly split
//      by the identity are exactly the 2^n coordinate subspaces.
//   X5 F is a ring homomorphism and F(r) = r^p for every r.
// ===========================================================================

pub type Mat = @linalg.Mat
pub type MPoly = @poly.MPoly

///| Data of a finite dimensional quotient R = F_p[x]/I: basis, structure
///| constants, Frobenius matrix.
pub struct FiniteRing {
  p : Int
  nvars : Int
  basis : Array[Array[Int]]
  mul : Array[Array[Array[Int]]]
  frob : Mat
  st : @ideal.IdealState
}

fn mat_col(m : Mat, j : Int) -> Array[Int] {
  let out : Array[Int] = []
  for i = 0; i < @linalg.mat_rows(m); i = i + 1 {
    out.push(m[i][j])
  }
  out
}

fn elem(p : Int, nvars : Int, e : Array[Int]) -> MPoly {
  @poly.mp_make(p, nvars, [(e, 1)])
}

///| Coordinates of a polynomial in the standard-monomial basis.
fn coords(st : @ideal.IdealState, basis : Array[Array[Int]], f : MPoly) -> Array[Int] {
  let nf = st.normal_form(f)
  let out : Array[Int] = Array::make(basis.length(), 0)
  for i = 0; i < nf.terms.length(); i = i + 1 {
    let (e, c) = nf.terms[i]
    let mut found = false
    for k = 0; k < basis.length(); k = k + 1 {
      if exp_eq(basis[k], e) {
        out[k] = @kernel.modpos(out[k] + c, st.p)
        found = true
      }
    }
    if !(found) {
      abort("[fsplit] normal form is not spanned by the basis")
    }
  }
  out
}

fn exp_eq(a : Array[Int], b : Array[Int]) -> Bool {
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

///| **A79 structure constants and Frobenius matrix**: build the finite ring data
///| from an ideal state whose quotient is finite.
pub fn finite_ring(st : @ideal.IdealState) -> FiniteRing {
  if !(st.finite) {
    abort("[fsplit] the quotient must be finite dimensional")
  }
  let p = st.p
  let nvars = st.nvars
  let basis = st.quotient_basis()
  let dim = basis.length()
  // build the structure constants with fresh rows (no aliasing)
  let mu : Array[Array[Array[Int]]] = []
  for a = 0; a < dim; a = a + 1 {
    let row : Array[Array[Int]] = []
    for j = 0; j < dim; j = j + 1 {
      let prod = elem(p, nvars, basis[a]).mul(elem(p, nvars, basis[j]))
      row.push(coords(st, basis, prod))
    }
    mu.push(row)
  }
  // mu[a][j][m] : coefficient of e_m in e_a * e_j
  let mul3 : Array[Array[Array[Int]]] = []
  for a = 0; a < dim; a = a + 1 {
    let row : Array[Array[Int]] = []
    for j = 0; j < dim; j = j + 1 {
      row.push(mu[a][j])
    }
    mul3.push(row)
  }
  let fr = @linalg.mat_make(dim, dim)
  for j = 0; j < dim; j = j + 1 {
    let mut pw = elem(p, nvars, basis[j])
    for _k = 1; _k < p; _k = _k + 1 {
      pw = pw.mul(elem(p, nvars, basis[j]))
    }
    let cs = coords(st, basis, pw)
    for i = 0; i < dim; i = i + 1 {
      fr[i][j] = cs[i]
    }
  }
  { p, nvars, basis, mul: mul3, frob: fr, st }
}

pub fn FiniteRing::dim(self : FiniteRing) -> Int {
  self.basis.length()
}

///| Element from coordinates.
pub fn FiniteRing::of_coords(self : FiniteRing, cs : Array[Int]) -> MPoly {
  let terms : Array[(Array[Int], Int)] = []
  for i = 0; i < self.dim(); i = i + 1 {
    if i < cs.length() && cs[i] != 0 {
      terms.push((self.basis[i], cs[i]))
    }
  }
  @poly.mp_make(self.p, self.nvars, terms)
}

///| Sum of two coordinate vectors.
pub fn FiniteRing::add_vec(self : FiniteRing, a : Array[Int], b : Array[Int]) -> Array[Int] {
  let dim = self.dim()
  let out : Array[Int] = Array::make(dim, 0)
  for i = 0; i < dim; i = i + 1 {
    let x = if i < a.length() {
      a[i]
    } else {
      0
    }
    let y = if i < b.length() {
      b[i]
    } else {
      0
    }
    out[i] = @kernel.addmod(x, y, self.p)
  }
  out
}

///| Product of two coordinate vectors, using the structure constants.
pub fn FiniteRing::mul_coords(self : FiniteRing, a : Array[Int], b : Array[Int]) -> Array[Int] {
  let dim = self.dim()
  let out : Array[Int] = Array::make(dim, 0)
  for i = 0; i < dim; i = i + 1 {
    for j = 0; j < dim; j = j + 1 {
      for m = 0; m < dim; m = m + 1 {
        out[m] = @kernel.addmod(out[m], a[i] * b[j] * self.mul[i][j][m], self.p)
      }
    }
  }
  out
}

///| p-th power map on coordinates (the Frobenius).
pub fn FiniteRing::frob_coords(self : FiniteRing, a : Array[Int]) -> Array[Int] {
  @linalg.mat_vec(self.frob, a, self.p)
}

///| All elements of R as coordinate vectors.
pub fn FiniteRing::elements(self : FiniteRing) -> Array[Array[Int]] {
  let dim = self.dim()
  let total = @arith.pow_int(self.p, dim)
  let out : Array[Array[Int]] = []
  for k = 0; k < total; k = k + 1 {
    let ds = @arith.digits(k, self.p)
    let v : Array[Int] = Array::make(dim, 0)
    for i = 0; i < dim; i = i + 1 {
      v[i] = if i < ds.length() {
        ds[i]
      } else {
        0
      }
    }
    out.push(v)
  }
  out
}

///| Is R reduced?  Equivalent to I = rad(I), tested by comparing the quotient
///| dimensions of I and its radical.
pub fn FiniteRing::is_reduced(self : FiniteRing) -> Bool {
  let rad = self.st.radical()
  rad.dim == self.dim()
}

///| **A81 reduced / perfect tests**.  Is R reduced, tested directly: no nonzero nilpotent, i.e. x^p = 0 => x = 0.
pub fn FiniteRing::is_reduced_bruteforce(self : FiniteRing) -> Bool {
  let els = self.elements()
  for i = 0; i < els.length(); i = i + 1 {
    let x = els[i]
    let mut zero = true
    for k = 0; k < self.dim(); k = k + 1 {
      if k < x.length() && x[k] != 0 {
        zero = false
      }
    }
    if zero {
      continue
    }
    if is_zero_vec(self.frob_coords(x)) {
      return false
    }
  }
  true
}

fn is_zero_vec(v : Array[Int]) -> Bool {
  for i = 0; i < v.length(); i = i + 1 {
    if v[i] != 0 {
      return false
    }
  }
  true
}

///| Is R perfect, i.e. is the Frobenius bijective?
pub fn FiniteRing::is_perfect(self : FiniteRing) -> Bool {
  @linalg.rank(self.frob, self.p) == self.dim()
}

///| **A80 affine splitting system**: solve the affine system for a Frobenius splitting.  Returns
///| (particular solution or None, dimension of the space of psi with
///| psi o F = 0 and psi R-linear).
pub fn FiniteRing::splitting_system(self : FiniteRing) -> (Array[Int]?, Int) {
  let dim = self.dim()
  let p = self.p
  let nu = dim * dim
  // rows of the augmented system
  let rows : Array[Array[Int]] = []
  // (A) phi(e_r^p) = e_r
  for r = 0; r < dim; r = r + 1 {
    let fcs = mat_col(self.frob, r)
    for i = 0; i < dim; i = i + 1 {
      let row : Array[Int] = Array::make(nu + 1, 0)
      for k = 0; k < dim; k = k + 1 {
        // sum_k (e_r^p)_k * M[i][k] = delta_{i,r}
        row[i * dim + k] = fcs[k]
      }
      row[nu] = if i == r {
        1
      } else {
        0
      }
      rows.push(row)
    }
  }
  // (B) phi(e_a^p * e_b) = e_a * phi(e_b)
  for a = 0; a < dim; a = a + 1 {
    let fcs = mat_col(self.frob, a)
    for b = 0; b < dim; b = b + 1 {
      // e_a^p * e_b
      let lhs_coords = self.mul_coords(fcs, unit_vec(dim, b))
      for i = 0; i < dim; i = i + 1 {
        let row : Array[Int] = Array::make(nu + 1, 0)
        for k = 0; k < dim; k = k + 1 {
          row[i * dim + k] = lhs_coords[k]
        }
        // minus sum_j mu[a][j][i] * M[j][b]
        for j = 0; j < dim; j = j + 1 {
          row[j * dim + b] = @kernel.modpos(row[j * dim + b] - self.mul[a][j][i], p)
        }
        rows.push(row)
      }
    }
  }
  // Gaussian elimination on the augmented matrix
  let (sol, consistent, nulldim) = solve_affine(rows, nu, p)
  if consistent {
    (Some(sol), nulldim)
  } else {
    (None, nulldim)
  }
}

fn unit_vec(n : Int, i : Int) -> Array[Int] {
  let v : Array[Int] = Array::make(n, 0)
  v[i] = 1
  v
}

///| Solve A x = b over F_p; returns (a particular solution, consistency,
///| dimension of the solution space of A x = 0).
fn solve_affine(rows : Array[Array[Int]], nu : Int, p : Int) -> (Array[Int], Bool, Int) {
  let a = rows.copy()
  let nr = a.length()
  let pivots : Array[Int] = []
  let mut row = 0
  for col = 0; col < nu; col = col + 1 {
    let mut sel = -1
    for r = row; r < nr; r = r + 1 {
      if a[r][col] != 0 {
        sel = r
      }
    }
    if sel < 0 {
      continue
    }
    let tmp = a[row]
    a[row] = a[sel]
    a[sel] = tmp
    let inv = @arith.invmod(a[row][col], p)
    for c = 0; c <= nu; c = c + 1 {
      a[row][c] = @kernel.modpos(a[row][c] * inv, p)
    }
    for r = 0; r < nr; r = r + 1 {
      if r == row || a[r][col] == 0 {
        continue
      }
      let f = a[r][col]
      for c = 0; c <= nu; c = c + 1 {
        a[r][c] = @kernel.modpos(a[r][c] - f * a[row][c], p)
      }
    }
    pivots.push(col)
    row = row + 1
  }
  // consistency: no row of the form 0 ... 0 | c with c != 0
  for r = row; r < nr; r = r + 1 {
    if a[r][nu] != 0 {
      return ([], false, nu - pivots.length())
    }
  }
  // particular solution: free variables = 0
  let sol : Array[Int] = Array::make(nu, 0)
  for i = 0; i < pivots.length(); i = i + 1 {
    sol[pivots[i]] = a[i][nu]
  }
  (sol, true, nu - pivots.length())
}

///| A splitting phi, as a matrix (columns = images of the basis), or None.
pub fn FiniteRing::splitting(self : FiniteRing) -> Mat? {
  let (sol, _) = self.splitting_system()
  match sol {
    None => None
    Some(xs) => {
      let dim = self.dim()
      let m = @linalg.mat_make(dim, dim)
      for i = 0; i < dim; i = i + 1 {
        for j = 0; j < dim; j = j + 1 {
          m[i][j] = xs[i * dim + j]
        }
      }
      Some(m)
    }
  }
}

///| Is R F-split?
pub fn FiniteRing::is_f_split(self : FiniteRing) -> Bool {
  match self.splitting() {
    Some(_) => true
    None => false
  }
}

///| Apply a splitting matrix to coordinates.
pub fn apply_phi(phi : Mat, v : Array[Int], p : Int) -> Array[Int] {
  @linalg.mat_vec(phi, v, p)
}

///| Is the ideal spanned by `gens` (as coordinate vectors) invariant under phi?
pub fn compatibly_split(phi : Mat, gens : Array[Array[Int]], p : Int) -> Bool {
  for i = 0; i < gens.length(); i = i + 1 {
    let img = apply_phi(phi, gens[i], p)
    // img must lie in the span of gens
    if !(in_span(img, gens, p)) {
      return false
    }
  }
  true
}

///| Membership in the F_p-span of a list of vectors.
pub fn in_span(v : Array[Int], gens : Array[Array[Int]], p : Int) -> Bool {
  let m : Mat = []
  for i = 0; i < gens.length(); i = i + 1 {
    m.push(gens[i].copy())
  }
  let r0 = @linalg.rank(m, p)
  m.push(v.copy())
  let r1 = @linalg.rank(m, p)
  r0 == r1
}
