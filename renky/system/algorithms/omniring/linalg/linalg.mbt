// ===========================================================================
//  Linear algebra: the engine behind homology, cohomology, F-splitting,
//  Berlekamp factorisation, ideal quotients and the tropical determinant.
//
//  Two coefficient worlds:
//    * `F_p`  (a field)  -> Gaussian elimination, rank, kernel, image, det.
//    * `Z`    (a PID)    -> Hermite and Smith normal forms, i.e. the
//                          structure theorem for finitely generated modules
//                          over a PID, which is what integral homology is.
// ===========================================================================

///| A matrix as an array of rows; every row has the same length.
pub type Mat = Array[Array[Int]]

pub fn mat_make(rows : Int, cols : Int) -> Mat {
  let out : Mat = []
  for _i = 0; _i < rows; _i = _i + 1 {
    out.push(Array::make(cols, 0))
  }
  out
}

pub fn mat_identity(n : Int) -> Mat {
  let m = mat_make(n, n)
  for i = 0; i < n; i = i + 1 {
    m[i][i] = 1
  }
  m
}

pub fn mat_rows(m : Mat) -> Int {
  m.length()
}

pub fn mat_cols(m : Mat) -> Int {
  if m.is_empty() {
    0
  } else {
    m[0].length()
  }
}

pub fn mat_copy(m : Mat) -> Mat {
  m.map(fn(r) { r.copy() })
}

pub fn mat_get(m : Mat, i : Int, j : Int) -> Int {
  m[i][j]
}

pub fn mat_set(m : Mat, i : Int, j : Int, v : Int) -> Unit {
  m[i][j] = v
}

pub fn mat_transpose(m : Mat) -> Mat {
  let out = mat_make(mat_cols(m), mat_rows(m))
  for i = 0; i < mat_rows(m); i = i + 1 {
    for j = 0; j < mat_cols(m); j = j + 1 {
      out[j][i] = m[i][j]
    }
  }
  out
}

pub fn mat_add(a : Mat, b : Mat, p : Int) -> Mat {
  let out = mat_make(mat_rows(a), mat_cols(a))
  for i = 0; i < mat_rows(a); i = i + 1 {
    for j = 0; j < mat_cols(a); j = j + 1 {
      out[i][j] = @kernel.addmod(a[i][j], b[i][j], p)
    }
  }
  out
}

pub fn mat_scale(a : Mat, c : Int, p : Int) -> Mat {
  let out = mat_make(mat_rows(a), mat_cols(a))
  for i = 0; i < mat_rows(a); i = i + 1 {
    for j = 0; j < mat_cols(a); j = j + 1 {
      out[i][j] = @kernel.mulmod(a[i][j], c, p)
    }
  }
  out
}

///| Matrix product over F_p.
pub fn mat_mul(a : Mat, b : Mat, p : Int) -> Mat {
  let n = mat_rows(a)
  let k = mat_cols(a)
  let m = mat_cols(b)
  let out = mat_make(n, m)
  for i = 0; i < n; i = i + 1 {
    for j = 0; j < m; j = j + 1 {
      let mut acc = 0
      for t = 0; t < k; t = t + 1 {
        acc = @kernel.addmod(acc, @kernel.mulmod(a[i][t], b[t][j], p), p)
      }
      out[i][j] = acc
    }
  }
  out
}

///| Matrix-vector product over F_p.
pub fn mat_vec(a : Mat, v : Array[Int], p : Int) -> Array[Int] {
  let n = mat_rows(a)
  let out : Array[Int] = Array::make(n, 0)
  for i = 0; i < n; i = i + 1 {
    let mut acc = 0
    for j = 0; j < mat_cols(a); j = j + 1 {
      acc = @kernel.addmod(acc, @kernel.mulmod(a[i][j], v[j], p), p)
    }
    out[i] = acc
  }
  out
}

///| **A27 reduced row echelon form over F_p.**  Returns the reduced matrix and
///| the list of pivot columns.
pub fn rref(m : Mat, p : Int) -> (Mat, Array[Int]) {
  let a = mat_copy(m)
  let cols = mat_cols(a)
  let rows = mat_rows(a)
  let pivots : Array[Int] = []
  let mut r = 0
  for c = 0; c < cols && r < rows; c = c + 1 {
    // find a pivot
    let mut piv = -1
    for i = r; i < rows; i = i + 1 {
      if a[i][c] != 0 && piv < 0 {
        piv = i
      }
    }
    if piv < 0 {
      continue
    }
    // swap
    if piv != r {
      let tmp = a[r]
      a[r] = a[piv]
      a[piv] = tmp
    }
    // normalise
    let inv = @arith.invmod(a[r][c], p)
    for j = 0; j < cols; j = j + 1 {
      a[r][j] = @kernel.mulmod(a[r][j], inv, p)
    }
    // eliminate
    for i = 0; i < rows; i = i + 1 {
      if i == r || a[i][c] == 0 {
        continue
      }
      let f = a[i][c]
      for j = 0; j < cols; j = j + 1 {
        a[i][j] = @kernel.submod(a[i][j], @kernel.mulmod(f, a[r][j], p), p)
      }
    }
    pivots.push(c)
    r = r + 1
  }
  (a, pivots)
}

///| Rank over F_p.
pub fn rank(m : Mat, p : Int) -> Int {
  rref(m, p).1.length()
}

///| **A28 kernel (null space) basis over F_p**: all `v` with `m*v == 0`.
pub fn nullspace(m : Mat, p : Int) -> Array[Array[Int]] {
  let (a, pivots) = rref(m, p)
  let cols = mat_cols(a)
  let is_pivot : Array[Bool] = Array::make(cols, false)
  for i = 0; i < pivots.length(); i = i + 1 {
    is_pivot[pivots[i]] = true
  }
  let basis : Array[Array[Int]] = []
  for c = 0; c < cols; c = c + 1 {
    if is_pivot[c] {
      continue
    }
    let v : Array[Int] = Array::make(cols, 0)
    v[c] = 1
    for r = 0; r < pivots.length(); r = r + 1 {
      if a[r][c] != 0 {
        v[pivots[r]] = @kernel.negmod(a[r][c], p)
      }
    }
    basis.push(v)
  }
  basis
}

///| Basis of the row space (image of the transpose) over F_p.
pub fn rowspace(m : Mat, p : Int) -> Array[Array[Int]] {
  let (a, pivots) = rref(m, p)
  let out : Array[Array[Int]] = []
  for i = 0; i < pivots.length(); i = i + 1 {
    out.push(a[i].copy())
  }
  out
}

///| Determinant over F_p by elimination.
pub fn det_mod(m : Mat, p : Int) -> Int {
  let a = mat_copy(m)
  let n = mat_rows(a)
  if n != mat_cols(a) {
    abort("det_mod: not square")
  }
  let mut sign = 1
  let mut acc = 1
  for c = 0; c < n; c = c + 1 {
    let mut piv = -1
    for i = c; i < n; i = i + 1 {
      if a[i][c] != 0 && piv < 0 {
        piv = i
      }
    }
    if piv < 0 {
      return 0
    }
    if piv != c {
      let tmp = a[c]
      a[c] = a[piv]
      a[piv] = tmp
      sign = -sign
    }
    let inv = @arith.invmod(a[c][c], p)
    acc = @kernel.mulmod(acc, a[c][c], p)
    for i = c + 1; i < n; i = i + 1 {
      if a[i][c] == 0 {
        continue
      }
      let f = @kernel.mulmod(a[i][c], inv, p)
      for j = c; j < n; j = j + 1 {
        a[i][j] = @kernel.submod(a[i][j], @kernel.mulmod(f, a[c][j], p), p)
      }
    }
  }
  if sign < 0 {
    @kernel.negmod(acc, p)
  } else {
    acc
  }
}

///| Solve `m * x = b` over F_p; returns `None` when inconsistent, otherwise one
///| particular solution together with a kernel basis.
pub fn solve_mod(m : Mat, b : Array[Int], p : Int) -> (Array[Int], Array[Array[Int]])? {
  let aug : Mat = []
  for i = 0; i < mat_rows(m); i = i + 1 {
    let row = m[i].copy()
    row.push(@kernel.modpos(b[i], p))
    aug.push(row)
  }
  let (a, pivots) = rref(aug, p)
  let cols = mat_cols(m)
  for i = 0; i < pivots.length(); i = i + 1 {
    if pivots[i] == cols {
      return None
    }
  }
  let x : Array[Int] = Array::make(cols, 0)
  for i = 0; i < pivots.length(); i = i + 1 {
    x[pivots[i]] = a[i][cols]
  }
  Some((x, nullspace(m, p)))
}

// ---------------------------------------------------------------------------
//  Integer matrices: Hermite and Smith normal forms
// ---------------------------------------------------------------------------

///| **A29 Hermite normal form** of an integer matrix (row style, upper
///| triangular with non-negative pivots).  Used for lattices: the group
///| completion of an affine semigroup and its saturation.
pub fn hnf(m : Mat) -> Mat {
  let a = mat_copy(m)
  let rows = mat_rows(a)
  let cols = mat_cols(a)
  let mut r = 0
  for c = 0; c < cols && r < rows; c = c + 1 {
    // pick the row with the smallest nonzero |entry| in column c, then reduce
    let mut done = false
    while !(done) {
      let mut piv = -1
      for i = r; i < rows; i = i + 1 {
        if a[i][c] != 0 {
          if piv < 0 || a[i][c].abs() < a[piv][c].abs() {
            piv = i
          }
        }
      }
      if piv < 0 {
        done = true
        continue
      }
      if piv != r {
        let tmp = a[r]
        a[r] = a[piv]
        a[piv] = tmp
      }
      let mut changed = false
      for i = r + 1; i < rows; i = i + 1 {
        if a[i][c] == 0 {
          continue
        }
        let q = a[i][c] / a[r][c]
        for j = 0; j < cols; j = j + 1 {
          a[i][j] = a[i][j] - q * a[r][j]
        }
        changed = true
      }
      if !(changed) {
        done = true
      }
    }
    if a[r][c] != 0 {
      if a[r][c] < 0 {
        for j = 0; j < cols; j = j + 1 {
          a[r][j] = -a[r][j]
        }
      }
      r = r + 1
    }
  }
  a
}

///| **A30 Smith normal form**: the diagonal of the Smith form of an integer
///| matrix, i.e. the invariant factors `d_1 | d_2 | ...`.  By the structure
///| theorem for modules over a PID this *is* the integral homology of a
///| chain complex: Betti numbers = the zero diagonal entries, torsion = the
///| non-unit nonzero ones.
pub fn snf(m : Mat) -> Array[Int] {
  let a = mat_copy(m)
  let rows = mat_rows(a)
  let cols = mat_cols(a)
  let n = Int::min(rows, cols)
  ignore(rows)
  let diag : Array[Int] = []
  for t = 0; t < n; t = t + 1 {
    // find the smallest nonzero |entry| in the submatrix [t.., t..] and move it
    let mut steps = 0
    while steps < 10000 {
      let mut pi = -1
      let mut pj = -1
      for i = t; i < rows; i = i + 1 {
        for j = t; j < cols; j = j + 1 {
          if a[i][j] != 0 {
            if pi < 0 || a[i][j].abs() < a[pi][pj].abs() {
              pi = i
              pj = j
            }
          }
        }
      }
      if pi < 0 {
        steps = 10000
        continue
      }
      if pi != t {
        for j = 0; j < cols; j = j + 1 {
          let tmp = a[t][j]
          a[t][j] = a[pi][j]
          a[pi][j] = tmp
        }
      }
      if pj != t {
        for i = 0; i < rows; i = i + 1 {
          let tmp = a[i][t]
          a[i][t] = a[i][pj]
          a[i][pj] = tmp
        }
      }
      let pivot = a[t][t]
      let mut cleared = true
      for i = t + 1; i < rows; i = i + 1 {
        if a[i][t] != 0 {
          let q = a[i][t] / pivot
          for j = t; j < cols; j = j + 1 {
            a[i][j] = a[i][j] - q * a[t][j]
          }
          if a[i][t] != 0 {
            cleared = false
          }
        }
      }
      for j = t + 1; j < cols; j = j + 1 {
        if a[t][j] != 0 {
          let q = a[t][j] / pivot
          for i = t; i < rows; i = i + 1 {
            a[i][j] = a[i][j] - q * a[i][t]
          }
          if a[t][j] != 0 {
            cleared = false
          }
        }
      }
      if cleared {
        steps = 10000
      } else {
        steps = steps + 1
      }
    }
    diag.push(if t < rows && t < cols {
      a[t][t].abs()
    } else {
      0
    })
  }
  // make the diagonal divide successively
  let mut changed = true
  while changed {
    changed = false
    for i = 0; i + 1 < diag.length(); i = i + 1 {
      let x = diag[i]
      let y = diag[i + 1]
      if x == 0 && y != 0 {
        diag[i] = y
        diag[i + 1] = 0
        changed = true
      } else if x != 0 && y != 0 && y % x != 0 {
        let g = @arith.gcd(x, y)
        let l = @arith.lcm(x, y)
        diag[i] = g
        diag[i + 1] = l
        changed = true
      }
    }
  }
  diag
}

///| Rank of an integer matrix (over Q), via the Smith form.
pub fn rank_z(m : Mat) -> Int {
  let d = snf(m)
  let mut r = 0
  for i = 0; i < d.length(); i = i + 1 {
    if d[i] != 0 {
      r = r + 1
    }
  }
  r
}
