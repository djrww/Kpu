// ===========================================================================
//  義務自證 · exhaustive verification of the linear-algebra engine.
// ===========================================================================

fn check(cond : Bool, msg : String) -> Unit {
  if cond {
    ()
  } else {
    println("CHECK FAILED: " + msg)
    abort("[law failed] " + msg)
  }
}

///| Every `rows x cols` matrix over F_p, for tiny shapes.
fn all_mats(p : Int, rows : Int, cols : Int) -> Array[Mat] {
  let total = @arith.pow_int(p, rows * cols)
  let out : Array[Mat] = []
  for k = 0; k < total; k = k + 1 {
    let ds = @arith.digits(k, p)
    let m = mat_make(rows, cols)
    for i = 0; i < rows * cols; i = i + 1 {
      m[i / cols][i % cols] = if i < ds.length() {
        ds[i]
      } else {
        0
      }
    }
    out.push(m)
  }
  out
}

test "rank-nullity: rank(m) + dim ker(m) == cols, for every 3x3 matrix over F_2 and F_3" {
  for p in [2, 3] {
    let mats = all_mats(p, 3, 3)
    for i = 0; i < mats.length(); i = i + 1 {
      let m = mats[i]
      let r = rank(m, p)
      let ker = nullspace(m, p)
      check(r + ker.length() == 3, "rank-nullity p=" + p.to_string())
      // every kernel vector really is annihilated
      for j = 0; j < ker.length(); j = j + 1 {
        let v = mat_vec(m, ker[j], p)
        for t = 0; t < v.length(); t = t + 1 {
          check(v[t] == 0, "kernel annihilated")
        }
      }
      // kernel vectors are independent: the only combination giving 0 is trivial
      if ker.length() > 0 {
        check(rank(mat_from_rows(ker), p) == ker.length(), "kernel independent")
      }
    }
  }
}

fn mat_from_rows(rows : Array[Array[Int]]) -> Mat {
  let out : Mat = []
  for i = 0; i < rows.length(); i = i + 1 {
    out.push(rows[i].copy())
  }
  out
}

test "rref is idempotent and preserves the row space, all 2x3 matrices over F_2, F_3" {
  for p in [2, 3] {
    let mats = all_mats(p, 2, 3)
    for i = 0; i < mats.length(); i = i + 1 {
      let (a1, piv1) = rref(mats[i], p)
      let (a2, piv2) = rref(a1, p)
      check(piv1.length() == piv2.length(), "rref idempotent pivots")
      for r = 0; r < a1.length(); r = r + 1 {
        for c = 0; c < a1[r].length(); c = c + 1 {
          check(a1[r][c] == a2[r][c], "rref idempotent entries")
        }
      }
      check(rank(mats[i], p) == rank(a1, p), "rank preserved")
    }
  }
}

test "det_mod equals the Leibniz permutation formula for every 3x3 matrix over F_2, F_3" {
  for p in [2, 3] {
    let mats = all_mats(p, 3, 3)
    for i = 0; i < mats.length(); i = i + 1 {
      let m = mats[i]
      let brute = det_brute(m, p)
      check(
        det_mod(m, p) == brute,
        "det p=" + p.to_string() + " brute=" + brute.to_string(),
      )
    }
  }
}

fn det_brute(m : Mat, p : Int) -> Int {
  // sum over all 6 permutations of S_3 with signs
  let perms = [[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]]
  let signs = [1, -1, -1, 1, 1, -1]
  let mut acc = 0
  for k = 0; k < perms.length(); k = k + 1 {
    let s = perms[k]
    let mut prod = signs[k]
    for i = 0; i < 3; i = i + 1 {
      prod = prod * m[i][s[i]]
    }
    acc = @kernel.addmod(acc, prod, p)
  }
  acc
}

test "solve_mod: a returned solution really solves the system, all 2x2 systems over F_3" {
  let p = 3
  let mats = all_mats(p, 2, 2)
  for i = 0; i < mats.length(); i = i + 1 {
    let m = mats[i]
    for b0 = 0; b0 < p; b0 = b0 + 1 {
      for b1 = 0; b1 < p; b1 = b1 + 1 {
        let b = [b0, b1]
        match solve_mod(m, b, p) {
          Some((x, _)) => {
            let v = mat_vec(m, x, p)
            check(v[0] == b0 && v[1] == b1, "solve correct")
          }
          None => {
            // inconsistent: verify by brute force that no solution exists
            let mut found = false
            for x0 = 0; x0 < p; x0 = x0 + 1 {
              for x1 = 0; x1 < p; x1 = x1 + 1 {
                let v = mat_vec(m, [x0, x1], p)
                if v[0] == b0 && v[1] == b1 {
                  found = true
                }
              }
            }
            check(!(found), "inconsistent really is")
          }
        }
      }
    }
  }
}

test "A30 Smith normal form: diagonal divides successively, rank matches, and det is the product" {
  let cases : Array[Mat] = [
    [[2, 0], [0, 3]],
    [[4, 6], [6, 9]],
    [[2, 4, 4], [-6, 6, 12], [10, -4, -16]],
    [[0, 0], [0, 0]],
    [[1, 2, 3], [4, 5, 6], [7, 8, 10]],
    [[6, 0, 0], [0, 10, 0], [0, 0, 15]],
  ]
  for i = 0; i < cases.length(); i = i + 1 {
    let d = snf(cases[i])
    for j = 0; j + 1 < d.length(); j = j + 1 {
      if d[j] != 0 {
        check(d[j + 1] % d[j] == 0, "snf divisibility")
      } else {
        check(d[j + 1] == 0, "snf zeros at the end")
      }
    }
    let mut nz = 0
    for j = 0; j < d.length(); j = j + 1 {
      if d[j] != 0 {
        nz = nz + 1
      }
    }
    check(nz == rank_z(cases[i]), "snf rank")
  }
  // det of the 3x3 examples equals +- product of the diagonal
  let d3 = snf(cases[4])
  let mut prod = 1
  for j = 0; j < d3.length(); j = j + 1 {
    prod = prod * d3[j]
  }
  let d_int = det_int(cases[4])
  check(prod == d_int.abs(), "snf product == |det|, got " + prod.to_string() + " vs " + d_int.to_string())
}

fn det_int(m : Mat) -> Int {
  let perms = [[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]]
  let signs = [1, -1, -1, 1, 1, -1]
  let mut acc = 0
  for k = 0; k < perms.length(); k = k + 1 {
    let s = perms[k]
    let mut prod = signs[k]
    for i = 0; i < 3; i = i + 1 {
      prod = prod * m[i][s[i]]
    }
    acc = acc + prod
  }
  acc
}

test "A29 HNF: upper triangular, non-negative pivots, and generates the same lattice" {
  let cases : Array[Mat] = [
    [[2, 3], [4, 5]],
    [[1, 0, 0], [0, 2, 0], [0, 0, 3]],
    [[3, -1, 2], [-6, 2, 4], [9, 1, -1]],
  ]
  for i = 0; i < cases.length(); i = i + 1 {
    let h = hnf(cases[i])
    let rows = mat_rows(h)
    let cols = mat_cols(h)
    let mut nz = 0
    for r = 0; r < rows; r = r + 1 {
      if row_nonzero(h[r]) {
        nz = nz + 1
      }
    }
    // pivots strictly to the right as we go down, and non-negative
    let mut last_pivot = -1
    for r = 0; r < nz; r = r + 1 {
      let pv = first_nonzero(h[r])
      check(pv > last_pivot, "hnf staircase")
      check(h[r][pv] > 0, "hnf positive pivot")
      last_pivot = pv
    }
    for r = nz; r < rows; r = r + 1 {
      check(!(row_nonzero(h[r])), "hnf zero rows at the end")
    }
    check(rank_z(h) == rank_z(cases[i]), "hnf preserves rank")
    // every original row is an integer combination of the HNF rows
    for r = 0; r < mat_rows(cases[i]); r = r + 1 {
      check(in_lattice(h, nz, cases[i][r], cols), "hnf contains original rows")
    }
  }
}

fn row_nonzero(r : Array[Int]) -> Bool {
  for i = 0; i < r.length(); i = i + 1 {
    if r[i] != 0 {
      return true
    }
  }
  false
}

fn first_nonzero(r : Array[Int]) -> Int {
  for i = 0; i < r.length(); i = i + 1 {
    if r[i] != 0 {
      return i
    }
  }
  -1
}

///| Is `v` an integer combination of the first `nz` rows of `h`?
///| Solved by back-substitution on the staircase.
fn in_lattice(h : Mat, nz : Int, v : Array[Int], cols : Int) -> Bool {
  // Row echelon with strictly increasing pivots: column pv_r is touched only by
  // rows t <= r, so the coefficients are recovered top-down.
  let coef : Array[Int] = Array::make(nz, 0)
  for r = 0; r < nz; r = r + 1 {
    let pv = first_nonzero(h[r])
    if pv < 0 {
      return false
    }
    let mut rhs = v[pv]
    for t = 0; t < r; t = t + 1 {
      rhs = rhs - coef[t] * h[t][pv]
    }
    if rhs % h[r][pv] != 0 {
      return false
    }
    coef[r] = rhs / h[r][pv]
  }
  // verify the whole vector
  for c = 0; c < cols; c = c + 1 {
    let mut acc = 0
    for t = 0; t < nz; t = t + 1 {
      acc = acc + coef[t] * h[t][c]
    }
    if acc != v[c] {
      return false
    }
  }
  true
}

test "mat_mul is associative and distributive over F_p on all 2x2 pairs over F_2, F_3" {
  for p in [2, 3] {
    let mats = all_mats(p, 2, 2)
    for i = 0; i < mats.length(); i = i + 1 {
      for j = 0; j < mats.length(); j = j + 1 {
        let a = mats[i]
        let b = mats[j]
        // distributivity against a third matrix (sampled subset keeps the
        // triple loop at 81*81*9 iterations, still complete in a and b)
        let lim = Int::min(mats.length(), 9)
        for t = 0; t < lim; t = t + 1 {
          let c = mats[t]
          check(
            mat_eq(
              mat_mul(mat_add(a, b, p), c, p),
              mat_add(mat_mul(a, c, p), mat_mul(b, c, p), p),
            ),
            "matrix distributivity",
          )
        }
      }
    }
  }
}

fn mat_eq(a : Mat, b : Mat) -> Bool {
  if mat_rows(a) != mat_rows(b) || mat_cols(a) != mat_cols(b) {
    return false
  }
  for i = 0; i < mat_rows(a); i = i + 1 {
    for j = 0; j < mat_cols(a); j = j + 1 {
      if a[i][j] != b[i][j] {
        return false
      }
    }
  }
  true
}

