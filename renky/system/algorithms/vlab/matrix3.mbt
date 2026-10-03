///|
// =====================================================================
// matrix3.mbt — 矩陣演算法（續）：QR、子空間、最小二乘、張量積、指數
// =====================================================================

/// M-30 QR 分解（Gram–Schmidt 於行上）：A = QR，Q 行標準正交
pub fn mat_qr(a : Mat) -> (Mat, Mat) {
  let m = a.rows
  let n = a.cols
  let q_cols : Array[Array[Double]] = []
  let r = mat_new(n, n)
  for j = 0; j < n; j = j + 1 {
    let mut v = mat_col(a, j)
    for i = 0; i < q_cols.length(); i = i + 1 {
      let c = euc_inner(q_cols[i], v)
      mat_set(r, i, j, c)
      v = euc_add(v, euc_scale(-c, q_cols[i]))
    }
    let nrm = ip_norm(euc_ip(m), v)
    mat_set(r, j, j, nrm)
    if nrm > 1.0e-12 {
      q_cols.push(euc_scale(1.0 / nrm, v))
    } else {
      q_cols.push(euc_zero(m))
    }
  }
  let q = mat_new(m, n)
  for j = 0; j < n; j = j + 1 {
    for i = 0; i < m; i = i + 1 {
      mat_set(q, i, j, q_cols[j][i])
    }
  }
  (q, r)
}

/// M-31 列空間基（A 的主元行）
pub fn mat_col_space_basis(a : Mat, eps : Double) -> Array[Array[Double]] {
  let pivs = mat_pivot_cols(a, eps)
  let out : Array[Array[Double]] = []
  for j in pivs {
    out.push(mat_col(a, j))
  }
  out
}

/// M-32 行空間基（RREF 非零列）
pub fn mat_row_space_basis(a : Mat, eps : Double) -> Array[Array[Double]] {
  let r = mat_rref(a, eps)
  let out : Array[Array[Double]] = []
  for i = 0; i < r.rows; i = i + 1 {
    let row = mat_row(r, i)
    let mut nonzero = false
    for x in row {
      if x.abs() > eps {
        nonzero = true
      }
    }
    if nonzero {
      out.push(row)
    }
  }
  out
}

/// M-33 零空間基（自由變量法）
pub fn mat_null_space_basis(a : Mat, eps : Double) -> Array[Array[Double]] {
  let r = mat_rref(a, eps)
  let m = r.rows
  let n = r.cols
  let pivs = mat_pivot_cols(a, eps)
  let is_piv = Array::make(n, false)
  for p in pivs {
    is_piv[p] = true
  }
  let out : Array[Array[Double]] = []
  for f = 0; f < n; f = f + 1 {
    if !is_piv[f] {
      let v = Array::make(n, 0.0)
      v[f] = 1.0
      for i = 0; i < pivs.length(); i = i + 1 {
        v[pivs[i]] = -mat_get(r, i, f)
      }
      let _ = m
      out.push(v)
    }
  }
  out
}

/// M-34 最小二乘解（法方程 AᵀAx = Aᵀb，經 Cholesky 或 LU）
pub fn least_squares(a : Mat, b : Array[Double]) -> Array[Double] {
  let at = mat_transpose(a)
  let ata = mat_mul(at, a)
  let atb = mat_vec(at, b)
  match mat_cholesky(ata) {
    Some(l) => {
      // 前代 Ly = atb
      let n = atb.length()
      let y = Array::make(n, 0.0)
      for i = 0; i < n; i = i + 1 {
        let mut s = atb[i]
        for j = 0; j < i; j = j + 1 {
          s = s - mat_get(l, i, j) * y[j]
        }
        y[i] = s / mat_get(l, i, i)
      }
      // 回代 Lᵀx = y
      let x = Array::make(n, 0.0)
      for i = n - 1; i >= 0; i = i - 1 {
        let mut s = y[i]
        for j = i + 1; j < n; j = j + 1 {
          s = s - mat_get(l, j, i) * x[j]
        }
        x[i] = s / mat_get(l, i, i)
      }
      x
    }
    None => {
      let (l, u, perm) = mat_lu(ata)
      let _ = l
      lu_solve(l, u, perm, atb)
    }
  }
}

/// M-35 外積 u vᵀ
pub fn outer_product(u : Array[Double], v : Array[Double]) -> Mat {
  let out = mat_new(u.length(), v.length())
  for i = 0; i < u.length(); i = i + 1 {
    for j = 0; j < v.length(); j = j + 1 {
      mat_set(out, i, j, u[i] * v[j])
    }
  }
  out
}

/// M-36 Kronecker 積 A⊗B
pub fn mat_kronecker(a : Mat, b : Mat) -> Mat {
  let out = mat_new(a.rows * b.rows, a.cols * b.cols)
  for i = 0; i < a.rows; i = i + 1 {
    for j = 0; j < a.cols; j = j + 1 {
      for k = 0; k < b.rows; k = k + 1 {
        for l = 0; l < b.cols; l = l + 1 {
          mat_set(
            out,
            i * b.rows + k,
            j * b.cols + l,
            mat_get(a, i, j) * mat_get(b, k, l),
          )
        }
      }
    }
  }
  out
}

/// M-37 對稱化：sym(A) = (A + Aᵀ)/2（到對稱矩陣空間之正交投影）
pub fn mat_symmetrize(a : Mat) -> Mat {
  let at = mat_transpose(a)
  mat_scale(0.5, mat_add(a, at))
}

/// M-38 反對稱化：skew(A) = (A − Aᵀ)/2
pub fn mat_skew(a : Mat) -> Mat {
  let at = mat_transpose(a)
  mat_scale(0.5, mat_sub(a, at))
}

/// M-39 對稱判定
pub fn mat_is_symmetric(a : Mat, eps : Double) -> Bool {
  mat_eq_eps(a, mat_transpose(a), eps)
}

/// M-40 矩陣指數 e^A（scaling-and-squaring + Taylor 20 項）
pub fn mat_exp(a : Mat) -> Mat {
  // 選 s 使 ‖A/2^s‖∞ ≤ 1/2
  let mut nrm = 0.0
  for i = 0; i < a.rows; i = i + 1 {
    let mut rs = 0.0
    for j = 0; j < a.cols; j = j + 1 {
      rs = rs + mat_get(a, i, j).abs()
    }
    if rs > nrm {
      nrm = rs
    }
  }
  let mut s = 0
  let mut scale = 1.0
  while nrm / scale > 0.5 {
    scale = scale * 2.0
    s = s + 1
  }
  let b = mat_scale(1.0 / scale, a)
  // Taylor: E = I + B + B²/2! + ...
  let n = a.rows
  let mut e = mat_eye(n)
  let mut term = mat_eye(n)
  for k = 1; k <= 20; k = k + 1 {
    term = mat_scale(1.0 / k.to_double(), mat_mul(term, b))
    e = mat_add(e, term)
  }
  // 平方 s 次
  for it_ = 0; it_ < s; it_ = it_ + 1 {
    e = mat_mul(e, e)
  }
  e
}

/// M-41 譜分解重構：由 (特徵值, Q) 重建 A ≈ Q diag(λ) Qᵀ
pub fn spectral_reconstruct(vals : Array[Double], q : Mat) -> Mat {
  let n = vals.length()
  let d = mat_new(n, n)
  for i = 0; i < n; i = i + 1 {
    mat_set(d, i, i, vals[i])
  }
  mat_mul(q, mat_mul(d, mat_transpose(q)))
}

/// M-42 Rayleigh 商 xᵀAx / xᵀx
pub fn rayleigh_quotient(a : Mat, x : Array[Double]) -> Double {
  let ax = mat_vec(a, x)
  euc_inner(x, ax) / euc_inner(x, x)
}

/// M-43 對角矩陣構造
pub fn mat_diag(d : Array[Double]) -> Mat {
  let n = d.length()
  let out = mat_new(n, n)
  for i = 0; i < n; i = i + 1 {
    mat_set(out, i, i, d[i])
  }
  out
}

/// M-44 取對角線
pub fn mat_diagonal(a : Mat) -> Array[Double] {
  let n = if a.rows < a.cols { a.rows } else { a.cols }
  let out : Array[Double] = []
  for i = 0; i < n; i = i + 1 {
    out.push(mat_get(a, i, i))
  }
  out
}
