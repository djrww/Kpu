///|
// =====================================================================
// matrix2.mbt — 矩陣演算法（續）：秩、逆、分解、特徵值、子空間
// =====================================================================

/// M-21 秩 = RREF 中主元數
pub fn mat_rank(a : Mat, eps : Double) -> Int {
  let r = mat_rref(a, eps)
  let mut rank = 0
  for i = 0; i < r.rows; i = i + 1 {
    let mut nonzero = false
    for j = 0; j < r.cols; j = j + 1 {
      if mat_get(r, i, j).abs() > eps {
        nonzero = true
      }
    }
    if nonzero {
      rank = rank + 1
    }
  }
  rank
}

/// M-22 主元行位置（RREF 中每列首個非零元之行號）
pub fn mat_pivot_cols(a : Mat, eps : Double) -> Array[Int] {
  let r = mat_rref(a, eps)
  let pivs : Array[Int] = []
  for i = 0; i < r.rows; i = i + 1 {
    let mut j = 0
    while j < r.cols && mat_get(r, i, j).abs() < eps {
      j = j + 1
    }
    if j < r.cols {
      pivs.push(j)
    }
  }
  pivs
}

/// M-23 逆矩陣（Gauss–Jordan；奇異時回 None）
pub fn mat_inv(a : Mat) -> Option[Mat] {
  let n = a.rows
  // 增廣 [A | I]
  let aug = mat_new(n, 2 * n)
  for i = 0; i < n; i = i + 1 {
    for j = 0; j < n; j = j + 1 {
      mat_set(aug, i, j, mat_get(a, i, j))
    }
    mat_set(aug, i, n + i, 1.0)
  }
  let r = mat_rref(aug, 1.0e-12)
  // 左半須為單位陣
  let mut ok = true
  for i = 0; i < n; i = i + 1 {
    for j = 0; j < n; j = j + 1 {
      let want = if i == j { 1.0 } else { 0.0 }
      if (mat_get(r, i, j) - want).abs() > 1.0e-8 {
        ok = false
      }
    }
  }
  if !ok {
    Option::None
  } else {
    let out = mat_new(n, n)
    for i = 0; i < n; i = i + 1 {
      for j = 0; j < n; j = j + 1 {
        mat_set(out, i, j, mat_get(r, i, n + j))
      }
    }
    Some(out)
  }
}

/// M-24 解 Ax = b（唯一解；奇異時回 None）
pub fn mat_solve(a : Mat, b : Array[Double]) -> Option[Array[Double]] {
  match mat_inv(a) {
    Some(inv) => Some(mat_vec(inv, b))
    None => None
  }
}

/// M-25 LU 分解（部分主元）：PA = LU，回傳 (L, U, perm)
pub fn mat_lu(a : Mat) -> (Mat, Mat, Array[Int]) {
  let n = a.rows
  let u = a.data.copy()
  let l = mat_eye(n)
  let perm : Array[Int] = []
  for i = 0; i < n; i = i + 1 {
    perm.push(i)
  }
  for col = 0; col < n - 1; col = col + 1 {
    let mut piv = col
    for r = col + 1; r < n; r = r + 1 {
      if u[r * n + col].abs() > u[piv * n + col].abs() {
        piv = r
      }
    }
    if piv != col {
      for j = 0; j < n; j = j + 1 {
        let t = u[col * n + j]
        u[col * n + j] = u[piv * n + j]
        u[piv * n + j] = t
      }
      let tp = perm[col]
      perm[col] = perm[piv]
      perm[piv] = tp
      // 交換 L 已算出的部分
      for j = 0; j < col; j = j + 1 {
        let t = mat_get(l, col, j)
        mat_set(l, col, j, mat_get(l, piv, j))
        mat_set(l, piv, j, t)
      }
    }
    if u[col * n + col].abs() < 1.0e-15 {
      // 奇異：保持現狀
      ()
    } else {
      for r = col + 1; r < n; r = r + 1 {
        let f = u[r * n + col] / u[col * n + col]
        mat_set(l, r, col, f)
        for j = col; j < n; j = j + 1 {
          u[r * n + j] = u[r * n + j] - f * u[col * n + j]
        }
      }
    }
  }
  ({ rows: n, cols: n, data: l.data }, { rows: n, cols: n, data: u }, perm)
}

/// M-26 由 LU 解方程：先解 Ly = Pb，再解 Ux = y
pub fn lu_solve(
  l : Mat,
  u : Mat,
  perm : Array[Int],
  b : Array[Double],
) -> Array[Double] {
  let n = b.length()
  // 前代
  let y = Array::make(n, 0.0)
  for i = 0; i < n; i = i + 1 {
    let mut s = b[perm[i]]
    for j = 0; j < i; j = j + 1 {
      s = s - mat_get(l, i, j) * y[j]
    }
    y[i] = s
  }
  // 回代
  let x = Array::make(n, 0.0)
  for i = n - 1; i >= 0; i = i - 1 {
    let mut s = y[i]
    for j = i + 1; j < n; j = j + 1 {
      s = s - mat_get(u, i, j) * x[j]
    }
    x[i] = s / mat_get(u, i, i)
  }
  x
}

/// M-27 Cholesky 分解 A = LLᵀ（A 須對稱正定；失敗回 None）
pub fn mat_cholesky(a : Mat) -> Option[Mat] {
  let n = a.rows
  let l = mat_new(n, n)
  let mut ok = true
  for i = 0; i < n; i = i + 1 {
    for j = 0; j < i + 1; j = j + 1 {
      let mut s = 0.0
      for k = 0; k < j; k = k + 1 {
        s = s + mat_get(l, i, k) * mat_get(l, j, k)
      }
      if i == j {
        let d = mat_get(a, i, i) - s
        if d <= 0.0 {
          ok = false
        } else {
          mat_set(l, i, j, d.sqrt())
        }
      } else {
        let diag = mat_get(l, j, j)
        if diag == 0.0 {
          ok = false
        } else {
          mat_set(l, i, j, (mat_get(a, i, j) - s) / diag)
        }
      }
    }
  }
  if ok {
    Some(l)
  } else {
    Option::None
  }
}

/// M-28 冪迭代：主特徵值與特徵向量（Rayleigh 商）
pub fn power_iteration(a : Mat, iters : Int) -> (Double, Array[Double]) {
  let n = a.rows
  let x : Array[Double] = []
  for i = 0; i < n; i = i + 1 {
    x.push(@math.sin(i.to_double() * 12.9898 + 4.1414) + 1.1)
  }
  // 正規化
  let mut nrm = 0.0
  for v in x {
    nrm = nrm + v * v
  }
  nrm = nrm.sqrt()
  for i = 0; i < n; i = i + 1 {
    x[i] = x[i] / nrm
  }
  let mut lam = 0.0
  for it_ = 0; it_ < iters; it_ = it_ + 1 {
    let ax = mat_vec(a, x)
    // Rayleigh 商
    let mut num = 0.0
    for i = 0; i < n; i = i + 1 {
      num = num + x[i] * ax[i]
    }
    lam = num
    let mut n2 = 0.0
    for v in ax {
      n2 = n2 + v * v
    }
    n2 = n2.sqrt()
    if n2 > 0.0 {
      for i = 0; i < n; i = i + 1 {
        x[i] = ax[i] / n2
      }
    }
  }
  (lam, x)
}

/// M-29 Jacobi 旋轉法求對稱矩陣特徵值
/// 回傳 (特徵值, Q)，Q 之行為對應標準正交特徵向量
pub fn jacobi_eigen(a : Mat, max_sweeps : Int) -> (Array[Double], Mat) {
  let n = a.rows
  let s = a.data.copy()
  let q = mat_eye(n)
  for it_ = 0; it_ < max_sweeps; it_ = it_ + 1 {
    // 離对角範數
    let mut off = 0.0
    for i = 0; i < n; i = i + 1 {
      for j = i + 1; j < n; j = j + 1 {
        off = off + s[i * n + j] * s[i * n + j]
      }
    }
    if off.sqrt() < 1.0e-13 {
      break
    }
    for p = 0; p < n - 1; p = p + 1 {
      for qq = p + 1; qq < n; qq = qq + 1 {
        let apq = s[p * n + qq]
        if apq.abs() < 1.0e-16 {
          continue
        }
        let app = s[p * n + p]
        let aqq = s[qq * n + qq]
        let theta = (aqq - app) / (2.0 * apq)
        let t =
          if theta >= 0.0 {
            1.0 / (theta + (theta * theta + 1.0).sqrt())
          } else {
            -1.0 / (-theta + (theta * theta + 1.0).sqrt())
          }
        let c = 1.0 / (t * t + 1.0).sqrt()
        let ss = t * c
        // 旋轉 S = Gᵀ S G
        for i = 0; i < n; i = i + 1 {
          let sip = s[i * n + p]
          let siq = s[i * n + qq]
          s[i * n + p] = c * sip - ss * siq
          s[i * n + qq] = ss * sip + c * siq
        }
        for j = 0; j < n; j = j + 1 {
          let spj = s[p * n + j]
          let sqj = s[qq * n + j]
          s[p * n + j] = c * spj - ss * sqj
          s[qq * n + j] = ss * spj + c * sqj
        }
        // 累積 Q = Q G
        for i = 0; i < n; i = i + 1 {
          let qip = q.data[i * n + p]
          let qiq = q.data[i * n + qq]
          q.data[i * n + p] = c * qip - ss * qiq
          q.data[i * n + qq] = ss * qip + c * qiq
        }
      }
    }
  }
  let vals : Array[Double] = []
  for i = 0; i < n; i = i + 1 {
    vals.push(s[i * n + i])
  }
  (vals, { rows: n, cols: n, data: q.data })
}
