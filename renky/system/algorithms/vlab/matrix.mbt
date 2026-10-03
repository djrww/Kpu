///|
// =====================================================================
// matrix.mbt — 矩陣空間 M_{m×n}(R)：代數運算、消去法、分解、特徵值
// 定理自證於 matrix_wbtest.mbt（det 乘法性、LU/QR/Cholesky 重構、
// 譜定理、Cayley–Hamilton 等）
// =====================================================================

/// 矩陣（列主序）
pub struct Mat {
  rows : Int
  cols : Int
  data : Array[Double]
} derive(Debug, Eq)

/// M-1 零矩陣
pub fn mat_new(rows : Int, cols : Int) -> Mat {
  { rows, cols, data: Array::make(rows * cols, 0.0) }
}

/// M-2 由列向量陣列構造
pub fn mat_from_rows(rs : Array[Array[Double]]) -> Mat {
  let rows = rs.length()
  let cols = if rows == 0 { 0 } else { rs[0].length() }
  let data : Array[Double] = []
  for i = 0; i < rows; i = i + 1 {
    for j = 0; j < cols; j = j + 1 {
      data.push(rs[i][j])
    }
  }
  { rows, cols, data }
}

/// M-3 取元素
pub fn mat_get(m : Mat, i : Int, j : Int) -> Double {
  m.data[i * m.cols + j]
}

/// M-4 設元素
pub fn mat_set(m : Mat, i : Int, j : Int, v : Double) -> Unit {
  m.data[i * m.cols + j] = v
}

/// M-5 第 i 列
pub fn mat_row(m : Mat, i : Int) -> Array[Double] {
  let out : Array[Double] = []
  for j = 0; j < m.cols; j = j + 1 {
    out.push(mat_get(m, i, j))
  }
  out
}

/// M-6 第 j 行
pub fn mat_col(m : Mat, j : Int) -> Array[Double] {
  let out : Array[Double] = []
  for i = 0; i < m.rows; i = i + 1 {
    out.push(mat_get(m, i, j))
  }
  out
}

/// M-7 加法
pub fn mat_add(a : Mat, b : Mat) -> Mat {
  let out = mat_new(a.rows, a.cols)
  for i = 0; i < a.rows * a.cols; i = i + 1 {
    out.data[i] = a.data[i] + b.data[i]
  }
  out
}

/// M-8 減法
pub fn mat_sub(a : Mat, b : Mat) -> Mat {
  let out = mat_new(a.rows, a.cols)
  for i = 0; i < a.rows * a.cols; i = i + 1 {
    out.data[i] = a.data[i] - b.data[i]
  }
  out
}

/// M-9 純量乘法
pub fn mat_scale(c : Double, a : Mat) -> Mat {
  let out = mat_new(a.rows, a.cols)
  for i = 0; i < a.rows * a.cols; i = i + 1 {
    out.data[i] = c * a.data[i]
  }
  out
}

/// M-10 矩陣乘法 (AB)ᵢₖ = Σⱼ AᵢⱼBⱼₖ
pub fn mat_mul(a : Mat, b : Mat) -> Mat {
  let out = mat_new(a.rows, b.cols)
  for i = 0; i < a.rows; i = i + 1 {
    for k = 0; k < b.cols; k = k + 1 {
      let mut s = 0.0
      for j = 0; j < a.cols; j = j + 1 {
        s = s + mat_get(a, i, j) * mat_get(b, j, k)
      }
      mat_set(out, i, k, s)
    }
  }
  out
}

/// M-11 矩陣-向量乘
pub fn mat_vec(a : Mat, x : Array[Double]) -> Array[Double] {
  let out = Array::make(a.rows, 0.0)
  for i = 0; i < a.rows; i = i + 1 {
    let mut s = 0.0
    for j = 0; j < a.cols; j = j + 1 {
      s = s + mat_get(a, i, j) * x[j]
    }
    out[i] = s
  }
  out
}

/// M-12 轉置
pub fn mat_transpose(a : Mat) -> Mat {
  let out = mat_new(a.cols, a.rows)
  for i = 0; i < a.rows; i = i + 1 {
    for j = 0; j < a.cols; j = j + 1 {
      mat_set(out, j, i, mat_get(a, i, j))
    }
  }
  out
}

/// M-13 單位矩陣
pub fn mat_eye(n : Int) -> Mat {
  let out = mat_new(n, n)
  for i = 0; i < n; i = i + 1 {
    mat_set(out, i, i, 1.0)
  }
  out
}

/// M-14 跡 tr(A) = Σ Aᵢᵢ
pub fn mat_trace(a : Mat) -> Double {
  let mut s = 0.0
  for i = 0; i < a.rows; i = i + 1 {
    s = s + mat_get(a, i, i)
  }
  s
}

/// M-15 Frobenius 內積 ⟨A,B⟩_F = Σ AᵢⱼBᵢⱼ = tr(AᵀB)
pub fn mat_frob_inner(a : Mat, b : Mat) -> Double {
  let mut s = 0.0
  for i = 0; i < a.rows * a.cols; i = i + 1 {
    s = s + a.data[i] * b.data[i]
  }
  s
}

/// M-16 Frobenius 範數
pub fn mat_frob_norm(a : Mat) -> Double {
  mat_frob_inner(a, a).sqrt()
}

/// M-17 近似相等
pub fn mat_eq_eps(a : Mat, b : Mat, eps : Double) -> Bool {
  if a.rows != b.rows || a.cols != b.cols {
    false
  } else {
    let mut ok = true
    for i = 0; i < a.rows * a.cols; i = i + 1 {
      if (a.data[i] - b.data[i]).abs() >= eps {
        ok = false
      }
    }
    ok
  }
}

/// M-18 扁平方陣行列式（部分主元高斯消去）——供內積空間體積公式使用
pub fn flat_det(data : Array[Double], n : Int) -> Double {
  let m : Mat = { rows: n, cols: n, data: data.copy() }
  mat_det(m)
}

/// M-19 行列式（部分主元高斯消去，O(n³)）
pub fn mat_det(a : Mat) -> Double {
  let n = a.rows
  if n != a.cols {
    0.0
  } else if n == 0 {
    1.0
  } else {
    let m = a.data.copy()
    let mut sign = 1.0
    let mut det = 1.0
    for col = 0; col < n; col = col + 1 {
      // 部分主元
      let mut piv = col
      for r = col + 1; r < n; r = r + 1 {
        if m[r * n + col].abs() > m[piv * n + col].abs() {
          piv = r
        }
      }
      if m[piv * n + col].abs() < 1.0e-14 {
        det = 0.0
        break
      }
      if piv != col {
        for j = 0; j < n; j = j + 1 {
          let t = m[col * n + j]
          m[col * n + j] = m[piv * n + j]
          m[piv * n + j] = t
        }
        sign = -sign
      }
      let d = m[col * n + col]
      det = det * d
      for r = col + 1; r < n; r = r + 1 {
        let f = m[r * n + col] / d
        for j = col; j < n; j = j + 1 {
          m[r * n + j] = m[r * n + j] - f * m[col * n + j]
        }
      }
    }
    sign * det
  }
}

/// M-20 約化列梯形（RREF，部分主元，容差 eps）
pub fn mat_rref(a : Mat, eps : Double) -> Mat {
  let m = a.data.copy()
  let rows = a.rows
  let cols = a.cols
  let mut r = 0
  let mut lead = 0
  while r < rows && lead < cols {
    // 覓主元
    let mut piv = r
    for k = r + 1; k < rows; k = k + 1 {
      if m[k * cols + lead].abs() > m[piv * cols + lead].abs() {
        piv = k
      }
    }
    if m[piv * cols + lead].abs() < eps {
      // 此行無主元：往右移
      lead = lead + 1
    } else {
      if piv != r {
        for j = 0; j < cols; j = j + 1 {
          let t = m[r * cols + j]
          m[r * cols + j] = m[piv * cols + j]
          m[piv * cols + j] = t
        }
      }
      let d = m[r * cols + lead]
      for j = 0; j < cols; j = j + 1 {
        m[r * cols + j] = m[r * cols + j] / d
      }
      for k = 0; k < rows; k = k + 1 {
        if k != r {
          let f = m[k * cols + lead]
          for j = 0; j < cols; j = j + 1 {
            m[k * cols + j] = m[k * cols + j] - f * m[r * cols + j]
          }
        }
      }
      r = r + 1
      lead = lead + 1
    }
  }
  { rows, cols, data: m }
}
