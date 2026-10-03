///|
// 矩陣空間定理之自證

test "matrix: 乘法與轉置性質" {
  let a = mat_from_rows([[1.0, 2.0], [3.0, 4.0]])
  let b = mat_from_rows([[0.0, 1.0], [1.0, 2.0]])
  // (AB)ᵀ = BᵀAᵀ
  assert_true(
    mat_eq_eps(
      mat_transpose(mat_mul(a, b)),
      mat_mul(mat_transpose(b), mat_transpose(a)),
      1.0e-12,
    ),
  )
  // 具體積
  assert_true(mat_eq_eps(mat_mul(a, b), mat_from_rows([[2.0, 5.0], [4.0, 11.0]]), 1.0e-12))
  // 單位元
  assert_true(mat_eq_eps(mat_mul(a, mat_eye(2)), a, 1.0e-12))
}

test "matrix: 行列式基本性質" {
  let a = mat_from_rows([[3.0, 1.0], [2.0, 4.0]])
  assert_true((mat_det(a) - 10.0).abs() < 1.0e-12)
  // 轉置不變
  assert_true((mat_det(mat_transpose(a)) - mat_det(a)).abs() < 1.0e-10)
  // 純量倍：det(cA) = c^n det(A)
  assert_true((mat_det(mat_scale(2.0, a)) - 4.0 * mat_det(a)).abs() < 1.0e-10)
  // 奇異矩陣
  let sing = mat_from_rows([[1.0, 2.0], [2.0, 4.0]])
  assert_true(mat_det(sing).abs() < 1.0e-10)
  // 3×3 具體值
  let b = mat_from_rows([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 10.0]])
  assert_true((mat_det(b) - (-3.0)).abs() < 1.0e-10)
}

test "matrix: RREF 具體形狀" {
  let a = mat_from_rows([[1.0, 2.0, 3.0], [2.0, 4.0, 7.0]])
  let r = mat_rref(a, 1.0e-12)
  assert_true(
    mat_eq_eps(r, mat_from_rows([[1.0, 2.0, 0.0], [0.0, 0.0, 1.0]]), 1.0e-9),
  )
  assert_eq(mat_rank(a, 1.0e-10), 2)
  let dep = mat_from_rows([[1.0, 2.0], [2.0, 4.0], [3.0, 6.0]])
  assert_eq(mat_rank(dep, 1.0e-10), 1)
}

test "matrix: 逆矩陣與解方程" {
  let a = mat_from_rows([[4.0, 7.0], [2.0, 6.0]])
  let inv = mat_inv(a).unwrap()
  assert_true(mat_eq_eps(mat_mul(inv, a), mat_eye(2), 1.0e-9))
  // 解 Ax = b
  let b = [1.0, 2.0]
  let x = mat_solve(a, b).unwrap()
  assert_true(euc_eq_eps(mat_vec(a, x), b, 1.0e-9))
  // 具體解：x = (-0.8, 0.6)
  assert_true(euc_eq_eps(x, [-0.8, 0.6], 1.0e-9))
  // 奇異矩陣無逆
  let sing = mat_from_rows([[1.0, 1.0], [1.0, 1.0]])
  match mat_inv(sing) {
    Some(_) => assert_true(false)
    None => ()
  }
}

test "matrix: LU 解方程與行列式" {
  let a = mat_from_rows([[2.0, 1.0, 1.0], [4.0, 3.0, 3.0], [8.0, 7.0, 9.0]])
  let (l, u, perm) = mat_lu(a)
  let b = [1.0, 1.0, 1.0]
  let x = lu_solve(l, u, perm, b)
  assert_true(euc_eq_eps(mat_vec(a, x), b, 1.0e-9))
  // det(U) 對角積 × 符號 = det(A)；驗證與直接計算一致
  let mut du = 1.0
  for i = 0; i < 3; i = i + 1 {
    du = du * mat_get(u, i, i)
  }
  let mut swaps = 0
  let visited = Array::make(3, false)
  for i = 0; i < 3; i = i + 1 {
    if !visited[i] {
      let mut j = i
      let mut clen = 0
      while !visited[j] {
        visited[j] = true
        j = perm[j]
        clen = clen + 1
      }
      swaps = swaps + clen - 1
    }
  }
  let sign = if swaps % 2 == 0 { 1.0 } else { -1.0 }
  assert_true((sign * du - mat_det(a)).abs() < 1.0e-9)
}

test "matrix: Cholesky 與求解" {
  // 對稱正定
  let a = mat_from_rows([[4.0, 2.0, 0.0], [2.0, 5.0, 1.0], [0.0, 1.0, 6.0]])
  let l = mat_cholesky(a).unwrap()
  assert_true(mat_eq_eps(mat_mul(l, mat_transpose(l)), a, 1.0e-9))
  // L 下三角
  assert_true(mat_get(l, 0, 1).abs() < 1.0e-12)
  assert_true(mat_get(l, 0, 2).abs() < 1.0e-12)
  assert_true(mat_get(l, 1, 2).abs() < 1.0e-12)
  // 非正定者失敗
  let bad = mat_from_rows([[1.0, 2.0], [2.0, 1.0]])
  match mat_cholesky(bad) {
    Some(_) => assert_true(false)
    None => ()
  }
}

test "matrix: QR 與最小二乘" {
  let a = mat_from_rows([[1.0, 1.0], [1.0, 2.0], [1.0, 3.0]])
  let (q, r) = mat_qr(a)
  assert_true(mat_eq_eps(mat_mul(q, r), a, 1.0e-9))
  // R 上三角
  assert_true(mat_get(r, 1, 0).abs() < 1.0e-9)
  // 最小二乘配線性數據須精確
  let b = [2.0, 3.0, 4.0] // 完美線性：斜率 1 截距 1
  let x = least_squares(a, b)
  assert_true(euc_eq_eps(x, [1.0, 1.0], 1.0e-8))
  // 非完美數據：殘差正交於列空間
  let b2 = [1.0, 2.0, 4.0]
  let x2 = least_squares(a, b2)
  let resid = euc_add(b2, euc_neg(mat_vec(a, x2)))
  assert_true(euc_inner(resid, mat_col(a, 0)).abs() < 1.0e-8)
  assert_true(euc_inner(resid, mat_col(a, 1)).abs() < 1.0e-8)
}

test "matrix: 對稱矩陣譜定理" {
  let a = mat_from_rows([[2.0, -1.0, 0.0], [-1.0, 2.0, -1.0], [0.0, -1.0, 2.0]])
  let (vals, q) = jacobi_eigen(a, 60)
  // 已知特徵值 2−√2, 2, 2+√2
  let sorted = vals.copy()
  sorted.sort()
  assert_true((sorted[0] - (2.0 - (2.0).sqrt())).abs() < 1.0e-8)
  assert_true((sorted[1] - 2.0).abs() < 1.0e-8)
  assert_true((sorted[2] - (2.0 + (2.0).sqrt())).abs() < 1.0e-8)
  // Q 正交
  assert_true(mat_eq_eps(mat_mul(mat_transpose(q), q), mat_eye(3), 1.0e-8))
  // AQ = QΛ（每行特徵方程）
  for j = 0; j < 3; j = j + 1 {
    let v = mat_col(q, j)
    let av = mat_vec(a, v)
    assert_true(euc_eq_eps(av, euc_scale(vals[j], v), 1.0e-7))
  }
  // 重構
  assert_true(mat_eq_eps(spectral_reconstruct(vals, q), a, 1.0e-8))
  // 跡 = 特徵值和、det = 特徵值積
  let mut vsum2 = 0.0
  let mut vprod = 1.0
  for v in vals {
    vsum2 = vsum2 + v
    vprod = vprod * v
  }
  assert_true((vsum2 - mat_trace(a)).abs() < 1.0e-8)
  assert_true((vprod - mat_det(a)).abs() < 1.0e-8)
}

test "matrix: 冪迭代收斂到主特徵值" {
  let a = mat_from_rows([[5.0, 1.0], [1.0, 2.0]])
  let (lam, v) = power_iteration(a, 300)
  // 主特徵值 = (7 + √(9+4))/2
  let expect = (7.0 + (13.0).sqrt()) / 2.0
  assert_true((lam - expect).abs() < 1.0e-6)
  assert_true((ip_norm(euc_ip(2), v) - 1.0).abs() < 1.0e-9)
}

test "matrix: 子空間定理（秩-零度）" {
  let a = mat_from_rows([[1.0, 2.0, 3.0, 4.0], [2.0, 4.0, 6.0, 9.0]])
  let rank = mat_rank(a, 1.0e-10)
  let nulls = mat_null_space_basis(a, 1.0e-10)
  assert_eq(rank + nulls.length(), 4) // 秩-零度定理
  // 零空間基確實被 A 消去
  for v in nulls {
    assert_true(euc_eq_eps(mat_vec(a, v), [0.0, 0.0], 1.0e-8))
  }
  // 列空間基生成 A 之各行
  let cols = mat_col_space_basis(a, 1.0e-10)
  assert_eq(cols.length(), rank)
  let rows = mat_row_space_basis(a, 1.0e-10)
  assert_eq(rows.length(), rank)
}

test "matrix: Cayley–Hamilton 定理" {
  // p(A) = 0，p 為特徵多項式；2×2：A² − tr(A)A + det(A)I = 0
  let a = mat_from_rows([[3.0, 1.0], [-2.0, 0.0]])
  let a2 = mat_mul(a, a)
  let ch = mat_add(
    mat_add(a2, mat_scale(-mat_trace(a), a)),
    mat_scale(mat_det(a), mat_eye(2)),
  )
  assert_true(mat_eq_eps(ch, mat_new(2, 2), 1.0e-9))
}

test "matrix: 矩陣指數性質" {
  // 對角陣之指數 = 逐元素指數
  let d = mat_diag([1.0, -1.0])
  let ed = mat_exp(d)
  assert_true((mat_get(ed, 0, 0) - @math.exp(1.0)).abs() < 1.0e-9)
  assert_true((mat_get(ed, 1, 1) - @math.exp(-1.0)).abs() < 1.0e-9)
  assert_true(mat_get(ed, 0, 1).abs() < 1.0e-9)
  // det(e^A) = e^{tr A}（於 2×2 反對稱陣：旋轉，det=1）
  let skew = mat_from_rows([[0.0, -0.7], [0.7, 0.0]])
  let esk = mat_exp(skew)
  assert_true((mat_det(esk) - 1.0).abs() < 1.0e-9)
  // 且為旋轉陣：正交
  assert_true(mat_eq_eps(mat_mul(esk, mat_transpose(esk)), mat_eye(2), 1.0e-9))
}

test "matrix: Kronecker 積定理" {
  let a = mat_from_rows([[1.0, 2.0], [3.0, 4.0]])
  let b = mat_from_rows([[0.0, 1.0], [1.0, 0.0]])
  let k = mat_kronecker(a, b)
  // det(A⊗B) = det(A)^n det(B)^m（m=n=2）
  assert_true((mat_det(k) - mat_det(a) * mat_det(a) * mat_det(b) * mat_det(b)).abs() < 1.0e-8)
  // (A⊗B)(x⊗y) = Ax ⊗ By
  let x = [1.0, -1.0]
  let y = [2.0, 3.0]
  let xy : Array[Double] = []
  for xi in x {
    for yi in y {
      xy.push(xi * yi)
    }
  }
  let ax = mat_vec(a, x)
  let by = mat_vec(b, y)
  let expect : Array[Double] = []
  for ai in ax {
    for bi in by {
      expect.push(ai * bi)
    }
  }
  assert_true(euc_eq_eps(mat_vec(k, xy), expect, 1.0e-9))
}

test "matrix: 對稱/反對稱分解正交性" {
  let a = mat_from_rows([[1.0, 5.0], [3.0, 2.0]])
  let s = mat_symmetrize(a)
  let k = mat_skew(a)
  // A = sym + skew
  assert_true(mat_eq_eps(mat_add(s, k), a, 1.0e-12))
  // sym ⊥ skew（Frobenius）
  assert_true(mat_frob_inner(s, k).abs() < 1.0e-12)
  // 冪等
  assert_true(mat_eq_eps(mat_symmetrize(s), s, 1.0e-12))
  assert_true(mat_eq_eps(mat_skew(k), k, 1.0e-12))
}

test "matrix: Frobenius 內積 = tr(AᵀB)" {
  let a = mat_from_rows([[1.0, 2.0], [3.0, 4.0]])
  let b = mat_from_rows([[5.0, 6.0], [7.0, 8.0]])
  let via_trace = mat_trace(mat_mul(mat_transpose(a), b))
  assert_true((mat_frob_inner(a, b) - via_trace).abs() < 1.0e-12)
  assert_true((mat_frob_inner(a, b) - 70.0).abs() < 1.0e-12)
}

test "matrix: 外積與跡" {
  let u = [1.0, 2.0, 3.0]
  let v = [4.0, 5.0, 6.0]
  let op = outer_product(u, v)
  assert_true((mat_trace(op) - euc_inner(u, v)).abs() < 1.0e-12)
  // op · x = u (v·x)
  let x = [1.0, 0.0, -1.0]
  assert_true(euc_eq_eps(mat_vec(op, x), euc_scale(euc_inner(v, x), u), 1.0e-12))
}
