///|
// 對稱矩陣空間與對角矩陣空間定理之自證

test "symm: 維數公式" {
  assert_eq(sym_dim(1), 1)
  assert_eq(sym_dim(2), 3)
  assert_eq(sym_dim(3), 6)
  assert_eq(sym_dim(4), 10)
  assert_eq(sym_basis(3).length(), sym_dim(3))
}

test "symm: pack/unpack 與標準正交基" {
  let n = 3
  let basis = sym_basis(n)
  // 標準正交
  for i = 0; i < basis.length(); i = i + 1 {
    for j = 0; j < basis.length(); j = j + 1 {
      let want = if i == j { 1.0 } else { 0.0 }
      assert_true((sym_inner(basis[i], basis[j]) - want).abs() < 1.0e-12)
    }
  }
  // 任意對稱矩陣可由基重構
  let v = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
  let m = sym_pack(n, v)
  let coords = ip_coords(sym_ip(n), m, basis)
  let recon = ip_from_coords(sym_ip(n), coords, basis)
  assert_true(mat_eq_eps(recon, m, 1.0e-10))
}

test "symm: 對稱化 = 正交投影（冪等 + 像 ⊥ 核）" {
  let a = mat_from_rows([[1.0, 4.0, 0.0], [2.0, 5.0, 1.0], [3.0, -1.0, 2.0]])
  let s = mat_symmetrize(a)
  let k = mat_skew(a)
  assert_true(mat_is_symmetric(s, 1.0e-12))
  // a = s + k
  assert_true(mat_eq_eps(mat_add(s, k), a, 1.0e-12))
  // s ⊥ k（Frobenius）
  assert_true(mat_frob_inner(s, k).abs() < 1.0e-12)
  // 冪等
  assert_true(mat_eq_eps(mat_symmetrize(s), s, 1.0e-12))
  assert_true(mat_eq_eps(mat_skew(k), k, 1.0e-12))
}

test "symm: 譜定理（特徵值全實 + 正交特徵基）" {
  let a = mat_from_rows([[4.0, 1.0, 0.0], [1.0, 3.0, 1.0], [0.0, 1.0, 2.0]])
  assert_true(mat_is_symmetric(a, 1.0e-12))
  assert_true(sym_spectral_error(a) < 1.0e-8)
  let (vals, q) = jacobi_eigen(a, 60)
  // 特徵方程
  for j = 0; j < 3; j = j + 1 {
    let v = mat_col(q, j)
    assert_true(euc_eq_eps(mat_vec(a, v), euc_scale(vals[j], v), 1.0e-7))
  }
  // Rayleigh 商範圍
  let x = [1.0, 1.0, 1.0]
  let rq = rayleigh_quotient(a, x)
  let mut lo = vals[0]
  let mut hi = vals[0]
  for v in vals {
    if v < lo {
      lo = v
    }
    if v > hi {
      hi = v
    }
  }
  assert_true(rq >= lo - 1.0e-8 && rq <= hi + 1.0e-8)
}

test "symm: 二次型與極值（在單位球上介於 λmin, λmax）" {
  let a = mat_from_rows([[2.0, 1.0], [1.0, 2.0]])
  let x = ip_normalize(euc_ip(2), [3.0, 1.0])
  let qv = sym_quadratic(a, x)
  assert_true(qv >= 1.0 - 1.0e-9 && qv <= 3.0 + 1.0e-9)
}

test "diag: 與 R^n 之代數同構" {
  let a = [1.0, 2.0, 3.0]
  let b = [4.0, 5.0, 6.0]
  // 乘法 ↔ 逐點積
  let ma = diag_to_mat(a)
  let mb = diag_to_mat(b)
  assert_true(mat_eq_eps(diag_to_mat(diag_mul(a, b)), mat_mul(ma, mb), 1.0e-12))
  // 加法
  assert_true(mat_eq_eps(diag_to_mat(diag_add(a, b)), mat_add(ma, mb), 1.0e-12))
  // 內積 = Frobenius
  assert_true((diag_inner(a, b) - mat_frob_inner(ma, mb)).abs() < 1.0e-12)
}

test "diag: 行列式、跡、逆" {
  let a = [2.0, -1.0, 3.0]
  assert_true((diag_det(a) - (-6.0)).abs() < 1.0e-12)
  assert_true((diag_trace(a) - 4.0).abs() < 1.0e-12)
  let inv = diag_inv(a).unwrap()
  assert_true(euc_eq_eps(diag_mul(a, inv), [1.0, 1.0, 1.0], 1.0e-12))
  // 含零者無逆
  match diag_inv([1.0, 0.0]) {
    Some(_) => assert_true(false)
    None => ()
  }
}

test "diag: 作用於向量與譜半徑" {
  let a = [2.0, -3.0, 1.0]
  let x = [1.0, 1.0, 1.0]
  assert_true(euc_eq_eps(diag_apply(a, x), a, 1.0e-12))
  assert_true((diag_spectral_radius(a) - 3.0).abs() < 1.0e-12)
  // 與 Jacobi 特徵值一致
  let (vals, _) = jacobi_eigen(diag_to_mat(a), 40)
  let mut mx = 0.0
  for v in vals {
    if v.abs() > mx {
      mx = v.abs()
    }
  }
  assert_true((mx - diag_spectral_radius(a)).abs() < 1.0e-8)
}

test "diag: 內積空間結構" {
  let ip = diag_ip(3)
  let a = [1.0, 2.0, 3.0]
  let b = [1.0, 0.0, 1.0]
  assert_true(((ip.inner)(a, b) - 4.0).abs() < 1.0e-12)
  let onb = ip_orthonormal(ip, [[1.0, 1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 2.0]], 1.0e-10)
  assert_eq(onb.length(), 3)
}
