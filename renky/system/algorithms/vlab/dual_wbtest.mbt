///|
// 對偶空間定理之自證

test "dual: 作用之線性" {
  let phi : DualVec = { coef: [1.0, -2.0, 3.0] }
  let u = [1.0, 1.0, 1.0]
  let v = [0.0, 2.0, 1.0]
  let a = 2.0
  let b = -1.0
  let lhs = dual_action(phi, euc_add(euc_scale(a, u), euc_scale(b, v)))
  let rhs = a * dual_action(phi, u) + b * dual_action(phi, v)
  assert_true((lhs - rhs).abs() < 1.0e-12)
}

test "dual: 對偶空間向量空間公理" {
  let ip = dual_ip(3)
  let phi : DualVec = { coef: [1.0, 0.0, 1.0] }
  let psi : DualVec = { coef: [0.0, 2.0, 1.0] }
  let s = (ip.sp.add)(phi, psi)
  assert_true(euc_eq_eps(s.coef, [1.0, 2.0, 2.0], 1.0e-12))
  let z = (ip.sp.add)(phi, (ip.sp.neg)(phi))
  assert_true(euc_eq_eps(z.coef, [0.0, 0.0, 0.0], 1.0e-12))
  // 內積 = 係數點積
  assert_true(((ip.inner)(phi, psi) - 1.0).abs() < 1.0e-12)
}

test "dual: 對偶基配對 = Kronecker delta" {
  let n = 4
  let phis = dual_basis(n)
  for i = 0; i < n; i = i + 1 {
    for j = 0; j < n; j = j + 1 {
      let want = if i == j { 1.0 } else { 0.0 }
      assert_true((dual_action(phis[i], evec(n, j)) - want).abs() < 1.0e-12)
    }
  }
  let m = dual_pairing_matrix(phis, n)
  assert_true(mat_eq_eps(m, mat_eye(n), 1.0e-12))
}

test "dual: Riesz 表示往返" {
  let v = [1.0, -2.0, 0.5]
  let phi = riesz_inverse(v)
  assert_true(euc_eq_eps(riesz_represent(phi), v, 1.0e-12))
  let w = [3.0, 3.0, 3.0]
  assert_true((dual_action(phi, w) - euc_inner(v, w)).abs() < 1.0e-12)
}

test "dual: 對偶映射與複合反序 (AB)* = B*A*" {
  let a = mat_from_rows([[1.0, 2.0], [3.0, 4.0]]) // R²→R²
  let b = mat_from_rows([[0.0, 1.0], [1.0, 0.0]])
  let lhs = dual_map(mat_mul(a, b))
  let rhs = mat_mul(dual_map(b), dual_map(a))
  assert_true(mat_eq_eps(lhs, rhs, 1.0e-12))
}

test "dual: 對偶映射之作用定義式" {
  let a = mat_from_rows([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]) // R³→R²
  let psi : DualVec = { coef: [1.0, -1.0] }
  let v = [1.0, 0.0, 2.0]
  let astar = dual_map(a) // 3×2
  let pulled : DualVec = { coef: mat_vec(astar, psi.coef) }
  assert_true((dual_action(pulled, v) - dual_action(psi, mat_vec(a, v))).abs() < 1.0e-12)
}

test "dual: 雙重對偶自然同構" {
  let v = [1.0, 2.0, 3.0]
  let phis = dual_basis(3)
  // J(v) 作用於 δᵢ 得 v 之分量
  assert_true(euc_eq_eps(double_dual_action(v, phis), v, 1.0e-12))
}

test "dual: 零化子維數定理 dim W° = n − dim W" {
  // W = span{e0, e1} ⊂ R³ 之零化子由 δ₂ 生成
  let w_gen = [evec(3, 0), evec(3, 1)]
  let phi0 : DualVec = { coef: [0.0, 0.0, 1.0] }
  let phi1 : DualVec = { coef: [0.0, 0.0, 5.0] }
  let bad : DualVec = { coef: [0.0, 1.0, 0.0] }
  assert_true(dual_in_annihilator(phi0, w_gen, 1.0e-12))
  assert_true(dual_in_annihilator(phi1, w_gen, 1.0e-12))
  assert_false(dual_in_annihilator(bad, w_gen, 1.0e-12))
}

test "dual: 維數相等 dim V = dim V*" {
  let n = 5
  assert_eq(dual_basis(n).length(), n)
  assert_eq(dual_ip(n).sp.zero.coef.length(), n)
}
