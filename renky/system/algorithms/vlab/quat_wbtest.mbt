///|
// 四元數空間定理之自證

test "quat: 分配律與結合律" {
  let a = quat_new(1.0, 2.0, -1.0, 0.0)
  let b = quat_new(0.0, 1.0, 1.0, 2.0)
  let c = quat_new(-1.0, 0.0, 2.0, 1.0)
  let lhs = quat_mul(a, quat_add(b, c))
  let rhs = quat_add(quat_mul(a, b), quat_mul(a, c))
  assert_true(ceq_q2(lhs, rhs))
  let lhs2 = quat_mul(quat_mul(a, b), c)
  let rhs2 = quat_mul(a, quat_mul(b, c))
  assert_true(ceq_q2(lhs2, rhs2))
  // 非交換！
  let ab = quat_mul(a, b)
  let ba = quat_mul(b, a)
  assert_false(ceq_q2(ab, ba))
}

test "quat: 共軛反序律" {
  let a = quat_new(1.0, 1.0, 1.0, 1.0)
  let b = quat_new(2.0, -1.0, 0.0, 3.0)
  assert_true(ceq_q2(quat_conj(quat_mul(a, b)), quat_mul(quat_conj(b), quat_conj(a))))
}

test "quat: q̄q = ‖q‖² 為實" {
  let q = quat_new(1.0, -2.0, 3.0, 0.5)
  let p = quat_mul(quat_conj(q), q)
  assert_true((p.w - quat_norm2(q)).abs() < 1.0e-12)
  assert_true(p.x.abs() < 1.0e-12 && p.y.abs() < 1.0e-12 && p.z.abs() < 1.0e-12)
}

test "quat: 逆元" {
  let q = quat_new(1.0, 2.0, -2.0, 1.0)
  assert_true(ceq_q2(quat_mul(q, quat_inv(q)), quat_real(1.0)))
  assert_true(ceq_q2(quat_mul(quat_inv(q), q), quat_real(1.0)))
}

test "quat: 純虛乘法 ↔ 叉積與點積" {
  // pq = −p·q + p×q（純虛四元數）
  let u = [1.0, 2.0, 3.0]
  let v = [-1.0, 0.0, 2.0]
  let p = quat_imag(u[0], u[1], u[2])
  let q = quat_imag(v[0], v[1], v[2])
  let pq = quat_mul(p, q)
  assert_true((pq.w + euc_inner(u, v)).abs() < 1.0e-12)
  let cr = euc_cross3(u, v)
  assert_true((pq.x - cr[0]).abs() < 1.0e-12)
  assert_true((pq.y - cr[1]).abs() < 1.0e-12)
  assert_true((pq.z - cr[2]).abs() < 1.0e-12)
}

test "quat: 旋轉複合 = 四元數乘法" {
  let q1 = quat_from_axis_angle([0.0, 0.0, 1.0], @math.PI / 2.0)
  let q2 = quat_from_axis_angle([1.0, 0.0, 0.0], @math.PI / 2.0)
  let v = [0.0, 1.0, 0.0]
  let r1 = quat_rotate_vec(q1, quat_rotate_vec(q2, v))
  let r2 = quat_rotate_vec(quat_compose(q2, q1), v)
  assert_true(euc_eq_eps(r1, r2, 1.0e-9))
}

test "quat: 旋轉是正交變換（保內積）" {
  let q = quat_from_axis_angle([1.0, 1.0, 1.0], 1.234)
  let u = [1.0, -2.0, 0.5]
  let v = [3.0, 0.0, 1.0]
  let ru = quat_rotate_vec(q, u)
  let rv = quat_rotate_vec(q, v)
  assert_true((euc_inner(ru, rv) - euc_inner(u, v)).abs() < 1.0e-9)
  // 軸不動
  let axis = [1.0, 1.0, 1.0]
  let raxis = quat_rotate_vec(q, axis)
  assert_true(euc_eq_eps(raxis, axis, 1.0e-9))
}

test "quat: H 之內積為實內積" {
  let a = quat_new(1.0, 2.0, 0.0, 1.0)
  let b = quat_new(0.0, 1.0, -1.0, 2.0)
  // 對稱
  assert_true((quat_inner(a, b) - quat_inner(b, a)).abs() < 1.0e-12)
  // 與 4 維點積一致
  assert_true(
    (quat_inner(a, b) - (a.w * b.w + a.x * b.x + a.y * b.y + a.z * b.z)).abs() < 1.0e-12,
  )
}

test "quat: H^n 作為實內積空間之公理" {
  let ip = hvec_ip(2)
  let u : HVec = {
    comps: [quat_new(1.0, 1.0, 0.0, 0.0), quat_new(0.0, 0.0, 1.0, 0.0)],
  }
  let v : HVec = {
    comps: [quat_new(0.0, 1.0, 1.0, 0.0), quat_new(2.0, 0.0, 0.0, 0.0)],
  }
  // 對稱、雙線性、正定
  assert_true(((ip.inner)(u, v) - (ip.inner)(v, u)).abs() < 1.0e-12)
  assert_true((ip.inner)(u, u) > 0.0)
  let s = (ip.sp.add)((ip.sp.scale)(2.0, u), v)
  assert_true(((ip.inner)(s, u) - 2.0 * (ip.inner)(u, u) - (ip.inner)(v, u)).abs() < 1.0e-12)
  // 實維數
  assert_eq(hvec_real_dim(2), 8)
}

test "quat: Gram–Schmidt 於 H²" {
  let vs : Array[HVec] = [
    { comps: [quat_real(1.0), quat_real(1.0)] },
    { comps: [quat_real(1.0), quat_imag(1.0, 0.0, 0.0)] },
    { comps: [quat_real(2.0), quat_new(1.0, 1.0, 0.0, 0.0)] }, // 依賴
  ]
  let onb = hvec_orthonormal(vs, 1.0e-10)
  assert_eq(onb.length(), 2) // 第三個被丟棄
  for i = 0; i < onb.length(); i = i + 1 {
    for j = 0; j < onb.length(); j = j + 1 {
      let want = if i == j { 1.0 } else { 0.0 }
      assert_true((hvec_inner(onb[i], onb[j]) - want).abs() < 1.0e-9)
    }
  }
}

test "quat: 單位四元數之範數為 1" {
  let q = quat_from_axis_angle([0.0, 1.0, 0.0], 0.77)
  assert_true((quat_norm(q) - 1.0).abs() < 1.0e-12)
  let u = quat_unit(quat_new(1.0, 2.0, 2.0, 0.0))
  assert_true((quat_norm(u) - 1.0).abs() < 1.0e-12)
}

/// 測試輔助：四元數近似相等
fn ceq_q2(a : Quat, b : Quat) -> Bool {
  (a.w - b.w).abs() < 1.0e-12 &&
    (a.x - b.x).abs() < 1.0e-12 &&
    (a.y - b.y).abs() < 1.0e-12 &&
    (a.z - b.z).abs() < 1.0e-12
}
