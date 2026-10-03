///|
// 空間組合定理之自證

test "product: 直和之向量空間公理" {
  let sp = vec_prod(euc_sp(2), euc_sp(3))
  let u = ([1.0, 1.0], [1.0, 0.0, 0.0])
  let v = ([0.0, 1.0], [0.0, 1.0, 0.0])
  let w = ([1.0, 0.0], [0.0, 0.0, 1.0])
  // 結合律
  let lhs = (sp.add)((sp.add)(u, v), w)
  let rhs = (sp.add)(u, (sp.add)(v, w))
  assert_true(euc_eq_eps(lhs.0, rhs.0, 1.0e-12))
  assert_true(euc_eq_eps(lhs.1, rhs.1, 1.0e-12))
  // 零元、負元
  let z = (sp.add)(u, (sp.neg)(u))
  assert_true(euc_eq_eps(z.0, [0.0, 0.0], 1.0e-12))
  assert_true(euc_eq_eps(z.1, [0.0, 0.0, 0.0], 1.0e-12))
}

test "product: 直和內積之勾股定理（ι₁V ⊥ ι₂W）" {
  let a = euc_ip(2)
  let b = euc_ip(2)
  let ip = ips_prod(a, b)
  let v = [3.0, 4.0]
  let w = [1.0, 2.0]
  let i1 = inj1(a, b, v)
  let i2 = inj2(a, b, w)
  assert_true((ip.inner)(i1, i2).abs() < 1.0e-12)
  let s = (ip.sp.add)(i1, i2)
  assert_true(
    ((ip.inner)(s, s) - (ip.inner)(i1, i1) - (ip.inner)(i2, i2)).abs() < 1.0e-12,
  )
  // 投影回分量
  assert_true(euc_eq_eps(prod_proj1(s), v, 1.0e-12))
  assert_true(euc_eq_eps(prod_proj2(s), w, 1.0e-12))
}

test "product: 加權直和" {
  let ip = ips_prod_weighted(euc_ip(1), euc_ip(1), 2.0, 3.0)
  let p = ([1.0], [1.0])
  assert_true(((ip.inner)(p, p) - 5.0).abs() < 1.0e-12)
}

test "product: 對角嵌入與投影" {
  let v = [1.0, -2.0]
  let d = diag_embed(v)
  assert_true(euc_eq_eps(d.0, v, 1.0e-12))
  assert_true(euc_eq_eps(d.1, v, 1.0e-12))
  // Δ 為等距嵌入（範數平方變兩倍）
  let ip = ips_prod(euc_ip(2), euc_ip(2))
  assert_true(
    ((ip.inner)(d, d) - 2.0 * euc_inner(v, v)).abs() < 1.0e-12,
  )
}

test "product: 映射直和之譜 = 譜之並" {
  let a = mat_from_rows([[2.0, 0.0], [0.0, 3.0]])
  let b = mat_from_rows([[5.0, 0.0], [0.0, 7.0]])
  let m = map_direct_sum(a, b)
  let (vals, _) = jacobi_eigen(m, 60)
  let sorted = vals.copy()
  sorted.sort()
  assert_true(euc_eq_eps(sorted, [2.0, 3.0, 5.0, 7.0], 1.0e-8))
  // det = det A · det B
  assert_true((mat_det(m) - mat_det(a) * mat_det(b)).abs() < 1.0e-9)
}

test "product: 迭代組合 (R² ⊕ R²) ⊕ R² ≅ Fun({0,1,2}, R²)" {
  let ip_left = ips_prod(ips_prod(euc_ip(2), euc_ip(2)), euc_ip(2))
  let ip_fn = fn_space_ip(3, [1.0, 1.0, 1.0], euc_ip(2))
  let v1 = [1.0, 2.0]
  let v2 = [3.0, 4.0]
  let v3 = [5.0, 6.0]
  let lhs = (ip_left.inner)(((v1, v2), v3), ((v1, v2), v3))
  let f : Array[Array[Double]] = [v1, v2, v3]
  let rhs = (ip_fn.inner)(f, f)
  assert_true((lhs - rhs).abs() < 1.0e-12)
}

test "product: Gram–Schmidt 於直和空間" {
  let ip = ips_prod(euc_ip(2), euc_ip(2))
  let vs : Array[(Array[Double], Array[Double])] = [
    ([1.0, 0.0], [0.0, 0.0]),
    ([1.0, 1.0], [0.0, 0.0]),
    ([0.0, 0.0], [1.0, 0.0]),
  ]
  let onb = ip_orthonormal(ip, vs, 1.0e-10)
  assert_eq(onb.length(), 3)
  for i = 0; i < onb.length(); i = i + 1 {
    for j = 0; j < onb.length(); j = j + 1 {
      let want = if i == j { 1.0 } else { 0.0 }
      assert_true(((ip.inner)(onb[i], onb[j]) - want).abs() < 1.0e-9)
    }
  }
}
