///|
// 有限集合上函數空間定理之自證

test "fspace: 加權內積公理" {
  let w = [1.0, 2.0, 3.0]
  let f = [1.0, 1.0, 1.0]
  let g = [1.0, 2.0, 3.0]
  assert_true((fn_inner_weighted(f, g, w) - (1.0 + 4.0 + 9.0)).abs() < 1.0e-12)
  // 對稱
  assert_true(
    (fn_inner_weighted(f, g, w) - fn_inner_weighted(g, f, w)).abs() < 1.0e-12,
  )
  // 正定
  assert_true(fn_inner_weighted(g, g, w) > 0.0)
}

test "fspace: delta 基之正交性（等權）" {
  let n = 4
  let basis = fn_delta_basis(n)
  for i = 0; i < n; i = i + 1 {
    for j = 0; j < n; j = j + 1 {
      let want = if i == j { 1.0 } else { 0.0 }
      assert_true((fn_inner(basis[i], basis[j]) - want).abs() < 1.0e-12)
    }
  }
}

test "fspace: 逐點積之代數律" {
  let f = [1.0, 2.0, 3.0]
  let g = [2.0, 0.0, -1.0]
  let h = [1.0, 1.0, 1.0]
  let lhs = fn_pointwise_mul(f, fn_pointwise_mul(g, h))
  let rhs = fn_pointwise_mul(fn_pointwise_mul(f, g), h)
  assert_true(euc_eq_eps(lhs, rhs, 1.0e-12))
  // 與純量積相容
  assert_true(euc_eq_eps(fn_pointwise_mul(f, h), f, 1.0e-12))
}

test "fspace: 差分算子與和式（離散微積分基本定理）" {
  let f = [1.0, 3.0, 2.0, 5.0]
  let df = fn_finite_diff(f)
  // 循環差分之和 = 0（伸縮和）
  assert_true(fn_sum(df).abs() < 1.0e-12)
  // 具體值
  assert_true(euc_eq_eps(df, [2.0, -1.0, 3.0, -4.0], 1.0e-12))
}

test "fspace: DFT 線性" {
  let f = [cnew(1.0, 0.0), cnew(0.0, 1.0), cnew(-1.0, 2.0), cnew(2.0, 1.0)]
  let g = [cnew(0.0, 1.0), cnew(1.0, 0.0), cnew(1.0, -1.0), cnew(0.0, 2.0)]
  let lhs = fn_dft(cvec_add(f, g))
  let rhs = cvec_add(fn_dft(f), fn_dft(g))
  for i = 0; i < lhs.length(); i = i + 1 {
    assert_true(ceq_eps(lhs[i], rhs[i], 1.0e-9))
  }
}

test "fspace: DFT 於常數函數 = n·δ₀" {
  let f : Array[C] = []
  for _i = 0; _i < 5; _i = _i + 1 {
    f.push(cone())
  }
  let fc = fn_dft(f)
  assert_true(ceq_eps(fc[0], cfrom_real(5.0), 1.0e-9))
  for k = 1; k < 5; k = k + 1 {
    assert_true(ceq_eps(fc[k], czero(), 1.0e-9))
  }
}

test "fspace: 循環卷積之交換律與結合律" {
  let f = [1.0, 2.0, 3.0, 4.0]
  let g = [0.0, 1.0, -1.0, 2.0]
  let h = [1.0, 0.0, 2.0, 1.0]
  assert_true(euc_eq_eps(fn_convolution(f, g), fn_convolution(g, f), 1.0e-10))
  let lhs = fn_convolution(fn_convolution(f, g), h)
  let rhs = fn_convolution(f, fn_convolution(g, h))
  assert_true(euc_eq_eps(lhs, rhs, 1.0e-10))
  // 與 delta 卷積 = 恒等
  let d = [1.0, 0.0, 0.0, 0.0]
  assert_true(euc_eq_eps(fn_convolution(f, d), f, 1.0e-12))
}

test "fspace: 卷積具體值（循環）" {
  let f = [1.0, 2.0]
  let g = [3.0, 4.0]
  // (f∗g)[0] = f0g0 + f1g_{-1 mod 2} = 3 + 8 = 11
  // (f∗g)[1] = f0g1 + f1g0 = 4 + 6 = 10
  assert_true(euc_eq_eps(fn_convolution(f, g), [11.0, 10.0], 1.0e-12))
}

test "fspace: Parseval 定理（一般情形）" {
  let f = [cnew(1.0, 2.0), cnew(-3.0, 0.0), cnew(2.0, -1.0)]
  assert_true(fn_parseval_residual(f).abs() < 1.0e-9)
}

test "fspace: DFT 矩陣為酉（Plancherel）" {
  // ⟨DFT f, DFT g⟩ = n⟨f, g⟩
  let f = [cnew(1.0, 0.0), cnew(0.0, 2.0), cnew(1.0, 1.0)]
  let g = [cnew(2.0, 1.0), cnew(1.0, 0.0), cnew(-1.0, 2.0)]
  let fc = fn_dft(f)
  let gc = fn_dft(g)
  let mut lhs = czero()
  for i = 0; i < fc.length(); i = i + 1 {
    lhs = cadd(lhs, cmul(cconj(fc[i]), gc[i]))
  }
  let mut rhs = czero()
  for i = 0; i < f.length(); i = i + 1 {
    rhs = cadd(rhs, cmul(cconj(f[i]), g[i]))
  }
  assert_true(ceq_eps(lhs, cscale(3.0, rhs), 1.0e-9))
}

test "fspace: Fun(X,V) 組合空間公理（V = 對稱矩陣）" {
  let w = [1.0, 1.0, 1.0]
  let ip = fn_space_ip(3, w, sym_ip(2))
  let e11 = sym_pack(2, [1.0, 0.0, 0.0])
  let e22 = sym_pack(2, [0.0, 0.0, 1.0])
  let f : Array[Mat] = [e11, e22, e11]
  let g : Array[Mat] = [e22, e11, e22]
  // 正交
  assert_true((ip.inner)(f, g).abs() < 1.0e-12)
  // 正定
  assert_true((ip.inner)(f, f) - 3.0 < 1.0e-12 && (ip.inner)(f, f) > 0.0)
  // 向量空間公理
  let s = (ip.sp.add)(f, g)
  let _ = s
  assert_true((ip.inner)((ip.sp.add)(f, (ip.sp.neg)(f)), f) == 0.0)
}

test "fspace: Fun(X,V) 與直和同構（內積一致）" {
  // Fun({0,1}, R²) ≅ R² ⊕ R²
  let w = [1.0, 1.0]
  let ipf = fn_space_ip(2, w, euc_ip(2))
  let ipd = ips_prod(euc_ip(2), euc_ip(2))
  let f = [[1.0, 2.0], [3.0, 4.0]]
  let g = [[5.0, 6.0], [7.0, 8.0]]
  let a = (ipf.inner)(f, g)
  let b = (ipd.inner)((f[0], f[1]), (g[0], g[1]))
  assert_true((a - b).abs() < 1.0e-12)
}
