///|
// 內積空間公理與定理之自證

test "ipspace: 內積公理（對稱、雙線性、正定）" {
  let ip = euc_ip(3)
  let u = [1.0, -2.0, 2.0]
  let v = [3.0, 0.0, 1.0]
  let w = [1.0, 1.0, 1.0]
  let a = 2.0
  let b = -0.5
  // 對稱
  assert_true(((ip.inner)(u, v) - (ip.inner)(v, u)).abs() < 1.0e-12)
  // 第一變量線性
  let lhs = (ip.inner)((ip.sp.add)((ip.sp.scale)(a, u), (ip.sp.scale)(b, v)), w)
  let rhs = a * (ip.inner)(u, w) + b * (ip.inner)(v, w)
  assert_true((lhs - rhs).abs() < 1.0e-12)
  // 正定
  assert_true((ip.inner)(u, u) > 0.0)
  assert_true((ip.inner)(ip.sp.zero, ip.sp.zero) == 0.0)
}

test "ipspace: Cauchy–Schwarz 不等式" {
  let ip = euc_ip(4)
  let u = [2.0, -1.0, 3.0, 0.5]
  let v = [1.0, 2.0, -2.0, 1.0]
  let c = (ip.inner)(u, v)
  assert_true(c.abs() <= ip_norm(ip, u) * ip_norm(ip, v) + 1.0e-12)
  // 等式情形：u, v 線性相關
  let v2 = (ip.sp.scale)(1.7, u)
  let c2 = (ip.inner)(u, v2)
  assert_true((c2.abs() - ip_norm(ip, u) * ip_norm(ip, v2)).abs() < 1.0e-10)
}

test "ipspace: 三角不等式" {
  let ip = euc_ip(3)
  let u = [1.0, 2.0, 3.0]
  let v = [-4.0, 1.0, 0.5]
  assert_true(
    ip_norm(ip, (ip.sp.add)(u, v)) <= ip_norm(ip, u) + ip_norm(ip, v) + 1.0e-12,
  )
}

test "ipspace: 餘弦範圍與角度" {
  let ip = euc_ip(3)
  let u = [1.0, 0.0, 0.0]
  let v = [0.0, 1.0, 0.0]
  assert_true((ip_cos(ip, u, v)).abs() < 1.0e-12)
  assert_true((ip_angle(ip, u, v) - @math.PI / 2.0).abs() < 1.0e-12)
  let w = [-1.0, 0.0, 0.0]
  assert_true((ip_angle(ip, u, w) - @math.PI).abs() < 1.0e-12)
  assert_true((ip_angle(ip, u, u)).abs() < 1.0e-12)
}

test "ipspace: 正交判定與單位化" {
  let ip = euc_ip(3)
  let u = [1.0, 1.0, 0.0]
  let v = [1.0, -1.0, 0.0]
  assert_true(ip_is_orthogonal(ip, u, v, 1.0e-12))
  let n = ip_normalize(ip, u)
  assert_true((ip_norm(ip, n) - 1.0).abs() < 1.0e-12)
  // 零向量單位化回零
  let z = ip_normalize(ip, ip.sp.zero)
  assert_true(ip_norm(ip, z) == 0.0)
}

test "ipspace: 投影正交性與冪等性" {
  let ip = euc_ip(3)
  let u = [3.0, 1.0, 2.0]
  let v = [1.0, 2.0, 1.0]
  let p = ip_proj(ip, u, v)
  let r = ip_rej(ip, u, v)
  // p + r = u
  assert_true(euc_eq_eps((ip.sp.add)(p, r), u, 1.0e-10))
  // p ⊥ r
  assert_true((ip.inner)(p, r).abs() < 1.0e-10)
  // 冪等
  assert_true(euc_eq_eps(ip_proj(ip, p, v), p, 1.0e-10))
  // 係數具體值
  let c = ip_proj_coeff(ip, u, v)
  assert_true((c - 7.0 / 6.0).abs() < 1.0e-12)
}

test "ipspace: Gram–Schmidt 保生成空間" {
  let ip = euc_ip(3)
  let vs : Array[Array[Double]] = [
    [1.0, 1.0, 0.0],
    [2.0, 1.0, 1.0],
    [3.0, 2.0, 1.0], // 依賴於前兩者
  ]
  let orth = ip_gram_schmidt(ip, vs, 1.0e-10)
  assert_eq(orth.length(), 2) // 第三個向量被丟棄
  let onb = ip_orthonormal(ip, vs, 1.0e-10)
  // 標準正交
  for i = 0; i < onb.length(); i = i + 1 {
    for j = 0; j < onb.length(); j = j + 1 {
      let want = if i == j { 1.0 } else { 0.0 }
      assert_true(((ip.inner)(onb[i], onb[j]) - want).abs() < 1.0e-9)
    }
  }
  // 投影回原向量 = 原向量（在生成空間內）
  for v in vs {
    let p = ip_proj_subspace(ip, v, onb)
    assert_true(euc_eq_eps(p, v, 1.0e-8))
  }
}

test "ipspace: Fourier 座標重構與 Parseval" {
  let ip = euc_ip(3)
  let v = [1.0, -2.0, 3.0]
  let onb = euc_basis(3)
  let cs = ip_coords(ip, v, onb)
  assert_true(euc_eq_eps(cs, v, 1.0e-12))
  assert_true(euc_eq_eps(ip_from_coords(ip, cs, onb), v, 1.0e-12))
  // Parseval：‖v‖² = Σ|⟨v,eᵢ⟩|²（完整基時殘差為 0）
  assert_true(ip_bessel_residual(ip, v, onb).abs() < 1.0e-10)
  // Bessel：不完整基時殘差 ≥ 0
  assert_true(ip_bessel_residual(ip, v, [evec(3, 0)]) > 0.0)
}

test "ipspace: 鏡射是等距且對合" {
  let ip = euc_ip(3)
  let n = ip_normalize(ip, [1.0, 2.0, -1.0])
  let v = [3.0, -1.0, 2.0]
  let rv = ip_reflect(ip, v, n)
  // 等距
  assert_true((ip_norm(ip, rv) - ip_norm(ip, v)).abs() < 1.0e-10)
  // 對合
  assert_true(euc_eq_eps(ip_reflect(ip, rv, n), v, 1.0e-10))
  // v − Rv 平行於 n
  let diff = vsub(ip.sp, v, rv)
  let c = ip_proj_coeff(ip, diff, n)
  assert_true(euc_eq_eps(diff, euc_scale(c, n), 1.0e-10))
}

test "ipspace: Gram 行列式與體積" {
  let ip = euc_ip(3)
  // 長方體 2×3×4
  let vs : Array[Array[Double]] = [
    [2.0, 0.0, 0.0],
    [0.0, 3.0, 0.0],
    [0.0, 0.0, 4.0],
  ]
  let g = ip_gram_matrix(ip, vs)
  assert_true(euc_eq_eps(g, [4.0, 0.0, 0.0, 0.0, 9.0, 0.0, 0.0, 0.0, 16.0], 1.0e-12))
  assert_true((ip_volume(ip, vs) - 24.0).abs() < 1.0e-9)
  // 線性相關 → 體積 0
  let dep : Array[Array[Double]] = [[1.0, 0.0, 0.0], [2.0, 0.0, 0.0]]
  assert_true(ip_volume(ip, dep).abs() < 1.0e-9)
}

test "ipspace: 到子空間距離" {
  let ip = euc_ip(3)
  let onb = [evec(3, 0), evec(3, 1)] // xy 平面
  let v = [1.0, 2.0, 5.0]
  assert_true((ip_dist_subspace(ip, v, onb) - 5.0).abs() < 1.0e-10)
  let p = ip_proj_subspace(ip, v, onb)
  assert_true(euc_eq_eps(p, [1.0, 2.0, 0.0], 1.0e-10))
}

test "ipspace: 加權內積尺度變換" {
  let ip = euc_ip(2)
  let ip2 = ips_scale(ip, 3.0)
  let u = [1.0, 2.0]
  let v = [3.0, 1.0]
  assert_true(((ip2.inner)(u, v) - 3.0 * (ip.inner)(u, v)).abs() < 1.0e-12)
  assert_true((ip_norm(ip2, u) - (3.0).sqrt() * ip_norm(ip, u)).abs() < 1.0e-12)
}

test "ipspace: 多項式空間內積精確值" {
  // ∫₀¹ x·x dx = 1/3
  let ip = poly_ip(2)
  assert_true(((ip.inner)([0.0, 1.0], [0.0, 1.0]) - 1.0 / 3.0).abs() < 1.0e-12)
  // ∫₀¹ 1·x dx = 1/2
  assert_true(((ip.inner)([1.0], [0.0, 1.0]) - 0.5).abs() < 1.0e-12)
}

test "ipspace: 對稱矩陣空間之內積 = Frobenius" {
  let ip = sym_ip(2)
  let a = sym_pack(2, [1.0, 2.0, 3.0])
  let b = sym_pack(2, [4.0, 5.0, 6.0])
  // tr(ab): a=[[1,2],[2,3]], b=[[4,5],[5,6]] → tr = 1*4+2*5+2*5+3*6 = 42
  assert_true(((ip.inner)(a, b) - 42.0).abs() < 1.0e-12)
}
