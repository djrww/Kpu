///|
// 多變量多項式定理之自證

test "mpoly: 加法與正規化" {
  let a = mpoly_add(mpoly_monomial([1, 0], 2.0), mpoly_monomial([0, 1], 3.0))
  let b = mpoly_add(mpoly_monomial([1, 0], -2.0), mpoly_monomial([0, 0], 1.0))
  let s = mpoly_add(a, b)
  // 應為 1 + 3y
  assert_true(euc_eq_eps([mpoly_eval(s, [0.0, 0.0])], [1.0], 1.0e-12))
  assert_true((mpoly_eval(s, [7.0, 2.0]) - 7.0).abs() < 1.0e-12)
  // 同類項合併
  let c = mpoly_add(mpoly_monomial([2, 1], 1.0), mpoly_monomial([2, 1], -1.0))
  assert_eq(c.terms.length(), 0)
}

test "mpoly: 乘法二項式定理 (x+y)³" {
  let x = mpoly_monomial([1, 0], 1.0)
  let y = mpoly_monomial([0, 1], 1.0)
  let s = mpoly_add(x, y)
  let mut p = mpoly_monomial([0, 0], 1.0)
  for _i = 0; _i < 3; _i = _i + 1 {
    p = mpoly_mul(p, s)
  }
  // 係數應為 1,3,3,1
  assert_eq(p.terms.length(), 4)
  let coefs : Array[Double] = []
  for t in p.terms {
    coefs.push(t.1)
  }
  assert_true(euc_eq_eps(coefs, [1.0, 3.0, 3.0, 1.0], 1.0e-12))
}

test "mpoly: 求值與代數一致" {
  let p = mpoly_add(
    mpoly_monomial([2, 0], 3.0),
    mpoly_add(mpoly_monomial([1, 1], -2.0), mpoly_monomial([0, 1], 5.0)),
  )
  // p(2,3) = 3·4 − 2·6 + 5·3 = 15
  assert_true((mpoly_eval(p, [2.0, 3.0]) - 15.0).abs() < 1.0e-12)
  // 乘積求值 = 求值之積
  let q = mpoly_monomial([1, 2], 1.0) // x y²
  assert_true(
    (mpoly_eval(mpoly_mul(p, q), [2.0, 3.0]) -
      mpoly_eval(p, [2.0, 3.0]) * mpoly_eval(q, [2.0, 3.0])).abs() < 1.0e-10,
  )
}

test "mpoly: 偏導乘積法則" {
  let p = mpoly_add(mpoly_monomial([1, 1], 1.0), mpoly_monomial([0, 2], 2.0))
  let q = mpoly_monomial([2, 0], 1.0)
  let lhs = mpoly_partial(mpoly_mul(p, q), 0)
  let rhs = mpoly_add(
    mpoly_mul(mpoly_partial(p, 0), q),
    mpoly_mul(p, mpoly_partial(q, 0)),
  )
  // 於隨機點比較
  assert_true(
    (mpoly_eval(lhs, [1.5, -0.5]) - mpoly_eval(rhs, [1.5, -0.5])).abs() < 1.0e-10,
  )
}

test "mpoly: Clairaut 定理（混合偏導交換）" {
  let p = mpoly_add(
    mpoly_monomial([3, 2], 2.0),
    mpoly_add(mpoly_monomial([1, 4], -1.0), mpoly_monomial([2, 2], 1.0)),
  )
  let dxy = mpoly_partial(mpoly_partial(p, 0), 1)
  let dyx = mpoly_partial(mpoly_partial(p, 1), 0)
  for pt in [[0.7, 1.3], [2.0, -1.0], [0.0, 0.0]] {
    assert_true(
      (mpoly_eval(dxy, pt) - mpoly_eval(dyx, pt)).abs() < 1.0e-10,
    )
  }
}

test "mpoly: 梯度與方向一致性" {
  // p = x² + y²; ∇p(1,2) = (2,4)
  let p = mpoly_add(mpoly_monomial([2, 0], 1.0), mpoly_monomial([0, 2], 1.0))
  let g = mpoly_gradient(p)
  assert_eq(g.length(), 2)
  assert_true((mpoly_eval(g[0], [1.0, 2.0]) - 2.0).abs() < 1.0e-12)
  assert_true((mpoly_eval(g[1], [1.0, 2.0]) - 4.0).abs() < 1.0e-12)
}

test "mpoly: 歐拉齊次函數定理" {
  // 齊次多項式次數 d：Σ xᵢ ∂p/∂xᵢ = d·p
  let p = mpoly_add(
    mpoly_monomial([3, 0], 1.0),
    mpoly_add(mpoly_monomial([2, 1], 2.0), mpoly_monomial([0, 3], -1.0)),
  )
  assert_eq(mpoly_total_degree(p), 3)
  let g = mpoly_gradient(p)
  let pt = [1.5, 2.0]
  let mut s = 0.0
  for i = 0; i < 2; i = i + 1 {
    s = s + pt[i] * mpoly_eval(g[i], pt)
  }
  assert_true((s - 3.0 * mpoly_eval(p, pt)).abs() < 1.0e-10)
}

test "mpoly: 內積正定與雙線性" {
  let p = mpoly_add(mpoly_monomial([1, 0], 2.0), mpoly_monomial([0, 1], 1.0))
  let q = mpoly_add(mpoly_monomial([1, 0], 1.0), mpoly_monomial([2, 0], 3.0))
  // ⟨p,q⟩ = 2·1 = 2（公共單項式只有 x）
  assert_true((mpoly_inner(p, q) - 2.0).abs() < 1.0e-12)
  assert_true(mpoly_inner(p, p) > 0.0)
  assert_true((mpoly_inner(p, q) - mpoly_inner(q, p)).abs() < 1.0e-12)
}

test "mpoly: 單項式計數公式" {
  assert_eq(monomial_count(2, 2), 6)
  assert_eq(monomial_count(3, 3), 20)
  assert_eq(monomial_count(1, 5), 6)
  assert_eq(monomial_count(5, 0), 1)
}

test "mpoly: 齊次化與回代" {
  let p = mpoly_add(
    mpoly_monomial([2, 1], 1.0),
    mpoly_add(mpoly_monomial([1, 0], 3.0), mpoly_monomial([0, 0], -2.0)),
  )
  let h = mpoly_homogenize(p)
  // 全部項總次數 = 3
  for t in h.terms {
    let mut s = 0
    for e in t.0 {
      s = s + e
    }
    assert_eq(s, 3)
  }
  // 於新變量 = 1 時回原值
  assert_true(
    (mpoly_eval(h, [2.0, 3.0, 1.0]) - mpoly_eval(p, [2.0, 3.0])).abs() < 1.0e-12,
  )
}

test "mpoly: 內積空間字典（Gram–Schmidt 可用）" {
  let ip = mpoly_ip(2)
  let basis : Array[MPoly] = [
    mpoly_monomial([0, 0], 1.0),
    mpoly_add(mpoly_monomial([0, 0], 1.0), mpoly_monomial([1, 0], 1.0)),
  ]
  let onb = ip_orthonormal(ip, basis, 1.0e-10)
  assert_eq(onb.length(), 2)
  assert_true((mpoly_inner(onb[0], onb[1])).abs() < 1.0e-9)
  assert_true((mpoly_inner(onb[0], onb[0]) - 1.0).abs() < 1.0e-9)
}
