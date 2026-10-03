///|
// 多項式空間定理之自證

test "poly: 環公理（分配律等）" {
  let p = [1.0, 2.0] // 1+2x
  let q = [-1.0, 0.0, 3.0] // -1+3x²
  let r = [2.0, 1.0, 1.0]
  // 分配律 p(q+r) = pq + pr
  let lhs = poly_mul(p, poly_add(q, r))
  let rhs = poly_add(poly_mul(p, q), poly_mul(p, r))
  assert_true(euc_eq_eps(lhs, rhs, 1.0e-12))
  // 結合律
  let lhs2 = poly_mul(poly_mul(p, q), r)
  let rhs2 = poly_mul(p, poly_mul(q, r))
  assert_true(euc_eq_eps(lhs2, rhs2, 1.0e-12))
}

test "poly: 導數乘積法則 (pq)' = p'q + pq'" {
  let p = [1.0, 1.0, 1.0]
  let q = [2.0, 0.0, -1.0, 3.0]
  let lhs = poly_deriv(poly_mul(p, q))
  let rhs = poly_add(poly_mul(poly_deriv(p), q), poly_mul(p, poly_deriv(q)))
  assert_true(euc_eq_eps(lhs, rhs, 1.0e-12))
}

test "poly: 積分線性與微積分基本定理" {
  let p = [3.0, -1.0, 2.0]
  let q = [1.0, 1.0]
  let a = 2.0
  let b = -1.0
  let lhs = poly_integral(poly_add(poly_scale(a, p), poly_scale(b, q)))
  let rhs = poly_add(poly_scale(a, poly_integral(p)), poly_scale(b, poly_integral(q)))
  assert_true(euc_eq_eps(lhs, rhs, 1.0e-12))
  // ∫₀¹ p = F(1) − F(0)
  let f = poly_integral(p)
  let expect = 3.0 - 0.5 + 2.0 / 3.0
  assert_true((poly_eval(f, 1.0) - expect).abs() < 1.0e-12)
}

test "poly: 內積對稱與正定" {
  let p = [1.0, 2.0, -1.0]
  let q = [0.0, 1.0, 1.0]
  assert_true((poly_inner(p, q) - poly_inner(q, p)).abs() < 1.0e-12)
  assert_true(poly_inner(p, p) > 0.0)
  // 具體值：∫₀¹ (1+2x−x²)² dx = 2 + 2/3 + 1/5
  assert_true((poly_inner(p, p) - (2.0 + 2.0 / 3.0 + 1.0 / 5.0)).abs() < 1.0e-12)
}

test "poly: 帶餘除法唯一性與次數" {
  let p = [1.0, 0.0, 0.0, 0.0, 1.0] // x⁴+1
  let q = [-1.0, 0.0, 1.0] // x²−1
  let (quot, rem) = poly_divmod(p, q)
  assert_true(euc_eq_eps(poly_add(poly_mul(q, quot), rem), p, 1.0e-9))
  assert_true(poly_degree(rem) < poly_degree(q))
  // x⁴+1 = (x²−1)(x²+1) + 2
  assert_true(euc_eq_eps(quot, [1.0, 0.0, 1.0], 1.0e-9))
  assert_true(euc_eq_eps(rem, [2.0], 1.0e-9))
}

test "poly: GCD 之 Bezout 性質（餘式為零）" {
  let g = [2.0, 1.0] // x+2
  let p = poly_mul(g, [1.0, 1.0])
  let q = poly_mul(g, [-1.0, 0.0, 1.0])
  let d = poly_gcd(p, q)
  // d 整除 p 與 q
  assert_true(euc_eq_eps(poly_mod(p, d), [0.0], 1.0e-7))
  assert_true(euc_eq_eps(poly_mod(q, d), [0.0], 1.0e-7))
  // 首一
  assert_true((poly_leading(d) - 1.0).abs() < 1.0e-7)
}

test "poly: 插值唯一性（Lagrange = Newton）" {
  let xs = [0.0, 1.0, 3.0, 4.0]
  let ys = [2.0, -1.0, 0.0, 5.0]
  let a = poly_lagrange(xs, ys)
  let b = poly_newton(xs, ys)
  assert_true(euc_eq_eps(a, b, 1.0e-8))
  assert_true(poly_degree(a) <= 3)
}

test "poly: 複合函數求值一致" {
  let p = [1.0, 0.0, 1.0] // x²+1
  let q = [1.0, 1.0] // x+1
  let pq = poly_compose(p, q) // (x+1)²+1
  let x = 2.5
  assert_true(
    (poly_eval(pq, x) - poly_eval(p, poly_eval(q, x))).abs() < 1.0e-10,
  )
}

test "poly: 由根構造與 Vieta" {
  let rs = [1.0, 2.0, 3.0]
  let p = poly_from_roots(rs) // x³−6x²+11x−6
  assert_true(euc_eq_eps(p, [-6.0, 11.0, -6.0, 1.0], 1.0e-12))
  // 根回代
  for r in rs {
    assert_true(poly_eval(p, r).abs() < 1.0e-12)
  }
}

test "poly: 二次公式（實根與複根）" {
  // x²−5x+6 = (x−2)(x−3)
  let (r1, r2) = poly_roots_quadratic(1.0, -5.0, 6.0)
  assert_true(ceq_eps(r1, cfrom_real(3.0), 1.0e-12))
  assert_true(ceq_eps(r2, cfrom_real(2.0), 1.0e-12))
  // x²+1 = 0 → ±i
  let (s1, s2) = poly_roots_quadratic(1.0, 0.0, 1.0)
  assert_true(ceq_eps(s1, ci(), 1.0e-12))
  assert_true(ceq_eps(s2, cneg(ci()), 1.0e-12))
}

test "poly: Durand–Kerner 與二次公式一致" {
  let p = [6.0, -11.0, 6.0, -1.0] // −(x−1)(x−2)(x−3)
  let roots = poly_roots_dk(p, 300)
  assert_eq(roots.length(), 3)
  let expected = [cfrom_real(1.0), cfrom_real(2.0), cfrom_real(3.0)]
  // 每個根都該在 expected 附近（排序不定）
  for e in expected {
    let mut found = false
    for r in roots {
      if ceq_eps(r, e, 1.0e-6) {
        found = true
      }
    }
    assert_true(found)
  }
}

test "poly: 移位 Legendre 前三項與標準公式一致" {
  let ls = shifted_legendre(2)
  // P₀ = 1 → 正規化後 1
  assert_true(euc_eq_eps(ls[0], [1.0], 1.0e-9))
  // P₁ = 2x−1 正規化：‖2x−1‖² = ∫(2x−1)² = 1/3 → ×√3
  let expect1 = [-1.0 * (3.0).sqrt(), 2.0 * (3.0).sqrt()]
  let ok1 = euc_eq_eps(ls[1], expect1, 1.0e-8) || euc_eq_eps(ls[1], poly_scale(-1.0, expect1), 1.0e-8)
  assert_true(ok1)
  // 正交
  assert_true(poly_inner(ls[0], ls[2]).abs() < 1.0e-8)
  assert_true(poly_inner(ls[1], ls[2]).abs() < 1.0e-8)
}
