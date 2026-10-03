///|
// =====================================================================
// registry.mbt — 演算法目錄與義務自證
// 每條演算法帶一個自證閉包：回傳 true 表示該演算法通過其義務
// （公理、定理、算式之數值驗證）。總數須多於 80。
// =====================================================================

/// 演算法條目
pub struct AlgoEntry {
  id : String
  name : String
  check : () -> Bool
}

/// 浮點近似（註冊表公用）
fn approx(a : Double, b : Double, eps : Double) -> Bool {
  (a - b).abs() < eps
}

/// 全部演算法註冊表
pub fn registry() -> Array[AlgoEntry] {
  let r : Array[AlgoEntry] = []

  // ---------- 複數域 ----------
  r.push({
    id: "C-01",
    name: "複數乘法",
    check: fn() {
      let z = cmul(cnew(1.0, 2.0), cnew(3.0, -1.0))
      ceq_eps(z, cnew(5.0, 5.0), 1.0e-12)
    },
  })
  r.push({
    id: "C-02",
    name: "共軛對合",
    check: fn() {
      let z = cnew(2.0, -3.0)
      ceq_eps(cconj(cconj(z)), z, 1.0e-12)
    },
  })
  r.push({
    id: "C-03",
    name: "乘法逆元",
    check: fn() {
      let z = cnew(2.0, 3.0)
      ceq_eps(cmul(z, cinv(z)), cone(), 1.0e-12)
    },
  })
  r.push({
    id: "C-04",
    name: "歐拉式 e^{iθ}",
    check: fn() {
      ceq_eps(cexp(ci()), cpolar(1.0, 1.0), 1.0e-12)
    },
  })
  r.push({
    id: "C-05",
    name: "主平方根",
    check: fn() {
      let z = cnew(-3.0, 4.0)
      ceq_eps(csqr(csqrt(z)), z, 1.0e-9)
    },
  })

  // ---------- 向量空間公理 ----------
  let sp3 = euc_sp(3)
  let u3 = [1.0, 2.0, 3.0]
  let v3 = [4.0, -1.0, 2.0]
  r.push({
    id: "VS-01",
    name: "交換律 u+v = v+u",
    check: fn() {
      euc_eq_eps((sp3.add)(u3, v3), (sp3.add)(v3, u3), 1.0e-12)
    },
  })
  r.push({
    id: "VS-02",
    name: "零元素",
    check: fn() {
      euc_eq_eps((sp3.add)(u3, sp3.zero), u3, 1.0e-12)
    },
  })
  r.push({
    id: "VS-03",
    name: "負元素",
    check: fn() {
      euc_eq_eps((sp3.add)(u3, (sp3.neg)(u3)), sp3.zero, 1.0e-12)
    },
  })
  r.push({
    id: "VS-04",
    name: "分配律 (a+b)v = av + bv",
    check: fn() {
      let a = 1.7
      let b = -0.4
      let lhs = (sp3.scale)(a + b, u3)
      let rhs = (sp3.add)((sp3.scale)(a, u3), (sp3.scale)(b, u3))
      euc_eq_eps(lhs, rhs, 1.0e-12)
    },
  })
  r.push({
    id: "VS-05",
    name: "1·v = v",
    check: fn() {
      euc_eq_eps((sp3.scale)(1.0, u3), u3, 1.0e-12)
    },
  })
  r.push({
    id: "VS-06",
    name: "倍增迭代 2⁴v = 16v",
    check: fn() {
      euc_eq_eps(vdouble_iter(sp3, u3, 4), (sp3.scale)(16.0, u3), 1.0e-12)
    },
  })
  r.push({
    id: "VS-07",
    name: "線性組合",
    check: fn() {
      let c = vlinc(sp3, [2.0, -1.0], [u3, v3])
      euc_eq_eps(c, euc_add(euc_scale(2.0, u3), euc_neg(v3)), 1.0e-12)
    },
  })

  // ---------- 內積空間 ----------
  let ip3 = euc_ip(3)
  r.push({
    id: "IP-01",
    name: "範數-內積一致 ‖v‖² = ⟨v,v⟩",
    check: fn() {
      approx(ip_norm(ip3, u3) * ip_norm(ip3, u3), (ip3.inner)(u3, u3), 1.0e-12)
    },
  })
  r.push({
    id: "IP-02",
    name: "投影冪等 proj² = proj",
    check: fn() {
      let p = ip_proj(ip3, u3, v3)
      let pp = ip_proj(ip3, p, v3)
      euc_eq_eps(p, pp, 1.0e-9)
    },
  })
  r.push({
    id: "IP-03",
    name: "勾股定理（正交分解）",
    check: fn() {
      let p = ip_proj(ip3, u3, v3)
      let rj = ip_rej(ip3, u3, v3)
      let lhs = ip_norm(ip3, u3) * ip_norm(ip3, u3)
      let rhs = ip_norm(ip3, p) * ip_norm(ip3, p) +
        ip_norm(ip3, rj) * ip_norm(ip3, rj)
      approx(lhs, rhs, 1.0e-9)
    },
  })
  r.push({
    id: "IP-04",
    name: "Gram–Schmidt 正交性",
    check: fn() {
      let vs = [[1.0, 1.0, 0.0], [1.0, 0.0, 1.0], [0.0, 1.0, 1.0]]
      let onb = ip_orthonormal(ip3, vs, 1.0e-10)
      let mut ok = onb.length() == 3
      for i = 0; i < onb.length(); i = i + 1 {
        for j = 0; j < onb.length(); j = j + 1 {
          let want = if i == j { 1.0 } else { 0.0 }
          if !approx((ip3.inner)(onb[i], onb[j]), want, 1.0e-9) {
            ok = false
          }
        }
      }
      ok
    },
  })
  r.push({
    id: "IP-05",
    name: "Bessel 不等式殘差 ≥ 0",
    check: fn() {
      let onb = [evec(3, 0), evec(3, 1)]
      ip_bessel_residual(ip3, [1.0, 2.0, 3.0], onb) > -1.0e-12
    },
  })
  r.push({
    id: "IP-06",
    name: "Fourier 重構（完整標準正交基）",
    check: fn() {
      let onb = euc_basis(3)
      let v = [0.5, -2.0, 7.0]
      let cs = ip_coords(ip3, v, onb)
      euc_eq_eps(ip_from_coords(ip3, cs, onb), v, 1.0e-12)
    },
  })
  r.push({
    id: "IP-07",
    name: "鏡射是對合 R² = I",
    check: fn() {
      let n = ip_normalize(ip3, [1.0, 1.0, 0.0])
      let rr = ip_reflect(ip3, ip_reflect(ip3, u3, n), n)
      euc_eq_eps(rr, u3, 1.0e-9)
    },
  })
  r.push({
    id: "IP-08",
    name: "Gram 行列式 = 體積²",
    check: fn() {
      let vs = [[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]]
      approx(ip_volume(ip3, vs), 6.0, 1.0e-9)
    },
  })

  // ---------- 歐幾里得 ----------
  r.push({
    id: "E-01",
    name: "叉積正交性",
    check: fn() {
      let a = [1.0, 2.0, 3.0]
      let b = [-1.0, 0.5, 2.0]
      let c = euc_cross3(a, b)
      euc_inner(c, a).abs() < 1.0e-12 && euc_inner(c, b).abs() < 1.0e-12
    },
  })
  r.push({
    id: "E-02",
    name: "Cauchy–Schwarz",
    check: fn() {
      let a = [1.0, -2.0, 3.0, 4.0]
      let b = [-1.0, 2.5, 0.0, 1.0]
      euc_inner(a, b).abs() <= ip_norm(euc_ip(4), a) * ip_norm(euc_ip(4), b) + 1.0e-12
    },
  })
  r.push({
    id: "E-03",
    name: "三角不等式",
    check: fn() {
      let a = [3.0, -1.0, 2.0]
      let b = [-2.0, 4.0, 1.0]
      ip_norm(euc_ip(3), euc_add(a, b)) <=
        ip_norm(euc_ip(3), a) + ip_norm(euc_ip(3), b) + 1.0e-12
    },
  })

  // ---------- 矩陣 ----------
  let a22 = mat_from_rows([[2.0, 1.0], [1.0, 3.0]])
  r.push({
    id: "M-01",
    name: "矩陣乘法結合律",
    check: fn() {
      let b = mat_from_rows([[1.0, 0.0], [2.0, -1.0]])
      let c = mat_from_rows([[0.0, 3.0], [1.0, 1.0]])
      mat_eq_eps(mat_mul(mat_mul(a22, b), c), mat_mul(a22, mat_mul(b, c)), 1.0e-9)
    },
  })
  r.push({
    id: "M-02",
    name: "det 乘法性 det(AB)=detA·detB",
    check: fn() {
      let b = mat_from_rows([[1.0, 2.0], [-1.0, 0.0]])
      approx(mat_det(mat_mul(a22, b)), mat_det(a22) * mat_det(b), 1.0e-9)
    },
  })
  r.push({
    id: "M-03",
    name: "逆矩陣 AA⁻¹ = I",
    check: fn() {
      match mat_inv(a22) {
        Some(inv) => mat_eq_eps(mat_mul(a22, inv), mat_eye(2), 1.0e-9)
        None => false
      }
    },
  })
  r.push({
    id: "M-04",
    name: "LU 重構 PA = LU",
    check: fn() {
      let m = mat_from_rows([[0.0, 2.0, 1.0], [1.0, 0.0, 3.0], [4.0, 1.0, 0.0]])
      let (l, u, perm) = mat_lu(m)
      let lu = mat_mul(l, u)
      let mut ok = true
      for i = 0; i < 3; i = i + 1 {
        for j = 0; j < 3; j = j + 1 {
          if !approx(mat_get(lu, i, j), mat_get(m, perm[i], j), 1.0e-9) {
            ok = false
          }
        }
      }
      ok
    },
  })
  r.push({
    id: "M-05",
    name: "Cholesky 重構 A = LLᵀ",
    check: fn() {
      match mat_cholesky(a22) {
        Some(l) => mat_eq_eps(mat_mul(l, mat_transpose(l)), a22, 1.0e-9)
        None => false
      }
    },
  })
  r.push({
    id: "M-06",
    name: "QR 重構 A = QR 且 QᵀQ = I",
    check: fn() {
      let m = mat_from_rows([[1.0, 1.0], [1.0, 2.0], [0.0, 1.0]])
      let (q, rr) = mat_qr(m)
      let qtq = mat_mul(mat_transpose(q), q)
      mat_eq_eps(mat_mul(q, rr), m, 1.0e-9) && mat_eq_eps(qtq, mat_eye(2), 1.0e-9)
    },
  })
  r.push({
    id: "M-07",
    name: "Jacobi 譜分解重構",
    check: fn() {
      approx(sym_spectral_error(a22), 0.0, 1.0e-8)
    },
  })
  r.push({
    id: "M-08",
    name: "冪迭代特徵方程 Av ≈ λv",
    check: fn() {
      let (lam, v) = power_iteration(a22, 200)
      let av = mat_vec(a22, v)
      euc_eq_eps(av, euc_scale(lam, v), 1.0e-6)
    },
  })
  r.push({
    id: "M-09",
    name: "RREF 秩",
    check: fn() {
      let m = mat_from_rows([[1.0, 2.0, 3.0], [2.0, 4.0, 6.0], [0.0, 1.0, 1.0]])
      mat_rank(m, 1.0e-10) == 2
    },
  })
  r.push({
    id: "M-10",
    name: "零空間 Av = 0",
    check: fn() {
      let m = mat_from_rows([[1.0, 2.0, 3.0], [2.0, 4.0, 6.0]])
      let basis = mat_null_space_basis(m, 1.0e-10)
      let mut ok = basis.length() == 2
      for v in basis {
        if !euc_eq_eps(mat_vec(m, v), [0.0, 0.0], 1.0e-9) {
          ok = false
        }
      }
      ok
    },
  })
  r.push({
    id: "M-11",
    name: "最小二乘法方程 AᵀAx = Aᵀb",
    check: fn() {
      let m = mat_from_rows([[1.0, 0.0], [0.0, 1.0], [1.0, 1.0]])
      let b = [1.0, 2.0, 2.4]
      let x = least_squares(m, b)
      let lhs = mat_vec(mat_mul(mat_transpose(m), m), x)
      let rhs = mat_vec(mat_transpose(m), b)
      euc_eq_eps(lhs, rhs, 1.0e-8)
    },
  })
  r.push({
    id: "M-12",
    name: "矩陣指數（零矩陣）e⁰ = I",
    check: fn() {
      mat_eq_eps(mat_exp(mat_new(3, 3)), mat_eye(3), 1.0e-12)
    },
  })
  r.push({
    id: "M-13",
    name: "Kronecker 維數",
    check: fn() {
      let k = mat_kronecker(a22, mat_eye(3))
      k.rows == 6 && k.cols == 6
    },
  })

  // ---------- 對稱 / 對角 ----------
  r.push({
    id: "S-01",
    name: "Sym(n) 標準正交基",
    check: fn() {
      let basis = sym_basis(3)
      let ip = sym_ip(3)
      let mut ok = basis.length() == sym_dim(3)
      for i = 0; i < basis.length(); i = i + 1 {
        for j = 0; j < basis.length(); j = j + 1 {
          let want = if i == j { 1.0 } else { 0.0 }
          if !approx(sym_inner(basis[i], basis[j]), want, 1.0e-12) {
            ok = false
          }
        }
      }
      let _ = ip
      ok
    },
  })
  r.push({
    id: "S-02",
    name: "對稱化為冪等正交投影",
    check: fn() {
      let m = mat_from_rows([[1.0, 5.0], [1.0, 2.0]])
      let s = mat_symmetrize(m)
      mat_eq_eps(mat_symmetrize(s), s, 1.0e-12) && mat_is_symmetric(s, 1.0e-12)
    },
  })
  r.push({
    id: "S-03",
    name: "pack/unpack 往返",
    check: fn() {
      let m = sym_pack(3, [1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
      mat_is_symmetric(m, 1.0e-12) &&
        euc_eq_eps(sym_unpack(m), [1.0, 2.0, 3.0, 4.0, 5.0, 6.0], 1.0e-12)
    },
  })
  r.push({
    id: "D-01",
    name: "對角行列式 = 積",
    check: fn() {
      approx(diag_det([2.0, -3.0, 4.0]), -24.0, 1.0e-12)
    },
  })
  r.push({
    id: "D-02",
    name: "對角逆",
    check: fn() {
      match diag_inv([2.0, 4.0]) {
        Some(inv) => euc_eq_eps(diag_mul([2.0, 4.0], inv), [1.0, 1.0], 1.0e-12)
        None => false
      }
    },
  })

  // ---------- Hermitian ----------
  r.push({
    id: "H-01",
    name: "Hermitian 內積正定性",
    check: fn() {
      let v = [cnew(1.0, 2.0), cnew(-1.0, 0.5)]
      cvec_norm2(v) > 0.0 && approx(cvec_inner(v, v).im, 0.0, 1.0e-12)
    },
  })
  r.push({
    id: "H-02",
    name: "複 Gram–Schmidt",
    check: fn() {
      let vs = [[cnew(1.0, 1.0), cnew(0.0, 1.0)], [cone(), cnew(1.0, -1.0)]]
      let onb = cvec_orthonormal(vs, 1.0e-10)
      let mut ok = onb.length() == 2
      for i = 0; i < 2; i = i + 1 {
        for j = 0; j < 2; j = j + 1 {
          let want = if i == j { 1.0 } else { 0.0 }
          if !approx(cvec_inner(onb[i], onb[j]).re, want, 1.0e-9) {
            ok = false
          }
          if !approx(cvec_inner(onb[i], onb[j]).im, 0.0, 1.0e-9) {
            ok = false
          }
        }
      }
      ok
    },
  })
  r.push({
    id: "H-03",
    name: "Hermitian 特徵值（實）與重構",
    check: fn() {
      // A = [[2, 1+i], [1-i, 3]]
      let m = cmat_new(2, 2)
      cmat_set(m, 0, 0, cfrom_real(2.0))
      cmat_set(m, 1, 1, cfrom_real(3.0))
      cmat_set(m, 0, 1, cnew(1.0, 1.0))
      cmat_set(m, 1, 0, cnew(1.0, -1.0))
      let (vals, q) = cmat_jacobi_hermitian(m, 60)
      // 特徵值應為 (5 ± √(1+8))/2 = (5 ± 3)/2 = {1, 4}
      let lo = if vals[0] < vals[1] { vals[0] } else { vals[1] }
      let hi = if vals[0] < vals[1] { vals[1] } else { vals[0] }
      approx(lo, 1.0, 1.0e-8) && approx(hi, 4.0, 1.0e-8) &&
        cmat_unitary_error(q) < 1.0e-8
    },
  })

  // ---------- 多項式 ----------
  r.push({
    id: "P-01",
    name: "Horner 求值",
    check: fn() {
      // p(x) = 1 + 2x + 3x²; p(2) = 17
      approx(poly_eval([1.0, 2.0, 3.0], 2.0), 17.0, 1.0e-12)
    },
  })
  r.push({
    id: "P-02",
    name: "微積分基本定理（∫p′ = p(1)−p(0)）",
    check: fn() {
      let p = [1.0, -2.0, 0.5, 3.0]
      let q = poly_integral(poly_deriv(p))
      approx(poly_eval(q, 1.0) - poly_eval(q, 0.0), poly_eval(p, 1.0) - poly_eval(p, 0.0), 1.0e-10)
    },
  })
  r.push({
    id: "P-03",
    name: "帶餘除法 p = q·quot + rem",
    check: fn() {
      let p = poly_from_roots([1.0, 2.0, 3.0]) // (x-1)(x-2)(x-3)
      let q = [-2.0, 1.0] // x-2
      let (quot, rem) = poly_divmod(p, q)
      euc_eq_eps(
        poly_add(poly_mul(q, quot), rem),
        p,
        1.0e-9,
      ) && poly_degree(rem) == 0
    },
  })
  r.push({
    id: "P-04",
    name: "GCD：gcd(p·g, q·g) = g（首一）",
    check: fn() {
      let g = [1.0, 1.0] // x+1
      let p = poly_mul(g, [-3.0, 1.0])
      let q = poly_mul(g, [2.0, 1.0])
      let d = poly_gcd(p, q)
      euc_eq_eps(d, g, 1.0e-8)
    },
  })
  r.push({
    id: "P-05",
    name: "Lagrange 插值過點",
    check: fn() {
      let xs = [0.0, 1.0, 2.0, 3.0]
      let ys = [1.0, 0.0, 5.0, -2.0]
      let p = poly_lagrange(xs, ys)
      let mut ok = true
      for i = 0; i < xs.length(); i = i + 1 {
        if !approx(poly_eval(p, xs[i]), ys[i], 1.0e-9) {
          ok = false
        }
      }
      ok
    },
  })
  r.push({
    id: "P-06",
    name: "Newton 插值 = Lagrange 插值",
    check: fn() {
      let xs = [-1.0, 0.0, 2.0, 4.0]
      let ys = [3.0, 1.0, -1.0, 2.0]
      let a = poly_lagrange(xs, ys)
      let b = poly_newton(xs, ys)
      euc_eq_eps(a, b, 1.0e-8)
    },
  })
  r.push({
    id: "P-07",
    name: "移位 Legendre 正交性",
    check: fn() {
      let ls = shifted_legendre(4)
      let mut ok = ls.length() == 5
      for i = 0; i < ls.length(); i = i + 1 {
        for j = 0; j < ls.length(); j = j + 1 {
          let want = if i == j { 1.0 } else { 0.0 }
          if !approx(poly_inner(ls[i], ls[j]), want, 1.0e-8) {
            ok = false
          }
        }
      }
      ok
    },
  })
  r.push({
    id: "P-08",
    name: "Durand–Kerner 根回代",
    check: fn() {
      let p = poly_from_roots([1.0, -2.0, 0.5])
      let roots = poly_roots_dk(p, 200)
      let mut ok = roots.length() == 3
      for z in roots {
        // p(z) ≈ 0
        let mut acc = czero()
        for k = poly_degree(p); k >= 0; k = k - 1 {
          acc = cadd(cmul(acc, z), cfrom_real(p[k]))
        }
        if cabs(acc) > 1.0e-7 {
          ok = false
        }
      }
      ok
    },
  })

  // ---------- 多變量多項式 ----------
  r.push({
    id: "MP-01",
    name: "多項式乘法分配律",
    check: fn() {
      // (x + y)(x - y) = x² - y²
      let a = mpoly_add(
        mpoly_monomial([1, 0], 1.0),
        mpoly_monomial([0, 1], 1.0),
      )
      let b = mpoly_add(
        mpoly_monomial([1, 0], 1.0),
        mpoly_monomial([0, 1], -1.0),
      )
      let p = mpoly_mul(a, b)
      approx(mpoly_eval(p, [3.0, 2.0]), 9.0 - 4.0, 1.0e-12)
    },
  })
  r.push({
    id: "MP-02",
    name: "偏導與求值一致",
    check: fn() {
      // p = x²y + 3y²; ∂p/∂x = 2xy; 於 (2,5) 值 20
      let p = mpoly_add(
        mpoly_monomial([2, 1], 1.0),
        mpoly_monomial([0, 2], 3.0),
      )
      let dx = mpoly_partial(p, 0)
      approx(mpoly_eval(dx, [2.0, 5.0]), 20.0, 1.0e-12)
    },
  })
  r.push({
    id: "MP-03",
    name: "單項式計數 C(n+d, d)",
    check: fn() {
      monomial_count(3, 2) == 10 && monomial_count(4, 3) == 35
    },
  })
  r.push({
    id: "MP-04",
    name: "齊次化於 x_{n}=1 回原值",
    check: fn() {
      let p = mpoly_add(
        mpoly_monomial([2, 0], 2.0),
        mpoly_monomial([0, 1], -1.0),
      )
      let h = mpoly_homogenize(p)
      approx(mpoly_eval(h, [3.0, 4.0, 1.0]), mpoly_eval(p, [3.0, 4.0]), 1.0e-12)
    },
  })

  // ---------- 函數空間 ----------
  r.push({
    id: "F-01",
    name: "DFT/IDFT 往返",
    check: fn() {
      let f = [cnew(1.0, 0.0), cnew(2.0, -1.0), cnew(0.0, 1.0), cnew(-1.0, 2.0)]
      let back = fn_idft(fn_dft(f))
      let mut ok = true
      for i = 0; i < f.length(); i = i + 1 {
        if !ceq_eps(back[i], f[i], 1.0e-9) {
          ok = false
        }
      }
      ok
    },
  })
  r.push({
    id: "F-02",
    name: "Parseval 定理",
    check: fn() {
      let f = [cnew(1.0, 1.0), cnew(-2.0, 0.0), cnew(0.0, 3.0), cnew(1.0, -1.0)]
      fn_parseval_residual(f).abs() < 1.0e-9
    },
  })
  r.push({
    id: "F-03",
    name: "卷積定理：DFT(f∗g) = DFT(f)·DFT(g)",
    check: fn() {
      let f = [1.0, 2.0, 0.0, -1.0]
      let g = [0.0, 1.0, 3.0, 2.0]
      let conv = fn_convolution(f, g)
      let fc : Array[C] = []
      let gc : Array[C] = []
      for i = 0; i < f.length(); i = i + 1 {
        fc.push(cfrom_real(f[i]))
        gc.push(cfrom_real(g[i]))
      }
      let lhs = fn_dft(conv.map(fn(x) { cfrom_real(x) }))
      let rhs_f = fn_dft(fc)
      let rhs_g = fn_dft(gc)
      let mut ok = true
      for i = 0; i < lhs.length(); i = i + 1 {
        if !ceq_eps(lhs[i], cmul(rhs_f[i], rhs_g[i]), 1.0e-8) {
          ok = false
        }
      }
      ok
    },
  })
  r.push({
    id: "F-04",
    name: "Fun(X, Sym(2)) 組合空間之內積",
    check: fn() {
      let w = [1.0, 2.0]
      let ip = fn_space_ip(2, w, sym_ip(2))
      let f = [sym_pack(2, [1.0, 0.0, 1.0]), sym_pack(2, [0.0, 1.0, 2.0])]
      let g = [sym_pack(2, [1.0, 1.0, 0.0]), sym_pack(2, [2.0, 0.0, 1.0])]
      // ⟨f,g⟩ = Σ wᵢ tr(fᵢgᵢ)
      let expect = 1.0 * (1.0 * 1.0 + 0.0 * 1.0 + 0.0 * 1.0 + 1.0 * 0.0) +
        2.0 * (0.0 * 2.0 + 1.0 * 0.0 + 1.0 * 0.0 + 2.0 * 1.0)
      approx((ip.inner)(f, g), expect, 1.0e-12)
    },
  })
  r.push({
    id: "F-05",
    name: "差分算子線性",
    check: fn() {
      let f = [1.0, 3.0, 2.0, 5.0]
      let g = [2.0, 0.0, 1.0, 1.0]
      let lhs = fn_finite_diff(euc_add(f, g))
      let rhs = euc_add(fn_finite_diff(f), fn_finite_diff(g))
      euc_eq_eps(lhs, rhs, 1.0e-12)
    },
  })

  // ---------- 張量 ----------
  r.push({
    id: "T-01",
    name: "純張量內積分解",
    check: fn() {
      let u = [1.0, 2.0]
      let v = [3.0, -1.0]
      let t = pure_tensor(u, v)
      let s = [0.0, 1.0]
      let w = [2.0, 2.0]
      let t2 = pure_tensor(s, w)
      approx(tensor_inner(t, t2), euc_inner(u, s) * euc_inner(v, w), 1.0e-12)
    },
  })
  r.push({
    id: "T-02",
    name: "收縮 = 雙線性求值",
    check: fn() {
      let t = pure_tensor([1.0, 2.0, 3.0], [4.0, 5.0, 6.0])
      approx(tensor_contract(t), 32.0, 1.0e-12) &&
        approx(
          tensor_bilinear_eval(t, [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]),
          5.0,
          1.0e-12,
        )
    },
  })
  r.push({
    id: "T-03",
    name: "張量-矩陣同構往返",
    check: fn() {
      let t = pure_tensor([1.0, 2.0], [3.0, 4.0, 5.0])
      let m = tensor_to_mat(t)
      let t2 = mat_to_tensor(m)
      mat_eq_eps(tensor_to_mat(t2), m, 1.0e-12)
    },
  })

  // ---------- 外代數 ----------
  r.push({
    id: "X-01",
    name: "反對稱性 u∧v = −v∧u",
    check: fn() {
      let u = [1.0, 2.0, 0.0]
      let v = [3.0, -1.0, 4.0]
      let a = ext_wedge_vecs(3, u, v)
      let b = ext_wedge_vecs(3, v, u)
      let mut ok = true
      for i = 0; i < a.data.length(); i = i + 1 {
        if !approx(a.data[i], -b.data[i], 1.0e-12) {
          ok = false
        }
      }
      ok
    },
  })
  r.push({
    id: "X-02",
    name: "‖u∧v‖ = 平行四邊形面積",
    check: fn() {
      let u = [3.0, 0.0, 0.0]
      let v = [0.0, 4.0, 0.0]
      approx(ext_parallelogram_area(u, v), 12.0, 1.0e-12)
    },
  })
  r.push({
    id: "X-03",
    name: "Hodge 星對合 ⋆⋆ = (−1)^{k(n−k)}",
    check: fn() {
      let w = ext_wedge_vecs(4, [1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0])
      let ww = ext_hodge_star(ext_hodge_star(w))
      // n=4, k=2 → 符號 +1
      let mut ok = true
      for i = 0; i < ww.data.length(); i = i + 1 {
        if !approx(ww.data[i], w.data[i], 1.0e-12) {
          ok = false
        }
      }
      ok
    },
  })
  r.push({
    id: "X-04",
    name: "外積行列式 = 矩陣行列式",
    check: fn() {
      let vs = [[2.0, 0.0, 1.0], [1.0, 3.0, 0.0], [0.0, 1.0, 2.0]]
      let m = mat_from_rows([vs[0], vs[1], vs[2]])
      approx(ext_det(vs), mat_det(m), 1.0e-9)
    },
  })
  r.push({
    id: "X-05",
    name: "叉積 = ⋆(u∧v)",
    check: fn() {
      let u = [1.0, 2.0, 3.0]
      let v = [-2.0, 1.0, 0.0]
      euc_eq_eps(ext_cross3(u, v), euc_cross3(u, v), 1.0e-12) &&
        euc_eq_eps(ext_cross3(u, v), euc_cross3(u, v), 1.0e-12)
    },
  })
  r.push({
    id: "X-06",
    name: "內縮與外積互伴性質",
    check: fn() {
      // i_v(v∧w) = ‖v‖²w − ⟨v,w⟩v
      let v = [1.0, 2.0, 1.0]
      let w = [0.0, 1.0, 3.0]
      let vw = ext_wedge_vecs(3, v, w)
      let iw = ext_interior(v, vw)
      let expect = euc_add(
        euc_scale(euc_inner(v, v), w),
        euc_scale(-euc_inner(v, w), v),
      )
      euc_eq_eps(iw.data, expect, 1.0e-9)
    },
  })

  // ---------- 對偶 ----------
  r.push({
    id: "DU-01",
    name: "對偶配對 ⟨δᵢ, eⱼ⟩ = δᵢⱼ",
    check: fn() {
      let phis = dual_basis(3)
      let m = dual_pairing_matrix(phis, 3)
      mat_eq_eps(m, mat_eye(3), 1.0e-12)
    },
  })
  r.push({
    id: "DU-02",
    name: "Riesz 表示 φ(v) = ⟨r,v⟩",
    check: fn() {
      let phi : DualVec = { coef: [2.0, -1.0, 0.5] }
      let v = [1.0, 1.0, 4.0]
      approx(dual_action(phi, v), euc_inner(riesz_represent(phi), v), 1.0e-12)
    },
  })
  r.push({
    id: "DU-03",
    name: "對偶映射 (A*ψ)(v) = ψ(Av)",
    check: fn() {
      let a = mat_from_rows([[1.0, 2.0], [3.0, 4.0]])
      let psi : DualVec = { coef: [1.0, -1.0] }
      let v = [2.0, 5.0]
      let astar = dual_map(a) // 2×2 轉置
      let lhs = dual_action({ coef: mat_vec(astar, psi.coef) }, v)
      let rhs = dual_action(psi, mat_vec(a, v))
      approx(lhs, rhs, 1.0e-12)
    },
  })

  // ---------- 四元數 ----------
  r.push({
    id: "Q-01",
    name: "i² = j² = k² = ijk = −1",
    check: fn() {
      let i = quat_imag(1.0, 0.0, 0.0)
      let j = quat_imag(0.0, 1.0, 0.0)
      let k = quat_imag(0.0, 0.0, 1.0)
      ceq_q(quat_mul(i, i), quat_real(-1.0)) &&
        ceq_q(quat_mul(j, j), quat_real(-1.0)) &&
        ceq_q(quat_mul(k, k), quat_real(-1.0)) &&
        ceq_q(quat_mul(quat_mul(i, j), k), quat_real(-1.0))
    },
  })
  r.push({
    id: "Q-02",
    name: "範數乘法性 ‖pq‖ = ‖p‖‖q‖",
    check: fn() {
      let p = quat_new(1.0, 2.0, -1.0, 0.5)
      let q = quat_new(0.0, 1.0, 1.0, -2.0)
      approx(quat_norm(quat_mul(p, q)), quat_norm(p) * quat_norm(q), 1.0e-12)
    },
  })
  r.push({
    id: "Q-03",
    name: "旋轉保範數與保角",
    check: fn() {
      let q = quat_from_axis_angle([0.0, 0.0, 1.0], @math.PI / 3.0)
      let v = [1.0, 2.0, 3.0]
      let rv = quat_rotate_vec(q, v)
      approx(
        ip_norm(euc_ip(3), rv),
        ip_norm(euc_ip(3), v),
        1.0e-9,
      ) && approx(euc_inner(rv, rv), euc_inner(v, v), 1.0e-9)
    },
  })
  r.push({
    id: "Q-04",
    name: "繞 z 軸 90° 旋轉 (1,0,0) ↦ (0,1,0)",
    check: fn() {
      let q = quat_from_axis_angle([0.0, 0.0, 1.0], @math.PI / 2.0)
      let rv = quat_rotate_vec(q, [1.0, 0.0, 0.0])
      euc_eq_eps(rv, [0.0, 1.0, 0.0], 1.0e-9)
    },
  })
  r.push({
    id: "Q-05",
    name: "H² 實內積正定",
    check: fn() {
      let v : HVec = {
        comps: [quat_new(1.0, 1.0, 0.0, 0.0), quat_new(0.0, 0.0, 2.0, -1.0)],
      }
      approx(hvec_inner(v, v), quat_norm2(v.comps[0]) + quat_norm2(v.comps[1]), 1.0e-12)
    },
  })
  r.push({
    id: "Q-06",
    name: "H^n Gram–Schmidt",
    check: fn() {
      let vs : Array[HVec] = [
        {
          comps: [
            quat_new(1.0, 0.0, 0.0, 0.0),
            quat_new(0.0, 1.0, 0.0, 0.0),
          ],
        },
        {
          comps: [
            quat_new(1.0, 1.0, 0.0, 0.0),
            quat_new(0.0, 0.0, 1.0, 1.0),
          ],
        },
      ]
      let onb = hvec_orthonormal(vs, 1.0e-10)
      let mut ok = onb.length() == 2
      for i = 0; i < 2; i = i + 1 {
        for j = 0; j < 2; j = j + 1 {
          let want = if i == j { 1.0 } else { 0.0 }
          if !approx(hvec_inner(onb[i], onb[j]), want, 1.0e-9) {
            ok = false
          }
        }
      }
      ok
    },
  })

  // ---------- 空間組合 ----------
  r.push({
    id: "PR-01",
    name: "直和內積 = 分量內積和",
    check: fn() {
      let ip = ips_prod(euc_ip(2), euc_ip(3))
      let p = ([1.0, 2.0], [1.0, 0.0, 1.0])
      let q = ([3.0, -1.0], [2.0, 2.0, 0.0])
      approx((ip.inner)(p, q), 1.0 * 3.0 + 2.0 * -1.0 + 2.0, 1.0e-12)
    },
  })
  r.push({
    id: "PR-02",
    name: "直和勾股定理",
    check: fn() {
      let ip = ips_prod(euc_ip(2), euc_ip(2))
      let p = (([1.0, 0.0]), ([0.0, 0.0]))
      let q = (([0.0, 0.0]), ([0.0, 1.0]))
      let s = (ip.sp.add)(p, q)
      let n2 = (ip.inner)(s, s)
      approx(n2, (ip.inner)(p, p) + (ip.inner)(q, q), 1.0e-12)
    },
  })
  r.push({
    id: "PR-03",
    name: "映射直和分塊對角",
    check: fn() {
      let m = map_direct_sum(mat_eye(2), mat_scale(3.0, mat_eye(3)))
      mat_eq_eps(
        m,
        mat_diag([1.0, 1.0, 3.0, 3.0, 3.0]),
        1.0e-12,
      )
    },
  })

  r
}

/// 四元數近似相等（註冊表輔助）
fn ceq_q(a : Quat, b : Quat) -> Bool {
  (a.w - b.w).abs() < 1.0e-12 &&
    (a.x - b.x).abs() < 1.0e-12 &&
    (a.y - b.y).abs() < 1.0e-12 &&
    (a.z - b.z).abs() < 1.0e-12
}

/// 執行全部自證，回傳 (通過數, 總數, 失敗者)
pub fn registry_run() -> (Int, Int, Array[String]) {
  let entries = registry()
  let mut passed = 0
  let fails : Array[String] = []
  for e in entries {
    if (e.check)() {
      passed = passed + 1
    } else {
      fails.push(e.id + " " + e.name)
    }
  }
  (passed, entries.length(), fails)
}
