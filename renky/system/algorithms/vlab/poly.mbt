///|
// =====================================================================
// poly.mbt — 單變量多項式空間 P_n（係數陣列，索引 = 次數）
// 內積：⟨p,q⟩ = ∫₀¹ p(x)q(x) dx（精確計算，非數值積分）
// =====================================================================

/// 去掉高次零係數
fn poly_trim(p : Array[Double]) -> Array[Double] {
  let mut n = p.length()
  while n > 1 && p[n - 1] == 0.0 {
    n = n - 1
  }
  let out : Array[Double] = []
  for i = 0; i < n; i = i + 1 {
    out.push(p[i])
  }
  out
}

/// P-1 零多項式
pub fn poly_zero() -> Array[Double] {
  [0.0]
}

/// P-2 常數多項式
pub fn poly_const(c : Double) -> Array[Double] {
  [c]
}

/// P-3 次數（零多項式定為 0）
pub fn poly_degree(p : Array[Double]) -> Int {
  poly_trim(p).length() - 1
}

/// P-4 首項係數
pub fn poly_leading(p : Array[Double]) -> Double {
  let t = poly_trim(p)
  t[t.length() - 1]
}

/// P-5 加法
pub fn poly_add(a : Array[Double], b : Array[Double]) -> Array[Double] {
  let n = if a.length() > b.length() { a.length() } else { b.length() }
  let out = Array::make(n, 0.0)
  for i = 0; i < a.length(); i = i + 1 {
    out[i] = out[i] + a[i]
  }
  for i = 0; i < b.length(); i = i + 1 {
    out[i] = out[i] + b[i]
  }
  poly_trim(out)
}

/// P-6 減法
pub fn poly_sub(a : Array[Double], b : Array[Double]) -> Array[Double] {
  poly_add(a, poly_scale(-1.0, b))
}

/// P-7 純量乘法
pub fn poly_scale(c : Double, p : Array[Double]) -> Array[Double] {
  let out : Array[Double] = []
  for x in p {
    out.push(c * x)
  }
  poly_trim(out)
}

/// P-8 乘法（卷積）
pub fn poly_mul(a : Array[Double], b : Array[Double]) -> Array[Double] {
  let out = Array::make(a.length() + b.length() - 1, 0.0)
  for i = 0; i < a.length(); i = i + 1 {
    for j = 0; j < b.length(); j = j + 1 {
      out[i + j] = out[i + j] + a[i] * b[j]
    }
  }
  poly_trim(out)
}

/// P-9 導數
pub fn poly_deriv(p : Array[Double]) -> Array[Double] {
  if p.length() <= 1 {
    poly_zero()
  } else {
    let out : Array[Double] = []
    for i = 1; i < p.length(); i = i + 1 {
      out.push(i.to_double() * p[i])
    }
    poly_trim(out)
  }
}

/// P-10 不定積分（常數項為 0）
pub fn poly_integral(p : Array[Double]) -> Array[Double] {
  let out = Array::make(p.length() + 1, 0.0)
  for i = 0; i < p.length(); i = i + 1 {
    out[i + 1] = p[i] / (i.to_double() + 1.0)
  }
  poly_trim(out)
}

/// P-11 Horner 求值
pub fn poly_eval(p : Array[Double], x : Double) -> Double {
  let mut acc = 0.0
  for i = p.length() - 1; i >= 0; i = i - 1 {
    acc = acc * x + p[i]
  }
  acc
}

/// P-12 L²[0,1] 內積（精確）：Σ aᵢbⱼ/(i+j+1)
pub fn poly_inner(a : Array[Double], b : Array[Double]) -> Double {
  let mut s = 0.0
  for i = 0; i < a.length(); i = i + 1 {
    for j = 0; j < b.length(); j = j + 1 {
      s = s + a[i] * b[j] / (i.to_double() + j.to_double() + 1.0)
    }
  }
  s
}

/// P-13 P_n 內積空間字典（次數 ≤ n）
pub fn poly_ip(n : Int) -> IPSp[Array[Double]] {
  ips_make(
    {
      zero: poly_zero(),
      add: poly_add,
      neg: fn(p) { poly_scale(-1.0, p) },
      scale: poly_scale,
    },
    poly_inner,
  )
}

/// P-14 移位 Legendre 多項式：Gram–Schmidt 於 {1,x,…,x^n}
pub fn shifted_legendre(n : Int) -> Array[Array[Double]] {
  let ip = poly_ip(n)
  let basis : Array[Array[Double]] = []
  for k = 0; k <= n; k = k + 1 {
    let m = Array::make(k + 1, 0.0)
    m[k] = 1.0
    basis.push(m)
  }
  ip_orthonormal(ip, basis, 1.0e-10)
}

/// P-15 帶餘除法：p = q·quot + rem（deg rem < deg q）
pub fn poly_divmod(
  p : Array[Double],
  q : Array[Double],
) -> (Array[Double], Array[Double]) {
  let mut rem = poly_trim(p.copy())
  let qt = poly_trim(q)
  let dq = poly_degree(qt)
  let lc = poly_leading(qt)
  let mut quot = poly_zero()
  while poly_degree(rem) >= dq && poly_leading(rem) != 0.0 {
    let d = poly_degree(rem) - dq
    let c = poly_leading(rem) / lc
    let term = Array::make(d + 1, 0.0)
    term[d] = c
    quot = poly_add(quot, term)
    rem = poly_sub(rem, poly_mul(term, qt))
  }
  (poly_trim(quot), poly_trim(rem))
}

/// P-16 餘式
pub fn poly_mod(p : Array[Double], q : Array[Double]) -> Array[Double] {
  poly_divmod(p, q).1
}

/// P-17 最大公因式（歐幾里得演算法，首一化）
pub fn poly_gcd(a : Array[Double], b : Array[Double]) -> Array[Double] {
  let mut x = poly_trim(a.copy())
  let mut y = poly_trim(b.copy())
  while poly_leading(y) != 0.0 {
    let r = poly_mod(x, y)
    x = y
    y = r
  }
  // 首一化
  let lc = poly_leading(x)
  poly_scale(1.0 / lc, x)
}

/// P-18 Lagrange 插值：過點 (xs[i], ys[i])
pub fn poly_lagrange(xs : Array[Double], ys : Array[Double]) -> Array[Double] {
  let n = xs.length()
  let mut out = poly_zero()
  for i = 0; i < n; i = i + 1 {
    // L_i(x) = Π_{j≠i} (x − x_j)/(x_i − x_j)
    let mut li = poly_const(1.0)
    let mut denom = 1.0
    for j = 0; j < n; j = j + 1 {
      if j != i {
        li = poly_mul(li, [-xs[j], 1.0])
        denom = denom * (xs[i] - xs[j])
      }
    }
    li = poly_scale(ys[i] / denom, li)
    out = poly_add(out, li)
  }
  out
}

/// P-19 Newton 插值（均差表）
pub fn poly_newton(xs : Array[Double], ys : Array[Double]) -> Array[Double] {
  let n = xs.length()
  let dd = ys.copy()
  let coeffs : Array[Double] = []
  for k = 0; k < n; k = k + 1 {
    coeffs.push(dd[k])
    for i = n - 1; i > k; i = i - 1 {
      dd[i] = (dd[i] - dd[i - 1]) / (xs[i] - xs[i - k - 1])
    }
  }
  // 展開 Π (x − x_i)
  let mut out = poly_zero()
  for k = 0; k < n; k = k + 1 {
    let mut basis = poly_const(1.0)
    for i = 0; i < k; i = i + 1 {
      basis = poly_mul(basis, [-xs[i], 1.0])
    }
    out = poly_add(out, poly_scale(coeffs[k], basis))
  }
  out
}

/// P-20 複合 p∘q
pub fn poly_compose(p : Array[Double], q : Array[Double]) -> Array[Double] {
  let mut acc = poly_zero()
  for i = p.length() - 1; i >= 0; i = i - 1 {
    acc = poly_add(poly_mul(acc, q), poly_const(p[i]))
  }
  acc
}

/// P-21 由根構造 Π(x − rᵢ)
pub fn poly_from_roots(rs : Array[Double]) -> Array[Double] {
  let mut out = poly_const(1.0)
  for r in rs {
    out = poly_mul(out, [-r, 1.0])
  }
  out
}

/// P-22 二次方程實/複根（回傳 (重數意義下) 兩複數根）
pub fn poly_roots_quadratic(a : Double, b : Double, c : Double) -> (C, C) {
  let disc = b * b - 4.0 * a * c
  if disc >= 0.0 {
    let sd = disc.sqrt()
    (cfrom_real((-b + sd) / (2.0 * a)), cfrom_real((-b - sd) / (2.0 * a)))
  } else {
    let sd = (-disc).sqrt()
    (
      cnew(-b / (2.0 * a), sd / (2.0 * a)),
      cnew(-b / (2.0 * a), -sd / (2.0 * a)),
    )
  }
}

/// P-23 Durand–Kerner 求全部複根（首一化後迭代）
pub fn poly_roots_dk(p : Array[Double], iters : Int) -> Array[C] {
  let t = poly_trim(p.copy())
  let n = poly_degree(t)
  let lc = poly_leading(t)
  let roots : Array[C] = []
  if n == 0 {
    return roots
  }
  // 初值：圓上分散點
  for k = 0; k < n; k = k + 1 {
    let ang = 2.0 * @math.PI * k.to_double() / n.to_double() + 0.4
    roots.push(cpolar(1.0 + 0.1 * k.to_double(), ang))
  }
  for _it = 0; _it < iters; _it = _it + 1 {
    for i = 0; i < n; i = i + 1 {
      // p(z_i)/lc
      let mut num = czero()
      for k = poly_degree(t); k >= 0; k = k - 1 {
        num = cadd(cmul(num, roots[i]), cfrom_real(t[k] / lc))
      }
      let mut den = cone()
      for j = 0; j < n; j = j + 1 {
        if j != i {
          den = cmul(den, csub(roots[i], roots[j]))
        }
      }
      roots[i] = csub(roots[i], cdiv(num, den))
    }
  }
  roots
}

/// P-24 多項式空間之維數（次數 ≤ n 者為 n+1）
pub fn poly_dim(n : Int) -> Int {
  n + 1
}
