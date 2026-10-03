///|
// =====================================================================
// fspace.mbt — 有限集合上函數空間 Fun(X, F)（X = {0,…,n−1}）
// 加權 ℓ² 內積、DFT、卷積、差分；可與任意內積空間 V 組合成
// Fun(X, V)（見 F-11）
// =====================================================================

/// F-1 加權 ℓ² 內積 Σ wᵢ f(i) g(i)
pub fn fn_inner_weighted(
  f : Array[Double],
  g : Array[Double],
  w : Array[Double],
) -> Double {
  let mut s = 0.0
  for i = 0; i < f.length(); i = i + 1 {
    s = s + w[i] * f[i] * g[i]
  }
  s
}

/// F-2 標準（等權）內積
pub fn fn_inner(f : Array[Double], g : Array[Double]) -> Double {
  let w = Array::make(f.length(), 1.0)
  fn_inner_weighted(f, g, w)
}

/// F-3 Fun(X, R) 內積空間字典（加權）
pub fn fn_ip(n : Int, w : Array[Double]) -> IPSp[Array[Double]] {
  ips_make(euc_sp(n), fn(a, b) { fn_inner_weighted(a, b, w) })
}

/// F-4 delta 基 {δ_0,…,δ_{n−1}}
pub fn fn_delta_basis(n : Int) -> Array[Array[Double]] {
  euc_basis(n)
}

/// F-5 逐點乘法（代數結構）
pub fn fn_pointwise_mul(f : Array[Double], g : Array[Double]) -> Array[Double] {
  euc_hadamard(f, g)
}

/// F-6 總和泛函
pub fn fn_sum(f : Array[Double]) -> Double {
  euc_sum(f)
}

/// F-7 循環前向差分 (Δf)(i) = f(i+1) − f(i)（指標模 n）
pub fn fn_finite_diff(f : Array[Double]) -> Array[Double] {
  let n = f.length()
  let out = Array::make(n, 0.0)
  for i = 0; i < n; i = i + 1 {
    out[i] = f[(i + 1) % n] - f[i]
  }
  out
}

/// F-8 離散傅里葉變換：F[k] = Σⱼ f[j] e^{−2πi·jk/n}
pub fn fn_dft(f : Array[C]) -> Array[C] {
  let n = f.length()
  let out : Array[C] = []
  for k = 0; k < n; k = k + 1 {
    let mut s = czero()
    for j = 0; j < n; j = j + 1 {
      let ang = -2.0 * @math.PI * j.to_double() * k.to_double() /
        n.to_double()
      s = cadd(s, cmul(f[j], cpolar(1.0, ang)))
    }
    out.push(s)
  }
  out
}

/// F-9 逆 DFT：f[j] = (1/n) Σₖ F[k] e^{2πi·jk/n}
pub fn fn_idft(fc : Array[C]) -> Array[C] {
  let n = fc.length()
  let out : Array[C] = []
  for j = 0; j < n; j = j + 1 {
    let mut s = czero()
    for k = 0; k < n; k = k + 1 {
      let ang = 2.0 * @math.PI * j.to_double() * k.to_double() /
        n.to_double()
      s = cadd(s, cmul(fc[k], cpolar(1.0, ang)))
    }
    out.push(cscale(1.0 / n.to_double(), s))
  }
  out
}

/// F-10 循環卷積 (f ∗ g)(k) = Σⱼ f[j] g(k−j mod n)
pub fn fn_convolution(f : Array[Double], g : Array[Double]) -> Array[Double] {
  let n = f.length()
  let out = Array::make(n, 0.0)
  for k = 0; k < n; k = k + 1 {
    let mut s = 0.0
    for j = 0; j < n; j = j + 1 {
      s = s + f[j] * g[(k - j + n) % n]
    }
    out[k] = s
  }
  out
}

/// F-11 函數空間組合子：Fun(X, V)，V 為任意內積空間
/// 內積 = Σᵢ wᵢ ⟨f(i), g(i)⟩_V（空間 × 空間之組合）
pub fn[V] fn_space_ip(
  n : Int,
  w : Array[Double],
  ip : IPSp[V],
) -> IPSp[Array[V]] {
  let zero : Array[V] = []
  for _i = 0; _i < n; _i = _i + 1 {
    zero.push(ip.sp.zero)
  }
  ips_make(
    {
      zero,
      add: fn(f, g) {
        let out : Array[V] = []
        for i = 0; i < f.length(); i = i + 1 {
          out.push((ip.sp.add)(f[i], g[i]))
        }
        out
      },
      neg: fn(f) {
        let out : Array[V] = []
        for v in f {
          out.push((ip.sp.neg)(v))
        }
        out
      },
      scale: fn(c, f) {
        let out : Array[V] = []
        for v in f {
          out.push((ip.sp.scale)(c, v))
        }
        out
      },
    },
    fn(f, g) {
      let mut s = 0.0
      for i = 0; i < f.length(); i = i + 1 {
        s = s + w[i] * (ip.inner)(f[i], g[i])
      }
      s
    },
  )
}

/// F-12 Parseval 殘差 ‖f‖² − (1/n)‖F‖²（對 DFT 必為 0）
pub fn fn_parseval_residual(f : Array[C]) -> Double {
  let n = f.length()
  let fc = fn_dft(f)
  let mut lhs = 0.0
  for x in f {
    lhs = lhs + cnorm2(x)
  }
  let mut rhs = 0.0
  for x in fc {
    rhs = rhs + cnorm2(x)
  }
  lhs - rhs / n.to_double()
}
