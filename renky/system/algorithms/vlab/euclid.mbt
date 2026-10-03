///|
// =====================================================================
// euclid.mbt — 標準歐幾里得空間 R^n（Array[Double] 實作）
// =====================================================================

/// E-1 零向量
pub fn euc_zero(n : Int) -> Array[Double] {
  Array::make(n, 0.0)
}

/// E-2 逐點加法
pub fn euc_add(a : Array[Double], b : Array[Double]) -> Array[Double] {
  let out = Array::make(a.length(), 0.0)
  for i = 0; i < a.length(); i = i + 1 {
    out[i] = a[i] + b[i]
  }
  out
}

/// E-3 負元
pub fn euc_neg(a : Array[Double]) -> Array[Double] {
  let out = Array::make(a.length(), 0.0)
  for i = 0; i < a.length(); i = i + 1 {
    out[i] = -a[i]
  }
  out
}

/// E-4 純量乘法
pub fn euc_scale(c : Double, a : Array[Double]) -> Array[Double] {
  let out = Array::make(a.length(), 0.0)
  for i = 0; i < a.length(); i = i + 1 {
    out[i] = c * a[i]
  }
  out
}

/// E-5 標準內積（點積）Σ aᵢbᵢ
pub fn euc_inner(a : Array[Double], b : Array[Double]) -> Double {
  let mut s = 0.0
  for i = 0; i < a.length(); i = i + 1 {
    s = s + a[i] * b[i]
  }
  s
}

/// E-6 R^n 向量空間字典
pub fn euc_sp(n : Int) -> VecSp[Array[Double]] {
  {
    zero: euc_zero(n),
    add: euc_add,
    neg: euc_neg,
    scale: euc_scale,
  }
}

/// E-7 R^n 內積空間字典
pub fn euc_ip(n : Int) -> IPSp[Array[Double]] {
  ips_make(euc_sp(n), euc_inner)
}

/// E-8 標準基 {e_0, ..., e_{n-1}}
pub fn euc_basis(n : Int) -> Array[Array[Double]] {
  let out : Array[Array[Double]] = []
  for i = 0; i < n; i = i + 1 {
    out.push(evec(n, i))
  }
  out
}

/// E-9 第 i 個標準基向量
pub fn evec(n : Int, i : Int) -> Array[Double] {
  let v = euc_zero(n)
  v[i] = 1.0
  v
}

/// E-10 三維叉積
pub fn euc_cross3(a : Array[Double], b : Array[Double]) -> Array[Double] {
  [
    a[1] * b[2] - a[2] * b[1],
    a[2] * b[0] - a[0] * b[2],
    a[0] * b[1] - a[1] * b[0],
  ]
}

/// E-11 Hadamard 逐點積
pub fn euc_hadamard(a : Array[Double], b : Array[Double]) -> Array[Double] {
  let out = Array::make(a.length(), 0.0)
  for i = 0; i < a.length(); i = i + 1 {
    out[i] = a[i] * b[i]
  }
  out
}

/// E-12 分量總和（線性泛函）
pub fn euc_sum(a : Array[Double]) -> Double {
  let mut s = 0.0
  for x in a {
    s = s + x
  }
  s
}

/// E-13 均值
pub fn euc_mean(a : Array[Double]) -> Double {
  euc_sum(a) / a.length().to_double()
}

/// E-14 近似相等
pub fn euc_eq_eps(a : Array[Double], b : Array[Double], eps : Double) -> Bool {
  if a.length() != b.length() {
    false
  } else {
    let mut ok = true
    for i = 0; i < a.length(); i = i + 1 {
      if (a[i] - b[i]).abs() >= eps {
        ok = false
      }
    }
    ok
  }
}

/// E-15 最大值分量
pub fn euc_max_component(a : Array[Double]) -> Double {
  let mut m = a[0]
  for x in a {
    if x > m {
      m = x
    }
  }
  m
}

/// E-16 最大絕對值分量（∞-範數）
pub fn euc_inf_norm(a : Array[Double]) -> Double {
  let mut m = a[0].abs()
  for x in a {
    if x.abs() > m {
      m = x.abs()
    }
  }
  m
}
