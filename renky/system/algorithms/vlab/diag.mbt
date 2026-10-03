///|
// =====================================================================
// diag.mbt — 對角矩陣空間 Diag(n) ≅ R^n（逐點運算）
// =====================================================================

/// D-1 對角空間維數
pub fn diag_dim(n : Int) -> Int {
  n
}

/// D-2 對角加法（逐點）
pub fn diag_add(a : Array[Double], b : Array[Double]) -> Array[Double] {
  euc_add(a, b)
}

/// D-3 對角純量乘法
pub fn diag_scale(c : Double, a : Array[Double]) -> Array[Double] {
  euc_scale(c, a)
}

/// D-4 對角矩陣乘法 ≙ 逐點積
pub fn diag_mul(a : Array[Double], b : Array[Double]) -> Array[Double] {
  euc_hadamard(a, b)
}

/// D-5 對角內積 ≙ Σ aᵢbᵢ（Frobenius 內積的限制）
pub fn diag_inner(a : Array[Double], b : Array[Double]) -> Double {
  euc_inner(a, b)
}

/// D-6 作用於向量：diag(a)·x
pub fn diag_apply(a : Array[Double], x : Array[Double]) -> Array[Double] {
  euc_hadamard(a, x)
}

/// D-7 行列式 = 對角線之積
pub fn diag_det(a : Array[Double]) -> Double {
  let mut p = 1.0
  for x in a {
    p = p * x
  }
  p
}

/// D-8 逆（無零分量時）
pub fn diag_inv(a : Array[Double]) -> Option[Array[Double]] {
  let mut ok = true
  for x in a {
    if x == 0.0 {
      ok = false
    }
  }
  if !ok {
    Option::None
  } else {
    let out : Array[Double] = []
    for x in a {
      out.push(1.0 / x)
    }
    Some(out)
  }
}

/// D-9 轉為完整矩陣（與 Diag(n) → M_n 之嵌入）
pub fn diag_to_mat(a : Array[Double]) -> Mat {
  mat_diag(a)
}

/// D-10 Diag(n) 內積空間字典
pub fn diag_ip(n : Int) -> IPSp[Array[Double]] {
  ips_make(euc_sp(n), diag_inner)
}

/// D-11 譜半徑（對角矩陣之算子範數）
pub fn diag_spectral_radius(a : Array[Double]) -> Double {
  euc_inf_norm(a)
}

/// D-12 跡 = 分量和
pub fn diag_trace(a : Array[Double]) -> Double {
  euc_sum(a)
}
