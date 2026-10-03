///|
// =====================================================================
// symm.mbt — 對稱矩陣空間 Sym(n) ⊂ M_n(R)，Frobenius 內積
// 維數 n(n+1)/2；正交投影 = 對稱化（定理自證於 symm_wbtest.mbt）
// =====================================================================

/// S-1 對稱空間維數 n(n+1)/2
pub fn sym_dim(n : Int) -> Int {
  n * (n + 1) / 2
}

/// S-2 上三角座標 (i≤j) → 線性索引
pub fn sym_index(n : Int, i : Int, j : Int) -> Int {
  let a = if i <= j { i } else { j }
  let b = if i <= j { j } else { i }
  // 前 a 行共有 a*n - a*(a-1)/2 個元素
  a * n - a * (a - 1) / 2 + (b - a)
}

/// S-3 向量打包為對稱矩陣（維數 n(n+1)/2 之座標）
pub fn sym_pack(n : Int, v : Array[Double]) -> Mat {
  let out = mat_new(n, n)
  for i = 0; i < n; i = i + 1 {
    for j = i; j < n; j = j + 1 {
      let x = v[sym_index(n, i, j)]
      mat_set(out, i, j, x)
      mat_set(out, j, i, x)
    }
  }
  out
}

/// S-4 對稱矩陣展開為向量（pack 之逆）
pub fn sym_unpack(m : Mat) -> Array[Double] {
  let n = m.rows
  let out : Array[Double] = []
  for i = 0; i < n; i = i + 1 {
    for j = i; j < n; j = j + 1 {
      out.push(mat_get(m, i, j))
    }
  }
  out
}

/// S-5 Sym(n) 加法
pub fn sym_add(a : Mat, b : Mat) -> Mat {
  mat_add(a, b)
}

/// S-6 Sym(n) 純量乘法
pub fn sym_scale(c : Double, a : Mat) -> Mat {
  mat_scale(c, a)
}

/// S-7 Sym(n) 內積（Frobenius 限制）
pub fn sym_inner(a : Mat, b : Mat) -> Double {
  mat_frob_inner(a, b)
}

/// S-8 Sym(n) 向量空間字典
pub fn sym_sp(n : Int) -> VecSp[Mat] {
  {
    zero: mat_new(n, n),
    add: sym_add,
    neg: fn(m) { mat_scale(-1.0, m) },
    scale: sym_scale,
  }
}

/// S-9 Sym(n) 內積空間字典
pub fn sym_ip(n : Int) -> IPSp[Mat] {
  ips_make(sym_sp(n), sym_inner)
}

/// S-10 Sym(n) 標準正交基：{E_ii} ∪ {(E_ij+E_ji)/√2}
pub fn sym_basis(n : Int) -> Array[Mat] {
  let out : Array[Mat] = []
  for i = 0; i < n; i = i + 1 {
    let m = mat_new(n, n)
    mat_set(m, i, i, 1.0)
    out.push(m)
  }
  for i = 0; i < n; i = i + 1 {
    for j = i + 1; j < n; j = j + 1 {
      let m = mat_new(n, n)
      let c = 1.0 / (2.0).sqrt()
      mat_set(m, i, j, c)
      mat_set(m, j, i, c)
      out.push(m)
    }
  }
  out
}

/// S-11 譜定理驗證用：對稱矩陣特徵分解重構誤差
pub fn sym_spectral_error(a : Mat) -> Double {
  let n = a.rows
  let (vals, q) = jacobi_eigen(a, 60)
  let recon = spectral_reconstruct(vals, q)
  let _ = n
  mat_frob_norm(mat_sub(a, recon))
}

/// S-12 二次型 q(x) = xᵀAx
pub fn sym_quadratic(a : Mat, x : Array[Double]) -> Double {
  euc_inner(x, mat_vec(a, x))
}
