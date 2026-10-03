///|
// =====================================================================
// tensor.mbt — 張量積空間 R^m ⊗ R^n（以 m×n 係數陣列表示）
// 與線性映射空間、矩陣空間之自然同構由測試自證
// =====================================================================

/// 張量（R^m ⊗ R^n 之元素）
pub struct Tensor {
  m : Int
  n : Int
  data : Array[Double]
} derive(Debug, Eq)

/// T-1 零張量
pub fn tensor_zero(m : Int, n : Int) -> Tensor {
  { m, n, data: Array::make(m * n, 0.0) }
}

/// T-2 純張量 u ⊗ v
pub fn pure_tensor(u : Array[Double], v : Array[Double]) -> Tensor {
  let out = tensor_zero(u.length(), v.length())
  for i = 0; i < u.length(); i = i + 1 {
    for j = 0; j < v.length(); j = j + 1 {
      out.data[i * v.length() + j] = u[i] * v[j]
    }
  }
  out
}

/// T-3 加法
pub fn tensor_add(a : Tensor, b : Tensor) -> Tensor {
  let out = tensor_zero(a.m, a.n)
  for i = 0; i < a.m * a.n; i = i + 1 {
    out.data[i] = a.data[i] + b.data[i]
  }
  out
}

/// T-4 純量乘法
pub fn tensor_scale(c : Double, a : Tensor) -> Tensor {
  let out = tensor_zero(a.m, a.n)
  for i = 0; i < a.m * a.n; i = i + 1 {
    out.data[i] = c * a.data[i]
  }
  out
}

/// T-5 誘導內積 ⟨Σ aᵢⱼ eᵢ⊗fⱼ, Σ bᵢⱼ eᵢ⊗fⱼ⟩ = Σ aᵢⱼbᵢⱼ
pub fn tensor_inner(a : Tensor, b : Tensor) -> Double {
  let mut s = 0.0
  for i = 0; i < a.m * a.n; i = i + 1 {
    s = s + a.data[i] * b.data[i]
  }
  s
}

/// T-6 內積空間字典（R^m ⊗ R^n）
pub fn tensor_ip(m : Int, n : Int) -> IPSp[Tensor] {
  ips_make(
    {
      zero: tensor_zero(m, n),
      add: tensor_add,
      neg: fn(t) { tensor_scale(-1.0, t) },
      scale: tensor_scale,
    },
    tensor_inner,
  )
}

/// T-7 與矩陣空間之同構（平坦化）
pub fn tensor_to_mat(t : Tensor) -> Mat {
  { rows: t.m, cols: t.n, data: t.data.copy() }
}

/// T-8 矩陣 → 張量（同構之逆）
pub fn mat_to_tensor(m : Mat) -> Tensor {
  { m: m.rows, n: m.cols, data: m.data.copy() }
}

/// T-9 全收縮（m = n 時）：Σᵢ tᵢᵢ ∈ R
pub fn tensor_contract(t : Tensor) -> Double {
  let k = if t.m < t.n { t.m } else { t.n }
  let mut s = 0.0
  for i = 0; i < k; i = i + 1 {
    s = s + t.data[i * t.n + i]
  }
  s
}

/// T-10 單腿收縮：以 φ ∈ (R^m)* 收縮第一腿，得 R^n 之向量
pub fn tensor_contract_first(t : Tensor, phi : Array[Double]) -> Array[Double] {
  let out = Array::make(t.n, 0.0)
  for j = 0; j < t.n; j = j + 1 {
    let mut s = 0.0
    for i = 0; i < t.m; i = i + 1 {
      s = s + phi[i] * t.data[i * t.n + j]
    }
    out[j] = s
  }
  out
}

/// T-11 單腿收縮：以 ψ ∈ (R^n)* 收縮第二腿，得 R^m 之向量
pub fn tensor_contract_second(t : Tensor, psi : Array[Double]) -> Array[Double] {
  let out = Array::make(t.m, 0.0)
  for i = 0; i < t.m; i = i + 1 {
    let mut s = 0.0
    for j = 0; j < t.n; j = j + 1 {
      s = s + t.data[i * t.n + j] * psi[j]
    }
    out[i] = s
  }
  out
}

/// T-12 雙線性求值：(φ ⊗ ψ)(T) = Σ tᵢⱼ φᵢ ψⱼ
pub fn tensor_bilinear_eval(
  t : Tensor,
  phi : Array[Double],
  psi : Array[Double],
) -> Double {
  let mut s = 0.0
  for i = 0; i < t.m; i = i + 1 {
    for j = 0; j < t.n; j = j + 1 {
      s = s + t.data[i * t.n + j] * phi[i] * psi[j]
    }
  }
  s
}

/// T-13 線性映射之張量積 A ⊗ B（Kronecker）
pub fn tensor_of_maps(a : Mat, b : Mat) -> Tensor {
  let k = mat_kronecker(a, b)
  mat_to_tensor(k)
}

/// T-14 張量秩上界估计（平坦化矩陣之秩）
pub fn tensor_flatten_rank(t : Tensor, eps : Double) -> Int {
  mat_rank(tensor_to_mat(t), eps)
}
