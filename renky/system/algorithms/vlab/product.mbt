///|
// =====================================================================
// product.mbt — 空間組合子：直和（積空間）、對角嵌入、加權組合
// 「空間 × 空間」之組合，可迭代套用（如 Fun(X, Sym(n)) ≅ 直和）
// =====================================================================

/// PR-1 向量空間直和 V ⊕ W
pub fn[V, W] vec_prod(a : VecSp[V], b : VecSp[W]) -> VecSp[(V, W)] {
  {
    zero: (a.zero, b.zero),
    add: fn(p, q) { ((a.add)(p.0, q.0), (b.add)(p.1, q.1)) },
    neg: fn(p) { ((a.neg)(p.0), (b.neg)(p.1)) },
    scale: fn(c, p) { ((a.scale)(c, p.0), (b.scale)(c, p.1)) },
  }
}

/// PR-2 內積空間直和：⟨(v₁,w₁),(v₂,w₂)⟩ = ⟨v₁,v₂⟩ + ⟨w₁,w₂⟩
pub fn[V, W] ips_prod(a : IPSp[V], b : IPSp[W]) -> IPSp[(V, W)] {
  ips_make(
    vec_prod(a.sp, b.sp),
    fn(p, q) { (a.inner)(p.0, q.0) + (b.inner)(p.1, q.1) },
  )
}

/// PR-3 加權直和：⟨,⟩ = α⟨,⟩_V + β⟨,⟩_W
pub fn[V, W] ips_prod_weighted(
  a : IPSp[V],
  b : IPSp[W],
  alpha : Double,
  beta : Double,
) -> IPSp[(V, W)] {
  ips_make(
    vec_prod(a.sp, b.sp),
    fn(p, q) { alpha * (a.inner)(p.0, q.0) + beta * (b.inner)(p.1, q.1) },
  )
}

/// PR-4 對角嵌入 Δ : V → V ⊕ V，v ↦ (v, v)
pub fn[V] diag_embed(v : V) -> (V, V) {
  (v, v)
}

/// PR-5 投影 π₁ : V ⊕ W → V
pub fn[V, W] prod_proj1(p : (V, W)) -> V {
  p.0
}

/// PR-6 投影 π₂ : V ⊕ W → W
pub fn[V, W] prod_proj2(p : (V, W)) -> W {
  p.1
}

/// PR-7 單射 ι₁ : V → V ⊕ W（w 取零）
pub fn[V, W] inj1(a : IPSp[V], b : IPSp[W], v : V) -> (V, W) {
  (v, b.sp.zero)
}

/// PR-8 單射 ι₂ : W → V ⊕ W
pub fn[V, W] inj2(a : IPSp[V], b : IPSp[W], w : W) -> (V, W) {
  (a.sp.zero, w)
}

/// PR-9 直和維數相加
pub fn prod_dim_sum(dv : Int, dw : Int) -> Int {
  dv + dw
}

/// PR-10 線性映射直和 A ⊕ B 之矩陣（分塊對角）
pub fn map_direct_sum(a : Mat, b : Mat) -> Mat {
  let out = mat_new(a.rows + b.rows, a.cols + b.cols)
  for i = 0; i < a.rows; i = i + 1 {
    for j = 0; j < a.cols; j = j + 1 {
      mat_set(out, i, j, mat_get(a, i, j))
    }
  }
  for i = 0; i < b.rows; i = i + 1 {
    for j = 0; j < b.cols; j = j + 1 {
      mat_set(out, a.rows + i, a.cols + j, mat_get(b, i, j))
    }
  }
  out
}
