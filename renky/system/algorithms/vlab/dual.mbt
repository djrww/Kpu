///|
// =====================================================================
// dual.mbt — 對偶空間 V*（有限維實向量空間之線性泛函空間）
// 以行向量表示；Riesz 表示、對偶映射、雙重對偶自證
// =====================================================================

/// 對偶向量（線性泛函，行向量表示）
pub struct DualVec {
  coef : Array[Double]
} derive(Debug, Eq)

/// D-1 零泛函
pub fn dual_zero(n : Int) -> DualVec {
  { coef: Array::make(n, 0.0) }
}

/// D-2 泛函作用 φ(v)
pub fn dual_action(phi : DualVec, v : Array[Double]) -> Double {
  euc_inner(phi.coef, v)
}

/// D-3 對偶加法
pub fn dual_add(a : DualVec, b : DualVec) -> DualVec {
  { coef: euc_add(a.coef, b.coef) }
}

/// D-4 對偶純量乘法
pub fn dual_scale(c : Double, phi : DualVec) -> DualVec {
  { coef: euc_scale(c, phi.coef) }
}

/// D-5 對偶內積（標準歐氏度量下 V ≅ V* 之誘導內積）
pub fn dual_inner(a : DualVec, b : DualVec) -> Double {
  euc_inner(a.coef, b.coef)
}

/// D-6 對偶空間字典（n 維）
pub fn dual_ip(n : Int) -> IPSp[DualVec] {
  ips_make(
    {
      zero: dual_zero(n),
      add: dual_add,
      neg: fn(p) { { coef: euc_neg(p.coef) } },
      scale: dual_scale,
    },
    dual_inner,
  )
}

/// D-7 對偶基（對標準基即 δᵢ = eᵢ*）
pub fn dual_basis(n : Int) -> Array[DualVec] {
  let out : Array[DualVec] = []
  for i = 0; i < n; i = i + 1 {
    out.push({ coef: evec(n, i) })
  }
  out
}

/// D-8 Riesz 表示定理（歐氏情形）：泛函 φ ↦ 表示向量 r，
/// 使得 φ(v) = ⟨r, v⟩。標準度量下 r = φ 之係數。
pub fn riesz_represent(phi : DualVec) -> Array[Double] {
  phi.coef.copy()
}

/// D-9 Riesz 逆：向量 ↦ 泛函 ⟨v, ·⟩
pub fn riesz_inverse(v : Array[Double]) -> DualVec {
  { coef: v.copy() }
}

/// D-10 對偶映射（拉回）：線性映射 A : V → W 誘導
/// A* : W* → V*，矩陣表示為轉置。輸入 A（m×n 行主序，
/// 作用 V=R^n → W=R^m），輸出 A* 之矩陣（n×m）
pub fn dual_map(a : Mat) -> Mat {
  mat_transpose(a)
}

/// D-11 對偶配對矩陣：Pᵢⱼ = φᵢ(eⱼ)（對偶基與基之配對）
pub fn dual_pairing_matrix(phis : Array[DualVec], n : Int) -> Mat {
  let out = mat_new(phis.length(), n)
  for i = 0; i < phis.length(); i = i + 1 {
    for j = 0; j < n; j = j + 1 {
      mat_set(out, i, j, dual_action(phis[i], evec(n, j)))
    }
  }
  out
}

/// D-12 雙重對偶嵌入 J : V → V**（有限維時為同構，自證 J = id）
pub fn double_dual_action(v : Array[Double], phis : Array[DualVec]) -> Array[Double] {
  let out : Array[Double] = []
  for phi in phis {
    out.push(dual_action(phi, v))
  }
  out
}

/// D-13 核之判定：φ ∈ (span S)° ⟺ ∀s∈S: φ(s)=0（回傳是否屬零化子）
pub fn dual_in_annihilator(
  phi : DualVec,
  span_gen : Array[Array[Double]],
  eps : Double,
) -> Bool {
  let mut ok = true
  for s in span_gen {
    if dual_action(phi, s).abs() > eps {
      ok = false
    }
  }
  ok
}
