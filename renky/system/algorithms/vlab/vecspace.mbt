///|
// =====================================================================
// vecspace.mbt — 向量空間公理結構與衍生演算法（實純量域）
// 公理（各實例自證於 vecspace_wbtest.mbt）：
//   (V,+) 交換群；純量乘法：1v=v, (ab)v=a(bv), a(u+v)=au+av, (a+b)v=av+bv
// 設計：以「字典記錄」攜帶運算，使空間可自由組合（直和、函數空間等）
// =====================================================================

/// 實向量空間字典
pub struct VecSp[V] {
  zero : V
  add : (V, V) -> V
  neg : (V) -> V
  scale : (Double, V) -> V
}

/// A-1 減法 u−v := u+(−v)
pub fn[V] vsub(sp : VecSp[V], a : V, b : V) -> V {
  (sp.add)(a, (sp.neg)(b))
}

/// A-2 線性組合 Σ cᵢvᵢ
pub fn[V] vlinc(sp : VecSp[V], cs : Array[Double], vs : Array[V]) -> V {
  let mut acc = sp.zero
  for i = 0; i < cs.length(); i = i + 1 {
    acc = (sp.add)(acc, (sp.scale)(cs[i], vs[i]))
  }
  acc
}

/// A-3 axpy: a·x + y（BLAS 基本算式）
pub fn[V] vaxpy(sp : VecSp[V], a : Double, x : V, y : V) -> V {
  (sp.add)((sp.scale)(a, x), y)
}

/// A-4 向量列求和
pub fn[V] vsum(sp : VecSp[V], vs : Array[V]) -> V {
  let mut acc = sp.zero
  for v in vs {
    acc = (sp.add)(acc, v)
  }
  acc
}

/// A-5 中點（仿射組合 ½a + ½b）
pub fn[V] vmidpoint(sp : VecSp[V], a : V, b : V) -> V {
  (sp.scale)(0.5, (sp.add)(a, b))
}

/// A-6 倍增迭代：2ⁿ·v（僅用加法實作純量乘法，驗證分配律用）
pub fn[V] vdouble_iter(sp : VecSp[V], v : V, n : Int) -> V {
  let mut acc = v
  for i = 0; i < n; i = i + 1 {
    let _ = i
    acc = (sp.add)(acc, acc)
  }
  acc
}

/// A-7 仿射組合 Σ cᵢvᵢ（Σcᵢ=1 時結果與原點無關之性質由測試驗證）
pub fn[V] vaffine(
  sp : VecSp[V],
  cs : Array[Double],
  vs : Array[V],
) -> V {
  vlinc(sp, cs, vs)
}

/// A-8 純量乘法之逆：(1/c)·v
pub fn[V] vunscale(sp : VecSp[V], c : Double, v : V) -> V {
  (sp.scale)(1.0 / c, v)
}
