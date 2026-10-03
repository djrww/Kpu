///|
// =====================================================================
// ipspace.mbt — 有限維內積空間公理結構與衍生演算法
// 公理（自證於 ipspace_wbtest.mbt）：
//   對稱性 ⟨u,v⟩=⟨v,u⟩、雙線性、正定性 ⟨v,v⟩≥0 且 =0 ⟺ v=0
// 衍生：範數、距離、角度、正交性、投影、Gram–Schmidt、Bessel 等
// =====================================================================

/// 內積空間字典 = 向量空間 + 內積
pub struct IPSp[V] {
  sp : VecSp[V]
  inner : (V, V) -> Double
}

/// 從向量空間字典與內積構造
pub fn[V] ips_make(sp : VecSp[V], inner : (V, V) -> Double) -> IPSp[V] {
  { sp, inner }
}

/// I-1 範數 ‖v‖ = √⟨v,v⟩
pub fn[V] ip_norm(ip : IPSp[V], v : V) -> Double {
  let g = (ip.inner)(v, v)
  if g < 0.0 {
    0.0
  } else {
    g.sqrt()
  }
}

/// I-2 距離 d(u,v) = ‖u−v‖
pub fn[V] ip_dist(ip : IPSp[V], a : V, b : V) -> Double {
  let d = vsub(ip.sp, a, b)
  ip_norm(ip, d)
}

/// I-3 餘弦角 cos θ = ⟨u,v⟩/(‖u‖‖v‖)（零向量時回 1）
pub fn[V] ip_cos(ip : IPSp[V], a : V, b : V) -> Double {
  let na = ip_norm(ip, a)
  let nb = ip_norm(ip, b)
  if na == 0.0 || nb == 0.0 {
    1.0
  } else {
    let c = (ip.inner)(a, b) / (na * nb)
    if c > 1.0 {
      1.0
    } else if c < -1.0 {
      -1.0
    } else {
      c
    }
  }
}

/// I-4 角度 θ = arccos(cos θ) ∈ [0, π]
pub fn[V] ip_angle(ip : IPSp[V], a : V, b : V) -> Double {
  @math.acos(ip_cos(ip, a, b))
}

/// I-5 正交判定 |⟨u,v⟩| < eps
pub fn[V] ip_is_orthogonal(ip : IPSp[V], a : V, b : V, eps : Double) -> Bool {
  (ip.inner)(a, b).abs() < eps
}

/// I-6 單位化 v/‖v‖（零向量回零）
pub fn[V] ip_normalize(ip : IPSp[V], v : V) -> V {
  let n = ip_norm(ip, v)
  if n == 0.0 {
    ip.sp.zero
  } else {
    (ip.sp.scale)(1.0 / n, v)
  }
}

/// I-7 投影係數 ⟨u,v⟩/⟨v,v⟩
pub fn[V] ip_proj_coeff(ip : IPSp[V], u : V, v : V) -> Double {
  (ip.inner)(u, v) / (ip.inner)(v, v)
}

/// I-8 正交投影 proj_v(u) = ⟨u,v⟩/⟨v,v⟩ · v
pub fn[V] ip_proj(ip : IPSp[V], u : V, v : V) -> V {
  (ip.sp.scale)(ip_proj_coeff(ip, u, v), v)
}

/// I-9 正交分量（拒絕）rej_v(u) = u − proj_v(u)
pub fn[V] ip_rej(ip : IPSp[V], u : V, v : V) -> V {
  vsub(ip.sp, u, ip_proj(ip, u, v))
}

/// I-10 Gram–Schmidt 正交化（丟棄線性相關者，容差 eps）
pub fn[V] ip_gram_schmidt(
  ip : IPSp[V],
  vs : Array[V],
  eps : Double,
) -> Array[V] {
  let out : Array[V] = []
  for v in vs {
    let mut w = v
    for u in out {
      let c = (ip.inner)(w, u) / (ip.inner)(u, u)
      w = vsub(ip.sp, w, (ip.sp.scale)(c, u))
    }
    if ip_norm(ip, w) > eps {
      out.push(w)
    }
  }
  out
}

/// I-11 Gram–Schmidt 標準正交化
pub fn[V] ip_orthonormal(ip : IPSp[V], vs : Array[V], eps : Double) -> Array[V] {
  let ortho = ip_gram_schmidt(ip, vs, eps)
  let out : Array[V] = []
  for v in ortho {
    out.push(ip_normalize(ip, v))
  }
  out
}

/// I-12 標準正交基下的座標 ⟨v, eᵢ⟩
pub fn[V] ip_coords(ip : IPSp[V], v : V, onb : Array[V]) -> Array[Double] {
  let out : Array[Double] = []
  for e in onb {
    out.push((ip.inner)(v, e))
  }
  out
}

/// I-13 由座標重構向量 Σ cᵢeᵢ（Fourier 重構）
pub fn[V] ip_from_coords(
  ip : IPSp[V],
  cs : Array[Double],
  onb : Array[V],
) -> V {
  vlinc(ip.sp, cs, onb)
}

/// I-14 子空間正交投影（輸入為子空間之標準正交基）
pub fn[V] ip_proj_subspace(ip : IPSp[V], v : V, onb : Array[V]) -> V {
  let cs = ip_coords(ip, v, onb)
  ip_from_coords(ip, cs, onb)
}

/// I-15 到子空間之距離 ‖v − proj v‖
pub fn[V] ip_dist_subspace(ip : IPSp[V], v : V, onb : Array[V]) -> Double {
  let p = ip_proj_subspace(ip, v, onb)
  ip_dist(ip, v, p)
}

/// I-16 鏡射（Householder 型）：沿單位法向量 n 反射 v ↦ v − 2⟨v,n⟩n
pub fn[V] ip_reflect(ip : IPSp[V], v : V, n_unit : V) -> V {
  let c = 2.0 * (ip.inner)(v, n_unit)
  vsub(ip.sp, v, (ip.sp.scale)(c, n_unit))
}

/// I-17 Gram 矩陣 Gᵢⱼ = ⟨vᵢ, vⱼ⟩（傳回扁平方陣，列主序）
pub fn[V] ip_gram_matrix(ip : IPSp[V], vs : Array[V]) -> Array[Double] {
  let n = vs.length()
  let g = Array::make(n * n, 0.0)
  for i = 0; i < n; i = i + 1 {
    for j = 0; j < n; j = j + 1 {
      g[i * n + j] = (ip.inner)(vs[i], vs[j])
    }
  }
  g
}

/// I-18 平行多面體體積 = √det(Gram)（依賴 matrix.mbt 之 det）
pub fn[V] ip_volume(ip : IPSp[V], vs : Array[V]) -> Double {
  let g = ip_gram_matrix(ip, vs)
  let d = flat_det(g, vs.length())
  if d < 0.0 {
    0.0
  } else {
    d.sqrt()
  }
}

/// I-19 Bessel 殘差 ‖v‖² − Σ⟨v,eᵢ⟩²（對標準正交系必 ≥ 0）
pub fn[V] ip_bessel_residual(ip : IPSp[V], v : V, onb : Array[V]) -> Double {
  let mut s = (ip.inner)(v, v)
  for e in onb {
    let c = (ip.inner)(v, e)
    s = s - c * c
  }
  s
}

/// I-20 加權內積尺度變換：⟨,⟩' = c·⟨,⟩
pub fn[V] ips_scale(ip : IPSp[V], c : Double) -> IPSp[V] {
  {
    sp: ip.sp,
    inner: fn(a, b) { c * (ip.inner)(a, b) },
  }
}
