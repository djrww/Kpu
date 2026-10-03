///|
// =====================================================================
// exterior.mbt — 外代數 Λ^k(R^n)
// 分量按遞增 k-組合索引；外積、Hodge 星、內縮、體積
// =====================================================================

/// k-向量（Λ^k(R^n) 之元素）
pub struct ExtVec {
  n : Int
  k : Int
  data : Array[Double] // 長度 C(n, k)，按組合字典序
} derive(Debug, Eq)

/// X-1 組合數 C(n, k)
pub fn comb_count(n : Int, k : Int) -> Int {
  if k < 0 || k > n {
    0
  } else {
    let mut num = 1.0
    let mut den = 1.0
    for i = 1; i <= k; i = i + 1 {
      num = num * (n - k + i).to_double()
      den = den * i.to_double()
    }
    (num / den + 0.5).trunc().to_int()
  }
}

/// 列舉全部遞增 k-組合（字典序）
pub fn comb_list(n : Int, k : Int) -> Array[Array[Int]] {
  let out : Array[Array[Int]] = []
  if k == 0 {
    out.push([])
    return out
  }
  let cur = Array::make(k, 0)
  fn rec(pos : Int, start : Int) -> Unit {
    if pos == k {
      out.push(cur.copy())
      return
    }
    for v = start; v <= n - (k - pos); v = v + 1 {
      cur[pos] = v
      rec(pos + 1, v + 1)
    }
  }
  rec(0, 0)
  out
}

/// X-2 組合 → 分量索引
pub fn comb_index(combs : Array[Array[Int]], c : Array[Int]) -> Int {
  let mut idx = -1
  for i = 0; i < combs.length(); i = i + 1 {
    let mut eq = true
    for j = 0; j < c.length(); j = j + 1 {
      if combs[i][j] != c[j] {
        eq = false
      }
    }
    if eq {
      idx = i
    }
  }
  idx
}

/// X-3 逆序數（排列奇偶性用）
pub fn inversion_count(p : Array[Int]) -> Int {
  let mut inv = 0
  for i = 0; i < p.length(); i = i + 1 {
    for j = i + 1; j < p.length(); j = j + 1 {
      if p[i] > p[j] {
        inv = inv + 1
      }
    }
  }
  inv
}

/// X-4 Λ^k 之零元素
pub fn ext_zero(n : Int, k : Int) -> ExtVec {
  { n, k, data: Array::make(comb_count(n, k), 0.0) }
}

/// X-5 向量 v ∈ R^n 嵌入 Λ^1
pub fn ext_from_vec(n : Int, v : Array[Double]) -> ExtVec {
  { n, k: 1, data: v.copy() }
}

/// X-6 外積 u ∧ v（兩個 1-向量 → 2-向量）
pub fn ext_wedge_vecs(n : Int, u : Array[Double], v : Array[Double]) -> ExtVec {
  let combs = comb_list(n, 2)
  let out = ext_zero(n, 2)
  for idx = 0; idx < combs.length(); idx = idx + 1 {
    let i = combs[idx][0]
    let j = combs[idx][1]
    out.data[idx] = u[i] * v[j] - u[j] * v[i]
  }
  out
}

/// X-7 一般外積：α ∈ Λ^p ∧ β ∈ Λ^q → Λ^{p+q}
pub fn ext_wedge(a : ExtVec, b : ExtVec) -> ExtVec {
  let n = a.n
  let p = a.k
  let q = b.k
  let ca = comb_list(n, p)
  let cb = comb_list(n, q)
  let cc = comb_list(n, p + q)
  let out = ext_zero(n, p + q)
  for ia = 0; ia < ca.length(); ia = ia + 1 {
    for ib = 0; ib < cb.length(); ib = ib + 1 {
      let coef = a.data[ia] * b.data[ib]
      if coef == 0.0 {
        continue
      }
      // 合併指數集；若相交則該項為 0
      let merged : Array[Int] = []
      let used = Array::make(n, false)
      let mut ok = true
      for x in ca[ia] {
        if used[x] {
          ok = false
        }
        used[x] = true
      }
      for x in cb[ib] {
        if used[x] {
          ok = false
        }
        used[x] = true
      }
      if !ok {
        continue
      }
      // 合併排列之奇偶性 = 交織符號
      let concat : Array[Int] = []
      for x in ca[ia] {
        concat.push(x)
      }
      for x in cb[ib] {
        concat.push(x)
      }
      for x = 0; x < n; x = x + 1 {
        if used[x] {
          merged.push(x)
        }
      }
      let sign = if inversion_count(concat) % 2 == 0 { 1.0 } else { -1.0 }
      let idx = comb_index(cc, merged)
      out.data[idx] = out.data[idx] + sign * coef
    }
  }
  out
}

/// X-8 Λ^k 之誘導內積（單形分量逐點積；對一般元素用 Gram 行列式性質）
pub fn ext_inner(a : ExtVec, b : ExtVec) -> Double {
  let mut s = 0.0
  for i = 0; i < a.data.length(); i = i + 1 {
    s = s + a.data[i] * b.data[i]
  }
  s
}

/// X-9 k-向量範數（= 對應平行體之體積，當為純量時）
pub fn ext_norm(a : ExtVec) -> Double {
  ext_inner(a, a).sqrt()
}

/// X-10 Hodge 星 ⋆ : Λ^k → Λ^{n−k}（標準定向）
pub fn ext_hodge_star(a : ExtVec) -> ExtVec {
  let n = a.n
  let k = a.k
  let ck = comb_list(n, k)
  let cnk = comb_list(n, n - k)
  let out = ext_zero(n, n - k)
  for idx = 0; idx < ck.length(); idx = idx + 1 {
    if a.data[idx] == 0.0 {
      continue
    }
    let used = Array::make(n, false)
    for x in ck[idx] {
      used[x] = true
    }
    let comp : Array[Int] = []
    for x = 0; x < n; x = x + 1 {
      if !used[x] {
        comp.push(x)
      }
    }
    // 符號 = (ck[idx] ++ comp) 之奇偶性
    let perm : Array[Int] = []
    for x in ck[idx] {
      perm.push(x)
    }
    for x in comp {
      perm.push(x)
    }
    let sign = if inversion_count(perm) % 2 == 0 { 1.0 } else { -1.0 }
    let j = comb_index(cnk, comp)
    out.data[j] = sign * a.data[idx]
  }
  out
}

/// X-11 內縮 i_v : Λ^k → Λ^{k−1}（與向量之內積）
pub fn ext_interior(v : Array[Double], a : ExtVec) -> ExtVec {
  let n = a.n
  let k = a.k
  if k == 0 {
    return ext_zero(n, 0)
  }
  let ck = comb_list(n, k)
  let ck1 = comb_list(n, k - 1)
  let out = ext_zero(n, k - 1)
  for idx = 0; idx < ck.length(); idx = idx + 1 {
    let coef = a.data[idx]
    if coef == 0.0 {
      continue
    }
    for pos = 0; pos < k; pos = pos + 1 {
      let removed : Array[Int] = []
      for j = 0; j < k; j = j + 1 {
        if j != pos {
          removed.push(ck[idx][j])
        }
      }
      let sign = if pos % 2 == 0 { 1.0 } else { -1.0 }
      let j = comb_index(ck1, removed)
      out.data[j] = out.data[j] + sign * v[ck[idx][pos]] * coef
    }
  }
  out
}

/// X-12 頂層外積判定行列式：v_1 ∧ … ∧ v_n = det(V)·e₁∧…∧e_n
pub fn ext_det(vectors : Array[Array[Double]]) -> Double {
  let n = vectors.length()
  if n == 0 {
    return 1.0
  }
  let mut acc = ext_from_vec(n, vectors[0])
  for i = 1; i < n; i = i + 1 {
    acc = ext_wedge(acc, ext_from_vec(n, vectors[i]))
  }
  acc.data[0]
}

/// X-13 R³ 叉積由外代數導出：u × v = ⋆(u ∧ v)
pub fn ext_cross3(u : Array[Double], v : Array[Double]) -> Array[Double] {
  let w = ext_wedge_vecs(3, u, v)
  let s = ext_hodge_star(w)
  s.data.copy()
}

/// X-14 Λ^k 之內積空間字典
pub fn ext_ip(n : Int, k : Int) -> IPSp[ExtVec] {
  ips_make(
    {
      zero: ext_zero(n, k),
      add: fn(a, b) {
        let out = ext_zero(n, k)
        for i = 0; i < a.data.length(); i = i + 1 {
          out.data[i] = a.data[i] + b.data[i]
        }
        out
      },
      neg: fn(a) {
        let out = ext_zero(n, k)
        for i = 0; i < a.data.length(); i = i + 1 {
          out.data[i] = -a.data[i]
        }
        out
      },
      scale: fn(c, a) {
        let out = ext_zero(n, k)
        for i = 0; i < a.data.length(); i = i + 1 {
          out.data[i] = c * a.data[i]
        }
        out
      },
    },
    ext_inner,
  )
}

/// X-15 兩向量張成平行四邊形之面積 = ‖u ∧ v‖
pub fn ext_parallelogram_area(u : Array[Double], v : Array[Double]) -> Double {
  ext_norm(ext_wedge_vecs(u.length(), u, v))
}
