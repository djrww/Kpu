///|
// =====================================================================
// hermitian.mbt — 複內積空間 C^n 與 Hermitian 矩陣空間
// 內積約定：⟨u,v⟩ = Σ conj(uᵢ)vᵢ（第一變量共軛線性）
// =====================================================================

/// H-1 C^n 零向量
pub fn cvec_zero(n : Int) -> Array[C] {
  Array::make(n, czero())
}

/// H-2 C^n 加法
pub fn cvec_add(a : Array[C], b : Array[C]) -> Array[C] {
  let out : Array[C] = []
  for i = 0; i < a.length(); i = i + 1 {
    out.push(cadd(a[i], b[i]))
  }
  out
}

/// H-3 C^n 負元
pub fn cvec_neg(a : Array[C]) -> Array[C] {
  let out : Array[C] = []
  for x in a {
    out.push(cneg(x))
  }
  out
}

/// H-4 C^n 複純量乘法
pub fn cvec_cscale(c : C, a : Array[C]) -> Array[C] {
  let out : Array[C] = []
  for x in a {
    out.push(cmul(c, x))
  }
  out
}

/// H-5 C^n 實純量乘法（視為實向量空間）
pub fn cvec_scale(s : Double, a : Array[C]) -> Array[C] {
  let out : Array[C] = []
  for x in a {
    out.push(cscale(s, x))
  }
  out
}

/// H-6 Hermitian 內積 ⟨u,v⟩ = Σ conj(uᵢ)vᵢ
pub fn cvec_inner(a : Array[C], b : Array[C]) -> C {
  let mut s = czero()
  for i = 0; i < a.length(); i = i + 1 {
    s = cadd(s, cmul(cconj(a[i]), b[i]))
  }
  s
}

/// H-7 實值二次範數 ⟨v,v⟩ ∈ R
pub fn cvec_norm2(v : Array[C]) -> Double {
  let g = cvec_inner(v, v)
  g.re
}

/// H-8 範數
pub fn cvec_norm(v : Array[C]) -> Double {
  cvec_norm2(v).sqrt()
}

/// H-9 C^n 之內積空間字典（視為實內積空間：取內積實部）
pub fn cvec_real_ip(n : Int) -> IPSp[Array[C]] {
  ips_make(
    {
      zero: cvec_zero(n),
      add: cvec_add,
      neg: cvec_neg,
      scale: cvec_scale,
    },
    fn(a, b) { cvec_inner(a, b).re },
  )
}

/// H-10 正交性（複）
pub fn cvec_is_orthogonal(a : Array[C], b : Array[C], eps : Double) -> Bool {
  cabs(cvec_inner(a, b)) < eps
}

/// H-11 複 Gram–Schmidt 標準正交化
pub fn cvec_orthonormal(vs : Array[Array[C]], eps : Double) -> Array[Array[C]] {
  let out : Array[Array[C]] = []
  for v in vs {
    let mut w = v
    for u in out {
      // c = ⟨u,w⟩/⟨u,u⟩; w -= c·u
      let c = cdiv(cvec_inner(u, w), cfrom_real(cvec_norm2(u)))
      w = cvec_add(w, cvec_cscale(cneg(c), u))
    }
    let nrm = cvec_norm(w)
    if nrm > eps {
      out.push(cvec_scale(1.0 / nrm, w))
    }
  }
  out
}

/// 複矩陣
pub struct CMat {
  rows : Int
  cols : Int
  data : Array[C]
}

/// H-12 複零矩陣
pub fn cmat_new(rows : Int, cols : Int) -> CMat {
  { rows, cols, data: Array::make(rows * cols, czero()) }
}

/// H-13 取元素
pub fn cmat_get(m : CMat, i : Int, j : Int) -> C {
  m.data[i * m.cols + j]
}

/// H-14 設元素
pub fn cmat_set(m : CMat, i : Int, j : Int, v : C) -> Unit {
  m.data[i * m.cols + j] = v
}

/// H-15 共軛轉置（dagger）
pub fn cmat_dagger(m : CMat) -> CMat {
  let out = cmat_new(m.cols, m.rows)
  for i = 0; i < m.rows; i = i + 1 {
    for j = 0; j < m.cols; j = j + 1 {
      cmat_set(out, j, i, cconj(cmat_get(m, i, j)))
    }
  }
  out
}

/// H-16 Hermitian 判定 A = A†
pub fn cmat_is_hermitian(m : CMat, eps : Double) -> Bool {
  let d = cmat_dagger(m)
  let mut ok = true
  for i = 0; i < m.rows * m.cols; i = i + 1 {
    if !ceq_eps(m.data[i], d.data[i], eps) {
      ok = false
    }
  }
  ok
}

/// H-17 複矩陣-向量乘
pub fn cmat_vec(m : CMat, x : Array[C]) -> Array[C] {
  let out : Array[C] = []
  for i = 0; i < m.rows; i = i + 1 {
    let mut s = czero()
    for j = 0; j < m.cols; j = j + 1 {
      s = cadd(s, cmul(cmat_get(m, i, j), x[j]))
    }
    out.push(s)
  }
  out
}

/// H-18 複矩陣乘法
pub fn cmat_mul(a : CMat, b : CMat) -> CMat {
  let out = cmat_new(a.rows, b.cols)
  for i = 0; i < a.rows; i = i + 1 {
    for k = 0; k < b.cols; k = k + 1 {
      let mut s = czero()
      for j = 0; j < a.cols; j = j + 1 {
        s = cadd(s, cmul(cmat_get(a, i, j), cmat_get(b, j, k)))
      }
      cmat_set(out, i, k, s)
    }
  }
  out
}

/// H-19 複 Rayleigh 商（Hermitian 時為實數）
pub fn cmat_rayleigh(m : CMat, x : Array[C]) -> Double {
  let ax = cmat_vec(m, x)
  cvec_inner(x, ax).re / cvec_norm2(x)
}

/// H-20 Hermitian 矩陣特徵值（複 Jacobi 旋轉，非對角元消去）
/// 回傳實特徵值陣列與酉矩陣 Q（行為特徵向量）
pub fn cmat_jacobi_hermitian(
  m : CMat,
  max_sweeps : Int,
) -> (Array[Double], CMat) {
  let n = m.rows
  let s = m.data.copy()
  let q = cmat_new(n, n)
  for i = 0; i < n; i = i + 1 {
    cmat_set(q, i, i, cone())
  }
  for _sw = 0; _sw < max_sweeps; _sw = _sw + 1 {
    let mut off = 0.0
    for i = 0; i < n; i = i + 1 {
      for j = i + 1; j < n; j = j + 1 {
        off = off + cnorm2(s[i * n + j])
      }
    }
    if off.sqrt() < 1.0e-13 {
      break
    }
    for p = 0; p < n - 1; p = p + 1 {
      for r = p + 1; r < n; r = r + 1 {
        let apq = s[p * n + r]
        if cabs(apq) < 1.0e-16 {
          continue
        }
        // G = P·R：先相位對角化 P[r,r]=e^{−iθ}，再實旋轉 R
        let rho = cabs(apq)
        let phase = cdiv(apq, cfrom_real(rho)) // e^{iθ}
        let phc = cconj(phase) // e^{−iθ}
        let app = s[p * n + p].re
        let aqq = s[r * n + r].re
        let tau = (aqq - app) / (2.0 * rho)
        // 歸零條件之最小根：t² − 2τt − 1 = 0
        let t =
          if tau >= 0.0 {
            -1.0 / (tau + (tau * tau + 1.0).sqrt())
          } else {
            1.0 / (-tau + (tau * tau + 1.0).sqrt())
          }
        let c = 1.0 / (t * t + 1.0).sqrt()
        let ss = t * c
        // S ← G† S G：先更新行（右乘 G）
        for i = 0; i < n; i = i + 1 {
          let sip = s[i * n + p]
          let siq = s[i * n + r]
          s[i * n + p] = cadd(cscale(c, sip), cmul(cscale(ss, phc), siq))
          s[i * n + r] = cadd(cscale(-ss, sip), cmul(cscale(c, phc), siq))
        }
        // 再更新列（左乘 G†）
        for j = 0; j < n; j = j + 1 {
          let spj = s[p * n + j]
          let sqj = s[r * n + j]
          s[p * n + j] = cadd(cscale(c, spj), cmul(cscale(ss, phase), sqj))
          s[r * n + j] = cadd(cscale(-ss, spj), cmul(cscale(c, phase), sqj))
        }
        // 累積 Q ← Q G（行變換）
        for i = 0; i < n; i = i + 1 {
          let qip = q.data[i * n + p]
          let qiq = q.data[i * n + r]
          q.data[i * n + p] = cadd(cscale(c, qip), cmul(cscale(ss, phc), qiq))
          q.data[i * n + r] = cadd(cscale(-ss, qip), cmul(cscale(c, phc), qiq))
        }
      }
    }
  }
  let vals : Array[Double] = []
  for i = 0; i < n; i = i + 1 {
    vals.push(s[i * n + i].re)
  }
  (vals, { rows: n, cols: n, data: q.data })
}

/// H-21 酉性誤差 ‖Q†Q − I‖_F
pub fn cmat_unitary_error(q : CMat) -> Double {
  let prod = cmat_mul(cmat_dagger(q), q)
  let mut err = 0.0
  for i = 0; i < prod.rows; i = i + 1 {
    for j = 0; j < prod.cols; j = j + 1 {
      let want = if i == j { 1.0 } else { 0.0 }
      let d = csub(cmat_get(prod, i, j), cfrom_real(want))
      err = err + cnorm2(d)
    }
  }
  err.sqrt()
}
