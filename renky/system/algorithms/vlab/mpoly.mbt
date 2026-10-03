///|
// =====================================================================
// mpoly.mbt — 多變量多項式空間 R[x_0,…,x_{n-1}]
// 表示：項之陣列 (指數向量, 係數)，按字典序正規化
// =====================================================================

/// 多變量多項式
pub struct MPoly {
  nvars : Int
  terms : Array[(Array[Int], Double)]
}

/// 指數向量字典序比較
fn expo_cmp(a : Array[Int], b : Array[Int]) -> Int {
  let mut i = 0
  while i < a.length() && a[i] == b[i] {
    i = i + 1
  }
  if i == a.length() {
    0
  } else if a[i] < b[i] {
    -1
  } else {
    1
  }
}

/// MP-1 正規化：排序、合併同類項、去零
pub fn mpoly_normalize(
  nvars : Int,
  terms : Array[(Array[Int], Double)],
) -> MPoly {
  let t = terms.copy()
  // 插入排序（按指數向量字典序）
  for i = 1; i < t.length(); i = i + 1 {
    let key = t[i]
    let mut j = i - 1
    while j >= 0 && expo_cmp(t[j].0, key.0) > 0 {
      t[j + 1] = t[j]
      j = j - 1
    }
    t[j + 1] = key
  }
  let out : Array[(Array[Int], Double)] = []
  for term in t {
    if term.1 != 0.0 {
      if out.length() > 0 && expo_cmp(out[out.length() - 1].0, term.0) == 0 {
        let last = out[out.length() - 1]
        out[out.length() - 1] = (last.0, last.1 + term.1)
      } else {
        out.push(term)
      }
    }
  }
  // 合併後可能產生零
  let final : Array[(Array[Int], Double)] = []
  for term in out {
    if term.1 != 0.0 {
      final.push(term)
    }
  }
  { nvars, terms: final }
}

/// MP-2 零多項式
pub fn mpoly_zero(nvars : Int) -> MPoly {
  { nvars, terms: [] }
}

/// MP-3 單項式 c·x^e
pub fn mpoly_monomial(e : Array[Int], c : Double) -> MPoly {
  mpoly_normalize(e.length(), [(e.copy(), c)])
}

/// MP-4 加法
pub fn mpoly_add(a : MPoly, b : MPoly) -> MPoly {
  let terms : Array[(Array[Int], Double)] = []
  for t in a.terms {
    terms.push(t)
  }
  for t in b.terms {
    terms.push(t)
  }
  mpoly_normalize(a.nvars, terms)
}

/// MP-5 純量乘法
pub fn mpoly_scale(c : Double, p : MPoly) -> MPoly {
  let terms : Array[(Array[Int], Double)] = []
  for t in p.terms {
    terms.push((t.0, c * t.1))
  }
  mpoly_normalize(p.nvars, terms)
}

/// MP-6 乘法
pub fn mpoly_mul(a : MPoly, b : MPoly) -> MPoly {
  let terms : Array[(Array[Int], Double)] = []
  for ta in a.terms {
    for tb in b.terms {
      let e : Array[Int] = []
      for i = 0; i < a.nvars; i = i + 1 {
        e.push(ta.0[i] + tb.0[i])
      }
      terms.push((e, ta.1 * tb.1))
    }
  }
  mpoly_normalize(a.nvars, terms)
}

/// MP-7 單項式求值
pub fn monomial_eval(e : Array[Int], x : Array[Double]) -> Double {
  let mut p = 1.0
  for i = 0; i < e.length(); i = i + 1 {
    for _k = 0; _k < e[i]; _k = _k + 1 {
      p = p * x[i]
    }
  }
  p
}

/// MP-8 多項式求值
pub fn mpoly_eval(p : MPoly, x : Array[Double]) -> Double {
  let mut s = 0.0
  for t in p.terms {
    s = s + t.1 * monomial_eval(t.0, x)
  }
  s
}

/// MP-9 偏導 ∂/∂x_i
pub fn mpoly_partial(p : MPoly, i : Int) -> MPoly {
  let terms : Array[(Array[Int], Double)] = []
  for t in p.terms {
    if t.0[i] > 0 {
      let e = t.0.copy()
      e[i] = e[i] - 1
      terms.push((e, t.1 * t.0[i].to_double()))
    }
  }
  mpoly_normalize(p.nvars, terms)
}

/// MP-10 梯度（全部偏導）
pub fn mpoly_gradient(p : MPoly) -> Array[MPoly] {
  let out : Array[MPoly] = []
  for i = 0; i < p.nvars; i = i + 1 {
    out.push(mpoly_partial(p, i))
  }
  out
}

/// MP-11 總次數
pub fn mpoly_total_degree(p : MPoly) -> Int {
  let mut d = 0
  for t in p.terms {
    let mut s = 0
    for e in t.0 {
      s = s + e
    }
    if s > d {
      d = s
    }
  }
  d
}

/// MP-12 係數內積（單項式基下之 ℓ² 內積）
pub fn mpoly_inner(a : MPoly, b : MPoly) -> Double {
  let mut s = 0.0
  let mut i = 0
  let mut j = 0
  while i < a.terms.length() && j < b.terms.length() {
    let c = expo_cmp(a.terms[i].0, b.terms[j].0)
    if c == 0 {
      s = s + a.terms[i].1 * b.terms[j].1
      i = i + 1
      j = j + 1
    } else if c < 0 {
      i = i + 1
    } else {
      j = j + 1
    }
  }
  s
}

/// MP-13 單項式計數：n 變量次數 ≤ d 之單項式數 = C(n+d, d)
pub fn monomial_count(n : Int, d : Int) -> Int {
  // C(n+d, d) 以乘法公式計算
  let mut num = 1.0
  let mut den = 1.0
  for k = 1; k <= d; k = k + 1 {
    num = num * (n + k).to_double()
    den = den * k.to_double()
  }
  (num / den + 0.5).trunc().to_int()
}

/// MP-14 多項式空間之內積字典（限制於次數 ≤ d 之子空間）
pub fn mpoly_ip(nvars : Int) -> IPSp[MPoly] {
  ips_make(
    {
      zero: mpoly_zero(nvars),
      add: mpoly_add,
      neg: fn(p) { mpoly_scale(-1.0, p) },
      scale: mpoly_scale,
    },
    mpoly_inner,
  )
}

/// MP-15 齊次化：引入新變量 x_n 使總次數齊次
pub fn mpoly_homogenize(p : MPoly) -> MPoly {
  let d = mpoly_total_degree(p)
  let terms : Array[(Array[Int], Double)] = []
  for t in p.terms {
    let mut s = 0
    for e in t.0 {
      s = s + e
    }
    let e = t.0.copy()
    e.push(d - s)
    terms.push((e, t.1))
  }
  mpoly_normalize(p.nvars + 1, terms)
}
