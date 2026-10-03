///|
// =====================================================================
// field.mbt — 純量域實作：實數 Double 與複數域 C
// 域公理（自證於 field_wbtest.mbt）：
//   加法交換群、乘法交換群(去0)、分配律、共軛對合、|zw|=|z||w|
// =====================================================================

/// 複數：re + im·i
pub struct C {
  re : Double
  im : Double
} derive(Debug, Eq)

/// 構造
pub fn cnew(re : Double, im : Double) -> C {
  { re, im }
}

/// 0
pub fn czero() -> C {
  { re: 0.0, im: 0.0 }
}

/// 1
pub fn cone() -> C {
  { re: 1.0, im: 0.0 }
}

/// i
pub fn ci() -> C {
  { re: 0.0, im: 1.0 }
}

/// 實數嵌入
pub fn cfrom_real(x : Double) -> C {
  { re: x, im: 0.0 }
}

/// 加法
pub fn cadd(a : C, b : C) -> C {
  { re: a.re + b.re, im: a.im + b.im }
}

/// 減法
pub fn csub(a : C, b : C) -> C {
  { re: a.re - b.re, im: a.im - b.im }
}

/// 負元
pub fn cneg(a : C) -> C {
  { re: -a.re, im: -a.im }
}

/// 實純量乘法
pub fn cscale(s : Double, a : C) -> C {
  { re: s * a.re, im: s * a.im }
}

/// 乘法 (a+bi)(c+di) = (ac−bd) + (ad+bc)i
pub fn cmul(a : C, b : C) -> C {
  { re: a.re * b.re - a.im * b.im, im: a.re * b.im + a.im * b.re }
}

/// 平方
pub fn csqr(a : C) -> C {
  cmul(a, a)
}

/// 共軛
pub fn cconj(a : C) -> C {
  { re: a.re, im: -a.im }
}

/// |z|²
pub fn cnorm2(a : C) -> Double {
  a.re * a.re + a.im * a.im
}

/// |z|
pub fn cabs(a : C) -> Double {
  cnorm2(a).sqrt()
}

/// 乘法逆元 z⁻¹ = z̄/|z|²
pub fn cinv(a : C) -> C {
  let d = cnorm2(a)
  { re: a.re / d, im: -a.im / d }
}

/// 除法
pub fn cdiv(a : C, b : C) -> C {
  cmul(a, cinv(b))
}

/// 極式 r·e^{iθ}
pub fn cpolar(r : Double, theta : Double) -> C {
  { re: r * @math.cos(theta), im: r * @math.sin(theta) }
}

/// 輻角主值
pub fn carg(a : C) -> Double {
  @math.atan2(a.im, a.re)
}

/// 指數函數 e^{a}
pub fn cexp(a : C) -> C {
  cscale(@math.exp(a.re), cpolar(1.0, a.im))
}

/// 主平方根
pub fn csqrt(a : C) -> C {
  let r = cabs(a)
  if r == 0.0 {
    czero()
  } else {
    let re = ((r + a.re) / 2.0).sqrt()
    let im = ((r - a.re) / 2.0).sqrt()
    if a.im >= 0.0 {
      { re, im }
    } else {
      { re, im: -im }
    }
  }
}

/// 近似相等
pub fn ceq_eps(a : C, b : C, eps : Double) -> Bool {
  (a.re - b.re).abs() < eps && (a.im - b.im).abs() < eps
}

/// 顯示
pub fn cstr(a : C) -> String {
  "\{a.re}+\{a.im}*i"
}
