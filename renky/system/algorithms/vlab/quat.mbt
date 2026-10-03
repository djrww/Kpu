///|
// =====================================================================
// quat.mbt — 四元數代數 H 與 H^n 作為實內積空間
// ⟨q,p⟩ = Re(q̄p)；‖q‖² = q̄q；旋轉由單位四元數共軛作用給出
// =====================================================================

/// 四元數 q = w + xi + yj + zk
pub struct Quat {
  w : Double
  x : Double
  y : Double
  z : Double
} derive(Debug, Eq)

/// Q-1 構造
pub fn quat_new(w : Double, x : Double, y : Double, z : Double) -> Quat {
  { w, x, y, z }
}

/// Q-2 零
pub fn quat_zero() -> Quat {
  { w: 0.0, x: 0.0, y: 0.0, z: 0.0 }
}

/// Q-3 純量四元數
pub fn quat_real(w : Double) -> Quat {
  { w, x: 0.0, y: 0.0, z: 0.0 }
}

/// Q-4 純虛四元數（對應 R³ 向量）
pub fn quat_imag(x : Double, y : Double, z : Double) -> Quat {
  { w: 0.0, x, y, z }
}

/// Q-5 加法
pub fn quat_add(a : Quat, b : Quat) -> Quat {
  { w: a.w + b.w, x: a.x + b.x, y: a.y + b.y, z: a.z + b.z }
}

/// Q-6 負元
pub fn quat_neg(a : Quat) -> Quat {
  { w: -a.w, x: -a.x, y: -a.y, z: -a.z }
}

/// Q-7 實純量乘法
pub fn quat_scale(c : Double, a : Quat) -> Quat {
  { w: c * a.w, x: c * a.x, y: c * a.y, z: c * a.z }
}

/// Q-8 Hamilton 乘法
pub fn quat_mul(a : Quat, b : Quat) -> Quat {
  {
    w: a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z,
    x: a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y,
    y: a.w * b.y - a.x * b.z + a.y * b.w + a.z * b.x,
    z: a.w * b.z + a.x * b.y - a.y * b.x + a.z * b.w,
  }
}

/// Q-9 共軛
pub fn quat_conj(a : Quat) -> Quat {
  { w: a.w, x: -a.x, y: -a.y, z: -a.z }
}

/// Q-10 範數平方
pub fn quat_norm2(a : Quat) -> Double {
  a.w * a.w + a.x * a.x + a.y * a.y + a.z * a.z
}

/// Q-11 範數
pub fn quat_norm(a : Quat) -> Double {
  quat_norm2(a).sqrt()
}

/// Q-12 逆 q⁻¹ = q̄/‖q‖²
pub fn quat_inv(a : Quat) -> Quat {
  let d = quat_norm2(a)
  { w: a.w / d, x: -a.x / d, y: -a.y / d, z: -a.z / d }
}

/// Q-13 單位化
pub fn quat_unit(a : Quat) -> Quat {
  quat_scale(1.0 / quat_norm(a), a)
}

/// Q-14 由軸角構造單位四元數（旋轉角 θ 繞單位軸 u）
pub fn quat_from_axis_angle(axis : Array[Double], theta : Double) -> Quat {
  let n = (axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]).sqrt()
  let s = @math.sin(theta / 2.0)
  {
    w: @math.cos(theta / 2.0),
    x: axis[0] / n * s,
    y: axis[1] / n * s,
    z: axis[2] / n * s,
  }
}

/// Q-15 旋轉 R³ 向量：v ↦ q v q̄（q 為單位四元數）
pub fn quat_rotate_vec(q : Quat, v : Array[Double]) -> Array[Double] {
  let p = quat_imag(v[0], v[1], v[2])
  let r = quat_mul(quat_mul(q, p), quat_conj(q))
  [r.x, r.y, r.z]
}

/// Q-16 H 上之內積 ⟨q,p⟩ = Re(q̄p)
pub fn quat_inner(a : Quat, b : Quat) -> Double {
  quat_mul(quat_conj(a), b).w
}

/// H^n 之元素
pub struct HVec {
  comps : Array[Quat]
} derive(Debug, Eq)

/// Q-17 H^n 零向量
pub fn hvec_zero(n : Int) -> HVec {
  { comps: Array::make(n, quat_zero()) }
}

/// Q-18 H^n 加法
pub fn hvec_add(a : HVec, b : HVec) -> HVec {
  let out : Array[Quat] = []
  for i = 0; i < a.comps.length(); i = i + 1 {
    out.push(quat_add(a.comps[i], b.comps[i]))
  }
  { comps: out }
}

/// Q-19 H^n 實純量乘法（H^n 視為實向量空間）
pub fn hvec_scale(c : Double, a : HVec) -> HVec {
  let out : Array[Quat] = []
  for q in a.comps {
    out.push(quat_scale(c, q))
  }
  { comps: out }
}

/// Q-20 H^n 之實內積 Σ Re(q̄ᵢ pᵢ)
pub fn hvec_inner(a : HVec, b : HVec) -> Double {
  let mut s = 0.0
  for i = 0; i < a.comps.length(); i = i + 1 {
    s = s + quat_inner(a.comps[i], b.comps[i])
  }
  s
}

/// Q-21 H^n 作為實內積空間之字典（實維數 4n）
pub fn hvec_ip(n : Int) -> IPSp[HVec] {
  ips_make(
    {
      zero: hvec_zero(n),
      add: hvec_add,
      neg: fn(v) { hvec_scale(-1.0, v) },
      scale: hvec_scale,
    },
    hvec_inner,
  )
}

/// Q-22 H^n 之實維數 = 4n
pub fn hvec_real_dim(n : Int) -> Int {
  4 * n
}

/// Q-23 Gram–Schmidt 於 H^n（經由實內積字典）
pub fn hvec_orthonormal(vs : Array[HVec], eps : Double) -> Array[HVec] {
  ip_orthonormal(hvec_ip(vs[0].comps.length()), vs, eps)
}

/// Q-24 四元數旋轉複合：q₂∘q₁ 對應四元數乘法（注意次序）
pub fn quat_compose(q1 : Quat, q2 : Quat) -> Quat {
  quat_mul(q2, q1)
}
