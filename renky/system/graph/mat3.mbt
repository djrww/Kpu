// 3x3 matrix. Naming: m<row><col>. Column vectors: v' = M * v.
// Uses: rotations, change of basis, normal matrix, tangent frames, 2D homogeneous affine
// transforms (UI / texture-space / map-space), homographies.

///|
#valtype
pub(all) struct Mat3 {
  m00 : Double
  m01 : Double
  m02 : Double
  m10 : Double
  m11 : Double
  m12 : Double
  m20 : Double
  m21 : Double
  m22 : Double
} derive(Eq, Debug)

///|
pub fn Mat3::identity() -> Mat3 {
  {
    m00: 1.0,
    m01: 0.0,
    m02: 0.0,
    m10: 0.0,
    m11: 1.0,
    m12: 0.0,
    m20: 0.0,
    m21: 0.0,
    m22: 1.0,
  }
}

///|
pub fn Mat3::from_cols(c0 : Vec3, c1 : Vec3, c2 : Vec3) -> Mat3 {
  {
    m00: c0.x,
    m01: c1.x,
    m02: c2.x,
    m10: c0.y,
    m11: c1.y,
    m12: c2.y,
    m20: c0.z,
    m21: c1.z,
    m22: c2.z,
  }
}

///|
pub fn Mat3::from_rows(r0 : Vec3, r1 : Vec3, r2 : Vec3) -> Mat3 {
  {
    m00: r0.x,
    m01: r0.y,
    m02: r0.z,
    m10: r1.x,
    m11: r1.y,
    m12: r1.z,
    m20: r2.x,
    m21: r2.y,
    m22: r2.z,
  }
}

///|
pub fn Mat3::col(self : Mat3, i : Int) -> Vec3 {
  match i {
    0 => { x: self.m00, y: self.m10, z: self.m20, }
    1 => { x: self.m01, y: self.m11, z: self.m21, }
    2 => { x: self.m02, y: self.m12, z: self.m22, }
    _ => abort("Mat3::col: index out of range")
  }
}

///|
pub fn Mat3::row(self : Mat3, i : Int) -> Vec3 {
  match i {
    0 => { x: self.m00, y: self.m01, z: self.m02, }
    1 => { x: self.m10, y: self.m11, z: self.m12, }
    2 => { x: self.m20, y: self.m21, z: self.m22, }
    _ => abort("Mat3::row: index out of range")
  }
}

///|
pub fn Mat3::diagonal(x : Double, y : Double, z : Double) -> Mat3 {
  { ..Mat3::identity(), m00: x, m11: y, m22: z, }
}

///|
/// Matrix product: `a * b` applies `b` first.
pub impl Mul for Mat3 with fn mul(a, b) {
  {
    m00: a.m00 * b.m00 + a.m01 * b.m10 + a.m02 * b.m20,
    m01: a.m00 * b.m01 + a.m01 * b.m11 + a.m02 * b.m21,
    m02: a.m00 * b.m02 + a.m01 * b.m12 + a.m02 * b.m22,
    m10: a.m10 * b.m00 + a.m11 * b.m10 + a.m12 * b.m20,
    m11: a.m10 * b.m01 + a.m11 * b.m11 + a.m12 * b.m21,
    m12: a.m10 * b.m02 + a.m11 * b.m12 + a.m12 * b.m22,
    m20: a.m20 * b.m00 + a.m21 * b.m10 + a.m22 * b.m20,
    m21: a.m20 * b.m01 + a.m21 * b.m11 + a.m22 * b.m21,
    m22: a.m20 * b.m02 + a.m21 * b.m12 + a.m22 * b.m22,
  }
}

///|
pub fn Mat3::mul_vec3(self : Mat3, v : Vec3) -> Vec3 {
  {
    x: self.m00 * v.x + self.m01 * v.y + self.m02 * v.z,
    y: self.m10 * v.x + self.m11 * v.y + self.m12 * v.z,
    z: self.m20 * v.x + self.m21 * v.y + self.m22 * v.z,
  }
}

///|
pub fn Mat3::transpose(self : Mat3) -> Mat3 {
  {
    m00: self.m00,
    m01: self.m10,
    m02: self.m20,
    m10: self.m01,
    m11: self.m11,
    m12: self.m21,
    m20: self.m02,
    m21: self.m12,
    m22: self.m22,
  }
}

///|
pub fn Mat3::determinant(self : Mat3) -> Double {
  self.m00 * (self.m11 * self.m22 - self.m12 * self.m21) -
  self.m01 * (self.m10 * self.m22 - self.m12 * self.m20) +
  self.m02 * (self.m10 * self.m21 - self.m11 * self.m20)
}

///|
/// General inverse; `None` when singular (det == 0) or not finite.
pub fn Mat3::inverse(self : Mat3) -> Mat3? {
  let c00 = self.m11 * self.m22 - self.m12 * self.m21
  let c01 = self.m12 * self.m20 - self.m10 * self.m22
  let c02 = self.m10 * self.m21 - self.m11 * self.m20
  let det = self.m00 * c00 + self.m01 * c01 + self.m02 * c02
  if det == 0.0 || det.is_nan() || det.is_inf() {
    return None
  }
  let i = 1.0 / det
  Some({
    m00: c00 * i,
    m01: (self.m02 * self.m21 - self.m01 * self.m22) * i,
    m02: (self.m01 * self.m12 - self.m02 * self.m11) * i,
    m10: c01 * i,
    m11: (self.m00 * self.m22 - self.m02 * self.m20) * i,
    m12: (self.m02 * self.m10 - self.m00 * self.m12) * i,
    m20: c02 * i,
    m21: (self.m01 * self.m20 - self.m00 * self.m21) * i,
    m22: (self.m00 * self.m11 - self.m01 * self.m10) * i,
  })
}

///|
/// (M^-1)^T, the matrix that transforms surface normals. `None` when singular.
pub fn Mat3::normal_matrix(self : Mat3) -> Mat3? {
  match self.inverse() {
    Some(i) => Some(i.transpose())
    None => None
  }
}

///|
pub fn Mat3::from_quat(q : Quat) -> Mat3 {
  let xx = q.x * q.x
  let yy = q.y * q.y
  let zz = q.z * q.z
  let xy = q.x * q.y
  let xz = q.x * q.z
  let yz = q.y * q.z
  let wx = q.w * q.x
  let wy = q.w * q.y
  let wz = q.w * q.z
  {
    m00: 1.0 - 2.0 * (yy + zz),
    m01: 2.0 * (xy - wz),
    m02: 2.0 * (xz + wy),
    m10: 2.0 * (xy + wz),
    m11: 1.0 - 2.0 * (xx + zz),
    m12: 2.0 * (yz - wx),
    m20: 2.0 * (xz - wy),
    m21: 2.0 * (yz + wx),
    m22: 1.0 - 2.0 * (xx + yy),
  }
}

///|
pub fn Mat3::rotation_x(a : Double) -> Mat3 {
  let c = @math.cos(a)
  let s = @math.sin(a)
  { ..Mat3::identity(), m11: c, m12: -s, m21: s, m22: c, }
}

///|
pub fn Mat3::rotation_y(a : Double) -> Mat3 {
  let c = @math.cos(a)
  let s = @math.sin(a)
  { ..Mat3::identity(), m00: c, m02: s, m20: -s, m22: c, }
}

///|
pub fn Mat3::rotation_z(a : Double) -> Mat3 {
  let c = @math.cos(a)
  let s = @math.sin(a)
  { ..Mat3::identity(), m00: c, m01: -s, m10: s, m11: c, }
}

///|
/// Rodrigues rotation about `axis` (normalised internally).
pub fn Mat3::from_axis_angle(axis : Vec3, angle : Double) -> Mat3 {
  match axis.try_normalize() {
    None => Mat3::identity()
    Some(a) => {
      let c = @math.cos(angle)
      let s = @math.sin(angle)
      let t = 1.0 - c
      {
        m00: c + a.x * a.x * t,
        m01: a.x * a.y * t - a.z * s,
        m02: a.x * a.z * t + a.y * s,
        m10: a.y * a.x * t + a.z * s,
        m11: c + a.y * a.y * t,
        m12: a.y * a.z * t - a.x * s,
        m20: a.z * a.x * t - a.y * s,
        m21: a.z * a.y * t + a.x * s,
        m22: c + a.z * a.z * t,
      }
    }
  }
}

///|
/// Orientation matrix (columns = right, up, back) whose local -Z is `forward` and whose
/// local +Y is closest to `up`. Degenerate inputs fall back to a valid orthonormal basis.
pub fn Mat3::look_rotation(forward : Vec3, up : Vec3) -> Mat3 {
  let z = (-forward).normalize()
  let z = if z.length_sq() == 0.0 { Vec3::unit_z() } else { z }
  let mut x = up.cross(z)
  if x.length_sq() < 1.0e-24 {
    x = z.any_perpendicular()
  } else {
    x = x.normalize()
  }
  let y = z.cross(x)
  Mat3::from_cols(x, y, z)
}

///|
pub fn Mat3::approx_eq(self : Mat3, o : Mat3, eps : Double) -> Bool {
  approx(self.m00, o.m00, eps) &&
  approx(self.m01, o.m01, eps) &&
  approx(self.m02, o.m02, eps) &&
  approx(self.m10, o.m10, eps) &&
  approx(self.m11, o.m11, eps) &&
  approx(self.m12, o.m12, eps) &&
  approx(self.m20, o.m20, eps) &&
  approx(self.m21, o.m21, eps) &&
  approx(self.m22, o.m22, eps)
}

///|
/// Column-major 9 values (GLSL `mat3` / WGSL `mat3x3` memory order, tightly packed).
pub fn Mat3::to_col_major(self : Mat3) -> FixedArray[Double] {
  [
    self.m00,
    self.m10,
    self.m20,
    self.m01,
    self.m11,
    self.m21,
    self.m02,
    self.m12,
    self.m22,
  ]
}

// ---------------------------------------------------------------------------
// 2D homogeneous affine / projective transforms (points are (x, y, 1)).
// ---------------------------------------------------------------------------

///|
pub fn Mat3::translation2d(tx : Double, ty : Double) -> Mat3 {
  { ..Mat3::identity(), m02: tx, m12: ty, }
}

///|
pub fn Mat3::scaling2d(sx : Double, sy : Double) -> Mat3 {
  { ..Mat3::identity(), m00: sx, m11: sy, }
}

///|
/// Counter-clockwise rotation about the origin (y-up frame).
pub fn Mat3::rotation2d(angle : Double) -> Mat3 {
  Mat3::rotation_z(angle)
}

///|
/// Translate * Rotate * Scale (scale applied first).
pub fn Mat3::from_trs2d(t : Vec2, angle : Double, s : Vec2) -> Mat3 {
  let c = @math.cos(angle)
  let sn = @math.sin(angle)
  {
    m00: c * s.x,
    m01: -sn * s.y,
    m02: t.x,
    m10: sn * s.x,
    m11: c * s.y,
    m12: t.y,
    m20: 0.0,
    m21: 0.0,
    m22: 1.0,
  }
}

///|
/// Affine point transform (assumes last row = 0 0 1).
pub fn Mat3::mul_point2(self : Mat3, p : Vec2) -> Vec2 {
  {
    x: self.m00 * p.x + self.m01 * p.y + self.m02,
    y: self.m10 * p.x + self.m11 * p.y + self.m12,
  }
}

///|
/// Direction transform (ignores translation).
pub fn Mat3::mul_dir2(self : Mat3, d : Vec2) -> Vec2 {
  { x: self.m00 * d.x + self.m01 * d.y, y: self.m10 * d.x + self.m11 * d.y, }
}

///|
/// Projective (homography) point transform with divide by w. `None` if w == 0.
pub fn Mat3::project_point2(self : Mat3, p : Vec2) -> Vec2? {
  let w = self.m20 * p.x + self.m21 * p.y + self.m22
  if w == 0.0 {
    return None
  }
  let iw = 1.0 / w
  Some({
    x: (self.m00 * p.x + self.m01 * p.y + self.m02) * iw,
    y: (self.m10 * p.x + self.m11 * p.y + self.m12) * iw,
  })
}
