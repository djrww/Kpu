// 4x4 matrix. Naming: m<row><col>. Column vectors: v' = M * v.  `a * b` applies `b` first.
// Typical chain:  clip = Projection * View * Model * local_position

///|
#valtype
pub(all) struct Mat4 {
  m00 : Double
  m01 : Double
  m02 : Double
  m03 : Double
  m10 : Double
  m11 : Double
  m12 : Double
  m13 : Double
  m20 : Double
  m21 : Double
  m22 : Double
  m23 : Double
  m30 : Double
  m31 : Double
  m32 : Double
  m33 : Double
} derive(Eq, Debug)

///|
pub fn Mat4::identity() -> Mat4 {
  {
    m00: 1.0,
    m01: 0.0,
    m02: 0.0,
    m03: 0.0,
    m10: 0.0,
    m11: 1.0,
    m12: 0.0,
    m13: 0.0,
    m20: 0.0,
    m21: 0.0,
    m22: 1.0,
    m23: 0.0,
    m30: 0.0,
    m31: 0.0,
    m32: 0.0,
    m33: 1.0,
  }
}

///|
pub fn Mat4::from_cols(c0 : Vec4, c1 : Vec4, c2 : Vec4, c3 : Vec4) -> Mat4 {
  {
    m00: c0.x,
    m01: c1.x,
    m02: c2.x,
    m03: c3.x,
    m10: c0.y,
    m11: c1.y,
    m12: c2.y,
    m13: c3.y,
    m20: c0.z,
    m21: c1.z,
    m22: c2.z,
    m23: c3.z,
    m30: c0.w,
    m31: c1.w,
    m32: c2.w,
    m33: c3.w,
  }
}

///|
pub fn Mat4::col(self : Mat4, i : Int) -> Vec4 {
  match i {
    0 => { x: self.m00, y: self.m10, z: self.m20, w: self.m30, }
    1 => { x: self.m01, y: self.m11, z: self.m21, w: self.m31, }
    2 => { x: self.m02, y: self.m12, z: self.m22, w: self.m32, }
    3 => { x: self.m03, y: self.m13, z: self.m23, w: self.m33, }
    _ => abort("Mat4::col: index out of range")
  }
}

///|
pub fn Mat4::row(self : Mat4, i : Int) -> Vec4 {
  match i {
    0 => { x: self.m00, y: self.m01, z: self.m02, w: self.m03, }
    1 => { x: self.m10, y: self.m11, z: self.m12, w: self.m13, }
    2 => { x: self.m20, y: self.m21, z: self.m22, w: self.m23, }
    3 => { x: self.m30, y: self.m31, z: self.m32, w: self.m33, }
    _ => abort("Mat4::row: index out of range")
  }
}

///|
/// Element (row, col).
pub fn Mat4::get(self : Mat4, r : Int, c : Int) -> Double {
  self.row(r).get(c)
}

///|
/// Matrix product: `a * b` applies `b` first.
pub impl Mul for Mat4 with fn mul(a, b) {
  {
    m00: a.m00 * b.m00 + a.m01 * b.m10 + a.m02 * b.m20 + a.m03 * b.m30,
    m01: a.m00 * b.m01 + a.m01 * b.m11 + a.m02 * b.m21 + a.m03 * b.m31,
    m02: a.m00 * b.m02 + a.m01 * b.m12 + a.m02 * b.m22 + a.m03 * b.m32,
    m03: a.m00 * b.m03 + a.m01 * b.m13 + a.m02 * b.m23 + a.m03 * b.m33,
    m10: a.m10 * b.m00 + a.m11 * b.m10 + a.m12 * b.m20 + a.m13 * b.m30,
    m11: a.m10 * b.m01 + a.m11 * b.m11 + a.m12 * b.m21 + a.m13 * b.m31,
    m12: a.m10 * b.m02 + a.m11 * b.m12 + a.m12 * b.m22 + a.m13 * b.m32,
    m13: a.m10 * b.m03 + a.m11 * b.m13 + a.m12 * b.m23 + a.m13 * b.m33,
    m20: a.m20 * b.m00 + a.m21 * b.m10 + a.m22 * b.m20 + a.m23 * b.m30,
    m21: a.m20 * b.m01 + a.m21 * b.m11 + a.m22 * b.m21 + a.m23 * b.m31,
    m22: a.m20 * b.m02 + a.m21 * b.m12 + a.m22 * b.m22 + a.m23 * b.m32,
    m23: a.m20 * b.m03 + a.m21 * b.m13 + a.m22 * b.m23 + a.m23 * b.m33,
    m30: a.m30 * b.m00 + a.m31 * b.m10 + a.m32 * b.m20 + a.m33 * b.m30,
    m31: a.m30 * b.m01 + a.m31 * b.m11 + a.m32 * b.m21 + a.m33 * b.m31,
    m32: a.m30 * b.m02 + a.m31 * b.m12 + a.m32 * b.m22 + a.m33 * b.m32,
    m33: a.m30 * b.m03 + a.m31 * b.m13 + a.m32 * b.m23 + a.m33 * b.m33,
  }
}

///|
pub fn Mat4::mul_vec4(self : Mat4, v : Vec4) -> Vec4 {
  {
    x: self.m00 * v.x + self.m01 * v.y + self.m02 * v.z + self.m03 * v.w,
    y: self.m10 * v.x + self.m11 * v.y + self.m12 * v.z + self.m13 * v.w,
    z: self.m20 * v.x + self.m21 * v.y + self.m22 * v.z + self.m23 * v.w,
    w: self.m30 * v.x + self.m31 * v.y + self.m32 * v.z + self.m33 * v.w,
  }
}

///|
/// Transform a point assuming an AFFINE matrix (last row 0 0 0 1): w stays 1, no divide.
/// Fastest path for model/view matrices.
pub fn Mat4::mul_point(self : Mat4, p : Vec3) -> Vec3 {
  {
    x: self.m00 * p.x + self.m01 * p.y + self.m02 * p.z + self.m03,
    y: self.m10 * p.x + self.m11 * p.y + self.m12 * p.z + self.m13,
    z: self.m20 * p.x + self.m21 * p.y + self.m22 * p.z + self.m23,
  }
}

///|
/// Transform a direction (w = 0): translation ignored. For normals use `normal_matrix`.
pub fn Mat4::mul_dir(self : Mat4, d : Vec3) -> Vec3 {
  {
    x: self.m00 * d.x + self.m01 * d.y + self.m02 * d.z,
    y: self.m10 * d.x + self.m11 * d.y + self.m12 * d.z,
    z: self.m20 * d.x + self.m21 * d.y + self.m22 * d.z,
  }
}

///|
/// Full homogeneous transform of a point (w = 1) followed by the perspective divide.
/// `None` if the resulting w == 0. NOTE: points behind the eye (w < 0) are returned
/// mirrored -- use `Camera::world_to_ndc` for visibility-aware projection.
pub fn Mat4::project_point(self : Mat4, p : Vec3) -> Vec3? {
  self.mul_vec4(p.to_vec4(1.0)).perspective_divide()
}

///|
pub fn Mat4::transpose(self : Mat4) -> Mat4 {
  {
    m00: self.m00,
    m01: self.m10,
    m02: self.m20,
    m03: self.m30,
    m10: self.m01,
    m11: self.m11,
    m12: self.m21,
    m13: self.m31,
    m20: self.m02,
    m21: self.m12,
    m22: self.m22,
    m23: self.m32,
    m30: self.m03,
    m31: self.m13,
    m32: self.m23,
    m33: self.m33,
  }
}

///|
pub fn Mat4::translation(t : Vec3) -> Mat4 {
  { ..Mat4::identity(), m03: t.x, m13: t.y, m23: t.z, }
}

///|
pub fn Mat4::scaling(s : Vec3) -> Mat4 {
  { ..Mat4::identity(), m00: s.x, m11: s.y, m22: s.z, }
}

///|
pub fn Mat4::from_mat3(r : Mat3) -> Mat4 {
  {
    m00: r.m00,
    m01: r.m01,
    m02: r.m02,
    m03: 0.0,
    m10: r.m10,
    m11: r.m11,
    m12: r.m12,
    m13: 0.0,
    m20: r.m20,
    m21: r.m21,
    m22: r.m22,
    m23: 0.0,
    m30: 0.0,
    m31: 0.0,
    m32: 0.0,
    m33: 1.0,
  }
}

///|
/// Upper-left 3x3 block (rotation * scale for an affine matrix).
pub fn Mat4::to_mat3(self : Mat4) -> Mat3 {
  {
    m00: self.m00,
    m01: self.m01,
    m02: self.m02,
    m10: self.m10,
    m11: self.m11,
    m12: self.m12,
    m20: self.m20,
    m21: self.m21,
    m22: self.m22,
  }
}

///|
pub fn Mat4::rotation(q : Quat) -> Mat4 {
  Mat4::from_mat3(Mat3::from_quat(q))
}

///|
pub fn Mat4::translation_part(self : Mat4) -> Vec3 {
  { x: self.m03, y: self.m13, z: self.m23, }
}

///|
/// Translate * Rotate * Scale (scale applied first). The standard model-matrix recipe.
pub fn Mat4::from_trs(t : Vec3, r : Quat, s : Vec3) -> Mat4 {
  let q = Mat3::from_quat(r)
  {
    m00: q.m00 * s.x,
    m01: q.m01 * s.y,
    m02: q.m02 * s.z,
    m03: t.x,
    m10: q.m10 * s.x,
    m11: q.m11 * s.y,
    m12: q.m12 * s.z,
    m13: t.y,
    m20: q.m20 * s.x,
    m21: q.m21 * s.y,
    m22: q.m22 * s.z,
    m23: t.z,
    m30: 0.0,
    m31: 0.0,
    m32: 0.0,
    m33: 1.0,
  }
}

///|
/// True when the last row is exactly (0, 0, 0, 1).
pub fn Mat4::is_affine(self : Mat4) -> Bool {
  self.m30 == 0.0 && self.m31 == 0.0 && self.m32 == 0.0 && self.m33 == 1.0
}

///|
pub fn Mat4::determinant(self : Mat4) -> Double {
  let s0 = self.m00 * self.m11 - self.m10 * self.m01
  let s1 = self.m00 * self.m12 - self.m10 * self.m02
  let s2 = self.m00 * self.m13 - self.m10 * self.m03
  let s3 = self.m01 * self.m12 - self.m11 * self.m02
  let s4 = self.m01 * self.m13 - self.m11 * self.m03
  let s5 = self.m02 * self.m13 - self.m12 * self.m03
  let c5 = self.m22 * self.m33 - self.m32 * self.m23
  let c4 = self.m21 * self.m33 - self.m31 * self.m23
  let c3 = self.m21 * self.m32 - self.m31 * self.m22
  let c2 = self.m20 * self.m33 - self.m30 * self.m23
  let c1 = self.m20 * self.m32 - self.m30 * self.m22
  let c0 = self.m20 * self.m31 - self.m30 * self.m21
  s0 * c5 - s1 * c4 + s2 * c3 + s3 * c2 - s4 * c1 + s5 * c0
}

///|
/// General 4x4 inverse (Laplace expansion on 2x2 minors). `None` if singular / not finite.
/// Prefer `affine_inverse` / `rigid_inverse` when the matrix structure is known.
pub fn Mat4::inverse(self : Mat4) -> Mat4? {
  let a = self
  let s0 = a.m00 * a.m11 - a.m10 * a.m01
  let s1 = a.m00 * a.m12 - a.m10 * a.m02
  let s2 = a.m00 * a.m13 - a.m10 * a.m03
  let s3 = a.m01 * a.m12 - a.m11 * a.m02
  let s4 = a.m01 * a.m13 - a.m11 * a.m03
  let s5 = a.m02 * a.m13 - a.m12 * a.m03
  let c5 = a.m22 * a.m33 - a.m32 * a.m23
  let c4 = a.m21 * a.m33 - a.m31 * a.m23
  let c3 = a.m21 * a.m32 - a.m31 * a.m22
  let c2 = a.m20 * a.m33 - a.m30 * a.m23
  let c1 = a.m20 * a.m32 - a.m30 * a.m22
  let c0 = a.m20 * a.m31 - a.m30 * a.m21
  let det = s0 * c5 - s1 * c4 + s2 * c3 + s3 * c2 - s4 * c1 + s5 * c0
  if det == 0.0 || det.is_nan() || det.is_inf() {
    return None
  }
  let i = 1.0 / det
  Some({
    m00: (a.m11 * c5 - a.m12 * c4 + a.m13 * c3) * i,
    m01: (-a.m01 * c5 + a.m02 * c4 - a.m03 * c3) * i,
    m02: (a.m31 * s5 - a.m32 * s4 + a.m33 * s3) * i,
    m03: (-a.m21 * s5 + a.m22 * s4 - a.m23 * s3) * i,
    m10: (-a.m10 * c5 + a.m12 * c2 - a.m13 * c1) * i,
    m11: (a.m00 * c5 - a.m02 * c2 + a.m03 * c1) * i,
    m12: (-a.m30 * s5 + a.m32 * s2 - a.m33 * s1) * i,
    m13: (a.m20 * s5 - a.m22 * s2 + a.m23 * s1) * i,
    m20: (a.m10 * c4 - a.m11 * c2 + a.m13 * c0) * i,
    m21: (-a.m00 * c4 + a.m01 * c2 - a.m03 * c0) * i,
    m22: (a.m30 * s4 - a.m31 * s2 + a.m33 * s0) * i,
    m23: (-a.m20 * s4 + a.m21 * s2 - a.m23 * s0) * i,
    m30: (-a.m10 * c3 + a.m11 * c1 - a.m12 * c0) * i,
    m31: (a.m00 * c3 - a.m01 * c1 + a.m02 * c0) * i,
    m32: (-a.m30 * s3 + a.m31 * s1 - a.m32 * s0) * i,
    m33: (a.m20 * s3 - a.m21 * s1 + a.m22 * s0) * i,
  })
}

///|
/// Inverse of an affine matrix (last row 0 0 0 1) via the 3x3 block: ~3x cheaper than
/// `inverse`. `None` if the 3x3 block is singular. The last row is NOT checked.
pub fn Mat4::affine_inverse(self : Mat4) -> Mat4? {
  match self.to_mat3().inverse() {
    None => None
    Some(i) => {
      let t = Vec3::new(self.m03, self.m13, self.m23)
      let it = i.mul_vec3(t)
      Some({
        m00: i.m00,
        m01: i.m01,
        m02: i.m02,
        m03: -it.x,
        m10: i.m10,
        m11: i.m11,
        m12: i.m12,
        m13: -it.y,
        m20: i.m20,
        m21: i.m21,
        m22: i.m22,
        m23: -it.z,
        m30: 0.0,
        m31: 0.0,
        m32: 0.0,
        m33: 1.0,
      })
    }
  }
}

///|
/// Inverse of a rigid transform (rotation + translation only): R^T, -R^T t. Cheapest and
/// exact; the caller guarantees orthonormality (e.g. camera pose, bone transforms).
pub fn Mat4::rigid_inverse(self : Mat4) -> Mat4 {
  let tx = self.m03
  let ty = self.m13
  let tz = self.m23
  {
    m00: self.m00,
    m01: self.m10,
    m02: self.m20,
    m03: -(self.m00 * tx + self.m10 * ty + self.m20 * tz),
    m10: self.m01,
    m11: self.m11,
    m12: self.m21,
    m13: -(self.m01 * tx + self.m11 * ty + self.m21 * tz),
    m20: self.m02,
    m21: self.m12,
    m22: self.m22,
    m23: -(self.m02 * tx + self.m12 * ty + self.m22 * tz),
    m30: 0.0,
    m31: 0.0,
    m32: 0.0,
    m33: 1.0,
  }
}

///|
/// Matrix for transforming surface normals of a model matrix: inverse-transpose of the 3x3
/// block. Correct under non-uniform scale and mirroring. `None` if singular.
/// Result is NOT normalised; normalise the transformed normal.
pub fn Mat4::normal_matrix(self : Mat4) -> Mat3? {
  self.to_mat3().normal_matrix()
}

///|
/// Split an affine TRS matrix into (translation, rotation, scale). A negative determinant
/// (mirroring) is reported as a negative X scale. Returns `None` for projective matrices,
/// zero scale, or when shear is present (columns not mutually orthogonal).
pub fn Mat4::decompose(self : Mat4) -> Transform? {
  if !self.is_affine() {
    return None
  }
  let c0 = Vec3::new(self.m00, self.m10, self.m20)
  let c1 = Vec3::new(self.m01, self.m11, self.m21)
  let c2 = Vec3::new(self.m02, self.m12, self.m22)
  let l0 = c0.length()
  let l1 = c1.length()
  let l2 = c2.length()
  if l0 == 0.0 || l1 == 0.0 || l2 == 0.0 {
    return None
  }
  // Orthogonality (no shear) test, relative to column lengths.
  let tol = 1.0e-7
  if c0.dot(c1).abs() > tol * l0 * l1 ||
    c0.dot(c2).abs() > tol * l0 * l2 ||
    c1.dot(c2).abs() > tol * l1 * l2 {
    return None
  }
  let det = self.to_mat3().determinant()
  let sx = if det < 0.0 { -l0 } else { l0 }
  let r = Mat3::from_cols(
    c0.scale(1.0 / sx),
    c1.scale(1.0 / l1),
    c2.scale(1.0 / l2),
  )
  Some({
    translation: self.translation_part(),
    rotation: Quat::from_mat3(r),
    scale: Vec3::new(sx, l1, l2),
  })
}

///|
/// View matrix (world -> view) of an eye at `eye` looking at `target`, `up` = approximate
/// up. Right-handed, camera looks down -Z, +Y up, +X right.
pub fn Mat4::look_at(eye : Vec3, target : Vec3, up : Vec3) -> Mat4 {
  let r = Mat3::look_rotation(target - eye, up)
  // camera-to-world = [r | eye]; view = its rigid inverse.
  Mat4::rigid_inverse({
    ..Mat4::from_mat3(r),
    m03: eye.x,
    m13: eye.y,
    m23: eye.z,
  })
}

///|
pub fn Mat4::approx_eq(self : Mat4, o : Mat4, eps : Double) -> Bool {
  approx(self.m00, o.m00, eps) &&
  approx(self.m01, o.m01, eps) &&
  approx(self.m02, o.m02, eps) &&
  approx(self.m03, o.m03, eps) &&
  approx(self.m10, o.m10, eps) &&
  approx(self.m11, o.m11, eps) &&
  approx(self.m12, o.m12, eps) &&
  approx(self.m13, o.m13, eps) &&
  approx(self.m20, o.m20, eps) &&
  approx(self.m21, o.m21, eps) &&
  approx(self.m22, o.m22, eps) &&
  approx(self.m23, o.m23, eps) &&
  approx(self.m30, o.m30, eps) &&
  approx(self.m31, o.m31, eps) &&
  approx(self.m32, o.m32, eps) &&
  approx(self.m33, o.m33, eps)
}

///|
/// 16 values in column-major order (GLSL / WGSL / Metal / Vulkan upload order).
pub fn Mat4::to_col_major(self : Mat4) -> FixedArray[Double] {
  [
    self.m00,
    self.m10,
    self.m20,
    self.m30,
    self.m01,
    self.m11,
    self.m21,
    self.m31,
    self.m02,
    self.m12,
    self.m22,
    self.m32,
    self.m03,
    self.m13,
    self.m23,
    self.m33,
  ]
}

///|
pub fn Mat4::from_col_major(a : FixedArray[Double]) -> Mat4 {
  if a.length() < 16 {
    abort("Mat4::from_col_major: need 16 values")
  }
  {
    m00: a[0],
    m10: a[1],
    m20: a[2],
    m30: a[3],
    m01: a[4],
    m11: a[5],
    m21: a[6],
    m31: a[7],
    m02: a[8],
    m12: a[9],
    m22: a[10],
    m32: a[11],
    m03: a[12],
    m13: a[13],
    m23: a[14],
    m33: a[15],
  }
}

///|
/// Column-major single-precision values for GPU upload.
pub fn Mat4::to_col_major_f32(self : Mat4) -> FixedArray[Float] {
  let d = self.to_col_major()
  let out : FixedArray[Float] = FixedArray::make(16, (0.0 : Float))
  for i in 0..<16 {
    out[i] = Float::from_double(d[i])
  }
  out
}

///|
/// Transform `count` points stored as packed xyz triples in `src` (affine matrix, w = 1)
/// writing to `dst` (may alias `src`). No per-point allocation.
pub fn Mat4::transform_points(
  self : Mat4,
  src : FixedArray[Double],
  dst : FixedArray[Double],
  count : Int,
) -> Unit {
  if src.length() < count * 3 || dst.length() < count * 3 {
    abort("Mat4::transform_points: buffer too small")
  }
  for i in 0..<count {
    let o = i * 3
    let x = src[o]
    let y = src[o + 1]
    let z = src[o + 2]
    dst[o] = self.m00 * x + self.m01 * y + self.m02 * z + self.m03
    dst[o + 1] = self.m10 * x + self.m11 * y + self.m12 * z + self.m13
    dst[o + 2] = self.m20 * x + self.m21 * y + self.m22 * z + self.m23
  }
}

///|
/// Single-precision variant of `transform_points`, math done in double then rounded once.
pub fn Mat4::transform_points_f32(
  self : Mat4,
  src : FixedArray[Float],
  dst : FixedArray[Float],
  count : Int,
) -> Unit {
  if src.length() < count * 3 || dst.length() < count * 3 {
    abort("Mat4::transform_points_f32: buffer too small")
  }
  for i in 0..<count {
    let o = i * 3
    let x = src[o].to_double()
    let y = src[o + 1].to_double()
    let z = src[o + 2].to_double()
    dst[o] = Float::from_double(
      self.m00 * x + self.m01 * y + self.m02 * z + self.m03,
    )
    dst[o + 1] = Float::from_double(
      self.m10 * x + self.m11 * y + self.m12 * z + self.m13,
    )
    dst[o + 2] = Float::from_double(
      self.m20 * x + self.m21 * y + self.m22 * z + self.m23,
    )
  }
}
