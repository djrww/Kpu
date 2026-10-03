// 3D vector: points, directions, normals, colors-as-xyz, euclidean coordinates in every space.

///|
#valtype
pub(all) struct Vec3 {
  x : Double
  y : Double
  z : Double
} derive(Eq, Debug)

///|
pub fn Vec3::new(x : Double, y : Double, z : Double) -> Vec3 {
  { x, y, z, }
}

///|
pub fn Vec3::zero() -> Vec3 {
  { x: 0.0, y: 0.0, z: 0.0, }
}

///|
pub fn Vec3::splat(v : Double) -> Vec3 {
  { x: v, y: v, z: v, }
}

///|
pub fn Vec3::unit_x() -> Vec3 {
  { x: 1.0, y: 0.0, z: 0.0, }
}

///|
pub fn Vec3::unit_y() -> Vec3 {
  { x: 0.0, y: 1.0, z: 0.0, }
}

///|
pub fn Vec3::unit_z() -> Vec3 {
  { x: 0.0, y: 0.0, z: 1.0, }
}

///|
pub impl Add for Vec3 with fn add(a, b) {
  { x: a.x + b.x, y: a.y + b.y, z: a.z + b.z, }
}

///|
pub impl Sub for Vec3 with fn sub(a, b) {
  { x: a.x - b.x, y: a.y - b.y, z: a.z - b.z, }
}

///|
pub impl Neg for Vec3 with fn neg(a) {
  { x: -a.x, y: -a.y, z: -a.z, }
}

///|
/// Component-wise product. For a scalar multiple use `scale`.
pub impl Mul for Vec3 with fn mul(a, b) {
  { x: a.x * b.x, y: a.y * b.y, z: a.z * b.z, }
}

///|
pub fn Vec3::scale(self : Vec3, s : Double) -> Vec3 {
  { x: self.x * s, y: self.y * s, z: self.z * s, }
}

///|
pub fn Vec3::dot(self : Vec3, o : Vec3) -> Double {
  self.x * o.x + self.y * o.y + self.z * o.z
}

///|
/// Right-handed cross product.
pub fn Vec3::cross(self : Vec3, o : Vec3) -> Vec3 {
  {
    x: self.y * o.z - self.z * o.y,
    y: self.z * o.x - self.x * o.z,
    z: self.x * o.y - self.y * o.x,
  }
}

///|
pub fn Vec3::length_sq(self : Vec3) -> Double {
  self.x * self.x + self.y * self.y + self.z * self.z
}

///|
pub fn Vec3::length(self : Vec3) -> Double {
  self.length_sq().sqrt()
}

///|
pub fn Vec3::distance(self : Vec3, o : Vec3) -> Double {
  (self - o).length()
}

///|
/// Unit vector, or `None` when the length is zero / not finite.
pub fn Vec3::try_normalize(self : Vec3) -> Vec3? {
  let l = self.length()
  if l > 0.0 && l.is_inf() == false {
    Some(self.scale(1.0 / l))
  } else {
    None
  }
}

///|
/// Unit vector; the zero vector is returned unchanged.
pub fn Vec3::normalize(self : Vec3) -> Vec3 {
  self.try_normalize().unwrap_or(self)
}

///|
pub fn Vec3::lerp(self : Vec3, o : Vec3, t : Double) -> Vec3 {
  {
    x: self.x + (o.x - self.x) * t,
    y: self.y + (o.y - self.y) * t,
    z: self.z + (o.z - self.z) * t,
  }
}

///|
pub fn Vec3::min(self : Vec3, o : Vec3) -> Vec3 {
  { x: self.x.min(o.x), y: self.y.min(o.y), z: self.z.min(o.z), }
}

///|
pub fn Vec3::max(self : Vec3, o : Vec3) -> Vec3 {
  { x: self.x.max(o.x), y: self.y.max(o.y), z: self.z.max(o.z), }
}

///|
pub fn Vec3::abs(self : Vec3) -> Vec3 {
  { x: self.x.abs(), y: self.y.abs(), z: self.z.abs(), }
}

///|
/// Component by index 0..2 (aborts otherwise).
pub fn Vec3::get(self : Vec3, i : Int) -> Double {
  match i {
    0 => self.x
    1 => self.y
    2 => self.z
    _ => abort("Vec3::get: index out of range")
  }
}

///|
/// Reflect `self` about unit normal `n`.
pub fn Vec3::reflect(self : Vec3, n : Vec3) -> Vec3 {
  self - n.scale(2.0 * self.dot(n))
}

///|
/// Component of `self` along `onto` (need not be unit). Zero if `onto` is zero.
pub fn Vec3::project_onto(self : Vec3, onto : Vec3) -> Vec3 {
  let d = onto.length_sq()
  if d == 0.0 {
    Vec3::zero()
  } else {
    onto.scale(self.dot(onto) / d)
  }
}

///|
/// Unsigned angle in [0, PI] between two vectors (robust for tiny angles).
pub fn Vec3::angle_to(self : Vec3, o : Vec3) -> Double {
  @math.atan2(self.cross(o).length(), self.dot(o))
}

///|
/// A vector perpendicular to `self` (unit length when `self` is non-zero).
pub fn Vec3::any_perpendicular(self : Vec3) -> Vec3 {
  let a = self.abs()
  let helper = if a.x <= a.y && a.x <= a.z {
    Vec3::unit_x()
  } else if a.y <= a.z {
    Vec3::unit_y()
  } else {
    Vec3::unit_z()
  }
  self.cross(helper).normalize()
}

///|
pub fn Vec3::xy(self : Vec3) -> Vec2 {
  { x: self.x, y: self.y, }
}

///|
pub fn Vec3::to_vec4(self : Vec3, w : Double) -> Vec4 {
  { x: self.x, y: self.y, z: self.z, w, }
}

///|
pub fn Vec3::approx_eq(self : Vec3, o : Vec3, eps : Double) -> Bool {
  approx(self.x, o.x, eps) &&
  approx(self.y, o.y, eps) &&
  approx(self.z, o.z, eps)
}

///|
pub fn Vec3::is_finite(self : Vec3) -> Bool {
  !(self.x.is_nan() ||
  self.x.is_inf() ||
  self.y.is_nan() ||
  self.y.is_inf() ||
  self.z.is_nan() ||
  self.z.is_inf())
}
