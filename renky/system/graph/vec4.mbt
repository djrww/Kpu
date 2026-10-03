// 4D vector: homogeneous / clip-space coordinates.

///|
#valtype
pub(all) struct Vec4 {
  x : Double
  y : Double
  z : Double
  w : Double
} derive(Eq, Debug)

///|
pub fn Vec4::new(x : Double, y : Double, z : Double, w : Double) -> Vec4 {
  { x, y, z, w, }
}

///|
pub fn Vec4::zero() -> Vec4 {
  { x: 0.0, y: 0.0, z: 0.0, w: 0.0, }
}

///|
pub impl Add for Vec4 with fn add(a, b) {
  { x: a.x + b.x, y: a.y + b.y, z: a.z + b.z, w: a.w + b.w, }
}

///|
pub impl Sub for Vec4 with fn sub(a, b) {
  { x: a.x - b.x, y: a.y - b.y, z: a.z - b.z, w: a.w - b.w, }
}

///|
pub impl Neg for Vec4 with fn neg(a) {
  { x: -a.x, y: -a.y, z: -a.z, w: -a.w, }
}

///|
pub fn Vec4::scale(self : Vec4, s : Double) -> Vec4 {
  { x: self.x * s, y: self.y * s, z: self.z * s, w: self.w * s, }
}

///|
pub fn Vec4::dot(self : Vec4, o : Vec4) -> Double {
  self.x * o.x + self.y * o.y + self.z * o.z + self.w * o.w
}

///|
pub fn Vec4::xyz(self : Vec4) -> Vec3 {
  { x: self.x, y: self.y, z: self.z, }
}

///|
pub fn Vec4::lerp(self : Vec4, o : Vec4, t : Double) -> Vec4 {
  {
    x: self.x + (o.x - self.x) * t,
    y: self.y + (o.y - self.y) * t,
    z: self.z + (o.z - self.z) * t,
    w: self.w + (o.w - self.w) * t,
  }
}

///|
/// Perspective divide: (x/w, y/w, z/w). `None` when w == 0 / non-finite result.
pub fn Vec4::perspective_divide(self : Vec4) -> Vec3? {
  if self.w == 0.0 {
    return None
  }
  let iw = 1.0 / self.w
  Some({ x: self.x * iw, y: self.y * iw, z: self.z * iw, })
}

///|
pub fn Vec4::get(self : Vec4, i : Int) -> Double {
  match i {
    0 => self.x
    1 => self.y
    2 => self.z
    3 => self.w
    _ => abort("Vec4::get: index out of range")
  }
}

///|
pub fn Vec4::approx_eq(self : Vec4, o : Vec4, eps : Double) -> Bool {
  approx(self.x, o.x, eps) &&
  approx(self.y, o.y, eps) &&
  approx(self.z, o.z, eps) &&
  approx(self.w, o.w, eps)
}
