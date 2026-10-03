// 2D vector. Also used for: UV / texel coords, pixel coords, NDC xy, polar->cartesian, map meters.

///|
#valtype
pub(all) struct Vec2 {
  x : Double
  y : Double
} derive(Eq, Debug)

///|
pub fn Vec2::new(x : Double, y : Double) -> Vec2 {
  { x, y, }
}

///|
pub fn Vec2::zero() -> Vec2 {
  { x: 0.0, y: 0.0, }
}

///|
pub fn Vec2::splat(v : Double) -> Vec2 {
  { x: v, y: v, }
}

///|
pub impl Add for Vec2 with fn add(a, b) {
  { x: a.x + b.x, y: a.y + b.y, }
}

///|
pub impl Sub for Vec2 with fn sub(a, b) {
  { x: a.x - b.x, y: a.y - b.y, }
}

///|
pub impl Neg for Vec2 with fn neg(a) {
  { x: -a.x, y: -a.y, }
}

///|
/// Component-wise product. For a scalar multiple use `scale`.
pub impl Mul for Vec2 with fn mul(a, b) {
  { x: a.x * b.x, y: a.y * b.y, }
}

///|
pub fn Vec2::scale(self : Vec2, s : Double) -> Vec2 {
  { x: self.x * s, y: self.y * s, }
}

///|
pub fn Vec2::dot(self : Vec2, o : Vec2) -> Double {
  self.x * o.x + self.y * o.y
}

///|
/// z-component of the 3D cross product (signed parallelogram area).
/// Positive when `o` is counter-clockwise from `self` (in a y-up frame).
pub fn Vec2::cross(self : Vec2, o : Vec2) -> Double {
  self.x * o.y - self.y * o.x
}

///|
pub fn Vec2::length_sq(self : Vec2) -> Double {
  self.x * self.x + self.y * self.y
}

///|
pub fn Vec2::length(self : Vec2) -> Double {
  @math.hypot(self.x, self.y)
}

///|
pub fn Vec2::distance(self : Vec2, o : Vec2) -> Double {
  (self - o).length()
}

///|
/// Unit vector, or `None` when the length is zero / not finite.
pub fn Vec2::try_normalize(self : Vec2) -> Vec2? {
  let l = self.length()
  if l > 0.0 && l.is_inf() == false {
    Some(self.scale(1.0 / l))
  } else {
    None
  }
}

///|
/// Unit vector; the zero vector is returned unchanged.
pub fn Vec2::normalize(self : Vec2) -> Vec2 {
  self.try_normalize().unwrap_or(self)
}

///|
/// 90 degrees counter-clockwise (y-up frame).
pub fn Vec2::perp(self : Vec2) -> Vec2 {
  { x: -self.y, y: self.x, }
}

///|
/// Rotate counter-clockwise by `angle` radians.
pub fn Vec2::rotate(self : Vec2, angle : Double) -> Vec2 {
  let c = @math.cos(angle)
  let s = @math.sin(angle)
  { x: self.x * c - self.y * s, y: self.x * s + self.y * c, }
}

///|
pub fn Vec2::lerp(self : Vec2, o : Vec2, t : Double) -> Vec2 {
  { x: self.x + (o.x - self.x) * t, y: self.y + (o.y - self.y) * t, }
}

///|
pub fn Vec2::min(self : Vec2, o : Vec2) -> Vec2 {
  { x: self.x.min(o.x), y: self.y.min(o.y), }
}

///|
pub fn Vec2::max(self : Vec2, o : Vec2) -> Vec2 {
  { x: self.x.max(o.x), y: self.y.max(o.y), }
}

///|
pub fn Vec2::abs(self : Vec2) -> Vec2 {
  { x: self.x.abs(), y: self.y.abs(), }
}

///|
/// Angle of the vector from +X, counter-clockwise, in (-PI, PI].
pub fn Vec2::angle(self : Vec2) -> Double {
  @math.atan2(self.y, self.x)
}

///|
pub fn Vec2::approx_eq(self : Vec2, o : Vec2, eps : Double) -> Bool {
  approx(self.x, o.x, eps) && approx(self.y, o.y, eps)
}

///|
pub fn Vec2::is_finite(self : Vec2) -> Bool {
  !(self.x.is_nan() || self.x.is_inf() || self.y.is_nan() || self.y.is_inf())
}
