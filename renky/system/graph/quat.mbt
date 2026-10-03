// Unit quaternion (Hamilton convention, x y z = vector part, w = scalar part).
// `a * b` applies `b` first, then `a` (same order as matrices: (a*b).rotate(v) == a.rotate(b.rotate(v))).
// Rotations are right-handed: positive angle about an axis is counter-clockwise looking
// down the axis toward the origin.

///|
#valtype
pub(all) struct Quat {
  x : Double
  y : Double
  z : Double
  w : Double
} derive(Eq, Debug)

///|
pub fn Quat::new(x : Double, y : Double, z : Double, w : Double) -> Quat {
  { x, y, z, w, }
}

///|
pub fn Quat::identity() -> Quat {
  { x: 0.0, y: 0.0, z: 0.0, w: 1.0, }
}

///|
pub impl Mul for Quat with fn mul(a, b) {
  {
    x: a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y,
    y: a.w * b.y - a.x * b.z + a.y * b.w + a.z * b.x,
    z: a.w * b.z + a.x * b.y - a.y * b.x + a.z * b.w,
    w: a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z,
  }
}

///|
pub fn Quat::conjugate(self : Quat) -> Quat {
  { x: -self.x, y: -self.y, z: -self.z, w: self.w, }
}

///|
pub fn Quat::dot(self : Quat, o : Quat) -> Double {
  self.x * o.x + self.y * o.y + self.z * o.z + self.w * o.w
}

///|
pub fn Quat::length(self : Quat) -> Double {
  self.dot(self).sqrt()
}

///|
/// Re-normalise (call after many accumulated multiplications). Zero -> identity.
pub fn Quat::normalize(self : Quat) -> Quat {
  let l = self.length()
  if l > 0.0 && l.is_inf() == false {
    let i = 1.0 / l
    { x: self.x * i, y: self.y * i, z: self.z * i, w: self.w * i, }
  } else {
    Quat::identity()
  }
}

///|
/// Inverse of a *unit* quaternion (= conjugate).
pub fn Quat::inverse(self : Quat) -> Quat {
  self.conjugate()
}

///|
pub fn Quat::negate(self : Quat) -> Quat {
  { x: -self.x, y: -self.y, z: -self.z, w: -self.w, }
}

///|
/// `axis` need not be unit length; a zero axis gives identity.
pub fn Quat::from_axis_angle(axis : Vec3, angle : Double) -> Quat {
  match axis.try_normalize() {
    None => Quat::identity()
    Some(a) => {
      let h = angle * 0.5
      let s = @math.sin(h)
      { x: a.x * s, y: a.y * s, z: a.z * s, w: @math.cos(h), }
    }
  }
}

///|
/// (unit axis, angle in [0, PI]). Identity -> (+X, 0).
pub fn Quat::to_axis_angle(self : Quat) -> (Vec3, Double) {
  let q = if self.w < 0.0 { self.negate() } else { self }
  let sv = Vec3::new(q.x, q.y, q.z)
  let s = sv.length()
  if s < 1.0e-15 {
    (Vec3::unit_x(), 0.0)
  } else {
    (sv.scale(1.0 / s), 2.0 * @math.atan2(s, q.w))
  }
}

///|
pub fn Quat::from_x(angle : Double) -> Quat {
  { x: @math.sin(angle * 0.5), y: 0.0, z: 0.0, w: @math.cos(angle * 0.5), }
}

///|
pub fn Quat::from_y(angle : Double) -> Quat {
  { x: 0.0, y: @math.sin(angle * 0.5), z: 0.0, w: @math.cos(angle * 0.5), }
}

///|
pub fn Quat::from_z(angle : Double) -> Quat {
  { x: 0.0, y: 0.0, z: @math.sin(angle * 0.5), w: @math.cos(angle * 0.5), }
}

///|
/// Rotate a vector (direction or point about the origin).
pub fn Quat::rotate(self : Quat, v : Vec3) -> Vec3 {
  // v' = v + w*t + u x t, with t = 2 (u x v)
  let tx = 2.0 * (self.y * v.z - self.z * v.y)
  let ty = 2.0 * (self.z * v.x - self.x * v.z)
  let tz = 2.0 * (self.x * v.y - self.y * v.x)
  {
    x: v.x + self.w * tx + (self.y * tz - self.z * ty),
    y: v.y + self.w * ty + (self.z * tx - self.x * tz),
    z: v.z + self.w * tz + (self.x * ty - self.y * tx),
  }
}

///|
/// Quaternion from a proper rotation matrix (Shepperd's method, numerically stable).
pub fn Quat::from_mat3(m : Mat3) -> Quat {
  let t = m.m00 + m.m11 + m.m22
  let q : Quat = if t > 0.0 {
    let s = (t + 1.0).sqrt() * 2.0
    {
      x: (m.m21 - m.m12) / s,
      y: (m.m02 - m.m20) / s,
      z: (m.m10 - m.m01) / s,
      w: 0.25 * s,
    }
  } else if m.m00 >= m.m11 && m.m00 >= m.m22 {
    let s = (1.0 + m.m00 - m.m11 - m.m22).sqrt() * 2.0
    {
      x: 0.25 * s,
      y: (m.m01 + m.m10) / s,
      z: (m.m02 + m.m20) / s,
      w: (m.m21 - m.m12) / s,
    }
  } else if m.m11 >= m.m22 {
    let s = (1.0 + m.m11 - m.m00 - m.m22).sqrt() * 2.0
    {
      x: (m.m01 + m.m10) / s,
      y: 0.25 * s,
      z: (m.m12 + m.m21) / s,
      w: (m.m02 - m.m20) / s,
    }
  } else {
    let s = (1.0 + m.m22 - m.m00 - m.m11).sqrt() * 2.0
    {
      x: (m.m02 + m.m20) / s,
      y: (m.m12 + m.m21) / s,
      z: 0.25 * s,
      w: (m.m10 - m.m01) / s,
    }
  }
  q.normalize()
}

///|
pub fn Quat::to_mat3(self : Quat) -> Mat3 {
  Mat3::from_quat(self)
}

///|
/// Shortest-arc rotation taking direction `a` onto direction `b` (inputs need not be unit).
pub fn Quat::from_to(a : Vec3, b : Vec3) -> Quat {
  let ua = a.normalize()
  let ub = b.normalize()
  let d = ua.dot(ub)
  if d < -1.0 + 1.0e-12 {
    // anti-parallel: any axis perpendicular to a
    let p = ua.any_perpendicular()
    { x: p.x, y: p.y, z: p.z, w: 0.0, }
  } else {
    let c = ua.cross(ub)
    Quat::new(c.x, c.y, c.z, 1.0 + d).normalize()
  }
}

///|
/// Orientation whose local -Z axis points along `forward` and whose local +Y is as close
/// as possible to `up` (the camera / OpenGL "look" convention). Degenerate `up` is handled.
pub fn Quat::look_rotation(forward : Vec3, up : Vec3) -> Quat {
  Quat::from_mat3(Mat3::look_rotation(forward, up))
}

///|
/// exp-map: rotation vector (axis * angle) -> quaternion. Stable for tiny angles.
pub fn Quat::from_rotation_vector(v : Vec3) -> Quat {
  let th = v.length()
  if th < 1.0e-8 {
    // sin(th/2)/th ~ 1/2 - th^2/48
    let k = 0.5 - th * th / 48.0
    Quat::new(v.x * k, v.y * k, v.z * k, 1.0 - th * th / 8.0).normalize()
  } else {
    let k = @math.sin(th * 0.5) / th
    { x: v.x * k, y: v.y * k, z: v.z * k, w: @math.cos(th * 0.5), }
  }
}

///|
/// log-map: quaternion -> rotation vector with angle in [0, PI] (shortest representation).
pub fn Quat::to_rotation_vector(self : Quat) -> Vec3 {
  let q = if self.w < 0.0 { self.negate() } else { self }
  let sv = Vec3::new(q.x, q.y, q.z)
  let s = sv.length()
  if s < 1.0e-12 {
    sv.scale(2.0 / q.w)
  } else {
    sv.scale(2.0 * @math.atan2(s, q.w) / s)
  }
}

///|
/// Advance orientation by a WORLD-space angular velocity `omega` (rad/s) over `dt` seconds.
pub fn Quat::integrate(self : Quat, omega : Vec3, dt : Double) -> Quat {
  (Quat::from_rotation_vector(omega.scale(dt)) * self).normalize()
}

///|
/// Normalised linear interpolation along the shortest arc.
pub fn Quat::nlerp(self : Quat, o : Quat, t : Double) -> Quat {
  let b = if self.dot(o) < 0.0 { o.negate() } else { o }
  Quat::new(
    self.x + (b.x - self.x) * t,
    self.y + (b.y - self.y) * t,
    self.z + (b.z - self.z) * t,
    self.w + (b.w - self.w) * t,
  ).normalize()
}

///|
/// Spherical interpolation along the shortest arc (constant angular speed).
pub fn Quat::slerp(self : Quat, o : Quat, t : Double) -> Quat {
  let mut d = self.dot(o)
  let b = if d < 0.0 {
    d = -d
    o.negate()
  } else {
    o
  }
  if d > 0.9995 {
    return self.nlerp(b, t)
  }
  let th = @math.acos(d)
  let s = @math.sin(th)
  let wa = @math.sin((1.0 - t) * th) / s
  let wb = @math.sin(t * th) / s
  {
    x: self.x * wa + b.x * wb,
    y: self.y * wa + b.y * wb,
    z: self.z * wa + b.z * wb,
    w: self.w * wa + b.w * wb,
  }
}

///|
/// Smallest rotation angle (radians, [0, PI]) between two orientations.
pub fn Quat::angle_to(self : Quat, o : Quat) -> Double {
  let d = self.dot(o).abs()
  2.0 * safe_acos(d)
}

///|
/// True when both quaternions represent the same rotation (q and -q are equal).
pub fn Quat::approx_eq(self : Quat, o : Quat, eps : Double) -> Bool {
  let s = if self.dot(o) < 0.0 { -1.0 } else { 1.0 }
  approx(self.x, o.x * s, eps) &&
  approx(self.y, o.y * s, eps) &&
  approx(self.z, o.z * s, eps) &&
  approx(self.w, o.w * s, eps)
}
