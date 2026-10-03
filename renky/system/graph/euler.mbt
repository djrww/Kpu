// Euler angles for all 12 axis sequences.
//
// Convention (INTRINSIC, "rotating axes", same as three.js / Eigen / scipy upper-case):
//   order ABC with angles (a, b, c)  =>  R = R_A(a) * R_B(b) * R_C(c)
// i.e. rotate about A by `a`, then about the *new* B axis by `b`, then about the *new* C
// axis by `c`. The same rotation is the EXTRINSIC (fixed axes) sequence CBA with angles
// (c, b, a).
//
// Tait-Bryan (3 distinct axes): XYZ XZY YXZ YZX ZXY ZYX   (yaw/pitch/roll families)
// Proper Euler (first = last):  XYX XZX YXY YZY ZXZ ZYZ

///|
pub(all) enum EulerOrder {
  XYZ
  XZY
  YXZ
  YZX
  ZXY
  ZYX
  XYX
  XZX
  YXY
  YZY
  ZXZ
  ZYZ
} derive(Eq, Debug)

///|
/// Axis indices (0=X, 1=Y, 2=Z) of the three rotations in order.
pub fn EulerOrder::axes(self : EulerOrder) -> (Int, Int, Int) {
  match self {
    XYZ => (0, 1, 2)
    XZY => (0, 2, 1)
    YXZ => (1, 0, 2)
    YZX => (1, 2, 0)
    ZXY => (2, 0, 1)
    ZYX => (2, 1, 0)
    XYX => (0, 1, 0)
    XZX => (0, 2, 0)
    YXY => (1, 0, 1)
    YZY => (1, 2, 1)
    ZXZ => (2, 0, 2)
    ZYZ => (2, 1, 2)
  }
}

///|
pub fn EulerOrder::is_proper(self : EulerOrder) -> Bool {
  let (a, _, c) = self.axes()
  a == c
}

///|
fn axis_quat(axis : Int, angle : Double) -> Quat {
  match axis {
    0 => Quat::from_x(angle)
    1 => Quat::from_y(angle)
    _ => Quat::from_z(angle)
  }
}

///|
/// Quaternion for R_A(a) * R_B(b) * R_C(c) with (A,B,C) = `order`.
pub fn Quat::from_euler(
  order : EulerOrder,
  a : Double,
  b : Double,
  c : Double,
) -> Quat {
  let (i, j, k) = order.axes()
  (axis_quat(i, a) * axis_quat(j, b) * axis_quat(k, c)).normalize()
}

///|
/// Matrix for R_A(a) * R_B(b) * R_C(c).
pub fn Mat3::from_euler(
  order : EulerOrder,
  a : Double,
  b : Double,
  c : Double,
) -> Mat3 {
  Mat3::from_quat(Quat::from_euler(order, a, b, c))
}

///|
fn m3_get(m : Mat3, r : Int, c : Int) -> Double {
  match (r, c) {
    (0, 0) => m.m00
    (0, 1) => m.m01
    (0, 2) => m.m02
    (1, 0) => m.m10
    (1, 1) => m.m11
    (1, 2) => m.m12
    (2, 0) => m.m20
    (2, 1) => m.m21
    _ => m.m22
  }
}

///|
/// Recover (a, b, c) such that `Mat3::from_euler(order, a, b, c)` reproduces `m`.
/// Ranges: Tait-Bryan: a, c in [-PI, PI], b in [-PI/2, PI/2]; proper: a, c in [-PI, PI],
/// b in [0, PI] (a ~ in (-PI, PI] up to the sign convention). At gimbal lock the
/// decomposition is not unique; one valid solution is returned (the matrix still round-trips).
pub fn Mat3::to_euler(
  self : Mat3,
  order : EulerOrder,
) -> (Double, Double, Double) {
  let (a0, a1, a2) = order.axes()
  let odd = if (a0 + 1) % 3 == a1 { 0 } else { 1 }
  let i = a0
  let j = (a0 + 1 + odd) % 3
  let k = (a0 + 2 - odd) % 3
  let mut r0 = 0.0
  let mut r1 = 0.0
  let mut r2 = 0.0
  if a0 == a2 {
    r0 = @math.atan2(m3_get(self, j, i), m3_get(self, k, i))
    let s2 = @math.hypot(m3_get(self, j, i), m3_get(self, k, i))
    if (odd == 1 && r0 < 0.0) || (odd == 0 && r0 > 0.0) {
      r0 = if r0 > 0.0 { r0 - PI } else { r0 + PI }
      r1 = -@math.atan2(s2, m3_get(self, i, i))
    } else {
      r1 = @math.atan2(s2, m3_get(self, i, i))
    }
    let s1 = @math.sin(r0)
    let c1 = @math.cos(r0)
    r2 = @math.atan2(
      c1 * m3_get(self, j, k) - s1 * m3_get(self, k, k),
      c1 * m3_get(self, j, j) - s1 * m3_get(self, k, j),
    )
  } else {
    r0 = @math.atan2(m3_get(self, j, k), m3_get(self, k, k))
    let c2 = @math.hypot(m3_get(self, i, i), m3_get(self, i, j))
    if (odd == 1 && r0 < 0.0) || (odd == 0 && r0 > 0.0) {
      r0 = if r0 > 0.0 { r0 - PI } else { r0 + PI }
      r1 = @math.atan2(-m3_get(self, i, k), -c2)
    } else {
      r1 = @math.atan2(-m3_get(self, i, k), c2)
    }
    let s1 = @math.sin(r0)
    let c1 = @math.cos(r0)
    r2 = @math.atan2(
      s1 * m3_get(self, k, i) - c1 * m3_get(self, j, i),
      c1 * m3_get(self, j, j) - s1 * m3_get(self, k, j),
    )
  }
  if odd == 0 {
    r0 = -r0
    r1 = -r1
    r2 = -r2
  }
  // Canonical branch: Tait-Bryan b in [-PI/2, PI/2]; proper b in [0, PI]; a, c in (-PI, PI].
  if a0 == a2 {
    if r1 < 0.0 {
      r0 = r0 + PI
      r1 = -r1
      r2 = r2 - PI
    }
  } else if r1 > HALF_PI || r1 < -HALF_PI {
    r0 = r0 + PI
    r1 = PI - r1
    r2 = r2 + PI
  }
  (wrap_pi(r0), wrap_pi(r1), wrap_pi(r2))
}

///|
pub fn Quat::to_euler(
  self : Quat,
  order : EulerOrder,
) -> (Double, Double, Double) {
  Mat3::from_quat(self).to_euler(order)
}
