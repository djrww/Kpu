// Axis conventions of other tools / engines / file formats and exact conversions to and from
// this library's world (right-handed, +X right, +Y up, -Z forward -- same as OpenGL, glTF,
// three.js, Maya, WebGPU/Vulkan world space).

///|
pub(all) enum Axis {
  PosX
  NegX
  PosY
  NegY
  PosZ
  NegZ
} derive(Eq, Debug)

///|
pub fn Axis::to_vec3(self : Axis) -> Vec3 {
  match self {
    PosX => { x: 1.0, y: 0.0, z: 0.0, }
    NegX => { x: -1.0, y: 0.0, z: 0.0, }
    PosY => { x: 0.0, y: 1.0, z: 0.0, }
    NegY => { x: 0.0, y: -1.0, z: 0.0, }
    PosZ => { x: 0.0, y: 0.0, z: 1.0, }
    NegZ => { x: 0.0, y: 0.0, z: -1.0, }
  }
}

///|
fn Axis::index(self : Axis) -> Int {
  match self {
    PosX | NegX => 0
    PosY | NegY => 1
    PosZ | NegZ => 2
  }
}

///|
/// Which world axes mean "right", "up" and "forward" (the direction a camera / character looks).
#valtype
pub(all) struct Convention {
  right : Axis
  up : Axis
  forward : Axis
} derive(Eq, Debug)

///|
/// This library / OpenGL / glTF / three.js / Maya: right +X, up +Y, forward -Z (right-handed).
pub fn Convention::native() -> Convention {
  { right: PosX, up: PosY, forward: NegZ, }
}

///|
/// Blender, 3ds Max, SketchUp, most CAD: right +X, up +Z, forward +Y (right-handed).
pub fn Convention::z_up_right_handed() -> Convention {
  { right: PosX, up: PosZ, forward: PosY, }
}

///|
/// Unity, Direct3D-style left-handed world: right +X, up +Y, forward +Z.
pub fn Convention::y_up_left_handed() -> Convention {
  { right: PosX, up: PosY, forward: PosZ, }
}

///|
/// Unreal Engine: forward +X, right +Y, up +Z (left-handed).
pub fn Convention::unreal() -> Convention {
  { right: PosY, up: PosZ, forward: PosX, }
}

///|
/// ROS REP-103 / Source engine style: forward +X, left +Y (right -Y), up +Z (right-handed).
pub fn Convention::ros() -> Convention {
  { right: NegY, up: PosZ, forward: PosX, }
}

///|
/// Valid when right/up/forward use three different coordinate axes.
pub fn Convention::is_valid(self : Convention) -> Bool {
  let r = self.right.index()
  let u = self.up.index()
  let f = self.forward.index()
  r != u && r != f && u != f
}

///|
/// Handedness of the world: right-handed when `right x up` points BACKWARD (= -forward).
pub fn Convention::is_right_handed(self : Convention) -> Bool {
  self.right.to_vec3().cross(self.up.to_vec3()).dot(self.forward.to_vec3()) <
  0.0
}

///|
fn Convention::basis(self : Convention) -> Mat3 {
  Mat3::from_cols(
    self.right.to_vec3(),
    self.up.to_vec3(),
    self.forward.to_vec3(),
  )
}

///|
/// Signed permutation matrix M such that `v_to = M * v_from`: right stays right, up stays up,
/// forward stays forward. det(M) = -1 when the two conventions differ in handedness (a mirror);
/// remember to flip triangle winding in that case.
pub fn Convention::basis_change(self : Convention, to : Convention) -> Mat3 {
  if !self.is_valid() || !to.is_valid() {
    abort(
      "Convention::basis_change: invalid convention (axes must be distinct)",
    )
  }
  to.basis() * self.basis().transpose()
}

///|
/// Re-express a point / direction given in `self` convention in `to` convention.
pub fn Convention::convert_vec3(
  self : Convention,
  to : Convention,
  v : Vec3,
) -> Vec3 {
  self.basis_change(to).mul_vec3(v)
}

///|
/// Re-express an orientation: the same physical rotation, described in the other axes.
/// Axis -> M * axis * det(M), angle unchanged (a mirror reverses the sense of rotation).
pub fn Convention::convert_quat(
  self : Convention,
  to : Convention,
  q : Quat,
) -> Quat {
  let m = self.basis_change(to)
  let det = m.determinant()
  let a = m.mul_vec3(Vec3::new(q.x, q.y, q.z)).scale(det)
  Quat::new(a.x, a.y, a.z, q.w)
}

///|
/// Re-express a 4x4 transform: M4 * T * M4^-1 (M4 is its own inverse transpose; orthogonal).
pub fn Convention::convert_mat4(
  self : Convention,
  to : Convention,
  t : Mat4,
) -> Mat4 {
  let m = Mat4::from_mat3(self.basis_change(to))
  m * t * m.transpose()
}
