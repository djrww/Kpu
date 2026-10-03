// Pose / TRS types.

///|
/// Translation * Rotation * Scale. `to_mat4` gives the model matrix.
/// A negative scale component represents mirroring.
#valtype
pub(all) struct Transform {
  translation : Vec3
  rotation : Quat
  scale : Vec3
} derive(Eq, Debug)

///|
pub fn Transform::identity() -> Transform {
  {
    translation: Vec3::zero(),
    rotation: Quat::identity(),
    scale: Vec3::splat(1.0),
  }
}

///|
pub fn Transform::new(t : Vec3, r : Quat, s : Vec3) -> Transform {
  { translation: t, rotation: r, scale: s, }
}

///|
pub fn Transform::to_mat4(self : Transform) -> Mat4 {
  Mat4::from_trs(self.translation, self.rotation, self.scale)
}

///|
/// local -> parent
pub fn Transform::mul_point(self : Transform, p : Vec3) -> Vec3 {
  self.rotation.rotate(p * self.scale) + self.translation
}

///|
pub fn Transform::mul_dir(self : Transform, d : Vec3) -> Vec3 {
  self.rotation.rotate(d * self.scale)
}

///|
/// Rigid body pose: rotation + translation only (no scale). Forms a group: composition and
/// inverse are exact. Use for camera poses, bones, physics bodies, frames.
#valtype
pub(all) struct Rigid {
  rotation : Quat
  translation : Vec3
} derive(Eq, Debug)

///|
pub fn Rigid::identity() -> Rigid {
  { rotation: Quat::identity(), translation: Vec3::zero(), }
}

///|
pub fn Rigid::new(rotation : Quat, translation : Vec3) -> Rigid {
  { rotation, translation, }
}

///|
/// `a * b` = apply `b` first, then `a` (same order as matrices).
pub impl Mul for Rigid with fn mul(a, b) {
  {
    rotation: (a.rotation * b.rotation).normalize(),
    translation: a.rotation.rotate(b.translation) + a.translation,
  }
}

///|
pub fn Rigid::inverse(self : Rigid) -> Rigid {
  let ri = self.rotation.conjugate()
  { rotation: ri, translation: -ri.rotate(self.translation), }
}

///|
pub fn Rigid::mul_point(self : Rigid, p : Vec3) -> Vec3 {
  self.rotation.rotate(p) + self.translation
}

///|
pub fn Rigid::mul_dir(self : Rigid, d : Vec3) -> Vec3 {
  self.rotation.rotate(d)
}

///|
pub fn Rigid::to_mat4(self : Rigid) -> Mat4 {
  Mat4::from_trs(self.translation, self.rotation, Vec3::splat(1.0))
}

///|
/// Extract from a rigid matrix (orthonormal 3x3, affine). Not validated.
pub fn Rigid::from_mat4(m : Mat4) -> Rigid {
  { rotation: Quat::from_mat3(m.to_mat3()), translation: m.translation_part(), }
}

///|
/// Pose (local -> world) of an object at `eye` whose local -Z looks at `target` with +Y
/// toward `up`. For a camera this is the camera-to-world transform.
pub fn Rigid::look_at(eye : Vec3, target : Vec3, up : Vec3) -> Rigid {
  { rotation: Quat::look_rotation(target - eye, up), translation: eye, }
}

///|
pub fn Rigid::lerp(self : Rigid, o : Rigid, t : Double) -> Rigid {
  {
    rotation: self.rotation.slerp(o.rotation, t),
    translation: self.translation.lerp(o.translation, t),
  }
}
