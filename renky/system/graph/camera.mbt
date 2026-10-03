// Camera: ties world, view, clip, NDC and window coordinates together.
//
//   local --Model--> world --View--> view --Projection--> clip --/w--> NDC --Viewport--> window

///|
pub(all) struct Camera {
  /// Camera -> world transform (position + orientation). Looks down its local -Z.
  pose : Rigid
  projection : Projection
  viewport : Viewport
} derive(Eq, Debug)

///|
pub fn Camera::new(
  pose : Rigid,
  projection : Projection,
  viewport : Viewport,
) -> Camera {
  { pose, projection, viewport, }
}

///|
pub fn Camera::look_at(
  eye : Vec3,
  target : Vec3,
  up : Vec3,
  projection : Projection,
  viewport : Viewport,
) -> Camera {
  { pose: Rigid::look_at(eye, target, up), projection, viewport, }
}

///|
/// Orbit / FPS style: `yaw`, `pitch` give the viewing direction as a `Spherical` direction
/// (yaw 0, pitch 0 looks down -Z; +yaw turns left about +Y; +pitch looks up). The eye sits
/// `distance` behind `target` along that direction.
pub fn Camera::orbit(
  target : Vec3,
  distance : Double,
  yaw : Double,
  pitch : Double,
  projection : Projection,
  viewport : Viewport,
) -> Camera {
  let f = Spherical::new(1.0, yaw, pitch).to_cartesian()
  let eye = target - f.scale(distance)
  { pose: Rigid::look_at(eye, target, Vec3::unit_y()), projection, viewport, }
}

///|
pub fn Camera::eye(self : Camera) -> Vec3 {
  self.pose.translation
}

///|
pub fn Camera::forward(self : Camera) -> Vec3 {
  self.pose.rotation.rotate(Vec3::new(0.0, 0.0, -1.0))
}

///|
pub fn Camera::right(self : Camera) -> Vec3 {
  self.pose.rotation.rotate(Vec3::unit_x())
}

///|
pub fn Camera::up(self : Camera) -> Vec3 {
  self.pose.rotation.rotate(Vec3::unit_y())
}

///|
/// World -> view matrix.
pub fn Camera::view_matrix(self : Camera) -> Mat4 {
  self.pose.to_mat4().rigid_inverse()
}

///|
pub fn Camera::proj_matrix(self : Camera) -> Mat4 {
  self.projection.to_mat4()
}

///|
/// World -> clip matrix.
pub fn Camera::view_proj(self : Camera) -> Mat4 {
  self.projection.to_mat4() * self.view_matrix()
}

///|
/// World -> clip with the eye moved to the origin ("camera-relative / floating origin"
/// rendering). Pair with model matrices whose translation has `eye` subtracted in DOUBLE
/// precision (`Mat4::relative_to`) before the single rounding to float: this removes the
/// jitter of huge world coordinates.
pub fn Camera::view_proj_rte(self : Camera) -> Mat4 {
  let inv_rot = Mat4::from_mat3(Mat3::from_quat(self.pose.rotation).transpose())
  self.projection.to_mat4() * inv_rot
}

///|
pub fn Camera::world_to_view(self : Camera, p : Vec3) -> Vec3 {
  self.pose.inverse().mul_point(p)
}

///|
pub fn Camera::view_to_world(self : Camera, p : Vec3) -> Vec3 {
  self.pose.mul_point(p)
}

///|
pub fn Camera::world_to_clip(self : Camera, p : Vec3) -> Vec4 {
  self.view_proj().mul_vec4(p.to_vec4(1.0))
}

///|
/// World -> NDC. `None` if behind the eye plane (perspective).
pub fn Camera::world_to_ndc(self : Camera, p : Vec3) -> Vec3? {
  self.projection.project_view(self.world_to_view(p))
}

///|
/// World -> window: (x px, y px, window depth). `None` if behind the camera.
pub fn Camera::world_to_screen(self : Camera, p : Vec3) -> Vec3? {
  match self.world_to_ndc(p) {
    Some(n) => Some(self.viewport.ndc_to_screen(n, self.projection.depth))
    None => None
  }
}

///|
/// Window position + window depth (as read back from the depth buffer) -> world.
pub fn Camera::screen_to_world(
  self : Camera,
  sx : Double,
  sy : Double,
  depth : Double,
) -> Vec3 {
  let ndc = self.viewport.screen_to_ndc(sx, sy, depth, self.projection.depth)
  self.pose.mul_point(self.projection.unproject_ndc(ndc))
}

///|
/// World point under window position (sx, sy) at `distance` in front of the camera plane.
pub fn Camera::screen_to_world_at_distance(
  self : Camera,
  sx : Double,
  sy : Double,
  distance : Double,
) -> Vec3 {
  let ray = self.pick_ray(sx, sy)
  match self.projection.kind {
    Perspective => {
      // point on the ray whose depth along the camera axis equals `distance`
      let cosang = ray.dir.dot(self.forward())
      ray.at(distance / cosang)
    }
    Orthographic =>
      // origin is on the near plane (depth = near)
      ray.at(distance - self.projection.near)
  }
}

///|
/// World-space picking ray through a window position (unit direction).
pub fn Camera::pick_ray(self : Camera, sx : Double, sy : Double) -> Ray {
  let ndc = self.viewport.screen_to_ndc_xy({ x: sx, y: sy, })
  let (o, d) = self.projection.view_ray(ndc.x, ndc.y)
  { origin: self.pose.mul_point(o), dir: self.pose.mul_dir(d), }
}

///|
pub fn Camera::frustum(self : Camera) -> Frustum {
  Frustum::from_view_proj(self.view_proj(), self.projection.depth)
}

///|
/// Change the viewport and keep the picture undistorted by updating the projection aspect.
pub fn Camera::resize(self : Camera, width : Double, height : Double) -> Camera {
  let vp = Viewport::new(self.viewport.x, self.viewport.y, width, height)
  {
    ..self,
    viewport: vp,
    projection: self.projection.with_aspect(width / height),
  }
}

///|
/// Model matrix expressed relative to `origin` (translation reduced in double precision).
/// Equivalent to `translation(-origin) * m` but without a full matrix multiply.
pub fn Mat4::relative_to(self : Mat4, origin : Vec3) -> Mat4 {
  {
    ..self,
    m03: self.m03 - origin.x,
    m13: self.m13 - origin.y,
    m23: self.m23 - origin.z,
  }
}
