// Projection: view space -> clip space -> NDC.
//
// CONVENTIONS (fixed for the whole library):
//   * World & view space: RIGHT-HANDED, +X right, +Y up, camera looks down -Z.
//   * Clip space: homogeneous (x, y, z, w). Perspective gives w = -z_view.
//   * NDC: x in [-1, 1] (left->right), y in [-1, 1] (bottom->top), z by `DepthMode`.

///|
/// How NDC depth is laid out.
///  - `ZeroToOne`:   near -> 0, far -> 1   (Vulkan / D3D / Metal / WebGPU native)
///  - `NegOneToOne`: near -> -1, far -> 1  (OpenGL default)
///  - `ReversedZ`:   near -> 1, far -> 0   (best float precision; use with a GREATER depth test
///    and clear depth to 0). With an infinite far plane this is the recommended modern setup.
pub(all) enum DepthMode {
  ZeroToOne
  NegOneToOne
  ReversedZ
} derive(Eq, Debug)

///|
pub(all) enum ProjKind {
  Perspective
  Orthographic
} derive(Eq, Debug)

///|
/// A projection volume. `left/right/bottom/top` are the extents in VIEW-space units measured
/// at the NEAR plane (perspective) or of the box (orthographic). `far` may be
/// `@double.infinity` for perspective projections.
#valtype
pub(all) struct Projection {
  kind : ProjKind
  left : Double
  right : Double
  bottom : Double
  top : Double
  near : Double
  far : Double
  depth : DepthMode
} derive(Eq, Debug)

///|
/// Symmetric perspective from vertical field of view (radians) and width/height aspect.
/// `far` may be `@double.infinity`.
pub fn Projection::perspective(
  fov_y : Double,
  aspect : Double,
  near : Double,
  far : Double,
  depth : DepthMode,
) -> Projection {
  if !(fov_y > 0.0 && fov_y < PI) {
    abort("Projection::perspective: fov_y must be in (0, PI)")
  }
  if !(aspect > 0.0) {
    abort("Projection::perspective: aspect must be > 0")
  }
  let top = near * @math.tan(fov_y * 0.5)
  Projection::frustum(-top * aspect, top * aspect, -top, top, near, far, depth)
}

///|
/// Off-centre perspective frustum (glFrustum): extents on the near plane. Used for
/// asymmetric FOV (VR eyes), tiled/sub-rect rendering, lens shift, portals.
pub fn Projection::frustum(
  left : Double,
  right : Double,
  bottom : Double,
  top : Double,
  near : Double,
  far : Double,
  depth : DepthMode,
) -> Projection {
  if !(near > 0.0) {
    abort("Projection::frustum: near must be > 0")
  }
  if !(far > near) {
    abort("Projection::frustum: far must be > near")
  }
  if !(right > left && top > bottom) {
    abort("Projection::frustum: need right > left and top > bottom")
  }
  { kind: Perspective, left, right, bottom, top, near, far, depth, }
}

///|
/// Asymmetric perspective from four half-angles (radians, all measured outward from the view
/// axis and normally positive): the VR/OpenXR "field of view" form.
pub fn Projection::perspective_fov(
  angle_left : Double,
  angle_right : Double,
  angle_down : Double,
  angle_up : Double,
  near : Double,
  far : Double,
  depth : DepthMode,
) -> Projection {
  Projection::frustum(
    -near * @math.tan(angle_left),
    near * @math.tan(angle_right),
    -near * @math.tan(angle_down),
    near * @math.tan(angle_up),
    near,
    far,
    depth,
  )
}

///|
pub fn Projection::orthographic(
  left : Double,
  right : Double,
  bottom : Double,
  top : Double,
  near : Double,
  far : Double,
  depth : DepthMode,
) -> Projection {
  if !(far > near) {
    abort("Projection::orthographic: far must be > near")
  }
  if far.is_inf() {
    abort("Projection::orthographic: far must be finite")
  }
  if !(right > left && top > bottom) {
    abort("Projection::orthographic: need right > left and top > bottom")
  }
  { kind: Orthographic, left, right, bottom, top, near, far, depth, }
}

///|
/// Centered orthographic box of the given width / height.
pub fn Projection::orthographic_centered(
  width : Double,
  height : Double,
  near : Double,
  far : Double,
  depth : DepthMode,
) -> Projection {
  Projection::orthographic(
    -width * 0.5,
    width * 0.5,
    -height * 0.5,
    height * 0.5,
    near,
    far,
    depth,
  )
}

///|
pub fn Projection::is_perspective(self : Projection) -> Bool {
  self.kind == Perspective
}

///|
pub fn Projection::is_infinite_far(self : Projection) -> Bool {
  self.far.is_inf()
}

///|
/// (A, B) such that  z_clip = A * z_view + B  (perspective: w = -z_view, so
/// z_ndc = (A z + B) / -z;  orthographic: w = 1, z_ndc = A z + B).
fn Projection::depth_coeffs(self : Projection) -> (Double, Double) {
  let n = self.near
  let f = self.far
  match self.kind {
    Perspective =>
      if f.is_inf() {
        match self.depth {
          ZeroToOne => (-1.0, -n)
          NegOneToOne => (-1.0, -2.0 * n)
          ReversedZ => (0.0, n)
        }
      } else {
        let d = f - n
        match self.depth {
          ZeroToOne => (-f / d, -f * n / d)
          NegOneToOne => (-(f + n) / d, -2.0 * f * n / d)
          ReversedZ => (n / d, n * f / d)
        }
      }
    Orthographic => {
      let d = f - n
      match self.depth {
        ZeroToOne => (-1.0 / d, -n / d)
        NegOneToOne => (-2.0 / d, -(f + n) / d)
        ReversedZ => (1.0 / d, f / d)
      }
    }
  }
}

///|
/// The 4x4 projection matrix (column vectors, v_clip = P * v_view).
pub fn Projection::to_mat4(self : Projection) -> Mat4 {
  let (a, b) = self.depth_coeffs()
  let rl = self.right - self.left
  let tb = self.top - self.bottom
  match self.kind {
    Perspective => {
      let n2 = 2.0 * self.near
      {
        m00: n2 / rl,
        m01: 0.0,
        m02: (self.right + self.left) / rl,
        m03: 0.0,
        m10: 0.0,
        m11: n2 / tb,
        m12: (self.top + self.bottom) / tb,
        m13: 0.0,
        m20: 0.0,
        m21: 0.0,
        m22: a,
        m23: b,
        m30: 0.0,
        m31: 0.0,
        m32: -1.0,
        m33: 0.0,
      }
    }
    Orthographic =>
      {
        m00: 2.0 / rl,
        m01: 0.0,
        m02: 0.0,
        m03: -(self.right + self.left) / rl,
        m10: 0.0,
        m11: 2.0 / tb,
        m12: 0.0,
        m13: -(self.top + self.bottom) / tb,
        m20: 0.0,
        m21: 0.0,
        m22: a,
        m23: b,
        m30: 0.0,
        m31: 0.0,
        m32: 0.0,
        m33: 1.0,
      }
  }
}

///|
/// Closed-form inverse (clip -> view), much cheaper and more accurate than a general
/// 4x4 inverse. The result of `inverse * clip` has a w component that must be divided out
/// for perspective (use `Mat4::project_point`).
pub fn Projection::inverse_mat4(self : Projection) -> Mat4 {
  let (a, b) = self.depth_coeffs()
  let rl = self.right - self.left
  let tb = self.top - self.bottom
  match self.kind {
    Perspective => {
      let m00 = 2.0 * self.near / rl
      let m02 = (self.right + self.left) / rl
      let m11 = 2.0 * self.near / tb
      let m12 = (self.top + self.bottom) / tb
      {
        m00: 1.0 / m00,
        m01: 0.0,
        m02: 0.0,
        m03: m02 / m00,
        m10: 0.0,
        m11: 1.0 / m11,
        m12: 0.0,
        m13: m12 / m11,
        m20: 0.0,
        m21: 0.0,
        m22: 0.0,
        m23: -1.0,
        m30: 0.0,
        m31: 0.0,
        m32: 1.0 / b,
        m33: a / b,
      }
    }
    Orthographic => {
      let m00 = 2.0 / rl
      let m03 = -(self.right + self.left) / rl
      let m11 = 2.0 / tb
      let m13 = -(self.top + self.bottom) / tb
      {
        m00: 1.0 / m00,
        m01: 0.0,
        m02: 0.0,
        m03: -m03 / m00,
        m10: 0.0,
        m11: 1.0 / m11,
        m12: 0.0,
        m13: -m13 / m11,
        m20: 0.0,
        m21: 0.0,
        m22: 1.0 / a,
        m23: -b / a,
        m30: 0.0,
        m31: 0.0,
        m32: 0.0,
        m33: 1.0,
      }
    }
  }
}

///|
/// NDC depth of a view-space z (negative in front of the camera).
pub fn Projection::ndc_depth(self : Projection, view_z : Double) -> Double {
  let (a, b) = self.depth_coeffs()
  match self.kind {
    Perspective => (a * view_z + b) / -view_z
    Orthographic => a * view_z + b
  }
}

///|
/// View-space z (negative in front) from an NDC depth value.
pub fn Projection::view_z_from_ndc(self : Projection, ndc_z : Double) -> Double {
  let (a, b) = self.depth_coeffs()
  match self.kind {
    Perspective => -b / (ndc_z + a)
    Orthographic => (ndc_z - b) / a
  }
}

///|
/// Distance in front of the camera plane (positive), from an NDC depth value. This is the
/// "linearized depth" used for fog, SSAO, soft particles, etc.
pub fn Projection::linear_depth(self : Projection, ndc_z : Double) -> Double {
  -self.view_z_from_ndc(ndc_z)
}

///|
/// View-space point -> NDC. `None` for perspective points on/behind the eye plane
/// (z_view >= 0), where the projection is not defined for rendering.
pub fn Projection::project_view(self : Projection, p : Vec3) -> Vec3? {
  let (a, b) = self.depth_coeffs()
  let rl = self.right - self.left
  let tb = self.top - self.bottom
  match self.kind {
    Perspective => {
      let w = -p.z
      if !(w > 0.0) {
        return None
      }
      let iw = 1.0 / w
      let n2 = 2.0 * self.near
      Some({
        x: (n2 / rl * p.x + (self.right + self.left) / rl * p.z) * iw,
        y: (n2 / tb * p.y + (self.top + self.bottom) / tb * p.z) * iw,
        z: (a * p.z + b) * iw,
      })
    }
    Orthographic =>
      Some({
        x: (2.0 * p.x - (self.right + self.left)) / rl,
        y: (2.0 * p.y - (self.top + self.bottom)) / tb,
        z: a * p.z + b,
      })
  }
}

///|
/// NDC -> view-space point (exact inverse of `project_view`).
pub fn Projection::unproject_ndc(self : Projection, ndc : Vec3) -> Vec3 {
  let zv = self.view_z_from_ndc(ndc.z)
  let rl = self.right - self.left
  let tb = self.top - self.bottom
  match self.kind {
    Perspective => {
      let m00 = 2.0 * self.near / rl
      let m02 = (self.right + self.left) / rl
      let m11 = 2.0 * self.near / tb
      let m12 = (self.top + self.bottom) / tb
      { x: -zv * (ndc.x + m02) / m00, y: -zv * (ndc.y + m12) / m11, z: zv, }
    }
    Orthographic =>
      {
        x: (ndc.x * rl + (self.right + self.left)) * 0.5,
        y: (ndc.y * tb + (self.top + self.bottom)) * 0.5,
        z: zv,
      }
  }
}

///|
/// View-space ray through an NDC xy position: (origin, unit direction).
///  - perspective: origin = eye (0,0,0)
///  - orthographic: origin lies on the near plane, direction = -Z
pub fn Projection::view_ray(
  self : Projection,
  ndc_x : Double,
  ndc_y : Double,
) -> (Vec3, Vec3) {
  let rl = self.right - self.left
  let tb = self.top - self.bottom
  match self.kind {
    Perspective => {
      let m00 = 2.0 * self.near / rl
      let m02 = (self.right + self.left) / rl
      let m11 = 2.0 * self.near / tb
      let m12 = (self.top + self.bottom) / tb
      let d = Vec3::new((ndc_x + m02) / m00, (ndc_y + m12) / m11, -1.0)
      (Vec3::zero(), d.normalize())
    }
    Orthographic =>
      (
        Vec3::new(
          (ndc_x * rl + (self.right + self.left)) * 0.5,
          (ndc_y * tb + (self.top + self.bottom)) * 0.5,
          -self.near,
        ),
        Vec3::new(0.0, 0.0, -1.0),
      )
  }
}

///|
/// Vertical field of view (radians) of a perspective projection (0 for orthographic).
pub fn Projection::fov_y(self : Projection) -> Double {
  match self.kind {
    Perspective =>
      @math.atan(self.top / self.near) + @math.atan(-self.bottom / self.near)
    Orthographic => 0.0
  }
}

///|
/// Width / height of the projection window.
pub fn Projection::aspect(self : Projection) -> Double {
  (self.right - self.left) / (self.top - self.bottom)
}

///|
/// Same vertical extent, new aspect ratio (keeps relative horizontal lens shift).
pub fn Projection::with_aspect(
  self : Projection,
  aspect : Double,
) -> Projection {
  let h = self.top - self.bottom
  let k = h * aspect / (self.right - self.left)
  { ..self, left: self.left * k, right: self.right * k, }
}

///|
pub fn Projection::with_depth(
  self : Projection,
  depth : DepthMode,
) -> Projection {
  { ..self, depth, }
}

///|
/// Convert an NDC z in any mode to the [0, 1] "window" range (0 = viewport.min_depth).
pub fn DepthMode::to_unit(self : DepthMode, ndc_z : Double) -> Double {
  match self {
    NegOneToOne => ndc_z * 0.5 + 0.5
    _ => ndc_z
  }
}

///|
/// Inverse of `to_unit`.
pub fn DepthMode::from_unit(self : DepthMode, unit_z : Double) -> Double {
  match self {
    NegOneToOne => unit_z * 2.0 - 1.0
    _ => unit_z
  }
}
