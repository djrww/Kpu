// Per-surface and per-primitive coordinate systems: orthonormal frames, tangent space
// (normal mapping), barycentric coordinates (rasterisation / interpolation).

///|
/// Orthonormal frame: tangent (T), bitangent (B), normal (N). Right-handed when
/// cross(T, B) == N, left-handed (mirrored UVs) when cross(T, B) == -N.
#valtype
pub(all) struct Frame {
  tangent : Vec3
  bitangent : Vec3
  normal : Vec3
} derive(Eq, Debug)

///|
/// Branchless orthonormal basis from a UNIT vector (Duff et al. 2017, "Building an
/// Orthonormal Basis, Revisited"): continuous except at n = (0, 0, -1) and free of
/// singularities for every other direction.
pub fn orthonormal_basis(n : Vec3) -> (Vec3, Vec3) {
  let s = sign_nz(n.z)
  let a = -1.0 / (s + n.z)
  let b = n.x * n.y * a
  (
    { x: 1.0 + s * n.x * n.x * a, y: s * b, z: -s * n.x, },
    { x: b, y: s + n.y * n.y * a, z: -n.y, },
  )
}

///|
/// Right-handed frame whose normal is `n` (normalised internally); T, B arbitrary but stable.
pub fn Frame::from_normal(n : Vec3) -> Frame {
  let nn = n.normalize()
  let (t, b) = orthonormal_basis(nn)
  { tangent: t, bitangent: b, normal: nn, }
}

///|
/// Local (T, B, N) coordinates -> parent space.
pub fn Frame::to_world(self : Frame, loc : Vec3) -> Vec3 {
  self.tangent.scale(loc.x) +
  self.bitangent.scale(loc.y) +
  self.normal.scale(loc.z)
}

///|
/// Parent space -> local (T, B, N) coordinates (the frame is orthonormal: transpose = inverse).
pub fn Frame::to_local(self : Frame, v : Vec3) -> Vec3 {
  { x: v.dot(self.tangent), y: v.dot(self.bitangent), z: v.dot(self.normal), }
}

///|
/// Matrix with columns T, B, N (local -> parent).
pub fn Frame::to_mat3(self : Frame) -> Mat3 {
  Mat3::from_cols(self.tangent, self.bitangent, self.normal)
}

///|
/// +1 for a right-handed frame, -1 for mirrored (the glTF `tangent.w` value).
pub fn Frame::handedness(self : Frame) -> Double {
  sign_nz(self.tangent.cross(self.bitangent).dot(self.normal))
}

///|
/// Tangent frame of a triangle from its positions and UVs, for normal mapping.
/// T points along increasing u (dP/du), B along increasing v (dP/dv, which is DOWN in the
/// top-left UV convention -- negate B, or `flip_v` the UVs, for a Y-up normal map).
/// T and B are orthonormalised against the supplied vertex normal `n`; the handedness of
/// the UV mapping is preserved in B. `None` for degenerate (zero-area UV) triangles.
pub fn tangent_frame(
  p0 : Vec3,
  p1 : Vec3,
  p2 : Vec3,
  uv0 : Vec2,
  uv1 : Vec2,
  uv2 : Vec2,
  n : Vec3,
) -> Frame? {
  let e1 = p1 - p0
  let e2 = p2 - p0
  let d1 = uv1 - uv0
  let d2 = uv2 - uv0
  let det = d1.x * d2.y - d2.x * d1.y
  if det == 0.0 {
    return None
  }
  let r = 1.0 / det
  let t_raw = (e1.scale(d2.y) - e2.scale(d1.y)).scale(r)
  let b_raw = (e2.scale(d1.x) - e1.scale(d2.x)).scale(r)
  let nn = n.normalize()
  let t = (t_raw - nn.scale(nn.dot(t_raw))).try_normalize()
  match t {
    None => None
    Some(t) => {
      let nxt = nn.cross(t)
      let h = sign_nz(nxt.dot(b_raw))
      Some({ tangent: t, bitangent: nxt.scale(h), normal: nn, })
    }
  }
}

///|
/// Twice the signed area of triangle (a, b, c): positive when counter-clockwise in a Y-UP frame
/// (NDC). In window space (Y down) the same triangle is negative.
pub fn triangle_area2(a : Vec2, b : Vec2, c : Vec2) -> Double {
  (b - a).cross(c - a)
}

///|
/// Edge function: > 0 when p is to the LEFT of the directed edge a -> b (Y-up frame).
pub fn edge_function(a : Vec2, b : Vec2, p : Vec2) -> Double {
  (b - a).cross(p - a)
}

///|
/// Barycentric coordinates (wa, wb, wc) of p in triangle (a, b, c); sum is 1; all >= 0 iff
/// p is inside. `None` for a degenerate triangle.
pub fn barycentric_2d(p : Vec2, a : Vec2, b : Vec2, c : Vec2) -> Vec3? {
  let area = triangle_area2(a, b, c)
  if area == 0.0 {
    return None
  }
  let inv = 1.0 / area
  let wa = edge_function(b, c, p) * inv
  let wb = edge_function(c, a, p) * inv
  Some({ x: wa, y: wb, z: 1.0 - wa - wb, })
}

///|
/// Barycentric coordinates of the projection of `p` onto the plane of triangle (a, b, c).
pub fn barycentric_3d(p : Vec3, a : Vec3, b : Vec3, c : Vec3) -> Vec3? {
  let v0 = b - a
  let v1 = c - a
  let v2 = p - a
  let d00 = v0.dot(v0)
  let d01 = v0.dot(v1)
  let d11 = v1.dot(v1)
  let d20 = v2.dot(v0)
  let d21 = v2.dot(v1)
  let den = d00 * d11 - d01 * d01
  if den == 0.0 {
    return None
  }
  let v = (d11 * d20 - d01 * d21) / den
  let w = (d00 * d21 - d01 * d20) / den
  Some({ x: 1.0 - v - w, y: v, z: w, })
}

///|
/// True when all three barycentric weights are >= -eps.
pub fn bary_inside(b : Vec3, eps : Double) -> Bool {
  b.x >= -eps && b.y >= -eps && b.z >= -eps
}

///|
/// Perspective-correct barycentrics. `screen_bary` are the weights measured in 2D window/NDC
/// space; `w0, w1, w2` are the clip-space w of the three vertices (= view-space depth for a
/// perspective projection). The result are the weights in the 3D triangle, to be used to
/// interpolate any attribute (UV, normal, colour, view-space position).
pub fn perspective_correct(
  screen_bary : Vec3,
  w0 : Double,
  w1 : Double,
  w2 : Double,
) -> Vec3 {
  let a = screen_bary.x / w0
  let b = screen_bary.y / w1
  let c = screen_bary.z / w2
  let s = 1.0 / (a + b + c)
  { x: a * s, y: b * s, z: c * s, }
}

///|
pub fn interpolate1(b : Vec3, a0 : Double, a1 : Double, a2 : Double) -> Double {
  b.x * a0 + b.y * a1 + b.z * a2
}

///|
pub fn interpolate2(b : Vec3, a0 : Vec2, a1 : Vec2, a2 : Vec2) -> Vec2 {
  a0.scale(b.x) + a1.scale(b.y) + a2.scale(b.z)
}

///|
pub fn interpolate3(b : Vec3, a0 : Vec3, a1 : Vec3, a2 : Vec3) -> Vec3 {
  a0.scale(b.x) + a1.scale(b.y) + a2.scale(b.z)
}
