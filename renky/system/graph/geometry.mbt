// Geometric primitives that live *in* the coordinate spaces: rays, planes, boxes, frusta.

///|
#valtype
pub(all) struct Ray {
  origin : Vec3
  /// Unit length for the intersection routines to return true distances.
  dir : Vec3
} derive(Eq, Debug)

///|
pub fn Ray::new(origin : Vec3, dir : Vec3) -> Ray {
  { origin, dir: dir.normalize(), }
}

///|
pub fn Ray::at(self : Ray, t : Double) -> Vec3 {
  self.origin + self.dir.scale(t)
}

///|
/// Apply a rigid / affine matrix to the ray (dir re-normalised; distances t scale by the
/// returned factor if the matrix scales).
pub fn Ray::transform(self : Ray, m : Mat4) -> Ray {
  { origin: m.mul_point(self.origin), dir: m.mul_dir(self.dir).normalize(), }
}

///|
/// Plane `dot(normal, p) + d = 0`; `normal` should be unit for `signed_distance` to be metric.
#valtype
pub(all) struct Plane {
  normal : Vec3
  d : Double
} derive(Eq, Debug)

///|
pub fn Plane::from_point_normal(point : Vec3, normal : Vec3) -> Plane {
  let n = normal.normalize()
  { normal: n, d: -n.dot(point), }
}

///|
/// Plane through three points (normal = (b-a) x (c-a), counter-clockwise => facing viewer).
pub fn Plane::from_points(a : Vec3, b : Vec3, c : Vec3) -> Plane {
  Plane::from_point_normal(a, (b - a).cross(c - a))
}

///|
pub fn Plane::signed_distance(self : Plane, p : Vec3) -> Double {
  self.normal.dot(p) + self.d
}

///|
/// Distance t >= 0 along the ray to the plane; `None` if parallel or behind the origin.
pub fn Ray::intersect_plane(self : Ray, pl : Plane) -> Double? {
  let denom = pl.normal.dot(self.dir)
  if denom == 0.0 {
    return None
  }
  let t = -(pl.normal.dot(self.origin) + pl.d) / denom
  if t >= 0.0 {
    Some(t)
  } else {
    None
  }
}

///|
/// Nearest distance t >= 0 to a sphere surface hit (entry, or exit if starting inside).
pub fn Ray::intersect_sphere(
  self : Ray,
  center : Vec3,
  radius : Double,
) -> Double? {
  let oc = self.origin - center
  let b = oc.dot(self.dir)
  let c = oc.dot(oc) - radius * radius
  let disc = b * b - c
  if disc < 0.0 {
    return None
  }
  let s = disc.sqrt()
  let t0 = -b - s
  let t1 = -b + s
  if t0 >= 0.0 {
    Some(t0)
  } else if t1 >= 0.0 {
    Some(t1)
  } else {
    None
  }
}

///|
/// Axis-aligned bounding box.
#valtype
pub(all) struct Aabb {
  min : Vec3
  max : Vec3
} derive(Eq, Debug)

///|
pub fn Aabb::new(min : Vec3, max : Vec3) -> Aabb {
  { min, max, }
}

///|
pub fn Aabb::from_center_half(center : Vec3, half : Vec3) -> Aabb {
  { min: center - half, max: center + half, }
}

///|
pub fn Aabb::center(self : Aabb) -> Vec3 {
  (self.min + self.max).scale(0.5)
}

///|
pub fn Aabb::half_extents(self : Aabb) -> Vec3 {
  (self.max - self.min).scale(0.5)
}

///|
pub fn Aabb::union(self : Aabb, o : Aabb) -> Aabb {
  { min: self.min.min(o.min), max: self.max.max(o.max), }
}

///|
pub fn Aabb::expand(self : Aabb, p : Vec3) -> Aabb {
  { min: self.min.min(p), max: self.max.max(p), }
}

///|
pub fn Aabb::contains_point(self : Aabb, p : Vec3) -> Bool {
  p.x >= self.min.x &&
  p.x <= self.max.x &&
  p.y >= self.min.y &&
  p.y <= self.max.y &&
  p.z >= self.min.z &&
  p.z <= self.max.z
}

///|
/// Corner `i` in 0..7 (bit0 -> x, bit1 -> y, bit2 -> z select max).
pub fn Aabb::corner(self : Aabb, i : Int) -> Vec3 {
  {
    x: if (i & 1) != 0 {
      self.max.x
    } else {
      self.min.x
    },
    y: if (i & 2) != 0 {
      self.max.y
    } else {
      self.min.y
    },
    z: if (i & 4) != 0 {
      self.max.z
    } else {
      self.min.z
    },
  }
}

///|
/// Exact world-space AABB of this box under an affine matrix (Arvo): cost of one 3x3
/// multiply, no corner enumeration.
pub fn Aabb::transform(self : Aabb, m : Mat4) -> Aabb {
  let c = self.center()
  let h = self.half_extents()
  let nc = m.mul_point(c)
  let nh : Vec3 = {
    x: m.m00.abs() * h.x + m.m01.abs() * h.y + m.m02.abs() * h.z,
    y: m.m10.abs() * h.x + m.m11.abs() * h.y + m.m12.abs() * h.z,
    z: m.m20.abs() * h.x + m.m21.abs() * h.y + m.m22.abs() * h.z,
  }
  { min: nc - nh, max: nc + nh, }
}

///|
/// Slab test. Returns (t_enter, t_exit) with t_exit >= max(t_enter, 0); t_enter is clamped to 0
/// when the origin is inside. `None` on a miss.
pub fn Ray::intersect_aabb(self : Ray, b : Aabb) -> (Double, Double)? {
  let mut t0 = 0.0
  let mut t1 = @double.infinity
  for i in 0..<3 {
    let o = self.origin.get(i)
    let d = self.dir.get(i)
    let lo = b.min.get(i)
    let hi = b.max.get(i)
    if d == 0.0 {
      if o < lo || o > hi {
        return None
      }
    } else {
      let inv = 1.0 / d
      let mut a = (lo - o) * inv
      let mut c = (hi - o) * inv
      if a > c {
        let tmp = a
        a = c
        c = tmp
      }
      if a > t0 {
        t0 = a
      }
      if c < t1 {
        t1 = c
      }
      if t0 > t1 {
        return None
      }
    }
  }
  Some((t0, t1))
}

///|
pub(all) enum Containment {
  Outside
  Intersecting
  Inside
} derive(Eq, Debug)

///|
/// The six planes of a view frustum with normals pointing INSIDE.
#valtype
pub(all) struct Frustum {
  left : Plane
  right : Plane
  bottom : Plane
  top : Plane
  near : Plane
  far : Plane
} derive(Eq, Debug)

///|
fn plane_from_vec4(v : Vec4) -> Plane {
  let n = Vec3::new(v.x, v.y, v.z)
  let l = n.length()
  if l > 0.0 {
    { normal: n.scale(1.0 / l), d: v.w / l, }
  } else {
    // degenerate (infinite far plane): constant, always satisfied when d >= 0
    { normal: n, d: v.w, }
  }
}

///|
/// Extract frustum planes from `clip = view_proj * world`; handles every `DepthMode`.
/// The planes are expressed in the space that `view_proj` takes as INPUT (world space when
/// `view_proj = P * V`; model space when `P * V * M`).
pub fn Frustum::from_view_proj(m : Mat4, depth : DepthMode) -> Frustum {
  let r0 = m.row(0)
  let r1 = m.row(1)
  let r2 = m.row(2)
  let r3 = m.row(3)
  let (near, far) = match depth {
    ZeroToOne => (r2, r3 - r2)
    NegOneToOne => (r3 + r2, r3 - r2)
    ReversedZ => (r3 - r2, r2)
  }
  {
    left: plane_from_vec4(r3 + r0),
    right: plane_from_vec4(r3 - r0),
    bottom: plane_from_vec4(r3 + r1),
    top: plane_from_vec4(r3 - r1),
    near: plane_from_vec4(near),
    far: plane_from_vec4(far),
  }
}

///|
pub fn Frustum::contains_point(self : Frustum, p : Vec3) -> Bool {
  self.left.signed_distance(p) >= 0.0 &&
  self.right.signed_distance(p) >= 0.0 &&
  self.bottom.signed_distance(p) >= 0.0 &&
  self.top.signed_distance(p) >= 0.0 &&
  self.near.signed_distance(p) >= 0.0 &&
  self.far.signed_distance(p) >= 0.0
}

///|
/// Conservative sphere test (never a false negative; rare false positives near frustum edges).
pub fn Frustum::intersects_sphere(self : Frustum, c : Vec3, r : Double) -> Bool {
  self.left.signed_distance(c) >= -r &&
  self.right.signed_distance(c) >= -r &&
  self.bottom.signed_distance(c) >= -r &&
  self.top.signed_distance(c) >= -r &&
  self.near.signed_distance(c) >= -r &&
  self.far.signed_distance(c) >= -r
}

///|
fn classify_plane(pl : Plane, b : Aabb) -> Int {
  // 0 = fully outside, 1 = straddles, 2 = fully inside
  let pv : Vec3 = {
    x: if pl.normal.x >= 0.0 {
      b.max.x
    } else {
      b.min.x
    },
    y: if pl.normal.y >= 0.0 {
      b.max.y
    } else {
      b.min.y
    },
    z: if pl.normal.z >= 0.0 {
      b.max.z
    } else {
      b.min.z
    },
  }
  if pl.signed_distance(pv) < 0.0 {
    return 0
  }
  let nv : Vec3 = {
    x: if pl.normal.x >= 0.0 {
      b.min.x
    } else {
      b.max.x
    },
    y: if pl.normal.y >= 0.0 {
      b.min.y
    } else {
      b.max.y
    },
    z: if pl.normal.z >= 0.0 {
      b.min.z
    } else {
      b.max.z
    },
  }
  if pl.signed_distance(nv) < 0.0 {
    1
  } else {
    2
  }
}

///|
/// Box classification. `Outside` is exact-safe (box is certainly invisible); `Inside` is
/// certain; `Intersecting` may include a few boxes that are actually just outside a corner.
pub fn Frustum::classify_aabb(self : Frustum, b : Aabb) -> Containment {
  let mut all_in = true
  for pl in [self.left, self.right, self.bottom, self.top, self.near, self.far] {
    match classify_plane(pl, b) {
      0 => return Outside
      1 => all_in = false
      _ => ()
    }
  }
  if all_in {
    Inside
  } else {
    Intersecting
  }
}
