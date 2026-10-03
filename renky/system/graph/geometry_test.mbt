// Rays, planes, boxes, frusta.

///|
test "Plane: construction, signed distance, orientation" {
  let a = @graph.Vec3::new(1.0, 0.0, 0.0)
  let b = @graph.Vec3::new(0.0, 1.0, 0.0)
  let c = @graph.Vec3::new(0.0, 0.0, 1.0)
  let p = @graph.Plane::from_points(a, b, c)
  let s = 1.0 / 3.0.sqrt()
  check_v3(
    p.normal,
    @graph.Vec3::new(s, s, s),
    1.0e-15,
    "normal (CCW => outward)",
  )
  check_close(p.signed_distance(a), 0.0, 1.0e-15, "a on plane")
  check_close(p.signed_distance(b), 0.0, 1.0e-15, "b on plane")
  check_close(p.signed_distance(c), 0.0, 1.0e-15, "c on plane")
  check_close(
    p.signed_distance(@graph.Vec3::zero()),
    -s,
    1.0e-15,
    "origin is behind",
  )
  check_close(
    p.signed_distance(@graph.Vec3::new(1.0, 1.0, 1.0)),
    3.0 * s - s,
    1.0e-15,
    "front",
  )
  let q = @graph.Plane::from_point_normal(
    @graph.Vec3::new(0.0, 5.0, 0.0),
    @graph.Vec3::new(0.0, 3.0, 0.0),
  )
  check_close(q.d, -5.0, 0.0, "d")
  check_close(
    q.signed_distance(@graph.Vec3::new(9.0, 7.0, -2.0)),
    2.0,
    0.0,
    "height above ground plane",
  )
}

///|
test "Ray: plane intersection (hit, behind, parallel), at()" {
  let ground = @graph.Plane::from_point_normal(
    @graph.Vec3::zero(),
    @graph.Vec3::unit_y(),
  )
  let r = @graph.Ray::new(
    @graph.Vec3::new(1.0, 10.0, 2.0),
    @graph.Vec3::new(0.0, -2.0, 0.0),
  )
  check_close(
    r.intersect_plane(ground).unwrap(),
    10.0,
    1.0e-15,
    "straight down",
  )
  check_v3(r.at(10.0), @graph.Vec3::new(1.0, 0.0, 2.0), 1.0e-15, "at")
  let slanted = @graph.Ray::new(
    @graph.Vec3::new(0.0, 3.0, 0.0),
    @graph.Vec3::new(1.0, -1.0, 0.0),
  )
  check_close(
    slanted.intersect_plane(ground).unwrap(),
    3.0 * 2.0.sqrt(),
    1.0e-14,
    "slanted distance",
  )
  // pointing away: intersection is behind the origin -> None
  check_true(
    @graph.Ray::new(@graph.Vec3::new(0.0, 3.0, 0.0), @graph.Vec3::unit_y()).intersect_plane(
      ground,
    )
    is None,
    "behind",
  )
  // very small positive distances are still hits
  let near = @graph.Ray::new(
    @graph.Vec3::new(0.0, 1.0e-3, 0.0),
    -@graph.Vec3::unit_y(),
  )
  check_close(near.intersect_plane(ground).unwrap(), 1.0e-3, 1.0e-12, "tiny t")
  // parallel
  check_true(
    @graph.Ray::new(@graph.Vec3::new(0.0, 3.0, 0.0), @graph.Vec3::unit_x()).intersect_plane(
      ground,
    )
    is None,
    "parallel",
  )
  // origin ON the plane -> t = 0
  check_close(
    @graph.Ray::new(@graph.Vec3::zero(), -@graph.Vec3::unit_y())
    .intersect_plane(ground)
    .unwrap(),
    0.0,
    0.0,
    "on plane",
  )
  // ray constructor normalises
  check_close(
    @graph.Ray::new(@graph.Vec3::zero(), @graph.Vec3::new(3.0, 4.0, 0.0)).dir.length(),
    1.0,
    1.0e-15,
    "unit dir",
  )
}

///|
test "Ray: sphere" {
  let c = @graph.Vec3::new(0.0, 0.0, -10.0)
  let r = @graph.Ray::new(@graph.Vec3::zero(), -@graph.Vec3::unit_z())
  check_close(r.intersect_sphere(c, 2.0).unwrap(), 8.0, 1.0e-14, "front hit")
  // starting inside: exit point
  let inside = @graph.Ray::new(c, -@graph.Vec3::unit_z())
  check_close(
    inside.intersect_sphere(c, 2.0).unwrap(),
    2.0,
    1.0e-14,
    "from inside",
  )
  // miss, and sphere behind
  check_true(
    @graph.Ray::new(@graph.Vec3::new(5.0, 0.0, 0.0), -@graph.Vec3::unit_z()).intersect_sphere(
      c, 2.0,
    )
    is None,
    "miss",
  )
  check_true(
    @graph.Ray::new(@graph.Vec3::zero(), @graph.Vec3::unit_z()).intersect_sphere(
      c, 2.0,
    )
    is None,
    "behind",
  )
  // tangent graze
  check_close(
    @graph.Ray::new(@graph.Vec3::new(2.0, 0.0, 0.0), -@graph.Vec3::unit_z())
    .intersect_sphere(c, 2.0)
    .unwrap(),
    10.0,
    1.0e-6,
    "tangent",
  )
  // random: hit point lies on the sphere
  let rng = Rng::new(40UL)
  for i in 0..<200 {
    let ctr = rng.vec3(-5.0, 5.0)
    let rad = rng.range(0.5, 3.0)
    let o = ctr + rng.unit_vec3().scale(rng.range(rad * 1.5, 12.0))
    let ray = @graph.Ray::new(o, ctr - o + rng.vec3(-0.3, 0.3))
    match ray.intersect_sphere(ctr, rad) {
      Some(t) =>
        check_close(ray.at(t).distance(ctr), rad, 1.0e-12, "on surface #\{i}")
      None => ()
    }
  }
}

///|
test "Aabb: transform is the exact bounding box of the transformed corners" {
  let rng = Rng::new(41UL)
  for i in 0..<300 {
    let c = rng.vec3(-5.0, 5.0)
    let h = rng.vec3(0.1, 3.0)
    let b = @graph.Aabb::from_center_half(c, h)
    let m = @graph.Mat4::from_trs(
      rng.vec3(-9.0, 9.0),
      rng.quat(),
      @graph.Vec3::new(
        rng.range(0.3, 3.0),
        -rng.range(0.3, 3.0),
        rng.range(0.3, 3.0),
      ),
    )
    let t = b.transform(m)
    // brute force over the 8 corners
    let mut mn = m.mul_point(b.corner(0))
    let mut mx = mn
    for k in 1..<8 {
      let p = m.mul_point(b.corner(k))
      mn = mn.min(p)
      mx = mx.max(p)
    }
    check_v3(t.min, mn, 1.0e-12, "min #\{i}")
    check_v3(t.max, mx, 1.0e-12, "max #\{i}")
  }
  let a = @graph.Aabb::new(@graph.Vec3::zero(), @graph.Vec3::new(1.0, 2.0, 3.0))
  check_v3(a.center(), @graph.Vec3::new(0.5, 1.0, 1.5), 0.0, "center")
  check_v3(a.half_extents(), @graph.Vec3::new(0.5, 1.0, 1.5), 0.0, "half")
  let u = a.union(
    @graph.Aabb::new(
      @graph.Vec3::new(-1.0, 1.0, 0.0),
      @graph.Vec3::new(0.5, 5.0, 1.0),
    ),
  )
  check_v3(u.min, @graph.Vec3::new(-1.0, 0.0, 0.0), 0.0, "union min")
  check_v3(u.max, @graph.Vec3::new(1.0, 5.0, 3.0), 0.0, "union max")
  let e = a.expand(@graph.Vec3::new(4.0, -1.0, 1.0))
  check_v3(e.min, @graph.Vec3::new(0.0, -1.0, 0.0), 0.0, "expand min")
  check_v3(e.max, @graph.Vec3::new(4.0, 2.0, 3.0), 0.0, "expand max")
  check_true(
    a.contains_point(@graph.Vec3::new(1.0, 2.0, 3.0)) &&
    !a.contains_point(@graph.Vec3::new(1.1, 0.0, 0.0)),
    "contains",
  )
  check_v3(a.corner(0), a.min, 0.0, "corner 0")
  check_v3(a.corner(7), a.max, 0.0, "corner 7")
  check_v3(
    a.corner(5),
    @graph.Vec3::new(1.0, 0.0, 3.0),
    0.0,
    "corner 5 (x and z max)",
  )
}

///|
test "Ray: AABB slab test vs dense sampling" {
  let b = @graph.Aabb::new(
    @graph.Vec3::new(-1.0, -2.0, -3.0),
    @graph.Vec3::new(2.0, 1.0, 0.0),
  )
  // straight through the middle along x
  let r = @graph.Ray::new(
    @graph.Vec3::new(-5.0, -0.5, -1.5),
    @graph.Vec3::unit_x(),
  )
  let (t0, t1) = r.intersect_aabb(b).unwrap()
  check_close(t0, 4.0, 1.0e-15, "enter")
  check_close(t1, 7.0, 1.0e-15, "exit")
  // origin inside -> enter clamped to 0
  let (i0, i1) = @graph.Ray::new(
      @graph.Vec3::new(0.0, 0.0, -1.0),
      @graph.Vec3::unit_y(),
    )
    .intersect_aabb(b)
    .unwrap()
  check_close(i0, 0.0, 0.0, "inside enter")
  check_close(i1, 1.0, 1.0e-15, "inside exit")
  // axis-parallel miss and pointing away
  check_true(
    @graph.Ray::new(@graph.Vec3::new(-5.0, 3.0, -1.5), @graph.Vec3::unit_x()).intersect_aabb(
      b,
    )
    is None,
    "parallel miss",
  )
  check_true(
    @graph.Ray::new(@graph.Vec3::new(-5.0, -0.5, -1.5), -@graph.Vec3::unit_x()).intersect_aabb(
      b,
    )
    is None,
    "away",
  )
  let rng = Rng::new(43UL)
  let mut hits = 0
  let mut misses = 0
  for i in 0..<400 {
    let o = rng.vec3(-6.0, 6.0)
    // aim near the box with some noise so that both hits and misses are common
    let ray = @graph.Ray::new(o, b.center() - o + rng.vec3(-2.5, 2.5))
    let res = ray.intersect_aabb(b)
    // brute-force by marching
    let mut found = false
    let mut t = 0.0
    while t < 30.0 {
      if b.contains_point(ray.at(t)) {
        found = true
        break
      }
      t += 0.002
    }
    match res {
      Some((a0, a1)) => {
        hits += 1
        check_true(a0 <= a1, "ordered #\{i}")
        // midpoint of the interval is inside the box
        check_true(
          b.contains_point(ray.at((a0 + a1) * 0.5)),
          "mid inside #\{i}",
        )
        check_true(found || a1 - a0 < 0.01, "marching found it #\{i}")
      }
      None => {
        misses += 1
        check_true(!found, "marching also misses #\{i}")
      }
    }
  }
  check_true(
    hits > 40 && misses > 40,
    "exercised hits=\{hits} misses=\{misses}",
  )
}

///|
test "Ray::transform equals transforming sampled points" {
  let rng = Rng::new(44UL)
  let m = @graph.Mat4::from_trs(
    rng.vec3(-9.0, 9.0),
    rng.quat(),
    @graph.Vec3::splat(1.0),
  )
  let r = @graph.Ray::new(rng.vec3(-3.0, 3.0), rng.unit_vec3())
  let rt = r.transform(m)
  for t in [0.0, 1.0, 7.5] {
    check_v3(rt.at(t), m.mul_point(r.at(t)), 1.0e-12, "point on ray")
  }
}

///|
test "Frustum planes face inward and sit at the right distance (near/far identity per depth mode)" {
  for m in all_modes() {
    let n = 0.5
    let f = 80.0
    let cam = @graph.Camera::look_at(
      @graph.Vec3::new(1.0, 2.0, 3.0),
      @graph.Vec3::new(4.0, 2.0, -9.0),
      @graph.Vec3::unit_y(),
      @graph.Projection::perspective(1.0, 1.5, n, f, m),
      @graph.Viewport::full(64.0, 64.0),
    )
    let fr = cam.frustum()
    let fwd = cam.forward()
    check_v3(fr.near.normal, fwd, 1.0e-12, "near normal faces forward")
    check_v3(fr.far.normal, -fwd, 1.0e-12, "far normal faces back")
    check_close(
      fr.near.signed_distance(cam.eye() + fwd.scale(n)),
      0.0,
      1.0e-12,
      "near plane at n",
    )
    check_close(
      fr.far.signed_distance(cam.eye() + fwd.scale(f)),
      0.0,
      1.0e-10,
      "far plane at f",
    )
    check_close(
      fr.near.signed_distance(cam.eye() + fwd.scale(2.0)),
      2.0 - n,
      1.0e-12,
      "metric distance",
    )
    // side planes pass through the eye
    for pl in [fr.left, fr.right, fr.bottom, fr.top] {
      check_close(
        pl.signed_distance(cam.eye()),
        0.0,
        1.0e-12,
        "side plane through eye",
      )
      check_close(pl.normal.length(), 1.0, 1.0e-14, "unit normal")
    }
    // the centre of the volume is inside every plane
    check_true(fr.contains_point(cam.eye() + fwd.scale(10.0)), "centre inside")
    // half-angle: the right plane is tilted by fov_x/2 from forward
    let half_x = @math.atan(@math.tan(0.5) * 1.5)
    check_close(
      fr.right.normal.dot(fwd),
      @math.sin(half_x),
      1.0e-12,
      "right plane tilt (inward normal leans forward)",
    )
  }
  // infinite far plane: never culls on distance
  let cam = @graph.Camera::look_at(
    @graph.Vec3::zero(),
    @graph.Vec3::new(0.0, 0.0, -1.0),
    @graph.Vec3::unit_y(),
    @graph.Projection::perspective(1.0, 1.0, 0.1, @double.infinity, ReversedZ),
    @graph.Viewport::full(8.0, 8.0),
  )
  check_true(
    cam.frustum().contains_point(@graph.Vec3::new(0.0, 0.0, -1.0e9)),
    "infinite far",
  )
  check_true(
    !cam.frustum().contains_point(@graph.Vec3::new(0.0, 0.0, -0.05)),
    "inside near distance",
  )
}
