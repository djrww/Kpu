// Projection / viewport / camera / frustum tests.
//   * numeric agreement with an independent numpy pipeline (golden_data_test.mbt)
//   * geometric invariants that must hold for ANY correct implementation

///|
fn all_modes() -> Array[@graph.DepthMode] {
  [@graph.DepthMode::ZeroToOne, NegOneToOne, ReversedZ]
}

///|
fn mode_of(i : Int) -> @graph.DepthMode {
  match i {
    0 => ZeroToOne
    1 => NegOneToOne
    _ => ReversedZ
  }
}

///|
fn gold_camera(g : FixedArray[Double], o : Int) -> @graph.Camera {
  let far = if g[o + 12] == 0.0 { @double.infinity } else { g[o + 12] }
  @graph.Camera::look_at(
    g_v3(g, o),
    g_v3(g, o + 3),
    g_v3(g, o + 6),
    @graph.Projection::perspective(
      g[o + 9],
      g[o + 10],
      g[o + 11],
      far,
      mode_of(g[o + 15].to_int()),
    ),
    @graph.Viewport::full(g[o + 13], g[o + 14]),
  )
}

///|
test "Camera pipeline == independent numpy pipeline (3 depth modes, finite & infinite far)" {
  let g = gold_pipe
  let k = 6
  let stride = 16 + k * 10
  for c in 0..<(g.length() / stride) {
    let o = c * stride
    let cam = gold_camera(g, o)
    let mode = mode_of(g[o + 15].to_int())
    for i in 0..<k {
      let q = o + 16 + i * 10
      let p = g_v3(g, q)
      let tag = "case \{c} pt \{i}"
      let want_s = g_v3(g, q + 3)
      let want_ndc = g_v3(g, q + 6)
      let want_vz = g[q + 9]
      // view space z
      check_close(cam.world_to_view(p).z, want_vz, 1.0e-10, "view z " + tag)
      // clip -> ndc through the 4x4 matrices
      let clip = cam.world_to_clip(p)
      match clip.perspective_divide() {
        Some(n) => check_v3(n, want_ndc, 1.0e-9, "ndc(matrix) " + tag)
        None => fail("w==0 " + tag)
      }
      // direct (matrix-free) projection
      match cam.world_to_ndc(p) {
        Some(n) => check_v3(n, want_ndc, 1.0e-9, "ndc(direct) " + tag)
        None => fail("None " + tag)
      }
      match cam.world_to_screen(p) {
        Some(s) => check_v3(s, want_s, 1.0e-8, "screen " + tag)
        None => fail("screen None " + tag)
      }
      // round trip
      let back = cam.screen_to_world(want_s.x, want_s.y, want_s.z)
      check_v3(back, p, 1.0e-6, "screen_to_world " + tag)
      // picking ray passes through the point
      let ray = cam.pick_ray(want_s.x, want_s.y)
      let t = (p - ray.origin).dot(ray.dir)
      check_true(t > 0.0, "point in front of ray " + tag)
      check_v3(ray.at(t), p, 1.0e-6, "ray hits point " + tag)
      // linearised depth
      check_close(
        cam.projection.linear_depth(want_ndc.z),
        -want_vz,
        1.0e-7,
        "linear depth " + tag,
      )
      check_close(
        cam.projection.ndc_depth(want_vz),
        want_ndc.z,
        1.0e-9,
        "ndc_depth " + tag,
      )
      let _ = mode
    }
  }
}

///|
test "Orthographic matrix and mapping == glOrtho / XMMatrixOrthographicOffCenterRH" {
  let g = gold_ortho
  let stride = 7 + 5 * 6
  for c in 0..<(g.length() / stride) {
    let o = c * stride
    let proj = @graph.Projection::orthographic(
      g[o],
      g[o + 1],
      g[o + 2],
      g[o + 3],
      g[o + 4],
      g[o + 5],
      mode_of(g[o + 6].to_int()),
    )
    let m = proj.to_mat4()
    for i in 0..<5 {
      let q = o + 7 + i * 6
      let pv = g_v3(g, q)
      let want = g_v3(g, q + 3)
      check_v3(m.mul_point(pv), want, 1.0e-12, "ortho matrix \{c}/\{i}")
      match proj.project_view(pv) {
        Some(n) => check_v3(n, want, 1.0e-12, "ortho project_view \{c}/\{i}")
        None => fail("none")
      }
      check_v3(
        proj.unproject_ndc(want),
        pv,
        1.0e-9,
        "ortho unproject \{c}/\{i}",
      )
      check_close(
        proj.view_z_from_ndc(want.z),
        pv.z,
        1.0e-9,
        "ortho view z \{c}/\{i}",
      )
    }
    check_true(m.m32 == 0.0 && m.m33 == 1.0, "affine last row")
  }
}

///|
test "known matrices: gluPerspective(90deg, 1, 1, 3) in all three depth layouts" {
  let p = @graph.HALF_PI
  let gl = @graph.Projection::perspective(p, 1.0, 1.0, 3.0, NegOneToOne).to_mat4()
  check_close(gl.m00, 1.0, 1.0e-15, "gl m00")
  check_close(gl.m11, 1.0, 1.0e-15, "gl m11")
  check_close(gl.m22, -2.0, 1.0e-15, "gl m22 = -(f+n)/(f-n)")
  check_close(gl.m23, -3.0, 1.0e-15, "gl m23 = -2fn/(f-n)")
  check_close(gl.m32, -1.0, 0.0, "gl m32")
  check_close(gl.m33, 0.0, 0.0, "gl m33")
  let d3 = @graph.Projection::perspective(p, 1.0, 1.0, 3.0, ZeroToOne).to_mat4()
  check_close(d3.m22, -1.5, 1.0e-15, "d3d m22 = -f/(f-n)")
  check_close(d3.m23, -1.5, 1.0e-15, "d3d m23 = -fn/(f-n)")
  let rz = @graph.Projection::perspective(p, 1.0, 1.0, 3.0, ReversedZ).to_mat4()
  check_close(rz.m22, 0.5, 1.0e-15, "rev m22 = n/(f-n)")
  check_close(rz.m23, 1.5, 1.0e-15, "rev m23 = nf/(f-n)")
  let ri = @graph.Projection::perspective(
    p,
    1.0,
    0.1,
    @double.infinity,
    ReversedZ,
  ).to_mat4()
  check_close(ri.m22, 0.0, 0.0, "reversed infinite m22")
  check_close(ri.m23, 0.1, 1.0e-15, "reversed infinite m23 = near")
  // fov 90 => at distance d the half-width is d
  check_v3(
    @graph.Projection::perspective(p, 1.0, 1.0, 3.0, ZeroToOne)
    .project_view(@graph.Vec3::new(2.0, -2.0, -2.0))
    .unwrap(),
    @graph.Vec3::new(1.0, -1.0, (-1.5 * -2.0 - 1.5) / 2.0),
    1.0e-15,
    "90deg edge",
  )
}

///|
test "depth: near/far planes land on the documented NDC depth for every mode" {
  let n = 0.25
  let f = 400.0
  for m in all_modes() {
    let pr = @graph.Projection::perspective(1.1, 1.7, n, f, m)
    let po = @graph.Projection::orthographic_centered(8.0, 6.0, n, f, m)
    let (dn, df) = match m {
      ZeroToOne => (0.0, 1.0)
      NegOneToOne => (-1.0, 1.0)
      ReversedZ => (1.0, 0.0)
    }
    for p in [pr, po] {
      check_close(p.ndc_depth(-n), dn, 1.0e-12, "near")
      check_close(p.ndc_depth(-f), df, 1.0e-12, "far")
      check_close(p.view_z_from_ndc(dn), -n, 1.0e-9, "near inv")
      check_close(p.view_z_from_ndc(df), -f, 1.0e-6, "far inv")
      // monotone in distance and bounded by [near, far] mapping
      let mut prev = p.ndc_depth(-n)
      let mut z = n
      for _ in 0..<50 {
        z = z * 1.15
        if z > f {
          break
        }
        let cur = p.ndc_depth(-z)
        check_true(
          if m == ReversedZ {
            cur < prev
          } else {
            cur > prev
          },
          "monotone",
        )
        prev = cur
      }
    }
  }
  // infinite far: depth -> far value, never crosses it
  let pi = @graph.Projection::perspective(
    1.0,
    1.0,
    0.1,
    @double.infinity,
    ReversedZ,
  )
  check_close(pi.ndc_depth(-0.1), 1.0, 1.0e-15, "inf near")
  check_close(pi.ndc_depth(-1.0e12), 1.0e-13, 1.0e-12, "inf far -> 0")
  check_true(pi.ndc_depth(-1.0e12) > 0.0, "positive")
  let pz = @graph.Projection::perspective(
    1.0,
    1.0,
    0.1,
    @double.infinity,
    ZeroToOne,
  )
  check_close(pz.ndc_depth(-1.0e9), 1.0, 1.0e-9, "inf far -> 1")
}

///|
test "Reversed-Z keeps far more depth resolution than standard Z in float32" {
  // Quantise ndc depth to float32 and count distinct values over [1m, 10km] on a 0.1..10000 range.
  let n = 0.1
  let f = 10000.0
  let std = @graph.Projection::perspective(1.0, 1.0, n, f, ZeroToOne)
  let rev = @graph.Projection::perspective(1.0, 1.0, n, f, ReversedZ)
  // Smallest distance step that still changes the float32 depth, at 5 km.
  let mut std_step = 0.0
  let mut rev_step = 0.0
  let z = 5000.0
  let mut d = 0.001
  while d < 1000.0 {
    let a = Float::from_double(std.ndc_depth(-z)).to_double()
    let b = Float::from_double(std.ndc_depth(-(z + d))).to_double()
    if a != b && std_step == 0.0 {
      std_step = d
    }
    let c = Float::from_double(rev.ndc_depth(-z)).to_double()
    let e = Float::from_double(rev.ndc_depth(-(z + d))).to_double()
    if c != e && rev_step == 0.0 {
      rev_step = d
    }
    d = d * 1.05
  }
  check_true(rev_step > 0.0 && std_step > 0.0, "both resolvable")
  check_true(
    rev_step * 8.0 < std_step,
    "reversed-Z resolves much finer steps at distance: rev=\{rev_step} std=\{std_step}",
  )
}

///|
test "Projection: closed-form inverse == general inverse; P * P^-1 = I; all kinds & modes" {
  let rng = Rng::new(31UL)
  for i in 0..<60 {
    let m = all_modes()[i % 3]
    let n = rng.range(0.05, 3.0)
    let f = n + rng.range(5.0, 3000.0)
    let projs : Array[@graph.Projection] = [
      @graph.Projection::perspective(
        rng.range(0.3, 2.5),
        rng.range(0.5, 3.0),
        n,
        f,
        m,
      ),
      @graph.Projection::perspective(
        rng.range(0.3, 2.5),
        rng.range(0.5, 3.0),
        n,
        @double.infinity,
        m,
      ),
      @graph.Projection::frustum(
        -rng.range(0.1, 2.0) * n,
        rng.range(0.1, 2.0) * n,
        -rng.range(0.1, 2.0) * n,
        rng.range(0.1, 2.0) * n,
        n,
        f,
        m,
      ),
      @graph.Projection::orthographic(
        -rng.range(1.0, 9.0),
        rng.range(1.0, 9.0),
        -rng.range(1.0, 9.0),
        rng.range(1.0, 9.0),
        n,
        f,
        m,
      ),
    ]
    for p in projs {
      let pm = p.to_mat4()
      let inv = p.inverse_mat4()
      check_m4(pm * inv, @graph.Mat4::identity(), 1.0e-9, "P*Pinv #\{i}")
      check_m4(inv * pm, @graph.Mat4::identity(), 1.0e-9, "Pinv*P #\{i}")
      match pm.inverse() {
        Some(g) => check_m4(inv, g, 1.0e-6, "closed form vs general #\{i}")
        None => fail("none")
      }
      // project/unproject round trip at random view points in the volume
      for _ in 0..<5 {
        let zv = -rng.range(n * 1.01, if f > 1.0e6 { 500.0 } else { f * 0.99 })
        let s = if p.is_perspective() { -zv / n } else { 1.0 }
        let pv = @graph.Vec3::new(
          rng.range(p.left, p.right) * s,
          rng.range(p.bottom, p.top) * s,
          zv,
        )
        match p.project_view(pv) {
          Some(ndc) => {
            check_v3(p.unproject_ndc(ndc), pv, 1.0e-7, "unproject #\{i}")
            check_v3(
              pm.project_point(pv).unwrap(),
              ndc,
              1.0e-9,
              "matrix vs direct #\{i}",
            )
            // |x|,|y| <= 1 inside the frustum side planes
            check_true(
              ndc.x.abs() <= 1.0 + 1.0e-9 && ndc.y.abs() <= 1.0 + 1.0e-9,
              "inside xy #\{i}",
            )
          }
          None => fail("project none")
        }
      }
    }
  }
}

///|
test "Frustum corners and asymmetric FOV map to the NDC cube / square" {
  let n = 0.5
  let f = 100.0
  for m in all_modes() {
    let p = @graph.Projection::frustum(-0.3, 0.7, -0.2, 0.4, n, f, m)
    // the near-plane corners (left,bottom),(right,top) in view space
    let c0 = p.project_view(@graph.Vec3::new(-0.3, -0.2, -n)).unwrap()
    let c1 = p.project_view(@graph.Vec3::new(0.7, 0.4, -n)).unwrap()
    check_close(c0.x, -1.0, 1.0e-14, "l")
    check_close(c0.y, -1.0, 1.0e-14, "b")
    check_close(c1.x, 1.0, 1.0e-14, "r")
    check_close(c1.y, 1.0, 1.0e-14, "t")
    // far-plane corners scale by f/n
    let s = f / n
    let f1 = p.project_view(@graph.Vec3::new(0.7 * s, 0.4 * s, -f)).unwrap()
    check_close(f1.x, 1.0, 1.0e-12, "far r")
    check_close(f1.y, 1.0, 1.0e-12, "far t")
  }
  // OpenXR-style half angles
  let a = @graph.Projection::perspective_fov(
    0.8,
    0.6,
    0.5,
    0.7,
    0.1,
    50.0,
    ZeroToOne,
  )
  let right_edge = @graph.Vec3::new(@math.tan(0.6) * 10.0, 0.0, -10.0)
  check_close(
    a.project_view(right_edge).unwrap().x,
    1.0,
    1.0e-13,
    "right half-angle",
  )
  let left_edge = @graph.Vec3::new(-@math.tan(0.8) * 10.0, 0.0, -10.0)
  check_close(
    a.project_view(left_edge).unwrap().x,
    -1.0,
    1.0e-13,
    "left half-angle",
  )
  let up_edge = @graph.Vec3::new(0.0, @math.tan(0.7) * 10.0, -10.0)
  check_close(a.project_view(up_edge).unwrap().y, 1.0, 1.0e-13, "up half-angle")
  let dn_edge = @graph.Vec3::new(0.0, -@math.tan(0.5) * 10.0, -10.0)
  check_close(
    a.project_view(dn_edge).unwrap().y,
    -1.0,
    1.0e-13,
    "down half-angle",
  )
  // fov_y / aspect accessors
  let s = @graph.Projection::perspective(1.2, 1.6, 0.1, 10.0, ZeroToOne)
  check_close(s.fov_y(), 1.2, 1.0e-14, "fov_y")
  check_close(s.aspect(), 1.6, 1.0e-14, "aspect")
  check_close(s.with_aspect(2.0).aspect(), 2.0, 1.0e-14, "with_aspect")
  check_close(
    s.with_aspect(2.0).fov_y(),
    1.2,
    1.0e-14,
    "with_aspect keeps fov_y",
  )
}

///|
test "projection of points on/behind the eye plane is rejected" {
  let p = @graph.Projection::perspective(1.0, 1.0, 0.1, 100.0, ZeroToOne)
  check_true(p.project_view(@graph.Vec3::new(0.0, 0.0, 0.0)) is None, "z=0")
  check_true(p.project_view(@graph.Vec3::new(1.0, 1.0, 2.0)) is None, "behind")
  check_true(
    p.project_view(@graph.Vec3::new(1.0, 1.0, -2.0)) is Some(_),
    "front",
  )
  // orthographic has no such restriction
  let o = @graph.Projection::orthographic_centered(
    4.0,
    4.0,
    -10.0,
    10.0,
    ZeroToOne,
  )
  check_true(
    o.project_view(@graph.Vec3::new(1.0, 1.0, 2.0)) is Some(_),
    "ortho behind ok",
  )
  let cam = @graph.Camera::look_at(
    @graph.Vec3::zero(),
    @graph.Vec3::new(0.0, 0.0, -1.0),
    @graph.Vec3::unit_y(),
    p,
    @graph.Viewport::full(100.0, 100.0),
  )
  check_true(
    cam.world_to_screen(@graph.Vec3::new(0.0, 0.0, 5.0)) is None,
    "camera: behind",
  )
  check_true(
    cam.world_to_screen(@graph.Vec3::new(0.0, 0.0, -5.0)) is Some(_),
    "camera: front",
  )
}

///|
test "Viewport: NDC corners, pixel centres, round trips, DPI" {
  let vp = @graph.Viewport::new(10.0, 20.0, 640.0, 480.0)
  // NDC (-1,+1) is the TOP-LEFT corner of the window, (+1,-1) the bottom-right
  let tl = vp.ndc_to_screen(@graph.Vec3::new(-1.0, 1.0, 0.0), ZeroToOne)
  let br = vp.ndc_to_screen(@graph.Vec3::new(1.0, -1.0, 1.0), ZeroToOne)
  check_v3(tl, @graph.Vec3::new(10.0, 20.0, 0.0), 0.0, "top-left")
  check_v3(br, @graph.Vec3::new(650.0, 500.0, 1.0), 0.0, "bottom-right")
  check_v3(
    vp.ndc_to_screen(@graph.Vec3::zero(), ZeroToOne),
    @graph.Vec3::new(330.0, 260.0, 0.0),
    0.0,
    "centre",
  )
  // GL depth range -1..1 maps to window 0..1
  check_close(
    vp.ndc_to_screen(@graph.Vec3::new(0.0, 0.0, -1.0), NegOneToOne).z,
    0.0,
    0.0,
    "gl near",
  )
  check_close(
    vp.ndc_to_screen(@graph.Vec3::new(0.0, 0.0, 1.0), NegOneToOne).z,
    1.0,
    0.0,
    "gl far",
  )
  // custom depth range
  let vd = { ..vp, min_depth: 0.25, max_depth: 0.75, }
  check_close(
    vd.ndc_to_screen(@graph.Vec3::new(0.0, 0.0, 0.5), ZeroToOne).z,
    0.5,
    1.0e-15,
    "depth range",
  )
  let rng = Rng::new(8UL)
  for m in all_modes() {
    for _ in 0..<50 {
      let n = @graph.Vec3::new(
        rng.range(-1.0, 1.0),
        rng.range(-1.0, 1.0),
        rng.range(0.0, 1.0),
      )
      let s = vd.ndc_to_screen(n, m)
      check_v3(vd.screen_to_ndc(s.x, s.y, s.z, m), n, 1.0e-13, "round trip")
    }
  }
  // pixel centres
  let (i, j) = @graph.pixel_of(@graph.pixel_center(3, 7))
  check_true(i == 3 && j == 7, "pixel_of(center)")
  let (a, b) = @graph.pixel_of(@graph.Vec2::new(-0.1, 4.999))
  check_true(a == -1 && b == 4, "pixel_of floor")
  check_true(
    vp.contains(@graph.Vec2::new(10.0, 20.0)) &&
    !vp.contains(@graph.Vec2::new(650.0, 100.0)),
    "contains half-open",
  )
  check_v2(
    vp.to_uv(@graph.Vec2::new(330.0, 260.0)),
    @graph.Vec2::new(0.5, 0.5),
    0.0,
    "to_uv",
  )
  check_v2(
    vp.from_uv(@graph.Vec2::new(1.0, 1.0)),
    @graph.Vec2::new(650.0, 500.0),
    0.0,
    "from_uv",
  )
  // DPI
  let hd = @graph.Viewport::from_logical(800.5, 600.0, 2.0)
  check_close(hd.width, 1601.0, 0.0, "physical width")
  check_close(hd.height, 1200.0, 0.0, "physical height")
  check_v2(
    @graph.logical_to_physical(@graph.Vec2::new(10.0, 5.0), 1.5),
    @graph.Vec2::new(15.0, 7.5),
    0.0,
    "l2p",
  )
  check_v2(
    @graph.physical_to_logical(@graph.Vec2::new(15.0, 7.5), 1.5),
    @graph.Vec2::new(10.0, 5.0),
    1.0e-15,
    "p2l",
  )
}

///|
test "Camera basis, view matrix, orbit and look_at" {
  let rng = Rng::new(77UL)
  let proj = @graph.Projection::perspective(1.0, 1.5, 0.1, 100.0, ZeroToOne)
  let vp = @graph.Viewport::full(300.0, 200.0)
  for i in 0..<100 {
    let eye = rng.vec3(-50.0, 50.0)
    let tgt = rng.vec3(-50.0, 50.0)
    let cam = @graph.Camera::look_at(eye, tgt, @graph.Vec3::unit_y(), proj, vp)
    check_v3(cam.forward(), (tgt - eye).normalize(), 1.0e-12, "forward #\{i}")
    check_close(cam.right().dot(cam.forward()), 0.0, 1.0e-12, "right ⟂ fwd")
    check_close(cam.up().dot(cam.forward()), 0.0, 1.0e-12, "up ⟂ fwd")
    // right-handed: right x up = back (= -forward)
    check_v3(
      cam.right().cross(cam.up()),
      -cam.forward(),
      1.0e-12,
      "right-handed",
    )
    check_m4(
      cam.view_matrix(),
      @graph.Mat4::look_at(eye, tgt, @graph.Vec3::unit_y()),
      1.0e-11,
      "view vs look_at #\{i}",
    )
    check_m4(
      cam.view_proj(),
      cam.proj_matrix() * cam.view_matrix(),
      1.0e-12,
      "view_proj",
    )
    // the target is on the optical axis => centre of the viewport
    match cam.world_to_screen(tgt) {
      Some(s) => {
        check_close(s.x, 150.0, 1.0e-7, "target x")
        check_close(s.y, 100.0, 1.0e-7, "target y")
      }
      None => fail("target not visible")
    }
  }
  // orbit: yaw 0 pitch 0 looks down -Z from +Z of the target
  let c0 = @graph.Camera::orbit(
    @graph.Vec3::new(1.0, 2.0, 3.0),
    10.0,
    0.0,
    0.0,
    proj,
    vp,
  )
  check_v3(c0.eye(), @graph.Vec3::new(1.0, 2.0, 13.0), 1.0e-12, "orbit eye")
  check_v3(
    c0.forward(),
    @graph.Vec3::new(0.0, 0.0, -1.0),
    1.0e-12,
    "orbit forward",
  )
  // positive pitch looks UP (eye goes below the target)
  let c1 = @graph.Camera::orbit(@graph.Vec3::zero(), 5.0, 0.0, 0.5, proj, vp)
  check_true(c1.forward().y > 0.0 && c1.eye().y < 0.0, "pitch up")
  // yaw by +90 deg turns forward from -Z to -X (right-hand rule about +Y)
  let c2 = @graph.Camera::orbit(
    @graph.Vec3::zero(),
    5.0,
    @graph.HALF_PI,
    0.0,
    proj,
    vp,
  )
  check_v3(c2.forward(), @graph.Vec3::new(-1.0, 0.0, 0.0), 1.0e-12, "yaw 90")
  check_v3(
    c2.pose.rotation.rotate(@graph.Vec3::new(0.0, 0.0, -1.0)),
    c2.forward(),
    1.0e-12,
    "yaw is Quat::from_y",
  )
}

///|
test "Camera: resize keeps aspect, screen_to_world_at_distance, ortho picking" {
  let proj = @graph.Projection::perspective(0.9, 1.0, 0.1, 500.0, ReversedZ)
  let cam = @graph.Camera::look_at(
    @graph.Vec3::new(3.0, 4.0, 5.0),
    @graph.Vec3::zero(),
    @graph.Vec3::unit_y(),
    proj,
    @graph.Viewport::full(100.0, 100.0),
  )
  let big = cam.resize(1920.0, 1080.0)
  check_close(big.projection.aspect(), 1920.0 / 1080.0, 1.0e-14, "aspect")
  check_close(big.projection.fov_y(), 0.9, 1.0e-14, "fov_y preserved")
  // a world point projected then re-placed at its distance returns the same point
  let rng = Rng::new(4UL)
  for _ in 0..<100 {
    let p = cam.eye() +
      cam.forward().scale(rng.range(1.0, 50.0)) +
      cam.right().scale(rng.range(-5.0, 5.0)) +
      cam.up().scale(rng.range(-5.0, 5.0))
    let s = cam.world_to_screen(p).unwrap()
    let dist = (p - cam.eye()).dot(cam.forward())
    check_v3(
      cam.screen_to_world_at_distance(s.x, s.y, dist),
      p,
      1.0e-9,
      "at distance",
    )
  }
  // ortho: rays are parallel with origin on the near plane
  let oc = @graph.Camera::look_at(
    @graph.Vec3::new(0.0, 0.0, 10.0),
    @graph.Vec3::zero(),
    @graph.Vec3::unit_y(),
    @graph.Projection::orthographic_centered(8.0, 6.0, 1.0, 50.0, ZeroToOne),
    @graph.Viewport::full(80.0, 60.0),
  )
  let r0 = oc.pick_ray(0.0, 0.0)
  let r1 = oc.pick_ray(80.0, 60.0)
  check_v3(r0.dir, r1.dir, 1.0e-15, "parallel")
  check_v3(
    r0.origin,
    @graph.Vec3::new(-4.0, 3.0, 9.0),
    1.0e-12,
    "ortho origin top-left on near plane",
  )
  check_v3(
    r1.origin,
    @graph.Vec3::new(4.0, -3.0, 9.0),
    1.0e-12,
    "ortho origin bottom-right",
  )
  let p = @graph.Vec3::new(1.0, -1.0, -3.0)
  let s = oc.world_to_screen(p).unwrap()
  check_v3(
    oc.screen_to_world_at_distance(s.x, s.y, 13.0),
    p,
    1.0e-12,
    "ortho at distance",
  )
}

///|
test "Frustum culling == clip-space inequalities (every depth mode, finite & infinite far)" {
  let rng = Rng::new(1234UL)
  let mut n_in = 0
  let mut n_out = 0
  for ci in 0..<30 {
    let mode = all_modes()[ci % 3]
    let far = if ci % 2 == 0 { 120.0 } else { @double.infinity }
    let proj = @graph.Projection::perspective(
      rng.range(0.5, 1.8),
      rng.range(0.7, 2.0),
      0.3,
      far,
      mode,
    )
    let cam = @graph.Camera::look_at(
      rng.vec3(-10.0, 10.0),
      rng.vec3(-10.0, 10.0),
      @graph.Vec3::unit_y(),
      proj,
      @graph.Viewport::full(64.0, 64.0),
    )
    let fr = cam.frustum()
    let vpm = cam.view_proj()
    for _ in 0..<400 {
      let p = rng.vec3(-60.0, 60.0)
      let c = vpm.mul_vec4(p.to_vec4(1.0))
      let (zlo, zhi) = match mode {
        NegOneToOne => (-c.w, c.w)
        _ => (0.0, c.w)
      }
      // slack of the clip-space test; skip razor-edge points
      let slack = [
        c.x + c.w,
        c.w - c.x,
        c.y + c.w,
        c.w - c.y,
        c.z - zlo,
        zhi - c.z,
      ]
      let mut mn = slack[0]
      for s in slack {
        mn = mn.min(s)
      }
      if mn.abs() < 1.0e-7 {
        continue
      }
      let expect = mn > 0.0
      check_true(
        fr.contains_point(p) == expect,
        "contains_point matches clip test (case \{ci}, mode index \{ci % 3})",
      )
      if expect {
        n_in += 1
      } else {
        n_out += 1
      }
      // sphere with radius 0 must agree, bigger sphere must never lose a visible point
      check_true(fr.intersects_sphere(p, 0.0) == expect, "sphere r=0")
      if expect {
        check_true(
          fr.intersects_sphere(p, 1.0),
          "sphere r=1 keeps visible point",
        )
      }
    }
  }
  check_true(
    n_in > 500 && n_out > 500,
    "test exercised both inside (\{n_in}) and outside (\{n_out}) points",
  )
}

///|
test "Frustum::classify_aabb: no false negatives, Inside implies all corners inside" {
  let rng = Rng::new(555UL)
  let proj = @graph.Projection::perspective(1.0, 1.4, 0.5, 80.0, ReversedZ)
  let cam = @graph.Camera::look_at(
    @graph.Vec3::new(0.0, 2.0, 10.0),
    @graph.Vec3::new(0.0, 0.0, 0.0),
    @graph.Vec3::unit_y(),
    proj,
    @graph.Viewport::full(64.0, 64.0),
  )
  let fr = cam.frustum()
  let counts = [0, 0, 0]
  for _ in 0..<2000 {
    let c = rng.vec3(-40.0, 40.0)
    let h = rng.vec3(0.1, 6.0)
    let b = @graph.Aabb::from_center_half(c, h)
    let cls = fr.classify_aabb(b)
    let mut any = false
    let mut all = true
    for i in 0..<8 {
      if fr.contains_point(b.corner(i)) {
        any = true
      } else {
        all = false
      }
    }
    match cls {
      Inside => {
        counts[2] += 1
        check_true(all, "Inside => every corner inside")
      }
      Outside => {
        counts[0] += 1
        check_true(!any, "Outside => no corner inside")
        // and no random point of the box is inside
        for _ in 0..<20 {
          let q = c + rng.vec3(-1.0, 1.0) * h
          check_true(!fr.contains_point(q), "Outside box has no visible point")
        }
      }
      Intersecting => counts[1] += 1
    }
    if any && !all {
      check_true(cls != Outside, "partially visible box is not Outside")
    }
  }
  check_true(
    counts[0] > 100 && counts[1] > 50 && counts[2] > 20,
    "all three classes exercised: \{counts[0]}/\{counts[1]}/\{counts[2]}",
  )
}

///|
test "camera-relative (floating origin) rendering removes float32 jitter at planetary scale" {
  // Camera 6,378,137 m from the origin (surface of the Earth), looking at geometry within 60 m.
  let eye = @graph.Vec3::new(6378137.0, 1234.5, -4321.25)
  let proj = @graph.Projection::perspective(
    1.0,
    1.0,
    0.1,
    @double.infinity,
    ReversedZ,
  )
  let cam = @graph.Camera::look_at(
    eye,
    eye + @graph.Vec3::new(-1.0, 0.1, 0.2),
    @graph.Vec3::unit_y(),
    proj,
    @graph.Viewport::full(1000.0, 1000.0),
  )
  let f32 = fn(x : Double) -> Double { Float::from_double(x).to_double() }
  let f32m = fn(m : @graph.Mat4) -> @graph.Mat4 {
    @graph.Mat4::from_col_major(m.to_col_major().map(f32))
  }
  let vp_naive = f32m(cam.view_proj())
  let vp_rte = f32m(cam.view_proj_rte())
  let rng = Rng::new(2UL)
  let mut worst_naive = 0.0
  let mut worst_rte = 0.0
  for _ in 0..<500 {
    // point 3..60 m in front of the camera
    let loc = cam.forward().scale(rng.range(3.0, 60.0)) +
      cam.right().scale(rng.range(-3.0, 3.0)) +
      cam.up().scale(rng.range(-3.0, 3.0))
    let world = eye + loc
    let truth = cam.world_to_screen(world).unwrap()
    // naive: world position stored as float32
    let pn = @graph.Vec3::new(f32(world.x), f32(world.y), f32(world.z))
    let sn = cam.viewport.ndc_to_screen(
      vp_naive.project_point(pn).unwrap(),
      ReversedZ,
    )
    // RTE: subtract the eye in double precision, THEN round to float32
    let rel = world - eye
    let pr = @graph.Vec3::new(f32(rel.x), f32(rel.y), f32(rel.z))
    let sr = cam.viewport.ndc_to_screen(
      vp_rte.project_point(pr).unwrap(),
      ReversedZ,
    )
    worst_naive = worst_naive.max(
      (sn.x - truth.x).abs().max((sn.y - truth.y).abs()),
    )
    worst_rte = worst_rte.max(
      (sr.x - truth.x).abs().max((sr.y - truth.y).abs()),
    )
  }
  check_true(
    worst_rte < 0.05,
    "RTE error in pixels should be tiny, got \{worst_rte}",
  )
  check_true(
    worst_naive > 20.0 * worst_rte && worst_naive > 1.0,
    "naive f32 world coords jitter visibly: naive=\{worst_naive}px rte=\{worst_rte}px",
  )
  // exactness in double: Mat4::relative_to composes with view_proj_rte to equal the full pipeline
  let model = @graph.Mat4::from_trs(
    eye + @graph.Vec3::new(5.0, -2.0, 7.0),
    @graph.Quat::from_y(0.3),
    @graph.Vec3::splat(2.0),
  )
  let p = @graph.Vec3::new(0.5, 0.25, -0.75)
  let full = cam.view_proj().project_point(model.mul_point(p)).unwrap()
  let rte = cam
    .view_proj_rte()
    .project_point(model.relative_to(eye).mul_point(p))
    .unwrap()
  check_v3(rte, full, 1.0e-9, "relative_to + view_proj_rte == full pipeline")
}
