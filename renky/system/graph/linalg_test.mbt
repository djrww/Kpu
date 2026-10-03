// Linear algebra checked against numpy / scipy reference values (golden_linalg_test.mbt)
// and against algebraic identities on random inputs.

///|
fn g_mat4(a : FixedArray[Double], off : Int) -> @graph.Mat4 {
  @graph.Mat4::from_cols(
    @graph.Vec4::new(a[off + 0], a[off + 4], a[off + 8], a[off + 12]),
    @graph.Vec4::new(a[off + 1], a[off + 5], a[off + 9], a[off + 13]),
    @graph.Vec4::new(a[off + 2], a[off + 6], a[off + 10], a[off + 14]),
    @graph.Vec4::new(a[off + 3], a[off + 7], a[off + 11], a[off + 15]),
  )
}

///|
fn g_mat3(a : FixedArray[Double], off : Int) -> @graph.Mat3 {
  @graph.Mat3::from_rows(
    @graph.Vec3::new(a[off], a[off + 1], a[off + 2]),
    @graph.Vec3::new(a[off + 3], a[off + 4], a[off + 5]),
    @graph.Vec3::new(a[off + 6], a[off + 7], a[off + 8]),
  )
}

///|
fn g_v3(a : FixedArray[Double], off : Int) -> @graph.Vec3 {
  @graph.Vec3::new(a[off], a[off + 1], a[off + 2])
}

///|
fn g_quat(a : FixedArray[Double], off : Int) -> @graph.Quat {
  @graph.Quat::new(a[off], a[off + 1], a[off + 2], a[off + 3])
}

///|
test "Mat4: product, determinant, inverse, mul_vec4 match numpy" {
  let g = gold_mat4
  let stride = 16 * 4 + 1 + 4 + 4
  for i in 0..<(g.length() / stride) {
    let o = i * stride
    let a = g_mat4(g, o)
    let b = g_mat4(g, o + 16)
    check_m4(a * b, g_mat4(g, o + 32), 1.0e-12, "A*B #\{i}")
    check_close(a.determinant(), g[o + 64], 1.0e-11, "det #\{i}")
    match a.inverse() {
      Some(inv) => check_m4(inv, g_mat4(g, o + 48), 1.0e-10, "inv #\{i}")
      None => fail("inverse returned None for invertible matrix #\{i}")
    }
    let v = @graph.Vec4::new(g[o + 65], g[o + 66], g[o + 67], g[o + 68])
    let want = @graph.Vec4::new(g[o + 69], g[o + 70], g[o + 71], g[o + 72])
    check_v4(a.mul_vec4(v), want, 1.0e-12, "A*v #\{i}")
  }
}

///|
test "Mat4: inverse of singular matrix is None" {
  // two identical rows -> det == 0 exactly
  let m : @graph.Mat4 = {
    m00: 1.0,
    m01: 2.0,
    m02: 3.0,
    m03: 4.0,
    m10: 1.0,
    m11: 2.0,
    m12: 3.0,
    m13: 4.0,
    m20: 0.0,
    m21: 1.0,
    m22: 0.0,
    m23: 1.0,
    m30: 5.0,
    m31: 0.0,
    m32: 0.0,
    m33: 1.0,
  }
  check_true(m.inverse() is None, "singular -> None")
  check_true(
    @graph.Mat3::diagonal(1.0, 0.0, 1.0).inverse() is None,
    "Mat3 singular -> None",
  )
}

///|
test "Mat4: TRS matrix, inverse paths, normal matrix, decompose vs numpy" {
  let g = gold_trs
  let stride = 3 + 4 + 3 + 16 + 16 + 9 + 3 + 3
  for i in 0..<(g.length() / stride) {
    let o = i * stride
    let t = g_v3(g, o)
    let q = g_quat(g, o + 3)
    let s = g_v3(g, o + 7)
    let m = @graph.Mat4::from_trs(t, q, s)
    check_m4(m, g_mat4(g, o + 10), 1.0e-12, "from_trs #\{i}")
    // affine_inverse and general inverse both match numpy
    match m.affine_inverse() {
      Some(inv) =>
        check_m4(inv, g_mat4(g, o + 26), 1.0e-11, "affine_inverse #\{i}")
      None => fail("affine_inverse None #\{i}")
    }
    match m.inverse() {
      Some(inv) => check_m4(inv, g_mat4(g, o + 26), 1.0e-11, "inverse #\{i}")
      None => fail("inverse None #\{i}")
    }
    match m.normal_matrix() {
      Some(n) => check_m3(n, g_mat3(g, o + 42), 1.0e-11, "normal matrix #\{i}")
      None => fail("normal None #\{i}")
    }
    let p = g_v3(g, o + 51)
    check_v3(m.mul_point(p), g_v3(g, o + 54), 1.0e-11, "mul_point #\{i}")
    // decompose -> recompose gives the same matrix (sign ambiguity of mirrored scale is fine)
    match m.decompose() {
      Some(tr) =>
        check_m4(tr.to_mat4(), m, 1.0e-11, "decompose roundtrip #\{i}")
      None => fail("decompose None #\{i}")
    }
  }
}

///|
test "Mat4::decompose rejects shear, projective and zero scale" {
  let shear = @graph.Mat4::identity()
  let shear = { ..shear, m01: 0.5, }
  check_true(shear.decompose() is None, "shear rejected")
  let p = @graph.Projection::perspective(1.0, 1.5, 0.1, 100.0, ZeroToOne).to_mat4()
  check_true(p.decompose() is None, "projective rejected")
  check_true(
    @graph.Mat4::scaling(@graph.Vec3::new(1.0, 0.0, 1.0)).decompose() is None,
    "zero scale rejected",
  )
}

///|
test "rigid_inverse == general inverse for rigid matrices; affine_inverse for TRS" {
  let rng = Rng::new(7UL)
  for i in 0..<200 {
    let rg = @graph.Mat4::from_trs(
      rng.vec3(-1000.0, 1000.0),
      rng.quat(),
      @graph.Vec3::splat(1.0),
    )
    let a = rg.rigid_inverse()
    match rg.inverse() {
      Some(b) => check_m4(a, b, 1.0e-12, "rigid vs general #\{i}")
      None => fail("none")
    }
    check_m4(rg * a, @graph.Mat4::identity(), 1.0e-12, "rigid * inv #\{i}")
  }
}

///|
test "Mat4: mul_point/mul_dir/project_point/transform_points agree with mul_vec4" {
  let rng = Rng::new(11UL)
  let m = @graph.Mat4::from_trs(
    rng.vec3(-9.0, 9.0),
    rng.quat(),
    @graph.Vec3::new(1.5, 0.7, 2.0),
  )
  let pr = @graph.Projection::perspective(1.0, 1.3, 0.5, 80.0, ZeroToOne).to_mat4() *
    m
  let n = 64
  let src : FixedArray[Double] = FixedArray::make(n * 3, 0.0)
  let pts : Array[@graph.Vec3] = []
  for i in 0..<n {
    let p = rng.vec3(-10.0, 10.0)
    pts.push(p)
    src[i * 3] = p.x
    src[i * 3 + 1] = p.y
    src[i * 3 + 2] = p.z
  }
  let dst : FixedArray[Double] = FixedArray::make(n * 3, 0.0)
  m.transform_points(src, dst, n)
  for i in 0..<n {
    let p = pts[i]
    let h = m.mul_vec4(p.to_vec4(1.0))
    check_v3(m.mul_point(p), h.xyz(), 1.0e-12, "mul_point")
    check_v3(g_v3(dst, i * 3), h.xyz(), 1.0e-12, "transform_points")
    check_close(h.w, 1.0, 0.0, "affine keeps w")
    let d = m.mul_vec4(p.to_vec4(0.0))
    check_v3(m.mul_dir(p), d.xyz(), 1.0e-12, "mul_dir")
    let c = pr.mul_vec4(p.to_vec4(1.0))
    match pr.project_point(p) {
      Some(q) => check_v3(q, c.xyz().scale(1.0 / c.w), 1.0e-12, "project_point")
      None => fail("project_point none")
    }
  }
  // in-place aliasing
  m.transform_points(src, src, n)
  for i in 0..<(n * 3) {
    check_close(src[i], dst[i], 1.0e-13, "in-place transform_points")
  }
}

///|
test "Mat4 column-major packing is GPU order and round-trips" {
  let g = gold_trs
  let m = g_mat4(g, 10)
  let cm = m.to_col_major()
  // column-major: element k = col (k/4), row (k%4)
  for k in 0..<16 {
    check_close(cm[k], m.get(k % 4, k / 4), 0.0, "cm[\{k}]")
  }
  check_m4(@graph.Mat4::from_col_major(cm), m, 0.0, "from_col_major")
  let f = m.to_col_major_f32()
  for k in 0..<16 {
    check_close(f[k].to_double(), cm[k], 1.0e-6, "f32[\{k}]")
  }
  // translation lives in elements 12..14 (GLSL convention)
  check_close(cm[12], m.m03, 0.0, "tx slot")
  check_close(cm[13], m.m13, 0.0, "ty slot")
  check_close(cm[14], m.m23, 0.0, "tz slot")
}

///|
test "Quat: product, rotate, slerp, matrix, rotvec match scipy" {
  let g = gold_quat
  let stride = 4 + 4 + 3 + 4 + 3 + 1 + 4 + 3 + 9
  for i in 0..<(g.length() / stride) {
    let o = i * stride
    let q1 = g_quat(g, o)
    let q2 = g_quat(g, o + 4)
    let v = g_v3(g, o + 8)
    check_quat(q1 * q2, g_quat(g, o + 11), 1.0e-13, "q1*q2 #\{i}")
    check_v3(q1.rotate(v), g_v3(g, o + 15), 1.0e-12, "rotate #\{i}")
    check_quat(
      q1.slerp(q2, g[o + 18]),
      g_quat(g, o + 19),
      1.0e-12,
      "slerp #\{i}",
    )
    check_v3(q1.to_rotation_vector(), g_v3(g, o + 23), 1.0e-12, "rotvec #\{i}")
    check_quat(
      @graph.Quat::from_rotation_vector(g_v3(g, o + 23)),
      q1,
      1.0e-12,
      "from rotvec #\{i}",
    )
    let m = g_mat3(g, o + 26)
    check_m3(q1.to_mat3(), m, 1.0e-13, "to_mat3 #\{i}")
    check_quat(@graph.Quat::from_mat3(m), q1, 1.0e-12, "from_mat3 #\{i}")
    // rotate via matrix == rotate via quat
    check_v3(
      m.mul_vec3(v),
      q1.rotate(v),
      1.0e-12,
      "matrix vs quat rotate #\{i}",
    )
  }
}

///|
test "Quat::from_mat3 covers all four Shepperd branches (180-degree flips)" {
  let axes = [
    @graph.Vec3::unit_x(),
    @graph.Vec3::unit_y(),
    @graph.Vec3::unit_z(),
    @graph.Vec3::new(1.0, 1.0, 0.0).normalize(),
    @graph.Vec3::new(0.0, 1.0, 1.0).normalize(),
  ]
  for a in axes {
    for ang in [0.0, 1.0, 3.0, @graph.PI, -@graph.PI, 3.14, 2.0] {
      let q = @graph.Quat::from_axis_angle(a, ang)
      let m = q.to_mat3()
      let q2 = @graph.Quat::from_mat3(m)
      check_m3(q2.to_mat3(), m, 1.0e-13, "branch roundtrip angle \{ang}")
      check_close(q2.length(), 1.0, 1.0e-14, "unit")
    }
  }
}

///|
test "Quat: axis-angle, from_to, look_rotation, integrate" {
  let rng = Rng::new(99UL)
  for i in 0..<200 {
    let axis = rng.unit_vec3()
    let ang = rng.range(0.0, @graph.PI)
    let q = @graph.Quat::from_axis_angle(axis, ang)
    let (a2, ang2) = q.to_axis_angle()
    check_close(ang2, ang, 1.0e-12, "angle #\{i}")
    check_v3(a2, axis, 1.0e-9, "axis #\{i}")
    // from_to maps a onto b
    let a = rng.unit_vec3()
    let b = rng.unit_vec3().scale(rng.range(0.1, 9.0))
    check_v3(
      @graph.Quat::from_to(a, b).rotate(a),
      b.normalize(),
      1.0e-12,
      "from_to #\{i}",
    )
    // look_rotation: local -Z -> forward, local +Y has non-negative dot with up
    let f = rng.unit_vec3()
    let up = rng.unit_vec3()
    let lr = @graph.Quat::look_rotation(f, up)
    check_v3(
      lr.rotate(@graph.Vec3::new(0.0, 0.0, -1.0)),
      f,
      1.0e-12,
      "look fwd #\{i}",
    )
    check_true(
      lr.rotate(@graph.Vec3::unit_y()).dot(up) >= -1.0e-12,
      "look up side #\{i}",
    )
    check_close(
      lr.rotate(@graph.Vec3::unit_x()).dot(f),
      0.0,
      1.0e-12,
      "look right perpendicular #\{i}",
    )
  }
  // antiparallel
  let q = @graph.Quat::from_to(
    @graph.Vec3::unit_x(),
    @graph.Vec3::new(-1.0, 0.0, 0.0),
  )
  check_v3(
    q.rotate(@graph.Vec3::unit_x()),
    @graph.Vec3::new(-1.0, 0.0, 0.0),
    1.0e-12,
    "antiparallel",
  )
  check_close(q.length(), 1.0, 1.0e-14, "antiparallel unit")
  // degenerate up (parallel to forward) must still give a valid rotation
  let d = @graph.Quat::look_rotation(
    @graph.Vec3::unit_y(),
    @graph.Vec3::unit_y(),
  )
  check_close(d.length(), 1.0, 1.0e-14, "degenerate unit")
  check_v3(
    d.rotate(@graph.Vec3::new(0.0, 0.0, -1.0)),
    @graph.Vec3::unit_y(),
    1.0e-12,
    "degenerate fwd",
  )
  // integrate: constant omega about y for 1s in 1000 steps == one exact rotation
  let mut o = @graph.Quat::identity()
  for _ in 0..<1000 {
    o = o.integrate(@graph.Vec3::new(0.0, 1.5, 0.0), 0.001)
  }
  check_quat(o, @graph.Quat::from_y(1.5), 1.0e-12, "integrate")
}

///|
test "Quat: slerp has constant angular speed and hits endpoints; shortest path" {
  let rng = Rng::new(5UL)
  for i in 0..<100 {
    let a = rng.quat()
    let b = rng.quat()
    check_quat(a.slerp(b, 0.0), a, 1.0e-12, "t=0")
    check_quat(a.slerp(b, 1.0), b, 1.0e-12, "t=1")
    let total = a.angle_to(b)
    check_true(total <= @graph.PI + 1.0e-12, "shortest")
    for t in [0.25, 0.5, 0.9] {
      check_close(
        a.angle_to(a.slerp(b, t)),
        total * t,
        1.0e-9,
        "uniform speed #\{i}",
      )
    }
  }
  // tiny rotation uses the nlerp branch but stays unit
  let a = @graph.Quat::identity()
  let b = @graph.Quat::from_x(1.0e-6)
  check_close(a.slerp(b, 0.5).length(), 1.0, 1.0e-15, "tiny unit")
  // q and -q are the same orientation: slerp must not travel the long way
  let c = @graph.Quat::from_y(0.3)
  check_quat(
    a.slerp(c.negate(), 0.5),
    @graph.Quat::from_y(0.15),
    1.0e-12,
    "neg shortest",
  )
}

///|
test "Euler: all 12 orders, intrinsic convention == scipy upper-case" {
  let g = gold_euler
  let orders = [
    @graph.EulerOrder::XYZ,
    XZY,
    YXZ,
    YZX,
    ZXY,
    ZYX,
    XYX,
    XZX,
    YXY,
    YZY,
    ZXZ,
    ZYZ,
  ]
  let stride = 3 + 9 + 3
  for oi in 0..<12 {
    let order = orders[oi]
    for k in 0..<8 {
      let o = (oi * 8 + k) * stride
      let a = g[o]
      let b = g[o + 1]
      let c = g[o + 2]
      let want = g_mat3(g, o + 3)
      let tag = "\{oi}/\{k}"
      check_m3(
        @graph.Mat3::from_euler(order, a, b, c),
        want,
        1.0e-13,
        "from_euler " + tag,
      )
      check_m3(
        @graph.Quat::from_euler(order, a, b, c).to_mat3(),
        want,
        1.0e-13,
        "quat from_euler " + tag,
      )
      // to_euler: reconstructs the matrix AND yields scipy's angles for non-singular input
      let (ra, rb, rc) = want.to_euler(order)
      check_m3(
        @graph.Mat3::from_euler(order, ra, rb, rc),
        want,
        1.0e-12,
        "to_euler roundtrip " + tag,
      )
      check_close(ra, g[o + 12], 1.0e-10, "angle a vs scipy " + tag)
      check_close(rb, g[o + 13], 1.0e-10, "angle b vs scipy " + tag)
      check_close(rc, g[o + 14], 1.0e-10, "angle c vs scipy " + tag)
    }
  }
}

///|
test "Euler: extrinsic sequence XYZ(a,b,c) == intrinsic ZYX(c,b,a)" {
  let g = gold_euler_ext
  let orders = [
    @graph.EulerOrder::XYZ,
    XZY,
    YXZ,
    YZX,
    ZXY,
    ZYX,
    XYX,
    XZX,
    YXY,
    YZY,
    ZXZ,
    ZYZ,
  ]
  let rev = [
    @graph.EulerOrder::ZYX,
    YZX,
    ZXY,
    XZY,
    YXZ,
    XYZ,
    XYX,
    XZX,
    YXY,
    YZY,
    ZXZ,
    ZYZ,
  ]
  let _ = orders
  for oi in 0..<12 {
    for k in 0..<4 {
      let o = (oi * 4 + k) * 12
      let a = g[o]
      let b = g[o + 1]
      let c = g[o + 2]
      check_m3(
        @graph.Mat3::from_euler(rev[oi], c, b, a),
        g_mat3(g, o + 3),
        1.0e-13,
        "extrinsic \{oi}/\{k}",
      )
    }
  }
}

///|
test "Euler: gimbal lock still round-trips the matrix" {
  let orders = [
    @graph.EulerOrder::XYZ,
    XZY,
    YXZ,
    YZX,
    ZXY,
    ZYX,
    XYX,
    XZX,
    YXY,
    YZY,
    ZXZ,
    ZYZ,
  ]
  for order in orders {
    let bs = if order.is_proper() {
      [0.0, @graph.PI]
    } else {
      [@graph.HALF_PI, -@graph.HALF_PI]
    }
    for b in bs {
      for (a, c) in [(0.3, 0.4), (-2.0, 1.1), (0.0, 0.0)] {
        let m = @graph.Mat3::from_euler(order, a, b, c)
        let (ra, rb, rc) = m.to_euler(order)
        check_m3(
          @graph.Mat3::from_euler(order, ra, rb, rc),
          m,
          1.0e-9,
          "gimbal b=\{b}",
        )
      }
    }
  }
}

///|
test "Mat4::look_at equals gluLookAt (numpy) and Rigid::look_at inverse" {
  let g = gold_lookat
  let stride = 9 + 16
  for i in 0..<(g.length() / stride) {
    let o = i * stride
    let eye = g_v3(g, o)
    let tgt = g_v3(g, o + 3)
    let up = g_v3(g, o + 6)
    let v = @graph.Mat4::look_at(eye, tgt, up)
    check_m4(v, g_mat4(g, o + 9), 1.0e-12, "look_at #\{i}")
    // eye maps to the origin, target maps onto the -Z axis
    check_v3(v.mul_point(eye), @graph.Vec3::zero(), 1.0e-12, "eye->origin")
    let t = v.mul_point(tgt)
    check_close(t.x, 0.0, 1.0e-11, "target x")
    check_close(t.y, 0.0, 1.0e-11, "target y")
    check_true(t.z < 0.0, "target in front (-Z)")
    let pose = @graph.Rigid::look_at(eye, tgt, up)
    check_m4(pose.to_mat4().rigid_inverse(), v, 1.0e-12, "pose vs view")
  }
}

///|
test "Rigid: group laws" {
  let rng = Rng::new(3UL)
  for i in 0..<100 {
    let a = @graph.Rigid::new(rng.quat(), rng.vec3(-100.0, 100.0))
    let b = @graph.Rigid::new(rng.quat(), rng.vec3(-100.0, 100.0))
    let p = rng.vec3(-5.0, 5.0)
    check_v3(
      (a * b).mul_point(p),
      a.mul_point(b.mul_point(p)),
      1.0e-12,
      "compose #\{i}",
    )
    check_v3((a * a.inverse()).mul_point(p), p, 1.0e-12, "inverse #\{i}")
    check_m4(
      (a * b).to_mat4(),
      a.to_mat4() * b.to_mat4(),
      1.0e-12,
      "matrix hom #\{i}",
    )
    check_m4(
      a.inverse().to_mat4(),
      a.to_mat4().rigid_inverse(),
      1.0e-12,
      "inverse mat #\{i}",
    )
    check_v3(a.lerp(b, 0.0).mul_point(p), a.mul_point(p), 1.0e-12, "lerp 0")
    check_v3(a.lerp(b, 1.0).mul_point(p), b.mul_point(p), 1.0e-12, "lerp 1")
  }
}

///|
test "Transform (TRS) point mapping equals its matrix" {
  let rng = Rng::new(21UL)
  for i in 0..<100 {
    let t = @graph.Transform::new(
      rng.vec3(-9.0, 9.0),
      rng.quat(),
      @graph.Vec3::new(
        rng.range(0.2, 3.0),
        -rng.range(0.2, 3.0),
        rng.range(0.2, 3.0),
      ),
    )
    let p = rng.vec3(-4.0, 4.0)
    check_v3(t.mul_point(p), t.to_mat4().mul_point(p), 1.0e-12, "point #\{i}")
    check_v3(t.mul_dir(p), t.to_mat4().mul_dir(p), 1.0e-12, "dir #\{i}")
  }
}

///|
test "Mat3 2D affine and homography" {
  let m = @graph.Mat3::from_trs2d(
    @graph.Vec2::new(10.0, -3.0),
    @graph.HALF_PI,
    @graph.Vec2::new(2.0, 3.0),
  )
  // scale first (2,3), rotate 90deg ccw -> (-3*1... ) then translate
  let p = m.mul_point2(@graph.Vec2::new(1.0, 1.0))
  check_v2(p, @graph.Vec2::new(10.0 - 3.0, -3.0 + 2.0), 1.0e-14, "trs2d")
  match m.inverse() {
    Some(i) =>
      check_v2(i.mul_point2(p), @graph.Vec2::new(1.0, 1.0), 1.0e-13, "inv2d")
    None => fail("none")
  }
  // pure scale+translate homography: w = 1
  check_true(
    @graph.Mat3::identity().project_point2(@graph.Vec2::new(3.0, 4.0))
    is Some(_),
    "id project",
  )
  // projective: x' = x / (1 + x)
  let h : @graph.Mat3 = { ..@graph.Mat3::identity(), m20: 1.0, }
  match h.project_point2(@graph.Vec2::new(3.0, 4.0)) {
    Some(q) => check_v2(q, @graph.Vec2::new(0.75, 1.0), 1.0e-15, "homography")
    None => fail("none")
  }
  check_true(
    h.project_point2(@graph.Vec2::new(-1.0, 0.0)) is None,
    "w==0 -> None",
  )
}

///|
test "Vec helpers" {
  let a = @graph.Vec3::new(1.0, 2.0, 3.0)
  let b = @graph.Vec3::new(-4.0, 0.5, 2.0)
  check_close(a.dot(b), -4.0 + 1.0 + 6.0, 0.0, "dot")
  check_v3(
    a.cross(b),
    @graph.Vec3::new(
      2.0 * 2.0 - 3.0 * 0.5,
      3.0 * -4.0 - 1.0 * 2.0,
      1.0 * 0.5 - 2.0 * -4.0,
    ),
    0.0,
    "cross",
  )
  check_true(@graph.Vec3::zero().try_normalize() is None, "zero normalize")
  check_close(a.normalize().length(), 1.0, 1.0e-15, "normalize")
  check_v3(
    @graph.Vec3::new(1.0, -1.0, 0.0).reflect(@graph.Vec3::unit_y()),
    @graph.Vec3::new(1.0, 1.0, 0.0),
    0.0,
    "reflect",
  )
  check_close(
    @graph.Vec3::unit_x().angle_to(@graph.Vec3::unit_y()),
    @graph.HALF_PI,
    1.0e-15,
    "angle",
  )
  check_close(
    @graph.Vec3::unit_x().angle_to(@graph.Vec3::new(1.0, 1.0e-9, 0.0)),
    1.0e-9,
    1.0e-18,
    "tiny angle precision",
  )
  let rng = Rng::new(1UL)
  for _ in 0..<100 {
    let v = rng.vec3(-3.0, 3.0)
    let p = v.any_perpendicular()
    check_close(p.dot(v), 0.0, 1.0e-12, "perp")
    check_close(p.length(), 1.0, 1.0e-12, "perp unit")
  }
  check_v2(
    @graph.Vec2::new(1.0, 0.0).rotate(@graph.HALF_PI),
    @graph.Vec2::new(0.0, 1.0),
    1.0e-15,
    "rot2",
  )
  check_close(
    @graph.Vec2::new(1.0, 0.0).cross(@graph.Vec2::new(0.0, 1.0)),
    1.0,
    0.0,
    "cross2",
  )
  check_close(@graph.Vec2::new(3.0, 4.0).length(), 5.0, 0.0, "len2")
  check_true(
    @graph.Vec4::new(1.0, 2.0, 3.0, 0.0).perspective_divide() is None,
    "w=0",
  )
  match @graph.Vec4::new(2.0, 4.0, 6.0, 2.0).perspective_divide() {
    Some(v) => check_v3(v, @graph.Vec3::new(1.0, 2.0, 3.0), 0.0, "divide")
    None => fail("none")
  }
}

///|
test "scalar: wrap_pi / wrap_tau ranges" {
  for k in -20..<21 {
    let a = Double::from_int(k) * 0.7
    let w = @graph.wrap_pi(a)
    check_true(w > -@graph.PI && w <= @graph.PI, "wrap_pi range \{a}")
    check_close(@math.cos(w), @math.cos(a), 1.0e-12, "wrap_pi same angle")
    check_close(@math.sin(w), @math.sin(a), 1.0e-12, "wrap_pi same angle")
    let t = @graph.wrap_tau(a)
    check_true(t >= 0.0 && t < @graph.TAU, "wrap_tau range \{a}")
  }
  check_close(@graph.wrap_pi(-@graph.PI), @graph.PI, 0.0, "-PI -> PI")
  check_close(@graph.wrap_pi(3.0 * @graph.PI), @graph.PI, 1.0e-12, "3PI -> PI")
  check_true(
    @graph.mod_floor(-1, 4) == 3 && @graph.mod_floor(9, 4) == 1,
    "mod_floor",
  )
}
