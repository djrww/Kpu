// Axis-convention conversion (Blender / Unity / Unreal / ROS / glTF) tests.

///|
fn presets() -> Array[(String, @graph.Convention)] {
  [
    ("native", @graph.Convention::native()),
    ("blender", @graph.Convention::z_up_right_handed()),
    ("unity", @graph.Convention::y_up_left_handed()),
    ("unreal", @graph.Convention::unreal()),
    ("ros", @graph.Convention::ros()),
  ]
}

///|
test "handedness of the presets" {
  check_true(@graph.Convention::native().is_right_handed(), "native RH")
  check_true(
    @graph.Convention::z_up_right_handed().is_right_handed(),
    "blender RH",
  )
  check_true(
    !@graph.Convention::y_up_left_handed().is_right_handed(),
    "unity LH",
  )
  check_true(!@graph.Convention::unreal().is_right_handed(), "unreal LH")
  check_true(@graph.Convention::ros().is_right_handed(), "ros RH")
  for (_, c) in presets() {
    check_true(c.is_valid(), "valid")
  }
  check_true(
    !({ right: PosX, up: NegX, forward: PosZ, } : @graph.Convention).is_valid(),
    "invalid: same axis twice",
  )
}

///|
test "well-known conversions into the native frame" {
  let nat = @graph.Convention::native()
  // Blender (x, y, z) Z-up -> glTF/OpenGL (x, z, -y)   [the glTF exporter's rule]
  check_v3(
    @graph.Convention::z_up_right_handed().convert_vec3(
      nat,
      @graph.Vec3::new(1.0, 2.0, 3.0),
    ),
    @graph.Vec3::new(1.0, 3.0, -2.0),
    0.0,
    "blender",
  )
  // Unity (x, y, z) -> (x, y, -z)
  check_v3(
    @graph.Convention::y_up_left_handed().convert_vec3(
      nat,
      @graph.Vec3::new(1.0, 2.0, 3.0),
    ),
    @graph.Vec3::new(1.0, 2.0, -3.0),
    0.0,
    "unity",
  )
  // Unreal: forward +X, right +Y, up +Z -> native: (x, y, z) => (y, z, -x)
  check_v3(
    @graph.Convention::unreal().convert_vec3(
      nat,
      @graph.Vec3::new(1.0, 2.0, 3.0),
    ),
    @graph.Vec3::new(2.0, 3.0, -1.0),
    0.0,
    "unreal",
  )
  // ROS: forward +X, left +Y, up +Z -> native: (x, y, z) => (-y, z, -x)
  check_v3(
    @graph.Convention::ros().convert_vec3(nat, @graph.Vec3::new(1.0, 2.0, 3.0)),
    @graph.Vec3::new(-2.0, 3.0, -1.0),
    0.0,
    "ros",
  )
  // identity conversion
  check_m3(nat.basis_change(nat), @graph.Mat3::identity(), 0.0, "identity")
}

///|
test "basis_change: orthogonal, det = handedness flip, composes, invertible, semantic roles preserved" {
  let ps = presets()
  for (n1, a) in ps {
    for (n2, b) in ps {
      let m = a.basis_change(b)
      check_m3(
        m * m.transpose(),
        @graph.Mat3::identity(),
        0.0,
        "orthogonal \{n1}->\{n2}",
      )
      let det = m.determinant()
      let want = if a.is_right_handed() == b.is_right_handed() {
        1.0
      } else {
        -1.0
      }
      check_close(det, want, 0.0, "det \{n1}->\{n2}")
      check_m3(m, b.basis_change(a).transpose(), 0.0, "inverse is transpose")
      // roles: a's right axis maps to b's right axis, etc.
      check_v3(m.mul_vec3(a.right.to_vec3()), b.right.to_vec3(), 0.0, "right")
      check_v3(m.mul_vec3(a.up.to_vec3()), b.up.to_vec3(), 0.0, "up")
      check_v3(
        m.mul_vec3(a.forward.to_vec3()),
        b.forward.to_vec3(),
        0.0,
        "forward",
      )
      for (_, c) in ps {
        check_m3(b.basis_change(c) * m, a.basis_change(c), 0.0, "compose")
      }
    }
  }
}

///|
test "rotations and 4x4 transforms convert consistently (also across a mirror)" {
  let rng = Rng::new(23UL)
  let ps = presets()
  for (n1, a) in ps {
    for (n2, b) in ps {
      for i in 0..<20 {
        let q = rng.quat()
        let v = rng.vec3(-3.0, 3.0)
        let cq = a.convert_quat(b, q)
        check_close(cq.length(), 1.0, 1.0e-14, "unit")
        // rotating then converting == converting then rotating by the converted quaternion
        check_v3(
          cq.rotate(a.convert_vec3(b, v)),
          a.convert_vec3(b, q.rotate(v)),
          1.0e-12,
          "quat \{n1}->\{n2} #\{i}",
        )
        let t = @graph.Mat4::from_trs(
          rng.vec3(-9.0, 9.0),
          q,
          @graph.Vec3::new(1.5, 0.5, 2.0),
        )
        let ct = a.convert_mat4(b, t)
        check_v3(
          ct.mul_point(a.convert_vec3(b, v)),
          a.convert_vec3(b, t.mul_point(v)),
          1.0e-12,
          "mat4 \{n1}->\{n2} #\{i}",
        )
        // a pure rotation matrix converts to the matrix of the converted quaternion
        check_m3(
          a.convert_mat4(b, @graph.Mat4::rotation(q)).to_mat3(),
          cq.to_mat3(),
          1.0e-12,
          "rot mat vs quat",
        )
      }
    }
  }
}

///|
test "mirror conversion flips triangle winding (signed area of projected faces)" {
  let tri = (
    @graph.Vec3::new(0.0, 0.0, 0.0),
    @graph.Vec3::new(1.0, 0.0, 0.0),
    @graph.Vec3::new(0.0, 1.0, 0.0),
  )
  let (a, b, c) = tri
  let n_native = (b - a).cross(c - a)
  let u = @graph.Convention::y_up_left_handed()
  let nat = @graph.Convention::native()
  let (ua, ub, uc) = (
    u.convert_vec3(nat, a),
    u.convert_vec3(nat, b),
    u.convert_vec3(nat, c),
  )
  let n_conv = (ub - ua).cross(uc - ua)
  // conversion Unity -> native is a mirror: the face normal computed in the new frame is the
  // MIRROR image of the converted old normal => opposite direction relative to "carried" normal
  let carried = u.convert_vec3(nat, n_native)
  check_true(
    n_conv.dot(carried) < 0.0,
    "winding flips across a handedness change",
  )
  // same-handedness conversion keeps it
  let bl = @graph.Convention::z_up_right_handed()
  let (ba, bb, bc) = (
    bl.convert_vec3(nat, a),
    bl.convert_vec3(nat, b),
    bl.convert_vec3(nat, c),
  )
  check_true(
    (bb - ba).cross(bc - ba).dot(bl.convert_vec3(nat, n_native)) > 0.0,
    "winding kept",
  )
}
