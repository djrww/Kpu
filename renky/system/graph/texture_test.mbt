// UV / texel / wrap / bilinear / mip / cube / equirect / octahedral mapping tests.

///|
test "uv <-> texel, centres, flip_v" {
  check_v2(
    @graph.uv_to_texel(@graph.Vec2::new(0.5, 0.25), 8, 4),
    @graph.Vec2::new(4.0, 1.0),
    0.0,
    "uv->texel",
  )
  check_v2(
    @graph.texel_to_uv(@graph.Vec2::new(4.0, 1.0), 8, 4),
    @graph.Vec2::new(0.5, 0.25),
    0.0,
    "texel->uv",
  )
  check_v2(
    @graph.texel_center_uv(0, 0, 4, 4),
    @graph.Vec2::new(0.125, 0.125),
    0.0,
    "centre uv",
  )
  check_v2(
    @graph.texel_center_uv(3, 1, 4, 2),
    @graph.Vec2::new(0.875, 0.75),
    0.0,
    "centre uv 2",
  )
  check_v2(
    @graph.flip_v(@graph.Vec2::new(0.2, 0.9)),
    @graph.Vec2::new(0.2, 1.0 - 0.9),
    1.0e-16,
    "flip_v",
  )
}

///|
test "wrap modes: coordinates and texel indices" {
  let cases : Array[(Double, Double, Double, Double)] = [
    // u, clamp, repeat, mirror
    (0.25, 0.25, 0.25, 0.25),
    (1.25, 1.0, 0.25, 0.75),
    (-0.25, 0.0, 0.75, 0.25),
    (2.25, 1.0, 0.25, 0.25),
    (-1.25, 0.0, 0.75, 0.75),
    (1.0, 1.0, 0.0, 1.0),
    (3.0, 1.0, 0.0, 1.0),
  ]
  for (u, c, r, m) in cases {
    check_close(@graph.wrap_coord(u, Clamp), c, 1.0e-15, "clamp \{u}")
    check_close(@graph.wrap_coord(u, Repeat), r, 1.0e-15, "repeat \{u}")
    check_close(@graph.wrap_coord(u, MirroredRepeat), m, 1.0e-15, "mirror \{u}")
  }
  // texel indices, n = 4
  let want_clamp = [0, 0, 0, 1, 2, 3, 3, 3]
  let want_repeat = [2, 3, 0, 1, 2, 3, 0, 1]
  let want_mirror = [1, 0, 0, 1, 2, 3, 3, 2]
  for k in 0..<8 {
    let i = k - 2
    check_true(@graph.wrap_texel(i, 4, Clamp) == want_clamp[k], "clamp i=\{i}")
    check_true(
      @graph.wrap_texel(i, 4, Repeat) == want_repeat[k],
      "repeat i=\{i}",
    )
    check_true(
      @graph.wrap_texel(i, 4, MirroredRepeat) == want_mirror[k],
      "mirror i=\{i}",
    )
  }
  // large negative indices stay in range
  for i in -1000..<1000 {
    for mode in [@graph.WrapMode::Clamp, Repeat, MirroredRepeat] {
      let w = @graph.wrap_texel(i, 7, mode)
      check_true(w >= 0 && w < 7, "range")
    }
  }
}

///|
test "bilinear footprint reproduces a linear image exactly (texel centres at (i+0.5)/w)" {
  let w = 16
  let h = 8
  let a = 3.0
  let b = -2.0
  let c = 5.0
  let tex = fn(i : Int, j : Int) -> Double {
    a * Double::from_int(i) + b * Double::from_int(j) + c
  }
  let rng = Rng::new(10UL)
  for _ in 0..<500 {
    // interior only (away from the clamp region where the linear ramp is truncated)
    let u = rng.range(
      0.5 / Double::from_int(w),
      1.0 - 0.5 / Double::from_int(w),
    )
    let v = rng.range(
      0.5 / Double::from_int(h),
      1.0 - 0.5 / Double::from_int(h),
    )
    let f = @graph.bilinear_footprint({ x: u, y: v, }, w, h, Clamp, Clamp)
    let s = (1.0 - f.fx) * (1.0 - f.fy) * tex(f.x0, f.y0) +
      f.fx * (1.0 - f.fy) * tex(f.x1, f.y0) +
      (1.0 - f.fx) * f.fy * tex(f.x0, f.y1) +
      f.fx * f.fy * tex(f.x1, f.y1)
    let x = u * Double::from_int(w) - 0.5
    let y = v * Double::from_int(h) - 0.5
    check_close(s, a * x + b * y + c, 1.0e-12, "linear reproduction")
    check_true(
      f.fx >= 0.0 && f.fx < 1.0 && f.fy >= 0.0 && f.fy < 1.0,
      "weights",
    )
  }
  // exactly on a texel centre: all weight on that texel
  let f = @graph.bilinear_footprint(
    @graph.texel_center_uv(5, 3, w, h),
    w,
    h,
    Clamp,
    Clamp,
  )
  check_true(
    f.x0 == 5 && f.y0 == 3 && f.fx == 0.0 && f.fy == 0.0,
    "centre sample",
  )
  // the edge: repeat wraps to the opposite side, clamp sticks
  let e0 = @graph.bilinear_footprint({ x: 0.0, y: 0.5, }, 4, 4, Repeat, Clamp)
  check_true(e0.x0 == 3 && e0.x1 == 0 && e0.fx == 0.5, "repeat edge")
  let e1 = @graph.bilinear_footprint({ x: 0.0, y: 0.5, }, 4, 4, Clamp, Clamp)
  check_true(e1.x0 == 0 && e1.x1 == 0, "clamp edge")
  let e2 = @graph.bilinear_footprint(
    { x: 0.0, y: 0.5, },
    4,
    4,
    MirroredRepeat,
    Clamp,
  )
  check_true(e2.x0 == 0 && e2.x1 == 0, "mirror edge")
}

///|
test "mip chain sizes and LOD from derivatives" {
  check_true(@graph.mip_level_count(1, 1) == 1, "1x1")
  check_true(@graph.mip_level_count(256, 64) == 9, "256x64")
  check_true(@graph.mip_level_count(5, 3) == 3, "5x3")
  check_true(@graph.mip_level_count(1024, 1024) == 11, "1024")
  check_true(@graph.mip_extent(256, 3) == 32, "extent")
  check_true(
    @graph.mip_extent(5, 2) == 1 && @graph.mip_extent(5, 9) == 1,
    "min 1",
  )
  // one texel per pixel -> lod 0; two texels per pixel -> lod 1; half -> lod -1
  check_close(
    @graph.mip_lod(
      @graph.Vec2::new(1.0 / 256.0, 0.0),
      @graph.Vec2::new(0.0, 1.0 / 256.0),
      256,
      256,
    ),
    0.0,
    1.0e-12,
    "lod0",
  )
  check_close(
    @graph.mip_lod(
      @graph.Vec2::new(2.0 / 256.0, 0.0),
      @graph.Vec2::new(0.0, 1.0 / 256.0),
      256,
      256,
    ),
    1.0,
    1.0e-12,
    "lod1",
  )
  check_close(
    @graph.mip_lod(
      @graph.Vec2::new(0.5 / 256.0, 0.0),
      @graph.Vec2::new(0.0, 0.25 / 256.0),
      256,
      256,
    ),
    -1.0,
    1.0e-12,
    "lod-1",
  )
  // non-square: derivative scaled by each axis size
  check_close(
    @graph.mip_lod(
      @graph.Vec2::new(1.0 / 64.0, 0.0),
      @graph.Vec2::new(0.0, 1.0 / 16.0),
      64,
      16,
    ),
    0.0,
    1.0e-12,
    "non-square",
  )
  check_true(
    @graph.mip_lod(@graph.Vec2::zero(), @graph.Vec2::zero(), 4, 4) < -1.0e300,
    "zero derivative",
  )
}

///|
test "cube map: OpenGL face table, hand-checked values" {
  let c = @graph.cube_from_direction(@graph.Vec3::new(1.0, 0.0, 0.0))
  check_true(c.face == @graph.CUBE_POS_X, "+X face")
  check_close(c.u, 0.5, 0.0, "+X u")
  check_close(c.v, 0.5, 0.0, "+X v")
  // +X: sc=-z, tc=-y: dir (1, 0.5, -0.5) -> s = (0.5+1)/2, t = (-0.5+1)/2
  let d = @graph.cube_from_direction(@graph.Vec3::new(1.0, 0.5, -0.5))
  check_true(d.face == 0, "face")
  check_close(d.u, 0.75, 1.0e-15, "u")
  check_close(d.v, 0.25, 1.0e-15, "v")
  // each axis -> its own face at the centre
  let axes = [
    @graph.Vec3::unit_x(),
    -@graph.Vec3::unit_x(),
    @graph.Vec3::unit_y(),
    -@graph.Vec3::unit_y(),
    @graph.Vec3::unit_z(),
    -@graph.Vec3::unit_z(),
  ]
  for f in 0..<6 {
    let cc = @graph.cube_from_direction(axes[f])
    check_true(cc.face == f, "axis face \{f}")
    check_close(cc.u, 0.5, 0.0, "centre u")
    check_close(cc.v, 0.5, 0.0, "centre v")
    check_v3(
      @graph.cube_to_direction(f, 0.5, 0.5),
      axes[f],
      1.0e-15,
      "centre dir",
    )
  }
  // +Y face (looking up): u -> +x, v -> +z (down in the image = towards +z)
  check_v3(
    @graph.cube_to_direction(2, 1.0, 0.5),
    @graph.Vec3::new(1.0, 1.0, 0.0).normalize(),
    1.0e-15,
    "+Y u",
  )
  check_v3(
    @graph.cube_to_direction(2, 0.5, 1.0),
    @graph.Vec3::new(0.0, 1.0, 1.0).normalize(),
    1.0e-15,
    "+Y v",
  )
  // side faces: v down means -y
  check_v3(
    @graph.cube_to_direction(4, 0.5, 1.0),
    @graph.Vec3::new(0.0, -1.0, 1.0).normalize(),
    1.0e-15,
    "+Z v down is -y",
  )
  // the zero vector doesn't crash
  check_true(@graph.cube_from_direction(@graph.Vec3::zero()).face == 0, "zero")
}

///|
test "cube map: round trip on random directions, ranges, edge continuity" {
  let rng = Rng::new(88UL)
  let seen = [0, 0, 0, 0, 0, 0]
  for i in 0..<3000 {
    let d = rng.unit_vec3().scale(rng.range(0.01, 100.0))
    let c = @graph.cube_from_direction(d)
    check_true(
      c.u >= 0.0 && c.u <= 1.0 && c.v >= 0.0 && c.v <= 1.0,
      "range #\{i}",
    )
    seen[c.face] += 1
    // the face is the one of the major axis
    let a = d.abs()
    let major = if a.x >= a.y && a.x >= a.z {
      0
    } else if a.y >= a.z {
      1
    } else {
      2
    }
    check_true(c.face / 2 == major, "major axis #\{i}")
    check_v3(
      @graph.cube_to_direction(c.face, c.u, c.v),
      d.normalize(),
      1.0e-12,
      "round trip #\{i}",
    )
  }
  for f in 0..<6 {
    check_true(seen[f] > 300, "face \{f} covered")
  }
  // continuity across every edge: two directions straddling a cube edge land on the
  // matching border of the two faces with the same position along the edge.
  let eps = 1.0e-7
  let a = @graph.cube_from_direction(@graph.Vec3::new(1.0, 0.2, 1.0 - eps)) // +X, near +Z edge
  let b = @graph.cube_from_direction(@graph.Vec3::new(1.0 - eps, 0.2, 1.0)) // +Z, near +X edge
  check_true(a.face == 0 && b.face == 4, "faces")
  check_close(a.u, 0.0, 1.0e-6, "+X meets +Z at u=0")
  check_close(b.u, 1.0, 1.0e-6, "+Z meets +X at u=1")
  check_close(a.v, b.v, 1.0e-6, "same v along the shared edge")
  let c = @graph.cube_from_direction(@graph.Vec3::new(0.3, 1.0, 1.0 - eps)) // +Y near +Z
  let d = @graph.cube_from_direction(@graph.Vec3::new(0.3, 1.0 - eps, 1.0)) // +Z near +Y
  check_true(c.face == 2 && d.face == 4, "faces 2")
  check_close(c.v, 1.0, 1.0e-6, "+Y bottom edge touches +Z")
  check_close(d.v, 0.0, 1.0e-6, "+Z top edge touches +Y")
  check_close(c.u, d.u, 1.0e-6, "same u along shared edge")
}

///|
test "equirectangular mapping" {
  let uv = @graph.equirect_from_direction(@graph.Vec3::new(0.0, 0.0, -1.0))
  check_v2(
    uv,
    @graph.Vec2::new(0.5, 0.5),
    1.0e-15,
    "forward is the image centre",
  )
  check_v2(
    @graph.equirect_from_direction(@graph.Vec3::unit_x()),
    @graph.Vec2::new(0.75, 0.5),
    1.0e-15,
    "+X right quarter",
  )
  check_v2(
    @graph.equirect_from_direction(-@graph.Vec3::unit_x()),
    @graph.Vec2::new(0.25, 0.5),
    1.0e-15,
    "-X left quarter",
  )
  check_close(
    @graph.equirect_from_direction(@graph.Vec3::unit_y()).y,
    0.0,
    1.0e-15,
    "up = top row",
  )
  check_close(
    @graph.equirect_from_direction(-@graph.Vec3::unit_y()).y,
    1.0,
    1.0e-15,
    "down = bottom row",
  )
  let back = @graph.equirect_from_direction(@graph.Vec3::unit_z())
  check_close(back.x.abs() - 0.5, 0.5, 1.0e-15, "behind = u at the seam")
  // consistent with the library's azimuth: u = 0.5 - azimuth / TAU
  let rng = Rng::new(6UL)
  for i in 0..<500 {
    let d = rng.unit_vec3()
    let e = @graph.equirect_from_direction(d)
    let s = @graph.Spherical::from_cartesian(d)
    if e.x > 1.0e-9 && e.x < 1.0 - 1.0e-9 {
      check_close(
        e.x,
        0.5 - s.azimuth / @graph.TAU,
        1.0e-12,
        "u vs azimuth #\{i}",
      )
    }
    check_close(
      e.y,
      0.5 - s.elevation / @graph.PI,
      1.0e-12,
      "v vs elevation #\{i}",
    )
    check_v3(@graph.equirect_to_direction(e), d, 1.0e-12, "round trip #\{i}")
  }
}

///|
test "octahedral mapping: round trip, domain, landmarks" {
  let rng = Rng::new(12UL)
  for i in 0..<3000 {
    let n = rng.unit_vec3()
    let e = @graph.oct_encode(n)
    check_true(
      e.x >= -1.0 && e.x <= 1.0 && e.y >= -1.0 && e.y <= 1.0,
      "domain #\{i}",
    )
    let d = @graph.oct_decode(e)
    check_close(d.length(), 1.0, 1.0e-14, "unit #\{i}")
    check_v3(d, n, 1.0e-12, "round trip #\{i}")
    // uv helpers invert each other
    check_v2(@graph.uv_to_oct(@graph.oct_to_uv(e)), e, 1.0e-15, "uv helpers")
  }
  check_v2(
    @graph.oct_encode(@graph.Vec3::unit_z()),
    @graph.Vec2::zero(),
    0.0,
    "+Z centre",
  )
  check_v2(
    @graph.oct_encode(@graph.Vec3::unit_x()),
    @graph.Vec2::new(1.0, 0.0),
    0.0,
    "+X",
  )
  check_v2(
    @graph.oct_encode(-@graph.Vec3::unit_x()),
    @graph.Vec2::new(-1.0, 0.0),
    0.0,
    "-X",
  )
  check_v2(
    @graph.oct_encode(@graph.Vec3::unit_y()),
    @graph.Vec2::new(0.0, 1.0),
    0.0,
    "+Y",
  )
  let nz = @graph.oct_encode(-@graph.Vec3::unit_z())
  check_true(nz.x.abs() == 1.0 && nz.y.abs() == 1.0, "-Z folds onto a corner")
  check_v3(@graph.oct_decode(nz), -@graph.Vec3::unit_z(), 1.0e-15, "-Z decode")
  // folded hemisphere: antipodal directions fold to mirrored points
  let p = @graph.Vec3::new(0.3, 0.5, 0.8).normalize()
  let q = @graph.Vec3::new(0.3, 0.5, -0.8).normalize()
  let ep = @graph.oct_encode(p)
  let eq = @graph.oct_encode(q)
  check_true(
    ep.x.abs() + ep.y.abs() <= 1.0 + 1.0e-15 &&
    eq.x.abs() + eq.y.abs() >= 1.0 - 1.0e-15,
    "upper inside diamond, lower outside",
  )
}
