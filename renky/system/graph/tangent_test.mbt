// Frames, tangent space, barycentrics, perspective-correct interpolation.

///|
fn mk_frame_checks(n : @graph.Vec3, tag : String) -> Unit raise {
  let (t, b) = @graph.orthonormal_basis(n)
  check_close(t.length(), 1.0, 1.0e-12, tag + " |t|")
  check_close(b.length(), 1.0, 1.0e-12, tag + " |b|")
  check_close(t.dot(b), 0.0, 1.0e-12, tag + " t.b")
  check_close(t.dot(n), 0.0, 1.0e-12, tag + " t.n")
  check_close(b.dot(n), 0.0, 1.0e-12, tag + " b.n")
  check_v3(t.cross(b), n, 1.0e-12, tag + " right-handed")
}

///|
test "orthonormal_basis: axes, poles, near-singular and random directions" {
  for
    a in [
      @graph.Vec3::unit_x(),
      -@graph.Vec3::unit_x(),
      @graph.Vec3::unit_y(),
      -@graph.Vec3::unit_y(),
      @graph.Vec3::unit_z(),
      -@graph.Vec3::unit_z(),
      @graph.Vec3::new(1.0e-9, 0.0, -1.0).normalize(),
      @graph.Vec3::new(0.0, 1.0e-9, 1.0).normalize(),
    ] {
    mk_frame_checks(a, "axis")
  }
  let rng = Rng::new(404UL)
  for i in 0..<2000 {
    mk_frame_checks(rng.unit_vec3(), "rand #\{i}")
  }
}

///|
test "Frame: to_world/to_local inverse, handedness, from_normal" {
  let rng = Rng::new(9UL)
  for i in 0..<200 {
    let f = @graph.Frame::from_normal(rng.vec3(-1.0, 1.0))
    check_close(f.handedness(), 1.0, 0.0, "right-handed")
    let v = rng.vec3(-5.0, 5.0)
    check_v3(f.to_world(f.to_local(v)), v, 1.0e-12, "roundtrip #\{i}")
    check_v3(
      f.to_world(@graph.Vec3::new(0.0, 0.0, 1.0)),
      f.normal,
      1.0e-15,
      "z -> n",
    )
    check_v3(f.to_mat3().mul_vec3(v), f.to_world(v), 1.0e-12, "mat3 #\{i}")
    check_close(f.to_mat3().determinant(), 1.0, 1.0e-12, "proper")
  }
  let mirrored : @graph.Frame = {
    tangent: @graph.Vec3::unit_x(),
    bitangent: @graph.Vec3::unit_z(),
    normal: @graph.Vec3::unit_y(),
  }
  check_close(mirrored.handedness(), -1.0, 0.0, "mirrored")
}

///|
test "tangent_frame follows dP/du and dP/dv of an affine UV mapping" {
  let rng = Rng::new(61UL)
  for i in 0..<500 {
    // random surface patch P(u,v) = O + U u + V v with non-orthogonal U, V
    let o = rng.vec3(-5.0, 5.0)
    let uu = rng.vec3(-3.0, 3.0)
    let vv = rng.vec3(-3.0, 3.0)
    let nrm = uu.cross(vv)
    if nrm.length() < 0.3 {
      continue
    }
    // random UV placement incl. mirrored (negative determinant) layouts
    let su = if i % 2 == 0 { 1.0 } else { -1.0 }
    let uv_of = fn(u : Double, v : Double) -> @graph.Vec2 {
      @graph.Vec2::new(su * u * 0.5 + 0.1, v * 0.25 + 0.3)
    }
    // P as a function of (u', v') in UV space: u = (u' - 0.1) / (su*0.5), v = (v' - 0.3) / 0.25
    let p_of = fn(u : Double, v : Double) -> @graph.Vec3 {
      o + uu.scale(u) + vv.scale(v)
    }
    let (a_u, a_v) = (0.2, 0.1)
    let (b_u, b_v) = (1.7, 0.3)
    let (c_u, c_v) = (0.5, 1.9)
    let n = nrm.normalize().scale(su)
    match
      @graph.tangent_frame(
        p_of(a_u, a_v),
        p_of(b_u, b_v),
        p_of(c_u, c_v),
        uv_of(a_u, a_v),
        uv_of(b_u, b_v),
        uv_of(c_u, c_v),
        n,
      ) {
      Some(f) => {
        check_close(f.tangent.length(), 1.0, 1.0e-12, "|t|")
        check_close(f.bitangent.length(), 1.0, 1.0e-12, "|b|")
        check_close(f.tangent.dot(f.normal), 0.0, 1.0e-12, "t.n")
        check_close(f.bitangent.dot(f.normal), 0.0, 1.0e-12, "b.n")
        check_close(f.tangent.dot(f.bitangent), 0.0, 1.0e-12, "t.b")
        // dP/du' and dP/dv' (derivatives w.r.t. UV coordinates)
        let dpdu = uu.scale(1.0 / (su * 0.5))
        let dpdv = vv.scale(1.0 / 0.25)
        check_true(f.tangent.dot(dpdu) > 0.0, "tangent along +u #\{i}")
        check_true(f.bitangent.dot(dpdv) > 0.0, "bitangent along +v #\{i}")
        // T = normalised projection of dP/du onto the tangent plane
        let proj = dpdu - f.normal.scale(f.normal.dot(dpdu))
        check_v3(f.tangent, proj.normalize(), 1.0e-9, "T value #\{i}")
        // handedness reflects the UV orientation relative to the supplied normal
        let expected_h = @graph.sign_nz(f.normal.cross(f.tangent).dot(dpdv))
        check_close(f.handedness(), expected_h, 0.0, "handedness #\{i}")
      }
      None => fail("unexpected None #\{i}")
    }
  }
  // degenerate UVs
  let z = @graph.Vec2::new(0.3, 0.3)
  check_true(
    @graph.tangent_frame(
      @graph.Vec3::zero(),
      @graph.Vec3::unit_x(),
      @graph.Vec3::unit_y(),
      z,
      z,
      z,
      @graph.Vec3::unit_z(),
    )
    is None,
    "degenerate -> None",
  )
}

///|
test "tangent_frame: normal-map vector (1,0,0) moves toward increasing u on the surface" {
  // quad in the XZ plane facing +Y, u along +X, v along +Z
  let p0 = @graph.Vec3::new(0.0, 0.0, 0.0)
  let p1 = @graph.Vec3::new(2.0, 0.0, 0.0)
  let p2 = @graph.Vec3::new(0.0, 0.0, 4.0)
  match
    @graph.tangent_frame(
      p0,
      p1,
      p2,
      @graph.Vec2::new(0.0, 0.0),
      @graph.Vec2::new(1.0, 0.0),
      @graph.Vec2::new(0.0, 1.0),
      @graph.Vec3::unit_y(),
    ) {
    Some(f) => {
      check_v3(
        f.to_world(@graph.Vec3::new(1.0, 0.0, 0.0)),
        @graph.Vec3::unit_x(),
        1.0e-15,
        "T=+x",
      )
      // v increases along +z (top-left UV convention: that is DOWN in the image), so B = +z
      check_v3(
        f.to_world(@graph.Vec3::new(0.0, 1.0, 0.0)),
        @graph.Vec3::unit_z(),
        1.0e-15,
        "B=+z",
      )
      check_v3(
        f.to_world(@graph.Vec3::new(0.0, 0.0, 1.0)),
        @graph.Vec3::unit_y(),
        1.0e-15,
        "N=+y",
      )
      // handedness: T x B = x cross z = -y  => opposite to the +Y normal => left-handed
      check_close(f.handedness(), -1.0, 0.0, "x cross z is -y, left-handed")
    }
    None => fail("None")
  }
}

///|
test "barycentric 2D: vertices, centroid, outside, reconstruction, degenerate" {
  let a = @graph.Vec2::new(1.0, 1.0)
  let b = @graph.Vec2::new(5.0, 2.0)
  let c = @graph.Vec2::new(2.0, 6.0)
  check_v3(
    @graph.barycentric_2d(a, a, b, c).unwrap(),
    @graph.Vec3::new(1.0, 0.0, 0.0),
    1.0e-15,
    "at a",
  )
  check_v3(
    @graph.barycentric_2d(b, a, b, c).unwrap(),
    @graph.Vec3::new(0.0, 1.0, 0.0),
    1.0e-15,
    "at b",
  )
  check_v3(
    @graph.barycentric_2d(c, a, b, c).unwrap(),
    @graph.Vec3::new(0.0, 0.0, 1.0),
    1.0e-15,
    "at c",
  )
  let g = (a + b + c).scale(1.0 / 3.0)
  check_v3(
    @graph.barycentric_2d(g, a, b, c).unwrap(),
    @graph.Vec3::splat(1.0 / 3.0),
    1.0e-15,
    "centroid",
  )
  let out = @graph.barycentric_2d(@graph.Vec2::new(0.0, 0.0), a, b, c).unwrap()
  check_true(!@graph.bary_inside(out, 0.0), "outside")
  check_close(out.x + out.y + out.z, 1.0, 1.0e-14, "sum 1")
  let rng = Rng::new(15UL)
  for i in 0..<500 {
    let p = @graph.Vec2::new(rng.range(-3.0, 9.0), rng.range(-3.0, 9.0))
    let w = @graph.barycentric_2d(p, a, b, c).unwrap()
    check_v2(
      a.scale(w.x) + b.scale(w.y) + c.scale(w.z),
      p,
      1.0e-12,
      "reconstruct #\{i}",
    )
    check_close(w.x + w.y + w.z, 1.0, 1.0e-13, "sum #\{i}")
  }
  // orientation: CCW (y-up) has positive area, swapping two vertices flips the sign
  check_true(
    @graph.triangle_area2(a, b, c) > 0.0 && @graph.triangle_area2(a, c, b) < 0.0,
    "area sign",
  )
  check_close(
    @graph.triangle_area2(
      @graph.Vec2::zero(),
      @graph.Vec2::new(2.0, 0.0),
      @graph.Vec2::new(0.0, 3.0),
    ),
    6.0,
    0.0,
    "area",
  )
  check_close(
    @graph.edge_function(
      @graph.Vec2::zero(),
      @graph.Vec2::new(1.0, 0.0),
      @graph.Vec2::new(0.5, 1.0),
    ),
    1.0,
    0.0,
    "left of edge is positive",
  )
  check_true(@graph.barycentric_2d(g, a, a, c) is None, "degenerate")
}

///|
test "barycentric 3D: reconstruction and plane projection" {
  let rng = Rng::new(16UL)
  for i in 0..<300 {
    let a = rng.vec3(-3.0, 3.0)
    let b = rng.vec3(-3.0, 3.0)
    let c = rng.vec3(-3.0, 3.0)
    if (b - a).cross(c - a).length() < 0.1 {
      continue
    }
    let lam = rng.vec3(0.0, 1.0)
    let lam = lam.scale(1.0 / (lam.x + lam.y + lam.z))
    let p = a.scale(lam.x) + b.scale(lam.y) + c.scale(lam.z)
    check_v3(
      @graph.barycentric_3d(p, a, b, c).unwrap(),
      lam,
      1.0e-10,
      "on-plane #\{i}",
    )
    // lifting the point off the plane does not change the weights
    let nrm = (b - a).cross(c - a).normalize()
    check_v3(
      @graph.barycentric_3d(p + nrm.scale(2.5), a, b, c).unwrap(),
      lam,
      1.0e-10,
      "off-plane #\{i}",
    )
  }
  check_true(
    @graph.barycentric_3d(
      @graph.Vec3::zero(),
      @graph.Vec3::zero(),
      @graph.Vec3::unit_x(),
      @graph.Vec3::unit_x().scale(2.0),
    )
    is None,
    "collinear",
  )
}

///|
test "perspective-correct interpolation recovers the true 3D weights" {
  let rng = Rng::new(17UL)
  let proj = @graph.Projection::perspective(1.2, 1.3, 0.1, 100.0, ZeroToOne)
  let mut worst_naive = 0.0
  for i in 0..<500 {
    // triangle in view space with strongly different depths
    let v0 = @graph.Vec3::new(
      rng.range(-2.0, 2.0),
      rng.range(-2.0, 2.0),
      -rng.range(1.0, 3.0),
    )
    let v1 = @graph.Vec3::new(
      rng.range(-2.0, 2.0),
      rng.range(-2.0, 2.0),
      -rng.range(6.0, 20.0),
    )
    let v2 = @graph.Vec3::new(
      rng.range(-2.0, 2.0),
      rng.range(-2.0, 2.0),
      -rng.range(2.0, 12.0),
    )
    let lam = rng.vec3(0.05, 1.0)
    let lam = lam.scale(1.0 / (lam.x + lam.y + lam.z))
    let p = v0.scale(lam.x) + v1.scale(lam.y) + v2.scale(lam.z)
    let n0 = proj.project_view(v0).unwrap()
    let n1 = proj.project_view(v1).unwrap()
    let n2 = proj.project_view(v2).unwrap()
    let np = proj.project_view(p).unwrap()
    match @graph.barycentric_2d(np.xy(), n0.xy(), n1.xy(), n2.xy()) {
      Some(sb) => {
        let pc = @graph.perspective_correct(sb, -v0.z, -v1.z, -v2.z)
        check_v3(pc, lam, 1.0e-9, "perspective-correct #\{i}")
        check_v3(
          @graph.interpolate3(pc, v0, v1, v2),
          p,
          1.0e-8,
          "interpolated position #\{i}",
        )
        // screen-space z (NDC depth) IS affine in screen space (the reason z-buffering works)
        check_close(
          @graph.interpolate1(sb, n0.z, n1.z, n2.z),
          np.z,
          1.0e-9,
          "ndc z affine #\{i}",
        )
        worst_naive = worst_naive.max((sb - lam).abs().x)
      }
      None => ()
    }
  }
  check_true(
    worst_naive > 0.05,
    "the test triangles really need the correction (naive error \{worst_naive})",
  )
  // attribute interpolation helpers
  let b = @graph.Vec3::new(0.2, 0.3, 0.5)
  check_close(
    @graph.interpolate1(b, 1.0, 2.0, 3.0),
    0.2 + 0.6 + 1.5,
    1.0e-15,
    "interp1",
  )
  check_v2(
    @graph.interpolate2(
      b,
      @graph.Vec2::new(1.0, 0.0),
      @graph.Vec2::new(0.0, 1.0),
      @graph.Vec2::new(1.0, 1.0),
    ),
    @graph.Vec2::new(0.7, 0.8),
    1.0e-15,
    "interp2",
  )
}
