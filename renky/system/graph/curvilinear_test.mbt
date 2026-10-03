// Polar / cylindrical / spherical coordinates.

///|
test "Polar" {
  check_v2(
    @graph.Polar::new(2.0, @graph.HALF_PI).to_cartesian(),
    @graph.Vec2::new(0.0, 2.0),
    1.0e-15,
    "90deg",
  )
  let p = @graph.Polar::from_cartesian(@graph.Vec2::new(-1.0, 0.0))
  check_close(p.radius, 1.0, 0.0, "r")
  check_close(p.angle, @graph.PI, 0.0, "angle pi")
  let z = @graph.Polar::from_cartesian(@graph.Vec2::zero())
  check_true(z.radius == 0.0 && z.angle == 0.0, "origin")
  let rng = Rng::new(50UL)
  for i in 0..<500 {
    let v = @graph.Vec2::new(rng.range(-9.0, 9.0), rng.range(-9.0, 9.0))
    let p = @graph.Polar::from_cartesian(v)
    check_true(
      p.radius >= 0.0 && p.angle > -@graph.PI - 1.0e-15 && p.angle <= @graph.PI,
      "ranges #\{i}",
    )
    check_v2(p.to_cartesian(), v, 1.0e-12, "round trip #\{i}")
  }
}

///|
test "Spherical: azimuth is yaw about +Y from -Z; elevation is above the XZ plane" {
  let s0 = @graph.Spherical::new(2.0, 0.0, 0.0).to_cartesian()
  check_v3(s0, @graph.Vec3::new(0.0, 0.0, -2.0), 1.0e-15, "az 0 = -Z")
  check_v3(
    @graph.Spherical::new(1.0, @graph.HALF_PI, 0.0).to_cartesian(),
    @graph.Vec3::new(-1.0, 0.0, 0.0),
    1.0e-15,
    "az 90 = -X",
  )
  check_v3(
    @graph.Spherical::new(1.0, -@graph.HALF_PI, 0.0).to_cartesian(),
    @graph.Vec3::unit_x(),
    1.0e-15,
    "az -90 = +X",
  )
  check_v3(
    @graph.Spherical::new(1.0, @graph.PI, 0.0).to_cartesian(),
    @graph.Vec3::unit_z(),
    1.0e-15,
    "az 180 = +Z",
  )
  check_v3(
    @graph.Spherical::new(3.0, 0.7, @graph.HALF_PI).to_cartesian(),
    @graph.Vec3::new(0.0, 3.0, 0.0),
    1.0e-15,
    "el 90 = +Y",
  )
  let rng = Rng::new(51UL)
  for i in 0..<500 {
    let az = rng.range(-@graph.PI + 1.0e-6, @graph.PI)
    let el = rng.range(-1.5, 1.5)
    let r = rng.range(0.1, 50.0)
    let v = @graph.Spherical::new(r, az, el).to_cartesian()
    check_close(v.length(), r, 1.0e-13, "radius #\{i}")
    // equals the quaternion yaw/pitch recipe
    let q = @graph.Quat::from_y(az) * @graph.Quat::from_x(el)
    check_v3(
      q.rotate(@graph.Vec3::new(0.0, 0.0, -r)),
      v,
      1.0e-12,
      "quaternion recipe #\{i}",
    )
    let b = @graph.Spherical::from_cartesian(v)
    check_close(b.radius, r, 1.0e-13, "r back")
    check_close(b.azimuth, az, 1.0e-12, "az back #\{i}")
    check_close(b.elevation, el, 1.0e-12, "el back #\{i}")
    let (yw, pt) = v.to_yaw_pitch()
    check_close(yw, az, 1.0e-12, "yaw")
    check_close(pt, el, 1.0e-12, "pitch")
    check_v3(
      @graph.Vec3::from_yaw_pitch(yw, pt).scale(r),
      v,
      1.0e-12,
      "from_yaw_pitch",
    )
  }
  // poles and origin
  let np = @graph.Spherical::from_cartesian(@graph.Vec3::new(0.0, 5.0, 0.0))
  check_close(np.elevation, @graph.HALF_PI, 1.0e-15, "pole el")
  check_close(np.azimuth, 0.0, 0.0, "pole az")
  let sp = @graph.Spherical::from_cartesian(@graph.Vec3::new(0.0, -5.0, 0.0))
  check_close(sp.elevation, -@graph.HALF_PI, 1.0e-15, "south pole el")
  let o = @graph.Spherical::from_cartesian(@graph.Vec3::zero())
  check_true(
    o.radius == 0.0 && o.azimuth == 0.0 && o.elevation == 0.0,
    "origin",
  )
  // anti-meridian
  check_close(
    @graph.Spherical::from_cartesian(@graph.Vec3::new(0.0, 0.0, 1.0)).azimuth.abs(),
    @graph.PI,
    1.0e-15,
    "azimuth pi",
  )
}

///|
test "Cylindrical" {
  check_v3(
    @graph.Cylindrical::new(2.0, 0.0, 3.0).to_cartesian(),
    @graph.Vec3::new(0.0, 3.0, -2.0),
    1.0e-15,
    "az 0",
  )
  check_v3(
    @graph.Cylindrical::new(2.0, @graph.HALF_PI, -1.0).to_cartesian(),
    @graph.Vec3::new(-2.0, -1.0, 0.0),
    1.0e-15,
    "az 90",
  )
  let rng = Rng::new(52UL)
  for i in 0..<500 {
    let v = rng.vec3(-9.0, 9.0)
    let c = @graph.Cylindrical::from_cartesian(v)
    check_true(c.radius >= 0.0, "radius")
    check_v3(c.to_cartesian(), v, 1.0e-12, "round trip #\{i}")
    // consistent with Spherical on the horizontal plane
    let s = @graph.Spherical::from_cartesian(v)
    check_close(c.azimuth, s.azimuth, 1.0e-12, "az agree")
  }
  let axis = @graph.Cylindrical::from_cartesian(@graph.Vec3::new(0.0, 4.0, 0.0))
  check_true(
    axis.radius == 0.0 && axis.azimuth == 0.0 && axis.height == 4.0,
    "on axis",
  )
}
