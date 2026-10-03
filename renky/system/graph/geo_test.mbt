// Geodesy + web-map tests: references are PROJ (pyproj) outputs and the OSM tile formulas.

///|
test "WGS84 constants" {
  let e = @graph.Ellipsoid::wgs84()
  check_close(e.b(), 6356752.314245179, 1.0e-15, "b")
  check_close(e.e2(), 0.00669437999014132, 1.0e-12, "e2")
  check_close(e.prime_vertical_radius(0.0), 6378137.0, 1.0e-15, "N(equator)")
  // polar radius of curvature c = a^2 / b
  check_close(
    e.prime_vertical_radius(@graph.HALF_PI),
    6399593.625758493,
    1.0e-12,
    "N(pole)",
  )
  check_close(e.meridional_radius(0.0), 6335439.327292, 1.0e-10, "M(equator)")
  check_close(
    e.meridional_radius(@graph.HALF_PI),
    6399593.625758493,
    1.0e-12,
    "M(pole)",
  )
  // 45 deg geodetic -> 44.8076 deg geocentric
  check_close(
    @graph.rad_to_deg(e.geocentric_latitude(@graph.deg_to_rad(45.0))),
    44.80757678,
    1.0e-9,
    "geocentric lat",
  )
  check_close(
    e.geodetic_latitude(e.geocentric_latitude(0.7)),
    0.7,
    1.0e-14,
    "lat roundtrip",
  )
  // metres per degree of latitude at the equator ~ 110574 m
  check_close(
    e.meridional_radius(0.0) * @graph.DEG2RAD,
    110574.2758217,
    1.0e-9,
    "m/deg",
  )
}

///|
test "Geodetic -> ECEF matches PROJ; ECEF -> Geodetic inverts it (poles, equator, space, underground)" {
  let g = gold_geo
  let e = @graph.Ellipsoid::wgs84()
  for i in 0..<(g.length() / 6) {
    let o = i * 6
    let geo = @graph.Geodetic::from_degrees(g[o], g[o + 1], g[o + 2])
    let want = g_v3(g, o + 3)
    let got = e.to_ecef(geo)
    // PROJ and this code agree to well below a micron
    check_close(got.x, want.x, 2.0e-14, "ecef x #\{i}")
    check_close(got.y, want.y, 2.0e-14, "ecef y #\{i}")
    check_close(got.z, want.z, 2.0e-14, "ecef z #\{i}")
    let back = e.from_ecef(want)
    check_close(back.lat, geo.lat, 1.0e-12, "lat #\{i}")
    if g[o].abs() < 89.99 {
      check_close(@graph.wrap_pi(back.lon - geo.lon), 0.0, 1.0e-12, "lon #\{i}")
    }
    check_close(back.height, geo.height, 1.0e-9, "height #\{i}")
  }
}

///|
test "ECEF special points: equator, pole, centre, axis conventions" {
  let e = @graph.Ellipsoid::wgs84()
  check_v3(
    e.to_ecef(@graph.Geodetic::new(0.0, 0.0, 0.0)),
    @graph.Vec3::new(6378137.0, 0.0, 0.0),
    1.0e-9,
    "lon 0",
  )
  check_v3(
    e.to_ecef(@graph.Geodetic::from_degrees(0.0, 90.0, 0.0)),
    @graph.Vec3::new(0.0, 6378137.0, 0.0),
    1.0e-6,
    "+Y is 90E",
  )
  check_v3(
    e.to_ecef(@graph.Geodetic::from_degrees(90.0, 0.0, 0.0)),
    @graph.Vec3::new(0.0, 0.0, e.b()),
    1.0e-6,
    "north pole +Z",
  )
  let np = e.from_ecef(@graph.Vec3::new(0.0, 0.0, e.b() + 1000.0))
  check_close(np.lat, @graph.HALF_PI, 1.0e-14, "pole lat")
  check_close(np.height, 1000.0, 1.0e-8, "pole height")
  let sp = e.from_ecef(@graph.Vec3::new(0.0, 0.0, -e.b()))
  check_close(sp.lat, -@graph.HALF_PI, 1.0e-14, "south pole lat")
  check_close(sp.height, 0.0, 1.0e-8, "south pole height")
  let c = e.from_ecef(@graph.Vec3::zero())
  check_close(c.height, -6378137.0, 1.0e-9, "centre")
  // sphere ellipsoid: ECEF is plain spherical coordinates
  let s = @graph.Ellipsoid::sphere(1000.0)
  check_v3(
    s.to_ecef(@graph.Geodetic::from_degrees(30.0, 60.0, 0.0)),
    @graph.Spherical::new(1000.0, 0.0, 0.0).to_cartesian().scale(0.0) +
    @graph.Vec3::new(
      1000.0 * @math.cos(@graph.deg_to_rad(30.0)) * 0.5,
      1000.0 *
      @math.cos(@graph.deg_to_rad(30.0)) *
      @math.sin(@graph.deg_to_rad(60.0)),
      500.0,
    ),
    1.0e-9,
    "sphere",
  )
  // heights stay consistent along the normal
  let g = @graph.Geodetic::from_degrees(40.0, -75.0, 0.0)
  let p0 = e.to_ecef(g)
  let p1 = e.to_ecef({ ..g, height: 123.0, })
  check_v3(p1 - p0, g.up_ecef().scale(123.0), 1.0e-8, "height along normal")
  // GRS80 differs from WGS84 by < 0.1 mm in height at the pole
  let gr = @graph.Ellipsoid::grs80()
  check_true((gr.b() - e.b()).abs() < 1.0e-3, "grs80 ~ wgs84")
}

///|
test "ENU/NED frames == PROJ topocentric; frames are orthonormal and right-handed" {
  let g = gold_enu
  let e = @graph.Ellipsoid::wgs84()
  for i in 0..<(g.length() / 9) {
    let o = i * 9
    let origin = @graph.Geodetic::from_degrees(g[o], g[o + 1], g[o + 2])
    let pt = @graph.Geodetic::from_degrees(g[o + 3], g[o + 4], g[o + 5])
    let want = g_v3(g, o + 6)
    let frame = @graph.LocalFrame::enu(e, origin)
    let ecef = e.to_ecef(pt)
    check_v3(frame.to_local(ecef), want, 1.0e-9, "ENU #\{i}")
    check_v3(frame.to_ecef(want), ecef, 1.0e-6, "ENU -> ECEF #\{i}")
    let ned = @graph.LocalFrame::ned(e, origin)
    check_v3(
      ned.to_local(ecef),
      @graph.Vec3::new(want.y, want.x, -want.z),
      1.0e-9,
      "NED #\{i}",
    )
    // axes orthonormal, E x N = U
    let a = frame.axes
    check_m3(a * a.transpose(), @graph.Mat3::identity(), 1.0e-14, "orthonormal")
    check_close(a.determinant(), 1.0, 1.0e-14, "proper rotation")
    check_v3(a.col(0).cross(a.col(1)), a.col(2), 1.0e-14, "E x N = U")
    check_close(ned.axes.determinant(), 1.0, 1.0e-14, "NED proper")
    // up is the ellipsoid normal: moving 1 m along up raises height by 1 m, lat/lon unchanged
    let up_pt = e.from_ecef(frame.to_ecef(@graph.Vec3::new(0.0, 0.0, 1.0)))
    check_close(up_pt.height, origin.height + 1.0, 1.0e-8, "up = +height")
    check_close(up_pt.lat, origin.lat, 1.0e-12, "up keeps lat")
    // moving 100 m north increases latitude by 100 / M
    let n_pt = e.from_ecef(frame.to_ecef(@graph.Vec3::new(0.0, 100.0, 0.0)))
    let dlat = n_pt.lat - origin.lat
    check_close(
      dlat,
      100.0 / (e.meridional_radius(origin.lat) + origin.height),
      1.0e-5,
      "north step",
    )
    check_true(
      frame
      .dir_to_local(frame.axes.col(1))
      .approx_eq(@graph.Vec3::unit_y(), 1.0e-14),
      "dir to local",
    )
  }
}

///|
test "render mapping: east -> +X, up -> +Y, north -> -Z; globe Y-up axes" {
  let e = @graph.Ellipsoid::wgs84()
  let g = @graph.Geodetic::from_degrees(22.3, 114.2, 5.0)
  let f = @graph.LocalFrame::enu(e, g)
  let to_r = f.ecef_to_render_enu()
  check_v3(to_r.mul_point(f.origin), @graph.Vec3::zero(), 1.0e-7, "origin")
  check_v3(
    to_r.mul_point(f.to_ecef(@graph.Vec3::new(1.0, 0.0, 0.0))),
    @graph.Vec3::unit_x(),
    1.0e-9,
    "east",
  )
  check_v3(
    to_r.mul_point(f.to_ecef(@graph.Vec3::new(0.0, 0.0, 1.0))),
    @graph.Vec3::unit_y(),
    1.0e-9,
    "up",
  )
  check_v3(
    to_r.mul_point(f.to_ecef(@graph.Vec3::new(0.0, 1.0, 0.0))),
    @graph.Vec3::new(0.0, 0.0, -1.0),
    1.0e-9,
    "north",
  )
  check_close(
    to_r.to_mat4().to_mat3().determinant(),
    1.0,
    1.0e-12,
    "proper rotation (no mirror)",
  )
  // a camera standing at the origin of the ENU frame looking north has forward = -Z in render space
  // and sees a point 100 m north dead-centre
  let cam = @graph.Camera::look_at(
    @graph.Vec3::new(0.0, 2.0, 0.0),
    @graph.Vec3::new(0.0, 2.0, -100.0),
    @graph.Vec3::unit_y(),
    @graph.Projection::perspective(1.0, 1.0, 0.1, 1000.0, ZeroToOne),
    @graph.Viewport::full(200.0, 200.0),
  )
  let world = to_r.mul_point(f.to_ecef(@graph.Vec3::new(0.0, 100.0, 2.0)))
  let s = cam.world_to_screen(world).unwrap()
  check_close(s.x, 100.0, 1.0e-6, "centre x")
  check_close(s.y, 100.0, 1.0e-6, "centre y")
  // globe axes: north pole is +Y, (lat 0, lon 90E) is -Z
  let np = @graph.ecef_to_yup(
    e.to_ecef(@graph.Geodetic::from_degrees(90.0, 0.0, 0.0)),
  )
  check_v3(np, @graph.Vec3::new(0.0, e.b(), 0.0), 1.0e-6, "pole +Y")
  let e90 = @graph.ecef_to_yup(@graph.Vec3::new(0.0, 1.0, 0.0))
  check_v3(e90, @graph.Vec3::new(0.0, 0.0, -1.0), 0.0, "90E -> -Z")
  check_v3(
    @graph.yup_to_ecef(@graph.ecef_to_yup(@graph.Vec3::new(1.0, 2.0, 3.0))),
    @graph.Vec3::new(1.0, 2.0, 3.0),
    0.0,
    "roundtrip",
  )
}

///|
test "Web Mercator matches PROJ EPSG:4326 -> EPSG:3857 and inverts" {
  let g = gold_merc
  for i in 0..<(g.length() / 4) {
    let o = i * 4
    let ll = @graph.LonLat::from_degrees(g[o], g[o + 1])
    let m = @graph.web_mercator_forward(ll)
    check_close(m.x, g[o + 2], 1.0e-13, "x #\{i}")
    check_close(m.y, g[o + 3], 1.0e-12, "y #\{i}")
    let back = @graph.web_mercator_inverse(@graph.Vec2::new(g[o + 2], g[o + 3]))
    check_close(back.lon, ll.lon, 1.0e-12, "lon #\{i}")
    check_close(back.lat, ll.lat, 1.0e-11, "lat #\{i}")
  }
  // clipping latitude: beyond 85.0511 deg the value saturates at the world edge
  let edge = @graph.web_mercator_forward(@graph.LonLat::from_degrees(0.0, 89.0)).y
  check_close(edge, 20037508.342789244, 1.0e-13, "clamped y")
  check_close(
    @graph.rad_to_deg(@graph.MERCATOR_MAX_LAT),
    85.0511287798066,
    1.0e-14,
    "max lat",
  )
}

///|
test "tile pyramid matches OSM formulas; quadkeys; ground resolution" {
  let g = gold_tile
  for i in 0..<(g.length() / 9) {
    let o = i * 9
    let ll = @graph.LonLat::from_degrees(g[o], g[o + 1])
    let z = g[o + 2].to_int()
    let (tx, ty) = @graph.lonlat_to_tile(ll, z, 256.0)
    check_true(
      tx == g[o + 3].to_int() && ty == g[o + 4].to_int(),
      "tile #\{i} z\{z}: got (\{tx},\{ty}) want (\{g[o + 3]},\{g[o + 4]})",
    )
    let nw = @graph.tile_to_lonlat(
      Double::from_int(tx),
      Double::from_int(ty),
      z,
    )
    check_close(@graph.rad_to_deg(nw.lon), g[o + 5], 1.0e-12, "nw lon #\{i}")
    check_close(@graph.rad_to_deg(nw.lat), g[o + 6], 1.0e-11, "nw lat #\{i}")
    let wp = @graph.lonlat_to_world_pixel(ll, Double::from_int(z), 256.0)
    check_close(wp.x, g[o + 7], 1.0e-12, "world px x #\{i}")
    check_close(wp.y, g[o + 8], 1.0e-11, "world px y #\{i}")
    let back = @graph.world_pixel_to_lonlat(wp, Double::from_int(z), 256.0)
    check_close(back.lon, ll.lon, 1.0e-11, "px->lon #\{i}")
    check_close(back.lat, ll.lat, 1.0e-10, "px->lat #\{i}")
    // meters <-> pixels agree with lon/lat <-> pixels
    let m = @graph.web_mercator_forward(ll)
    check_v2(
      @graph.mercator_to_world_pixel(m, Double::from_int(z), 256.0),
      wp,
      1.0e-12,
      "merc->px #\{i}",
    )
    check_v2(
      @graph.world_pixel_to_mercator(wp, Double::from_int(z), 256.0),
      m,
      1.0e-12,
      "px->merc #\{i}",
    )
    // the point lies inside its tile's bounds
    let (a, b) = @graph.tile_bounds(tx, ty, z)
    check_true(
      ll.lon >= a.lon - 1.0e-12 &&
      ll.lon <= b.lon + 1.0e-12 &&
      ll.lat <= a.lat + 1.0e-12 &&
      ll.lat >= b.lat - 1.0e-12,
      "inside bounds #\{i}",
    )
    // quadkey round trip
    let key = @graph.tile_to_quadkey(tx, ty, z)
    check_true(key.length() == z, "quadkey length")
    match @graph.quadkey_to_tile(key) {
      Some((qx, qy, qz)) =>
        check_true(qx == tx && qy == ty && qz == z, "quadkey roundtrip #\{i}")
      None => fail("quadkey parse")
    }
  }
  // documented example of the Bing Maps tile system: tile (3, 5) at level 3 -> "213"
  check_true(@graph.tile_to_quadkey(3, 5, 3) == "213", "bing example")
  check_true(@graph.quadkey_to_tile("213") is Some((3, 5, 3)), "bing parse")
  check_true(@graph.quadkey_to_tile("2x3") is None, "bad key")
  check_true(@graph.tile_to_quadkey(0, 0, 0) == "", "z0 key")
  // z0 / equator: 156543.03392804097 m per 256 px pixel
  check_close(
    @graph.ground_resolution(0.0, 0.0, 256.0),
    156543.03392804097,
    1.0e-14,
    "z0 resolution",
  )
  check_close(
    @graph.ground_resolution(0.0, 1.0, 256.0),
    156543.03392804097 / 2.0,
    1.0e-14,
    "z1 resolution",
  )
  check_close(
    @graph.ground_resolution(@graph.deg_to_rad(60.0), 3.0, 256.0),
    156543.03392804097 / 8.0 * 0.5,
    1.0e-13,
    "60 deg",
  )
  // origin of the pyramid
  let (ox, oy) = @graph.lonlat_to_tile(
    @graph.LonLat::from_degrees(-179.9, 85.0),
    5,
    256.0,
  )
  check_true(ox == 0 && oy == 0, "NW tile")
  let (sx, sy) = @graph.lonlat_to_tile(
    @graph.LonLat::from_degrees(179.9, -85.0),
    5,
    256.0,
  )
  check_true(sx == 31 && sy == 31, "SE tile")
  // known: zoom 1, lon 0 lat 0 is the centre corner => x = 1 (east side by floor), y = 1
  let (cx, cy) = @graph.lonlat_to_tile(
    @graph.LonLat::from_degrees(0.0, 0.0),
    1,
    256.0,
  )
  check_true(cx == 1 && cy == 1, "centre corner")
}
