// Geodetic coordinates on an oblate ellipsoid and the Earth-centred frames built on them.
//
//   Geodetic (lat, lon, height)  <->  ECEF (x, y, z)  <->  local tangent frame (ENU / NED)
//                                                    <->  render space (Y-up, -Z = north)
//
// ECEF: origin at the ellipsoid centre, +X through (lat 0, lon 0), +Y through (lat 0, lon +90 deg
// east), +Z through the north pole. Right-handed. Units: metres. Angles: radians.
// Latitude is GEODETIC (angle between the ellipsoid normal and the equatorial plane).

///|
#valtype
pub(all) struct Ellipsoid {
  /// semi-major (equatorial) radius, metres
  a : Double
  /// flattening  f = (a - b) / a
  f : Double
} derive(Eq, Debug)

///|
/// WGS 84 (GPS, Web maps, glTF/3D-Tiles ECEF).
pub fn Ellipsoid::wgs84() -> Ellipsoid {
  { a: 6378137.0, f: 1.0 / 298.257223563, }
}

///|
/// GRS 80 (NAD83, ETRS89, CGCS2000 differs slightly in f).
pub fn Ellipsoid::grs80() -> Ellipsoid {
  { a: 6378137.0, f: 1.0 / 298.257222101, }
}

///|
/// A perfect sphere (f = 0) -- the model used by simple globes and Web Mercator.
pub fn Ellipsoid::sphere(radius : Double) -> Ellipsoid {
  { a: radius, f: 0.0, }
}

///|
/// semi-minor (polar) radius
pub fn Ellipsoid::b(self : Ellipsoid) -> Double {
  self.a * (1.0 - self.f)
}

///|
/// first eccentricity squared
pub fn Ellipsoid::e2(self : Ellipsoid) -> Double {
  self.f * (2.0 - self.f)
}

///|
/// Radius of curvature in the prime vertical, N(lat): distance from the surface to the polar
/// axis along the normal.
pub fn Ellipsoid::prime_vertical_radius(
  self : Ellipsoid,
  lat : Double,
) -> Double {
  let s = @math.sin(lat)
  self.a / (1.0 - self.e2() * s * s).sqrt()
}

///|
/// Meridional radius of curvature M(lat) (north-south): metres per radian of latitude.
pub fn Ellipsoid::meridional_radius(self : Ellipsoid, lat : Double) -> Double {
  let e2 = self.e2()
  let s = @math.sin(lat)
  let w = 1.0 - e2 * s * s
  self.a * (1.0 - e2) / (w * w.sqrt())
}

///|
/// Geodetic position: latitude / longitude in RADIANS, height above the ellipsoid in metres.
#valtype
pub(all) struct Geodetic {
  lat : Double
  lon : Double
  height : Double
} derive(Eq, Debug)

///|
pub fn Geodetic::new(lat : Double, lon : Double, height : Double) -> Geodetic {
  { lat, lon, height, }
}

///|
pub fn Geodetic::from_degrees(
  lat_deg : Double,
  lon_deg : Double,
  height : Double,
) -> Geodetic {
  { lat: deg_to_rad(lat_deg), lon: deg_to_rad(lon_deg), height, }
}

///|
/// Geodetic -> ECEF (closed form).
pub fn Ellipsoid::to_ecef(self : Ellipsoid, g : Geodetic) -> Vec3 {
  let sl = @math.sin(g.lat)
  let cl = @math.cos(g.lat)
  let n = self.a / (1.0 - self.e2() * sl * sl).sqrt()
  {
    x: (n + g.height) * cl * @math.cos(g.lon),
    y: (n + g.height) * cl * @math.sin(g.lon),
    z: (n * (1.0 - self.e2()) + g.height) * sl,
  }
}

///|
/// ECEF -> Geodetic (Bowring start + fixed-point iteration on latitude; ~1e-6 m or better for any point from deep
/// underground to beyond geostationary orbit; poles and the equator are handled). The Earth's
/// exact centre has no unique latitude: (0, 0, -a) is returned.
pub fn Ellipsoid::from_ecef(self : Ellipsoid, p : Vec3) -> Geodetic {
  let e2 = self.e2()
  let rho = @math.hypot(p.x, p.y)
  if rho == 0.0 && p.z == 0.0 {
    return { lat: 0.0, lon: 0.0, height: -self.a, }
  }
  let lon = if rho == 0.0 { 0.0 } else { @math.atan2(p.y, p.x) }
  // Bowring's closed-form first guess (~1e-10 rad), then fixed-point polish: 2-3 sweeps, not ~10.
  let b = self.a * (1.0 - e2).sqrt()
  let beta = @math.atan2(p.z * self.a, rho * b)
  let sb = @math.sin(beta)
  let cb = @math.cos(beta)
  let ep2 = e2 / (1.0 - e2)
  let mut lat = @math.atan2(
    p.z + ep2 * b * sb * sb * sb,
    rho - e2 * self.a * cb * cb * cb,
  )
  for _ in 0..<50 {
    let s = @math.sin(lat)
    let n = self.a / (1.0 - e2 * s * s).sqrt()
    let nl = @math.atan2(p.z + e2 * n * s, rho)
    let done = (nl - lat).abs() < 1.0e-13
    lat = nl
    if done {
      break
    }
  }
  let s = @math.sin(lat)
  let c = @math.cos(lat)
  // numerically stable at the poles as well as the equator
  let h = rho * c + p.z * s - self.a * (1.0 - e2 * s * s).sqrt()
  { lat, lon, height: h, }
}

///|
/// Geocentric latitude (angle from the equatorial plane to the position vector) for a
/// point ON the ellipsoid surface at geodetic latitude `lat`.
pub fn Ellipsoid::geocentric_latitude(self : Ellipsoid, lat : Double) -> Double {
  @math.atan((1.0 - self.e2()) * @math.tan(lat))
}

///|
/// Inverse of `geocentric_latitude`.
pub fn Ellipsoid::geodetic_latitude(
  self : Ellipsoid,
  geocentric : Double,
) -> Double {
  @math.atan(@math.tan(geocentric) / (1.0 - self.e2()))
}

///|
/// Unit ellipsoid normal ("up") in ECEF at a geodetic position (independent of height).
pub fn Geodetic::up_ecef(self : Geodetic) -> Vec3 {
  let c = @math.cos(self.lat)
  {
    x: c * @math.cos(self.lon),
    y: c * @math.sin(self.lon),
    z: @math.sin(self.lat),
  }
}

///|
pub fn Geodetic::east_ecef(self : Geodetic) -> Vec3 {
  { x: -@math.sin(self.lon), y: @math.cos(self.lon), z: 0.0, }
}

///|
pub fn Geodetic::north_ecef(self : Geodetic) -> Vec3 {
  let sl = @math.sin(self.lat)
  {
    x: -sl * @math.cos(self.lon),
    y: -sl * @math.sin(self.lon),
    z: @math.cos(self.lat),
  }
}

///|
/// A local Cartesian frame anchored at `origin` (ECEF metres) whose orthonormal axes
/// (columns of `axes`, expressed in ECEF) define the local coordinates.
#valtype
pub(all) struct LocalFrame {
  origin : Vec3
  axes : Mat3
} derive(Eq, Debug)

///|
/// East-North-Up tangent frame at geodetic position `g` (local x = east, y = north, z = up).
pub fn LocalFrame::enu(e : Ellipsoid, g : Geodetic) -> LocalFrame {
  {
    origin: e.to_ecef(g),
    axes: Mat3::from_cols(g.east_ecef(), g.north_ecef(), g.up_ecef()),
  }
}

///|
/// North-East-Down tangent frame (aerospace convention: x = north, y = east, z = down).
pub fn LocalFrame::ned(e : Ellipsoid, g : Geodetic) -> LocalFrame {
  {
    origin: e.to_ecef(g),
    axes: Mat3::from_cols(g.north_ecef(), g.east_ecef(), -g.up_ecef()),
  }
}

///|
/// ECEF -> local frame coordinates.
pub fn LocalFrame::to_local(self : LocalFrame, ecef : Vec3) -> Vec3 {
  self.axes.transpose().mul_vec3(ecef - self.origin)
}

///|
/// Local frame coordinates -> ECEF.
pub fn LocalFrame::to_ecef(self : LocalFrame, loc : Vec3) -> Vec3 {
  self.axes.mul_vec3(loc) + self.origin
}

///|
/// Rotate a free vector (e.g. velocity, direction) from ECEF to local axes.
pub fn LocalFrame::dir_to_local(self : LocalFrame, ecef_dir : Vec3) -> Vec3 {
  self.axes.transpose().mul_vec3(ecef_dir)
}

///|
/// Rigid transform ECEF -> RENDER space for a scene centred on an ENU frame: east -> +X,
/// up -> +Y, north -> -Z (the library's Y-up, -Z-forward world). Use the returned pose's
/// `to_mat4()` as a view/world matrix for geospatial content, or `inverse()` to go back.
pub fn LocalFrame::ecef_to_render_enu(self : LocalFrame) -> Rigid {
  // render = Q * (axes^T (p - origin)),  Q maps (e, n, u) -> (e, u, -n)
  let q = Mat3::from_rows(
    Vec3::new(1.0, 0.0, 0.0),
    Vec3::new(0.0, 0.0, 1.0),
    Vec3::new(0.0, -1.0, 0.0),
  )
  let r = q * self.axes.transpose()
  { rotation: Quat::from_mat3(r), translation: -r.mul_vec3(self.origin), }
}

///|
/// Globe-centred render axes: ECEF (x, y, z) -> render (x, z, -y), i.e. north pole = +Y.
/// A proper rotation (determinant +1).
pub fn ecef_to_yup(p : Vec3) -> Vec3 {
  { x: p.x, y: p.z, z: -p.y, }
}

///|
pub fn yup_to_ecef(p : Vec3) -> Vec3 {
  { x: p.x, y: -p.z, z: p.y, }
}
