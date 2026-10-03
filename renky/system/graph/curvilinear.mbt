// Curvilinear coordinates: polar, cylindrical, spherical.
//
// All three are expressed in the library's world frame (right-handed, +Y up, -Z forward) with
// ONE azimuth convention so they compose with rotations:
//     azimuth  = yaw angle about +Y by the right-hand rule, measured from the -Z (forward)
//                axis: azimuth 0 -> -Z, +PI/2 -> -X, PI -> +Z, -PI/2 -> +X.
//                (i.e. direction = Quat::from_y(azimuth).rotate(-Z))
//     elevation = angle above the XZ ground plane, in [-PI/2, PI/2] (+PI/2 = +Y).

///|
/// 2D polar coordinates: `angle` counter-clockwise from +X (y-up frame).
#valtype
pub(all) struct Polar {
  radius : Double
  angle : Double
} derive(Eq, Debug)

///|
pub fn Polar::new(radius : Double, angle : Double) -> Polar {
  { radius, angle, }
}

///|
pub fn Polar::to_cartesian(self : Polar) -> Vec2 {
  {
    x: self.radius * @math.cos(self.angle),
    y: self.radius * @math.sin(self.angle),
  }
}

///|
/// Radius >= 0, angle in (-PI, PI]; the origin maps to angle 0.
pub fn Polar::from_cartesian(p : Vec2) -> Polar {
  {
    radius: p.length(),
    angle: if p.x == 0.0 && p.y == 0.0 {
      0.0
    } else {
      @math.atan2(p.y, p.x)
    },
  }
}

///|
/// Cylindrical coordinates about the +Y axis.
#valtype
pub(all) struct Cylindrical {
  radius : Double
  azimuth : Double
  height : Double
} derive(Eq, Debug)

///|
pub fn Cylindrical::new(
  radius : Double,
  azimuth : Double,
  height : Double,
) -> Cylindrical {
  { radius, azimuth, height, }
}

///|
pub fn Cylindrical::to_cartesian(self : Cylindrical) -> Vec3 {
  {
    x: -self.radius * @math.sin(self.azimuth),
    y: self.height,
    z: -self.radius * @math.cos(self.azimuth),
  }
}

///|
pub fn Cylindrical::from_cartesian(p : Vec3) -> Cylindrical {
  let r = @math.hypot(p.x, p.z)
  {
    radius: r,
    azimuth: if r == 0.0 {
      0.0
    } else {
      @math.atan2(-p.x, -p.z)
    },
    height: p.y,
  }
}

///|
/// Spherical coordinates (radius, azimuth, elevation) -- see the module comment.
#valtype
pub(all) struct Spherical {
  radius : Double
  azimuth : Double
  elevation : Double
} derive(Eq, Debug)

///|
pub fn Spherical::new(
  radius : Double,
  azimuth : Double,
  elevation : Double,
) -> Spherical {
  { radius, azimuth, elevation, }
}

///|
pub fn Spherical::to_cartesian(self : Spherical) -> Vec3 {
  let ce = @math.cos(self.elevation)
  {
    x: -self.radius * @math.sin(self.azimuth) * ce,
    y: self.radius * @math.sin(self.elevation),
    z: -self.radius * @math.cos(self.azimuth) * ce,
  }
}

///|
/// Radius >= 0, azimuth in (-PI, PI], elevation in [-PI/2, PI/2]. At the poles (and the
/// origin) azimuth is reported as 0; at the origin elevation is also 0.
pub fn Spherical::from_cartesian(p : Vec3) -> Spherical {
  let r = p.length()
  let h = @math.hypot(p.x, p.z)
  {
    radius: r,
    azimuth: if h == 0.0 {
      0.0
    } else {
      @math.atan2(-p.x, -p.z)
    },
    elevation: if r == 0.0 {
      0.0
    } else {
      @math.atan2(p.y, h)
    },
  }
}

///|
/// Unit direction for yaw (azimuth) and pitch (elevation).
pub fn Vec3::from_yaw_pitch(yaw : Double, pitch : Double) -> Vec3 {
  Spherical::new(1.0, yaw, pitch).to_cartesian()
}

///|
/// (yaw, pitch) of a direction (need not be unit).
pub fn Vec3::to_yaw_pitch(self : Vec3) -> (Double, Double) {
  let s = Spherical::from_cartesian(self)
  (s.azimuth, s.elevation)
}
