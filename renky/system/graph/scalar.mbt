// Scalar helpers shared by every coordinate-system module.
//
// UNITS: every angle in this library is in RADIANS.

///|
pub const PI : Double = 3.141592653589793

///|
pub const TAU : Double = 6.283185307179586

///|
pub const HALF_PI : Double = 1.5707963267948966

///|
pub const DEG2RAD : Double = 0.017453292519943295

///|
pub const RAD2DEG : Double = 57.29577951308232

///|
pub fn deg_to_rad(d : Double) -> Double {
  d * DEG2RAD
}

///|
pub fn rad_to_deg(r : Double) -> Double {
  r * RAD2DEG
}

///|
pub fn clamp(x : Double, lo : Double, hi : Double) -> Double {
  if x < lo {
    lo
  } else if x > hi {
    hi
  } else {
    x
  }
}

///|
pub fn saturate(x : Double) -> Double {
  clamp(x, 0.0, 1.0)
}

///|
pub fn lerp(a : Double, b : Double, t : Double) -> Double {
  a + (b - a) * t
}

///|
/// |a - b| <= eps (absolute).
pub fn approx(a : Double, b : Double, eps : Double) -> Bool {
  (a - b).abs() <= eps
}

///|
/// asin that tolerates |x| slightly above 1 from rounding error.
pub fn safe_asin(x : Double) -> Double {
  @math.asin(clamp(x, -1.0, 1.0))
}

///|
/// acos that tolerates |x| slightly above 1 from rounding error.
pub fn safe_acos(x : Double) -> Double {
  @math.acos(clamp(x, -1.0, 1.0))
}

///|
/// Wrap an angle into (-PI, PI].
pub fn wrap_pi(a : Double) -> Double {
  if a > -PI && a <= PI {
    return a
  }
  let r = a - TAU * @math.floor(a / TAU + 0.5)
  // floor(a/TAU + 0.5) puts r in [-PI, PI); move -PI to +PI for a half-open (-PI, PI].
  if r <= -PI {
    r + TAU
  } else {
    r
  }
}

///|
/// Wrap an angle into [0, TAU).
pub fn wrap_tau(a : Double) -> Double {
  let r = a - TAU * @math.floor(a / TAU)
  if r >= TAU {
    0.0
  } else {
    r
  }
}

///|
/// +1.0 for x >= 0 (including +0 and -0), -1.0 otherwise. Never returns 0.
pub fn sign_nz(x : Double) -> Double {
  if x >= 0.0 {
    1.0
  } else {
    -1.0
  }
}

///|
/// Floor division remainder that is always in [0, n) for n > 0.
pub fn mod_floor(a : Int, n : Int) -> Int {
  let r = a % n
  if r < 0 {
    r + n
  } else {
    r
  }
}
