// Shared helpers for all blackbox tests of this package.

///|
/// Mixed absolute / relative tolerance: |got-want| <= tol * (1 + |want|).
fn check_close(
  got : Double,
  want : Double,
  tol : Double,
  what : String,
) -> Unit raise {
  let d = (got - want).abs()
  if !(d <= tol * (1.0 + want.abs())) {
    fail("\{what}: got \{got}, want \{want}, |diff|=\{d}, tol=\{tol}")
  }
}

///|
fn check_v2(
  got : @graph.Vec2,
  want : @graph.Vec2,
  tol : Double,
  what : String,
) -> Unit raise {
  check_close(got.x, want.x, tol, what + ".x")
  check_close(got.y, want.y, tol, what + ".y")
}

///|
fn check_v3(
  got : @graph.Vec3,
  want : @graph.Vec3,
  tol : Double,
  what : String,
) -> Unit raise {
  check_close(got.x, want.x, tol, what + ".x")
  check_close(got.y, want.y, tol, what + ".y")
  check_close(got.z, want.z, tol, what + ".z")
}

///|
fn check_v4(
  got : @graph.Vec4,
  want : @graph.Vec4,
  tol : Double,
  what : String,
) -> Unit raise {
  check_close(got.x, want.x, tol, what + ".x")
  check_close(got.y, want.y, tol, what + ".y")
  check_close(got.z, want.z, tol, what + ".z")
  check_close(got.w, want.w, tol, what + ".w")
}

///|
fn check_m3(
  got : @graph.Mat3,
  want : @graph.Mat3,
  tol : Double,
  what : String,
) -> Unit raise {
  for r in 0..<3 {
    for c in 0..<3 {
      check_close(
        got.row(r).get(c),
        want.row(r).get(c),
        tol,
        "\{what}[\{r}][\{c}]",
      )
    }
  }
}

///|
fn check_m4(
  got : @graph.Mat4,
  want : @graph.Mat4,
  tol : Double,
  what : String,
) -> Unit raise {
  for r in 0..<4 {
    for c in 0..<4 {
      check_close(got.get(r, c), want.get(r, c), tol, "\{what}[\{r}][\{c}]")
    }
  }
}

///|
/// Rotation equality up to q / -q.
fn check_quat(
  got : @graph.Quat,
  want : @graph.Quat,
  tol : Double,
  what : String,
) -> Unit raise {
  let s = if got.dot(want) < 0.0 { -1.0 } else { 1.0 }
  check_close(got.x, want.x * s, tol, what + ".x")
  check_close(got.y, want.y * s, tol, what + ".y")
  check_close(got.z, want.z * s, tol, what + ".z")
  check_close(got.w, want.w * s, tol, what + ".w")
}

///|
fn check_true(c : Bool, what : String) -> Unit raise {
  if !c {
    fail("assertion failed: " + what)
  }
}

///|
/// Deterministic SplitMix64 generator for property tests.
struct Rng {
  mut state : UInt64
}

///|
fn Rng::new(seed : UInt64) -> Rng {
  { state: seed, }
}

///|
fn Rng::next_u64(self : Rng) -> UInt64 {
  self.state = self.state + 0x9E3779B97F4A7C15UL
  let mut z = self.state
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9UL
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBUL
  z ^ (z >> 31)
}

///|
/// Uniform in [0, 1).
fn Rng::unit(self : Rng) -> Double {
  Double::convert_uint64(self.next_u64() >> 11) * 1.1102230246251565e-16
}

///|
fn Rng::range(self : Rng, lo : Double, hi : Double) -> Double {
  lo + (hi - lo) * self.unit()
}

///|
fn Rng::vec3(self : Rng, lo : Double, hi : Double) -> @graph.Vec3 {
  @graph.Vec3::new(self.range(lo, hi), self.range(lo, hi), self.range(lo, hi))
}

///|
/// Uniformly distributed random rotation (Shoemake).
fn Rng::quat(self : Rng) -> @graph.Quat {
  let u1 = self.unit()
  let u2 = self.unit() * @graph.TAU
  let u3 = self.unit() * @graph.TAU
  let a = (1.0 - u1).sqrt()
  let b = u1.sqrt()
  @graph.Quat::new(
    a * @math.sin(u2),
    a * @math.cos(u2),
    b * @math.sin(u3),
    b * @math.cos(u3),
  )
}

///|
fn Rng::unit_vec3(self : Rng) -> @graph.Vec3 {
  let z = self.range(-1.0, 1.0)
  let t = self.range(0.0, @graph.TAU)
  let r = (1.0 - z * z).sqrt()
  @graph.Vec3::new(r * @math.cos(t), r * @math.sin(t), z)
}
