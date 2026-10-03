// Micro-benchmarks of the hot paths.   moon run cmd/benchmark --target native --release

///|
fn main {
  let b = @bench.Bench()
  let m1 = @graph.Mat4::from_trs(
    @graph.Vec3::new(1.0, 2.0, 3.0),
    @graph.Quat::from_axis_angle(@graph.Vec3::new(1.0, 2.0, 3.0), 0.7),
    @graph.Vec3::new(1.0, 2.0, 0.5),
  )
  let m2 = @graph.Mat4::look_at(
    @graph.Vec3::new(3.0, 4.0, 5.0),
    @graph.Vec3::zero(),
    @graph.Vec3::unit_y(),
  )
  let proj = @graph.Projection::perspective(
    1.0,
    1.6,
    0.1,
    @double.infinity,
    ReversedZ,
  )
  let cam = @graph.Camera::look_at(
    @graph.Vec3::new(3.0, 4.0, 5.0),
    @graph.Vec3::zero(),
    @graph.Vec3::unit_y(),
    proj,
    @graph.Viewport::full(1920.0, 1080.0),
  )
  let q = @graph.Quat::from_axis_angle(@graph.Vec3::new(1.0, 2.0, 3.0), 0.7)
  let p = @graph.Vec3::new(0.3, -0.2, -4.0)
  let g = @graph.Geodetic::new(0.5, 1.0, 1200.0)
  let wgs = @graph.Ellipsoid::wgs84()
  let e = wgs.to_ecef(g)
  b.bench(name="Mat4 * Mat4", fn() { b.keep(m1 * m2) })
  b.bench(name="Mat4::inverse", fn() { b.keep(m1.inverse()) })
  b.bench(name="Mat4::mul_point", fn() { b.keep(m1.mul_point(p)) })
  b.bench(name="Quat::rotate", fn() { b.keep(q.rotate(p)) })
  b.bench(name="Camera::world_to_screen", fn() {
    b.keep(cam.world_to_screen(p))
  })
  b.bench(name="Camera::pick_ray", fn() { b.keep(cam.pick_ray(960.5, 540.5)) })
  b.bench(name="Geodetic -> ECEF", fn() { b.keep(wgs.to_ecef(g)) })
  b.bench(name="ECEF -> Geodetic", fn() { b.keep(wgs.from_ecef(e)) })
  b.bench(name="cube_from_direction", fn() {
    b.keep(@graph.cube_from_direction(p))
  })
  println(b.dump_summaries())
}
