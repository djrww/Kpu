// END-TO-END: a tiny reference rasteriser (test code only) built purely from the library's
// coordinate systems is compared, pixel by pixel, against ANALYTIC ray/plane intersections.
// This is what proves the spaces compose: world -> view -> clip -> NDC -> window and back,
// barycentrics, perspective-correct interpolation, depth modes and picking rays.

///|
struct Quad {
  o : @graph.Vec3
  u : @graph.Vec3 // edge along s
  v : @graph.Vec3 // edge along t
}

///|
fn Quad::corner(self : Quad, s : Double, t : Double) -> @graph.Vec3 {
  self.o + self.u.scale(s) + self.v.scale(t)
}

///|
/// Analytic hit of a camera ray with the quad: (distance t along ray, s, t) or None.
fn Quad::analytic(self : Quad, ray : @graph.Ray) -> (Double, Double, Double)? {
  let pl = @graph.Plane::from_point_normal(self.o, self.u.cross(self.v))
  let denom = pl.normal.dot(ray.dir)
  if denom == 0.0 {
    return None
  }
  let t = -(pl.normal.dot(ray.origin) + pl.d) / denom
  if t < 0.0 {
    return None
  }
  let h = ray.at(t)
  // solve  d = u*s + v*r  (u, v need not be orthogonal) with the 2x2 Gram system
  let d = h - self.o
  let uu = self.u.dot(self.u)
  let uv = self.u.dot(self.v)
  let vv = self.v.dot(self.v)
  let det = uu * vv - uv * uv
  let s = (d.dot(self.u) * vv - d.dot(self.v) * uv) / det
  let r = (d.dot(self.v) * uu - d.dot(self.u) * uv) / det
  Some((t, s, r))
}

///|
/// Rasterise `quad` (two triangles) into callbacks-free arrays: returns per-pixel
/// (covered, world position, uv, window depth); everything computed through the library.
fn raster_quad(
  cam : @graph.Camera,
  q : Quad,
  wpx : Int,
  hpx : Int,
) -> Array[(Bool, @graph.Vec3, @graph.Vec2, Double)] {
  let corners = [
    q.corner(0.0, 0.0),
    q.corner(1.0, 0.0),
    q.corner(1.0, 1.0),
    q.corner(0.0, 1.0),
  ]
  let uvs = [
    @graph.Vec2::new(0.0, 0.0),
    @graph.Vec2::new(1.0, 0.0),
    @graph.Vec2::new(1.0, 1.0),
    @graph.Vec2::new(0.0, 1.0),
  ]
  let tris = [(0, 1, 2), (0, 2, 3)]
  let out : Array[(Bool, @graph.Vec3, @graph.Vec2, Double)] = []
  // pre-transform vertices once
  let scr : Array[@graph.Vec3] = []
  let ws : Array[Double] = []
  for c in corners {
    scr.push(cam.world_to_screen(c).unwrap())
    ws.push(cam.world_to_clip(c).w)
  }
  for py in 0..<hpx {
    for px in 0..<wpx {
      let p = @graph.pixel_center(px, py)
      let mut hit = (false, @graph.Vec3::zero(), @graph.Vec2::zero(), 0.0)
      for (i0, i1, i2) in tris {
        if hit.0 {
          continue
        }
        match
          @graph.barycentric_2d(p, scr[i0].xy(), scr[i1].xy(), scr[i2].xy()) {
          Some(b) =>
            if @graph.bary_inside(b, 0.0) {
              let pc = @graph.perspective_correct(b, ws[i0], ws[i1], ws[i2])
              hit = (
                true,
                @graph.interpolate3(pc, corners[i0], corners[i1], corners[i2]),
                @graph.interpolate2(pc, uvs[i0], uvs[i1], uvs[i2]),
                // window depth is affine in screen space: interpolate with the 2D weights
                @graph.interpolate1(b, scr[i0].z, scr[i1].z, scr[i2].z),
              )
            }
          None => ()
        }
      }
      out.push(hit)
    }
  }
  out
}

///|
fn run_scene(
  cam : @graph.Camera,
  q : Quad,
  w : Int,
  h : Int,
  tag : String,
) -> Int raise {
  let img = raster_quad(cam, q, w, h)
  let mut covered = 0
  let mut checked = 0
  for py in 0..<h {
    for px in 0..<w {
      let (cov, wp, uv, depth) = img[py * w + px]
      let c = @graph.pixel_center(px, py)
      let ray = cam.pick_ray(c.x, c.y)
      match q.analytic(ray) {
        Some((t, s, r)) => {
          // skip pixels whose centre is within 1e-7 of a quad edge or a triangle diagonal
          let margin = [s, 1.0 - s, r, 1.0 - r, (s - r).abs()]
          let mut m = margin[0]
          for x in margin {
            m = m.min(x)
          }
          if m < 1.0e-7 && m > -1.0e-7 {
            continue
          }
          let inside = s > 0.0 && s < 1.0 && r > 0.0 && r < 1.0
          check_true(
            cov == inside,
            tag + ": coverage disagrees with analytic ray test at (\{px},\{py})",
          )
          if inside {
            covered += 1
            checked += 1
            check_v3(
              wp,
              ray.at(t),
              1.0e-8,
              tag + ": interpolated world position vs ray hit",
            )
            check_close(uv.x, s, 1.0e-9, tag + ": u")
            check_close(uv.y, r, 1.0e-9, tag + ": v")
            let sd = cam.world_to_screen(ray.at(t)).unwrap().z
            check_close(depth, sd, 1.0e-9, tag + ": window depth")
            // depth buffer value converts back to the true distance in front of the camera
            let ndc_z = cam.projection.depth.from_unit(depth)
            let lin = cam.projection.linear_depth(ndc_z)
            check_close(
              lin,
              (ray.at(t) - cam.eye()).dot(cam.forward()),
              1.0e-6,
              tag + ": linear depth",
            )
          }
        }
        None => check_true(!cov, tag + ": covered but ray misses plane")
      }
    }
  }
  check_true(
    checked > w * h / 10,
    tag + ": enough pixels covered to be meaningful (\{checked})",
  )
  covered
}

///|
test "reference rasteriser == analytic ray casting (perspective: all depth modes, infinite far)" {
  let q = {
    o: @graph.Vec3::new(-3.0, -1.0, -4.0),
    u: @graph.Vec3::new(6.0, 0.0, -2.0),
    v: @graph.Vec3::new(0.0, 3.0, -8.0),
  }
  let mut n = 0
  for m in all_modes() {
    for far in [100.0, @double.infinity] {
      let cam = @graph.Camera::look_at(
        @graph.Vec3::new(0.5, 0.3, -0.5),
        @graph.Vec3::new(0.0, 0.3, -8.0),
        @graph.Vec3::unit_y(),
        @graph.Projection::perspective(1.3, 4.0 / 3.0, 0.1, far, m),
        @graph.Viewport::full(80.0, 60.0),
      )
      n += run_scene(cam, q, 80, 60, "persp mode \{n}")
    }
  }
  check_true(n > 0, "something rendered")
}

///|
test "reference rasteriser == analytic ray casting (orthographic, off-centre frustum, sub-viewport)" {
  let q = {
    o: @graph.Vec3::new(-2.0, -2.0, -5.0),
    u: @graph.Vec3::new(4.0, 0.0, 1.0),
    v: @graph.Vec3::new(0.0, 4.0, 0.0),
  }
  let cam_o = @graph.Camera::look_at(
    @graph.Vec3::new(0.0, 0.0, 5.0),
    @graph.Vec3::zero(),
    @graph.Vec3::unit_y(),
    @graph.Projection::orthographic(
      -3.0,
      3.0,
      -2.0,
      2.0,
      0.5,
      30.0,
      NegOneToOne,
    ),
    @graph.Viewport::full(60.0, 40.0),
  )
  let _ = run_scene(cam_o, q, 60, 40, "ortho")
  // lens-shifted perspective (asymmetric frustum)
  let cam_a = @graph.Camera::look_at(
    @graph.Vec3::new(0.0, 0.0, 6.0),
    @graph.Vec3::zero(),
    @graph.Vec3::unit_y(),
    @graph.Projection::frustum(-0.2, 0.5, -0.3, 0.1, 0.5, 50.0, ReversedZ),
    @graph.Viewport::full(70.0, 50.0),
  )
  let _ = run_scene(cam_a, q, 70, 50, "asymmetric")
}

///|
test "z-buffer: nearest surface wins in every depth mode (ties to analytic distances)" {
  // two crossing quads at different depths, the closer one varies per pixel
  let qa : Quad = {
    o: @graph.Vec3::new(-4.0, -3.0, -6.0),
    u: @graph.Vec3::new(8.0, 0.0, -3.0),
    v: @graph.Vec3::new(0.0, 6.0, 0.0),
  }
  let qb : Quad = {
    o: @graph.Vec3::new(-3.0, -3.0, -12.0),
    u: @graph.Vec3::new(6.0, 0.0, 9.0),
    v: @graph.Vec3::new(0.0, 6.0, 0.0),
  }
  let w = 64
  let h = 48
  for m in all_modes() {
    let cam = @graph.Camera::look_at(
      @graph.Vec3::new(0.0, 0.0, 4.0),
      @graph.Vec3::new(0.0, 0.0, -8.0),
      @graph.Vec3::unit_y(),
      @graph.Projection::perspective(1.0, 4.0 / 3.0, 0.1, @double.infinity, m),
      @graph.Viewport::full(Double::from_int(w), Double::from_int(h)),
    )
    let ia = raster_quad(cam, qa, w, h)
    let ib = raster_quad(cam, qb, w, h)
    let mut a_wins = 0
    let mut b_wins = 0
    for py in 0..<h {
      for px in 0..<w {
        let (ca, _, _, da) = ia[py * w + px]
        let (cb, _, _, db) = ib[py * w + px]
        if ca && cb {
          // depth test as the GPU would do it
          let a_nearer_by_buffer = match m {
            ReversedZ => da > db
            _ => da < db
          }
          let c = @graph.pixel_center(px, py)
          let ray = cam.pick_ray(c.x, c.y)
          match (qa.analytic(ray), qb.analytic(ray)) {
            (Some((ta, _, _)), Some((tb, _, _))) =>
              if (ta - tb).abs() > 1.0e-6 {
                check_true(
                  a_nearer_by_buffer == (ta < tb),
                  "z-test picks the geometrically nearer surface (mode \{m == ReversedZ})",
                )
                if ta < tb {
                  a_wins += 1
                } else {
                  b_wins += 1
                }
              }
            _ => fail("rasteriser and analytic disagree on overlap")
          }
        }
      }
    }
    check_true(
      a_wins > 50 && b_wins > 50,
      "both surfaces win somewhere (a=\{a_wins}, b=\{b_wins})",
    )
  }
}

///|
test "model -> world -> view -> clip -> ndc -> window chain, step by step" {
  // Every intermediate space is checked against a hand calculation.
  let model = @graph.Mat4::from_trs(
    @graph.Vec3::new(10.0, 0.0, -20.0),
    @graph.Quat::from_y(@graph.HALF_PI),
    @graph.Vec3::splat(2.0),
  )
  let loc = @graph.Vec3::new(1.0, 0.5, 0.0)
  // scale 2 -> (2, 1, 0); rotate +90deg about Y: x->-z => (0, 1, -2); translate -> (10, 1, -22)
  let world = model.mul_point(loc)
  check_v3(world, @graph.Vec3::new(10.0, 1.0, -22.0), 1.0e-14, "world")
  let cam = @graph.Camera::look_at(
    @graph.Vec3::new(10.0, 1.0, 0.0),
    world,
    @graph.Vec3::unit_y(),
    @graph.Projection::perspective(@graph.HALF_PI, 2.0, 1.0, 101.0, ZeroToOne),
    @graph.Viewport::full(200.0, 100.0),
  )
  // camera looks straight down -Z at the point 22 m away
  let view = cam.world_to_view(world)
  check_v3(view, @graph.Vec3::new(0.0, 0.0, -22.0), 1.0e-13, "view")
  let clip = cam.world_to_clip(world)
  // fov 90 (f=1), aspect 2: clip.x = 0.5*0, clip.y = 0, clip.z = (-101/100)(-22) - 101/100, w = 22
  check_v4(
    clip,
    @graph.Vec4::new(0.0, 0.0, 1.01 * 22.0 - 1.01, 22.0),
    1.0e-12,
    "clip",
  )
  let ndc = cam.world_to_ndc(world).unwrap()
  check_v3(
    ndc,
    @graph.Vec3::new(0.0, 0.0, (1.01 * 22.0 - 1.01) / 22.0),
    1.0e-12,
    "ndc",
  )
  let scr = cam.world_to_screen(world).unwrap()
  check_v3(scr, @graph.Vec3::new(100.0, 50.0, ndc.z), 1.0e-12, "window")
  // a point 22 m ahead and 22 m to the right is at the right edge (tan(45deg) * 22 = 22)
  let edge = cam
    .world_to_screen(@graph.Vec3::new(10.0 + 22.0 * 2.0, 1.0, -22.0))
    .unwrap()
  check_close(
    edge.x,
    200.0,
    1.0e-10,
    "right edge (aspect 2 => half-width = 2 * half-height)",
  )
  let top = cam
    .world_to_screen(@graph.Vec3::new(10.0, 1.0 + 22.0, -22.0))
    .unwrap()
  check_close(top.y, 0.0, 1.0e-10, "top edge: window y grows downward")
  // everything inverts
  check_v3(cam.screen_to_world(scr.x, scr.y, scr.z), world, 1.0e-9, "inverse")
}
