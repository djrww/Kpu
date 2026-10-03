// Software-rendered demo built ONLY from the coordinate-system library: prints a PPM image.
//   moon run cmd/demo --target native > demo.ppm
// What it exercises: look_at / orbit camera, reversed-Z infinite projection, viewport transform,
// frustum culling, TRS model matrices + normal matrix, perspective-correct barycentric
// interpolation, picking rays for the sky (cube-map lookup) and ground plane, axis gizmo.

///|
struct Framebuffer {
  w : Int
  h : Int
  color : FixedArray[Int] // packed 0xRRGGBB
  depth : FixedArray[Double] // reversed-Z: larger = nearer, cleared to 0
}

///|
fn Framebuffer::new(w : Int, h : Int) -> Framebuffer {
  {
    w,
    h,
    color: FixedArray::make(w * h, 0),
    depth: FixedArray::make(w * h, 0.0),
  }
}

///|
fn pack(r : Double, g : Double, b : Double) -> Int {
  let c = fn(x : Double) -> Int {
    (@graph.clamp(x, 0.0, 1.0) * 255.0 + 0.5).to_int()
  }
  (c(r) << 16) | (c(g) << 8) | c(b)
}

///|
fn checker(u : Double, v : Double, n : Double) -> Bool {
  let iu = @math.floor(u * n).to_int()
  let iv = @math.floor(v * n).to_int()
  ((iu + iv) & 1) == 0
}

///|
/// Sky colour from a world direction through the cube-map coordinate system: every face gets
/// its own hue and a UV grid, so seams and orientation are visible.
fn sky(dir : @graph.Vec3) -> Int {
  let c = @graph.cube_from_direction(dir)
  let hues : Array[(Double, Double, Double)] = [
    (0.85, 0.35, 0.35),
    (0.55, 0.20, 0.20),
    (0.35, 0.80, 0.40),
    (0.15, 0.40, 0.20),
    (0.35, 0.45, 0.90),
    (0.18, 0.22, 0.50),
  ]
  let (r, g, b) = hues[c.face]
  let line = c.u * 8.0 - @math.floor(c.u * 8.0) < 0.04 ||
    c.v * 8.0 - @math.floor(c.v * 8.0) < 0.04
  let k = if line { 1.0 } else { 0.55 }
  // darker toward the horizon-less bottom hemisphere
  pack(r * k, g * k, b * k)
}

///|
fn draw_triangle(
  fb : Framebuffer,
  cam : @graph.Camera,
  p : Array[@graph.Vec3], // world positions
  uv : Array[@graph.Vec2],
  normal_world : @graph.Vec3,
  base : (Double, Double, Double),
) -> Unit {
  let scr : Array[@graph.Vec3] = []
  let ws : Array[Double] = []
  for v in p {
    match cam.world_to_screen(v) {
      Some(s) => {
        scr.push(s)
        ws.push(cam.world_to_clip(v).w)
      }
      None => return // behind camera: a real renderer would clip; the demo keeps cubes in front
    }
  }
  let minx = @math.floor(scr[0].x.min(scr[1].x).min(scr[2].x))
    .to_int()
    .clamp(min=0, max=fb.w - 1)
  let maxx = @math.floor(scr[0].x.max(scr[1].x).max(scr[2].x))
    .to_int()
    .clamp(min=0, max=fb.w - 1)
  let miny = @math.floor(scr[0].y.min(scr[1].y).min(scr[2].y))
    .to_int()
    .clamp(min=0, max=fb.h - 1)
  let maxy = @math.floor(scr[0].y.max(scr[1].y).max(scr[2].y))
    .to_int()
    .clamp(min=0, max=fb.h - 1)
  let light = @graph.Vec3::new(-0.4, 0.8, 0.5).normalize()
  let lam = @graph.Vec3::new(0.0, 0.0, 0.0)
  let _ = lam
  let (br, bg, bb) = base
  let shade = 0.25 + 0.75 * @graph.saturate(normal_world.dot(light))
  for y in miny..<=maxy {
    for x in minx..<=maxx {
      let c = @graph.pixel_center(x, y)
      match @graph.barycentric_2d(c, scr[0].xy(), scr[1].xy(), scr[2].xy()) {
        Some(b) =>
          if @graph.bary_inside(b, 0.0) {
            let z = @graph.interpolate1(b, scr[0].z, scr[1].z, scr[2].z)
            let i = y * fb.w + x
            if z > fb.depth[i] {
              let pc = @graph.perspective_correct(b, ws[0], ws[1], ws[2])
              let t = @graph.interpolate2(pc, uv[0], uv[1], uv[2])
              let k = if checker(t.x, t.y, 4.0) { 1.0 } else { 0.55 }
              fb.depth[i] = z
              fb.color[i] = pack(br * shade * k, bg * shade * k, bb * shade * k)
            }
          }
        None => ()
      }
    }
  }
}

///|
/// Draw the six faces of a unit cube [-0.5, 0.5]^3 under model matrix `m`.
fn draw_cube(
  fb : Framebuffer,
  cam : @graph.Camera,
  m : @graph.Mat4,
  base : (Double, Double, Double),
) -> Unit {
  // frustum culling with the world-space AABB of the transformed cube
  let box = @graph.Aabb::from_center_half(
    @graph.Vec3::zero(),
    @graph.Vec3::splat(0.5),
  ).transform(m)
  if cam.frustum().classify_aabb(box) == Outside {
    return
  }
  // normals transform with the inverse transpose (correct under non-uniform scale)
  let nm = m.normal_matrix().unwrap()
  // face: normal, u axis, v axis (right-handed), 4 corners
  let faces : Array[(@graph.Vec3, @graph.Vec3, @graph.Vec3)] = [
    (@graph.Vec3::unit_x(), -@graph.Vec3::unit_z(), @graph.Vec3::unit_y()),
    (-@graph.Vec3::unit_x(), @graph.Vec3::unit_z(), @graph.Vec3::unit_y()),
    (@graph.Vec3::unit_y(), @graph.Vec3::unit_x(), -@graph.Vec3::unit_z()),
    (-@graph.Vec3::unit_y(), @graph.Vec3::unit_x(), @graph.Vec3::unit_z()),
    (@graph.Vec3::unit_z(), @graph.Vec3::unit_x(), @graph.Vec3::unit_y()),
    (-@graph.Vec3::unit_z(), -@graph.Vec3::unit_x(), @graph.Vec3::unit_y()),
  ]
  for (n, u, v) in faces {
    let c = n.scale(0.5)
    let q = [
      c - u.scale(0.5) - v.scale(0.5),
      c + u.scale(0.5) - v.scale(0.5),
      c + u.scale(0.5) + v.scale(0.5),
      c - u.scale(0.5) + v.scale(0.5),
    ]
    let wq = q.map(fn(p) { m.mul_point(p) })
    let nw = nm.mul_vec3(n).normalize()
    // back-face culling in world space
    if nw.dot(wq[0] - cam.eye()) >= 0.0 {
      continue
    }
    let uvs = [
      @graph.Vec2::new(0.0, 0.0),
      @graph.Vec2::new(1.0, 0.0),
      @graph.Vec2::new(1.0, 1.0),
      @graph.Vec2::new(0.0, 1.0),
    ]
    draw_triangle(
      fb,
      cam,
      [wq[0], wq[1], wq[2]],
      [uvs[0], uvs[1], uvs[2]],
      nw,
      base,
    )
    draw_triangle(
      fb,
      cam,
      [wq[0], wq[2], wq[3]],
      [uvs[0], uvs[2], uvs[3]],
      nw,
      base,
    )
  }
}

///|
fn draw_line(
  fb : Framebuffer,
  a : @graph.Vec3,
  b : @graph.Vec3,
  color : Int,
) -> Unit {
  let steps = @math.hypot(b.x - a.x, b.y - a.y).to_int() + 1
  for i in 0..<=steps {
    let t = Double::from_int(i) / Double::from_int(steps)
    let x = @math.floor(a.x + (b.x - a.x) * t).to_int()
    let y = @math.floor(a.y + (b.y - a.y) * t).to_int()
    if x >= 0 && x < fb.w && y >= 0 && y < fb.h {
      fb.color[y * fb.w + x] = color
    }
  }
}

///|
fn main {
  let w = 480
  let h = 300
  let fb = Framebuffer::new(w, h)
  let proj = @graph.Projection::perspective(
    1.0,
    Double::from_int(w) / Double::from_int(h),
    0.1,
    @double.infinity,
    ReversedZ,
  )
  let cam = @graph.Camera::orbit(
    @graph.Vec3::new(0.0, 0.6, 0.0),
    7.5,
    -0.55, // yaw: look toward +X-ish (negative yaw turns right)
    -0.28, // pitch: look slightly down
    proj,
    @graph.Viewport::full(Double::from_int(w), Double::from_int(h)),
  )
  // 1. background + ground by picking rays through every pixel centre
  let ground = @graph.Plane::from_point_normal(
    @graph.Vec3::zero(),
    @graph.Vec3::unit_y(),
  )
  for y in 0..<h {
    for x in 0..<w {
      let c = @graph.pixel_center(x, y)
      let ray = cam.pick_ray(c.x, c.y)
      match ray.intersect_plane(ground) {
        Some(t) => {
          let hit = ray.at(t)
          let fade = 1.0 / (1.0 + t * 0.03)
          let k = if checker(hit.x * 0.5, hit.z * 0.5, 1.0) {
            0.75
          } else {
            0.45
          }
          fb.color[y * w + x] = pack(
            0.55 * k * fade + 0.1,
            0.6 * k * fade + 0.12,
            0.55 * k * fade + 0.15,
          )
          // write ground depth so cubes are depth-tested against it
          let s = cam.world_to_screen(hit).unwrap()
          fb.depth[y * w + x] = s.z
        }
        None => fb.color[y * w + x] = sky(ray.dir)
      }
    }
  }
  // 2. cubes with different TRS (one mirrored scale to prove the normal matrix is right)
  let cubes : Array[(@graph.Mat4, (Double, Double, Double))] = [
    (
      @graph.Mat4::from_trs(
        @graph.Vec3::new(0.0, 0.5, 0.0),
        @graph.Quat::from_y(0.5),
        @graph.Vec3::splat(1.0),
      ),
      (0.95, 0.85, 0.25),
    ),
    (
      @graph.Mat4::from_trs(
        @graph.Vec3::new(2.4, 0.4, -1.5),
        @graph.Quat::from_axis_angle(@graph.Vec3::new(1.0, 1.0, 0.0), 0.6),
        @graph.Vec3::new(0.8, 0.8, 0.8),
      ),
      (0.9, 0.35, 0.35),
    ),
    (
      @graph.Mat4::from_trs(
        @graph.Vec3::new(-2.2, 1.0, -1.0),
        @graph.Quat::from_y(-0.3),
        @graph.Vec3::new(1.0, 2.0, 1.0),
      ),
      (0.35, 0.65, 0.95),
    ),
    (
      @graph.Mat4::from_trs(
        @graph.Vec3::new(1.2, 0.3, 2.0),
        @graph.Quat::from_y(0.9),
        @graph.Vec3::new(0.6, 0.6, 0.6),
      ),
      (0.5, 0.9, 0.5),
    ),
    (
      @graph.Mat4::from_trs(
        @graph.Vec3::new(0.0, 0.5, -5.0),
        @graph.Quat::identity(),
        @graph.Vec3::new(1.0, 1.0, 1.0),
      ),
      (0.9, 0.9, 0.9),
    ),
    // behind the camera / outside the view: must be frustum-culled
    (
      @graph.Mat4::translation(@graph.Vec3::new(0.0, 0.5, 30.0)),
      (1.0, 0.0, 1.0),
    ),
  ]
  for (m, col) in cubes {
    draw_cube(fb, cam, m, col)
  }
  // 3. world axis gizmo from the origin: +X red, +Y green, +Z blue
  let o = cam.world_to_screen(@graph.Vec3::zero()).unwrap()
  for
    (axis, col) in [
      (@graph.Vec3::unit_x(), 0xFF2020),
      (@graph.Vec3::unit_y(), 0x20FF20),
      (@graph.Vec3::unit_z(), 0x3030FF),
    ] {
    match cam.world_to_screen(axis.scale(2.5)) {
      Some(e) => draw_line(fb, o, e, col)
      None => ()
    }
  }
  // 4. output PPM (P3)
  let sb = StringBuilder()
  sb.write_string("P3\n\{w} \{h}\n255\n")
  for i in 0..<(w * h) {
    let c = fb.color[i]
    sb.write_string("\{(c >> 16) & 255} \{(c >> 8) & 255} \{c & 255}\n")
  }
  println(sb.to_string())
}
