// Viewport: NDC <-> window ("screen" / framebuffer pixel) coordinates.
//
// Window space: origin at the TOP-LEFT, +x right, +y DOWN, units = pixels. Pixel (i, j) covers
// [i, i+1) x [j, j+1) and its centre is (i + 0.5, j + 0.5) (D3D10+/Vulkan/Metal/WebGPU rule).
// Window depth = min_depth + unit_depth * (max_depth - min_depth).

///|
#valtype
pub(all) struct Viewport {
  x : Double
  y : Double
  width : Double
  height : Double
  min_depth : Double
  max_depth : Double
} derive(Eq, Debug)

///|
pub fn Viewport::new(
  x : Double,
  y : Double,
  width : Double,
  height : Double,
) -> Viewport {
  if !(width > 0.0 && height > 0.0) {
    abort("Viewport::new: width and height must be > 0")
  }
  { x, y, width, height, min_depth: 0.0, max_depth: 1.0, }
}

///|
pub fn Viewport::full(width : Double, height : Double) -> Viewport {
  Viewport::new(0.0, 0.0, width, height)
}

///|
/// Viewport for a window given in LOGICAL (DPI-independent) units; the pixel size is the
/// logical size times `scale`, rounded to whole device pixels.
pub fn Viewport::from_logical(
  logical_width : Double,
  logical_height : Double,
  scale : Double,
) -> Viewport {
  Viewport::full(
    @math.round(logical_width * scale),
    @math.round(logical_height * scale),
  )
}

///|
pub fn Viewport::aspect(self : Viewport) -> Double {
  self.width / self.height
}

///|
/// NDC (x, y, z) -> window coordinates (x px, y px, depth). `mode` tells how ndc z maps to
/// [0,1] before min/max depth are applied.
pub fn Viewport::ndc_to_screen(
  self : Viewport,
  ndc : Vec3,
  mode : DepthMode,
) -> Vec3 {
  let u = mode.to_unit(ndc.z)
  {
    x: self.x + (ndc.x * 0.5 + 0.5) * self.width,
    y: self.y + (0.5 - ndc.y * 0.5) * self.height,
    z: self.min_depth + u * (self.max_depth - self.min_depth),
  }
}

///|
/// Window coordinates + window depth -> NDC.
pub fn Viewport::screen_to_ndc(
  self : Viewport,
  sx : Double,
  sy : Double,
  depth : Double,
  mode : DepthMode,
) -> Vec3 {
  let u = (depth - self.min_depth) / (self.max_depth - self.min_depth)
  {
    x: (sx - self.x) / self.width * 2.0 - 1.0,
    y: 1.0 - (sy - self.y) / self.height * 2.0,
    z: mode.from_unit(u),
  }
}

///|
/// Window xy -> NDC xy only.
pub fn Viewport::screen_to_ndc_xy(self : Viewport, p : Vec2) -> Vec2 {
  {
    x: (p.x - self.x) / self.width * 2.0 - 1.0,
    y: 1.0 - (p.y - self.y) / self.height * 2.0,
  }
}

///|
pub fn Viewport::ndc_to_screen_xy(self : Viewport, ndc : Vec2) -> Vec2 {
  {
    x: self.x + (ndc.x * 0.5 + 0.5) * self.width,
    y: self.y + (0.5 - ndc.y * 0.5) * self.height,
  }
}

///|
/// Window position -> normalised [0,1]^2 coordinates inside the viewport (origin top-left;
/// the same layout as texture UVs of a full-screen image).
pub fn Viewport::to_uv(self : Viewport, p : Vec2) -> Vec2 {
  { x: (p.x - self.x) / self.width, y: (p.y - self.y) / self.height, }
}

///|
pub fn Viewport::from_uv(self : Viewport, uv : Vec2) -> Vec2 {
  { x: self.x + uv.x * self.width, y: self.y + uv.y * self.height, }
}

///|
/// Does the window position fall inside the viewport rectangle?
pub fn Viewport::contains(self : Viewport, p : Vec2) -> Bool {
  p.x >= self.x &&
  p.x < self.x + self.width &&
  p.y >= self.y &&
  p.y < self.y + self.height
}

///|
/// Centre of pixel (ix, iy) in window coordinates.
pub fn pixel_center(ix : Int, iy : Int) -> Vec2 {
  { x: Double::from_int(ix) + 0.5, y: Double::from_int(iy) + 0.5, }
}

///|
/// Integer pixel containing a window position.
pub fn pixel_of(p : Vec2) -> (Int, Int) {
  (@math.floor(p.x).to_int(), @math.floor(p.y).to_int())
}

///|
pub fn logical_to_physical(p : Vec2, scale : Double) -> Vec2 {
  p.scale(scale)
}

///|
pub fn physical_to_logical(p : Vec2, scale : Double) -> Vec2 {
  p.scale(1.0 / scale)
}
