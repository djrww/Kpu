// Texture-side coordinate systems.
//
// UV: origin at the TOP-LEFT of the image, u -> right, v -> DOWN (WebGPU / Vulkan / Metal / D3D
// / glTF layout). Texel (i, j) covers uv [i/w, (i+1)/w) x [j/h, (j+1)/h); its centre is
// ((i + 0.5)/w, (j + 0.5)/h). Use `flip_v` for OpenGL's bottom-left origin images.

///|
/// UV -> continuous texel coordinates (texel i spans [i, i+1)). Pixel centre = i + 0.5.
pub fn uv_to_texel(uv : Vec2, width : Int, height : Int) -> Vec2 {
  { x: uv.x * Double::from_int(width), y: uv.y * Double::from_int(height), }
}

///|
pub fn texel_to_uv(t : Vec2, width : Int, height : Int) -> Vec2 {
  { x: t.x / Double::from_int(width), y: t.y / Double::from_int(height), }
}

///|
/// UV of the CENTRE of texel (i, j).
pub fn texel_center_uv(i : Int, j : Int, width : Int, height : Int) -> Vec2 {
  {
    x: (Double::from_int(i) + 0.5) / Double::from_int(width),
    y: (Double::from_int(j) + 0.5) / Double::from_int(height),
  }
}

///|
/// Switch between top-left (v down) and bottom-left (v up) conventions.
pub fn flip_v(uv : Vec2) -> Vec2 {
  { x: uv.x, y: 1.0 - uv.y, }
}

///|
pub(all) enum WrapMode {
  Clamp
  Repeat
  MirroredRepeat
} derive(Eq, Debug)

///|
/// Wrap a normalised coordinate into [0, 1].
pub fn wrap_coord(u : Double, mode : WrapMode) -> Double {
  match mode {
    Clamp => clamp(u, 0.0, 1.0)
    Repeat => u - @math.floor(u)
    MirroredRepeat => {
      let t = u - 2.0 * @math.floor(u * 0.5)
      if t > 1.0 {
        2.0 - t
      } else {
        t
      }
    }
  }
}

///|
/// Wrap an integer texel index into [0, n).
pub fn wrap_texel(i : Int, n : Int, mode : WrapMode) -> Int {
  match mode {
    Clamp => i.clamp(min=0, max=n - 1)
    Repeat => mod_floor(i, n)
    MirroredRepeat => {
      let m = mod_floor(i, 2 * n)
      if m >= n {
        2 * n - 1 - m
      } else {
        m
      }
    }
  }
}

///|
/// The 2x2 texel neighbourhood and weights of a bilinear fetch at `uv`.
/// Result texel = (1-fx)(1-fy) T[x0,y0] + fx(1-fy) T[x1,y0] + (1-fx)fy T[x0,y1] + fx fy T[x1,y1].
#valtype
pub(all) struct Bilinear {
  x0 : Int
  y0 : Int
  x1 : Int
  y1 : Int
  fx : Double
  fy : Double
} derive(Eq, Debug)

///|
pub fn bilinear_footprint(
  uv : Vec2,
  width : Int,
  height : Int,
  wrap_u : WrapMode,
  wrap_v : WrapMode,
) -> Bilinear {
  let x = uv.x * Double::from_int(width) - 0.5
  let y = uv.y * Double::from_int(height) - 0.5
  let fxl = @math.floor(x)
  let fyl = @math.floor(y)
  let ix = fxl.to_int()
  let iy = fyl.to_int()
  {
    x0: wrap_texel(ix, width, wrap_u),
    y0: wrap_texel(iy, height, wrap_v),
    x1: wrap_texel(ix + 1, width, wrap_u),
    y1: wrap_texel(iy + 1, height, wrap_v),
    fx: x - fxl,
    fy: y - fyl,
  }
}

///|
/// Size of a mip level (>= 1).
pub fn mip_extent(size : Int, level : Int) -> Int {
  let s = size >> level
  if s < 1 {
    1
  } else {
    s
  }
}

///|
/// Number of mip levels of a full chain down to 1x1.
pub fn mip_level_count(width : Int, height : Int) -> Int {
  let mut m = if width > height { width } else { height }
  let mut n = 1
  while m > 1 {
    m = m >> 1
    n += 1
  }
  n
}

///|
/// Level of detail from screen-space UV derivatives (what `textureSample` computes):
/// lod = log2(max(|d(uv * size)/dx|, |d(uv * size)/dy|)). May be negative (magnification).
pub fn mip_lod(
  duv_dx : Vec2,
  duv_dy : Vec2,
  width : Int,
  height : Int,
) -> Double {
  let sz : Vec2 = { x: Double::from_int(width), y: Double::from_int(height), }
  let lx = (duv_dx * sz).length()
  let ly = (duv_dy * sz).length()
  let m = lx.max(ly)
  if m <= 0.0 {
    -@double.infinity
  } else {
    @math.ln(m) / @math.ln(2.0)
  }
}

// ---------------------------------------------------------------------------
// Cube maps
// ---------------------------------------------------------------------------

///|
/// Cube face indices (same order as OpenGL / D3D / Vulkan array layers).
pub const CUBE_POS_X : Int = 0

///|
pub const CUBE_NEG_X : Int = 1

///|
pub const CUBE_POS_Y : Int = 2

///|
pub const CUBE_NEG_Y : Int = 3

///|
pub const CUBE_POS_Z : Int = 4

///|
pub const CUBE_NEG_Z : Int = 5

///|
/// Face + coordinates inside that face image (u right, v down, both in [0, 1]).
#valtype
pub(all) struct CubeCoord {
  face : Int
  u : Double
  v : Double
} derive(Eq, Debug)

///|
/// Direction -> (face, u, v) following the OpenGL cube-map selection table (also used by
/// D3D and Vulkan). `dir` need not be unit; zero gives face 0 centre. Ties between equal
/// major axes resolve X, then Y, then Z.
pub fn cube_from_direction(dir : Vec3) -> CubeCoord {
  let ax = dir.x.abs()
  let ay = dir.y.abs()
  let az = dir.z.abs()
  let mut face = 0
  let mut sc = 0.0
  let mut tc = 0.0
  let mut ma = 0.0
  if ax >= ay && ax >= az {
    ma = ax
    if dir.x >= 0.0 {
      face = CUBE_POS_X
      sc = -dir.z
      tc = -dir.y
    } else {
      face = CUBE_NEG_X
      sc = dir.z
      tc = -dir.y
    }
  } else if ay >= az {
    ma = ay
    if dir.y >= 0.0 {
      face = CUBE_POS_Y
      sc = dir.x
      tc = dir.z
    } else {
      face = CUBE_NEG_Y
      sc = dir.x
      tc = -dir.z
    }
  } else {
    ma = az
    if dir.z >= 0.0 {
      face = CUBE_POS_Z
      sc = dir.x
      tc = -dir.y
    } else {
      face = CUBE_NEG_Z
      sc = -dir.x
      tc = -dir.y
    }
  }
  if ma == 0.0 {
    return { face: 0, u: 0.5, v: 0.5, }
  }
  { face, u: (sc / ma + 1.0) * 0.5, v: (tc / ma + 1.0) * 0.5, }
}

///|
/// (face, u, v) -> UNIT direction through that point of the face (exact inverse of
/// `cube_from_direction`).
pub fn cube_to_direction(face : Int, u : Double, v : Double) -> Vec3 {
  let sc = 2.0 * u - 1.0
  let tc = 2.0 * v - 1.0
  let d : Vec3 = match face {
    0 => { x: 1.0, y: -tc, z: -sc, }
    1 => { x: -1.0, y: -tc, z: sc, }
    2 => { x: sc, y: 1.0, z: tc, }
    3 => { x: sc, y: -1.0, z: -tc, }
    4 => { x: sc, y: -tc, z: 1.0, }
    5 => { x: -sc, y: -tc, z: -1.0, }
    _ => abort("cube_to_direction: face must be 0..5")
  }
  d.normalize()
}

// ---------------------------------------------------------------------------
// Equirectangular (lat/long) panoramas
// ---------------------------------------------------------------------------

///|
/// Direction -> equirectangular UV. The image centre (0.5, 0.5) looks down -Z (forward);
/// u increases toward +X (turning right), v increases DOWNWARD (v = 0 is straight up).
pub fn equirect_from_direction(dir : Vec3) -> Vec2 {
  let lon = @math.atan2(dir.x, -dir.z)
  let lat = @math.atan2(dir.y, @math.hypot(dir.x, dir.z))
  { x: 0.5 + lon / TAU, y: 0.5 - lat / PI, }
}

///|
pub fn equirect_to_direction(uv : Vec2) -> Vec3 {
  let lon = (uv.x - 0.5) * TAU
  let lat = (0.5 - uv.y) * PI
  let c = @math.cos(lat)
  { x: @math.sin(lon) * c, y: @math.sin(lat), z: -@math.cos(lon) * c, }
}

// ---------------------------------------------------------------------------
// Octahedral mapping (Cigolle et al. 2014): unit vector <-> [-1, 1]^2, uniform and cheap.
// Used for normal compression, probes, environment atlases.
// ---------------------------------------------------------------------------

///|
pub fn oct_encode(n : Vec3) -> Vec2 {
  let l1 = n.x.abs() + n.y.abs() + n.z.abs()
  let px = n.x / l1
  let py = n.y / l1
  if n.z < 0.0 {
    { x: (1.0 - py.abs()) * sign_nz(px), y: (1.0 - px.abs()) * sign_nz(py), }
  } else {
    { x: px, y: py, }
  }
}

///|
pub fn oct_decode(e : Vec2) -> Vec3 {
  let z = 1.0 - e.x.abs() - e.y.abs()
  let v : Vec3 = if z < 0.0 {
    {
      x: (1.0 - e.y.abs()) * sign_nz(e.x),
      y: (1.0 - e.x.abs()) * sign_nz(e.y),
      z,
    }
  } else {
    { x: e.x, y: e.y, z, }
  }
  v.normalize()
}

///|
/// [-1,1]^2 octahedral coordinate -> [0,1]^2 texture UV.
pub fn oct_to_uv(e : Vec2) -> Vec2 {
  { x: e.x * 0.5 + 0.5, y: e.y * 0.5 + 0.5, }
}

///|
pub fn uv_to_oct(uv : Vec2) -> Vec2 {
  { x: uv.x * 2.0 - 1.0, y: uv.y * 2.0 - 1.0, }
}
