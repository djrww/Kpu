// Web Mercator (EPSG:3857) and the slippy-map tile pyramid.
//
// Angles in radians (longitude, latitude), like every other API here. Tile / world-pixel
// coordinates use the web-map layout: origin at the NORTH-WEST corner, +x east, +y SOUTH.

///|
/// Radius of the sphere used by EPSG:3857 (WGS84 semi-major axis), metres.
pub const MERCATOR_RADIUS : Double = 6378137.0

///|
/// Latitude beyond which Web Mercator is clipped: atan(sinh(PI)) = 85.0511287798 deg.
pub const MERCATOR_MAX_LAT : Double = 1.4844222297453324

///|
#valtype
pub(all) struct LonLat {
  lon : Double
  lat : Double
} derive(Eq, Debug)

///|
pub fn LonLat::new(lon : Double, lat : Double) -> LonLat {
  { lon, lat, }
}

///|
pub fn LonLat::from_degrees(lon_deg : Double, lat_deg : Double) -> LonLat {
  { lon: deg_to_rad(lon_deg), lat: deg_to_rad(lat_deg), }
}

///|
/// lon/lat -> EPSG:3857 metres (x east, y NORTH). Latitude is clamped to +-MERCATOR_MAX_LAT.
pub fn web_mercator_forward(p : LonLat) -> Vec2 {
  let lat = clamp(p.lat, -MERCATOR_MAX_LAT, MERCATOR_MAX_LAT)
  {
    x: MERCATOR_RADIUS * p.lon,
    y: MERCATOR_RADIUS * @math.atanh(@math.sin(lat)),
  }
}

///|
/// EPSG:3857 metres -> lon/lat.
pub fn web_mercator_inverse(m : Vec2) -> LonLat {
  {
    lon: m.x / MERCATOR_RADIUS,
    lat: @math.atan(@math.sinh(m.y / MERCATOR_RADIUS)),
  }
}

///|
/// Side of the whole-world square in pixels at (possibly fractional) `zoom`.
pub fn world_size_px(zoom : Double, tile_size : Double) -> Double {
  tile_size * @math.pow(2.0, zoom)
}

///|
/// lon/lat -> world pixel coordinates (origin NW corner, y down).
pub fn lonlat_to_world_pixel(
  p : LonLat,
  zoom : Double,
  tile_size : Double,
) -> Vec2 {
  let s = world_size_px(zoom, tile_size)
  let lat = clamp(p.lat, -MERCATOR_MAX_LAT, MERCATOR_MAX_LAT)
  {
    x: (p.lon + PI) / TAU * s,
    y: (0.5 - @math.atanh(@math.sin(lat)) / TAU) * s,
  }
}

///|
pub fn world_pixel_to_lonlat(
  px : Vec2,
  zoom : Double,
  tile_size : Double,
) -> LonLat {
  let s = world_size_px(zoom, tile_size)
  {
    lon: px.x / s * TAU - PI,
    lat: @math.atan(@math.sinh((0.5 - px.y / s) * TAU)),
  }
}

///|
/// EPSG:3857 metres -> world pixels (and back); handy for mixing meter-based and tile-based math.
pub fn mercator_to_world_pixel(
  m : Vec2,
  zoom : Double,
  tile_size : Double,
) -> Vec2 {
  let s = world_size_px(zoom, tile_size)
  let half = PI * MERCATOR_RADIUS
  { x: (m.x + half) / (2.0 * half) * s, y: (half - m.y) / (2.0 * half) * s, }
}

///|
pub fn world_pixel_to_mercator(
  px : Vec2,
  zoom : Double,
  tile_size : Double,
) -> Vec2 {
  let s = world_size_px(zoom, tile_size)
  let half = PI * MERCATOR_RADIUS
  { x: px.x / s * 2.0 * half - half, y: half - px.y / s * 2.0 * half, }
}

///|
/// Ground resolution: metres per pixel at a latitude (true ground distance, not Mercator units).
pub fn ground_resolution(
  lat : Double,
  zoom : Double,
  tile_size : Double,
) -> Double {
  @math.cos(lat) * TAU * MERCATOR_RADIUS / world_size_px(zoom, tile_size)
}

///|
/// Integer tile (x, y) containing a position at integer zoom (clamped to the pyramid).
pub fn lonlat_to_tile(p : LonLat, zoom : Int, tile_size : Double) -> (Int, Int) {
  let px = lonlat_to_world_pixel(p, Double::from_int(zoom), tile_size)
  let n = 1 << zoom
  let tx = @math.floor(px.x / tile_size).to_int()
  let ty = @math.floor(px.y / tile_size).to_int()
  (tx.clamp(min=0, max=n - 1), ty.clamp(min=0, max=n - 1))
}

///|
/// North-west corner of tile (x, y) at `zoom` (x, y may be fractional for sub-tile positions).
pub fn tile_to_lonlat(x : Double, y : Double, zoom : Int) -> LonLat {
  let n = Double::from_int(1 << zoom)
  {
    lon: x / n * TAU - PI,
    lat: @math.atan(@math.sinh(PI * (1.0 - 2.0 * y / n))),
  }
}

///|
/// (west, south, east, north) lon/lat bounds of a tile.
pub fn tile_bounds(x : Int, y : Int, zoom : Int) -> (LonLat, LonLat) {
  let nw = tile_to_lonlat(Double::from_int(x), Double::from_int(y), zoom)
  let se = tile_to_lonlat(
    Double::from_int(x + 1),
    Double::from_int(y + 1),
    zoom,
  )
  (nw, se)
}

///|
/// Bing-style quadkey of a tile ("" at zoom 0).
pub fn tile_to_quadkey(x : Int, y : Int, zoom : Int) -> String {
  let sb = StringBuilder()
  let mut i = zoom
  while i > 0 {
    let mask = 1 << (i - 1)
    let mut digit = 0
    if (x & mask) != 0 {
      digit += 1
    }
    if (y & mask) != 0 {
      digit += 2
    }
    sb.write_string(digit.to_string())
    i -= 1
  }
  sb.to_string()
}

///|
/// Parse a quadkey back into (x, y, zoom). `None` for any character other than 0-3.
pub fn quadkey_to_tile(key : String) -> (Int, Int, Int)? {
  let mut x = 0
  let mut y = 0
  let mut zoom = 0
  for c in key {
    x = x << 1
    y = y << 1
    zoom += 1
    match c {
      '0' => ()
      '1' => x = x | 1
      '2' => y = y | 1
      '3' => {
        x = x | 1
        y = y | 1
      }
      _ => return None
    }
  }
  Some((x, y, zoom))
}
