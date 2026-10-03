#!/usr/bin/env python3
"""Mutation testing: inject a deliberate bug into the library, run the suite, and require it to FAIL.
A surviving mutant would mean the tests cannot see that bug.  Run: python3 tools/mutation_check.py"""
import shutil, subprocess, os, sys, tempfile, re
PKG_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # .../system/graph
MODROOT = os.path.dirname(os.path.dirname(PKG_DIR))                       # module root (moon.mod)
PKG = "ky678/renky/system/graph"
MUTANTS = [
 ("quat.mbt", "x: v.x + self.w * tx + (self.y * tz - self.z * ty),", "x: v.x + self.w * tx - (self.y * tz - self.z * ty),", "Quat::rotate sign error"),
 ("quat.mbt", "m.m21 - m.m12) / s,\n      y: (m.m02 - m.m20) / s,", "m.m12 - m.m21) / s,\n      y: (m.m02 - m.m20) / s,", "Quat::from_mat3 sign"),
 ("quat.mbt", "let wa = @math.sin((1.0 - t) * th) / s", "let wa = @math.sin(t * th) / s", "slerp weight swap"),
 ("mat3.mbt", "m01: 2.0 * (xy - wz),", "m01: 2.0 * (xy + wz),", "Mat3::from_quat sign"),
 ("mat4.mbt", "m01: a.m00 * b.m01 + a.m01 * b.m11 + a.m02 * b.m21 + a.m03 * b.m31,", "m01: a.m00 * b.m01 + a.m01 * b.m11 + a.m02 * b.m21 + a.m03 * b.m30,", "Mat4 mul typo"),
 ("mat4.mbt", "m13: (a.m20 * s5 - a.m22 * s2 + a.m23 * s1) * i,", "m13: (a.m20 * s5 - a.m22 * s2 - a.m23 * s1) * i,", "Mat4 inverse cofactor sign"),
 ("mat4.mbt", "m03: -(self.m00 * tx + self.m10 * ty + self.m20 * tz),", "m03: -(self.m00 * tx + self.m01 * ty + self.m02 * tz),", "rigid_inverse translation"),
 ("mat4.mbt", "let sx = if det < 0.0 { -l0 } else { l0 }", "let sx = l0", "decompose ignores mirroring"),
 ("mat3.mbt", "Some(i) => Some(i.transpose())", "Some(i) => Some(i)", "normal matrix missing transpose"),
 ("mat3.mbt", "let y = z.cross(x)\n  Mat3::from_cols(x, y, z)", "let y = x.cross(z)\n  Mat3::from_cols(x, y, z)", "look_rotation handedness"),
 ("euler.mbt", "let odd = if (a0 + 1) % 3 == a1 { 0 } else { 1 }", "let odd = if (a0 + 2) % 3 == a1 { 0 } else { 1 }", "euler parity"),
 ("euler.mbt", "(i, j, k) = order.axes()\n  (axis_quat(i, a) * axis_quat(j, b) * axis_quat(k, c)).normalize()", "(i, j, k) = order.axes()\n  (axis_quat(k, c) * axis_quat(j, b) * axis_quat(i, a)).normalize()", "euler composition order"),
 ("projection.mbt", "ZeroToOne => (-f / d, -f * n / d)", "ZeroToOne => (-f / d, f * n / d)", "Z01 B sign"),
 ("projection.mbt", "ReversedZ => (n / d, n * f / d)", "ReversedZ => (-n / d, n * f / d)", "reversed A sign"),
 ("projection.mbt", "NegOneToOne => (-(f + n) / d, -2.0 * f * n / d)", "NegOneToOne => (-(f - n) / d, -2.0 * f * n / d)", "GL A"),
 ("projection.mbt", "m02: (self.right + self.left) / rl,\n        m03: 0.0,", "m02: (self.right - self.left) / rl,\n        m03: 0.0,", "off-centre frustum m02"),
 ("projection.mbt", "x: -zv * (ndc.x + m02) / m00,", "x: zv * (ndc.x + m02) / m00,", "unproject x sign"),
 ("projection.mbt", "Perspective => -b / (ndc_z + a)", "Perspective => b / (ndc_z + a)", "view_z_from_ndc"),
 ("viewport.mbt", "y: self.y + (0.5 - ndc.y * 0.5) * self.height,\n    z: self.min_depth", "y: self.y + (0.5 + ndc.y * 0.5) * self.height,\n    z: self.min_depth", "viewport y not flipped"),
 ("projection.mbt", "NegOneToOne => ndc_z * 0.5 + 0.5", "NegOneToOne => ndc_z", "GL depth to unit"),
 ("camera.mbt", "let inv_rot = Mat4::from_mat3(Mat3::from_quat(self.pose.rotation).transpose())", "let inv_rot = Mat4::from_mat3(Mat3::from_quat(self.pose.rotation))", "RTE view uses R not R^T"),
 ("camera.mbt", "let eye = target - f.scale(distance)", "let eye = target + f.scale(distance)", "orbit direction"),
 ("geometry.mbt", "ReversedZ => (r3 - r2, r2)", "ReversedZ => (r2, r3 - r2)", "frustum reversed-Z near/far"),
 ("geometry.mbt", "left: plane_from_vec4(r3 + r0),", "left: plane_from_vec4(r3 - r0),", "frustum left plane"),
 ("geometry.mbt", "m.m00.abs() * h.x + m.m01.abs() * h.y + m.m02.abs() * h.z", "m.m00 * h.x + m.m01 * h.y + m.m02 * h.z", "Aabb::transform without abs"),
 ("transform.mbt", "translation: a.rotation.rotate(b.translation) + a.translation,", "translation: b.rotation.rotate(a.translation) + b.translation,", "Rigid compose order"),
 ("geodetic.mbt", "z: (n * (1.0 - self.e2()) + g.height) * sl,", "z: (n + g.height) * sl,", "ECEF z uses N not N(1-e2)"),
 ("geodetic.mbt", "let nl = @math.atan2(p.z + e2 * n * s, rho)", "let nl = @math.atan2(p.z, rho)", "from_ecef iteration"),
 ("geodetic.mbt", "axes: Mat3::from_cols(g.east_ecef(), g.north_ecef(), g.up_ecef()),", "axes: Mat3::from_cols(g.north_ecef(), g.east_ecef(), g.up_ecef()),", "ENU axis order"),
 ("geodetic.mbt", "Vec3::new(0.0, -1.0, 0.0),\n  )\n  let r = q", "Vec3::new(0.0, 1.0, 0.0),\n  )\n  let r = q", "render mapping north sign"),
 ("mercator.mbt", "y: MERCATOR_RADIUS * @math.atanh(@math.sin(lat))", "y: MERCATOR_RADIUS * @math.sin(lat)", "mercator y"),
 ("mercator.mbt", "y: (0.5 - @math.atanh(@math.sin(lat)) / TAU) * s,", "y: (0.5 + @math.atanh(@math.sin(lat)) / TAU) * s,", "tile y direction"),
 ("mercator.mbt", "'1' => x = x | 1", "'1' => y = y | 1", "quadkey digit"),
 ("texture.mbt", "      sc = -dir.z\n      tc = -dir.y\n    } else {\n      face = CUBE_NEG_X", "      sc = dir.z\n      tc = -dir.y\n    } else {\n      face = CUBE_NEG_X", "cube +X sc"),
 ("texture.mbt", "3 => { x: sc, y: -1.0, z: -tc }", "3 => { x: sc, y: -1.0, z: tc }", "cube_to_direction -Y"),
 ("texture.mbt", "let lon = @math.atan2(dir.x, -dir.z)", "let lon = @math.atan2(-dir.x, -dir.z)", "equirect handedness"),
 ("texture.mbt", "2 * n - 1 - m", "2 * n - m", "mirrored repeat off by one"),
 ("texture.mbt", "let x = uv.x * Double::from_int(width) - 0.5", "let x = uv.x * Double::from_int(width)", "bilinear half-texel"),
 ("texture.mbt", "{ x: (1.0 - e.y.abs()) * sign_nz(e.x), y: (1.0 - e.x.abs()) * sign_nz(e.y), z }", "{ x: (1.0 - e.y.abs()) * sign_nz(e.y), y: (1.0 - e.x.abs()) * sign_nz(e.x), z }", "oct decode fold"),
 ("tangent.mbt", "Some({ tangent: t, bitangent: nxt.scale(h), normal: nn })", "Some({ tangent: t, bitangent: nxt, normal: nn })", "tangent handedness dropped"),
 ("tangent.mbt", "let a = screen_bary.x / w0", "let a = screen_bary.x * w0", "perspective-correct inverted"),
 ("tangent.mbt", "let wa = edge_function(b, c, p) * inv", "let wa = edge_function(a, c, p) * inv", "barycentric edge"),
 ("convention.mbt", "let a = m.mul_vec3(Vec3::new(q.x, q.y, q.z)).scale(det)", "let a = m.mul_vec3(Vec3::new(q.x, q.y, q.z))", "quat conversion ignores mirror"),
 ("convention.mbt", "to.basis() * self.basis().transpose()", "self.basis().transpose() * to.basis()", "basis_change order"),
 ("curvilinear.mbt", "azimuth: if h == 0.0 { 0.0 } else { @math.atan2(-p.x, -p.z) },", "azimuth: if h == 0.0 { 0.0 } else { @math.atan2(p.x, -p.z) },", "spherical azimuth sign"),
 ("geometry.mbt", "if t >= 0.0 {\n    Some(t)", "if t > 5.0 {\n    Some(t)", "ray-plane t threshold"),
 ("scalar.mbt", "let r = a - TAU * @math.floor(a / TAU + 0.5)", "let r = a - TAU * @math.floor(a / TAU)", "wrap_pi"),
 ("geometry.mbt", "let t = -(pl.normal.dot(self.origin) + pl.d) / denom", "let t = (pl.normal.dot(self.origin) + pl.d) / denom", "ray-plane distance sign"),
 ("geometry.mbt", "x: if (i & 1) != 0 { self.max.x } else { self.min.x },", "x: if (i & 2) != 0 { self.max.x } else { self.min.x },", "Aabb::corner bit mapping"),
 ("geometry.mbt", "Plane::from_point_normal(a, (b - a).cross(c - a))", "Plane::from_point_normal(a, (c - a).cross(b - a))", "Plane::from_points winding"),
 ("geometry.mbt", "if a > t0 {\n        t0 = a\n      }", "if a < t0 {\n        t0 = a\n      }", "slab t0 update"),
 ("geometry.mbt", "if pl.signed_distance(pv) < 0.0 {\n    return 0", "if pl.signed_distance(pv) < -1000.0 {\n    return 0", "classify_aabb never Outside"),
 ("transform.mbt", "{ rotation: ri, translation: -ri.rotate(self.translation) }", "{ rotation: ri, translation: ri.rotate(self.translation) }", "Rigid inverse translation"),
 ("camera.mbt", "{ origin: self.pose.mul_point(o), dir: self.pose.mul_dir(d) }", "{ origin: self.pose.mul_point(o), dir: d }", "pick_ray not rotated to world"),
 ("camera.mbt", "ray.at(distance / cosang)", "ray.at(distance)", "screen_to_world_at_distance perspective"),
 ("projection.mbt", "ZeroToOne => (-1.0, -n)", "ZeroToOne => (-1.0, n)", "infinite-far Z01"),
 ("projection.mbt", "ZeroToOne => (-1.0 / d, -n / d)", "ZeroToOne => (-1.0 / d, n / d)", "ortho Z01 B"),
 ("projection.mbt", "ReversedZ => (1.0 / d, f / d)", "ReversedZ => (1.0 / d, n / d)", "ortho reversed B"),
 ("projection.mbt", "m11: 1.0 / m11,\n        m12: 0.0,\n        m13: m12 / m11,", "m11: 1.0 / m11,\n        m12: 0.0,\n        m13: -m12 / m11,", "inverse_mat4 y shift"),
 ("mat3.mbt", "m01: (self.m02 * self.m21 - self.m01 * self.m22) * i,", "m01: (self.m01 * self.m22 - self.m02 * self.m21) * i,", "Mat3 inverse cofactor"),
 ("quat.mbt", "sv.scale(2.0 * @math.atan2(s, q.w) / s)", "sv.scale(@math.atan2(s, q.w) / s)", "rotation vector scale"),
 ("quat.mbt", "{ x: p.x, y: p.y, z: p.z, w: 0.0 }", "{ x: p.x, y: p.y, z: p.z, w: 1.0 }", "from_to antiparallel"),
 ("geodetic.mbt", "self.axes.transpose().mul_vec3(ecef - self.origin)", "self.axes.mul_vec3(ecef - self.origin)", "to_local not transposed"),
 ("geodetic.mbt", "axes: Mat3::from_cols(g.north_ecef(), g.east_ecef(), -g.up_ecef()),", "axes: Mat3::from_cols(g.north_ecef(), g.east_ecef(), g.up_ecef()),", "NED down sign"),
 ("geodetic.mbt", "@math.atan((1.0 - self.e2()) * @math.tan(lat))", "@math.atan((1.0 - self.f) * @math.tan(lat))", "geocentric latitude uses f"),
 ("mercator.mbt", "@math.cos(lat) * TAU * MERCATOR_RADIUS / world_size_px(zoom, tile_size)", "TAU * MERCATOR_RADIUS / world_size_px(zoom, tile_size)", "ground resolution cos"),
 ("texture.mbt", "{ x: uv.x * Double::from_int(width), y: uv.y * Double::from_int(height) }", "{ x: uv.x * Double::from_int(height), y: uv.y * Double::from_int(width) }", "uv_to_texel axes"),
 ("texture.mbt", "let t = u - 2.0 * @math.floor(u * 0.5)", "let t = u - @math.floor(u)", "mirrored coord"),
 ("curvilinear.mbt", "z: -self.radius * @math.cos(self.azimuth) * ce,", "z: self.radius * @math.cos(self.azimuth) * ce,", "spherical z sign"),
 ("viewport.mbt", "x: (p.x - self.x) / self.width * 2.0 - 1.0,\n    y: 1.0 - (p.y - self.y) / self.height * 2.0,\n  }\n}\n\n///|\npub fn Viewport::ndc_to_screen_xy", "x: (p.x - self.x) / self.width * 2.0 - 1.0,\n    y: (p.y - self.y) / self.height * 2.0 - 1.0,\n  }\n}\n\n///|\npub fn Viewport::ndc_to_screen_xy", "screen_to_ndc_xy y"),
 ("euler.mbt", "r1 = PI - r1", "r1 = -r1", "euler canonicalisation"),
]

def make_regex(pat):
    """Whitespace-insensitive pattern that also tolerates the trailing commas `moon fmt` inserts."""
    out = []
    for ch in pat:
        if ch.isspace():
            continue
        if ch in "})":
            out.append(r"(?:\s*,)?")
        out.append(r"\s*" + re.escape(ch))
    return re.compile("".join(out)[3:])  # drop the leading \s*
def run(cwd):
    p = subprocess.run(["moon", "test", "-p", PKG], cwd=cwd, capture_output=True, text=True, timeout=600)
    return p.returncode, p.stdout + p.stderr
def main():
    env_path = os.path.expanduser("~/.moon/bin")
    os.environ["PATH"] = env_path + os.pathsep + os.environ["PATH"]
    tmp = tempfile.mkdtemp(prefix="mut_")
    work = os.path.join(tmp, "w")
    shutil.copytree(MODROOT, work, ignore=shutil.ignore_patterns("_build", ".git", "target", ".mooncakes"))
    rc, out = run(work)
    if rc != 0:
        print("baseline must pass first:\n", out[-1500:]); sys.exit(2)
    killed = 0; survived = []
    for (f, old, new, desc) in MUTANTS:
        path = os.path.join(work, "system", "graph", f)
        src = open(path).read()
        rx = make_regex(old)
        found = list(rx.finditer(src))
        if len(found) != 1:
            print(f"!! pattern for '{desc}' matches {len(found)} times in {f}"); survived.append(desc + " (BAD PATTERN)"); continue
        mm = found[0]
        open(path, "w").write(src[:mm.start()] + new + src[mm.end():])
        rc, out = run(work)
        open(path, "w").write(src)
        if rc != 0:
            killed += 1
            # name of first failing test
            first = [l for l in out.splitlines() if "failed" in l and "test" in l][:1]
            print(f"KILLED   {desc:45s} <- {first[0][:110] if first else 'compile error?'}")
        else:
            survived.append(desc); print(f"SURVIVED {desc}")
    print(f"\n{killed}/{len(MUTANTS)} mutants killed")
    if survived:
        print("SURVIVORS:", survived); sys.exit(1)
if __name__ == "__main__":
    main()
