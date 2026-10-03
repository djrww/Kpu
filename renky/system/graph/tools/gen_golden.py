#!/usr/bin/env python3
"""Generates golden_*_test.mbt: reference values computed by INDEPENDENT software
(numpy / scipy / pyproj / mpmath), never by the code under test.
Run:  python3 tools/gen_golden.py   (needs numpy, scipy, pyproj, mpmath)"""
import numpy as np
from scipy.spatial.transform import Rotation as R, Slerp

rng = np.random.default_rng(20261002)
OUT = []

def fl(x):
    return repr(float(x))

def emit(name, arr):
    arr = np.asarray(arr, dtype=float).ravel()
    OUT.append(f"///|\nlet gold_{name[2:].lower()} : FixedArray[Double] = [\n  " +
               ",\n  ".join(", ".join(fl(v) for v in arr[i:i+4]) for i in range(0, len(arr), 4)) + ",\n]\n")

def rowmajor(m): return np.asarray(m).ravel()

# ---- 1. generic 4x4: A, B, A@B, inv(A), det(A), v, A@v  (stride 16*4+1+4+4)
rec = []
N_MAT = 24
for _ in range(N_MAT):
    A = rng.uniform(-2, 2, (4, 4)) + np.eye(4) * rng.uniform(-3, 3)
    B = rng.uniform(-2, 2, (4, 4))
    v = rng.uniform(-3, 3, 4)
    rec += [*rowmajor(A), *rowmajor(B), *rowmajor(A @ B), *rowmajor(np.linalg.inv(A)),
            np.linalg.det(A), *v, *(A @ v)]
emit("G_MAT4", rec)
# stride:
MAT_STRIDE = 16*4 + 1 + 4 + 4

# ---- 2. TRS: t(3) q(4 xyzw) s(3) | M(16 rowmajor) | inv(16) | normal matrix(9) | probe p(3) M@p(3)
rec = []
N_TRS = 24
for _ in range(N_TRS):
    t = rng.uniform(-50, 50, 3)
    q = R.random(random_state=int(rng.integers(1 << 30)))
    s = rng.uniform(0.2, 4, 3) * rng.choice([1, 1, 1, -1], 3)
    # keep at most reflections allowed; generic
    M = np.eye(4); M[:3, :3] = q.as_matrix() @ np.diag(s); M[:3, 3] = t
    p = rng.uniform(-5, 5, 3)
    N = np.linalg.inv(M[:3, :3]).T
    rec += [*t, *q.as_quat(), *s, *rowmajor(M), *rowmajor(np.linalg.inv(M)), *rowmajor(N), *p, *(M[:3, :3] @ p + t)]
emit("G_TRS", rec)

# ---- 3. quats: q1(4) q2(4) v(3) | (q1*q2)(4) | q1*v (3) | slerp t | slerp(4) | rotvec q1 (3) | mat q1 (9 rowmajor)
rec = []
N_Q = 32
for _ in range(N_Q):
    r1 = R.random(random_state=int(rng.integers(1 << 30)))
    r2 = R.random(random_state=int(rng.integers(1 << 30)))
    v = rng.uniform(-4, 4, 3)
    t = rng.uniform(0, 1)
    sl = Slerp([0, 1], R.concatenate([r1, r2]))([t])[0]
    rec += [*r1.as_quat(), *r2.as_quat(), *v, *(r1 * r2).as_quat(), *r1.apply(v), t, *sl.as_quat(),
            *r1.as_rotvec(), *rowmajor(r1.as_matrix())]
emit("G_QUAT", rec)

# ---- 4. euler: for each of 12 seqs, 8 cases: angles(3) | matrix(9)  [intrinsic == scipy upper-case]
SEQS = ["XYZ","XZY","YXZ","YZX","ZXY","ZYX","XYX","XZX","YXY","YZY","ZXZ","ZYZ"]
rec = []
for seq in SEQS:
    for _ in range(8):
        proper = seq[0] == seq[2]
        a = rng.uniform(-3.0, 3.0)
        c = rng.uniform(-3.0, 3.0)
        b = rng.uniform(0.15, 2.9) if proper else rng.uniform(-1.4, 1.4)
        r = R.from_euler(seq, [a, b, c])
        rec += [a, b, c, *rowmajor(r.as_matrix())]
        # scipy's own decomposition (angles) for cross-check of to_euler on generic input
        rec += [*r.as_euler(seq)]
emit("G_EULER", rec)
# also extrinsic cross-check: lower-case seq == reversed intrinsic
rec = []
for seq in SEQS:
    for _ in range(4):
        proper = seq[0] == seq[2]
        a = rng.uniform(-3.0, 3.0); c = rng.uniform(-3.0, 3.0)
        b = rng.uniform(0.15, 2.9) if proper else rng.uniform(-1.4, 1.4)
        r = R.from_euler(seq.lower(), [a, b, c])  # extrinsic: first about fixed seq[0] by a ...
        rec += [a, b, c, *rowmajor(r.as_matrix())]
emit("G_EULER_EXT", rec)

# ---- 5. look_at (gluLookAt formula, numpy): eye(3) target(3) up(3) | view(16 rowmajor)
def glu_look_at(eye, tgt, up):
    f = tgt - eye; f /= np.linalg.norm(f)
    s = np.cross(f, up); s /= np.linalg.norm(s)
    u = np.cross(s, f)
    M = np.eye(4); M[0, :3] = s; M[1, :3] = u; M[2, :3] = -f
    T = np.eye(4); T[:3, 3] = -eye
    return M @ T
rec = []
for _ in range(16):
    eye = rng.uniform(-20, 20, 3); tgt = rng.uniform(-20, 20, 3); up = rng.uniform(-1, 1, 3)
    rec += [*eye, *tgt, *up, *rowmajor(glu_look_at(eye, tgt, up))]
emit("G_LOOKAT", rec)


# =====================================================================================
# Independent pipeline oracle (numpy; textbook gluPerspective / XMMatrixPerspectiveFovRH /
# swapped-plane reversed-Z formulas; D3D viewport transform).
# case: eye(3) target(3) up(3) fov aspect near far(0=inf) W H mode | K * [p(3) -> sx sy depth ndc(3) viewz]
# =====================================================================================
def persp_matrix(fov, aspect, n, f, mode):
    t = 1.0 / np.tan(fov / 2)
    P = np.zeros((4, 4)); P[0, 0] = t / aspect; P[1, 1] = t; P[3, 2] = -1
    if mode == 1:   # OpenGL, gluPerspective
        if f is None: P[2, 2] = -1; P[2, 3] = -2 * n
        else:
            P[2, 2] = (f + n) / (n - f); P[2, 3] = 2 * f * n / (n - f)
    elif mode == 0:  # D3D / Vulkan / Metal RH zero-to-one (XMMatrixPerspectiveFovRH)
        if f is None: P[2, 2] = -1; P[2, 3] = -n
        else:
            P[2, 2] = f / (n - f); P[2, 3] = f * n / (n - f)
    else:            # reversed-Z: D3D formula with near and far swapped
        if f is None: P[2, 2] = 0; P[2, 3] = n
        else:
            P[2, 2] = n / (f - n); P[2, 3] = n * f / (f - n)
    return P

K = 6
rec = []
rng2 = np.random.default_rng(777)
for case in range(36):
    mode = case % 3
    inf = (case // 3) % 2 == 1
    eye = rng2.uniform(-30, 30, 3); tgt = rng2.uniform(-5, 5, 3); up = np.array([0, 1, 0.0]) + rng2.uniform(-.3, .3, 3)
    fov = rng2.uniform(0.4, 1.9); aspect = rng2.uniform(0.6, 2.2); n = rng2.uniform(0.05, 2.0); f = None if inf else rng2.uniform(50, 5000)
    W = float(rng2.integers(200, 2000)); H = float(rng2.integers(200, 1500))
    V = glu_look_at(eye, tgt, up)
    P = persp_matrix(fov, aspect, n, f, mode)
    rec += [*eye, *tgt, *up, fov, aspect, n, 0.0 if f is None else f, W, H, mode]
    got = 0
    while got < K:
        # random point in front of camera: pick via view space
        pv = np.array([rng2.uniform(-1, 1) * 5, rng2.uniform(-1, 1) * 5, -rng2.uniform(n * 1.5 + 0.1, 40)])
        pw = (np.linalg.inv(V) @ np.append(pv, 1))[:3]
        c = P @ V @ np.append(pw, 1)
        ndc = c[:3] / c[3]
        if abs(ndc[0]) > 1.5 or abs(ndc[1]) > 1.5: continue
        sx = (ndc[0] + 1) * W / 2; sy = (1 - ndc[1]) * H / 2
        d = ndc[2] if mode != 1 else ndc[2] * 0.5 + 0.5   # window depth in [0,1]
        rec += [*pw, sx, sy, d, *ndc, pv[2]]
        got += 1
emit("G_PIPE", rec)

# ---- orthographic: l r b t n f mode | pv(3) -> ndc(3) | formulas from glOrtho / XMMatrixOrthographicOffCenterRH
rec = []
for case in range(18):
    mode = case % 3
    l = rng2.uniform(-12, -1); r = rng2.uniform(1, 12); b = rng2.uniform(-9, -1); t = rng2.uniform(1, 9)
    n = rng2.uniform(-1, 2); f = n + rng2.uniform(5, 300)
    O = np.zeros((4, 4)); O[0, 0] = 2 / (r - l); O[1, 1] = 2 / (t - b); O[3, 3] = 1
    O[0, 3] = -(r + l) / (r - l); O[1, 3] = -(t + b) / (t - b)
    if mode == 1:   # glOrtho
        O[2, 2] = -2 / (f - n); O[2, 3] = -(f + n) / (f - n)
    elif mode == 0:  # XMMatrixOrthographicOffCenterRH
        O[2, 2] = 1 / (n - f); O[2, 3] = n / (n - f)
    else:            # reversed: swap n,f
        O[2, 2] = 1 / (f - n); O[2, 3] = f / (f - n)
    rec += [l, r, b, t, n, f, mode]
    for _ in range(5):
        pv = np.array([rng2.uniform(l, r), rng2.uniform(b, t), -rng2.uniform(n, f)])
        c = O @ np.append(pv, 1)
        rec += [*pv, *c[:3]]
emit("G_ORTHO", rec)

# =====================================================================================
# Geodesy via PROJ (pyproj)
# =====================================================================================
from pyproj import Transformer
geo2ecef = Transformer.from_crs("EPSG:4979", "EPSG:4978", always_xy=True)
ecef2geo = Transformer.from_crs("EPSG:4978", "EPSG:4979", always_xy=True)
ll2merc = Transformer.from_crs("EPSG:4326", "EPSG:3857", always_xy=True)
merc2ll = Transformer.from_crs("EPSG:3857", "EPSG:4326", always_xy=True)
# lat lon h | x y z
pts = [(0, 0, 0), (45, 0, 0), (89.9999, 10, 0), (-89.9999, -170, 50), (0, 90, 0), (0, 180, 100),
       (22.3193, 114.1694, 8.0), (27.9881, 86.9250, 8848.86), (-33.8688, 151.2093, 58),
       (51.4778, -0.0015, 46), (35.3606, 138.7274, 3776), (-90 + 1e-7, 0, 0), (10, 20, 35786000.0),
       (60, -120, -5000), (-12, 77, -8000)]
for _ in range(40):
    pts.append((rng2.uniform(-89.5, 89.5), rng2.uniform(-180, 180), rng2.uniform(-10000, 2.0e6)))
rec = []
for lat, lon, h in pts:
    x, y, z = geo2ecef.transform(lon, lat, h)
    rec += [lat, lon, h, x, y, z]
emit("G_GEO", rec)

# ENU via PROJ topocentric: origin lat lon h ; point lat lon h ; => e n u
rec = []
for _ in range(30):
    o = (rng2.uniform(-80, 80), rng2.uniform(-180, 180), rng2.uniform(0, 3000))
    p = (o[0] + rng2.uniform(-0.2, 0.2), o[1] + rng2.uniform(-0.2, 0.2), rng2.uniform(-100, 12000))
    tr = Transformer.from_pipeline(
        f"+proj=pipeline +step +proj=cart +ellps=WGS84 +step +proj=topocentric +ellps=WGS84 +lat_0={o[0]} +lon_0={o[1]} +h_0={o[2]}")
    e, n_, u = tr.transform(p[1], p[0], p[2])
    # NED order = (n, e, -u)
    rec += [*o, *p, e, n_, u]
emit("G_ENU", rec)

# mercator: lon lat | mx my   (degrees)
rec = []
mpts = [(0, 0), (180, 0), (-180, 0), (0, 85.0511287798066), (0, -85.0511287798066), (114.1694, 22.3193), (-0.1276, 51.5072)]
for _ in range(40): mpts.append((rng2.uniform(-180, 180), rng2.uniform(-85, 85)))
for lon, lat in mpts:
    x, y = ll2merc.transform(lon, lat)
    rec += [lon, lat, x, y]
emit("G_MERC", rec)

# tiles (OSM wiki formulas): lon lat zoom | xtile ytile | nw lon lat of that tile | pixel x y at tile 256
import math
def osm_tile(lat, lon, z):
    n = 2 ** z
    x = int(math.floor((lon + 180.0) / 360.0 * n))
    lr = math.radians(lat)
    y = int(math.floor((1.0 - math.asinh(math.tan(lr)) / math.pi) / 2.0 * n))
    return max(0, min(n - 1, x)), max(0, min(n - 1, y))
def tile_nw(x, y, z):
    n = 2 ** z
    lon = x / n * 360.0 - 180.0
    lat = math.degrees(math.atan(math.sinh(math.pi * (1 - 2 * y / n))))
    return lon, lat
rec = []
for _ in range(60):
    z = int(rng2.integers(0, 20))
    lon = rng2.uniform(-179.99, 179.99); lat = rng2.uniform(-85.0, 85.0)
    x, y = osm_tile(lat, lon, z)
    nl, nt = tile_nw(x, y, z)
    # fractional world pixel (tile 256)
    wx = (lon + 180) / 360 * 256 * 2 ** z
    wy = (1 - math.asinh(math.tan(math.radians(lat))) / math.pi) / 2 * 256 * 2 ** z
    rec += [lon, lat, z, x, y, nl, nt, wx, wy]
emit("G_TILE", rec)

if __name__ == "__main__":
    import sys
    header = "// GENERATED by tools/gen_golden.py from numpy/scipy/pyproj/mpmath -- DO NOT EDIT.\n\n"
    open("golden_data_test.mbt", "w").write(header + "\n".join(OUT))
    print("ok", len("".join(OUT)))
