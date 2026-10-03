# ky678/renky/system/graph — coordinate systems for a renderer (MoonBit, zero dependencies)

Only the *coordinate-system* layer of a renderer: every space, every transform between them, and the
conventions that glue them together. No third-party packages (only `moonbitlang/core`).

## Conventions (one place, no surprises)

| Item | Choice |
|---|---|
| Handedness / axes | **Right-handed, Y-up**; camera looks down **-Z**, +X right, +Y up |
| Scalars / angles | `Double`, radians (`Float` only for GPU upload helpers) |
| Matrices | column vectors `v' = M v`, fields `mRC`, `Mat4::to_col_major` is GL/WebGPU/Vulkan layout |
| Quaternion | Hamilton, `a * b` applies `b` first |
| Euler | all 12 orders, intrinsic `R = R_A(a) R_B(b) R_C(c)`, canonical branch on extraction |
| Clip / NDC depth | `DepthMode`: `NegOneToOne` (OpenGL), `ZeroToOne` (D3D/Vulkan/Metal/WebGPU), `ReversedZ`; finite **or infinite** far plane |
| Viewport / pixels | origin top-left, y down, pixel centres at +0.5 |
| Cube map | faces +X,-X,+Y,-Y,+Z,-Z (OpenGL table), v down; equirect `u = 0.5 + atan2(x,-z)/2π` |
| Spherical | azimuth 0 = -Z, +π/2 = -X; elevation from the horizon |
| Geodesy | WGS-84; ECEF → local ENU/NED → render frame (east=+X, up=+Y, north=-Z) |

## Spaces and what converts between them

```
 model ──Rigid/Mat4 TRS──▶ world ──view (Camera)──▶ view ──Projection──▶ clip ──÷w──▶ NDC ──Viewport──▶ screen
   ▲                          ▲                                                                           │
   └─ tangent / TBN           └─ convention (Z-up, Y-down, left-handed … ⇄ native) ◀── pick_ray / unproject ┘
 geodetic (lat,lon,h) ⇄ ECEF ⇄ ENU/NED ⇄ render world      Web-Mercator ⇄ tile / pixel / quadkey
 direction ⇄ cube face (u,v) ⇄ equirect ⇄ octahedral        texel ⇄ uv, wrap modes, bilinear footprint
 cartesian ⇄ spherical / cylindrical / polar                camera-relative (floating-origin) view matrix
```

Files: `vec*.mbt quat.mbt mat3.mbt mat4.mbt euler.mbt transform.mbt` (algebra) ·
`projection.mbt viewport.mbt camera.mbt geometry.mbt` (pipeline, rays, planes, frusta, AABB) ·
`geodetic.mbt mercator.mbt curvilinear.mbt` (earth & curvilinear systems) ·
`texture.mbt tangent.mbt` (UV, cube, octahedral, TBN, barycentrics) · `convention.mbt` (other engines' frames).

## Run

From the module root (`renky/`):

```
moon test -p ky678/renky/system/graph      # 75 tests; also --target wasm | js | native
python3 system/graph/tools/mutation_check.py   # injects 71 bugs; every one must be caught
moon run system/graph/cmd/demo --target native > demo.ppm      # software-rendered scene (docs/demo.png)
moon run system/graph/cmd/benchmark --target native --release  # micro-benchmarks
```

Use it: `moon add ky678/renky`, then `import { "ky678/renky/system/graph" }` and refer to it as `@graph`.

## How the tests avoid being decorative

* **External oracles** (`tools/gen_golden.py` → `golden_data_test.mbt`): numpy textbook pipelines
  (gluPerspective, XMMatrixPerspectiveFovRH, glOrtho, reversed-Z, D3D viewport), scipy
  (`Rotation`, Euler, rotvec), **PROJ/pyproj** (EPSG:4979→4978, EPSG:4326→3857, `topocentric` ENU),
  PROJ `Geod` for metres/degree, OSM tile formulas, Bing quadkey example.
* **Independent re-derivations**: a reference triangle rasteriser is compared pixel-by-pixel with
  analytic ray casting through `pick_ray`; frustum culling vs clip-space inequalities; closed-form inverses vs general inverse.
* **Mutation testing**: `tools/mutation_check.py` (71 mutants — sign flips, swapped operands, dropped
  terms) — all killed. It found and led to fixing three real gaps (no tests for AABB transform, ray-plane, frustum near/far).
