# renky

MoonBit module **`ky678/renky`** — layout:

```
renky/
└─ system/
   ├─ graph/             ← coordinate systems for a renderer (right-handed, Y-up): see system/graph/README.md
   └─ algorithms/
      └─ omniring/        ← the axiom book (公理書): sources, tests, proofs, docs
         ├─ kernel/ select/            A01–A07   proof-enabled (moon prove, Why3 + Z3)
         ├─ arith/ linalg/ ff/ poly/   A08–A34   finite-field skeleton + linear algebra
         ├─ ideal/                     A35–A43   live ideal state (即時理想狀態)
         ├─ heat/ tropical/            A44–A51   Maslov dequantization → (min,+) semiring
         ├─ graded/ semigroup/         A52–A62   graded rings, affine semigroup rings
         ├─ padic/ smooth/             A63–A71   p-adic integers, ring of smooth germs
         ├─ homology/ cohom/           A72–A78   homology + cohomology rings
         ├─ fsplit/                    A79–A81   perfect F-split rings
         ├─ classgroup/                A82–A86   ideal class group correspondence ring
         ├─ loop/                      A87–A89   the closed computational chain (閉環)
         ├─ cmd/report/                          ledger printer (main package)
         └─ docs/{AXIOMS,THEOREMS,ALGORITHMS}.md
```

Package paths are `ky678/renky/system/algorithms/omniring/<pkg>`, e.g.

```moonbit
import { "ky678/renky/system/algorithms/omniring/loop" }

fn main {
  println(@loop.run_all(2).to_string())   // closed = true
}
```

## Verification ledger

| layer | command | result |
|---|---|---|
| machine-checked proofs | `moon prove` | 2 of 2 packages proved, **9 goals proved** |
| executable laws | `moon test` | **122 passed, 0 failed** |
| build | `moon check` | 0 errors |
| algorithms | `docs/ALGORITHMS.md` | **A01–A89** (89 > 18), no third-party dependencies |

See [`system/algorithms/omniring/README.md`](system/algorithms/omniring/README.md) for the full
axiom book (命案對帳表, 套件地圖, 閉環鏈條, 義務自證的兩層), and
[`system/algorithms/omniring/docs/`](system/algorithms/omniring/docs/) for
公理/定義/算式 (`AXIOMS.md`), 定理帳本 (`THEOREMS.md`) and 演算法登記冊 (`ALGORITHMS.md`).

## Reproduce

```bash
cd renky
moon check && moon test && moon prove
moon run system/algorithms/omniring/cmd/report
```

License: Apache-2.0.

## system/graph — coordinate systems for rendering (added in 0.1.2)

Package path `ky678/renky/system/graph`. Linear algebra (Vec/Quat/Mat/Euler/Rigid), the full
model→world→view→clip→NDC→screen pipeline (OpenGL / D3D / reversed-Z, finite or infinite far),
picking rays, frusta, WGS-84 geodesy (ECEF/ENU/NED), Web Mercator tiles, cube/equirect/octahedral
mappings, tangent frames and convention conversion. 75 tests (numpy/scipy/PROJ oracles) plus a
71-mutant mutation check. Details: [system/graph/README.md](system/graph/README.md).

```
moon test -p ky678/renky/system/graph
```
