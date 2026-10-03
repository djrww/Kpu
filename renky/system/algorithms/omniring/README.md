# omniring — 公理書（Axiom Book）

> **MoonBit** 實作：以**有限域**為地基的定義骨架 + **即時理想狀態**（live ideal state），
> 並把指定的 **8 個代數結構**各自「**轉熱**」（Maslov 去量子化 → 熱帶半環）與「**同調**」
> （鏈複形 → 同調／上同調），串成一條**終極運算力閉環鏈條**。
>
> **無第 3 方依賴**、**含定理／定義／算式**、**義務自證**（機器證明 + 窮盡驗證）、
> **演算法 89 條（> 18）**。

```
moon check      →  0 errors
moon test       →  Total tests: 122, passed: 122, failed: 0
moon prove      →  2 of 2 packages proved, 9 goals proved
moon run system/algorithms/omniring/cmd/report  →  印出公開帳本（closed = true）
```

---

## 1. 命案（驗收條件）逐條對帳

| 命案 | 交付 | 證據 |
|---|---|---|
| 用 **MoonBit** 實作 | 全專案 MoonBit（`moon.mod` / `moon.pkg` / `.mbt` / `.mbtp`），`preferred_target = wasm-gc` | 18 個套件 |
| **有限域**的定義骨架 | `kernel`（模運算，機器證明）、`arith`、`ff`（`F_p`、`F_{p^n}`、不可約、因式分解、離散對數）、`linalg`（`F_p` 上 rref/nullspace、`Z` 上 HNF/SNF） | `docs/AXIOMS.md` §0–§1 |
| **即時理想狀態** | `ideal.IdealState`：每次 `add(g)` 後**精確重算** Gröbner 基、零維判準、標準單項式基、`dim`、Krull 維度、Hilbert 函數 | `docs/AXIOMS.md` §2 |
| **8 個結構** | ① 仿射半群環 ② 上同調環 ③ 熱帶半環 ④ 光滑函數環 ⑤ 完美的 F-分裂環 ⑥ 分次環 ⑦ p 進整數環 ⑧ 理想類群對應環 | `semigroup/ cohom/ tropical+heat/ smooth/ fsplit/ graded/ padic/ classgroup/` |
| **轉熱** | Maslov 去量子化 `φ_h(x) = −h ln x`、`softmin_h → min`、Legendre `F = ⟨E⟩ − hS`、殘餘熵 `ln g`；閉環中另以 **J-adic 賦值** `v : R → (N∪{∞}, min, +)` 把每個環熱帶化 | `heat/`, `tropical/`, `loop/` A87–A88 |
| **同調** | 通用鏈複形引擎（`F_p` 用高斯消去、`Z` 用 Smith 標準形）、Koszul、張量積（Künneth）、映射錐、圖複形；上同調環（杯積、Poincaré 對偶、de Rham DG 代数） | `homology/`, `cohom/` |
| **終極運算力閉環鏈條** | `loop.run_all(p)`：結構 → 有限域資料 → 熱帶化 → 同調 → 不變量 → 約化 `R ↦ R/rad(0)`（冪等）→ 不動點；`closed = true` | `loop/`（L1–L6） |
| **定理** | 90+ 條定理／法律，分組 T/G/S/P/M/H/C/X/K/L | `docs/THEOREMS.md` |
| **定義** | 每個結構的定義（含結構常數、賦值、約化、對應） | `docs/AXIOMS.md` |
| **算式** | 每條定義旁的封閉算式（Hilbert–Serre、Ehrhart、Legendre、Kronecker、genus `2^{t−1}`、Faà di Bruno、Newton 多邊形、`χ = Σ(−1)^i β_i` …） | `docs/AXIOMS.md` |
| **義務自證** | (a) **Formal**：`moon prove` 9 goals（Why3 + Z3）；(b) **Exhaustive**：有限域上**完全枚舉**（例：`fsplit` 對 288 個環驗證 `F-split ⟺ reduced ⟺ perfect`；`classgroup` 對所有三元組驗證結合律）；(c) **Structural**：多條獨立路徑互驗 | `moon test` / `moon prove` |
| **無第 3 方依賴** | 只用 MoonBit 隨附的 `moonbitlang/core`（自動注入）；`Double` 缺 exp/log ⇒ 自行以級數實作（A44/A45）；Z3 僅為**證明器**（非程式依賴） | `moon.mod` 無 `deps` |
| **演算法 > 18 條** | **89 條**，A01–A89，全部自行實作、逐條標注於原始碼 docstring | `docs/ALGORITHMS.md`（由原始碼 grep 產生） |

---

## 2. 套件地圖

```
omniring/
├─ kernel/      A01–A05  modpos addmod submod mulmod negmod     ← proof-enabled（Formal）
├─ select/      A06–A07  argmin bsearch                         ← proof-enabled（Formal）
├─ arith/       A08–A17  gcd egcd valuation powmod invmod Miller–Rabin Pollard-rho
│                        jacobi Tonelli–Shanks crt divisors totient mobius
├─ linalg/      A27–A30  rref nullspace hnf snf  (+ rank, det_mod, solve_mod, mat_* )
├─ ff/          A18–A26  F_p[x] 運算、gcd、powmod、BSGS dlog、Ben-Or、不可約搜尋、
│                        ddf、Berlekamp edf、完全因式分解；F_{p^n} 表示
├─ poly/        A31–A34  多變數 MPoly、mp_divmod、S-polynomial、Buchberger、標準單項式
├─ ideal/       A35–A43  **即時理想狀態**：Krull 維度、成員+證書、Hilbert 函數、簇、
│                        Seidenberg 根、極小多項式/squarefree、極大性、乘法矩陣、Rabinowitsch
├─ heat/        A44–A47  expd logd softmin/softmax、partition/free_energy/gibbs（轉熱）
├─ tropical/    A48–A51  Hungarian 熱帶行列式、min-plus 矩陣積與 Floyd–Warshall 閉包、
│                        Newton 細分、平衡條件
├─ graded/      A52–A56  Hilbert 級數/多項式、Veronese、理想冪、adic 分次、外(Koszul)代數
├─ semigroup/   A57–A62  半群成員 DP、HNF 格成員、正規性、凸包+膨脹點計數(Ehrhart)、
│                        Newton 二項式基插值、環面理想
├─ padic/       A63–A67  p 進數字運算與 Newton 逆、Teichmüller、Hensel、exp/log、Newton 多邊形
├─ smooth/      A68–A71  jet 合成(Faà di Bruno)、Hasse 除法冪導子、jet 熱流、Lagrange 插值
├─ homology/    A72–A76  鏈複形同調(F_p 與 Z)、張量積(Künneth)、映射錐、Koszul、圖複形+並查集
├─ cohom/       A77–A78  分次交換杯積 + Poincaré 配對；jet 代數的 de Rham DG 代数（同調）
├─ fsplit/      A79–A81  結構常數與 Frobenius 矩陣、仿射分裂系統求解、reduced/perfect 判定
├─ classgroup/  A82–A86  Gauss 約化、二次序運算+理想乘法+型↔理想對應、Gauss 合成、
│                        Kronecker 符號+genus 特徵、群代數 F_p[Cl(D)]（對應環）
├─ loop/        A87–A89  J-adic 賦值（熱帶化）、Maslov 橋、**閉環執行器** run_all
└─ cmd/report/           主程式：印出公開帳本
```

---

## 3. 閉環鏈條（終極運算力）

```
     ①有限域資料            ②熱帶化(轉熱)          ③同調               ④不變量
 X ─────────────▶ R=F_p[x]/I ─────────▶ (R, v_J, min,+) ─────▶ Koszul 0→R─·a→R→0 ─────▶ inv
 (八結構)          dim, H(d), kdim,       v(x+y) ≥ min(vx,vy)     H_1 = ann(a)          100·Loewy
                   reduced, F-split,      v(xy) ≥ v(x)+v(y)       H_0 = R/aR              + dim
                   Loewy                  softmin_h(v) → min v    χ(C) = χ(H)
                                                                                    │
                     ⑤ 約化 R ↦ R/rad(0)（閉包算子，冪等）◀──────────────────────────┘
                        reduced ⟺ F-split ⟺ Loewy = 1  ⇒ 不動點，closed = true
```

`moon run system/algorithms/omniring/cmd/report` 的實際輸出（節錄，`p = 2`）：

```
closed loop over F_2
  1. 仿射半群環 | F_p[x,y,z]/(I_A + (x^2,y^2,z^2)) | dim=6 H=[1,3,2,0,…] kdim=0 reduced=false Fsplit=false Loewy=3 beta=[3,3] chi=0
  2. 上同調環   | H^*(CP^2) = F_p[h]/(h^3)         | dim=3 H=[1,1,1,0,…] kdim=0 reduced=false Fsplit=false Loewy=3 beta=[1,1] chi=0
  3. 熱帶半環   | F_p[x]/(x^p - x) = F_p^p         | dim=2 H=[1,1,0,…]   kdim=0 reduced=true  Fsplit=true  Loewy=1 beta=[1,1] chi=0
  4. 光滑函數環 | jet algebra F_p[t]/(t^p)         | dim=2 H=[1,1,0,…]   kdim=0 reduced=false Fsplit=false Loewy=2 beta=[1,1] chi=0
  5. 完美的F分裂環 | F_{p^2}                        | dim=2 H=[1,1,0,…]   kdim=0 reduced=true  Fsplit=true  Loewy=1 beta=[0,0] chi=0
  6. 分次環     | F_p[x,y]/(x^2,y^2)               | dim=4 H=[1,2,1,0,…] kdim=0 reduced=false Fsplit=false Loewy=3 beta=[2,2] chi=0
  7. 理想類群對應環 | F_2[Cl(-23)] = F_p[x]/(x^3-1) | dim=3 H=[1,1,1,0,…] kdim=0 reduced=true  Fsplit=true  Loewy=1 beta=[0,0] chi=0
  [ok]  …每結構 8 條交叉檢查（Σ_d H(d)=dim R、χ(C)=χ(H)、reduced⟺F-split⟺Loewy1、
          熱帶賦值公理、Maslov 橋、rank–nullity、kdim=0）…
  [ok]   p進: valuation axioms (Ostrowski) / Teichmuller / exp-log / Hensel / p-torsion
  [ok]   reduction R -> R/rad(0) is idempotent
  [ok]   上同調環 = H^*(CP^2): the regraded Hilbert function = the Betti numbers, chi = 3
  [ok]   上同調環: Poincare duality and associativity
  [ok]   光滑函數環 = J_(p-1): t^p = 0, t^(p-1) != 0 and dim = p
  [ok]   光滑函數環: jet exp/log are inverse and d/dt obeys Leibniz
  closed = true

  h(-3) = 1   (1,1,1)        genera = 1        h(-23) = 3  (1,1,6)(2,-1,3)(2,1,3)   genera = 1
  h(-15) = 2  (1,1,4)(2,1,2) genera = 2        h(-47) = 5  …                          genera = 1
  h(-20) = 2  (1,0,5)(2,2,3) genera = 2        h(-71) = 7  …                          genera = 1
  h(-24) = 2  (1,0,6)(2,0,3) genera = 2        h(-163) = 1 (1,1,41)                   genera = 1

  p=2 F_p[x]/(x^2)      : dim=2 reduced=false F-split=false perfect=false
  p=2 F_p[x]/(x^2-x)    : dim=2 reduced=true  F-split=true  perfect=true
```

---

## 4. 義務自證的兩層

**(a) Formal — `moon prove`（Why3 + Z3）**
`kernel`、`select` 啟用 `options("proof-enabled": true)`，以 `where { proof_require / proof_ensure }`
與 `.mbtp` 謂詞契約寫成，迴圈帶 `proof_invariant` / `proof_decrease`。9 個目標全數 discharged：

```
ky678/renky/system/algorithms/omniring/kernel   Succeeded: 6 goals proved
ky678/renky/system/algorithms/omniring/select   Succeeded: 3 goals proved
Summary: 2 of 2 packages proved, 9 goals proved
```

**(b) Exhaustive / Structural — `moon test`（122 laws）**
凡屬**有限**的命題一律**完全枚舉**，不抽樣：

- `fsplit` X2：`F_p[x,y]/(x^a,y^b,M)` 的**所有** `M ⊆ box`（`p ∈ {2,3}`，`(a,b) ∈ {(2,2),(2,3),(3,2)}`，288 個環）
  ⇒ `F-split ⟺ reduced ⟺ perfect ⟺ rank(F) = dim`。
- `fsplit` X3：枚舉**所有** `dim²` 個矩陣，直接驗證分裂公理，`#分裂 = p^{nullity}`。
- `classgroup` K2：所有三元組的結合律；K6：所有類 × 所有 genus 特徵。
- `homology` H3：**4 個頂點的全部 64 個圖**；`arith`：Fermat 小定理對所有質數 `p ≤ 100` 的每個 `a`。
- `ideal` T6：有限域 Nullstellensatz 對所有 `deg ≤ 2` 的 `f`。

**數學上被測試抓出的陷阱**（每条都曾导致真实失败，现已成为定理的一部分）：
`Array::make(n, Array::make(...))` 会共享同一个内层数组；分次交换律与 Leibniz **只对齊次元素**成立；
`d/dt` 下降到 `F_p[t]/(t^n)` **若且唯若** `p | n`；`softmin_h(a,a,h) = a − h ln 2`（有限 `h` 下不冪等）；
胞腔維度 ≠ 指標個數（須以 `hull_dim` 分類）；平衡貢獻已含格長度（勿再乘 `gcd`）；
`F_p[x,y]/(x^4)` 是**無限維**（測試理想必須界住每個變數）；型 ↔ 理想對應中 `b` 的移位必須同步重算 `c`。

---

## 5. 環境與指令

```bash
# 工具鏈（沙箱重置後需重裝；注意路徑是 nightly，core 解壓時「不要」strip）
export MOON_HOME=$HOME/.local/moonbit
curl -sL https://cli.moonbitlang.com/binaries/nightly/moonbit-linux-x86_64.tar.gz \
  | tar xz --strip-components=1 -C "$MOON_HOME"
curl -sL https://cli.moonbitlang.com/cores/core-nightly.tar.gz | tar xz -C "$MOON_HOME/lib"
chmod +x "$MOON_HOME"/bin/*
moon -C "$MOON_HOME/lib/core" bundle --warn-list -a --all
moon -C "$MOON_HOME/lib/core" bundle --warn-list -a --target wasm-gc --quiet

# 證明器（僅 moon prove 需要）
curl -sL -o /tmp/z3.zip https://github.com/Z3Prover/z3/releases/download/z3-5.1.0/z3-5.1.0-x64-glibc-2.39.zip
unzip -q /tmp/z3.zip -d /tmp/z3x && cp /tmp/z3x/*/bin/z3 "$HOME/.local/z3/"

export PATH="$MOON_HOME/bin:$HOME/.local/z3:$PATH"
cd renky   # module ky678/renky; sources under system/algorithms/omniring
moon check            # 0 errors
moon test             # 122 passed
moon prove            # 9 goals proved
moon run system/algorithms/omniring/cmd/report   # 印出帳本
```

---

## 6. 文件

| 檔案 | 內容 |
|---|---|
| `docs/AXIOMS.md` | **公理與定義（含算式）**：底層模運算公理、有限域骨架、即時理想狀態、8 大結構的定義與封閉算式、閉環定義 |
| `docs/THEOREMS.md` | **定理帳本**：Formal 9 條 + 各套件法律（T/G/S/P/M/H/C/X/K/L 分組），逐條標注證明方式 |
| `docs/ALGORITHMS.md` | **演算法登記冊 A01–A89**（由原始碼標籤自動產生，與程式碼一一對應） |
