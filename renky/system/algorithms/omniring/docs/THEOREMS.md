# THEOREMS · 定理帳本（義務自證）

> 每一條定理都標注**證明方式**：
> - **Formal** — `moon prove`（Why3 + Z3）機器 discharged 的驗證條件。
> - **Exhaustive** — 在**整個有限域上完全枚舉**驗證（無抽樣、無隨機）。
> - **Structural** — 由已證引理組合而成（在測試中以多條獨立路徑互驗）。
>
> 執行：`moon test`（122 條法律，全數通過）與 `moon prove`（9 個目標，全數 proved）。

---

## 0. Formal layer（9 goals, `moon prove` → *2 of 2 packages proved, 9 goals proved*）

| # | 目標 | 敘述 | 方式 |
|---|---|---|---|
| 1 | `kernel.modpos` | `0 ≤ r < m ∧ r ≡ x (mod m)` | Formal |
| 2 | `kernel.addmod` | `addmod(a,b,m) ≡ a+b`，且在 `[0,m)` | Formal |
| 3 | `kernel.submod` | `submod(a,b,m) ≡ a−b`，且在 `[0,m)` | Formal |
| 4 | `kernel.mulmod` | `mulmod(a,b,m) ≡ a·b`，且在 `[0,m)` | Formal |
| 5 | `kernel.negmod` | `negmod(a,m) ≡ −a`，且在 `[0,m)` | Formal |
| 6 | `kernel.cong_refl` | `∀ x m. cong(x,x,m)` | Formal |
| 7 | `select.argmin` | `∀ i ∈ [lo,hi). xs[argmin] ≤ xs[i]` | Formal |
| 8 | `select.bsearch` | `Some(i) ⇒ xs[i] = x`；`None ⇒ x ∉ xs` | Formal |
| 9 | `select.argmin_is_min` | argmin 的結果確實是最小值 | Formal |

---

## 1. `arith` — 18 laws（`arith/laws_wbtest.mbt`）

| 定理 | 敘述 | 方式 |
|---|---|---|
| gcd 正確性 | `gcd(a,b)` = 暴力最大公因數，`[1,60]²` 全枚舉 | Exhaustive |
| Bézout | `a·x + b·y = gcd(a,b)`，`[-40,40]²` 全枚舉 | Exhaustive |
| `lcm·gcd = |ab|` | `[1,40]²` | Exhaustive |
| powmod | = 朴素重複乘法，`p ∈ {5,7,13,101}`, `e ≤ 40` | Exhaustive |
| **Fermat 小定理** | `a^p ≡ a (mod p)`，所有 `a ∈ Z/p`、所有質數 `p ≤ 100` | Exhaustive |
| invmod | `a·a^{-1} ≡ 1`，`Z/p` 的每個單位、`p ≤ 60` | Exhaustive |
| Miller–Rabin | 與試除法一致，所有 `n ≤ 2000` | Exhaustive |
| 因數分解 | 質數冪之積還原 `n`，所有 `n ≤ 600` | Exhaustive |
| 賦值 | `p^v ‖ n`，`n ≤ 200`, `p ∈ {2,3,5,7}` | Exhaustive |
| digits 往返 | `from_digits(digits(n,b)) = n`，`n ≤ 500`, `b ∈ [2,12]` | Exhaustive |
| totient | `φ(n)` = 互質計數，`n ≤ 120` | Exhaustive |
| **Möbius 反演** | `Σ_{d|n} μ(d) = [n=1]` 且 `φ(n) = Σ_{d|n} μ(d)·(n/d)`，`n ≤ 120` | Exhaustive |
| Jacobi／Legendre | `(a/p)` 與 Euler 判準 `a^{(p-1)/2}` 一致 | Exhaustive |
| Tonelli–Shanks | `sqrt_mod(a,p)² ≡ a`，對所有 QR；非 QR 回傳 `-1` | Exhaustive |
| CRT | 中國剩餘定理對所有互質模對 | Exhaustive |
| 同餘保持 | `cong` 對 `+ − ×` 保持（機器證明之補） | Exhaustive |

## 2. `linalg` — 7 laws

| 定理 | 敘述 | 方式 |
|---|---|---|
| **rank–nullity** | `rank(m) + dim ker(m) = cols(m)`，所有 `3×3` 矩陣 over `F_2, F_3` | Exhaustive |
| rref | 冪等且保持列空間，所有 `2×3` over `F_2,F_3` | Exhaustive |
| det | `det_mod` = Leibniz 排列公式，所有 `3×3` over `F_2,F_3` | Exhaustive |
| solve | 回傳解確實滿足方程組，所有 `2×2` over `F_3` | Exhaustive |
| **Smith 標準形** | `d_i | d_{i+1}`、秩相符、`det` = 對角積；已知矩陣對照 | Exhaustive + Structural |
| **Hermite 標準形** | 上三角、樞軸非負、生成同一格 | Exhaustive |
| 矩陣乘法 | 結合律與分配律，所有 `2×2` over `F_2,F_3` | Exhaustive |

## 3. `ff` — 18 laws

| 定理 | 敘述 | 方式 |
|---|---|---|
| `F_p` 環公理 | 所有三元組，`p ∈ {2,3,5,7}` | Exhaustive |
| Frobenius on `F_p` | `a^p = a`，所有 `a`、所有質數 `p ≤ 100` | Exhaustive |
| **`F_p^*` 循環** | 生成元存在，且恰有 `φ(p−1)` 個，`p ≤ 60` | Exhaustive |
| `F_p[x]` 除法 | `f = qg + r`, `deg r < deg g`，所有次數 `< 4` over `F_2,F_3` | Exhaustive |
| Bézout（`F_p[x]`） | `sf + tg = gcd(f,g)`，所有次數 `< 4` 之對 | Exhaustive |
| **Ben-Or** | 與暴力不可約判定一致，所有 monic `f`，`deg ≤ 4` over `F_2,F_3` | Exhaustive |
| **Gauss 計數** | 次數 `d` 的 monic 不可約數 `= (1/d)Σ_{k|d} μ(d/k)p^k` | Exhaustive |
| **有限域基本定理** | `x^{p^n} − x = ∏`（所有次數整除 `n` 的 monic 不可約） | Structural |
| 完全因式分解 | 正確，所有 monic `deg ≤ 4` over `F_2`、`deg ≤ 3` over `F_3` | Exhaustive |
| ddf + edf | 重現 squarefree 部分 | Exhaustive |
| `F_{p^n}` 環公理 | 所有三元組：`F_4, F_9, F_25` | Exhaustive |
| Freshman's dream | `(a+b)^p = a^p + b^p`；`F^n = id` | Exhaustive |
| 離散對數 | BSGS 與暴力枚舉一致 | Exhaustive |
| 範數 | `N(c·a) = c^n N(a)`；`N(ab) = N(a)N(b)` | Exhaustive |
| 極小多項式 | 共軛之積；`F_{p^n}` 上每個元素的極小多項式整除 `x^{p^n} − x` | Structural |

## 4. `ideal` — 13 laws（T1–T8）

| 定理 | 敘述 | 方式 |
|---|---|---|
| **T1** | `f ∈ I ⟺ nf(f) = 0`，並附證書 `f = Σ q_i g_i`（驗證後代入） | Structural |
| **T2** Buchberger | 所得基的每個 S-多項式對該基歸約為 0 | Exhaustive |
| **T3** | 標準單項式在 `R/I` 中線性獨立（且張成）⇒ 是 `F_p`-基 | Exhaustive |
| **T4** | `dim_{F_p} R/I = Σ_d H(d)` | Structural |
| **T5** | Krull 維度：零維 ⟺ 商有限維（一組理想目錄上） | Exhaustive |
| **T6** 有限域 Nullstellensatz | `I(V(I)) = rad(I + (x_i^p − x_i))`，所有 `deg ≤ 2` 之 `f`；且 `f ∈ I ⟺ f 在 V(I) 上消失` | Exhaustive |
| **T7** Rabinowitsch | 根成員關係與根狀態一致，`I = (x²,y²)` over `F_2,F_3` | Exhaustive |
| **T8** | `I` 極大 ⟺ `R/I` 無零因子（⟺ 是域） | Exhaustive |
| 環公理／求值同態 | `F_p[x,y]` 是交換環；`ev_a : F_p[x,y] → F_p` 是環同態（所有點） | Exhaustive |
| 除法演算法 | `f = Σ q_i g_i + r`，`r` 無可約項 | Exhaustive |
| 正規形典型性 | `nf` 不依賴基的順序 | Exhaustive |
| **即時性** | 逐一加入生成元 = 一次加入全部（增量 Gröbner 與批次一致） | Exhaustive |

## 5. `heat` — 6 laws（轉熱層）

| 定理 | 敘述 | 方式 |
|---|---|---|
| `log ∘ exp = id` | 在整個格點上 | Exhaustive |
| `exp` 是同態 | `(R,+) → (R_{>0},×)`，所有格點對 | Exhaustive |
| `φ_h` 是 semifield 同構 | 熱加法 ↔ 普通加法、熱乘法 ↔ 普通乘法 | Exhaustive |
| 熱加法律 | 交換、結合、平移不變、`h→0` 時冪等 | Exhaustive |
| **零溫極限** | `softmin_h → min`，且 `0 ≤ min − softmin_h ≤ h·ln 2`；`softmin_h` 隨 `h` 遞減而單調遞增 | Exhaustive |
| **Legendre 關係** | `F = ⟨E⟩ − hS`；`F ≤ min E`；`lim_{h→0}F = min E`；`S(h→0) = ln g`；Gibbs `argmax` = 熱帶 `argmin` | Exhaustive |

## 6. `tropical` — 9 laws

| 定理 | 敘述 | 方式 |
|---|---|---|
| `T = Z ∪ {∞}` 是冪等交換半環 | `(min,+)` 全部公理，所有三元組 | Exhaustive |
| 求值是半環同態 | `TPoly → T`，對 box 內每個 `w` | Exhaustive |
| 熱帶多項式律 | 加法／乘法／冪等在全族上 | Exhaustive |
| 熱帶超曲面 | `x ⊕ y ⊕ 1` 的重數 `≥ 2` 恰在三條射線上 | Exhaustive |
| **A48 熱帶行列式** | Hungarian = 所有排列之最小值（`n ≤ 5`） | Exhaustive |
| Hungarian vs 暴力 | 隨機矩陣 `n = 4,5,6`，且回傳一個排列 | Exhaustive |
| min-plus 閉包 | Floyd–Warshall = 路徑枚舉 | Exhaustive |
| **A50/A51 Newton 細分與平衡** | 胞腔維度、`Σ len(E)·ν_E = 0`（一組雙變熱帶曲線目錄） | Exhaustive |
| Maslov 橋 | `h → 0` 的單調收斂至熱帶極限 | Exhaustive |

## 7. `graded` — 5 laws（G1–G7）

| 定理 | 敘述 | 方式 |
|---|---|---|
| **G1/G2/G3** | `R/I` 的結構常數是分次的、結合的、交換的、含幺的 | Exhaustive |
| **G4 Hilbert–Serre** | Hilbert 函數最終是次數 `kdim − 1` 的多項式（三種獨立算法互驗） | Structural |
| **G5 Veronese** | `H_{A^{(k)}}(d) = H_A(kd)`，且子代數律成立 | Exhaustive |
| **G6 adic 分次** | `Σ_k dim(J^k/J^{k+1}) = dim R`，各分片非負 | Exhaustive |
| **G7 外（Koszul）代數** | 分次交換、結合、`a∧a = 0`、`dim Λ^d = C(n,d)` | Exhaustive |

## 8. `semigroup` — 5 laws（S1–S6）

| 定理 | 敘述 | 方式 |
|---|---|---|
| **S1/S2** | 半群成員（DP）與格成員（HNF）皆與窮盡搜尋一致 | Exhaustive |
| **S3 正規性** | 古典例子的真／假判定 | Exhaustive |
| **S4** | `k[S]` 的 Hilbert 函數；標準單項式 ↔ 半群元素雙射（經 `toric_ideal`） | Structural |
| **S5 Ehrhart** | 多項式性、首項 `= 2!·vol`、Macdonald 互反 vs `dilated_interior` | Exhaustive |
| **S6 環面理想** | 二元式在單項式映射下消失；格關係 `∈ I_A` | Exhaustive |

## 9. `padic` — 7 laws（P1–P8）

| 定理 | 敘述 | 方式 |
|---|---|---|
| **P1** | `Z/p^k` 是交換環，且截斷是環同態（截斷相容性） | Exhaustive |
| **P2/P3** | `v(xy) = v(x)+v(y)`；`v(x+y) ≥ min`；等號情形；Ostrowski 超度量 | Exhaustive |
| **P4** | `Z_p` 是 DVR：有限個元素生成的理想 `= (p^{min v})` | Exhaustive |
| **P5** | Teichmüller：`ω(a)^{p−1}=1`、`ω(a) ≡ a (mod p)`、`ω(ab)=ω(a)ω(b)`、單位分裂 | Exhaustive |
| **P6** | Hensel：mod `p` 的單根唯一提升到 `Z/p^k` 的根 | Exhaustive（含唯一性） |
| **P7** | `exp`／`log` 在極大理想上互逆且為同態（`p = 3`） | Exhaustive |
| **P8** | Newton 多邊形斜率 = 根的賦值；`Σ dx = deg`；`Σ dy = v(a_n) − v(a_0)`；`x²−p` 斜率 `−1/2` | Exhaustive |

## 10. `smooth` — 7 laws（M1–M6）

| 定理 | 敘述 | 方式 |
|---|---|---|
| **M1** | `J_k` 是交換局部 `F_p`-代數，極大理想冪零（`t^{k+1}=0`）；單位 ⟺ `a_0 ≠ 0` | Exhaustive |
| **M2** | Kock–Lawvereux 與每個點的 Taylor 展開 | Exhaustive |
| **M3** | 導子：Leibniz、線性、`[d/dt, E] = d/dt`、Jacobi 恆等式 | Exhaustive |
| **M4** | jet 合成 = 代入；結合、單位、鏈式法則、Faà di Bruno | Exhaustive |
| **M5** | jet 熱流：指數律、閉式、熱方程 `∂_s u = Δu` | Exhaustive |
| **M6** | 每個 `F_p → F_p` 的映射都是次數 `< p` 的多項式，且唯一（`p ∈ {2,3}`） | Exhaustive |
| jet exp/log | 在 unipotent／nilpotent jets 上互逆 | Exhaustive |

## 11. `homology` — 7 laws（H1–H6）

| 定理 | 敘述 | 方式 |
|---|---|---|
| **H1** | `χ(C) = χ(H)`，隨機 `d² = 0` 複形；每步 rank–nullity | Structural |
| **H2** | `F_p[x]/(x^n)` 上 `x` 的 Koszul 同調：`H_0 = H_1 = F_p` | Exhaustive |
| **H3** | 圖同調：`β_0` = 連通分量、`β_1` = 圈秩，**4 頂點的全部 64 個圖** | Exhaustive |
| **H4 Künneth** | `β_k(A ⊗ B) = Σ_{i+j=k} β_i(A)β_j(B)`；`χ` 相乘 | Exhaustive |
| **H5** | `cone(id)` 零調；`cone(0) = C ⊕ C[−1]`；長正合序列交替和為 0 | Structural |
| **H6** | `Z`-扭轉由 SNF 讀出（`[2]`、`[6]` 之 CRT、`Z/2 ⊕ Z`、三項複形）；分裂短正合序列中 `χ`／Betti 可加 | Exhaustive |
| SNF 不變因子 | 已知矩陣對照 | Exhaustive |

## 12. `cohom` — 4 laws（C1–C6）

| 定理 | 敘述 | 方式 |
|---|---|---|
| **C1/C2/C3** | `S^n, T², CP², RP²(F_2), S²×S²` 的環律（結合、次數可加、分次交換、單位）、Poincaré 對偶（配對矩陣之秩）、對稱 Betti 數；`χ = 0,0,0,0,3,1,4`；`a²=0`、`a∪b=[T²]`、`b∪a=−a∪b`、`h²=[CP²]`、`x²=[RP²]` | Exhaustive |
| **C4** | jet 代數的 de Rham DG 代数（`p | n` 時 `d/dt` 下降）：`d² = 0`、楔積結合／雙線性／分次交換／含幺、`dt∧dt = 0`、**齊次元素**上的分次 Leibniz | Exhaustive |
| **C5** | DG 代數之上同調是分次環：cocycles 對 `∧` 封閉、coboundary × **exact** = coboundary、`dr_well_defined(p,n) ⟺ p | n`、`dt = d(t)` 是 exact、`t^{n−1}dt` 生成 `H¹` 且非 exact、`[1]∪[t^{n−1}dt] = [t^{n−1}dt]` | Exhaustive |
| **C6** | `dr_cohomology_dims`（導子矩陣之秩）`= dr_cohomology_count`（計數公式），`p ∈ {2,3,5}`, `n ≤ 6`；`χ(Ω) = χ(H) = 0`；獨立枚舉：`#closed 0-forms = p^{h_0}`、`#nonzero exact 1-forms = p^{rank} − 1` | Exhaustive（兩條獨立路徑） |

## 13. `fsplit` — 5 laws（X1–X6）

| 定理 | 敘述 | 方式 |
|---|---|---|
| **X1** | 分裂條件 `(i)(ii)` 是 `φ` 矩陣元素的**仿射線性系統**，以高斯消去判定 ⇒ F-splitness 可判定 | Structural |
| **X2** | 對**全部** `F_p[x,y]/(x^a, y^b, M)`（`M` 任取 box 內單項式子集，`p ∈ {2,3}`, `(a,b) ∈ {(2,2),(2,3),(3,2)}`，共 288 個環）：`F-split ⟺ reduced ⟺ perfect ⟺ rank(F) = dim R`；兩種 reducedness 判法一致 | Exhaustive |
| **X3** | 分裂集合是 torsor：`φ₀ + {ψ : ψ R-線性, ψ∘F = 0}`，且 `#分裂 = p^{nullity}`（與**所有矩陣**的暴力枚舉一致） | Exhaustive |
| **X4** | `F_p^n = F_p[x]/(x^n − x)`：`F = id`、唯一分裂 `φ = id`（nullity 0）、理想恰有 `2^n` 個、每個理想由冪等元生成且被 `φ` 相容分裂 | Exhaustive |
| **X5** | Frobenius 是環同態，且 `F(a) = a^p`（矩陣路徑 vs 重複乘法路徑） | Exhaustive |
| **X6** | `F_{p²}`（`p ∈ {2,3,5}`，由 A23 找到的不可約二次式）：reduced、perfect、唯一分裂且 `φ = F = F^{-1}`、`φ(F(x)) = x` 逐點成立 | Exhaustive |

## 14. `classgroup` — 7 laws（K1–K7）

| 定理 | 敘述 | 方式 |
|---|---|---|
| **K1 Gauss 約化** | 保判別式、冪等、`SL_2(Z)`-不變（`b → b+2ak`、`(a,b,c) → (c,−b,a)`）；每個 primitive 型約化後落在 reduced 清單中（`1 ≤ a ≤ 6`, `|b| ≤ 14`, `1 ≤ c ≤ 20` 全枚舉）；同類的 reduced 型相等（代表元唯一） | Exhaustive |
| **K2** | `Cl(D)` 在合成下是有限交換群：封閉、交換、結合（所有三元組）、單位 = 主型、`(a,b,c)^{-1} = (a,−b,c)`，`D ∈ {−15,−20,−23,−24,−31,−47}` | Exhaustive |
| **K3 類數** | `h(−3)=h(−4)=h(−7)=h(−8)=h(−11)=h(−19)=h(−163)=1`、`h(−15)=h(−20)=h(−24)=h(−40)=2`、`h(−23)=h(−31)=3`、`h(−47)=5`、`h(−71)=7` | Exhaustive（古典值對照） |
| **K4 Lagrange** | `g^{h} = 1`；每個類的階整除 `h`；階極小；`⟨g⟩` 恰有 `ord(g)` 個元素 | Exhaustive |
| **K5 對應** | `N(I(f)) = a`；`N(IJ) = N(I)N(J)`（格的行列式）；`q` ramified ⟺ `q | D` ⟺ `(D/q) = 0`；split ⟺ `(D/q) = 1` ⟺ `b² ≡ D (mod 4q)` 可解；對 split 的 `q`，`(q,b,(b²−D)/4q)` 是 primitive、判別式為 `D`、落在類群中，其共軛是其逆，且階整除 `h` | Exhaustive |
| **K6 genus 理論** | discriminant divisor 有 `2^t − 1` 個；genus 數 `= 2^{t−1}`；每個 genus 特徵取值 `{±1}`、在類上**良定義**（不同代表元同值）、是群同態；所有 genus 特徵之積平凡；各 genus 大小相等 | Exhaustive |
| **K7 對應環** | `F_p[Cl(D)]` 是交換結合含幺環（所有元素對／三元組）；分次：`[g][h] = [gh]`，類群表每列是排列；群元素皆為單位；augmentation 是環同態 | Exhaustive |

## 15. `loop` — 4 laws（L1–L6）：閉環

| 定理 | 敘述 | 方式 |
|---|---|---|
| **L1** | 對 `p = 2, 3`，`run_all(p)` 的**每一條**交叉檢查通過 ⇒ `closed = true`（含八結構之 `Σ_d H(d) = dim R`、`χ(C) = χ(H)`、`reduced ⟺ F-split ⟺ Loewy = 1`、熱帶賦值公理、Maslov 橋、rank–nullity、`kdim = 0`、p-adic 賦值／Teichmüller／exp-log／Hensel／`p`-扭轉、`H^*(CP²)` 之 PD 與結合律、jet 代數 `t^p = 0`） | Structural（整合） |
| **L2** | 每個 stage 都是有限維 Artin `F_p`-代數，且各不變量彼此一致 | Exhaustive |
| **L3** | 八大結構全部出現在鏈條中 | Structural |
| **L4 冪等** | `R ↦ R/rad(0)` 是閉包算子：`R = F_3[x,y]/(x³,y²)` 非 reduced、Loewy `> 1`；約化後 reduced、Loewy `= 1`、F-split；再約化不變（`dim`、Loewy、`χ`、Hilbert 函數、旗標、不變量全部穩定）⇒ 不動點 | Exhaustive |
| **L5 Maschke** | `F_p[Cl(−23)]`：`p = 2` 半單（reduced、Loewy 1）；`p = 3 | h` 非 reduced、Loewy `= v_p(h)+1 = 2`；`F_5[Cl(−47)]` 同理 | Exhaustive |
| **L6 跨引擎一致** | 截斷環面環之 `dim` = 標準單項式數 = 有界次數之半群元素數（獨立枚舉）；`F_p^p` 半單（`(min,+)` 冪等性的影子）；熱帶半環律；Maslov 橋之界與單調性與殘餘熵；jet 代數 Loewy `= p` | Exhaustive |

---

## 總計

| 層 | 數量 | 指令 |
|---|---|---|
| Formal（Why3 + Z3） | **9 goals proved**（2/2 packages） | `moon prove` |
| Executable（窮盡／結構） | **122 tests passed, 0 failed** | `moon test` |
| Algorithms（自製、無第三方） | **89 條（A01–A89）** | `docs/ALGORITHMS.md` |
| 結構 | 8/8 + 有限域骨架 + 即時理想狀態 + 閉環 | `loop.run_all(p)` |
