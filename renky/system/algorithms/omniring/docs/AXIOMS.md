# AXIOMS · 公理與定義（含算式）

> 本文件是 `omniring` 的**公理層**。每一條公理／定義都對應到程式碼中的一個型別或函式，
> 且都带有**義務自證**（`moon prove` 機器證明 或 `moon test` 窮盡驗證）。
> 對應的定理清單見 `THEOREMS.md`，演算法登記冊見 `ALGORITHMS.md`。

---

## 0. 底層公理： residue arithmetic（`kernel`, `select`）

這一層是**唯一啟用機器證明**（`options("proof-enabled": true)` + `.mbtp` 契約）的層，
因為它是一切模運算的地基。Why3 + Z3  discharged 9 個目標。

| 公理 | 算式 | 契約 |
|---|---|---|
| **A01** `modpos` | `0 ≤ modpos(x,m) < m` 且 `modpos(x,m) ≡ x (mod m)` | `proof_ensure: result => 0 <= result && result < m && cong(result, x, m)` |
| **A02** `addmod` | `addmod(a,b,m) = modpos(a+b,m)` | 同上，含 `cong` 保持 |
| **A03** `submod` | `submod(a,b,m) = modpos(a-b,m)` | 同上 |
| **A04** `mulmod` | `mulmod(a,b,m) = modpos(a*b,m)` | 同上 |
| **A05** `negmod` | `negmod(a,m) = modpos(-a,m)` | `addmod(a, negmod(a,m), m) == 0` |
| **A06** `argmin` | `∀ i ∈ [lo,hi). xs[argmin] ≤ xs[i]` | 迴圈不變式 + `proof_decrease` |
| **A07** `bsearch` | 有序陣列上 `bsearch(x) = Some(i) ⟺ xs[i] == x` | 不變式 `xs[lo] ≤ x < xs[hi]` 之補 |

**同餘公理（`cong_refl`）**：`∀ x m. cong(x, x, m)`（自反性）——已機器證明；
對稱／傳遞／加乘保持（`cong_sym/add/mul`）在 Z3 下是非線性死路，改由**窮盡層**驗證
（見 `THEOREMS.md` §1）。

---

## 1. 有限域骨架（`arith`, `ff`, `linalg`）

**定義 1.1（F_p）** 設 `p` 為質數。`F_p = Z/pZ`，運算由 A02–A05 給出。
*義務*：環公理（結合、交換、分配、單位、逆元）對**所有**三元組 `(a,b,c) ∈ F_p³` 窮盡成立。

**定義 1.2（Frobenius）** `F : F_{p^n} → F_{p^n}`, `F(a) = a^p`。
*算式*：freshman's dream `(a+b)^p = a^p + b^p`；`F^n = id`；`a^{p^n} = a`。

**定義 1.3（F_{p^n}）** `F_{p^n} = F_p[x]/(f)`，`f` 為 `n` 次不可約多項式（A23 搜尋、A22 Ben-Or 驗證）。
*算式*：`x^{p^n} - x = ∏_{d | n} ∏_{f 不可約, deg f = d} f`（有限域基本定理）；
不可約多項式個數 `N_d = (1/d) Σ_{k|d} μ(d/k) p^k`（Gauss）。

**定義 1.4（線性代數）** `rref`（A27）、`nullspace`（A28）在 `F_p` 上；`hnf`（A29）、`snf`（A30）在 `Z` 上。
*算式*：rank–nullity `rank(m) + dim ker(m) = cols(m)`；
Smith 標準形 `diag(d_1,…,d_r)` 且 `d_i | d_{i+1}`，即 `Z^m/(im d)` 的不變因子。

---

## 2. 即時理想狀態（`poly`, `ideal`）—— 實時 ideal state

**定義 2.1（`IdealState`）** 一個**活的**（live）零維理想狀態：

```
IdealState = { p, nvars, gens（追加式原始生成元）, gb（Gröbner 基）,
               dim, kdim, finite, var_bounds, basis（標準單項式）, … }
```

**即時性公理**：每次 `st.add(g)` 之後，**所有**導出不變量被*精確*重算（增量 Buchberger，A33），
不使用任何次數启发式：
- `finite ⟺ ∀ v ∃ k. x_v^k ∈ LM(I)`（零維判準，精確）
- `basis` = 標準單項式（不被任何 `LM(g)` 整除者），是 `R/I` 的 `F_p`-基（T3）
- `dim = |basis| = Σ_d H(d)`（T4）
- `kdim` = Krull 維度，由初始理想算出（T5，A35）

**定義 2.2（單項式序）** grlex：先比總次數 `total_deg`，再字典序。
**定義 2.3（Gröbner 基）** `G` 是 `I` 的 Gröbner 基 ⟺ `∀ f,g ∈ G. S(f,g) →_G 0`（Buchberger 判準，T2）。
**定義 2.4（正規形）** `nf(f)` = 對 `G` 反覆 `mp_divmod`（A31）至不可約；`f ∈ I ⟺ nf(f) = 0`（T1），
並附**證書** `f = Σ q_i g_i`（A36 `membership_witness`）。
**定義 2.5（Hilbert 函數）** `H(d) = dim_{F_p} (R/I)_d`（A37）。
**定義 2.6（簇）** `V(I)(F_p) = { a ∈ F_p^n : f(a) = 0 ∀ f ∈ I }`（A38，因 `F_p^n` 有限故可完全枚舉）。
**定義 2.7（根）** `rad(I) = { f : f^m ∈ I }`；零維時用 Seidenberg 引理（A39）：
`rad(I) = I + (sf(minpoly_v(x_v)))_v`，其中 `minpoly_v` 是 `I ∩ F_p[x_v]` 的單調生成元（A40）。
**定義 2.8（域方程）** `(x_0^p - x_0, …, x_{n-1}^p - x_{n-1})`（A38/A39）。
*算式（有限域 Nullstellensatz，T6）*：`I(V(I)) = rad(I + (x_i^p - x_i))`。
**定義 2.9（極大）** `I` 極大 ⟺ `R/I` 是域 ⟺ `R/I` 無零因子（T8，A41）。
**定義 2.10（乘法矩陣）** `m_f : R/I → R/I` 在基 `basis` 上的矩陣（A42）；
`det(m_f) = 0 ⟺ f` 是零因子；`f` 的極小多項式 = `m_f` 的極小多項式。
**定義 2.11（Rabinowitsch）** `f ∈ rad(I) ⟺ 1 ∈ (I, 1 - y f) ⊂ R[y]`（T7，A43）。

---

## 3. 八大結構：定義與算式

### 3.1 仿射半群環（`semigroup`）

**定義** `S = N A ⊂ N^d`（`A` 的列為生成元），`k[S] = ⊕_{s ∈ S} k·x^s ≅ k[x_1..x_m]/I_A`。
*算式*：**環面理想** `I_A = ( x^u - x^v : A u = A v, u ≠ v )`（A62）。
**定義（格完成）** `L = Z S`，`S` 正規（saturated）⟺ `∀ v ∈ L ∩ cone(S) ∃ k ≥ 1. kv ∈ S ⇒ v ∈ S`（A59）。
*算式（Ehrhart）*：`L_P(k) = |kP ∩ Z^d|` 是 `k` 的 `d` 次多項式，首項 `= vol(P)`，
且 `L_P(-k) = (-1)^d |int(kP) ∩ Z^d|`（Macdonald 互反，A60/A61）。
*算式（Hilbert）*：`H_{k[S]}(k) = |S ∩ {deg = k}|`，且標準單項式 ↔ 半群元素雙射（S4）。

### 3.2 上同調環（`cohom`）

**定義** `H^*(X; F_p) = ⊕_k H^k`，杯積 `∪ : H^i × H^j → H^{i+j}`（A77）。
*算式*：**分次交換律** `a ∪ b = (-1)^{|a||b|} b ∪ a`；結合律；單位 `[X]`。
*算式（Poincaré 對偶）*：配對矩陣 `H^k × H^{n-k} → H^n ≅ F_p` 的秩 = `dim H^k`，
故 `β_k = β_{n-k}`；`χ(X) = Σ_k (-1)^k β_k`。
**模型**：`S^n = F_p[x]/(x²), |x| = n`；`T² = Λ(a,b), |a|=|b|=1`；`CP² = F_p[h]/(h³), |h|=2`；
`RP²(F_2) = F_2[x]/(x³), |x|=1`；`S²×S² = F_p[a,b]/(a²,b²)`。
**定義（de Rham DG 代数）** `Ω^*(J) = J ⊕ J·dt`，`J = F_p[t]/(t^n)`，
`d(f + g dt) = f' dt`，楔積 `(f + g dt) ∧ (f' + g' dt) = (ff' + (fg' + gf')dt)` 且 `dt ∧ dt = 0`（A78）。
*算式*：`d² = 0`；分次 Leibniz `d(α∧β) = dα∧β + (-1)^{|α|} α∧dβ`；
`d/dt` 下降到 `J` **若且唯若** `p | n`（因 `d(t^n) = n t^{n-1}`）；
`dim H⁰ = #{ i < n : p | i }`，`dim H¹ = n - rank(d)`，`χ(Ω) = χ(H) = 0`。

### 3.3 熱帶半環（`tropical`, `heat`）

**定義** `T = Z ∪ {∞}`，`a ⊕ b = min(a,b)`，`a ⊗ b = a + b`；`∞` 是加法單位（零元）。
*義務*：冪等 `a ⊕ a = a`、交換、結合、`⊗` 對 `⊕` 分配 —— 對**所有**三元組窮盡驗證。
**定義（熱帶多項式）** `f = ⨁_a c_a ⊗ x^{⟨a,w⟩}`；`f(w) = min_a (c_a + ⟨a,w⟩)`（A48 求值同態）。
**定義（熱帶超曲面）** `V(f) = { w : f 的最小值至少被兩項達到 }`，重數 = 達到最小值之項數。
*算式（熱帶行列式）*：`trop_det(M) = min_σ Σ_i M[i][σ(i)]`（A48 Hungarian 演算法 = 暴力枚舉，n ≤ 6）。
*算式（Newton 細分）*：`V(f)` 的胞腔 = 下凸包 `hull{(a, c_a)}` 的垂直投影面（A50）；
**平衡條件** `Σ_E len(E)·ν_E = 0`（A51）。
**定義（Maslov 去量子化 / 轉熱）** `φ_h(x) = -h ln x`，
`x ⊕_h y = -h ln(e^{-x/h} + e^{-y/h})`（A46 softmin），`x ⊗_h y = x + y`。
*算式*：`φ_h` 是 semifield 同構；`lim_{h→0} (x ⊕_h y) = min(x,y)`；
界 `-h ln n ≤ softmin_h(E) - min E ≤ 0`；殘餘熵 `S(h→0) = ln g`（`g` = 最小值重數）；
`softmin_h` 隨 `h` 遞減而**單調遞增**。
*算式（Legendre／Gibbs，A47）*：`Z(h) = Σ e^{-E_i/h}`，`F = -h ln Z = ⟨E⟩ - h S`，
`F ≤ min E`，`lim_{h→0} F = min E`，Gibbs 態的 `argmax` = 熱帶 `argmin`。

### 3.4 光滑函數環（`smooth`）

**定義（jet 代數）** `J_k = F_p[t]/(t^{k+1})`：光滑函數芽的 `k` 階截斷（合成律 = 代入）。
*算式*：`t^{k+1} = 0`；`f` 可逆 ⟺ `f(0) ≠ 0`；非單位皆冪零；剩余同態 `J_k → F_p`。
**定義（Kock–Lawvereux）** `∀ f ∈ J_1 ∃! (a, d). f = a + d·ε`（`ε² = 0`）—— 微分算子的公理化。
*算式（Taylor）*：`f(x + ε) = Σ_n (D^n f)(x) ε^n / n!`，其中 `D^n` 是 **Hasse 除法冪導子**（A69）：
`D^n(Σ a_i t^i) = Σ C(i+n, n) a_{i+n} t^i`（特徵 `p` 下仍有意義）。
*算式（導子）*：Leibniz `d(fg) = f'g + fg'`；`[d/dt, E] = d/dt`（`E = t·d/dt`）；Jacobi 恆等式。
*算式（Faà di Bruno，A68）*：`(f∘g)^{(n)} = Σ_{k} f^{(k)}(g) · B_{n,k}(g', g'', …)`。
*算式（熱流，A70）*：`u(s) = e^{sΔ} u_0 = Σ_m s^m Δ^m u_0/m!`，`∂_s u = Δu`，`u(s+r) = e^{rΔ}u(s)`。
*算式（Lagrange，A71）*：每個映射 `F_p → F_p` 都是次數 `< p` 的多項式，且此多項式唯一。

### 3.5 完美的 F-分裂環（`fsplit`）

**設定** `R = F_p[x_1..x_n]/I`，`dim_{F_p} R < ∞`；`B = (e_0 = 1, e_1, …)` 標準單項式基（A34），
結構常數 `μ[a][j][m]`：`e_a e_j = Σ_m μ[a][j][m] e_m`（A79）。

**定義（Frobenius）** `F(r) = r^p`，是環同態；其矩陣的第 `i` 欄 = `e_i^p` 的座標。
**定義（F-分裂）** `φ : F_*R → R` 是 `R`-線性且 `φ ∘ F = id_R`。以 `r·s := r^p s` 表示 `F_*R` 的
`R`-作用，則 `R`-線性等價於
```
   (i)  φ(r^p)   = r                ∀ r ∈ R      （分裂條件，含 φ(1) = 1）
   (ii) φ(r^p s) = r · φ(s)         ∀ r, s ∈ R   （R-線性）
```
`R` 稱為 **F-split** 若此系統有解。
**定義（完美）** `F` 為雙射（有限維時 ⟺ 滿射 ⟺ `rank(F) = dim R`）。
**定義（相容分裂理想）** `J ≤ R` 被 `φ` 相容分裂 ⟺ `φ(J) ⊆ J`。

*算式（X1，可判定性）*：把 `φ` 的矩陣元素 `M[i][j]` 當未知數，(i)(ii) 是 `F_p` 上的
**仿射線性系統**（`dim²` 個未知數，`dim² + dim³` 個方程），以 A27 高斯消去判定相容性：
```
   (i)  Σ_k (e_r^p)_k · M[i][k] = δ_{i,r}
   (ii) Σ_k (e_a^p e_b)_k · M[i][k] − Σ_j μ[a][j][i] · M[j][b] = 0
```
*算式（X3，torsor）*：分裂集合 = `φ₀ + V`，`V = { ψ : ψ R-線性, ψ∘F = 0 }`，`|V| = p^{nullity}`。
*算式（X2）*：`F-split ⟺ reduced ⟺ perfect ⟺ rank(F) = dim R`（有限維 `F_p`-代數）。
*算式（X4）*：`R = F_p^n` 時 `φ = id` 唯一（nullity 0），理想恰有 `2^n` 個，皆由冪等元生成。
*算式（X6）*：`R = F_{p²}` 時 `F² = id`，故唯一分裂 `φ = F = F^{-1}`。

### 3.6 分次環（`graded`）

**定義** `A = ⊕_{d ≥ 0} A_d`，`A_i · A_j ⊆ A_{i+j}`；以結構常數 `mult[i][j][k]` 表示
`e_i · e_j = Σ_k mult[i][j][k] e_k`（`deg e_i` 已知）。
*義務*：次數保持 `mult[i][j][k] ≠ 0 ⇒ deg_k = deg_i + deg_j`；結合律
`Σ_k mult[i][j][k] mult[k][t][m] = Σ_k mult[j][t][k] mult[i][k][m]`；
分次交換律 `mult[i][j][k] = (-1)^{deg_i deg_j} mult[j][i][k]`；雙線性。
*算式（Hilbert 級數／Hilbert–Serre，A52/G4）*：`H_A(t) = Σ_d H(d) t^d = Q(t)/(1-t)^{kdim}`，
`H(d)` 對 `d ≫ 0` 是次數 `kdim - 1` 的多項式（以有限差分偵測，並與三種獨立算法互驗）。
*算式（Veronese，A53/G5）*：`H_{A^{(k)}}(d) = H_A(kd)`。
*算式（adic 分次，A55/G6）*：`Σ_k dim(J^k/J^{k+1}) = dim R`。
*算式（外代数，A56/G7）*：`Λ(V)` 上 `e_i ∧ e_j = -e_j ∧ e_i`，`e_i ∧ e_i = 0`，`dim Λ^d = C(n,d)`。

### 3.7 p 進整數環（`padic`）

**定義** `Z_p = lim← Z/p^k`，以 `prec` 位數字截斷表示（A63）；截斷是同態（P1）。
**定義（賦值）** `v_p(x) = max{ k : p^k | x }`，`v_p(0) = ∞`（A10）。
*算式*：`v(xy) = v(x) + v(y)`；**超度量不等式** `v(x+y) ≥ min(v x, v y)`，
且 `v(x) ≠ v(y) ⇒ v(x+y) = min`（等號情形）；`|x|_p = p^{-v(x)}`（Ostrowski）。
**定義（DVR）** `Z_p` 是離散賦值環：每個理想形如 `(p^{min v})`（P4）。
*算式（Teichmüller，A64）*：`ω(a) = lim_n a^{p^n}`，`ω(a)^{p-1} = 1`，`ω(a) ≡ a (mod p)`，
`ω(ab) = ω(a)ω(b)`；單位分裂 `Z_p^* ≅ μ_{p-1} × (1 + pZ_p)`（P5）。
*算式（Hensel，A65）*：`f(r_0) ≡ 0, f'(r_0) ≢ 0 (mod p) ⇒ ∃! r ∈ Z_p. f(r) = 0 ∧ r ≡ r_0`，
Newton 迭代 `x ← x - f(x)/f'(x)`（P6，含窮盡唯一性驗證）。
*算式（exp/log，A66）*：`exp(x) = Σ x^n/n!`（`v(x) ≥ 1, p ≥ 3` 收斂），`log(1+u) = Σ (-1)^{n+1} u^n/n`，
互逆且為同態（P7）。
*算式（Newton 多邊形，A67）*：`{(i, v(a_i))}` 的下凸包斜率 `= -`(根的賦值)；
`Σ dx = deg f`，`Σ dy = v(a_n) - v(a_0)`；例：`x² - p` 的斜率 `-1/2`（P8）。

### 3.8 理想類群對應環（`classgroup`）

**定義（判別式）** `D < 0`，`D ≡ 0,1 (mod 4)`；二次序 `O_D = Z[ω]`，
```
   ω = (D + √D)/2,     ω² = D·ω − D(D−1)/4,
   (x + yω)(x' + y'ω) = (xx' − yy'·D(D−1)/4) + (xy' + yx' + yy'D)·ω      （A83）
```
**定義（二元二次型）** `f = (a,b,c) = a x² + b xy + c y²`，`disc(f) = b² − 4ac`；
primitive ⟺ `gcd(a,b,c) = 1`；reduced ⟺ `|b| ≤ a ≤ c` 且（`|b| = a` 或 `a = c` 時 `b ≥ 0`）。
*算式（Gauss 約化，A82）*：`b ← b + 2ak`（同時 `c ← (b²−D)/4a`）與 `(a,b,c) ← (c,−b,a)`
（皆屬 `SL_2(Z)`）反覆施行至 reduced；`D < 0` 時 reduced 代表元**唯一**（K1）。
**定義（對應）** `(a,b,c) ↔ I(f) = Z·a + Z·((-b+√D)/2)`，在 `ω`-座標下基為
`(a, 0), ((-b-D)/2 mod a, 1)`；`N(I) = |det(基矩陣)|`（K5）。
**定義（合成）** `I(f·g) = I(f)·I(g)`：四個基向量乘積的 `Z`-span，經 HNF（A29）回到基（A84）。
*算式（由基讀回型）*：對基 `(α, β)`，
```
   f(x,y) = N(αx + βy)/N(I),   a = N(α)/N(I),   b = −Tr(α·β̄)/N(I),   c = N(β)/N(I)
```
（改變 `SL_2(Z)` 基得到等價型，故 HNF 給出典型類，再約化得代表元。）
*算式（類數，K3）*：`h(D) = #{reduced primitive forms of discriminant D}`，
`a ≤ √(|D|/3)`，`b² ≡ D (mod 4a)`。
*算式（Kronecker／genus，A85/K6）*：`(d/n)` 由 Jacobi 符號（A15）與 `(d/2) = (-1)^{(d²-1)/8}` 組合；
質數 `q ∤ D`：`q` ramified ⟺ `q | D`，split ⟺ `(D/q) = 1`，inert ⟺ `(D/q) = -1`；
genus 數 `= 2^{t-1}`（`t` = `D` 的不同質因數個數），discriminant divisor 有 `2^t − 1` 個，
所有 genus 特徵之積為平凡特徵，各 genus 大小相等。
**定義（對應環）** `F_p[Cl(D)]`：以類群表 `table[i][j]` 為分次結構的**群代數**（A86），
是 `Cl(D)`-分次的交換結合含幺環；augmentation `ε(Σ c_g[g]) = Σ c_g` 是環同態（K7）。
*算式（Maschke）*：`F_p[Cl(D)]` reduced（半單）⟺ `p ∤ h(D)`；`p | h(D)` 時 Loewy 長 `= v_p(h)+1`。

---

## 4. 終極運算力閉環鏈條（`loop`）

**定義（閉環）** 對每個結構 `X`，鏈條為

```
   X  ──①──▶  R = F_p[x]/I（有限維）  ──②──▶  熱帶化  ──③──▶  同調  ──④──▶  不變量
                                                                              │
                    ⑤ 約化 R ↦ R/rad(0)（閉包算子，冪等）◀──────────────────────┘
```

- **① 有限域資料**：`dim R`、Hilbert 函數 `H(d)`（`Σ_d H(d) = dim R`）、Krull 維度、
  reduced（`I = rad I`）、F-split、Loewy 長度。
- **② 熱帶化（轉熱）**：`J = rad(I)/I` 之 **J-adic 賦值**（A87）
  `v(x) = max{ k : x ∈ J^k }`（`v(0) = ∞`），是取值於熱帶半環 `(N ∪ {∞}, min, +)` 的 Krull 賦值：
  ```
     v(x + y) ≥ min(v x, v y),      v(xy) ≥ v(x) + v(y),      v(x) = ∞ ⟺ x = 0
  ```
  再以 **Maslov 橋**（A88）驗證 `softmin_h(v(·)) → min v`（`h → 0`）。
- **③ 同調**：Koszul 複形 `0 → R --·a→ R → 0`（A75），`H_1 = ann(a)`、`H_0 = R/aR`，
  `χ(C) = χ(H)`，並以 rank–nullity 交叉驗證 `dim H_i = dim R − rank(a·)`。
- **④ 不變量**：`inv = 100·Loewy + dim`（可比較、可重現）。
- **⑤ 閉包**：`R ↦ R/rad(0)` 冪等；`reduced ⟺ F-split ⟺ Loewy = 1`，
  故約化後的環是鏈條的**不動點**（L1/L4）。

**八結構的具體實例**（`run_all(p)`）：

| # | 結構 | 環 `R` | `p=2` 的 `dim / Loewy / reduced` |
|---|---|---|---|
| 1 | 仿射半群環 | `F_p[x,y,z]/(I_A + (x²,y²,z²))`, `A = ⟨(2,0),(0,2),(1,1)⟩` | 6 / 3 / false |
| 2 | 上同調環 | `H^*(CP²) = F_p[h]/(h³)` | 3 / 3 / false |
| 3 | 熱帶半環 | `F_p[x]/(x^p − x) = F_p^p`（冪等元豐富 = `(min,+)` 的影子） | 2 / 1 / true |
| 4 | 光滑函數環 | jet 代數 `F_p[t]/(t^p)` | 2 / 2 / false |
| 5 | 完美的 F-分裂環 | `F_{p²}` | 2 / 1 / true |
| 6 | 分次環 | `F_p[x,y]/(x²,y²)` | 4 / 3 / false |
| 7 | p 進整數環 | `Z/p^k`（賦值、Teichmüller、exp/log、Hensel、`0 → Z/p^k --p→ Z/p^k → 0` 之 `H_0 = F_p`） | k / k / false |
| 8 | 理想類群對應環 | `F_p[Cl(−23)] ≅ F_p[x]/(x³−1)` | 3 / 1 / true |

跨層閉環驗證（L1）：`②` 的 Hilbert 函數經**次數加倍重分次**後 `=` `H^*(CP²)` 的 Betti 數；
`④` 的 `χ` 由三條獨立路徑（複形秩、同調秩、rank–nullity）一致；
`⑤` 的約化對所有結構冪等 ⇒ **`closed = true`**。
