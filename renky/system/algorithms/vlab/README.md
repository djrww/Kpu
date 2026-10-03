# vlab — 向量空間與有限維內積空間之組合式演算法庫（MoonBit）

以 **MoonBit**（v0.1.20260920 工具鏈）實作向量空間與有限維內積空間之公理化框架，
將十二類具體空間以「字典記錄」統一，並提供空間組合子（直和、函數空間），
使任意空間可迭代組合。全庫 **228 條編號演算法**，其中 **82 條帶可執行義務自證**，
另有 **152 項定理級白盒測試**。零第三方依賴（僅用語言內建與官方工具鏈）。

## 覆蓋之具體空間

| 空間 | 檔案 | 字典 |
|---|---|---|
| 標準歐幾里得空間 Rⁿ | `euclid.mbt` | `euc_ip(n)` |
| 矩陣空間 M_{m×n}（Frobenius） | `matrix.mbt` 等 | 經 `mat_frob_inner` |
| 對稱矩陣空間 Sym(n) | `symm.mbt` | `sym_ip(n)` |
| 對角矩陣空間 Diag(n) ≅ Rⁿ | `diag.mbt` | `diag_ip(n)` |
| 複空間 Cⁿ（Hermitian 內積） | `hermitian.mbt` | `cvec_real_ip(n)` |
| Hermitian 矩陣空間 | `hermitian.mbt` | `cmat_jacobi_hermitian` |
| 多項式空間 Pₙ（L²[0,1] 精確內積） | `poly.mbt` | `poly_ip(n)` |
| 多變量多項式空間 R[x₁,…,xₙ] | `mpoly.mbt` | `mpoly_ip(n)` |
| 有限集合上函數空間 Fun(X,F)（加權 ℓ²、DFT） | `fspace.mbt` | `fn_ip(n,w)` |
| 張量空間 R^m ⊗ R^n | `tensor.mbt` | `tensor_ip(m,n)` |
| 外代數 Λ^k(Rⁿ) | `exterior.mbt` | `ext_ip(n,k)` |
| 對偶空間 V* | `dual.mbt` | `dual_ip(n)` |
| 四元數空間 Hⁿ（實內積 ⟨q,p⟩=Re(q̄p)） | `quat.mbt` | `hvec_ip(n)` |

## 空間組合

- `ips_prod(a, b)` / `ips_prod_weighted` — 直和 V ⊕ W（可迭代）
- `fn_space_ip(n, w, ip)` — 函數空間 Fun(X, V)：V 為**任意**內積空間，
  例如 `Fun(X, Sym(2))`、`Fun(X, R²⊕R²)`
- `map_direct_sum` — 線性映射直和；`tensor_of_maps` — 映射張量積

## 衍生演算法摘錄（完整見 `ALGORITHMS.md`）

- **通用內積演算法**：範數、角度、正交投影、鏡射、Gram–Schmidt、
  Fourier 座標重構、Bessel 殘差、Gram 行列式、平行多面體體積
- **矩陣**：行列式（部分主元）、RREF、秩、逆、LU、Cholesky、QR、
  Jacobi 對角化、冪迭代、列/行/零空間基、最小二乘、Kronecker 積、矩陣指數
- **多項式**：Horner、導數/積分、帶餘除法、GCD（歐幾里得）、
  Lagrange/Newton 插值、複合、Durand–Kerner 全複根、移位 Legendre
- **函數空間**：循環差分、DFT/IDFT、循環卷積、Parseval
- **外代數**：外積、Hodge 星、內縮、外積行列式、叉積導出、面積/體積
- **四元數**：Hamilton 乘法、軸角↔四元數、三維旋轉、Hⁿ 實 Gram–Schmidt

## 運行

```bash
moon test          # 152 項定理/公理自證測試
moon run cmd/main  # 執行 82 條演算法義務自證並輸出報告
```

## 義務自證結構

1. **公理層**（`*_wbtest.mbt`）：每個具體實例驗證向量空間八條公理 / 內積三條公理；
2. **定理層**：Cauchy–Schwarz、三角不等式、勾股、Parseval、Cayley–Hamilton、
   譜定理、秩-零度、det 乘法性、Riesz、Bezout、Clairaut、歐拉齊次函數定理等；
3. **註冊表層**（`registry.mbt`）：82 條演算法各帶自證閉包，`registry_run()` 回傳
   (通過數, 總數, 失敗列表)，由 `registry_wbtest` 強制全部通過。

## 突變測試（自研，零依賴）

```bash
python3 tools/mutate.py        # MAX_MUTANTS=400 固定種子抽樣
cat tools/mutation_report.json # 結果
```

突變算子：`+↔−`、`*↔/`、`<↔<=`、`>↔>=`、`&&↔||`、數值字面量邊界。
分類：**UNVIABLE**（moon check 拒絕）／**KILLED**（moon test 失敗）／**SURVIVED**。
擊殺率 = KILLED/(KILLED+SURVIVED)。registry.mbt 屬測試預言（oracle），按慣例排除。
