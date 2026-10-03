# 演算法目錄（共 260 條）

本庫以 MoonBit 實作，零第三方依賴。每條演算法於源碼中以編號標記；
其中 82 條帶有可執行義務自證（`registry.mbt`），另有 152 項定理級白盒測試。


## 向量空間衍生（vecspace.mbt）

- **A-1** — 減法 u−v := u+(−v)  
  <sub>vecspace.mbt</sub>
- **A-2** — 線性組合 Σ cᵢvᵢ  
  <sub>vecspace.mbt</sub>
- **A-3** — axpy: a·x + y（BLAS 基本算式）  
  <sub>vecspace.mbt</sub>
- **A-4** — 向量列求和  
  <sub>vecspace.mbt</sub>
- **A-5** — 中點（仿射組合 ½a + ½b）  
  <sub>vecspace.mbt</sub>
- **A-6** — 倍增迭代：2ⁿ·v（僅用加法實作純量乘法，驗證分配律用）  
  <sub>vecspace.mbt</sub>
- **A-7** — 仿射組合 Σ cᵢvᵢ（Σcᵢ=1 時結果與原點無關之性質由測試驗證）  
  <sub>vecspace.mbt</sub>
- **A-8** — 純量乘法之逆：(1/c)·v  
  <sub>vecspace.mbt</sub>

## 內積空間衍生（ipspace.mbt）

- **I-1** — 範數 ‖v‖ = √⟨v,v⟩  
  <sub>ipspace.mbt</sub>
- **I-2** — 距離 d(u,v) = ‖u−v‖  
  <sub>ipspace.mbt</sub>
- **I-3** — 餘弦角 cos θ = ⟨u,v⟩/(‖u‖‖v‖)（零向量時回 1）  
  <sub>ipspace.mbt</sub>
- **I-4** — 角度 θ = arccos(cos θ) ∈ [0, π]  
  <sub>ipspace.mbt</sub>
- **I-5** — 正交判定 |⟨u,v⟩| < eps  
  <sub>ipspace.mbt</sub>
- **I-6** — 單位化 v/‖v‖（零向量回零）  
  <sub>ipspace.mbt</sub>
- **I-7** — 投影係數 ⟨u,v⟩/⟨v,v⟩  
  <sub>ipspace.mbt</sub>
- **I-8** — 正交投影 proj_v(u) = ⟨u,v⟩/⟨v,v⟩ · v  
  <sub>ipspace.mbt</sub>
- **I-9** — 正交分量（拒絕）rej_v(u) = u − proj_v(u)  
  <sub>ipspace.mbt</sub>
- **I-10** — Gram–Schmidt 正交化（丟棄線性相關者，容差 eps）  
  <sub>ipspace.mbt</sub>
- **I-11** — Gram–Schmidt 標準正交化  
  <sub>ipspace.mbt</sub>
- **I-12** — 標準正交基下的座標 ⟨v, eᵢ⟩  
  <sub>ipspace.mbt</sub>
- **I-13** — 由座標重構向量 Σ cᵢeᵢ（Fourier 重構）  
  <sub>ipspace.mbt</sub>
- **I-14** — 子空間正交投影（輸入為子空間之標準正交基）  
  <sub>ipspace.mbt</sub>
- **I-15** — 到子空間之距離 ‖v − proj v‖  
  <sub>ipspace.mbt</sub>
- **I-16** — 鏡射（Householder 型）：沿單位法向量 n 反射 v ↦ v − 2⟨v,n⟩n  
  <sub>ipspace.mbt</sub>
- **I-17** — Gram 矩陣 Gᵢⱼ = ⟨vᵢ, vⱼ⟩（傳回扁平方陣，列主序）  
  <sub>ipspace.mbt</sub>
- **I-18** — 平行多面體體積 = √det(Gram)（依賴 matrix.mbt 之 det）  
  <sub>ipspace.mbt</sub>
- **I-19** — Bessel 殘差 ‖v‖² − Σ⟨v,eᵢ⟩²（對標準正交系必 ≥ 0）  
  <sub>ipspace.mbt</sub>
- **I-20** — 加權內積尺度變換：⟨,⟩' = c·⟨,⟩  
  <sub>ipspace.mbt</sub>

## 標準歐幾里得空間 R^n（euclid.mbt）

- **E-1** — 零向量  
  <sub>euclid.mbt</sub>
- **E-2** — 逐點加法  
  <sub>euclid.mbt</sub>
- **E-3** — 負元  
  <sub>euclid.mbt</sub>
- **E-4** — 純量乘法  
  <sub>euclid.mbt</sub>
- **E-5** — 標準內積（點積）Σ aᵢbᵢ  
  <sub>euclid.mbt</sub>
- **E-6** — R^n 向量空間字典  
  <sub>euclid.mbt</sub>
- **E-7** — R^n 內積空間字典  
  <sub>euclid.mbt</sub>
- **E-8** — 標準基 {e_0, ..., e_{n-1}}  
  <sub>euclid.mbt</sub>
- **E-9** — 第 i 個標準基向量  
  <sub>euclid.mbt</sub>
- **E-10** — 三維叉積  
  <sub>euclid.mbt</sub>
- **E-11** — Hadamard 逐點積  
  <sub>euclid.mbt</sub>
- **E-12** — 分量總和（線性泛函）  
  <sub>euclid.mbt</sub>
- **E-13** — 均值  
  <sub>euclid.mbt</sub>
- **E-14** — 近似相等  
  <sub>euclid.mbt</sub>
- **E-15** — 最大值分量  
  <sub>euclid.mbt</sub>
- **E-16** — 最大絕對值分量（∞-範數）  
  <sub>euclid.mbt</sub>

## 矩陣空間 M_{m×n}（matrix*.mbt）

- **M-1** — 零矩陣  
  <sub>matrix.mbt</sub>
- **M-2** — 由列向量陣列構造  
  <sub>matrix.mbt</sub>
- **M-3** — 取元素  
  <sub>matrix.mbt</sub>
- **M-4** — 設元素  
  <sub>matrix.mbt</sub>
- **M-5** — 第 i 列  
  <sub>matrix.mbt</sub>
- **M-6** — 第 j 行  
  <sub>matrix.mbt</sub>
- **M-7** — 加法  
  <sub>matrix.mbt</sub>
- **M-8** — 減法  
  <sub>matrix.mbt</sub>
- **M-9** — 純量乘法  
  <sub>matrix.mbt</sub>
- **M-10** — 矩陣乘法 (AB)ᵢₖ = Σⱼ AᵢⱼBⱼₖ  
  <sub>matrix.mbt</sub>
- **M-11** — 矩陣-向量乘  
  <sub>matrix.mbt</sub>
- **M-12** — 轉置  
  <sub>matrix.mbt</sub>
- **M-13** — 單位矩陣  
  <sub>matrix.mbt</sub>
- **M-14** — 跡 tr(A) = Σ Aᵢᵢ  
  <sub>matrix.mbt</sub>
- **M-15** — Frobenius 內積 ⟨A,B⟩_F = Σ AᵢⱼBᵢⱼ = tr(AᵀB)  
  <sub>matrix.mbt</sub>
- **M-16** — Frobenius 範數  
  <sub>matrix.mbt</sub>
- **M-17** — 近似相等  
  <sub>matrix.mbt</sub>
- **M-18** — 扁平方陣行列式（部分主元高斯消去）——供內積空間體積公式使用  
  <sub>matrix.mbt</sub>
- **M-19** — 行列式（部分主元高斯消去，O(n³)）  
  <sub>matrix.mbt</sub>
- **M-20** — 約化列梯形（RREF，部分主元，容差 eps）  
  <sub>matrix.mbt</sub>
- **M-21** — 秩 = RREF 中主元數  
  <sub>matrix2.mbt</sub>
- **M-22** — 主元行位置（RREF 中每列首個非零元之行號）  
  <sub>matrix2.mbt</sub>
- **M-23** — 逆矩陣（Gauss–Jordan；奇異時回 None）  
  <sub>matrix2.mbt</sub>
- **M-24** — 解 Ax = b（唯一解；奇異時回 None）  
  <sub>matrix2.mbt</sub>
- **M-25** — LU 分解（部分主元）：PA = LU，回傳 (L, U, perm)  
  <sub>matrix2.mbt</sub>
- **M-26** — 由 LU 解方程：先解 Ly = Pb，再解 Ux = y  
  <sub>matrix2.mbt</sub>
- **M-27** — Cholesky 分解 A = LLᵀ（A 須對稱正定；失敗回 None）  
  <sub>matrix2.mbt</sub>
- **M-28** — 冪迭代：主特徵值與特徵向量（Rayleigh 商）  
  <sub>matrix2.mbt</sub>
- **M-29** — Jacobi 旋轉法求對稱矩陣特徵值  
  <sub>matrix2.mbt</sub>
- **M-30** — QR 分解（Gram–Schmidt 於行上）：A = QR，Q 行標準正交  
  <sub>matrix3.mbt</sub>
- **M-31** — 列空間基（A 的主元行）  
  <sub>matrix3.mbt</sub>
- **M-32** — 行空間基（RREF 非零列）  
  <sub>matrix3.mbt</sub>
- **M-33** — 零空間基（自由變量法）  
  <sub>matrix3.mbt</sub>
- **M-34** — 最小二乘解（法方程 AᵀAx = Aᵀb，經 Cholesky 或 LU）  
  <sub>matrix3.mbt</sub>
- **M-35** — 外積 u vᵀ  
  <sub>matrix3.mbt</sub>
- **M-36** — Kronecker 積 A⊗B  
  <sub>matrix3.mbt</sub>
- **M-37** — 對稱化：sym(A) = (A + Aᵀ)/2（到對稱矩陣空間之正交投影）  
  <sub>matrix3.mbt</sub>
- **M-38** — 反對稱化：skew(A) = (A − Aᵀ)/2  
  <sub>matrix3.mbt</sub>
- **M-39** — 對稱判定  
  <sub>matrix3.mbt</sub>
- **M-40** — 矩陣指數 e^A（scaling-and-squaring + Taylor 20 項）  
  <sub>matrix3.mbt</sub>
- **M-41** — 譜分解重構：由 (特徵值, Q) 重建 A ≈ Q diag(λ) Qᵀ  
  <sub>matrix3.mbt</sub>
- **M-42** — Rayleigh 商 xᵀAx / xᵀx  
  <sub>matrix3.mbt</sub>
- **M-43** — 對角矩陣構造  
  <sub>matrix3.mbt</sub>
- **M-44** — 取對角線  
  <sub>matrix3.mbt</sub>

## 對稱矩陣空間 Sym(n)（symm.mbt）

- **S-1** — 對稱空間維數 n(n+1)/2  
  <sub>symm.mbt</sub>
- **S-2** — 上三角座標 (i≤j) → 線性索引  
  <sub>symm.mbt</sub>
- **S-3** — 向量打包為對稱矩陣（維數 n(n+1)/2 之座標）  
  <sub>symm.mbt</sub>
- **S-4** — 對稱矩陣展開為向量（pack 之逆）  
  <sub>symm.mbt</sub>
- **S-5** — Sym(n) 加法  
  <sub>symm.mbt</sub>
- **S-6** — Sym(n) 純量乘法  
  <sub>symm.mbt</sub>
- **S-7** — Sym(n) 內積（Frobenius 限制）  
  <sub>symm.mbt</sub>
- **S-8** — Sym(n) 向量空間字典  
  <sub>symm.mbt</sub>
- **S-9** — Sym(n) 內積空間字典  
  <sub>symm.mbt</sub>
- **S-10** — Sym(n) 標準正交基：{E_ii} ∪ {(E_ij+E_ji)/√2}  
  <sub>symm.mbt</sub>
- **S-11** — 譜定理驗證用：對稱矩陣特徵分解重構誤差  
  <sub>symm.mbt</sub>
- **S-12** — 二次型 q(x) = xᵀAx  
  <sub>symm.mbt</sub>

## 對角矩陣空間 / 對偶空間（diag.mbt / dual.mbt）

- **D-1** — 對角空間維數  
  <sub>diag.mbt</sub>
- **D-2** — 對角加法（逐點）  
  <sub>diag.mbt</sub>
- **D-3** — 對角純量乘法  
  <sub>diag.mbt</sub>
- **D-4** — 對角矩陣乘法 ≙ 逐點積  
  <sub>diag.mbt</sub>
- **D-5** — 對角內積 ≙ Σ aᵢbᵢ（Frobenius 內積的限制）  
  <sub>diag.mbt</sub>
- **D-6** — 作用於向量：diag(a)·x  
  <sub>diag.mbt</sub>
- **D-7** — 行列式 = 對角線之積  
  <sub>diag.mbt</sub>
- **D-8** — 逆（無零分量時）  
  <sub>diag.mbt</sub>
- **D-9** — 轉為完整矩陣（與 Diag(n) → M_n 之嵌入）  
  <sub>diag.mbt</sub>
- **D-10** — Diag(n) 內積空間字典  
  <sub>diag.mbt</sub>
- **D-11** — 譜半徑（對角矩陣之算子範數）  
  <sub>diag.mbt</sub>
- **D-12** — 跡 = 分量和  
  <sub>diag.mbt</sub>
- **D-1** — 零泛函  
  <sub>dual.mbt</sub>
- **D-2** — 泛函作用 φ(v)  
  <sub>dual.mbt</sub>
- **D-3** — 對偶加法  
  <sub>dual.mbt</sub>
- **D-4** — 對偶純量乘法  
  <sub>dual.mbt</sub>
- **D-5** — 對偶內積（標準歐氏度量下 V ≅ V* 之誘導內積）  
  <sub>dual.mbt</sub>
- **D-6** — 對偶空間字典（n 維）  
  <sub>dual.mbt</sub>
- **D-7** — 對偶基（對標準基即 δᵢ = eᵢ*）  
  <sub>dual.mbt</sub>
- **D-8** — Riesz 表示定理（歐氏情形）：泛函 φ ↦ 表示向量 r，  
  <sub>dual.mbt</sub>
- **D-9** — Riesz 逆：向量 ↦ 泛函 ⟨v, ·⟩  
  <sub>dual.mbt</sub>
- **D-10** — 對偶映射（拉回）：線性映射 A : V → W 誘導  
  <sub>dual.mbt</sub>
- **D-11** — 對偶配對矩陣：Pᵢⱼ = φᵢ(eⱼ)（對偶基與基之配對）  
  <sub>dual.mbt</sub>
- **D-12** — 雙重對偶嵌入 J : V → V**（有限維時為同構，自證 J = id）  
  <sub>dual.mbt</sub>
- **D-13** — 核之判定：φ ∈ (span S)° ⟺ ∀s∈S: φ(s)=0（回傳是否屬零化子）  
  <sub>dual.mbt</sub>

## 複空間 C^n 與 Hermitian（hermitian.mbt）

- **H-1** — C^n 零向量  
  <sub>hermitian.mbt</sub>
- **H-2** — C^n 加法  
  <sub>hermitian.mbt</sub>
- **H-3** — C^n 負元  
  <sub>hermitian.mbt</sub>
- **H-4** — C^n 複純量乘法  
  <sub>hermitian.mbt</sub>
- **H-5** — C^n 實純量乘法（視為實向量空間）  
  <sub>hermitian.mbt</sub>
- **H-6** — Hermitian 內積 ⟨u,v⟩ = Σ conj(uᵢ)vᵢ  
  <sub>hermitian.mbt</sub>
- **H-7** — 實值二次範數 ⟨v,v⟩ ∈ R  
  <sub>hermitian.mbt</sub>
- **H-8** — 範數  
  <sub>hermitian.mbt</sub>
- **H-9** — C^n 之內積空間字典（視為實內積空間：取內積實部）  
  <sub>hermitian.mbt</sub>
- **H-10** — 正交性（複）  
  <sub>hermitian.mbt</sub>
- **H-11** — 複 Gram–Schmidt 標準正交化  
  <sub>hermitian.mbt</sub>
- **H-12** — 複零矩陣  
  <sub>hermitian.mbt</sub>
- **H-13** — 取元素  
  <sub>hermitian.mbt</sub>
- **H-14** — 設元素  
  <sub>hermitian.mbt</sub>
- **H-15** — 共軛轉置（dagger）  
  <sub>hermitian.mbt</sub>
- **H-16** — Hermitian 判定 A = A†  
  <sub>hermitian.mbt</sub>
- **H-17** — 複矩陣-向量乘  
  <sub>hermitian.mbt</sub>
- **H-18** — 複矩陣乘法  
  <sub>hermitian.mbt</sub>
- **H-19** — 複 Rayleigh 商（Hermitian 時為實數）  
  <sub>hermitian.mbt</sub>
- **H-20** — Hermitian 矩陣特徵值（複 Jacobi 旋轉，非對角元消去）  
  <sub>hermitian.mbt</sub>
- **H-21** — 酉性誤差 ‖Q†Q − I‖_F  
  <sub>hermitian.mbt</sub>

## 多項式空間 P_n（poly.mbt）

- **P-1** — 零多項式  
  <sub>poly.mbt</sub>
- **P-2** — 常數多項式  
  <sub>poly.mbt</sub>
- **P-3** — 次數（零多項式定為 0）  
  <sub>poly.mbt</sub>
- **P-4** — 首項係數  
  <sub>poly.mbt</sub>
- **P-5** — 加法  
  <sub>poly.mbt</sub>
- **P-6** — 減法  
  <sub>poly.mbt</sub>
- **P-7** — 純量乘法  
  <sub>poly.mbt</sub>
- **P-8** — 乘法（卷積）  
  <sub>poly.mbt</sub>
- **P-9** — 導數  
  <sub>poly.mbt</sub>
- **P-10** — 不定積分（常數項為 0）  
  <sub>poly.mbt</sub>
- **P-11** — Horner 求值  
  <sub>poly.mbt</sub>
- **P-12** — L²[0,1] 內積（精確）：Σ aᵢbⱼ/(i+j+1)  
  <sub>poly.mbt</sub>
- **P-13** — P_n 內積空間字典（次數 ≤ n）  
  <sub>poly.mbt</sub>
- **P-14** — 移位 Legendre 多項式：Gram–Schmidt 於 {1,x,…,x^n}  
  <sub>poly.mbt</sub>
- **P-15** — 帶餘除法：p = q·quot + rem（deg rem < deg q）  
  <sub>poly.mbt</sub>
- **P-16** — 餘式  
  <sub>poly.mbt</sub>
- **P-17** — 最大公因式（歐幾里得演算法，首一化）  
  <sub>poly.mbt</sub>
- **P-18** — Lagrange 插值：過點 (xs[i], ys[i])  
  <sub>poly.mbt</sub>
- **P-19** — Newton 插值（均差表）  
  <sub>poly.mbt</sub>
- **P-20** — 複合 p∘q  
  <sub>poly.mbt</sub>
- **P-21** — 由根構造 Π(x − rᵢ)  
  <sub>poly.mbt</sub>
- **P-22** — 二次方程實/複根（回傳 (重數意義下) 兩複數根）  
  <sub>poly.mbt</sub>
- **P-23** — Durand–Kerner 求全部複根（首一化後迭代）  
  <sub>poly.mbt</sub>
- **P-24** — 多項式空間之維數（次數 ≤ n 者為 n+1）  
  <sub>poly.mbt</sub>

## 多變量多項式空間（mpoly.mbt）

- **MP-1** — 正規化：排序、合併同類項、去零  
  <sub>mpoly.mbt</sub>
- **MP-2** — 零多項式  
  <sub>mpoly.mbt</sub>
- **MP-3** — 單項式 c·x^e  
  <sub>mpoly.mbt</sub>
- **MP-4** — 加法  
  <sub>mpoly.mbt</sub>
- **MP-5** — 純量乘法  
  <sub>mpoly.mbt</sub>
- **MP-6** — 乘法  
  <sub>mpoly.mbt</sub>
- **MP-7** — 單項式求值  
  <sub>mpoly.mbt</sub>
- **MP-8** — 多項式求值  
  <sub>mpoly.mbt</sub>
- **MP-9** — 偏導 ∂/∂x_i  
  <sub>mpoly.mbt</sub>
- **MP-10** — 梯度（全部偏導）  
  <sub>mpoly.mbt</sub>
- **MP-11** — 總次數  
  <sub>mpoly.mbt</sub>
- **MP-12** — 係數內積（單項式基下之 ℓ² 內積）  
  <sub>mpoly.mbt</sub>
- **MP-13** — 單項式計數：n 變量次數 ≤ d 之單項式數 = C(n+d, d)  
  <sub>mpoly.mbt</sub>
- **MP-14** — 多項式空間之內積字典（限制於次數 ≤ d 之子空間）  
  <sub>mpoly.mbt</sub>
- **MP-15** — 齊次化：引入新變量 x_n 使總次數齊次  
  <sub>mpoly.mbt</sub>

## 有限集合上函數空間（fspace.mbt）

- **F-1** — 加權 ℓ² 內積 Σ wᵢ f(i) g(i)  
  <sub>fspace.mbt</sub>
- **F-2** — 標準（等權）內積  
  <sub>fspace.mbt</sub>
- **F-3** — Fun(X, R) 內積空間字典（加權）  
  <sub>fspace.mbt</sub>
- **F-4** — delta 基 {δ_0,…,δ_{n−1}}  
  <sub>fspace.mbt</sub>
- **F-5** — 逐點乘法（代數結構）  
  <sub>fspace.mbt</sub>
- **F-6** — 總和泛函  
  <sub>fspace.mbt</sub>
- **F-7** — 循環前向差分 (Δf)(i) = f(i+1) − f(i)（指標模 n）  
  <sub>fspace.mbt</sub>
- **F-8** — 離散傅里葉變換：F[k] = Σⱼ f[j] e^{−2πi·jk/n}  
  <sub>fspace.mbt</sub>
- **F-9** — 逆 DFT：f[j] = (1/n) Σₖ F[k] e^{2πi·jk/n}  
  <sub>fspace.mbt</sub>
- **F-10** — 循環卷積 (f ∗ g)(k) = Σⱼ f[j] g(k−j mod n)  
  <sub>fspace.mbt</sub>
- **F-11** — 函數空間組合子：Fun(X, V)，V 為任意內積空間  
  <sub>fspace.mbt</sub>
- **F-12** — Parseval 殘差 ‖f‖² − (1/n)‖F‖²（對 DFT 必為 0）  
  <sub>fspace.mbt</sub>

## 張量空間 R^m ⊗ R^n（tensor.mbt）

- **T-1** — 零張量  
  <sub>tensor.mbt</sub>
- **T-2** — 純張量 u ⊗ v  
  <sub>tensor.mbt</sub>
- **T-3** — 加法  
  <sub>tensor.mbt</sub>
- **T-4** — 純量乘法  
  <sub>tensor.mbt</sub>
- **T-5** — 誘導內積 ⟨Σ aᵢⱼ eᵢ⊗fⱼ, Σ bᵢⱼ eᵢ⊗fⱼ⟩ = Σ aᵢⱼbᵢⱼ  
  <sub>tensor.mbt</sub>
- **T-6** — 內積空間字典（R^m ⊗ R^n）  
  <sub>tensor.mbt</sub>
- **T-7** — 與矩陣空間之同構（平坦化）  
  <sub>tensor.mbt</sub>
- **T-8** — 矩陣 → 張量（同構之逆）  
  <sub>tensor.mbt</sub>
- **T-9** — 全收縮（m = n 時）：Σᵢ tᵢᵢ ∈ R  
  <sub>tensor.mbt</sub>
- **T-10** — 單腿收縮：以 φ ∈ (R^m)* 收縮第一腿，得 R^n 之向量  
  <sub>tensor.mbt</sub>
- **T-11** — 單腿收縮：以 ψ ∈ (R^n)* 收縮第二腿，得 R^m 之向量  
  <sub>tensor.mbt</sub>
- **T-12** — 雙線性求值：(φ ⊗ ψ)(T) = Σ tᵢⱼ φᵢ ψⱼ  
  <sub>tensor.mbt</sub>
- **T-13** — 線性映射之張量積 A ⊗ B（Kronecker）  
  <sub>tensor.mbt</sub>
- **T-14** — 張量秩上界估计（平坦化矩陣之秩）  
  <sub>tensor.mbt</sub>

## 外代數 Λ^k(R^n)（exterior.mbt）

- **X-1** — 組合數 C(n, k)  
  <sub>exterior.mbt</sub>
- **X-2** — 組合 → 分量索引  
  <sub>exterior.mbt</sub>
- **X-3** — 逆序數（排列奇偶性用）  
  <sub>exterior.mbt</sub>
- **X-4** — Λ^k 之零元素  
  <sub>exterior.mbt</sub>
- **X-5** — 向量 v ∈ R^n 嵌入 Λ^1  
  <sub>exterior.mbt</sub>
- **X-6** — 外積 u ∧ v（兩個 1-向量 → 2-向量）  
  <sub>exterior.mbt</sub>
- **X-7** — 一般外積：α ∈ Λ^p ∧ β ∈ Λ^q → Λ^{p+q}  
  <sub>exterior.mbt</sub>
- **X-8** — Λ^k 之誘導內積（單形分量逐點積；對一般元素用 Gram 行列式性質）  
  <sub>exterior.mbt</sub>
- **X-9** — k-向量範數（= 對應平行體之體積，當為純量時）  
  <sub>exterior.mbt</sub>
- **X-10** — Hodge 星 ⋆ : Λ^k → Λ^{n−k}（標準定向）  
  <sub>exterior.mbt</sub>
- **X-11** — 內縮 i_v : Λ^k → Λ^{k−1}（與向量之內積）  
  <sub>exterior.mbt</sub>
- **X-12** — 頂層外積判定行列式：v_1 ∧ … ∧ v_n = det(V)·e₁∧…∧e_n  
  <sub>exterior.mbt</sub>
- **X-13** — R³ 叉積由外代數導出：u × v = ⋆(u ∧ v)  
  <sub>exterior.mbt</sub>
- **X-14** — Λ^k 之內積空間字典  
  <sub>exterior.mbt</sub>
- **X-15** — 兩向量張成平行四邊形之面積 = ‖u ∧ v‖  
  <sub>exterior.mbt</sub>

## 四元數空間 H 與 H^n（quat.mbt）

- **Q-1** — 構造  
  <sub>quat.mbt</sub>
- **Q-2** — 零  
  <sub>quat.mbt</sub>
- **Q-3** — 純量四元數  
  <sub>quat.mbt</sub>
- **Q-4** — 純虛四元數（對應 R³ 向量）  
  <sub>quat.mbt</sub>
- **Q-5** — 加法  
  <sub>quat.mbt</sub>
- **Q-6** — 負元  
  <sub>quat.mbt</sub>
- **Q-7** — 實純量乘法  
  <sub>quat.mbt</sub>
- **Q-8** — Hamilton 乘法  
  <sub>quat.mbt</sub>
- **Q-9** — 共軛  
  <sub>quat.mbt</sub>
- **Q-10** — 範數平方  
  <sub>quat.mbt</sub>
- **Q-11** — 範數  
  <sub>quat.mbt</sub>
- **Q-12** — 逆 q⁻¹ = q̄/‖q‖²  
  <sub>quat.mbt</sub>
- **Q-13** — 單位化  
  <sub>quat.mbt</sub>
- **Q-14** — 由軸角構造單位四元數（旋轉角 θ 繞單位軸 u）  
  <sub>quat.mbt</sub>
- **Q-15** — 旋轉 R³ 向量：v ↦ q v q̄（q 為單位四元數）  
  <sub>quat.mbt</sub>
- **Q-16** — H 上之內積 ⟨q,p⟩ = Re(q̄p)  
  <sub>quat.mbt</sub>
- **Q-17** — H^n 零向量  
  <sub>quat.mbt</sub>
- **Q-18** — H^n 加法  
  <sub>quat.mbt</sub>
- **Q-19** — H^n 實純量乘法（H^n 視為實向量空間）  
  <sub>quat.mbt</sub>
- **Q-20** — H^n 之實內積 Σ Re(q̄ᵢ pᵢ)  
  <sub>quat.mbt</sub>
- **Q-21** — H^n 作為實內積空間之字典（實維數 4n）  
  <sub>quat.mbt</sub>
- **Q-22** — H^n 之實維數 = 4n  
  <sub>quat.mbt</sub>
- **Q-23** — Gram–Schmidt 於 H^n（經由實內積字典）  
  <sub>quat.mbt</sub>
- **Q-24** — 四元數旋轉複合：q₂∘q₁ 對應四元數乘法（注意次序）  
  <sub>quat.mbt</sub>

## 空間組合子：直和等（product.mbt）

- **PR-1** — 向量空間直和 V ⊕ W  
  <sub>product.mbt</sub>
- **PR-2** — 內積空間直和：⟨(v₁,w₁),(v₂,w₂)⟩ = ⟨v₁,v₂⟩ + ⟨w₁,w₂⟩  
  <sub>product.mbt</sub>
- **PR-3** — 加權直和：⟨,⟩ = α⟨,⟩_V + β⟨,⟩_W  
  <sub>product.mbt</sub>
- **PR-4** — 對角嵌入 Δ : V → V ⊕ V，v ↦ (v, v)  
  <sub>product.mbt</sub>
- **PR-5** — 投影 π₁ : V ⊕ W → V  
  <sub>product.mbt</sub>
- **PR-6** — 投影 π₂ : V ⊕ W → W  
  <sub>product.mbt</sub>
- **PR-7** — 單射 ι₁ : V → V ⊕ W（w 取零）  
  <sub>product.mbt</sub>
- **PR-8** — 單射 ι₂ : W → V ⊕ W  
  <sub>product.mbt</sub>
- **PR-9** — 直和維數相加  
  <sub>product.mbt</sub>
- **PR-10** — 線性映射直和 A ⊕ B 之矩陣（分塊對角）  
  <sub>product.mbt</sub>
