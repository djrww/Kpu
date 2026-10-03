# ALGORITHMS · 演算法登記冊

本文件由原始碼中的 `**Axx ...**` 標籤自動產生（`grep` 抽取），與程式碼一一對應。
全部演算法皆為本專案自行實作，**無任何第 3 方依賴**（僅使用 MoonBit 隨附之 `moonbitlang/core`）。

登記總數：**89 條**（命案要求 > 18 條 ✓）。

| 編號 | 演算法 | 位置 |
|---|---|---|
| `A01` | modpos — canonical representative of a residue class | `kernel/mod.mbt:2` |
| `A02` | addmod — addition in Z/p | `kernel/mod.mbt:3` |
| `A03` | submod — subtraction in Z/p | `kernel/mod.mbt:4` |
| `A04` | mulmod — multiplication in Z/p | `kernel/mod.mbt:5` |
| `A05` | negmod — additive inverse in Z/p | `kernel/mod.mbt:6` |
| `A06` | argmin — verified arg-minimum | `select/select.mbt:2` |
| `A07` | bsearch — verified binary search | `select/select.mbt:3` |
| `A08` | gcd — Euclid | `arith/div.mbt:23` |
| `A09` | egcd — extended gcd / Bézout | `arith/div.mbt:36` |
| `A10` | valuation v_p(n) | `arith/div.mbt:79` |
| `A11` | powmod — square-and-multiply | `arith/prime.mbt:6` |
| `A12` | invmod — modular inverse | `arith/prime.mbt:25` |
| `A13` | is_prime — Miller–Rabin (deterministic bases) | `arith/prime.mbt:41` |
| `A14` | factorize — trial division + Pollard rho | `arith/prime.mbt:104` |
| `A15` | jacobi — Jacobi symbol | `arith/prime.mbt:170` |
| `A16` | sqrt_mod — Tonelli–Shanks | `arith/prime.mbt:219` |
| `A17` | crt — Chinese remainder theorem | `arith/prime.mbt:268` |
| `A18` | F_p[x] arithmetic | `ff/polyf.mbt:118` |
| `A19` | fp_gcd — Euclid in F_p[x] | `ff/polyf.mbt:261` |
| `A20` | fp_powmod | `ff/polyf.mbt:301` |
| `A21` | fp_dlog — baby-step giant-step | `ff/fpn.mbt:307` |
| `A22` | fp_irreducible — Ben-Or test | `ff/irred.mbt:25` |
| `A23` | fp_find_irreducible — search | `ff/irred.mbt:52` |
| `A24` | fp_ddf — distinct-degree factorisation | `ff/irred.mbt:111` |
| `A25` | fp_edf — Berlekamp equal-degree factorisation | `ff/irred.mbt:137` |
| `A26` | fp_factor — full factorisation | `ff/irred.mbt:220` |
| `A27` | rref — reduced row echelon form over F_p | `linalg/linalg.mbt:117` |
| `A28` | nullspace — kernel basis over F_p | `linalg/linalg.mbt:168` |
| `A29` | hnf — Hermite normal form over Z | `linalg/linalg.mbt:274` |
| `A30` | snf — Smith normal form over Z | `linalg/linalg.mbt:330` |
| `A31` | mp_divmod — multivariate division with remainder | `poly/mpoly.mbt:281` |
| `A32` | mp_spoly — S-polynomial | `poly/mpoly.mbt:324` |
| `A33` | mp_groebner — Buchberger (incremental) | `poly/mpoly.mbt:340` |
| `A34` | mp_standard_monomials — standard monomials of a Gröbner basis | `poly/mpoly.mbt:472` |
| `A35` | mp_krull_dim — Krull dimension from the initial ideal | `ideal/ideal.mbt:246` |
| `A36` | in_ideal / membership_witness — membership test and certificate | `ideal/ideal.mbt:228` |
| `A37` | hilbert_fn — Hilbert function of R/I | `ideal/ideal.mbt:299` |
| `A38` | variety — all F_p-rational points of V(I) | `ideal/ideal.mbt:321` |
| `A39` | radical — Seidenberg radical (+ field equations) | `ideal/ideal.mbt:373` |
| `A40` | univariate_minpoly / fp_squarefree — minimal polynomial and squarefree part | `ideal/ideal.mbt:394` |
| `A41` | is_maximal — maximality test | `ideal/ideal.mbt:462` |
| `A42` | mul_matrix — matrix of multiplication on the quotient basis | `ideal/ideal.mbt:496` |
| `A43` | in_radical — Rabinowitsch trick | `ideal/ideal.mbt:516` |
| `A44` | expd — real exponential (Taylor + range reduction) | `heat/heat.mbt:37` |
| `A45` | logd — real logarithm (atanh series) | `heat/heat.mbt:83` |
| `A46` | softmin_h / softmax_h — Maslov soft min/max | `heat/heat.mbt:112` |
| `A47` | partition / free_energy / gibbs — Gibbs state and Legendre transform | `heat/heat.mbt:181` |
| `A48` | hungarian / trop_det — tropical determinant by assignment | `tropical/tropical.mbt:240` |
| `A49` | trop_mat_mul / trop_closure — min-plus product and Floyd–Warshall | `tropical/tropical.mbt:358` |
| `A50` | newton_cells — Newton subdivision | `tropical/tropical.mbt:407` |
| `A51` | newton_balancing / hull_dim — balancing condition | `tropical/tropical.mbt:701` |
| `A52` | hilbert_series / hilbert_polynomial — Hilbert series and finite-difference polynomial | `graded/graded.mbt:110` |
| `A53` | veronese — Veronese subalgebra | `graded/graded.mbt:152` |
| `A54` | ideal_power — power of an ideal | `graded/graded.mbt:217` |
| `A55` | adic grading / associated graded dimensions | `graded/graded.mbt:248` |
| `A56` | exterior_algebra — Koszul (exterior) algebra with signs | `graded/graded.mbt:275` |
| `A57` | semigroup membership DP + enumeration | `semigroup/semigroup.mbt:71` |
| `A58` | in_lattice — HNF lattice membership | `semigroup/semigroup.mbt:176` |
| `A59` | is_normal — normality / saturation test | `semigroup/semigroup.mbt:220` |
| `A60` | convex_hull / dilated_points / dilated_interior — hull and Ehrhart counting | `semigroup/semigroup.mbt:313` |
| `A61` | newton_coeffs / newton_eval — binomial-basis interpolation | `semigroup/semigroup.mbt:490` |
| `A62` | toric_ideal — toric ideal by exhaustive binomial search | `semigroup/semigroup.mbt:531` |
| `A63` | Padic::inv — Newton inverse | `padic/padic.mbt:198` |
| `A64` | teichmuller — Teichmüller lift | `padic/padic.mbt:243` |
| `A65` | hensel_lift — Hensel/Newton lifting | `padic/padic.mbt:301` |
| `A66` | padic_exp / padic_log — p-adic exponential and logarithm | `padic/padic.mbt:341` |
| `A67` | newton_polygon / newton_slopes — Newton polygon | `padic/padic.mbt:428` |
| `A68` | Jet::compose — substitution and Faà di Bruno | `smooth/smooth.mbt:145` |
| `A69` | Jet::hasse — Hasse divided-power derivatives | `smooth/smooth.mbt:188` |
| `A70` | jet_heat_flow — heat flow on jets | `smooth/smooth.mbt:293` |
| `A71` | lagrange_interpolate — every map F_p → F_p is a polynomial | `smooth/smooth.mbt:327` |
| `A72` | Complex::homology — homology over F_p (rref) and over Z (SNF) | `homology/homology.mbt:136` |
| `A73` | tensor_product — Künneth tensor product of complexes | `homology/homology.mbt:254` |
| `A74` | mapping_cone — mapping cone and the long exact sequence | `homology/homology.mbt:388` |
| `A75` | koszul_two_term — Koszul complex | `homology/homology.mbt:473` |
| `A76` | graph_complex / cycle_rank — simplicial complex of a graph, union-find | `homology/homology.mbt:491` |
| `A77` | CohomRing — graded-commutative cup product and Poincaré pairing | `cohom/cohom.mbt:44` |
| `A78` | DRForm / dr_cohomology_dims — de Rham DG algebra of a jet algebra | `cohom/cohom.mbt:282` |
| `A79` | finite_ring — structure constants and Frobenius matrix | `fsplit/fsplit.mbt:97` |
| `A80` | splitting_system — affine system for a Frobenius splitting | `fsplit/fsplit.mbt:256` |
| `A81` | is_reduced / is_reduced_bruteforce / is_perfect — reduced and perfect tests | `fsplit/fsplit.mbt:221` |
| `A82` | bqf_reduce — Gauss reduction of binary quadratic forms | `classgroup/classgroup.mbt:94` |
| `A83` | omega_mul / ideal_mul / form_to_ideal / ideal_to_form — quadratic-order arithmetic and the form↔ideal correspondence | `classgroup/classgroup.mbt:199` |
| `A84` | bqf_compose — Gauss composition via ideal multiplication | `classgroup/classgroup.mbt:270` |
| `A85` | kronecker / genus_char — Kronecker symbol and genus characters | `classgroup/classgroup.mbt:338` |
| `A86` | class_ring — the group algebra F_p[Cl(D)] | `classgroup/classgroup.mbt:490` |
| `A87` | jadic_filtration / valuation_of — J-adic valuation of an Artin ring | `loop/loop.mbt:216` |
| `A88` | maslov_bridge_ok — soft-min → tropical-min consistency | `loop/loop.mbt:308` |
| `A89` | run_all — the closed-loop pipeline runner | `loop/loop.mbt:662` |
