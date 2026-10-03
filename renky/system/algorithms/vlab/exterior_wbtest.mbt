///|
// 外代數定理之自證

test "exterior: 組合計數與列舉" {
  assert_eq(comb_count(4, 2), 6)
  assert_eq(comb_count(5, 3), 10)
  assert_eq(comb_count(6, 0), 1)
  assert_eq(comb_count(3, 4), 0)
  assert_eq(comb_list(4, 2).length(), 6)
  assert_eq(comb_list(3, 0).length(), 1)
}

test "exterior: 反對稱性與交替性" {
  let u = [1.0, 2.0, 3.0, 0.0]
  let v = [-1.0, 0.0, 2.0, 1.0]
  let a = ext_wedge_vecs(4, u, v)
  let b = ext_wedge_vecs(4, v, u)
  for i = 0; i < a.data.length(); i = i + 1 {
    assert_true((a.data[i] + b.data[i]).abs() < 1.0e-12)
  }
  // u ∧ u = 0
  let z = ext_wedge_vecs(4, u, u)
  for x in z.data {
    assert_true(x.abs() < 1.0e-12)
  }
}

test "exterior: 結合律 (u∧v)∧w = u∧(v∧w)" {
  let n = 4
  let u = ext_from_vec(n, [1.0, 0.0, 1.0, 0.0])
  let v = ext_from_vec(n, [0.0, 1.0, 1.0, 1.0])
  let w = ext_from_vec(n, [1.0, 1.0, 0.0, 2.0])
  let lhs = ext_wedge(ext_wedge(u, v), w)
  let rhs = ext_wedge(u, ext_wedge(v, w))
  for i = 0; i < lhs.data.length(); i = i + 1 {
    assert_true((lhs.data[i] - rhs.data[i]).abs() < 1.0e-12)
  }
}

test "exterior: 分次交換律 α∧β = (−1)^{pq} β∧α" {
  let n = 4
  let a = ext_wedge_vecs(n, [1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0]) // 2-向量
  let b = ext_from_vec(n, [0.0, 0.0, 1.0, 1.0]) // 1-向量
  let ab = ext_wedge(a, b)
  let ba = ext_wedge(b, a)
  // p=2, q=1 → 符號 +1
  for i = 0; i < ab.data.length(); i = i + 1 {
    assert_true((ab.data[i] - ba.data[i]).abs() < 1.0e-12)
  }
  // 兩個 1-向量：符號 −1
  let c = ext_from_vec(n, [1.0, 1.0, 0.0, 0.0])
  let d = ext_from_vec(n, [0.0, 1.0, 1.0, 0.0])
  let cd = ext_wedge(c, d)
  let dc = ext_wedge(d, c)
  for i = 0; i < cd.data.length(); i = i + 1 {
    assert_true((cd.data[i] + dc.data[i]).abs() < 1.0e-12)
  }
}

test "exterior: 誘導內積與面積" {
  let u = [1.0, 0.0, 0.0]
  let v = [0.0, 1.0, 0.0]
  let w = ext_wedge_vecs(3, u, v)
  assert_true((ext_inner(w, w) - 1.0).abs() < 1.0e-12)
  // 一般情形：‖u∧v‖² = ‖u‖²‖v‖² − ⟨u,v⟩²（Lagrange 恆等式）
  let u2 = [1.0, 2.0, 2.0]
  let v2 = [2.0, -1.0, 1.0]
  let lhs = ext_norm(ext_wedge_vecs(3, u2, v2))
  let rhs = (
    euc_inner(u2, u2) * euc_inner(v2, v2) - euc_inner(u2, v2) * euc_inner(u2, v2),
  ).sqrt()
  assert_true((lhs - rhs).abs() < 1.0e-10)
}

test "exterior: Hodge 星性質" {
  // ⋆1 = e₁∧…∧e_n
  let one = ext_zero(3, 0)
  one.data[0] = 1.0
  let s = ext_hodge_star(one)
  assert_eq(s.k, 3)
  assert_true((s.data[0] - 1.0).abs() < 1.0e-12)
  // ⋆⋆ = (−1)^{k(n−k)}
  let n = 3
  let w = ext_wedge_vecs(n, [1.0, 2.0, 0.0], [0.0, 1.0, 3.0]) // 2-向量
  let ww = ext_hodge_star(ext_hodge_star(w))
  // k=2, n=3 → (−1)^2 = +1
  for i = 0; i < w.data.length(); i = i + 1 {
    assert_true((ww.data[i] - w.data[i]).abs() < 1.0e-12)
  }
  // k=1, n=3 → (−1)^2 = +1
  let v = ext_from_vec(3, [1.0, -2.0, 3.0])
  let vv = ext_hodge_star(ext_hodge_star(v))
  for i = 0; i < v.data.length(); i = i + 1 {
    assert_true((vv.data[i] - v.data[i]).abs() < 1.0e-12)
  }
}

test "exterior: 內縮恆等式 i_v(v∧w) = ‖v‖²w − ⟨v,w⟩v" {
  let v = [1.0, 2.0, -1.0]
  let w = [3.0, 0.0, 1.0]
  let lhs = ext_interior(v, ext_wedge_vecs(3, v, w))
  let rhs = euc_add(euc_scale(euc_inner(v, v), w), euc_scale(-euc_inner(v, w), v))
  assert_true(euc_eq_eps(lhs.data, rhs, 1.0e-10))
}

test "exterior: 外積行列式等於矩陣行列式（4×4）" {
  let vs : Array[Array[Double]] = [
    [1.0, 2.0, 0.0, 1.0],
    [0.0, 1.0, 1.0, 0.0],
    [2.0, 0.0, 1.0, 1.0],
    [1.0, 1.0, 1.0, 2.0],
  ]
  let m = mat_from_rows([vs[0], vs[1], vs[2], vs[3]])
  assert_true((ext_det(vs) - mat_det(m)).abs() < 1.0e-9)
}

test "exterior: 叉積由外代數導出且滿足右手系" {
  let i = [1.0, 0.0, 0.0]
  let j = [0.0, 1.0, 0.0]
  let k = [0.0, 0.0, 1.0]
  assert_true(euc_eq_eps(ext_cross3(i, j), k, 1.0e-12))
  assert_true(euc_eq_eps(ext_cross3(j, k), i, 1.0e-12))
  assert_true(euc_eq_eps(ext_cross3(k, i), j, 1.0e-12))
  assert_true(euc_eq_eps(ext_cross3(j, i), euc_neg(k), 1.0e-12))
  // 與直接公式一致
  let u = [1.0, -2.0, 3.0]
  let v = [2.0, 1.0, -1.0]
  assert_true(euc_eq_eps(ext_cross3(u, v), euc_cross3(u, v), 1.0e-12))
}

test "exterior: 體積之多重線性" {
  let u = [1.0, 0.0, 0.0]
  let v = [0.0, 2.0, 0.0]
  assert_true((ext_parallelogram_area(u, v) - 2.0).abs() < 1.0e-12)
  // 純量倍：面積乘 |c|
  assert_true((ext_parallelogram_area(euc_scale(3.0, u), v) - 6.0).abs() < 1.0e-12)
  // 剪切不變：v ↦ v + cu 面積不變
  assert_true(
    (ext_parallelogram_area(u, euc_add(v, euc_scale(5.0, u))) - 2.0).abs() < 1.0e-10,
  )
}
