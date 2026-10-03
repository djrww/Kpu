///|
// 複內積空間與 Hermitian 矩陣定理之自證

test "hermitian: 內積共軛對稱 ⟨u,v⟩ = conj⟨v,u⟩" {
  let u = [cnew(1.0, 2.0), cnew(-1.0, 0.0)]
  let v = [cnew(0.0, 1.0), cnew(3.0, -2.0)]
  assert_true(ceq_eps(cvec_inner(u, v), cconj(cvec_inner(v, u)), 1.0e-12))
}

test "hermitian: 第二變量線性" {
  let u = [cnew(1.0, 1.0)]
  let v = [cnew(2.0, 0.0)]
  let w = [cnew(0.0, 3.0)]
  let a = cnew(1.0, -1.0)
  let b = cnew(0.5, 2.0)
  let lhs = cvec_inner(u, cvec_add(cvec_cscale(a, v), cvec_cscale(b, w)))
  // ⟨u,v⟩ = Σ conj(u)v ⇒ 第二變量線性
  let rhs2 = cadd(cmul(a, cvec_inner(u, v)), cmul(b, cvec_inner(u, w)))
  assert_true(ceq_eps(lhs, rhs2, 1.0e-12))
}

test "hermitian: 正定性" {
  let v = [cnew(1.0, 1.0), cnew(-2.0, 3.0)]
  assert_true(cvec_norm2(v) > 0.0)
  assert_true(cvec_norm2(cvec_zero(2)) == 0.0)
  // ⟨v,v⟩ 必為實
  assert_true(cvec_inner(v, v).im.abs() < 1.0e-12)
}

test "hermitian: 範數具體值" {
  let v = [cnew(3.0, 4.0)] // |3+4i| = 5
  assert_true((cvec_norm(v) - 5.0).abs() < 1.0e-12)
}

test "hermitian: 複勾股定理" {
  let u = [cnew(1.0, 1.0), cnew(0.0, 0.0)]
  let v = [cnew(0.0, 0.0), cnew(1.0, -1.0)]
  // 正交
  assert_true(cvec_is_orthogonal(u, v, 1.0e-12))
  let s = cvec_add(u, v)
  assert_true(
    (cvec_norm2(s) - cvec_norm2(u) - cvec_norm2(v)).abs() < 1.0e-12,
  )
}

test "hermitian: 複 Gram–Schmidt 生成空間不變" {
  let vs : Array[Array[C]] = [
    [cone(), cnew(1.0, 1.0)],
    [cnew(1.0, -1.0), cnew(2.0, 1.0)],
  ]
  let onb = cvec_orthonormal(vs, 1.0e-10)
  assert_eq(onb.length(), 2)
  // 標準正交
  for i = 0; i < 2; i = i + 1 {
    for j = 0; j < 2; j = j + 1 {
      let g = cvec_inner(onb[i], onb[j])
      let want = if i == j { 1.0 } else { 0.0 }
      assert_true((g.re - want).abs() < 1.0e-9)
      assert_true(g.im.abs() < 1.0e-9)
    }
  }
}

test "hermitian: dagger 對合與反序律" {
  let a = cmat_new(2, 2)
  cmat_set(a, 0, 0, cnew(1.0, 1.0))
  cmat_set(a, 0, 1, cnew(2.0, 0.0))
  cmat_set(a, 1, 0, cnew(0.0, -1.0))
  cmat_set(a, 1, 1, cnew(3.0, 2.0))
  let b = cmat_new(2, 2)
  cmat_set(b, 0, 0, cnew(0.0, 1.0))
  cmat_set(b, 0, 1, cnew(1.0, 1.0))
  cmat_set(b, 1, 0, cnew(-2.0, 0.0))
  cmat_set(b, 1, 1, cnew(1.0, 0.0))
  // (A†)† = A
  let dd = cmat_dagger(cmat_dagger(a))
  let mut ok = true
  for i = 0; i < 4; i = i + 1 {
    if !ceq_eps(dd.data[i], a.data[i], 1.0e-12) {
      ok = false
    }
  }
  assert_true(ok)
  // (AB)† = B†A†
  let lhs = cmat_dagger(cmat_mul(a, b))
  let rhs = cmat_mul(cmat_dagger(b), cmat_dagger(a))
  let mut ok2 = true
  for i = 0; i < 4; i = i + 1 {
    if !ceq_eps(lhs.data[i], rhs.data[i], 1.0e-12) {
      ok2 = false
    }
  }
  assert_true(ok2)
}

test "hermitian: Hermitian 判定" {
  let h = cmat_new(2, 2)
  cmat_set(h, 0, 0, cfrom_real(1.0))
  cmat_set(h, 1, 1, cfrom_real(2.0))
  cmat_set(h, 0, 1, cnew(3.0, 4.0))
  cmat_set(h, 1, 0, cnew(3.0, -4.0))
  assert_true(cmat_is_hermitian(h, 1.0e-12))
  cmat_set(h, 1, 0, cnew(3.0, 4.0)) // 破壞
  assert_false(cmat_is_hermitian(h, 1.0e-12))
}

test "hermitian: Hermitian 特徵值為實且特徵向量酉" {
  // 3×3 Hermitian
  let h = cmat_new(3, 3)
  cmat_set(h, 0, 0, cfrom_real(2.0))
  cmat_set(h, 1, 1, cfrom_real(2.0))
  cmat_set(h, 2, 2, cfrom_real(2.0))
  cmat_set(h, 0, 1, cnew(1.0, 1.0))
  cmat_set(h, 1, 0, cnew(1.0, -1.0))
  cmat_set(h, 1, 2, cnew(0.0, -1.0))
  cmat_set(h, 2, 1, cnew(0.0, 1.0))
  let (vals, q) = cmat_jacobi_hermitian(h, 80)
  assert_true(cmat_unitary_error(q) < 1.0e-8)
  // 特徵方程 H q_j = λ_j q_j
  for j = 0; j < 3; j = j + 1 {
    let v = mat_col_to_cvec(q, j)
    let hv = cmat_vec(h, v)
    let lv = cvec_scale(vals[j], v)
    let mut ok = true
    for i = 0; i < 3; i = i + 1 {
      if !ceq_eps(hv[i], lv[i], 1.0e-7) {
        ok = false
      }
    }
    assert_true(ok)
  }
  // 跡不變
  let mut vsum2 = 0.0
  for v in vals {
    vsum2 = vsum2 + v
  }
  assert_true((vsum2 - 6.0).abs() < 1.0e-8)
}

/// 輔助：CMat 之行取出為 Array[C]
fn mat_col_to_cvec(m : CMat, j : Int) -> Array[C] {
  let out : Array[C] = []
  for i = 0; i < m.rows; i = i + 1 {
    out.push(cmat_get(m, i, j))
  }
  out
}

test "hermitian: Rayleigh 商範圍 [λmin, λmax]" {
  let h = cmat_new(2, 2)
  cmat_set(h, 0, 0, cfrom_real(1.0))
  cmat_set(h, 1, 1, cfrom_real(3.0))
  cmat_set(h, 0, 1, cnew(0.0, 1.0))
  cmat_set(h, 1, 0, cnew(0.0, -1.0))
  let (vals, _) = cmat_jacobi_hermitian(h, 60)
  let lo = if vals[0] < vals[1] { vals[0] } else { vals[1] }
  let hi = if vals[0] < vals[1] { vals[1] } else { vals[0] }
  let x = [cnew(1.0, 2.0), cnew(-1.0, 0.0)]
  let rq = cmat_rayleigh(h, x)
  assert_true(rq >= lo - 1.0e-8 && rq <= hi + 1.0e-8)
}
