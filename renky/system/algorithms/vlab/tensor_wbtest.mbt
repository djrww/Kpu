///|
// 張量空間定理之自證

test "tensor: 雙線性展開" {
  // (u1+u2)⊗v = u1⊗v + u2⊗v
  let u1 = [1.0, 0.0]
  let u2 = [0.0, 2.0]
  let v = [3.0, 4.0, 5.0]
  let lhs = pure_tensor(euc_add(u1, u2), v)
  let rhs = tensor_add(pure_tensor(u1, v), pure_tensor(u2, v))
  assert_true((tensor_inner(lhs, lhs) > 0.0))
  assert_eq(lhs.data, rhs.data)
  // 純量結合 (cu)⊗v = u⊗(cv) = c(u⊗v)
  let c = 2.5
  let a1 = pure_tensor(euc_scale(c, u1), v)
  let a2 = pure_tensor(u1, euc_scale(c, v))
  let a3 = tensor_scale(c, pure_tensor(u1, v))
  assert_eq(a1.data, a2.data)
  assert_eq(a2.data, a3.data)
}

test "tensor: 內積公理" {
  let ip = tensor_ip(2, 3)
  let a = pure_tensor([1.0, 2.0], [1.0, 0.0, 1.0])
  let b = pure_tensor([0.0, 1.0], [2.0, 2.0, 0.0])
  assert_true(((ip.inner)(a, b) - (ip.inner)(b, a)).abs() < 1.0e-12)
  assert_true((ip.inner)(a, a) > 0.0)
  // 具體值：⟨u⊗v, s⊗w⟩ = ⟨u,s⟩⟨v,w⟩ = 2·2 = 4
  assert_true(((ip.inner)(a, b) - 4.0).abs() < 1.0e-12)
}

test "tensor: 收縮之線性與循環性" {
  let t1 = pure_tensor([1.0, 2.0], [3.0, 4.0])
  let t2 = pure_tensor([1.0, 0.0], [1.0, 1.0])
  let t = tensor_add(t1, t2)
  assert_true((tensor_contract(t) - tensor_contract(t1) - tensor_contract(t2)).abs() < 1.0e-12)
  // tr = Σ t_ii = (1·3+2·4) + (1·1+0·1) = 11 + 1 = 12
  assert_true((tensor_contract(t) - 12.0).abs() < 1.0e-12)
}

test "tensor: 單腿收縮 = 雙線性求值" {
  let t = pure_tensor([1.0, 2.0, 3.0], [4.0, 5.0])
  let phi = [1.0, 0.0, -1.0]
  let psi = [2.0, 1.0]
  let v1 = tensor_contract_first(t, phi)
  assert_true((euc_inner(v1, psi) - tensor_bilinear_eval(t, phi, psi)).abs() < 1.0e-12)
  let v2 = tensor_contract_second(t, psi)
  assert_true((euc_inner(phi, v2) - tensor_bilinear_eval(t, phi, psi)).abs() < 1.0e-12)
}

test "tensor: 同構保持內積（與矩陣空間）" {
  let ip = tensor_ip(2, 2)
  let a = tensor_add(
    pure_tensor([1.0, 1.0], [1.0, 0.0]),
    pure_tensor([0.0, 2.0], [1.0, 1.0]),
  )
  let m = tensor_to_mat(a)
  assert_true((mat_frob_inner(m, m) - tensor_inner(a, a)).abs() < 1.0e-12)
  let a2 = mat_to_tensor(m)
  assert_eq(a.data, a2.data)
}

test "tensor: 映射張量積 = Kronecker" {
  let a = mat_from_rows([[1.0, 2.0], [3.0, 4.0]])
  let b = mat_from_rows([[0.0, 1.0], [1.0, 0.0]])
  let t = tensor_of_maps(a, b)
  assert_eq(t.m, 4)
  assert_eq(t.n, 4)
  // (A⊗B)ᵀ = Aᵀ⊗Bᵀ 之驗證
  let t2 = tensor_of_maps(mat_transpose(a), mat_transpose(b))
  assert_true(mat_eq_eps(tensor_to_mat(t2), mat_transpose(tensor_to_mat(t)), 1.0e-12))
}

test "tensor: 平坦化秩 ≤ min(m,n) 且純張量秩 = 1" {
  let t = pure_tensor([1.0, 2.0, 3.0], [4.0, 5.0])
  assert_eq(tensor_flatten_rank(t, 1.0e-10), 1)
  let t2 = tensor_add(
    pure_tensor([1.0, 0.0, 0.0], [1.0, 0.0]),
    pure_tensor([0.0, 1.0, 0.0], [0.0, 1.0]),
  )
  assert_eq(tensor_flatten_rank(t2, 1.0e-10), 2)
}

test "tensor: Gram–Schmidt 於張量空間" {
  let ip = tensor_ip(2, 2)
  let vs : Array[Tensor] = [
    pure_tensor([1.0, 0.0], [1.0, 0.0]),
    tensor_add(
      pure_tensor([1.0, 0.0], [1.0, 0.0]),
      pure_tensor([0.0, 1.0], [0.0, 1.0]),
    ),
  ]
  let onb = ip_orthonormal(ip, vs, 1.0e-10)
  assert_eq(onb.length(), 2)
  assert_true((ip.inner)(onb[0], onb[1]).abs() < 1.0e-9)
  assert_true(((ip.inner)(onb[0], onb[0]) - 1.0).abs() < 1.0e-9)
}
