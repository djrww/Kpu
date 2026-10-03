///|
// 向量空間公理之自證（於多個具體實例上驗證）

test "vecspace R³: 全部八條公理" {
  let sp = euc_sp(3)
  let u = [1.0, -2.0, 0.5]
  let v = [0.0, 3.0, 1.0]
  let w = [-1.0, 1.0, 2.0]
  let a = 1.5
  let b = -2.0
  // 交換律
  assert_true(euc_eq_eps((sp.add)(u, v), (sp.add)(v, u), 1.0e-12))
  // 結合律
  assert_true(
    euc_eq_eps((sp.add)((sp.add)(u, v), w), (sp.add)(u, (sp.add)(v, w)), 1.0e-12),
  )
  // 零元
  assert_true(euc_eq_eps((sp.add)(u, sp.zero), u, 1.0e-12))
  // 負元
  assert_true(euc_eq_eps((sp.add)(u, (sp.neg)(u)), sp.zero, 1.0e-12))
  // 純量乘法結合
  assert_true(
    euc_eq_eps((sp.scale)(a * b, u), (sp.scale)(a, (sp.scale)(b, u)), 1.0e-12),
  )
  // 單位純量
  assert_true(euc_eq_eps((sp.scale)(1.0, u), u, 1.0e-12))
  // 分配律一
  assert_true(
    euc_eq_eps(
      (sp.scale)(a, (sp.add)(u, v)),
      (sp.add)((sp.scale)(a, u), (sp.scale)(a, v)),
      1.0e-12,
    ),
  )
  // 分配律二
  assert_true(
    euc_eq_eps(
      (sp.scale)(a + b, u),
      (sp.add)((sp.scale)(a, u), (sp.scale)(b, u)),
      1.0e-12,
    ),
  )
}

test "vecspace: 減法定義 u−v = u+(−v)" {
  let sp = euc_sp(3)
  let u = [2.0, 2.0, 2.0]
  let v = [1.0, 4.0, 0.0]
  assert_true(
    euc_eq_eps(vsub(sp, u, v), (sp.add)(u, (sp.neg)(v)), 1.0e-12),
  )
  assert_true(euc_eq_eps(vsub(sp, u, v), [1.0, -2.0, 2.0], 1.0e-12))
}

test "vecspace: axpy" {
  let sp = euc_sp(2)
  let x = [1.0, 1.0]
  let y = [0.5, -1.5]
  assert_true(euc_eq_eps(vaxpy(sp, 3.0, x, y), [3.5, 1.5], 1.0e-12))
}

test "vecspace: 線性組合具體值" {
  let sp = euc_sp(2)
  let r = vlinc(sp, [2.0, -1.0, 0.5], [[1.0, 0.0], [0.0, 1.0], [2.0, 2.0]])
  assert_true(euc_eq_eps(r, [3.0, 0.0], 1.0e-12))
}

test "vecspace: 中點" {
  let sp = euc_sp(2)
  assert_true(euc_eq_eps(vmidpoint(sp, [0.0, 0.0], [4.0, 6.0]), [2.0, 3.0], 1.0e-12))
}

test "vecspace: 多項式空間亦為向量空間" {
  let sp = poly_ip(3).sp
  let p = [1.0, 1.0]
  let q = [-1.0, 0.0, 2.0]
  assert_true(
    euc_eq_eps((sp.add)(p, q), poly_add(p, q), 1.0e-12),
  )
  assert_true(
    euc_eq_eps((sp.scale)(2.0, q), poly_scale(2.0, q), 1.0e-12),
  )
  assert_true(euc_eq_eps((sp.add)(p, (sp.neg)(p)), poly_zero(), 1.0e-12))
}

test "vecspace: 函數空間（加權）亦為向量空間" {
  let ip = fn_ip(4, [1.0, 1.0, 2.0, 2.0])
  let sp = ip.sp
  let f = [1.0, 2.0, 3.0, 4.0]
  let g = [4.0, 3.0, 2.0, 1.0]
  assert_true(euc_eq_eps((sp.add)(f, g), [5.0, 5.0, 5.0, 5.0], 1.0e-12))
  assert_true(euc_eq_eps((sp.add)(f, (sp.neg)(f)), sp.zero, 1.0e-12))
}

test "vecspace: 直和空間公理" {
  let sp = vec_prod(euc_sp(2), poly_ip(2).sp)
  let u = ([1.0, 1.0], [1.0, 2.0])
  let v = ([0.0, 2.0], [0.0, -1.0])
  let s = (sp.add)(u, v)
  assert_true(euc_eq_eps(s.0, [1.0, 3.0], 1.0e-12))
  assert_true(euc_eq_eps(s.1, [1.0, 1.0], 1.0e-12))
  let z = (sp.add)(u, (sp.neg)(u))
  assert_true(euc_eq_eps(z.0, [0.0, 0.0], 1.0e-12))
}

test "vecspace: vsum 與 vdouble_iter" {
  let sp = euc_sp(2)
  let s = vsum(sp, [[1.0, 0.0], [0.0, 1.0], [1.0, 1.0]])
  assert_true(euc_eq_eps(s, [2.0, 2.0], 1.0e-12))
  assert_true(euc_eq_eps(vdouble_iter(sp, [1.0, -1.0], 3), [8.0, -8.0], 1.0e-12))
}
