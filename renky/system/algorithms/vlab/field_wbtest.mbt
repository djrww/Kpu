///|
// 複數域公理之自證

fn capprox(a : C, b : C) -> Bool {
  ceq_eps(a, b, 1.0e-12)
}

test "field: 加法交換律與結合律" {
  let a = cnew(1.0, -2.0)
  let b = cnew(3.0, 4.0)
  let c = cnew(-1.0, 0.5)
  assert_true(capprox(cadd(a, b), cadd(b, a)))
  assert_true(capprox(cadd(cadd(a, b), c), cadd(a, cadd(b, c))))
}

test "field: 加法零元與負元" {
  let a = cnew(2.5, 3.5)
  assert_true(capprox(cadd(a, czero()), a))
  assert_true(capprox(cadd(a, cneg(a)), czero()))
}

test "field: 乘法結合律與交換律" {
  let a = cnew(1.0, 2.0)
  let b = cnew(-2.0, 1.0)
  let c = cnew(0.5, -3.0)
  assert_true(capprox(cmul(cmul(a, b), c), cmul(a, cmul(b, c))))
  assert_true(capprox(cmul(a, b), cmul(b, a)))
}

test "field: 分配律" {
  let a = cnew(1.0, 2.0)
  let b = cnew(3.0, -1.0)
  let c = cnew(-2.0, 2.0)
  assert_true(capprox(cmul(a, cadd(b, c)), cadd(cmul(a, b), cmul(a, c))))
}

test "field: 乘法單位元與逆元" {
  let a = cnew(-1.0, 3.0)
  assert_true(capprox(cmul(a, cone()), a))
  assert_true(capprox(cmul(a, cinv(a)), cone()))
}

test "field: 共軛性質" {
  let a = cnew(2.0, -5.0)
  let b = cnew(1.0, 1.0)
  assert_true(capprox(cconj(cadd(a, b)), cadd(cconj(a), cconj(b))))
  assert_true(capprox(cconj(cmul(a, b)), cmul(cconj(a), cconj(b))))
  assert_true(capprox(cmul(a, cconj(a)), cfrom_real(cnorm2(a))))
}

test "field: 模長乘法性" {
  let a = cnew(3.0, 4.0)
  let b = cnew(-1.0, 2.0)
  assert_true((cabs(cmul(a, b)) - cabs(a) * cabs(b)).abs() < 1.0e-12)
  assert_true((cabs(a) - 5.0).abs() < 1.0e-12)
}

test "field: 除法與逆元一致" {
  let a = cnew(1.0, 1.0)
  let b = cnew(2.0, -1.0)
  assert_true(capprox(cdiv(a, b), cmul(a, cinv(b))))
  assert_true(capprox(cdiv(b, b), cone()))
}

test "field: 極式與歐拉公式" {
  let z = cpolar(2.0, @math.PI / 4.0)
  assert_true((cabs(z) - 2.0).abs() < 1.0e-12)
  assert_true((carg(z) - @math.PI / 4.0).abs() < 1.0e-12)
  // e^{iπ} = −1
  assert_true(
    ceq_eps(cexp(cscale(@math.PI, ci())), cnew(-1.0, 0.0), 1.0e-10),
  )
}

test "field: 平方根與指數互逆" {
  let z = cnew(2.0, -3.0)
  let w = csqrt(z)
  assert_true(ceq_eps(csqr(w), z, 1.0e-9))
  // exp(log) 循環性：e^{iθ} 之模長為 1
  let e = cexp(cnew(0.0, 1.234))
  assert_true((cabs(e) - 1.0).abs() < 1.0e-12)
}

test "field: 實數嵌入保運算" {
  let x = 2.5
  let y = -1.5
  assert_true(capprox(cadd(cfrom_real(x), cfrom_real(y)), cfrom_real(x + y)))
  assert_true(capprox(cmul(cfrom_real(x), cfrom_real(y)), cfrom_real(x * y)))
}
