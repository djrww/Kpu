// ===========================================================================
//  義務自證 · layer 2 — EXHAUSTIVE proofs for the finite-field skeleton.
//
//  Each test states a theorem and then verifies it on the *whole* of the
//  finite structure it quantifies over.  Since the domains are finite and
//  enumerated completely, these are not statistical tests: they are complete
//  decision procedures for the stated sentences.
// ===========================================================================

fn check(cond : Bool, msg : String) -> Unit {
  if cond {
    ()
  } else {
    println("CHECK FAILED: " + msg)
    abort("[law failed] " + msg)
  }
}

///| Every polynomial over F_p of degree < `len` (all p^len of them).
fn all_polys(p : Int, len : Int) -> Array[FPoly] {
  let out : Array[FPoly] = []
  let total = @arith.pow_int(p, len)
  for k = 0; k < total; k = k + 1 {
    let ds = @arith.digits(k, p)
    let c : Array[Int] = Array::make(len, 0)
    for i = 0; i < len && i < ds.length(); i = i + 1 {
      c[i] = ds[i]
    }
    out.push(fpoly(p, c))
  }
  out
}

// ---------------------------------------------------------------------------
//  1. F_p is a field
// ---------------------------------------------------------------------------

test "F_p ring axioms hold for every triple of elements, p in {2,3,5,7}" {
  let ps = [2, 3, 5, 7]
  for k = 0; k < ps.length(); k = k + 1 {
    let p = ps[k]
    let els = fp_elements(p)
    let one = Fp::one(p)
    let zero = Fp::zero(p)
    for i = 0; i < els.length(); i = i + 1 {
      let a = els[i]
      check(a.add(zero).equal(a), "additive identity")
      check(a.mul(one).equal(a), "multiplicative identity")
      check(a.add(a.neg()).equal(zero), "additive inverse")
      if !(a.is_zero()) {
        check(a.mul(a.inv()).equal(one), "multiplicative inverse")
      }
      for j = 0; j < els.length(); j = j + 1 {
        let b = els[j]
        check(a.add(b).equal(b.add(a)), "add comm")
        check(a.mul(b).equal(b.mul(a)), "mul comm")
        for t = 0; t < els.length(); t = t + 1 {
          let c = els[t]
          check(a.add(b).add(c).equal(a.add(b.add(c))), "add assoc")
          check(a.mul(b).mul(c).equal(a.mul(b.mul(c))), "mul assoc")
          check(
            a.mul(b.add(c)).equal(a.mul(b).add(a.mul(c))),
            "distributivity p=" + p.to_string(),
          )
        }
      }
    }
  }
}

test "Frobenius on F_p is the identity: a^p == a for every a, every prime p <= 100" {
  for p = 2; p <= 100; p = p + 1 {
    if !(@arith.is_prime(p)) {
      continue
    }
    for a = 0; a < p; a = a + 1 {
      let x = fp_mk(p, a)
      check(x.frobenius().equal(x), "frobenius id p=" + p.to_string())
      check(x.pow(p - 1).equal(Fp::one(p)) || x.is_zero(), "Lagrange")
    }
  }
}

test "F_p^* is cyclic: a generator exists and there are exactly phi(p-1) of them, p <= 60" {
  for p = 2; p <= 60; p = p + 1 {
    if !(@arith.is_prime(p)) {
      continue
    }
    let els = fp_elements(p)
    let mut gens = 0
    for i = 1; i < els.length(); i = i + 1 {
      let a = els[i]
      // order of a is p-1 iff a^((p-1)/q) != 1 for every prime q | p-1
      let f = @arith.factorize(p - 1)
      let mut is_gen = true
      for j = 0; j < f.length(); j = j + 1 {
        if a.pow((p - 1) / f[j].0).equal(Fp::one(p)) {
          is_gen = false
        }
      }
      if is_gen {
        gens = gens + 1
        // the powers of a hit every nonzero element exactly once
        let seen : Array[Bool] = Array::make(p, false)
        let mut cur = Fp::one(p)
        for k = 0; k < p - 1; k = k + 1 {
          check(!(seen[cur.repr()]), "generator powers distinct")
          seen[cur.repr()] = true
          cur = cur.mul(a)
        }
        check(cur.equal(Fp::one(p)), "generator closes")
      }
    }
    check(gens == @arith.totient(p - 1), "phi(p-1) generators, p=" + p.to_string())
  }
}

// ---------------------------------------------------------------------------
//  2. F_p[x]: Euclidean domain
// ---------------------------------------------------------------------------

test "F_p[x] division: f == q*g + r with deg r < deg g, all f,g of degree < 4 over F_2, F_3" {
  let ps = [2, 3]
  for k = 0; k < ps.length(); k = k + 1 {
    let p = ps[k]
    let fs = all_polys(p, 4)
    for i = 0; i < fs.length(); i = i + 1 {
      for j = 0; j < fs.length(); j = j + 1 {
        let g = fs[j]
        if g.is_zero() {
          continue
        }
        let (q, r) = fs[i].divmod(g)
        check(
          q.mul(g).add(r).equal(fs[i]),
          "divmod identity p=" + p.to_string(),
        )
        check(r.degree() < g.degree(), "remainder degree")
      }
    }
  }
}

test "Bezout in F_p[x]: s*f + t*g == gcd(f,g), all pairs of degree < 4 over F_2, F_3" {
  let ps = [2, 3]
  for k = 0; k < ps.length(); k = k + 1 {
    let p = ps[k]
    let fs = all_polys(p, 4)
    for i = 0; i < fs.length(); i = i + 1 {
      for j = 0; j < fs.length(); j = j + 1 {
        if fs[i].is_zero() && fs[j].is_zero() {
          continue
        }
        let (g, s, t) = fp_egcd(fs[i], fs[j])
        check(
          s.mul(fs[i]).add(t.mul(fs[j])).equal(g),
          "poly bezout p=" + p.to_string(),
        )
        if !(g.is_zero()) {
          check(fs[i].rem(g).is_zero(), "gcd divides f")
          check(fs[j].rem(g).is_zero(), "gcd divides g")
        }
      }
    }
  }
}

// ---------------------------------------------------------------------------
//  3. Irreducibles: Ben-Or, Gauss, and the fundamental theorem
// ---------------------------------------------------------------------------

test "Ben-Or agrees with brute-force irreducibility for every monic f of degree <= 4 over F_2, F_3" {
  let ps = [2, 3]
  for k = 0; k < ps.length(); k = k + 1 {
    let p = ps[k]
    for d = 1; d <= 4; d = d + 1 {
      let fs = all_polys(p, d)
      for i = 0; i < fs.length(); i = i + 1 {
        let f = fs[i].add(fp_monom(p, d, 1))
        // brute force: f reducible iff some monic g of degree in [1,d/2] divides it
        let mut brute_irr = true
        for e = 1; e <= d / 2; e = e + 1 {
          let gs = all_polys(p, e)
          for t = 0; t < gs.length(); t = t + 1 {
            let g = gs[t].add(fp_monom(p, e, 1))
            if f.rem(g).is_zero() {
              brute_irr = false
            }
          }
        }
        check(
          fp_is_irreducible(f) == brute_irr,
          "ben-or p=" + p.to_string() + " d=" + d.to_string() + " f=" + f.to_string(),
        )
      }
    }
  }
}

test "Gauss: number of monic irreducibles of degree d == (1/d) sum mu(d/k) p^k" {
  let ps = [2, 3, 5]
  for k = 0; k < ps.length(); k = k + 1 {
    let p = ps[k]
    for d = 1; d <= 4; d = d + 1 {
      let enumerated = fp_all_irred(p, d).length()
      check(
        enumerated == fp_count_monic_irred(p, d),
        "gauss p=" + p.to_string() + " d=" + d.to_string(),
      )
    }
  }
}

test "Fundamental theorem: x^(p^n) - x == product of all monic irreducibles of degree | n" {
  let cases = [(2, 1), (2, 2), (2, 3), (2, 4), (3, 1), (3, 2), (3, 3), (5, 1), (5, 2)]
  for k = 0; k < cases.length(); k = k + 1 {
    let p = cases[k].0
    let n = cases[k].1
    let x = fp_x(p)
    let big = fp_monom(p, @arith.pow_int(p, n), 1).sub(x)
    let mut prod = fp_one(p)
    for d = 1; d <= n; d = d + 1 {
      if n % d != 0 {
        continue
      }
      let irrs = fp_all_irred(p, d)
      for i = 0; i < irrs.length(); i = i + 1 {
        prod = prod.mul(irrs[i])
      }
    }
    check(big.equal(prod), "x^(p^n)-x factorisation p=" + p.to_string() + " n=" + n.to_string())
    check(big.degree() == @arith.pow_int(p, n), "degree of x^(p^n)-x")
  }
}

test "A26 factorisation is correct for every monic f of degree <= 4 over F_2 and degree <= 3 over F_3" {
  let p2 = 2
  for d = 1; d <= 4; d = d + 1 {
    let fs = all_polys(p2, d)
    for i = 0; i < fs.length(); i = i + 1 {
      let f = fs[i].add(fp_monom(p2, d, 1))
      let fac = fp_factor(f)
      let mut prod = fp_one(p2)
      for j = 0; j < fac.length(); j = j + 1 {
        check(fp_is_irreducible(fac[j].0), "factor irreducible " + f.to_string())
        check(fac[j].1 >= 1, "multiplicity positive")
        for _k = 0; _k < fac[j].1; _k = _k + 1 {
          prod = prod.mul(fac[j].0)
        }
      }
      check(prod.equal(f), "factor product " + f.to_string())
    }
  }
  let p3 = 3
  for d = 1; d <= 3; d = d + 1 {
    let fs = all_polys(p3, d)
    for i = 0; i < fs.length(); i = i + 1 {
      let f = fs[i].add(fp_monom(p3, d, 1))
      let fac = fp_factor(f)
      let mut prod = fp_one(p3)
      for j = 0; j < fac.length(); j = j + 1 {
        check(fp_is_irreducible(fac[j].0), "factor irreducible3 " + f.to_string())
        for _k = 0; _k < fac[j].1; _k = _k + 1 {
          prod = prod.mul(fac[j].0)
        }
      }
      check(prod.equal(f), "factor product3 " + f.to_string())
    }
  }
}

test "A24/A25 ddf then edf reproduces the squarefree part, monic f of degree <= 4 over F_2" {
  let p = 2
  for d = 1; d <= 4; d = d + 1 {
    let fs = all_polys(p, d)
    for i = 0; i < fs.length(); i = i + 1 {
      let f = fs[i].add(fp_monom(p, d, 1))
      let g = f.gcd(f.deriv())
      if !(g.is_zero()) && g.degree() > 0 {
        continue // not squarefree, handled by the previous test
      }
      let parts = fp_ddf(f)
      let mut prod = fp_one(p)
      let mut total = 0
      for j = 0; j < parts.length(); j = j + 1 {
        let irr = fp_edf(parts[j].1, parts[j].0)
        for t = 0; t < irr.length(); t = t + 1 {
          check(irr[t].degree() == parts[j].0, "edf degrees")
          check(fp_is_irreducible(irr[t]), "edf irreducible")
          prod = prod.mul(irr[t])
          total = total + 1
        }
      }
      check(prod.equal(f), "ddf+edf product " + f.to_string())
      check(total == fp_factor(f).length(), "ddf+edf count")
    }
  }
}

// ---------------------------------------------------------------------------
//  4. F_{p^n} is a field, Frobenius generates Gal(F_{p^n}/F_p)
// ---------------------------------------------------------------------------

test "F_{p^n} axioms hold for every triple of elements: F_4, F_9, F_25" {
  let cases = [(2, 2), (3, 2), (5, 2)]
  for k = 0; k < cases.length(); k = k + 1 {
    let fld = fpn_create(cases[k].0, cases[k].1)
    let els = fld.elements()
    let one = fld.one()
    let zero = fld.zero()
    check(els.length() == fld.size(), "field size")
    // the defining relation: irr(gen) == 0
    let g = fld.gen()
    let mut val = zero
    let mut pw = one
    for i = 0; i <= fld.irr.degree(); i = i + 1 {
      val = val.add(pw.mul(fld.lift(fld.irr.coeff(i))))
      pw = pw.mul(g)
    }
    check(val.is_zero(), "irr(gen) == 0")
    for i = 0; i < els.length(); i = i + 1 {
      let a = els[i]
      check(a.add(zero).equal(a), "Fpn add id")
      check(a.mul(one).equal(a), "Fpn mul id")
      check(a.add(a.neg()).equal(zero), "Fpn add inv")
      if !(a.is_zero()) {
        check(a.mul(a.inv()).equal(one), "Fpn mul inv")
      }
      for j = 0; j < els.length(); j = j + 1 {
        let b = els[j]
        check(a.mul(b).equal(b.mul(a)), "Fpn comm")
        for t = 0; t < els.length(); t = t + 1 {
          let c = els[t]
          check(a.add(b).add(c).equal(a.add(b.add(c))), "Fpn add assoc")
          check(a.mul(b).mul(c).equal(a.mul(b.mul(c))), "Fpn mul assoc")
          check(a.mul(b.add(c)).equal(a.mul(b).add(a.mul(c))), "Fpn distrib")
        }
      }
    }
  }
}

test "Freshman's dream and Frobenius order in F_{p^n}: (a+b)^p == a^p+b^p, F^n == id" {
  let cases = [(2, 3), (3, 2), (2, 4), (5, 2), (2, 2), (3, 3)]
  for k = 0; k < cases.length(); k = k + 1 {
    let fld = fpn_create(cases[k].0, cases[k].1)
    let els = fld.elements()
    for i = 0; i < els.length(); i = i + 1 {
      let a = els[i]
      check(a.frobenius_k(fld.n).equal(a), "F^n == id")
      if fld.n > 1 {
        check(!(a.frobenius_k(fld.n - 1).equal(a)) || in_subfield(a, fld), "F^(n-1) fixes only F_p")
      }
      for j = 0; j < els.length(); j = j + 1 {
        let b = els[j]
        check(
          a.add(b).frobenius().equal(a.frobenius().add(b.frobenius())),
          "freshman dream",
        )
        check(
          a.mul(b).frobenius().equal(a.frobenius().mul(b.frobenius())),
          "frobenius multiplicative",
        )
      }
    }
  }
}

fn in_subfield(a : Fpn, fld : Fpn) -> Bool {
  // a lies in the prime field iff all coefficients except the constant vanish
  for i = 1; i < fld.n; i = i + 1 {
    if a.coeffs()[i] != 0 {
      return false
    }
  }
  true
}

test "Galois fixed fields: |{x : F^d(x) == x}| == p^gcd(d,n)" {
  let cases = [(2, 4), (3, 2), (2, 6), (2, 3), (5, 2)]
  for k = 0; k < cases.length(); k = k + 1 {
    let p = cases[k].0
    let n = cases[k].1
    let fld = fpn_create(p, n)
    if fld.size() > 200 {
      continue
    }
    let els = fld.elements()
    for d = 1; d <= n; d = d + 1 {
      let mut cnt = 0
      for i = 0; i < els.length(); i = i + 1 {
        if els[i].frobenius_k(d).equal(els[i]) {
          cnt = cnt + 1
        }
      }
      check(
        cnt == @arith.pow_int(p, @arith.gcd(d, n)),
        "fixed field p=" + p.to_string() + " n=" + n.to_string() + " d=" + d.to_string(),
      )
    }
  }
}

test "Trace is F_p-linear and Norm is multiplicative, verified on all of F_8, F_9, F_16" {
  let cases = [(2, 3), (3, 2), (2, 4)]
  for k = 0; k < cases.length(); k = k + 1 {
    let p = cases[k].0
    let fld = fpn_create(p, cases[k].1)
    let els = fld.elements()
    for i = 0; i < els.length(); i = i + 1 {
      let a = els[i]
      for j = 0; j < els.length(); j = j + 1 {
        let b = els[j]
        check(
          @kernel.addmod(a.trace(), b.trace(), p) == a.add(b).trace(),
          "trace additive",
        )
        check(
          @kernel.mulmod(a.norm(), b.norm(), p) == a.mul(b).norm(),
          "norm multiplicative",
        )
      }
      for c = 0; c < p; c = c + 1 {
        check(
          @kernel.mulmod(c, a.trace(), p) == a.mul(fld.lift(c)).trace(),
          "trace F_p-linear",
        )
        check(
          @kernel.mulmod(@arith.powmod(c, fld.n, p), a.norm(), p) ==
            a.mul(fld.lift(c)).norm(),
          "norm scaling N(c*a) == c^n * N(a)",
        )
      }
      // Norm as the (p^n-1)/(p-1)-th power
      if !(a.is_zero()) {
        let e = (fld.size() - 1) / (p - 1)
        check(a.pow(e).coeffs()[0] == a.norm(), "N(a) == a^((p^n-1)/(p-1))")
      }
    }
  }
}

test "Wedderburn: F_{p^n}^* is cyclic of order p^n-1 with phi(p^n-1) generators" {
  let cases = [(2, 2), (2, 3), (3, 2), (2, 4), (5, 2), (2, 5), (3, 3)]
  for k = 0; k < cases.length(); k = k + 1 {
    let fld = fpn_create(cases[k].0, cases[k].1)
    let q = fld.size() - 1
    let g = fld.primitive_element()
    check(g.order() == q, "primitive element order")
    // count all generators
    let els = fld.elements()
    let mut cnt = 0
    for i = 0; i < els.length(); i = i + 1 {
      if !(els[i].is_zero()) && els[i].order() == q {
        cnt = cnt + 1
      }
    }
    check(cnt == @arith.totient(q), "number of generators = phi(q)")
    // discrete log round trip with the verified binary search inside
    for e = 0; e < q; e = e + 1 {
      let h = g.pow(e)
      check(fld.dlog(g, h) == e, "dlog round trip e=" + e.to_string())
    }
  }
}

test "Minimal polynomial of x over F_p has degree equal to the Frobenius orbit and divides the irreducible" {
  let cases = [(2, 3), (2, 4), (3, 2), (3, 3), (5, 2)]
  for k = 0; k < cases.length(); k = k + 1 {
    let fld = fpn_create(cases[k].0, cases[k].1)
    let els = fld.elements()
    for i = 0; i < els.length(); i = i + 1 {
      let a = els[i]
      let conj = a.conjugates()
      // orbit size = number of *distinct* conjugates
      let mut orbit = 0
      for t = 0; t < conj.length(); t = t + 1 {
        let mut fresh = true
        for u = 0; u < t; u = u + 1 {
          if conj[u].equal(conj[t]) {
            fresh = false
          }
        }
        if fresh {
          orbit = orbit + 1
        }
      }
      let mp = a.minpoly()
      check(mp.degree() == orbit, "minpoly degree " + a.to_string())
      check(mp.is_monic(), "minpoly monic")
      check(fp_is_irreducible(mp), "minpoly irreducible")
      // a is a root of its minimal polynomial
      let mut val = fld.zero()
      let mut pw = fld.one()
      for j = 0; j <= mp.degree(); j = j + 1 {
        val = val.add(pw.mul(fld.lift(mp.coeff(j))))
        pw = pw.mul(a)
      }
      check(val.is_zero(), "minpoly vanishes at a")
      // every minimal polynomial divides x^(p^n) - x
      let big = fp_monom(fld.p, fld.size(), 1).sub(fp_x(fld.p))
      check(
        big.rem(mp).is_zero(),
        "minpoly divides x^(p^n)-x: " + mp.to_string(),
      )
      // and when a generates the field over F_p, its minimal polynomial is one
      // of the monic irreducibles of degree n enumerated independently
      if orbit == fld.n {
        let irrs = fp_all_irred(fld.p, fld.n)
        let mut found = false
        for t = 0; t < irrs.length(); t = t + 1 {
          if irrs[t].equal(mp) {
            found = true
          }
        }
        check(found, "minpoly is an enumerated irreducible")
      }
      // Vieta: sum of roots == -coeff_{d-1}
      let mut sum = fld.zero()
      for t = 0; t < conj.length(); t = t + 1 {
        if t < orbit {
          sum = sum.add(conj[t])
        }
      }
      check(
        sum.coeffs()[0] == @kernel.negmod(mp.coeff(mp.degree() - 1), fld.p),
        "viete trace",
      )
    }
  }
}

test "Subfield lattice: F_{p^d} <= F_{p^n} iff d | n (via fixed field of F^d)" {
  let p = 2
  let n = 6
  let fld = fpn_create(p, n)
  for d = 1; d <= n; d = d + 1 {
    let mut cnt = 0
    for k = 0; k < fld.size(); k = k + 1 {
      let a = fpn_from_index(fld, k)
      if a.frobenius_k(d).equal(a) {
        cnt = cnt + 1
      }
    }
    let expect = if n % d == 0 {
      @arith.pow_int(p, d)
    } else {
      @arith.pow_int(p, @arith.gcd(d, n))
    }
    check(cnt == expect, "subfield d=" + d.to_string())
  }
}

test "F_2[x]/(x^2+x+1) == F_4 concretely, and every element satisfies a^(2^2) == a" {
  let f = fpoly(2, [1, 1, 1])
  check(fp_is_irreducible(f), "x^2+x+1 irreducible over F_2")
  let fld = fpn_of(2, 2, f)
  let els = fld.elements()
  check(els.length() == 4, "four elements")
  for i = 0; i < els.length(); i = i + 1 {
    check(els[i].pow(4).equal(els[i]), "a^4 == a")
    check(els[i].frobenius().frobenius().equal(els[i]), "F^2 == id")
  }
  let g = fld.gen()
  check(g.mul(g).equal(g.add(fld.one())), "w^2 == w + 1 in F_4")
}
