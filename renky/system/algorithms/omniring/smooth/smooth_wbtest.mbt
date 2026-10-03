// ===========================================================================
//  義務自證 · M1–M6 for the smooth (jet / Weil algebra) layer.  Everything is
//  checked by complete enumeration of the finite jet algebras involved.
// ===========================================================================

fn check(cond : Bool, msg : String) -> Unit {
  if cond {
    ()
  } else {
    println("CHECK FAILED: " + msg)
    abort("[smooth law failed] " + msg)
  }
}

///| Equality of jets up to degree k (the top coefficients of a truncated
///| algebra lose information, so identities of differential operators are
///| compared one or two degrees below the truncation).
fn equal_upto(a : Jet, b : Jet, k : Int) -> Bool {
  for i = 0; i <= Int::min(k, Int::min(a.order, b.order)); i = i + 1 {
    if a.coeffs[i] != b.coeffs[i] {
      return false
    }
  }
  true
}

///| Every jet of the given shape (p^(order+1) of them).
fn all_jets(p : Int, order : Int) -> Array[Jet] {
  let total = @arith.pow_int(p, order + 1)
  let out : Array[Jet] = []
  for k = 0; k < total; k = k + 1 {
    out.push(jet_of_coeffs(p, order, @arith.digits(k, p)))
  }
  out
}

///| Independent definition of composition: sum_i f_i * g^i.
fn compose_direct(f : Jet, g : Jet) -> Jet {
  let order = Int::min(f.order, g.order)
  let mut acc = jet_zero(f.p, order)
  let mut pw = jet_const(f.p, order, 1)
  let gg = jet_of_coeffs(f.p, order, g.coeffs)
  for i = 0; i <= order; i = i + 1 {
    if i <= f.order {
      acc = acc.add(pw.scale(f.coeffs[i]))
    }
    pw = pw.mul(gg)
  }
  acc
}

test "M1 J_k is a commutative local F_p-algebra with nilpotent maximal ideal" {
  for p in [2, 3] {
    let order = 2
    let els = all_jets(p, order)
    let zero = jet_zero(p, order)
    let one = jet_const(p, order, 1)
    let t = jet_t(p, order)
    // t^{order+1} = 0
    let mut pw = one
    for _i = 0; _i <= order; _i = _i + 1 {
      pw = pw.mul(t)
    }
    check(pw.is_zero(), "t^(k+1) = 0")
    for i = 0; i < els.length(); i = i + 1 {
      let a = els[i]
      check(a.add(zero).equal(a), "additive identity")
      check(a.mul(one).equal(a), "multiplicative identity")
      check(a.add(a.neg()).equal(zero), "additive inverse")
      // residue map J_k -> F_p, a |-> a_0, is a ring homomorphism
      check(a.is_unit() == (a.coeffs[0] != 0), "units are exactly a_0 != 0")
      if a.is_unit() {
        check(a.mul(a.inv()).equal(one), "inverse")
      } else {
        // every non-unit is nilpotent: a^{order+1} = 0
        let mut q = one
        for _j = 0; _j <= order; _j = _j + 1 {
          q = q.mul(a)
        }
        check(q.is_zero(), "non-units are nilpotent")
      }
      for j = 0; j < els.length(); j = j + 1 {
        let b = els[j]
        check(a.add(b).equal(b.add(a)), "commutative addition")
        check(a.mul(b).equal(b.mul(a)), "commutative multiplication")
        check(
          a.mul(b).coeffs[0] == (a.coeffs[0] * b.coeffs[0]) % p,
          "residue map multiplicative",
        )
        check(
          a.add(b).coeffs[0] == (a.coeffs[0] + b.coeffs[0]) % p,
          "residue map additive",
        )
        for t2 = 0; t2 < els.length(); t2 = t2 + 1 {
          let c = els[t2]
          check(a.add(b).add(c).equal(a.add(b.add(c))), "associative addition")
          check(a.mul(b).mul(c).equal(a.mul(b.mul(c))), "associative multiplication")
          check(a.mul(b.add(c)).equal(a.mul(b).add(a.mul(c))), "distributivity")
        }
      }
    }
  }
}

test "M2 Kock-Lawvereux and the Taylor expansion at every point" {
  let p = 3
  // dual numbers: order 1
  let duals = all_jets(p, 1)
  for i = 0; i < duals.length(); i = i + 1 {
    let f = duals[i]
    for a = 0; a < p; a = a + 1 {
      let g = jet_of_coeffs(p, 1, [a, 1])
      let lhs = f.compose(g)
      let rhs = jet_of_coeffs(p, 1, [f.eval_at(a), f.deriv().eval_at(a)])
      check(lhs.equal(rhs), "M2 f(a+eps) = f(a) + f'(a) eps")
    }
  }
  // higher order: f(a+t) = sum_i D^{(i)}f(a) t^i
  let order = 3
  let els = all_jets(p, order)
  for i = 0; i < els.length(); i = i + 1 {
    let f = els[i]
    for a = 0; a < p; a = a + 1 {
      let g = jet_of_coeffs(p, order, [a, 1])
      let lhs = f.compose(g)
      let cs : Array[Int] = []
      for k = 0; k <= order; k = k + 1 {
        cs.push(f.hasse(k).eval_at(a))
      }
      let rhs = jet_of_coeffs(p, order, cs)
      check(lhs.equal(rhs), "M2 Taylor expansion at a, f=" + f.to_string())
    }
    // hasse(1) is the ordinary derivative
    check(f.hasse(1).equal(f.deriv()), "D^(1) = d/dt")
  }
}

test "M3 derivations: Leibniz, linearity, [d/dt, E] = d/dt, Jacobi identity" {
  let p = 3
  let order = 4
  let els = all_jets(p, order)
  let t = jet_t(p, order)
  let t2 = t.mul(t)
  for i = 0; i < els.length(); i = i + 1 {
    let f = els[i]
    // exact identities (no products involved)
    check(
      f.deriv().euler().sub(f.euler().deriv()).equal(f.deriv().neg()),
      "M3 [d/dt,E] = d/dt",
    )
    check(f.hasse(1).equal(f.deriv()), "D^(1) = d/dt")
    // Jacobi identity for D1 = d/dt, D2 = E, D3 = t^2 d/dt, expanded fully,
    // compared two degrees below the truncation
    let d1f = f.deriv()
    let d2f = f.euler()
    let d3f = f.deriv().mul(t2)
    let a1 = d2f.deriv().mul(t2).sub(d3f.euler())
    let a1p = d1f.euler().deriv().mul(t2).sub(d1f.deriv().mul(t2).euler())
    let term1 = a1.deriv().sub(a1p)
    let a2 = d3f.deriv().sub(d1f.deriv().mul(t2))
    let a2p = d2f.deriv().mul(t2).deriv().sub(d2f.deriv().deriv().mul(t2))
    let term2 = a2.euler().sub(a2p)
    let a3 = d1f.euler().sub(d2f.deriv())
    let a3p = d3f.deriv().euler().sub(d3f.euler().deriv())
    let term3 = a3.deriv().mul(t2).sub(a3p)
    check(
      equal_upto(term1.add(term2).add(term3), jet_zero(p, order), order - 2),
      "M3 Jacobi identity, f=" + f.to_string(),
    )
    for j = 0; j < els.length(); j = j + 1 {
      let g = els[j]
      check(f.add(g).deriv().equal(f.deriv().add(g.deriv())), "M3 linearity")
      check(
        equal_upto(
          f.mul(g).deriv(),
          f.deriv().mul(g).add(f.mul(g.deriv())),
          order - 1,
        ),
        "M3 Leibniz d/dt",
      )
      check(
        equal_upto(
          f.mul(g).euler(),
          f.euler().mul(g).add(f.mul(g.euler())),
          order - 1,
        ),
        "M3 Leibniz Euler",
      )
      check(
        equal_upto(
          f.mul(g).deriv().mul(t2),
          f.deriv().mul(t2).mul(g).add(f.mul(g.deriv().mul(t2))),
          order - 1,
        ),
        "M3 Leibniz t^2 d/dt",
      )
    }
  }
}

test "M4 composition of jets: substitution, associativity, identity, chain rule, Faà di Bruno" {
  let p = 2
  let order = 3
  let els = all_jets(p, order)
  let ident = jet_t(p, order)
  // augmented jets (constant term 0): substitution is a monoid on these
  let aug : Array[Jet] = []
  for i = 0; i < els.length(); i = i + 1 {
    if els[i].coeffs[0] == 0 {
      aug.push(els[i])
    }
  }
  for i = 0; i < els.length(); i = i + 1 {
    let f = els[i]
    check(f.compose(ident).equal(f), "M4 f o t = f")
    for j = 0; j < els.length(); j = j + 1 {
      let g = els[j]
      check(
        f.compose(g).equal(compose_direct(f, g)),
        "M4 composition = substitution, f=" + f.to_string() + " g=" + g.to_string(),
      )
      check(
        equal_upto(
          f.compose(g).deriv(),
          f.deriv().compose(g).mul(g.deriv()),
          order - 1,
        ),
        "M4 chain rule",
      )
    }
    for j = 0; j < aug.length(); j = j + 1 {
      let g = aug[j]
      check(ident.compose(g).equal(g), "M4 t o g = g")
      for k = 0; k < aug.length(); k = k + 1 {
        let hh = aug[k]
        check(
          f.compose(g).compose(hh).equal(f.compose(g.compose(hh))),
          "M4 associativity of composition on augmented jets",
        )
      }
    }
  }
  // Faà di Bruno with partial Bell polynomials, for g with g_0 = 0
  let p2 = 3
  let ord = 3
  let fs = all_jets(p2, ord)
  let gs : Array[Jet] = []
  for i = 0; i < fs.length(); i = i + 1 {
    if fs[i].coeffs[0] == 0 {
      gs.push(fs[i])
    }
  }
  for i = 0; i < fs.length(); i = i + 1 {
    let f = fs[i]
    for j = 0; j < gs.length(); j = j + 1 {
      let g = gs[j]
      let lhs = f.compose(g)
      let cs : Array[Int] = []
      for n = 0; n <= ord; n = n + 1 {
        let mut acc = 0
        for m = 0; m <= n; m = m + 1 {
          if m > f.order {
            continue
          }
          acc = acc + f.coeffs[m] * bell_poly(g.coeffs, n, m, p2)
        }
        cs.push(acc)
      }
      let rhs = jet_of_coeffs(p2, ord, cs)
      check(lhs.equal(rhs), "M4 Faa di Bruno, f=" + f.to_string() + " g=" + g.to_string())
    }
  }
}

test "M5 heat flow on jets: exponential law, closed form, and the heat equation" {
  let p = 5
  let order = 4
  let els = all_jets(p, order)
  for i = 0; i < els.length(); i = i + 1 {
    let u0 = els[i]
    check(jet_heat_flow(u0, 0).equal(u0), "M5 u(0) = u0")
    // the Laplacian is the second derivative
    check(u0.laplacian().equal(u0.deriv().deriv()), "M5 Delta = d^2/dt^2")
    for s = 0; s < p; s = s + 1 {
      let us = jet_heat_flow(u0, s)
      // exponential law
      for r = 0; r < p; r = r + 1 {
        check(
          jet_heat_flow(u0, @kernel.modpos(s + r, p)).equal(jet_heat_flow(us, r)),
          "M5 exponential law, u0=" + u0.to_string(),
        )
      }
      // closed form: u(s)_i = sum_m (s^m/m!) * prod_{j=1..m}(i+2j)(i+2j-1) * a_{i+2m}
      for idx = 0; idx <= order; idx = idx + 1 {
        let mut want = 0
        let mut m = 0
        while idx + 2 * m <= order {
          let mut prod = 1
          for j = 1; j <= m; j = j + 1 {
            prod = prod * (idx + 2 * j) * (idx + 2 * j - 1)
          }
          let pw = pow_mod(s, m, p)
          let finv = factinv_mod(m, p)
          want = @kernel.modpos(want + pw * finv * prod * u0.coeffs[idx + 2 * m], p)
          m = m + 1
        }
        check(
          us.coeffs[idx] == want,
          "M5 closed form, u0=" + u0.to_string() + " i=" + idx.to_string(),
        )
      }
      // heat equation, with d/ds u extracted *independently* by interpolating
      // the polynomial r |-> u(s+r) and reading off its linear coefficient
      let ds : Array[Int] = []
      for idx = 0; idx <= order; idx = idx + 1 {
        let vals : Array[Int] = []
        for r = 0; r < p; r = r + 1 {
          vals.push(jet_heat_flow(u0, @kernel.modpos(s + r, p)).coeffs[idx])
        }
        ds.push(interpolate(p, vals)[1])
      }
      let du = jet_of_coeffs(p, order, ds)
      check(du.equal(us.laplacian()), "M5 heat equation d/ds u = Delta u, u0=" + u0.to_string())
    }
  }
}

///| Partial ordinary Bell polynomial B_{n,m}(g_1, g_2, ...) =
///| sum over compositions of n into m positive parts of the product of g's.
fn bell_poly(g : Array[Int], n : Int, m : Int, p : Int) -> Int {
  let comps = compositions_pos(n, m)
  let mut acc = 0
  for c = 0; c < comps.length(); c = c + 1 {
    let mut prod = 1
    for i = 0; i < comps[c].length(); i = i + 1 {
      let k = comps[c][i]
      prod = prod * (if k < g.length() {
        g[k]
      } else {
        0
      })
    }
    acc = acc + prod
  }
  @kernel.modpos(acc, p)
}

fn compositions_pos(n : Int, m : Int) -> Array[Array[Int]] {
  if m == 0 {
    return if n == 0 {
      [[]]
    } else {
      []
    }
  }
  if n < m {
    return []
  }
  let out : Array[Array[Int]] = []
  for first = 1; first <= n - m + 1; first = first + 1 {
    let sub = compositions_pos(n - first, m - 1)
    for i = 0; i < sub.length(); i = i + 1 {
      out.push(Array::add([first], sub[i]))
    }
  }
  out
}

fn factinv_mod(m : Int, p : Int) -> Int {
  let mut acc = 1
  for k = 1; k <= m; k = k + 1 {
    acc = @kernel.modpos(acc * @arith.invmod(k, p), p)
  }
  acc
}

fn pow_mod(b : Int, e : Int, p : Int) -> Int {
  let mut acc = 1
  for _i = 0; _i < e; _i = _i + 1 {
    acc = @kernel.modpos(acc * b, p)
  }
  acc
}

test "M6 every function F_p -> F_p is a polynomial of degree < p, uniquely" {
  for p in [2, 3] {
    let total = @arith.pow_int(p, p)
    for f = 0; f < total; f = f + 1 {
      let values = @arith.digits(f, p)
      let vs : Array[Int] = Array::make(p, 0)
      for i = 0; i < p; i = i + 1 {
        vs[i] = if i < values.length() {
          values[i]
        } else {
          0
        }
      }
      let coeffs = interpolate(p, vs)
      for a = 0; a < p; a = a + 1 {
        check(
          poly_eval_mod(coeffs, a, p) == vs[a],
          "M6 interpolation reproduces the function",
        )
      }
      // uniqueness: re-interpolating the evaluations gives the same polynomial
      let again : Array[Int] = []
      for a = 0; a < p; a = a + 1 {
        again.push(poly_eval_mod(coeffs, a, p))
      }
      let coeffs2 = interpolate(p, again)
      for i = 0; i < p; i = i + 1 {
        check(coeffs2[i] == coeffs[i], "M6 uniqueness of the interpolant")
      }
    }
    // the ring of functions has p^p elements and F_p[t]/(t^p - t) models it:
    // t^p - t vanishes on all of F_p
    let mut vanishes = true
    for a = 0; a < p; a = a + 1 {
      if @kernel.modpos(pow_mod(a, p, p) - a, p) != 0 {
        vanishes = false
      }
    }
    check(vanishes, "M6 t^p - t vanishes on F_p (Freshman's dream)")
  }
}

test "jet exp and log are mutually inverse on the unipotent/nilpotent jets" {
  let p = 5
  let order = 3
  let els = all_jets(p, order)
  for i = 0; i < els.length(); i = i + 1 {
    let f = els[i]
    if f.coeffs[0] == 0 {
      check(f.jet_exp().jet_log().equal(f), "log(exp(x)) = x")
      check(f.jet_exp().coeffs[0] == 1, "exp is unipotent")
    }
    if f.coeffs[0] == 1 {
      check(f.jet_log().jet_exp().equal(f), "exp(log(y)) = y")
    }
    if f.is_unit() {
      check(f.mul(f.inv()).equal(jet_const(p, order, 1)), "unit inverse")
    }
  }
}
