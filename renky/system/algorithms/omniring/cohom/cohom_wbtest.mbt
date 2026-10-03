// ===========================================================================
//  義務自證 · C1–C6 for cohomology rings.  The algebraic laws of the models are
//  checked on every basis triple, Poincaré duality by the rank of the pairing
//  matrix, and the de Rham DG algebra by complete enumeration of all forms.
// ===========================================================================

fn check(cond : Bool, msg : String) -> Unit {
  if cond {
    ()
  } else {
    println("CHECK FAILED: " + msg)
    abort("[cohom law failed] " + msg)
  }
}

test "C1/C2/C3 the cohomology rings of S^n, T^2, CP^2, RP^2, S^2 x S^2" {
  for p in [2, 3, 5] {
    let models : Array[CohomRing] = [
      sphere(p, 1),
      sphere(p, 2),
      sphere(p, 3),
      torus2(p),
      cp2(p),
      s2xs2(p),
    ]
    for i = 0; i < models.length(); i = i + 1 {
      let r = models[i]
      check(r.degree_ok(), "C1 degrees are additive, " + r.to_string())
      check(r.associativity_ok(), "C1 the cup product is associative, " + r.to_string())
      check(
        r.graded_commutative_ok(),
        "C1 graded commutativity, " + r.to_string(),
      )
      // the unit acts as the identity on every basis element
      for k = 0; k < r.cup.length(); k = k + 1 {
        check(r.cup[0][k][k] == 1, "C1 unit on the left, " + r.to_string())
        check(r.cup[k][0][k] == 1, "C1 unit on the right, " + r.to_string())
      }
      // products land in degree <= top
      for a = 0; a < r.cup.length(); a = a + 1 {
        for b = 0; b < r.cup.length(); b = b + 1 {
          for c = 0; c < r.cup.length(); c = c + 1 {
            if r.cup[a][b][c] != 0 {
              check(r.deg(a) + r.deg(b) <= r.top, "C1 products stay below the top degree")
            }
          }
        }
      }
      check(
        r.poincare_duality_ok(),
        "C2 Poincare duality, " + r.to_string(),
      )
      check(r.dims[r.top] == 1, "C2 one-dimensional top degree")
      check(r.dims[0] == 1, "C2 one-dimensional H^0")
      for k = 0; k <= r.top; k = k + 1 {
        check(r.dims[k] == r.dims[r.top - k], "C2 symmetric Betti numbers")
      }
    }
    // C3 Euler characteristics
    check(sphere(p, 1).euler() == 0, "C3 chi(S^1) = 0")
    check(sphere(p, 2).euler() == 2, "C3 chi(S^2) = 2")
    check(sphere(p, 3).euler() == 0, "C3 chi(S^3) = 0")
    check(torus2(p).euler() == 0, "C3 chi(T^2) = 0")
    check(cp2(p).euler() == 3, "C3 chi(CP^2) = 3")
    check(s2xs2(p).euler() == 4, "C3 chi(S^2 x S^2) = 4")
    // the Betti numbers of T^2 are the binomial coefficients
    check(torus2(p).dims[0] == 1 && torus2(p).dims[1] == 2 && torus2(p).dims[2] == 1, "C3 H^*(T^2)")
    // the fundamental class of T^2 is a cup b, and a cup a = 0
    let t2 = torus2(p)
    check(t2.cup[1][1][3] == 0, "C3 odd classes square to zero")
    check(t2.cup[1][2][3] == 1, "C3 a cup b is the fundamental class")
    check(
      t2.cup[2][1][3] == @kernel.modpos(-1, p),
      "C3 b cup a = -(a cup b)",
    )
    // CP^2: h^2 is the fundamental class
    let c2 = cp2(p)
    check(c2.cup[1][1][2] == 1, "C3 h^2 = [CP^2]")
  }
  // RP^2 with F_2 coefficients: x^2 is the fundamental class and signs vanish
  let rp = rp2(2)
  check(rp.associativity_ok() && rp.degree_ok(), "C1 RP^2 laws")
  check(rp.graded_commutative_ok(), "C1 RP^2 is commutative over F_2")
  check(rp.poincare_duality_ok(), "C2 Poincare duality for RP^2 over F_2")
  check(rp.euler() == 1, "C3 chi(RP^2) = 1")
  check(rp.cup[1][1][2] == 1, "C3 x^2 = [RP^2]")
}

test "C4 the de Rham DG algebra of the jet algebra: d^2 = 0 and graded Leibniz" {
  // d/dt descends to the jet algebra F_p[t]/(t^n) exactly when p | n, since
  // d(t^n) = n t^{n-1} must vanish in the quotient; only those shapes are DG
  // algebras and only those are tested here (the criterion itself is checked
  // in C6).
  let shapes = [(2, 2), (2, 4), (3, 3)]
  for sh = 0; sh < shapes.length(); sh = sh + 1 {
    let p = shapes[sh].0
    let n = shapes[sh].1
    let forms = all_forms(p, n)
    let hom = homogeneous_forms(p, n)
    let dt = dr_form(p, n, [], [1])
    let one = dr_form(p, n, [1], [])
    check(dt.wedge(dt).is_zero(), "C4 dt^dt = 0")
    for i = 0; i < forms.length(); i = i + 1 {
      let a = forms[i]
      check(a.deriv().deriv().is_zero(), "C4 d^2 = 0")
      check(a.wedge(a.deriv()).deg() <= 1, "C4 forms have degree <= 1")
      check(one.wedge(a).sub(a).is_zero() && a.wedge(one).sub(a).is_zero(), "C4 unit")
    }
    // the graded laws need homogeneous elements (a degree is only defined there)
    for i = 0; i < hom.length(); i = i + 1 {
      let a = hom[i]
      for j = 0; j < hom.length(); j = j + 1 {
        let b = hom[j]
        let sign = if (a.deg() * b.deg()) % 2 == 0 {
          1
        } else {
          -1
        }
        check(
          a.wedge(b).sub(b.wedge(a).wedge(dr_form(p, n, [sign], []))).is_zero(),
          "C4 graded commutativity of the wedge",
        )
        let lhs = a.wedge(b).deriv()
        let sgn0 = if a.deg() % 2 == 0 {
          1
        } else {
          p - 1
        }
        let rhs = a.deriv().wedge(b).add(a.wedge(b.deriv()).wedge(dr_form(p, n, [sgn0], [])))
        check(
          lhs.sub(rhs).is_zero(),
          "C4 graded Leibniz, p=" + p.to_string() + " n=" + n.to_string() + " a=" + a.to_string() +
            " b=" + b.to_string() + " lhs=" + lhs.to_string() + " rhs=" + rhs.to_string(),
        )
        for k = 0; k < hom.length(); k = k + 1 {
          let c = hom[k]
          check(
            a.wedge(b).wedge(c).sub(a.wedge(b.wedge(c))).is_zero(),
            "C4 associativity of the wedge",
          )
          check(
            a.wedge(b.add(c)).sub(a.wedge(b).add(a.wedge(c))).is_zero(),
            "C4 bilinearity of the wedge",
          )
        }
      }
    }
  }
}

///| The homogeneous forms: the 0-forms (g = 0) and the 1-forms (f = 0).
fn homogeneous_forms(p : Int, n : Int) -> Array[DRForm] {
  let out : Array[DRForm] = []
  let total = @arith.pow_int(p, n)
  for k = 0; k < total; k = k + 1 {
    let ds = @arith.digits(k, p)
    out.push(dr_form(p, n, ds, []))
    out.push(dr_form(p, n, [], ds))
  }
  out
}

///| The image of d, enumerated once.
fn exact_set(p : Int, n : Int) -> Array[DRForm] {
  let forms = all_forms(p, n)
  let out : Array[DRForm] = []
  for i = 0; i < forms.length(); i = i + 1 {
    out.push(forms[i].deriv())
  }
  out
}

///| Is the form exact, i.e. in the image of d?
fn is_exact_in(w : DRForm, xs : Array[DRForm]) -> Bool {
  for i = 0; i < xs.length(); i = i + 1 {
    if xs[i].sub(w).is_zero() {
      return true
    }
  }
  false
}

test "C5 cohomology of a DG algebra is a graded ring" {
  let shapes = [(2, 2), (3, 3)]
  for sh = 0; sh < shapes.length(); sh = sh + 1 {
    let p = shapes[sh].0
    let n = shapes[sh].1
    let hom = homogeneous_forms(p, n)
    let xs = exact_set(p, n)
    let closed : Array[DRForm] = []
    for i = 0; i < hom.length(); i = i + 1 {
      if hom[i].deriv().is_zero() {
        closed.push(hom[i])
      }
    }
    let exact : Array[DRForm] = []
    for i = 0; i < hom.length(); i = i + 1 {
      exact.push(hom[i].deriv())
    }
    for i = 0; i < closed.length(); i = i + 1 {
      for j = 0; j < closed.length(); j = j + 1 {
        check(
          closed[i].wedge(closed[j]).deriv().is_zero(),
          "C5 cocycles are closed under the product",
        )
        // the product on cohomology is graded commutative
        let a = closed[i]
        let b = closed[j]
        let sign = if (a.deg() * b.deg()) % 2 == 0 {
          1
        } else {
          -1
        }
        let diff = a.wedge(b).sub(b.wedge(a).wedge(dr_form(p, n, [sign], [])))
        check(is_exact_in(diff, xs), "C5 H^* is graded commutative")
      }
      for j = 0; j < exact.length(); j = j + 1 {
        check(is_exact_in(exact[j].wedge(closed[i]), xs), "C5 coboundary * cocycle is exact")
        check(is_exact_in(closed[i].wedge(exact[j]), xs), "C5 cocycle * coboundary is exact")
      }
    }
    // H^0 and H^1 of the jet algebra are one dimensional when n = p
    let (h0, h1) = dr_cohomology_dims(p, n)
    check(
      h0 == dr_cohomology_count(p, n).0 && h1 == dr_cohomology_count(p, n).1,
      "C5 dimensions agree with the count",
    )
    check(h0 == 1, "C5 H^0 = F_p (only constants are closed), n = p")
    let one = dr_form(p, n, [1], [])
    let dt = dr_form(p, n, [], [1])
    check(one.wedge(dt).sub(dt).is_zero(), "C5 [1] cup [dt] = [dt]")
    check(dt.wedge(dt).is_zero(), "C5 [dt]^2 = 0")
    // dt = d(t) is exact, but t^{n-1} dt generates H^1 and is not exact
    check(is_exact_in(dt, xs), "C5 dt = d(t) is exact")
    let top = dr_form(p, n, [], Array::add(Array::make(n - 1, 0), [1]))
    check(!(is_exact_in(top, xs)), "C5 t^(n-1) dt generates H^1 and is not exact")
    check(top.wedge(top).is_zero(), "C5 the top class squares to zero")
    check(one.wedge(top).sub(top).is_zero(), "C5 [1] cup [t^(n-1) dt] = [t^(n-1) dt]")
  }
}

test "C6 de Rham dimensions: linear algebra agrees with counting, and chi(Omega) = chi(H)" {
  for p in [2, 3, 5] {
    for n = 1; n <= 6; n = n + 1 {
      let by_linalg = dr_cohomology_dims(p, n)
      let by_count = dr_cohomology_count(p, n)
      check(
        by_linalg.0 == by_count.0 && by_linalg.1 == by_count.1,
        "C6 dimensions, p=" + p.to_string() + " n=" + n.to_string(),
      )
      // chi(Omega) = n - n = 0 = chi(H) = h0 - h1
      check(by_linalg.0 - by_linalg.1 == 0, "C6 chi(Omega) = chi(H) = 0")
      // d/dt descends to F_p[t]/(t^n) exactly when p | n
      check(
        dr_well_defined(p, n) == (n % p == 0),
        "C6 the de Rham differential is well defined iff p | n",
      )
      // independent count by enumerating closed and exact forms (small cases)
      if n % p == 0 && @arith.pow_int(p, 2 * n) <= 1024 {
        let forms = all_forms(p, n)
        let xs = exact_set(p, n)
        let mut closed0 = 0
        let mut exact1 = 0
        for i = 0; i < forms.length(); i = i + 1 {
          let w = forms[i]
          // a 0-form is one with g = 0
          let mut is0 = true
          for k = 0; k < n; k = k + 1 {
            if w.g[k] != 0 {
              is0 = false
            }
          }
          if is0 && w.deriv().is_zero() {
            closed0 = closed0 + 1
          }
          let mut is1 = false
          for k = 0; k < n; k = k + 1 {
            if w.g[k] != 0 {
              is1 = true
            }
          }
          if is1 && is_exact_in(w, xs) {
            exact1 = exact1 + 1
          }
        }
        check(
          closed0 == @arith.pow_int(p, by_linalg.0),
          "C6 #closed 0-forms = p^h0, p=" + p.to_string() + " n=" + n.to_string(),
        )
        check(
          exact1 == @arith.pow_int(p, n - by_linalg.1) - 1,
          "C6 #nonzero exact 1-forms = p^rank(d) - 1",
        )
      }
    }
  }
}
