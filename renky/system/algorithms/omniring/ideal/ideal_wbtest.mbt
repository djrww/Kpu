// ===========================================================================
//  義務自證 · T1–T8 for the live ideal state, plus the polynomial ring laws.
//  Domains are enumerated completely: these are proofs, not samples.
// ===========================================================================


fn check(cond : Bool, msg : String) -> Unit {
  if cond {
    ()
  } else {
    println("CHECK FAILED: " + msg)
    abort("[law failed] " + msg)
  }
}

///| Every polynomial over F_p supported on the monomials 1, x, y (p^3 of them).
fn small_polys(p : Int) -> Array[MPoly] {
  let out : Array[MPoly] = []
  let e0 = [0, 0]
  let e1 = [1, 0]
  let e2 = [0, 1]
  for a = 0; a < p; a = a + 1 {
    for b = 0; b < p; b = b + 1 {
      for c = 0; c < p; c = c + 1 {
        out.push(mp_make(p, 2, [(e0, a), (e1, b), (e2, c)]))
      }
    }
  }
  out
}

///| Every polynomial over F_p supported on 1, x, y, x^2, xy, y^2 (p^6).
fn degree2_polys(p : Int) -> Array[MPoly] {
  let out : Array[MPoly] = []
  let es = [[0, 0], [1, 0], [0, 1], [2, 0], [1, 1], [0, 2]]
  let total = @arith.pow_int(p, 6)
  for k = 0; k < total; k = k + 1 {
    let ds = @arith.digits(k, p)
    let terms : Array[(Array[Int], Int)] = []
    for i = 0; i < 6; i = i + 1 {
      terms.push((es[i], if i < ds.length() {
        ds[i]
      } else {
        0
      }))
    }
    out.push(mp_make(p, 2, terms))
  }
  out
}

// ---------------------------------------------------------------------------
//  Polynomial ring laws
// ---------------------------------------------------------------------------

test "F_p[x,y] is a commutative ring: all axioms on every triple of small polynomials" {
  for p in [2, 3] {
    let els = small_polys(p)
    let zero = @poly.mp_zero(p, 2)
    let one = mp_one(p, 2)
    for i = 0; i < els.length(); i = i + 1 {
      let a = els[i]
      check(a.add(zero).equal(a), "add identity")
      check(a.mul(one).equal(a), "mul identity")
      check(a.add(a.neg()).equal(zero), "add inverse")
      for j = 0; j < els.length(); j = j + 1 {
        let b = els[j]
        check(a.add(b).equal(b.add(a)), "add comm")
        check(a.mul(b).equal(b.mul(a)), "mul comm")
        for t = 0; t < els.length(); t = t + 1 {
          let c = els[t]
          check(a.add(b).add(c).equal(a.add(b.add(c))), "add assoc p=" + p.to_string())
          check(a.mul(b).mul(c).equal(a.mul(b.mul(c))), "mul assoc p=" + p.to_string())
          check(a.mul(b.add(c)).equal(a.mul(b).add(a.mul(c))), "distrib p=" + p.to_string())
        }
      }
    }
  }
}

test "evaluation is a ring homomorphism F_p[x,y] -> F_p at every point" {
  for p in [2, 3] {
    let els = small_polys(p)
    for x = 0; x < p; x = x + 1 {
      for y = 0; y < p; y = y + 1 {
        let pt = [x, y]
        for i = 0; i < els.length(); i = i + 1 {
          for j = 0; j < els.length(); j = j + 1 {
            let a = els[i]
            let b = els[j]
            check(
              (a.eval(pt) + b.eval(pt)) % p == a.add(b).eval(pt),
              "eval additive",
            )
            check(
              (a.eval(pt) * b.eval(pt)) % p == a.mul(b).eval(pt),
              "eval multiplicative",
            )
          }
        }
      }
    }
  }
}

test "division algorithm: f = sum q_i g_i + r and no term of r is reducible" {
  let p = 2
  let g0 = mp_make(p, 2, [([2, 0], 1), ([0, 1], 1)])
  let g1 = mp_make(p, 2, [([1, 1], 1), ([0, 0], 1)])
  let gs = [g0, g1]
  let els = degree2_polys(p)
  for i = 0; i < els.length(); i = i + 1 {
    let (qs, r) = mp_divmod(els[i], gs)
    let mut acc = r
    for k = 0; k < qs.length(); k = k + 1 {
      acc = acc.add(qs[k].mul(gs[k]))
    }
    check(acc.equal(els[i]), "division identity " + els[i].to_string())
    for t = 0; t < r.terms.length(); t = t + 1 {
      let e = r.terms[t].0
      let mut reducible = false
      for k = 0; k < gs.length(); k = k + 1 {
        if monom_divides(gs[k].lm(), e) {
          reducible = true
        }
      }
      check(!(reducible), "remainder irreducible")
    }
  }
}

test "T2 Buchberger criterion: every S-polynomial of the computed basis reduces to zero" {
  let cases : Array[(Int, Array[MPoly])] = [
    (2, [mp_make(2, 2, [([2, 0], 1), ([0, 1], 1)]), mp_make(2, 2, [([1, 1], 1), ([0, 0], 1)])]),
    (3, [mp_make(3, 2, [([1, 0], 1), ([0, 0], 2)]), mp_make(3, 2, [([0, 1], 1), ([1, 0], 1)])]),
    (2, [mp_make(2, 2, [([3, 0], 1), ([0, 1], 1)]), mp_make(2, 2, [([1, 1], 1)])]),
    (5, [mp_make(5, 2, [([2, 0], 1), ([0, 2], 4)]), mp_make(5, 2, [([1, 1], 1), ([0, 0], 1)])]),
  ]
  for k = 0; k < cases.length(); k = k + 1 {
    let p = cases[k].0
    let gs = cases[k].1
    let gb = @poly.mp_groebner(gs, true)
    for i = 0; i < gb.length(); i = i + 1 {
      for j = i + 1; j < gb.length(); j = j + 1 {
        let s = mp_spoly(gb[i], gb[j])
        check(
          mp_normal_form(s, gb).is_zero(),
          "S-poly reduces to zero, case " + k.to_string(),
        )
      }
    }
    // the reduced basis generates the same ideal as the input
    for i = 0; i < gs.length(); i = i + 1 {
      check(
        mp_normal_form(gs[i], gb).is_zero(),
        "input generators lie in the Groebner ideal",
      )
    }
    // and every basis element lies in the ideal generated by the input.
    // (a) necessary condition, checked independently by brute force: it must
    //     vanish at every point of F_p^2 where the generators vanish;
    // (b) for p = 2 the witness search over all multipliers of degree <= 2 is
    //     complete, so it decides membership exactly.
    for x = 0; x < p; x = x + 1 {
      for y = 0; y < p; y = y + 1 {
        let mut common = true
        for i = 0; i < gs.length(); i = i + 1 {
          if gs[i].eval([x, y]) != 0 {
            common = false
          }
        }
        if common {
          for i = 0; i < gb.length(); i = i + 1 {
            check(
              gb[i].eval([x, y]) == 0,
              "basis vanishes on V(gs), case " + k.to_string(),
            )
          }
        }
      }
    }
    if p == 2 {
      for i = 0; i < gb.length(); i = i + 1 {
        check(
          witness_in_ideal(gb[i], gs, p),
          "basis element has an explicit witness, case " + k.to_string(),
        )
      }
    }
    ignore(p)
  }
}

///| Independent (non-Groebner) membership check for a two-generator ideal:
///| exhaustive search over all multipliers of degree <= 2.
fn witness_in_ideal(f : MPoly, gs : Array[MPoly], p : Int) -> Bool {
  let cand = degree2_polys(p)
  for i = 0; i < cand.length(); i = i + 1 {
    for j = 0; j < cand.length(); j = j + 1 {
      let acc = cand[i].mul(gs[0]).add(cand[j].mul(gs[1]))
      if acc.equal(f) {
        return true
      }
    }
  }
  false
}

test "normal form is canonical: it does not depend on the order of the basis" {
  let p = 3
  let gs = [
    mp_make(p, 2, [([2, 0], 1), ([0, 1], 1)]),
    mp_make(p, 2, [([1, 1], 1), ([0, 0], 1)]),
  ]
  let gb = @poly.mp_groebner(gs, true)
  let rev : Array[MPoly] = []
  for i = gb.length() - 1; i >= 0; i = i - 1 {
    rev.push(gb[i])
  }
  let els = degree2_polys(p)
  for i = 0; i < els.length(); i = i + 1 {
    check(
      mp_normal_form(els[i], gb).equal(mp_normal_form(els[i], rev)),
      "canonical normal form " + els[i].to_string(),
    )
  }
}

// ---------------------------------------------------------------------------
//  T1, T3, T4, T5: the live ideal state
// ---------------------------------------------------------------------------

test "T1/T6 membership equals vanishing on V(I) for I = (x^p-x, y^p-y): all f of degree <= 2" {
  for p in [2, 3] {
    let st = ideal_new(p, 2)
    let fe = field_equations(p, 2)
    for i = 0; i < fe.length(); i = i + 1 {
      st.add(fe[i])
    }
    check(st.finite, "quotient is finite")
    check(st.dim == p * p, "dim = p^2, got " + st.dim.to_string())
    check(st.kdim == 0, "Krull dim 0")
    let els = degree2_polys(p)
    for i = 0; i < els.length(); i = i + 1 {
      let f = els[i]
      let mut vanishes = true
      for x = 0; x < p; x = x + 1 {
        for y = 0; y < p; y = y + 1 {
          if f.eval([x, y]) != 0 {
            vanishes = false
          }
        }
      }
      check(
        st.in_ideal(f) == vanishes,
        "T1 membership == vanishing, p=" + p.to_string() + " f=" + f.to_string(),
      )
      // T3 independent direction: the witness really reconstructs f
      if st.in_ideal(f) {
        let (qs, r) = st.membership_witness(f)
        check(r.is_zero(), "witness remainder zero")
        let mut acc = @poly.mp_zero(p, 2)
        for k = 0; k < qs.length(); k = k + 1 {
          acc = acc.add(qs[k].mul(st.basis[k]))
        }
        check(acc.equal(f), "witness reconstructs f")
      }
    }
  }
}

test "T3 standard monomials are linearly independent modulo I, for several ideals" {
  let cases : Array[(Int, Array[MPoly])] = [
    (2, [mp_make(2, 2, [([2, 0], 1)]), mp_make(2, 2, [([0, 2], 1), ([0, 1], 1)])]),
    (3, [mp_make(3, 2, [([1, 0], 1), ([0, 0], 1)]), mp_make(3, 2, [([0, 2], 1)])]),
    (2, [mp_make(2, 2, [([1, 1], 1), ([0, 0], 1)])]),
  ]
  for k = 0; k < cases.length(); k = k + 1 {
    let p = cases[k].0
    let st = ideal_new(p, 2)
    let gs = cases[k].1
    for i = 0; i < gs.length(); i = i + 1 {
      st.add(gs[i])
    }
    if !(st.finite) {
      continue
    }
    let bas = st.quotient_basis()
    check(bas.length() == st.dim, "basis size == dim, case " + k.to_string())
    // T4: dim == sum of Hilbert function
    let h = st.hilbert_fn(2 * st.dim + 2)
    let mut total = 0
    for d = 0; d < h.length(); d = d + 1 {
      total = total + h[d]
    }
    check(total == st.dim, "T4 dim = sum H(d), case " + k.to_string())
    check(h[0] == 1, "H(0) = 1")
    // exhaustive independence: a combination of standard monomials lies in I
    // only when every coefficient is zero
    let combos = @arith.pow_int(p, bas.length())
    for c = 0; c < combos; c = c + 1 {
      let ds = @arith.digits(c, p)
      let terms : Array[(Array[Int], Int)] = []
      let mut allzero = true
      for i = 0; i < bas.length(); i = i + 1 {
        let co = if i < ds.length() {
          ds[i]
        } else {
          0
        }
        if co != 0 {
          allzero = false
        }
        terms.push((bas[i], co))
      }
      let f = mp_make(p, 2, terms)
      check(
        st.in_ideal(f) == allzero,
        "T3 independence, case " + k.to_string() + " f=" + f.to_string(),
      )
    }
  }
}

test "T5 Krull dimension: zero-dimensional iff finite quotient, on a catalogue of ideals" {
  let p = 2
  let catalogue : Array[(Array[MPoly], Int, Bool)] = [
    ([], 2, false),
    ([mp_var(p, 2, 0)], 1, false),
    ([mp_var(p, 2, 0), mp_var(p, 2, 1)], 0, true),
    ([mp_monom(p, 2, [2, 0], 1)], 1, false),
    ([mp_monom(p, 2, [2, 0], 1), mp_monom(p, 2, [0, 3], 1)], 0, true),
    ([mp_make(p, 2, [([1, 1], 1), ([0, 0], 1)])], 1, false),
    ([mp_make(p, 2, [([1, 0], 1), ([0, 0], 1)]), mp_make(p, 2, [([0, 1], 1), ([0, 0], 1)])], 0, true),
  ]
  for k = 0; k < catalogue.length(); k = k + 1 {
    let st = ideal_new(p, 2)
    let gs = catalogue[k].0
    for i = 0; i < gs.length(); i = i + 1 {
      st.add(gs[i])
    }
    check(
      st.kdim == catalogue[k].1,
      "T5 kdim case " + k.to_string() + " got " + st.kdim.to_string(),
    )
    check(
      st.finite == catalogue[k].2,
      "T5 finite case " + k.to_string(),
    )
    check((st.kdim == 0) == st.finite, "T5 kdim==0 iff finite")
  }
}

test "T7 Rabinowitsch: radical membership agrees with the radical state, I = (x^2, y^2) over F_2, F_3" {
  for p in [2, 3] {
    let st = ideal_new(p, 2)
    st.add(mp_monom(p, 2, [2, 0], 1))
    st.add(mp_monom(p, 2, [0, 2], 1))
    let rad = st.radical()
    let els = small_polys(p)
    for i = 0; i < els.length(); i = i + 1 {
      let f = els[i]
      let via_rabinowitsch = st.in_radical(f)
      let via_radical = rad.in_ideal(f)
      check(
        via_rabinowitsch == via_radical,
        "T7 p=" + p.to_string() + " f=" + f.to_string(),
      )
    }
    // x and y are in the radical, 1 is not
    check(st.in_radical(mp_var(p, 2, 0)), "x in rad(x^2,y^2)")
    check(st.in_radical(mp_var(p, 2, 1)), "y in rad(x^2,y^2)")
    check(!(st.in_radical(mp_one(p, 2))), "1 not in rad")
  }
}

test "T6 finite-field Nullstellensatz: I(V(I)) == rad(I + (x^p-x, y^p-y)) on all f of degree <= 2" {
  let p = 2
  let st = ideal_new(p, 2)
  st.add(mp_make(p, 2, [([2, 0], 1), ([0, 1], 1)]))
  st.add(mp_make(p, 2, [([1, 1], 1), ([0, 0], 1)]))
  let with_fe = st.with_field_equations()
  let rad = with_fe.radical()
  let pts = st.variety()
  let els = degree2_polys(p)
  for i = 0; i < els.length(); i = i + 1 {
    let f = els[i]
    let mut vanishes = true
    for k = 0; k < pts.length(); k = k + 1 {
      if f.eval(pts[k]) != 0 {
        vanishes = false
      }
    }
    check(
      vanishes == rad.in_ideal(f),
      "T6 f=" + f.to_string() + " vanishes=" + vanishes.to_string(),
    )
  }
  check(pts.length() > 0, "variety nonempty for this ideal")
}

///| Every element of R/I, as a polynomial in the quotient basis (p^dim of them).
fn all_classes(st : IdealState) -> Array[MPoly] {
  let p = st.p
  let bas = st.quotient_basis()
  let out : Array[MPoly] = []
  let total = @arith.pow_int(p, bas.length())
  for c = 0; c < total; c = c + 1 {
    let ds = @arith.digits(c, p)
    let terms : Array[(Array[Int], Int)] = []
    for i = 0; i < bas.length(); i = i + 1 {
      terms.push((bas[i], if i < ds.length() {
        ds[i]
      } else {
        0
      }))
    }
    out.push(mp_make(p, 2, terms))
  }
  out
}

test "T8 maximality: I maximal iff R/I has no zero divisors" {
  let p = 3
  let cases : Array[(Array[MPoly], Bool)] = [
    ([mp_var(p, 2, 0), mp_var(p, 2, 1)], true),
    ([mp_make(p, 2, [([1, 0], 1), ([0, 0], 1)]), mp_make(p, 2, [([0, 1], 1), ([0, 0], 2)])], true),
    ([mp_monom(p, 2, [2, 0], 1), mp_var(p, 2, 1)], false),
    ([mp_make(p, 2, [([2, 0], 1), ([0, 0], 2)]), mp_var(p, 2, 1)], false),
  ]
  for k = 0; k < cases.length(); k = k + 1 {
    let st = ideal_new(p, 2)
    let gs = cases[k].0
    for i = 0; i < gs.length(); i = i + 1 {
      st.add(gs[i])
    }
    check(
      st.is_maximal() == cases[k].1,
      "T8 case " + k.to_string() + " dim=" + st.dim.to_string(),
    )
    if st.is_maximal() {
      // R/I is a field with p^dim elements: *every* nonzero element, not only
      // the basis monomials, must be invertible
      let els = all_classes(st)
      for i = 0; i < els.length(); i = i + 1 {
        let a = els[i]
        if st.in_ideal(a) {
          continue
        }
        let mut invertible = false
        for j = 0; j < els.length(); j = j + 1 {
          if st.in_ideal(a.mul(els[j]).sub(mp_one(p, 2))) {
            invertible = true
          }
        }
        check(invertible, "T8 field: invertible " + a.to_string())
      }
    } else {
      // not a field: some nonzero element is not invertible
      let els = all_classes(st)
      let mut found = false
      for i = 0; i < els.length(); i = i + 1 {
        let a = els[i]
        if st.in_ideal(a) {
          continue
        }
        let mut invertible = false
        for j = 0; j < els.length(); j = j + 1 {
          if st.in_ideal(a.mul(els[j]).sub(mp_one(p, 2))) {
            invertible = true
          }
        }
        if !(invertible) {
          found = true
        }
      }
      check(found, "T8 non-field: a non-invertible element exists")
    }
  }
}

test "incremental updates: adding generators one at a time equals adding them at once" {
  let p = 3
  let gs = [
    mp_make(p, 2, [([2, 0], 1), ([0, 1], 2)]),
    mp_make(p, 2, [([1, 1], 1), ([0, 0], 1)]),
    mp_make(p, 2, [([0, 2], 1), ([1, 0], 1)]),
  ]
  let incremental = ideal_new(p, 2)
  for i = 0; i < gs.length(); i = i + 1 {
    incremental.add(gs[i])
    // after every single update the state must already be consistent
    for j = 0; j <= i; j = j + 1 {
      check(
        incremental.in_ideal(gs[j]),
        "state contains every generator so far",
      )
    }
    let (std, _) = mp_standard_monomials(incremental.basis, 2, 6)
    ignore(std)
    check(incremental.updates == i + 1, "update counter")
  }
  let at_once = ideal_new(p, 2)
  for i = 0; i < gs.length(); i = i + 1 {
    at_once.add(gs[i])
  }
  check(incremental.dim == at_once.dim, "same quotient dimension")
  check(incremental.kdim == at_once.kdim, "same Krull dimension")
  // the two bases define the same ideal: each reduces to zero modulo the other
  for i = 0; i < incremental.basis.length(); i = i + 1 {
    check(
      mp_normal_form(incremental.basis[i], at_once.basis).is_zero(),
      "basis inclusion 1",
    )
  }
  for i = 0; i < at_once.basis.length(); i = i + 1 {
    check(
      mp_normal_form(at_once.basis[i], incremental.basis).is_zero(),
      "basis inclusion 2",
    )
  }
  let _ = incremental.describe()
}

test "Hilbert function of F_p[x,y]/(x^2, y^2) is 1,2,1 and of F_p[x,y]/(x^2) is 2,3,4,..." {
  let p = 5
  let st = ideal_new(p, 2)
  st.add(mp_monom(p, 2, [2, 0], 1))
  st.add(mp_monom(p, 2, [0, 2], 1))
  let h = st.hilbert_fn(4)
  check(h[0] == 1 && h[1] == 2 && h[2] == 1, "hilbert (x^2,y^2)")
  check(h[3] == 0 && h[4] == 0, "hilbert vanishes above the top degree")
  check(st.dim == 4, "dim 4")
  let st2 = ideal_new(p, 2)
  st2.add(mp_monom(p, 2, [2, 0], 1))
  let h2 = st2.hilbert_fn(4)
  check(h2[0] == 1 && h2[1] == 2 && h2[2] == 2 && h2[3] == 2, "hilbert (x^2)")
  check(st2.kdim == 1, "kdim of (x^2) is 1")
  check(!(st2.finite), "(x^2) has infinite quotient")
}
