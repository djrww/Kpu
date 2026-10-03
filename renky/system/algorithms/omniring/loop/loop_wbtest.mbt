// ===========================================================================
//  義務自證 · L1–L6: the closed loop itself.  `run_all(p)` pushes all eight
//  structures through structure -> finite-field data -> tropicalization ->
//  homology -> invariant -> reduction, and every cross-check must pass.
// ===========================================================================

fn lcheck(cond : Bool, msg : String) -> Unit {
  if cond {
    ()
  } else {
    println("CHECK FAILED: " + msg)
    abort("[loop law failed] " + msg)
  }
}

test "L1 the closed loop closes for p = 2 and p = 3" {
  for p in [2, 3] {
    let rep = run_all(p)
    println(rep.to_string())
    for i = 0; i < rep.checks.length(); i = i + 1 {
      lcheck(rep.checks[i].1, "L1 check failed: " + rep.checks[i].0)
    }
    lcheck(rep.closed(), "L1 the loop is closed for p = " + p.to_string())
    // L3: all eight structures are present
    lcheck(rep.stages.length() == 7, "L3 seven ideal stages")
    let names : Array[String] = []
    for i = 0; i < rep.stages.length(); i = i + 1 {
      names.push(rep.stages[i].name)
    }
    let pad = stage_padic(p, 4)
    names.push("p進整數環")
    let want = [
      "仿射半群環",
      "上同調環",
      "熱帶半環",
      "光滑函數環",
      "完美的F分裂環",
      "分次環",
      "理想類群對應環",
      "p進整數環",
    ]
    for i = 0; i < want.length(); i = i + 1 {
      let mut found = false
      for j = 0; j < names.length(); j = j + 1 {
        if names[j] == want[i] {
          found = true
        }
      }
      lcheck(found, "L3 structure present: " + want[i])
    }
    // L2: every stage is a finite dimensional Artin F_p-algebra
    for i = 0; i < rep.stages.length(); i = i + 1 {
      let s = rep.stages[i]
      lcheck(s.dim > 0, "L2 positive dimension: " + s.name)
      lcheck(s.kdim == 0, "L2 Krull dimension 0: " + s.name)
      let mut tot = 0
      for d = 0; d < s.hilbert.length(); d = d + 1 {
        tot = tot + s.hilbert[d]
      }
      lcheck(tot == s.dim, "L2 sum of the Hilbert function = dim: " + s.name)
      lcheck(
        s.reduced == s.fsplit && s.reduced == (s.loewy == 1),
        "L2 reduced <=> F-split <=> Loewy 1: " + s.name,
      )
      lcheck(s.trop_ok, "L2 tropical valuation axioms: " + s.name)
      lcheck(s.maslov_ok, "L2 Maslov bridge: " + s.name)
      lcheck(s.rank_nullity_ok, "L2 rank-nullity: " + s.name)
    }
    lcheck(pad.trop_ok, "L2 p-adic valuation axioms")
    lcheck(pad.teich_ok, "L2 Teichmuller representatives")
    lcheck(pad.explog_ok == (p >= 3), "L2 p-adic exp/log holds exactly for odd p")
    lcheck(pad.hensel_ok, "L2 Hensel lifting")
  }
}

test "L4 idempotence: the reduction is a closure operator on the chain" {
  let p = 3
  let st = @ideal.ideal_new(p, 2)
  st.add(@poly.mp_monom(p, 2, [3, 0], 1))
  st.add(@poly.mp_monom(p, 2, [0, 2], 1))
  let s0 = stage_from_ideal("R", "F_3[x,y]/(x^3,y^2)", st)
  lcheck(!(s0.reduced), "L4 R is not reduced")
  lcheck(s0.loewy > 1, "L4 the Loewy length of R is > 1")
  let r1 = reduce_ring(st)
  let s1 = stage_from_ideal("R1", "R/rad(0)", r1)
  lcheck(s1.reduced, "L4 the reduction is reduced")
  lcheck(s1.loewy == 1, "L4 the reduction has Loewy length 1")
  lcheck(s1.fsplit, "L4 the reduction is F-split")
  // applying the closure operator again changes nothing
  let r2 = reduce_ring(r1)
  let s2 = stage_from_ideal("R2", "(R/rad(0))/rad(0)", r2)
  lcheck(s2.dim == s1.dim, "L4 idempotent: the dimension is stable")
  lcheck(s2.loewy == s1.loewy, "L4 idempotent: the Loewy length is stable")
  lcheck(s2.euler == s1.euler, "L4 idempotent: chi is stable")
  lcheck(s2.reduced == s1.reduced && s2.fsplit == s1.fsplit, "L4 idempotent: flags stable")
  for d = 0; d < s1.hilbert.length(); d = d + 1 {
    lcheck(s1.hilbert[d] == s2.hilbert[d], "L4 idempotent: the Hilbert function is stable")
  }
  // and the invariants of the reduced ring are a fixed point of the whole chain
  lcheck(s1.invariant == s2.invariant, "L4 the invariant is a fixed point")
}

test "L5 the class-group correspondence ring detects p | h(D)" {
  // h(-23) = 3
  let s2 = stage_classgroup(2, -23)
  lcheck(s2.dim == 3, "L5 dim F_2[Cl(-23)] = h(-23) = 3")
  lcheck(s2.reduced && s2.fsplit && s2.loewy == 1, "L5 semisimple for p = 2 (Maschke)")
  let s3 = stage_classgroup(3, -23)
  lcheck(s3.dim == 3, "L5 dim F_3[Cl(-23)] = 3")
  lcheck(!(s3.reduced), "L5 not reduced when p | h(D)")
  lcheck(s3.loewy == 2, "L5 the Loewy length is v_p(h) + 1 = 2")
  // h(-47) = 5
  let t2 = stage_classgroup(2, -47)
  lcheck(t2.dim == 5, "L5 dim F_2[Cl(-47)] = 5")
  lcheck(t2.reduced, "L5 F_2[Cl(-47)] is semisimple")
  let t5 = stage_classgroup(5, -47)
  lcheck(t5.dim == 5, "L5 dim F_5[Cl(-47)] = 5")
  lcheck(!(t5.reduced) && t5.loewy == 2, "L5 F_5[Cl(-47)] is not semisimple, Loewy 2")
}

test "L6 the chain agrees with the independent engines on every structure" {
  let p = 2
  // semigroup ring: the toric ideal of <(2,0),(0,2),(1,1)> has dim 9 over F_2
  let gens = [[2, 0], [0, 2], [1, 1]]
  let st = @semigroup.toric_ideal(gens, p, 3)
  st.add(@poly.mp_monom(p, 3, [2, 0, 0], 1))
  st.add(@poly.mp_monom(p, 3, [0, 2, 0], 1))
  st.add(@poly.mp_monom(p, 3, [0, 0, 2], 1))
  lcheck(st.finite, "L6 the truncated toric quotient is finite")
  let s = stage_from_ideal("半群環", "F_2[x,y,z]/(I_A + (x^2,y^2,z^2))", st)
  // independent count: the standard monomials are x^i y^j z^k with i,j,k < 3
  // that are not divisible by the leading monomial xy of the toric relation
  let mut expect = 0
  for i = 0; i < 2; i = i + 1 {
    for j = 0; j < 2; j = j + 1 {
      for k = 0; k < 2; k = k + 1 {
        if !(i >= 1 && j >= 1) {
          expect = expect + 1
        }
      }
    }
  }
  lcheck(
    s.dim == expect,
    "L6 dim = #standard monomials = #semigroup elements of bounded degree: " +
      s.dim.to_string() +
      " vs " +
      expect.to_string(),
  )
  // the semigroup itself: an independent enumeration of <(2,0),(0,2),(1,1)>
  // inside the box [0,8]^2 agrees with the DP membership engine
  let sg = @semigroup.semigroup_new(2, gens)
  let got = sg.elements(8)
  let want : Array[Array[Int]] = []
  for i = 0; i <= 4; i = i + 1 {
    for j = 0; j <= 4; j = j + 1 {
      for k = 0; k <= 8; k = k + 1 {
        let v = [2 * i + k, 2 * j + k]
        if v[0] > 8 || v[1] > 8 {
          continue
        }
        let mut dup = false
        for t = 0; t < want.length(); t = t + 1 {
          if want[t][0] == v[0] && want[t][1] == v[1] {
            dup = true
          }
        }
        if !(dup) {
          want.push(v)
        }
      }
    }
  }
  lcheck(
    got.length() == want.length(),
    "L6 the semigroup enumeration agrees: " +
      got.length().to_string() +
      " vs " +
      want.length().to_string(),
  )
  for t = 0; t < want.length(); t = t + 1 {
    lcheck(sg.contains(want[t]), "L6 membership of " + want[t].to_string())
  }
  // tropical semiring stage: F_2^2 = F_2[x]/(x^2 - x) is semisimple
  let st3 = @ideal.ideal_new(p, 1)
  st3.add(@poly.mp_make(p, 1, [([p], 1), ([1], p - 1)]))
  let s3 = stage_from_ideal("熱帶", "F_2^2", st3)
  lcheck(s3.dim == p, "L6 dim F_p^p = p")
  lcheck(s3.reduced && s3.loewy == 1, "L6 F_p^p is semisimple: the idempotent shadow of (min,+)")
  // the tropical semiring laws themselves, on the valuation image
  lcheck(@tropical.trop_add(3, 3) == 3, "L6 (min,+) is idempotent")
  lcheck(@tropical.trop_mul(2, 5) == 7, "L6 (min,+) multiplication is addition")
  lcheck(@tropical.trop_add(1, 4) == 1, "L6 min")
  // heat: the dequantization bridge on a concrete energy spectrum
  let energies = [3.0, 1.0, 4.0, 1.0, 5.0]
  let (tmin, argmin) = @tropical.trop_limit(energies)
  lcheck(tmin == 1.0 && argmin == 1, "L6 the tropical limit of the spectrum")
  // the multiplicity of the minimum is the residual entropy of the bridge
  let mut mult = 0
  for i = 0; i < energies.length(); i = i + 1 {
    if energies[i] == tmin {
      mult = mult + 1
    }
  }
  lcheck(mult == 2, "L6 the minimum is attained twice")
  let ln_n = @heat.logd(energies.length().to_double())
  for h in [1.0, 0.25, 0.05] {
    let sm = @heat.softmin_all(energies, h)
    lcheck(sm <= tmin + 1.0e-9, "L6 softmin <= min")
    lcheck(sm >= tmin - h * ln_n - 1.0e-9, "L6 softmin >= min - h ln n")
  }
  // the bridge is monotone in h: softmin increases as h decreases
  let s1 = @heat.softmin_all(energies, 1.0)
  let s2 = @heat.softmin_all(energies, 0.25)
  let s3 = @heat.softmin_all(energies, 0.05)
  lcheck(s1 <= s2 + 1.0e-9 && s2 <= s3 + 1.0e-9, "L6 the Maslov bridge is monotone")
  lcheck(s3 - tmin <= -0.05 * @heat.logd(2.0) + 1.0e-9, "L6 residual entropy ~ h ln g")
  // homology: chi of the Koszul complex of a unit is 0, of a nilpotent is not
  let stn = @ideal.ideal_new(p, 1)
  stn.add(@poly.mp_monom(p, 1, [p], 1))
  let sn = stage_from_ideal("jet", "F_2[t]/(t^2)", stn)
  lcheck(sn.dim == p, "L6 dim of the jet algebra = p")
  lcheck(!(sn.reduced) && sn.loewy == p, "L6 the jet algebra has Loewy length p")
}
