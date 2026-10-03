// ===========================================================================
//  義務自證 · S1–S6 for affine semigroup rings.  Membership and lattice
//  questions are cross-checked against complete brute-force searches.
// ===========================================================================

fn check(cond : Bool, msg : String) -> Unit {
  if cond {
    ()
  } else {
    println("CHECK FAILED: " + msg)
    abort("[semigroup law failed] " + msg)
  }
}

///| Brute force: v in N*gens with coefficients in [0,lim]?
fn brute_semigroup(gens : Array[Array[Int]], v : Array[Int], lim : Int) -> Bool {
  let m = gens.length()
  let n = v.length()
  for c0 = 0; c0 <= lim; c0 = c0 + 1 {
    for c1 = 0; c1 <= lim; c1 = c1 + 1 {
      for c2 = 0; c2 <= lim; c2 = c2 + 1 {
        if m == 2 && c2 > 0 {
          continue
        }
        let cs = [c0, c1, c2]
        let mut ok = true
        for j = 0; j < n; j = j + 1 {
          let mut s = 0
          for i = 0; i < m; i = i + 1 {
            s = s + cs[i] * gens[i][j]
          }
          if s != v[j] {
            ok = false
          }
        }
        if ok {
          return true
        }
      }
    }
  }
  false
}

///| Brute force: v in Z*gens with coefficients in [-lim,lim]?
fn brute_lattice(gens : Array[Array[Int]], v : Array[Int], lim : Int) -> Bool {
  let m = gens.length()
  let n = v.length()
  for c0 = -lim; c0 <= lim; c0 = c0 + 1 {
    for c1 = -lim; c1 <= lim; c1 = c1 + 1 {
      for c2 = -lim; c2 <= lim; c2 = c2 + 1 {
        if m == 2 && c2 != 0 {
          continue
        }
        let cs = [c0, c1, c2]
        let mut ok = true
        for j = 0; j < n; j = j + 1 {
          let mut s = 0
          for i = 0; i < m; i = i + 1 {
            s = s + cs[i] * gens[i][j]
          }
          if s != v[j] {
            ok = false
          }
        }
        if ok {
          return true
        }
      }
    }
  }
  false
}

test "S1/S2 semigroup and lattice membership agree with exhaustive search" {
  let cases : Array[Array[Array[Int]]] = [
    [[2, 0], [0, 2], [1, 1]],
    [[1, 0], [0, 1]],
    [[3, 0], [0, 3], [1, 1]],
    [[2, 3], [3, 2], [1, 1]],
  ]
  for c = 0; c < cases.length(); c = c + 1 {
    let s = semigroup_new(2, cases[c])
    for x = 0; x <= 6; x = x + 1 {
      for y = 0; y <= 6; y = y + 1 {
        let v = [x, y]
        check(
          s.contains(v) == brute_semigroup(cases[c], v, 6),
          "S1 membership, case " + c.to_string() + " v=" + v.to_string(),
        )
        check(
          s.in_lattice(v) == brute_lattice(cases[c], v, 6),
          "S2 lattice membership, case " + c.to_string() + " v=" + v.to_string(),
        )
        // S subset L
        if s.contains(v) {
          check(s.in_lattice(v), "S is contained in its group completion")
        }
      }
    }
    // the enumeration inside a box matches pointwise membership
    let els = s.elements(4)
    for x = 0; x <= 4; x = x + 1 {
      for y = 0; y <= 4; y = y + 1 {
        let mut listed = false
        for i = 0; i < els.length(); i = i + 1 {
          if els[i][0] == x && els[i][1] == y {
            listed = true
          }
        }
        check(listed == s.contains([x, y]), "S1 enumeration agrees with membership")
      }
    }
  }
}

test "S3 normality: the classical examples" {
  let normal = semigroup_new(2, [[1, 0], [0, 1]])
  check(normal.is_normal(6, 4), "N^2 is normal")
  let normal2 = semigroup_new(2, [[1, 0], [1, 1], [0, 1]])
  check(normal2.is_normal(6, 4), "the cone of the unit square is normal")
  let bad = semigroup_new(2, [[2, 0], [0, 2], [1, 1]])
  check(!(bad.is_normal(6, 4)), "<(2,0),(0,2),(1,1)> is not normal")
  check(bad.contains([2, 0]), "(2,0) in S")
  check(!(bad.contains([1, 0])), "(1,0) not in S")
  check(bad.contains([2, 0]) && !(bad.contains([1, 0])), "witness of non-normality")
  let bad3 = semigroup_new(2, [[3, 0], [0, 3], [1, 1]])
  check(!(bad3.is_normal(6, 6)), "<(3,0),(0,3),(1,1)> is not normal")
}

test "S4 Hilbert function of k[S] and the standard-monomial/semigroup bijection" {
  let gens = [[2, 0], [0, 2], [1, 1]]
  let s = semigroup_new(2, gens)
  let w = [1, 1]
  let h = s.hilbert(w, 4, 8)
  // independent count: elements of S with <w,e> = d
  for d = 0; d <= 8; d = d + 1 {
    let mut cnt = 0
    for x = 0; x <= 4; x = x + 1 {
      for y = 0; y <= 4; y = y + 1 {
        if x + y == d && s.contains([x, y]) {
          cnt = cnt + 1
        }
      }
    }
    check(h[d] == cnt, "S4 Hilbert function at d=" + d.to_string())
  }
  check(h[0] == 1 && h[1] == 0 && h[2] == 3, "S4 hand values 1,0,3")
  // the quotient of F_p[x_1..x_m] by the toric ideal has standard monomials in
  // bijection with the semigroup elements reachable with |u| <= D
  let p = 3
  let dd = 3
  let st = toric_ideal(gens, p, dd)
  let (std, _) = @poly.mp_standard_monomials(st.basis, gens.length(), dd)
  let mut total = 0
  let exps : Array[Array[Int]] = []
  for dg = 0; dg <= dd; dg = dg + 1 {
    total = total + binom3(dg + 2)
  }
  ignore(total)
  for dg = 0; dg <= dd; dg = dg + 1 {
    let cs = compositions_pub(dg, gens.length())
    for i = 0; i < cs.length(); i = i + 1 {
      exps.push(cs[i])
    }
  }
  let images : Array[Array[Int]] = []
  for i = 0; i < exps.length(); i = i + 1 {
    let img = monomial_map_exponent(gens, exps[i])
    let mut fresh = true
    for j = 0; j < images.length(); j = j + 1 {
      if vec_eq(images[j], img) {
        fresh = false
      }
    }
    if fresh {
      images.push(img)
    }
  }
  check(
    std.length() == images.length(),
    "S4 #standard monomials = #semigroup elements, got " +
      std.length().to_string() +
      " vs " +
      images.length().to_string(),
  )
  // and the map from standard monomials to semigroup elements is injective
  let seen : Array[Array[Int]] = []
  for i = 0; i < std.length(); i = i + 1 {
    let img = monomial_map_exponent(gens, std[i])
    for j = 0; j < seen.length(); j = j + 1 {
      check(!(vec_eq(seen[j], img)), "S4 injectivity of the standard-monomial map")
    }
    seen.push(img)
  }
}

fn binom3(n : Int) -> Int {
  n * (n + 1) / 2
}

fn vec_eq(a : Array[Int], b : Array[Int]) -> Bool {
  if a.length() != b.length() {
    return false
  }
  for i = 0; i < a.length(); i = i + 1 {
    if a[i] != b[i] {
      return false
    }
  }
  true
}

fn compositions_pub(k : Int, parts : Int) -> Array[Array[Int]] {
  if parts == 0 {
    return if k == 0 {
      [[]]
    } else {
      []
    }
  }
  let out : Array[Array[Int]] = []
  for first = 0; first <= k; first = first + 1 {
    let sub = compositions_pub(k - first, parts - 1)
    for i = 0; i < sub.length(); i = i + 1 {
      out.push(Array::add([first], sub[i]))
    }
  }
  out
}

test "S6 toric ideal: binomials vanish under the monomial map and lattice relations lie in it" {
  let gens = [[2, 0], [0, 2], [1, 1]]
  let p = 2
  let st = toric_ideal(gens, p, 3)
  // (a) every generator of the ideal maps to zero: its two monomials have the
  //     same image under x_i -> t^{a_i}
  for i = 0; i < st.gens.length(); i = i + 1 {
    let f = st.gens[i]
    check(f.terms.length() == 2, "binomial has two terms")
    let a = monomial_map_exponent(gens, f.terms[0].0)
    let b = monomial_map_exponent(gens, f.terms[1].0)
    check(vec_eq(a, b), "S6 the binomial vanishes under the monomial map")
    check(f.terms[0].1 + f.terms[1].1 == 0 || f.terms[0].1 == f.terms[1].1, "coefficients")
  }
  check(st.gens.length() > 0, "the toric ideal is nonzero")
  // (b) every lattice relation u with A u = 0 and |u_+|,|u_-| <= 4 gives a
  //     binomial that lies in the ideal built from degree <= 3 relations
  let lim = 4
  for u0 = -lim; u0 <= lim; u0 = u0 + 1 {
    for u1 = -lim; u1 <= lim; u1 = u1 + 1 {
      for u2 = -lim; u2 <= lim; u2 = u2 + 1 {
        let img = [
          2 * u0 + u2,
          2 * u1 + u2,
        ]
        if img[0] != 0 || img[1] != 0 {
          continue
        }
        let plus = [Int::max(0, u0), Int::max(0, u1), Int::max(0, u2)]
        let minus = [Int::max(0, -u0), Int::max(0, -u1), Int::max(0, -u2)]
        let dp = plus[0] + plus[1] + plus[2]
        let dm = minus[0] + minus[1] + minus[2]
        if dp > lim || dm > lim || (dp == 0 && dm == 0) {
          continue
        }
        let bin = @poly.mp_monom(p, 3, plus, 1).sub(@poly.mp_monom(p, 3, minus, 1))
        check(
          st.in_ideal(bin),
          "S6 lattice binomial in the toric ideal: " + bin.to_string(),
        )
      }
    }
  }
}

test "S5 Ehrhart: polynomiality, leading coefficient = volume, and reciprocity" {
  let triangle = [[0, 0], [1, 0], [0, 1]]
  let square = [[0, 0], [1, 0], [0, 1], [1, 1]]
  let big = [[0, 0], [2, 0], [0, 3]]
  let cases : Array[(Array[Array[Int]], Array[Int])] = [
    (triangle, [1, 3, 6, 10, 15, 21]),
    (square, [1, 4, 9, 16, 25, 36]),
    (big, [1, 7, 19, 37, 61, 91]),
  ]
  for c = 0; c < cases.length(); c = c + 1 {
    let pts = cases[c].0
    let want = cases[c].1
    let lf = ehrhart_fn(pts, 5)
    for k = 0; k <= 5; k = k + 1 {
      check(
        lf[k] == want[k],
        "S5 L_P(" + k.to_string() + ") case " + c.to_string() + " got " + lf[k].to_string(),
      )
    }
    // polynomiality: L_P(k) is a polynomial of degree n = 2, so the third
    // finite difference vanishes
    let coeffs = newton_coeffs(lf)
    check(coeffs.length() == 6, "six samples")
    for k = 0; k <= 5; k = k + 1 {
      check(newton_eval(coeffs, k) == lf[k], "interpolation agrees")
    }
    check(coeffs[3] == 0 && coeffs[4] == 0 && coeffs[5] == 0, "degree <= 2")
    // leading coefficient c_2 satisfies c_2 / 2! = volume(P)
    let area2 = if c == 0 {
      1
    } else if c == 1 {
      2
    } else {
      6
    }
    check(coeffs[2] == area2, "c_2 = 2! * volume, case " + c.to_string())
    // Ehrhart-Macdonald reciprocity: L_P(-k) = (-1)^n * #interior(k*P)
    for k = 1; k <= 4; k = k + 1 {
      let interior = dilated_interior(pts, k).length()
      let lhs = newton_eval(coeffs, -k)
      check(
        lhs == interior,
        "S5 reciprocity case " +
          c.to_string() +
          " k=" +
          k.to_string() +
          " got " +
          lhs.to_string() +
          " interior " +
          interior.to_string(),
      )
    }
  }
}
