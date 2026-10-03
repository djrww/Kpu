// ===========================================================================
//  義務自證 · K1–K7 for the ideal class group and its correspondence ring.
//  Every law is checked by complete enumeration of the (finite) class group.
// ===========================================================================

fn check(cond : Bool, msg : String) -> Unit {
  if cond {
    ()
  } else {
    println("CHECK FAILED: " + msg)
    abort("[classgroup law failed] " + msg)
  }
}

let test_discs : Array[Int] = [-15, -20, -23, -24, -31, -47]

test "K1 Gauss reduction: idempotent, discriminant preserving, unique per class" {
  // every primitive form with small coefficients reduces into the reduced list
  for a = 1; a <= 6; a = a + 1 {
    for b = -14; b <= 14; b = b + 1 {
      // c is determined by requiring a negative discriminant
      for c = 1; c <= 20; c = c + 1 {
        let f = bqf(a, b, c)
        let d = f.disc()
        if d >= 0 || @kernel.modpos(d, 4) > 1 {
          continue
        }
        if !(f.is_primitive()) {
          continue
        }
        let r = bqf_reduce(f)
        check(r.disc() == d, "K1 reduction preserves the discriminant")
        check(r.is_reduced(), "K1 the result is reduced, f=" + f.to_string())
        check(bqf_reduce(r).equal(r), "K1 reduction is idempotent")
        check(r.a > 0, "K1 positive definite")
        // it must be one of the enumerated reduced forms of that discriminant
        let list = reduced_forms(d)
        let mut found = false
        for i = 0; i < list.length(); i = i + 1 {
          if list[i].equal(r) {
            found = true
          }
        }
        check(found, "K1 reduction lands in the reduced list, f=" + f.to_string())
        // reduction is invariant under the SL_2(Z) move (a,b,c) -> (c,-b,a)
        let sw = bqf(c, -b, a)
        check(bqf_reduce(sw).equal(r), "K1 reduction is SL_2(Z) invariant")
        // and under b -> b + 2ak
        for k = -2; k <= 2; k = k + 1 {
          let nb = b + 2 * a * k
          let nc = (nb * nb - d) / (4 * a)
          check(
            bqf_reduce(bqf(a, nb, nc)).equal(r),
            "K1 reduction is invariant under translation of b",
          )
        }
      }
    }
  }
  // reduced forms of the same class are equal: uniqueness of the representative
  for di = 0; di < test_discs.length(); di = di + 1 {
    let d = test_discs[di]
    let cl = reduced_forms(d)
    for i = 0; i < cl.length(); i = i + 1 {
      for j = 0; j < cl.length(); j = j + 1 {
        check(
          bqf_equivalent(d, cl[i], cl[j]) == cl[i].equal(cl[j]),
          "K1 uniqueness of the reduced representative, D=" + d.to_string(),
        )
      }
    }
  }
}

test "K2 Cl(D) is a finite abelian group under composition" {
  for di = 0; di < test_discs.length(); di = di + 1 {
    let d = test_discs[di]
    let cl = reduced_forms(d)
    let one = bqf_principal(d)
    let mut idx = -1
    for i = 0; i < cl.length(); i = i + 1 {
      if cl[i].equal(one) {
        idx = i
      }
    }
    check(idx >= 0, "K2 the principal form is reduced, D=" + d.to_string())
    for i = 0; i < cl.length(); i = i + 1 {
      // identity
      check(
        bqf_compose(d, one, cl[i]).equal(cl[i]) &&
          bqf_compose(d, cl[i], one).equal(cl[i]),
        "K2 identity, D=" + d.to_string(),
      )
      // inverse
      check(
        is_principal(d, bqf_compose(d, cl[i], cl[i].inverse())),
        "K2 inverse, D=" + d.to_string(),
      )
      check(
        bqf_reduce(cl[i].inverse()).disc() == d,
        "K2 the inverse has the same discriminant",
      )
      for j = 0; j < cl.length(); j = j + 1 {
        // closure
        let prod = bqf_compose(d, cl[i], cl[j])
        check(prod.disc() == d, "K2 closure: the discriminant is preserved")
        let mut found = false
        for k = 0; k < cl.length(); k = k + 1 {
          if cl[k].equal(prod) {
            found = true
          }
        }
        check(found, "K2 closure, D=" + d.to_string())
        // commutativity
        check(
          prod.equal(bqf_compose(d, cl[j], cl[i])),
          "K2 commutativity, D=" + d.to_string(),
        )
        for t = 0; t < cl.length(); t = t + 1 {
          // associativity
          check(
            bqf_compose(d, bqf_compose(d, cl[i], cl[j]), cl[t]).equal(
              bqf_compose(d, cl[i], bqf_compose(d, cl[j], cl[t])),
            ),
            "K2 associativity, D=" + d.to_string(),
          )
        }
      }
    }
  }
}

test "K3 the classical class numbers" {
  let known = [
    (-3, 1),
    (-4, 1),
    (-7, 1),
    (-8, 1),
    (-11, 1),
    (-15, 2),
    (-19, 1),
    (-20, 2),
    (-23, 3),
    (-24, 2),
    (-31, 3),
    (-40, 2),
    (-47, 5),
    (-71, 7),
    (-163, 1),
  ]
  for i = 0; i < known.length(); i = i + 1 {
    let (d, h) = known[i]
    check(
      class_number(d) == h,
      "K3 h(" + d.to_string() + ") = " + h.to_string() + ", got " + class_number(d).to_string(),
    )
  }
  // D with h = 1 are exactly the Euler idylls among the listed ones
  check(class_number(-163) == 1, "K3 Euler: h(-163) = 1")
}

test "K4 Lagrange: g^h = 1 and the order of every class divides h" {
  for di = 0; di < test_discs.length(); di = di + 1 {
    let d = test_discs[di]
    let cl = reduced_forms(d)
    let h = cl.length()
    for i = 0; i < cl.length(); i = i + 1 {
      check(is_principal(d, bqf_pow(d, cl[i], h)), "K4 g^h = 1, D=" + d.to_string())
      let ord = class_order(d, cl[i])
      check(ord >= 1 && ord <= h, "K4 the order is in [1,h]")
      check(h % ord == 0, "K4 Lagrange: the order divides h, D=" + d.to_string())
      check(is_principal(d, bqf_pow(d, cl[i], ord)), "K4 g^ord = 1")
      if ord > 1 {
        check(
          !(is_principal(d, bqf_pow(d, cl[i], ord - 1))),
          "K4 the order is minimal",
        )
      }
      // the cyclic subgroup generated by cl[i] has exactly `ord` elements
      let seen : Array[BQF] = []
      for k = 0; k < ord; k = k + 1 {
        let e = bqf_pow(d, cl[i], k)
        let mut dup = false
        for t = 0; t < seen.length(); t = t + 1 {
          if seen[t].equal(e) {
            dup = true
          }
        }
        check(!(dup), "K4 the powers are distinct")
        seen.push(e)
      }
      check(seen.length() == ord, "K4 the subgroup has `ord` elements")
    }
  }
}

///| The norm of the ideal of a form = |det| of its Z-basis.
fn ideal_norm(d : Int, f : BQF) -> Int {
  let b = form_to_ideal(d, f)
  (b[0][0] * b[1][1] - b[0][1] * b[1][0]).abs()
}

test "K5 correspondence: ideal norms multiply, and primes split by the Kronecker symbol" {
  for di = 0; di < test_discs.length(); di = di + 1 {
    let d = test_discs[di]
    let cl = reduced_forms(d)
    // the norm of the ideal of a reduced form is its first coefficient
    for i = 0; i < cl.length(); i = i + 1 {
      check(
        ideal_norm(d, cl[i]) == form_norm(cl[i]),
        "K5 N(I(f)) = a, D=" + d.to_string(),
      )
    }
    // multiplicativity of the ideal norm at the lattice level
    for i = 0; i < cl.length(); i = i + 1 {
      for j = 0; j < cl.length(); j = j + 1 {
        let prod = ideal_mul(d, form_to_ideal(d, cl[i]), form_to_ideal(d, cl[j]))
        let det = (prod[0][0] * prod[1][1] - prod[0][1] * prod[1][0]).abs()
        check(
          det == ideal_norm(d, cl[i]) * ideal_norm(d, cl[j]),
          "K5 N(IJ) = N(I)N(J), D=" + d.to_string(),
        )
      }
    }
    // prime decomposition: ramified / split / inert
    for q in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41] {
      let kr = kronecker(d, q)
      let divides = d % q == 0
      check((kr == 0) == divides, "K5 q ramifies iff q | D, q=" + q.to_string())
      // b^2 = D (mod 4q) is solvable iff q splits
      let mut solvable = false
      let mut witness = 0
      for b = 0; b < 4 * q; b = b + 1 {
        if @kernel.modpos(b * b - d, 4 * q) == 0 {
          solvable = true
          witness = b
          break
        }
      }
      check(
        solvable == (kr == 1 || divides),
        "K5 splitting criterion, D=" + d.to_string() + " q=" + q.to_string(),
      )
      if kr == 1 {
        // the form above q is primitive of discriminant D, and its conjugate is
        // its inverse in the class group
        let bb = if witness > 2 * q {
          witness - 4 * q
        } else {
          witness
        }
        let num = bb * bb - d
        check(num % (4 * q) == 0, "K5 integrality of c")
        let f = bqf(q, bb, num / (4 * q))
        check(f.disc() == d, "K5 the form above q has discriminant D")
        check(f.is_primitive(), "K5 the form above q is primitive")
        let r = bqf_reduce(f)
        let mut found = false
        for i = 0; i < cl.length(); i = i + 1 {
          if cl[i].equal(r) {
            found = true
          }
        }
        check(found, "K5 the form above q lies in the class group")
        check(
          is_principal(d, bqf_compose(d, f, f.inverse())),
          "K5 the conjugate prime is the inverse",
        )
        // the class above q has order dividing h
        check(
          is_principal(d, bqf_pow(d, r, cl.length())),
          "K5 Lagrange for the class above q",
        )
      }
    }
  }
}

test "K6 genus theory: 2^(t-1) genera, characters are homomorphisms, equal fibre sizes" {
  for di = 0; di < test_discs.length(); di = di + 1 {
    let d = test_discs[di]
    let cl = reduced_forms(d)
    let t = @arith.factorize(d.abs()).length()
    let discs = genus_discriminants(d)
    check(
      discs.length() == @arith.pow_int(2, t) - 1,
      "K6 2^t - 1 discriminant divisors, D=" + d.to_string(),
    )
    check(
      genus_count(d) == @arith.pow_int(2, t - 1),
      "K6 2^(t-1) genera, D=" + d.to_string(),
    )
    // each character is a group homomorphism into {+-1} and is well defined on
    // the class (independent of the representative)
    for e = 0; e < discs.length(); e = e + 1 {
      let dd = discs[e]
      for i = 0; i < cl.length(); i = i + 1 {
        let v = genus_char(dd, cl[i], d)
        check(v == 1 || v == -1, "K6 the character takes values in {+-1}")
        // well defined: the same value on equivalent non-reduced representatives
        for k = -2; k <= 2; k = k + 1 {
          let a = cl[i].a
          let b = cl[i].b + 2 * a * k
          let c = (b * b - d) / (4 * a)
          check(
            genus_char(dd, bqf(a, b, c), d) == v,
            "K6 the genus character is well defined on the class",
          )
        }
        for j = 0; j < cl.length(); j = j + 1 {
          check(
            genus_char(dd, bqf_compose(d, cl[i], cl[j]), d) ==
              genus_char(dd, cl[i], d) * genus_char(dd, cl[j], d),
            "K6 the genus character is a homomorphism, D=" + d.to_string(),
          )
        }
      }
    }
    // the product of all genus characters (with d = D included) is trivial
    for i = 0; i < cl.length(); i = i + 1 {
      let mut prod = 1
      for e = 0; e < discs.length(); e = e + 1 {
        prod = prod * genus_char(discs[e], cl[i], d)
      }
      check(
        prod == 1,
        "K6 the product of all genus characters is trivial, D=" +
          d.to_string() +
          " class=" +
          cl[i].to_string() +
          " discs=" +
          fmt_ints(discs) +
          " prod=" +
          prod.to_string(),
      )
    }
    // the genera (fibres of the character map) all have the same size
    let keys : Array[String] = []
    let sizes : Array[Int] = []
    for i = 0; i < cl.length(); i = i + 1 {
      let parts : Array[String] = []
      for e = 0; e < discs.length(); e = e + 1 {
        parts.push(genus_char(discs[e], cl[i], d).to_string())
      }
      let key = parts.join(",")
      let mut found = -1
      for t2 = 0; t2 < keys.length(); t2 = t2 + 1 {
        if keys[t2] == key {
          found = t2
        }
      }
      if found < 0 {
        keys.push(key)
        sizes.push(1)
      } else {
        sizes[found] = sizes[found] + 1
      }
    }
    check(
      keys.length() == @arith.pow_int(2, t - 1),
      "K6 the number of genera is 2^(t-1), got " + keys.length().to_string(),
    )
    for i = 0; i < sizes.length(); i = i + 1 {
      check(
        sizes[i] * keys.length() == cl.length(),
        "K6 all genera have the same size, D=" + d.to_string(),
      )
    }
  }
}

test "K7 the correspondence ring F_p[Cl(D)] is a graded commutative ring" {
  for p in [2, 3] {
    for di = 0; di < 3; di = di + 1 {
      let d = test_discs[di]
      let r = class_ring(p, d)
      let els = r.elements()
      let one = r.unit_vec()
      for i = 0; i < els.length(); i = i + 1 {
        let x = els[i]
        for j = 0; j < els.length(); j = j + 1 {
          let y = els[j]
          check(
            vec_eq(r.mul_vec(x, y), r.mul_vec(y, x)),
            "K7 commutative, D=" + d.to_string(),
          )
          check(
            vec_eq(r.mul_vec(x, one), x),
            "K7 unit, D=" + d.to_string(),
          )
          for t = 0; t < els.length(); t = t + 1 {
            let z = els[t]
            check(
              vec_eq(r.mul_vec(r.mul_vec(x, y), z), r.mul_vec(x, r.mul_vec(y, z))),
              "K7 associative, D=" + d.to_string(),
            )
            check(
              vec_eq(
                r.mul_vec(x, r.add_vec(y, z)),
                r.add_vec(r.mul_vec(x, y), r.mul_vec(x, z)),
              ),
              "K7 distributive, D=" + d.to_string(),
            )
          }
        }
        // the augmentation is a ring homomorphism
        for j = 0; j < els.length(); j = j + 1 {
          check(
            r.augmentation(r.mul_vec(x, els[j])) ==
              @kernel.modpos(r.augmentation(x) * r.augmentation(els[j]), p),
            "K7 the augmentation is multiplicative",
          )
        }
        check(
          r.augmentation(r.add_vec(x, els[0])) ==
            @kernel.modpos(r.augmentation(x) + r.augmentation(els[0]), p),
          "K7 the augmentation is additive",
        )
      }
      // grading: the basis elements multiply as the group, and every row of the
      // table is a permutation (the group law)
      let n = r.dim()
      for i = 0; i < n; i = i + 1 {
        let seen : Array[Bool] = Array::make(n, false)
        for j = 0; j < n; j = j + 1 {
          let k = r.table[i][j]
          check(!(seen[k]), "K7 each row of the table is a permutation")
          seen[k] = true
          // [g_i] * [g_j] = [g_i g_j]
          let ei : Array[Int] = Array::make(n, 0)
          ei[i] = 1
          let ej : Array[Int] = Array::make(n, 0)
          ej[j] = 1
          let prod = r.mul_vec(ei, ej)
          let ek : Array[Int] = Array::make(n, 0)
          ek[k] = 1
          check(vec_eq(prod, ek), "K7 the ring is graded by Cl(D)")
        }
      }
      // every basis element is a unit in the ring
      for i = 0; i < n; i = i + 1 {
        let ei : Array[Int] = Array::make(n, 0)
        ei[i] = 1
        let inv : Array[Int] = Array::make(n, 0)
        // find the index j with table[i][j] = one
        for j = 0; j < n; j = j + 1 {
          if r.table[i][j] == r.one {
            inv[j] = 1
          }
        }
        check(vec_eq(r.mul_vec(ei, inv), one), "K7 group elements are units")
      }
    }
  }
}

fn fmt_ints(xs : Array[Int]) -> String {
  let parts : Array[String] = []
  for i = 0; i < xs.length(); i = i + 1 {
    parts.push(xs[i].to_string())
  }
  "[" + parts.join(",") + "]"
}

fn vec_eq(a : Array[Int], b : Array[Int]) -> Bool {
  let n = Int::max(a.length(), b.length())
  for i = 0; i < n; i = i + 1 {
    let x = if i < a.length() {
      a[i]
    } else {
      0
    }
    let y = if i < b.length() {
      b[i]
    } else {
      0
    }
    if x != y {
      return false
    }
  }
  true
}
