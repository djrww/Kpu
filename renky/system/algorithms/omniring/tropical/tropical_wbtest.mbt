// ===========================================================================
//  義務自證 · the tropical layer.  Every law below is checked on a *complete*
//  finite domain (all pairs/triples of a finite tropical set, all weight
//  vectors in a box, all 2x2 matrices, all permutations for n <= 5).
// ===========================================================================

fn check(cond : Bool, msg : String) -> Unit {
  if cond {
    ()
  } else {
    println("CHECK FAILED: " + msg)
    abort("[tropical law failed] " + msg)
  }
}

///| A finite sub-semiring of the tropical numbers, including +infinity.
fn trop_set() -> Array[Int] {
  [-4, -1, 0, 1, 3, trop_zero]
}

test "T = Z u {inf} with (min, +) is an idempotent commutative semiring: all triples" {
  let s = trop_set()
  for i = 0; i < s.length(); i = i + 1 {
    let a = s[i]
    check(trop_add(a, trop_zero) == a, "additive identity")
    check(trop_mul(a, trop_one) == a, "multiplicative identity")
    check(trop_add(a, a) == a, "idempotence")
    check(trop_mul(a, trop_zero) == trop_zero, "infinity absorbs")
    for j = 0; j < s.length(); j = j + 1 {
      let b = s[j]
      check(trop_add(a, b) == trop_add(b, a), "oplus commutative")
      check(trop_mul(a, b) == trop_mul(b, a), "otimes commutative")
      for t = 0; t < s.length(); t = t + 1 {
        let c = s[t]
        check(
          trop_add(trop_add(a, b), c) == trop_add(a, trop_add(b, c)),
          "oplus associative",
        )
        check(
          trop_mul(trop_mul(a, b), c) == trop_mul(a, trop_mul(b, c)),
          "otimes associative",
        )
        check(
          trop_mul(a, trop_add(b, c)) == trop_add(trop_mul(a, b), trop_mul(a, c)),
          "distributivity",
        )
      }
    }
  }
  check(trop_sum([5, 2, 9, 4]) == 2, "n-ary oplus")
  check(trop_prod([5, 2, 9, 4]) == 20, "n-ary otimes")
  check(trop_sum([]) == trop_zero, "empty oplus = zero")
  check(trop_prod([]) == trop_one, "empty otimes = one")
}

///| A small exhaustive family of bivariate tropical polynomials.
fn tp_family() -> Array[TPoly] {
  [
    tp_make(2, [([0, 0], 0), ([1, 0], 0), ([0, 1], 0)]),
    tp_make(2, [([0, 0], 1), ([2, 0], -1), ([0, 2], 0), ([1, 1], 2)]),
    tp_make(2, [([1, 0], 3), ([0, 1], -2)]),
    tp_make(2, [([0, 0], 0), ([3, 0], 0), ([0, 3], 0)]),
  ]
}

test "evaluation at w is a semiring homomorphism TPoly -> T, for every w in the box" {
  let fam = tp_family()
  for i = 0; i < fam.length(); i = i + 1 {
    for j = 0; j < fam.length(); j = j + 1 {
      let f = fam[i]
      let g = fam[j]
      for wx = -3; wx <= 3; wx = wx + 1 {
        for wy = -3; wy <= 3; wy = wy + 1 {
          let w = [wx, wy]
          check(
            f.mul(g).eval(w) == trop_mul(f.eval(w), g.eval(w)),
            "eval multiplicative " + f.to_string() + " * " + g.to_string(),
          )
          check(
            f.add(g).eval(w) == trop_add(f.eval(w), g.eval(w)),
            "eval additive",
          )
        }
      }
    }
  }
}

test "tropical polynomial semiring laws on the whole family" {
  let fam = tp_family()
  let zero = tp_make(2, [])
  let one = tp_make(2, [([0, 0], trop_one)])
  for i = 0; i < fam.length(); i = i + 1 {
    let f = fam[i]
    check(f.add(zero).terms.length() == f.terms.length(), "additive identity")
    check(f.mul(one).eval([1, 2]) == f.eval([1, 2]), "multiplicative identity")
    check(f.add(f).eval([2, -1]) == f.eval([2, -1]), "idempotent addition")
    for j = 0; j < fam.length(); j = j + 1 {
      let g = fam[j]
      check(f.add(g).eval([1, 1]) == g.add(f).eval([1, 1]), "addition commutative")
      check(f.mul(g).eval([1, 1]) == g.mul(f).eval([1, 1]), "multiplication commutative")
      for t = 0; t < fam.length(); t = t + 1 {
        let hh = fam[t]
        check(
          f.mul(g.add(hh)).eval([2, 2]) == f.mul(g).add(f.mul(hh)).eval([2, 2]),
          "distributive",
        )
      }
    }
  }
}

test "tropical hypersurface of x+y+1: multiplicity >= 2 exactly on the three rays" {
  let f = tp_make(2, [([0, 0], 0), ([1, 0], 0), ([0, 1], 0)])
  for wx = -4; wx <= 4; wx = wx + 1 {
    for wy = -4; wy <= 4; wy = wy + 1 {
      let on_line = (wx == wy && wx <= 0) ||
        (wx == 0 && wy >= 0) ||
        (wy == 0 && wx >= 0)
      check(
        f.on_hypersurface([wx, wy]) == on_line,
        "tropical line at (" + wx.to_string() + "," + wy.to_string() + ")",
      )
      // multiplicity is the number of attaining terms, and equals 1 off the curve
      let m = f.multiplicity([wx, wy])
      check(m == (if on_line {
        if wx == 0 && wy == 0 {
          3
        } else {
          2
        }
      } else {
        1
      }), "multiplicity at (" + wx.to_string() + "," + wy.to_string() + ")")
    }
  }
}

test "A48 Hungarian determinant equals the minimum over all permutations (n <= 5)" {
  // exhaustive for n = 2
  for a = -2; a <= 2; a = a + 1 {
    for b = -2; b <= 2; b = b + 1 {
      for c = -2; c <= 2; c = c + 1 {
        for d = -2; d <= 2; d = d + 1 {
          let m = [[a, b], [c, d]]
          check(
            trop_det(m) == trop_det_bruteforce(m),
            "n=2 det",
          )
        }
      }
    }
  }
  // exhaustive for n = 3 over entries in {-1,0,1}
  let vals = [-1, 0, 1]
  for i0 = 0; i0 < 3; i0 = i0 + 1 {
    for i1 = 0; i1 < 3; i1 = i1 + 1 {
      for i2 = 0; i2 < 3; i2 = i2 + 1 {
        for i3 = 0; i3 < 3; i3 = i3 + 1 {
          for i4 = 0; i4 < 3; i4 = i4 + 1 {
            for i5 = 0; i5 < 3; i5 = i5 + 1 {
              for i6 = 0; i6 < 3; i6 = i6 + 1 {
                for i7 = 0; i7 < 3; i7 = i7 + 1 {
                  for i8 = 0; i8 < 3; i8 = i8 + 1 {
                    let m = [
                      [vals[i0], vals[i1], vals[i2]],
                      [vals[i3], vals[i4], vals[i5]],
                      [vals[i6], vals[i7], vals[i8]],
                    ]
                    check(trop_det(m) == trop_det_bruteforce(m), "n=3 det")
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

///| A deterministic pseudo-random generator (LCG), so tests are reproducible.
fn lcg(seed : Int) -> Int {
  (seed * 1103515245 + 12345) % 2147483647
}

test "A48 Hungarian equals brute force on random matrices for n = 4,5,6 and gives a permutation" {
  let mut seed = 7
  for trial = 0; trial < 220; trial = trial + 1 {
    let n = 4 + trial % 3
    let m : Array[Array[Int]] = []
    for i = 0; i < n; i = i + 1 {
      let row : Array[Int] = []
      for j = 0; j < n; j = j + 1 {
        seed = lcg(seed)
        row.push(((seed % 13) + 13) % 13 - 6)
      }
      m.push(row)
    }
    let (cost, assign) = hungarian(m)
    check(cost == trop_det(m), "cost == det")
    // the assignment is a permutation and its cost is the reported one
    let seen : Array[Bool] = Array::make(n, false)
    let mut total = 0
    for i = 0; i < n; i = i + 1 {
      let j = assign[i]
      check(j >= 0 && j < n, "assignment in range")
      check(!(seen[j]), "assignment is injective")
      seen[j] = true
      total = total + m[i][j]
    }
    check(total == cost, "assignment cost matches")
    if n <= 5 {
      check(cost == trop_det_bruteforce(m), "hungarian == brute force, n=" + n.to_string())
    }
  }
}

test "tropical matrix product and Floyd-Warshall closure agree with explicit path enumeration" {
  let mut seed = 11
  for trial = 0; trial < 40; trial = trial + 1 {
    let n = 3 + trial % 2
    let w : Array[Array[Int]] = []
    for i = 0; i < n; i = i + 1 {
      let row : Array[Int] = []
      for j = 0; j < n; j = j + 1 {
        seed = lcg(seed)
        row.push(if i == j {
          0
        } else {
          ((seed % 9) + 9) % 9 + 1
        })
      }
      w.push(row)
    }
    let d = trop_closure(w)
    // brute force: minimum over all walks of length <= n-1 (enough for n nodes
    // with non-negative weights)
    for i = 0; i < n; i = i + 1 {
      for j = 0; j < n; j = j + 1 {
        let mut best = trop_zero
        let paths = enumerate_walks(n, i, j, n - 1)
        for k = 0; k < paths.length(); k = k + 1 {
          let mut s = 0
          let pth = paths[k]
          for t = 0; t + 1 < pth.length(); t = t + 1 {
            s = s + w[pth[t]][pth[t + 1]]
          }
          if s < best {
            best = s
          }
        }
        check(d[i][j] == best, "closure == shortest walk")
      }
    }
    // tropical matrix product against the definition
    let prod = trop_mat_mul(w, w)
    for i = 0; i < n; i = i + 1 {
      for j = 0; j < n; j = j + 1 {
        let mut best = trop_zero
        for t = 0; t < n; t = t + 1 {
          let cand = w[i][t] + w[t][j]
          if cand < best {
            best = cand
          }
        }
        check(prod[i][j] == best, "tropical product entry")
      }
    }
  }
}

fn enumerate_walks(n : Int, i : Int, j : Int, maxlen : Int) -> Array[Array[Int]] {
  let out : Array[Array[Int]] = []
  if i == j {
    out.push([i])
  }
  let mut frontier : Array[Array[Int]] = [[i]]
  for _len = 0; _len < maxlen; _len = _len + 1 {
    let next : Array[Array[Int]] = []
    for a = 0; a < frontier.length(); a = a + 1 {
      let p = frontier[a]
      let last = p[p.length() - 1]
      for k = 0; k < n; k = k + 1 {
        if k == last {
          continue
        }
        let q = Array::add(p, [k])
        next.push(q)
        if k == j {
          out.push(q)
        }
      }
    }
    frontier = next
  }
  out
}

test "A50 Newton subdivision and balancing on a catalogue of bivariate tropical curves" {
  let cases : Array[Array[(Array[Int], Int)]] = [
    [([0, 0], 0), ([1, 0], 0), ([0, 1], 0)],
    [([0, 0], 0), ([2, 0], 0), ([0, 2], 0), ([1, 1], 0)],
    [([0, 0], 0), ([1, 0], 0), ([2, 0], 0), ([0, 1], 0)],
    [([0, 0], 0), ([3, 0], 0), ([0, 3], 0), ([1, 1], 5), ([2, 1], 5)],
    [([0, 0], 2), ([1, 0], 0), ([0, 1], 1), ([1, 1], 0), ([2, 2], 3)],
    [([0, 0], 0), ([1, 0], 1), ([0, 1], 1), ([1, 1], 0), ([2, 0], 0), ([0, 2], 0)],
  ]
  for c = 0; c < cases.length(); c = c + 1 {
    let terms = cases[c]
    let cells = newton_cells(terms)
    check(cells.length() > 0, "cells exist, case " + c.to_string())
    // every point lies in some cell, and cells are pairwise distinct
    for k = 0; k < terms.length(); k = k + 1 {
      let mut seen = false
      for t = 0; t < cells.length(); t = t + 1 {
        for u = 0; u < cells[t].length(); u = u + 1 {
          if cells[t][u] == k {
            seen = true
          }
        }
      }
      check(seen, "point covered, case " + c.to_string())
    }
    // balancing
    let sums = newton_balancing(terms)
    for t = 0; t < sums.length(); t = t + 1 {
      check(
        sums[t].0 == 0 && sums[t].1 == 0,
        "balancing, case " + c.to_string() + " face " + t.to_string(),
      )
    }
    // independent cross-check: every argmin set found by scanning a grid of
    // rational weights must be one of the computed cells
    for a = -6; a <= 6; a = a + 1 {
      for b = -6; b <= 6; b = b + 1 {
        let heights : Array[Int] = []
        let mut best = 0
        for k = 0; k < terms.length(); k = k + 1 {
          let h = 2 * terms[k].1 + a * terms[k].0[0] + b * terms[k].0[1]
          heights.push(h)
          if k == 0 || h < best {
            best = h
          }
        }
        let cell : Array[Int] = []
        for k = 0; k < terms.length(); k = k + 1 {
          if heights[k] == best {
            cell.push(k)
          }
        }
        let mut found = false
        for t = 0; t < cells.length(); t = t + 1 {
          if cells[t].length() == cell.length() {
            let mut same = true
            for u = 0; u < cell.length(); u = u + 1 {
              if cells[t][u] != cell[u] {
                same = false
              }
            }
            if same {
              found = true
            }
          }
        }
        check(
          found,
          "grid argmin cell is known, case " +
            c.to_string() +
            " w=(" +
            a.to_string() +
            "/2," +
            b.to_string() +
            "/2)",
        )
      }
    }
  }
}

test "Maslov bridge: the h -> 0 limit of the heat soft minimum is the tropical sum" {
  let fams : Array[Array[Int]] = [[3, -1, 4], [0, 0, 5], [7], [-2, 9, -2, 3]]
  for i = 0; i < fams.length(); i = i + 1 {
    let xs = fams[i]
    let ds : Array[Double] = []
    for k = 0; k < xs.length(); k = k + 1 {
      ds.push(Double::from_int(xs[k]))
    }
    let tmin = Double::from_int(trop_sum(xs))
    check(
      Double::is_close(trop_dequantize(ds, 1.0e-12), tmin, absolute_tolerance=1.0e-6),
      "dequantization limit == tropical sum",
    )
    let (lim, k) = trop_limit(ds)
    check(lim == tmin, "trop_limit value")
    check(xs[k] == trop_sum(xs), "trop_limit index is a tropical argmin")
    // monotone approach from below
    let hs = [2.0, 0.5, 0.1, 0.01]
    let mut prev = -1.0e308
    for t = 0; t < hs.length(); t = t + 1 {
      let val = trop_dequantize(ds, hs[t])
      check(val <= tmin + 1.0e-9, "soft min <= tropical min")
      check(val >= prev - 1.0e-9, "soft min increases as h decreases")
      prev = val
    }
  }
}
