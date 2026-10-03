// ===========================================================================
//  義務自證 · the heat (Maslov) layer, verified on complete finite grids.
// ===========================================================================

fn close(a : Double, b : Double, tol : Double) -> Bool {
  Double::is_close(a, b, absolute_tolerance=tol)
}

fn grid() -> Array[Double] {
  [-4.0, -2.5, -1.0, -0.25, 0.0, 0.5, 1.0, 2.0, 3.5, 7.0]
}

test "exp and log are mutually inverse on the whole grid" {
  for i = 0; i < grid().length(); i = i + 1 {
    let x = grid()[i]
    check(close(logd(expd(x)), x, 1.0e-9), "log(exp(x)) = x at " + x.to_string())
    let y = expd(x)
    check(close(expd(logd(y)), y, 1.0e-6), "exp(log(y)) = y at " + y.to_string())
  }
  check(close(expd(0.0), 1.0, 1.0e-12), "exp(0) = 1")
  check(close(logd(1.0), 0.0, 1.0e-12), "log(1) = 0")
  check(close(expd(1.0), 2.718281828459045, 1.0e-9), "exp(1) = e")
  check(close(logd(2.0), ln2, 1.0e-12), "log(2) = ln2")
}

test "exp is a homomorphism (R,+) -> (R>0,*) on every pair of the grid" {
  let g = grid()
  for i = 0; i < g.length(); i = i + 1 {
    for j = 0; j < g.length(); j = j + 1 {
      check(
        close(expd(g[i] + g[j]), expd(g[i]) * expd(g[j]), 1.0e-6),
        "exp(a+b) = exp a * exp b",
      )
      let a = expd(g[i])
      let b = expd(g[j])
      check(close(logd(a * b), logd(a) + logd(b), 1.0e-9), "log(ab) = log a + log b")
    }
  }
}

test "phi_h is a semifield isomorphism: heat addition becomes ordinary addition" {
  let g = [-3.0, -1.0, 0.0, 0.7, 2.0, 5.0]
  let hs = [2.0, 1.0, 0.5, 0.1, 0.01]
  for t = 0; t < hs.length(); t = t + 1 {
    let h = hs[t]
    for i = 0; i < g.length(); i = i + 1 {
      for j = 0; j < g.length(); j = j + 1 {
        let a = g[i]
        let b = g[j]
        // phi(a (+)_h b) = phi(a) + phi(b)
        check(
          close(phi(heat_add(a, b, h), h), phi(a, h) + phi(b, h), 1.0e-6),
          "phi preserves heat addition, h=" + h.to_string(),
        )
        // phi(a (*)_h b) = phi(a) * phi(b)
        check(
          close(phi(heat_mul(a, b), h), phi(a, h) * phi(b, h), 1.0e-6),
          "phi preserves heat multiplication, h=" + h.to_string(),
        )
        // phi^{-1} really inverts phi
        check(close(phi_inv(phi(a, h), h), a, 1.0e-6), "phi_inv(phi(a)) = a")
      }
    }
  }
}

test "heat addition is commutative, associative, translation invariant and idempotent" {
  let g = [-3.0, -1.0, 0.0, 0.7, 2.0, 5.0]
  for t = 0; t < 3; t = t + 1 {
    let h = [1.0, 0.3, 0.05][t]
    for i = 0; i < g.length(); i = i + 1 {
      let a = g[i]
      // (+)_h is idempotent only in the h -> 0 limit; at finite h the exact
      // identity is a (+)_h a = a - h*log 2
      check(close(heat_add(a, a, h), a - h * ln2, 1.0e-9), "a (+)_h a = a - h ln2")
      for j = 0; j < g.length(); j = j + 1 {
        let b = g[j]
        check(close(heat_add(a, b, h), heat_add(b, a, h), 1.0e-12), "commutativity")
        // translation invariance = distributivity of (*)_h over (+)_h
        check(
          close(
            heat_mul(a, heat_add(b, g[(j + 1) % g.length()], h)),
            heat_add(heat_mul(a, b), heat_mul(a, g[(j + 1) % g.length()]), h),
            1.0e-9,
          ),
          "distributivity",
        )
        for k = 0; k < g.length(); k = k + 1 {
          let c = g[k]
          check(
            close(
              heat_add(heat_add(a, b, h), c, h),
              heat_add(a, heat_add(b, c, h), h),
              1.0e-6,
            ),
            "associativity",
          )
        }
      }
    }
  }
}

test "zero temperature limit: softmin -> min with the bound 0 <= gap <= h*ln2" {
  let g = grid()
  let hs = [4.0, 1.0, 0.3, 0.1, 0.03, 0.001]
  for i = 0; i < g.length(); i = i + 1 {
    for j = 0; j < g.length(); j = j + 1 {
      let a = g[i]
      let b = g[j]
      let mn = Double::min(a, b)
      let mut prev = -1.0e308
      for t = 0; t < hs.length(); t = t + 1 {
        let h = hs[t]
        let s = softmin_h(a, b, h)
        let gap = s - mn
        check(gap <= 1.0e-9, "softmin <= min")
        check(gap >= -(h * ln2) - 1.0e-9, "softmin >= min - h ln 2, h=" + h.to_string())
        check(s >= prev - 1.0e-9, "softmin increases as h decreases")
        prev = s
      }
      check(close(softmin_h(a, b, 1.0e-12), mn, 1.0e-9), "h -> 0 limit is min")
    }
  }
}

test "free energy, mean energy and entropy: Legendre relation and the h -> 0 limit" {
  let families : Array[Array[Double]] = [
    [0.0, 1.0],
    [3.0, -1.0, 2.0],
    [0.0, 0.0, 5.0, 5.0],
    [-2.0, 0.5, 0.5, 4.0, 9.0],
  ]
  let hs = [2.0, 1.0, 0.3, 0.1, 0.01, 0.001]
  for f = 0; f < families.length(); f = f + 1 {
    let e = families[f]
    let (mn, k) = to_tropical(e)
    let mut prev = -1.0e308
    for t = 0; t < hs.length(); t = t + 1 {
      let h = hs[t]
      let fe = free_energy(e, h)
      // F = <E> - h S  (Legendre relation)
      let (lhs, rhs) = legendre_check(e, h)
      check(close(lhs, rhs, 1.0e-6), "Legendre relation, family " + f.to_string())
      check(close(fe, rhs, 1.0e-9), "free energy = n-ary softmin")
      check(fe <= mn + 1.0e-9, "F <= min E")
      check(fe >= prev - 1.0e-9, "F is increasing in h")
      prev = fe
      check(close(fe, softmin_all(e, h), 1.0e-12), "F = softmin_all")
    }
    check(close(free_energy(e, 1.0e-12), mn, 1.0e-9), "F -> min E")
    check(close(mean_energy(e, 1.0e-12), mn, 1.0e-6), "<E> -> min E")
    // residual entropy: S -> log(g) where g is the degeneracy of the ground
    // state (g = 1 gives S -> 0, the usual third-law statement)
    let mut g = 0
    for i = 0; i < e.length(); i = i + 1 {
      if e[i] == mn {
        g = g + 1
      }
    }
    check(
      close(entropy(e, 1.0e-12), logd(Double::from_int(g)), 1.0e-6),
      "S -> log(degeneracy), family " + f.to_string(),
    )
    // the Gibbs state concentrates on the tropical argmin
    let ps = gibbs(e, 1.0e-9)
    let mut best = 0
    for i = 1; i < ps.length(); i = i + 1 {
      if ps[i] > ps[best] {
        best = i
      }
    }
    check(e[best] == mn, "ground state is the tropical minimum")
    check(best == k || e[k] == e[best], "tropical argmin agrees with the Gibbs argmax")
    // entropies are non-negative and bounded by log n
    check(entropy(e, 1.0) >= -1.0e-9, "S >= 0")
    check(entropy(e, 1.0) <= logd(Double::from_int(e.length())) + 1.0e-9, "S <= log n")
  }
}

fn check(cond : Bool, msg : String) -> Unit {
  if cond {
    ()
  } else {
    println("CHECK FAILED: " + msg)
    abort("[heat law failed] " + msg)
  }
}
