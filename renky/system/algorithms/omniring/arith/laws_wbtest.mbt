// ===========================================================================
//  義務自證 (obligatory self-proof), layer 2: EXHAUSTIVE verification.
//
//  Everything below is a *complete* check, not a sample: the quantified
//  variables range over the entire finite domain named in the test title.
//  For a finite ring a universally quantified identity is decidable by
//  enumeration, so these tests are proofs of the stated theorems for those
//  rings.
// ===========================================================================

fn check(cond : Bool, msg : String) -> Unit {
  if cond {
    ()
  } else {
    abort("[law failed] " + msg)
  }
}

fn naive_powmod(b : Int, e : Int, p : Int) -> Int {
  let mut acc = @kernel.modpos(1, p)
  for _i = 0; _i < e; _i = _i + 1 {
    acc = @kernel.mulmod(acc, b, p)
  }
  acc
}

fn brute_gcd(a : Int, b : Int) -> Int {
  let mut best = 1
  for d = 1; d <= a; d = d + 1 {
    if a % d == 0 && b % d == 0 {
      best = d
    }
  }
  best
}

test "A08 gcd equals the brute-force greatest common divisor on [1,60]^2" {
  for a = 1; a <= 60; a = a + 1 {
    for b = 1; b <= 60; b = b + 1 {
      check(gcd(a, b) == brute_gcd(a, b), "gcd " + a.to_string() + " " + b.to_string())
    }
  }
}

test "A09 extended gcd: Bezout identity a*x+b*y == gcd(a,b) on [-40,40]^2" {
  for a = -40; a <= 40; a = a + 1 {
    for b = -40; b <= 40; b = b + 1 {
      if a == 0 && b == 0 {
        continue
      }
      let z = egcd(a, b)
      check(a * z.x + b * z.y == z.g, "bezout " + a.to_string() + " " + b.to_string())
      check(z.g == gcd(a, b), "gcd agree " + a.to_string() + " " + b.to_string())
      check(z.g >= 0, "gcd nonneg")
    }
  }
}

test "lcm * gcd == |a*b| on [1,40]^2" {
  for a = 1; a <= 40; a = a + 1 {
    for b = 1; b <= 40; b = b + 1 {
      check(lcm(a, b) * gcd(a, b) == a * b, "lcm " + a.to_string())
    }
  }
}

test "A11 powmod equals naive repeated multiplication, p in {5,7,13,101}, b<p, e<=40" {
  let primes = [5, 7, 13, 101]
  for k = 0; k < primes.length(); k = k + 1 {
    let p = primes[k]
    for b = 0; b < p; b = b + 1 {
      for e = 0; e <= 40; e = e + 1 {
        check(
          powmod(b, e, p) == naive_powmod(b, e, p),
          "powmod " + b.to_string() + "^" + e.to_string() + " mod " + p.to_string(),
        )
      }
    }
  }
}

test "Fermat's little theorem: a^p == a for every a in Z/p, every prime p <= 100" {
  for p = 2; p <= 100; p = p + 1 {
    if !(is_prime(p)) {
      continue
    }
    for a = 0; a < p; a = a + 1 {
      check(powmod(a, p, p) == a, "fermat " + a.to_string() + " p=" + p.to_string())
    }
  }
}

test "A12 invmod: a * a^-1 == 1 for every unit of Z/p, every prime p <= 60" {
  for p = 2; p <= 60; p = p + 1 {
    if !(is_prime(p)) {
      continue
    }
    for a = 1; a < p; a = a + 1 {
      let inv = invmod(a, p)
      check(
        @kernel.mulmod(a, inv, p) == 1,
        "invmod " + a.to_string() + " p=" + p.to_string(),
      )
      check(
        invmod_fermat(a, p) == inv,
        "fermat inverse agrees " + a.to_string() + " p=" + p.to_string(),
      )
    }
  }
}

test "A13 Miller-Rabin agrees with trial division on every n <= 2000" {
  for n = 0; n <= 2000; n = n + 1 {
    let mut trial = n >= 2
    let mut d = 2
    while d * d <= n {
      if n % d == 0 {
        trial = false
        d = n
      }
      d = d + 1
    }
    check(is_prime(n) == trial, "is_prime " + n.to_string())
  }
}

test "A14 factorization: product of prime powers recovers n, for every n <= 600" {
  for n = 2; n <= 600; n = n + 1 {
    let f = factorize(n)
    let mut prod = 1
    for i = 0; i < f.length(); i = i + 1 {
      check(is_prime(f[i].0), "factor prime " + n.to_string())
      check(f[i].1 >= 1, "factor exponent")
      let mut pw = 1
      for _k = 0; _k < f[i].1; _k = _k + 1 {
        pw = pw * f[i].0
      }
      prod = prod * pw
      if i > 0 {
        check(f[i - 1].0 < f[i].0, "factors sorted " + n.to_string())
      }
    }
    check(prod == n, "factorization product " + n.to_string())
  }
}

test "A10 valuation: p^v divides n but p^(v+1) does not, n <= 200, p in {2,3,5,7}" {
  let ps = [2, 3, 5, 7]
  for k = 0; k < ps.length(); k = k + 1 {
    let p = ps[k]
    for n = 1; n <= 200; n = n + 1 {
      let v = valuation(n, p)
      let mut pw = 1
      for _i = 0; _i < v; _i = _i + 1 {
        pw = pw * p
      }
      check(n % pw == 0, "valuation divides " + n.to_string())
      check(n % (pw * p) != 0, "valuation maximal " + n.to_string())
    }
  }
}

test "digits/from_digits round trip for every n <= 500 in bases 2..12" {
  for b = 2; b <= 12; b = b + 1 {
    for n = 0; n <= 500; n = n + 1 {
      check(from_digits(digits(n, b), b) == n, "digits " + n.to_string())
    }
  }
}

test "totient counts coprimes exactly, for every n <= 120" {
  for n = 1; n <= 120; n = n + 1 {
    let mut c = 0
    for k = 1; k <= n; k = k + 1 {
      if brute_gcd(k, n) == 1 {
        c = c + 1
      }
    }
    check(totient(n) == c, "totient " + n.to_string())
  }
}

test "Möbius inversion: sum_{d|n} mu(d) == [n == 1] and phi(n) == sum mu(d)*(n/d), n <= 120" {
  for n = 1; n <= 120; n = n + 1 {
    let ds = divisors(n)
    let mut s = 0
    let mut t = 0
    for i = 0; i < ds.length(); i = i + 1 {
      let d = ds[i]
      check(n % d == 0, "divisor divides")
      s = s + mobius(d)
      t = t + mobius(d) * (n / d)
    }
    let expect = if n == 1 {
      1
    } else {
      0
    }
    check(s == expect, "mobius sum " + n.to_string())
    check(t == totient(n), "mobius totient " + n.to_string())
  }
}

test "A15 Jacobi symbol is multiplicative in the numerator, all odd n <= 51" {
  for n = 3; n <= 51; n = n + 2 {
    for a = -n; a <= 2 * n; a = a + 1 {
      for b = -n; b <= 2 * n; b = b + 1 {
        let lhs = jacobi(a * b, n)
        let rhs = jacobi(a, n) * jacobi(b, n)
        check(lhs == rhs, "jacobi mult n=" + n.to_string())
      }
    }
  }
}

test "Euler criterion: a^((p-1)/2) == (a/p) for every a in Z/p, every prime p <= 100" {
  for p = 3; p <= 100; p = p + 1 {
    if !(is_prime(p)) {
      continue
    }
    for a = 0; a < p; a = a + 1 {
      let e = euler_criterion(a, p)
      let l = legendre(a, p)
      let expect = if l == 0 {
        0
      } else if l == 1 {
        1
      } else {
        p - 1
      }
      check(e == expect, "euler " + a.to_string() + " p=" + p.to_string())
      // exactly (p-1)/2 nonzero squares
    }
  }
}

test "A16 Tonelli-Shanks: r^2 == a for every residue a of every prime p <= 200" {
  for p = 3; p <= 200; p = p + 1 {
    if !(is_prime(p)) {
      continue
    }
    let mut squares = 0
    for a = 0; a < p; a = a + 1 {
      let r = sqrt_mod(a, p)
      if is_qr(a, p) {
        squares = squares + 1
        check(r >= 0, "sqrt exists " + a.to_string() + " p=" + p.to_string())
        check(
          @kernel.mulmod(r, r, p) == @kernel.modpos(a, p),
          "sqrt correct " + a.to_string() + " p=" + p.to_string(),
        )
      } else {
        check(r == -1, "non-residue " + a.to_string() + " p=" + p.to_string())
      }
    }
    check(squares == (p + 1) / 2, "half are squares p=" + p.to_string())
  }
}

test "A17 CRT: unique solution in [0, m1*m2) for all coprime pairs m1,m2 <= 18" {
  for m1 = 2; m1 <= 18; m1 = m1 + 1 {
    for m2 = 2; m2 <= 18; m2 = m2 + 1 {
      if gcd(m1, m2) != 1 {
        continue
      }
      for r1 = 0; r1 < m1; r1 = r1 + 1 {
        for r2 = 0; r2 < m2; r2 = r2 + 1 {
          let x = crt(r1, m1, r2, m2)
          check(0 <= x && x < m1 * m2, "crt range")
          check(x % m1 == r1, "crt m1")
          check(x % m2 == r2, "crt m2")
        }
      }
    }
  }
}

test "A14 Pollard rho returns a proper factor of every composite n <= 300" {
  for n = 4; n <= 300; n = n + 1 {
    if is_prime(n) {
      continue
    }
    let f = pollard_rho(n)
    check(f > 1 && f < n, "rho proper " + n.to_string())
    check(n % f == 0, "rho divides " + n.to_string())
  }
}

fn prime_trial(n : Int) -> Bool {
  if n < 2 {
    return false
  }
  let mut d = 2
  while d * d <= n {
    if n % d == 0 {
      return false
    }
    d = d + 1
  }
  true
}

test "next_prime enumerates exactly the primes, and skips none, up to 500" {
  let mut prev = 1
  let mut count = 0
  let mut q = next_prime(prev)
  while q <= 500 {
    check(is_prime(q), "next prime is prime")
    check(q > prev, "next prime grows")
    let mut between = prev + 1
    while between < q {
      check(!(is_prime(between)), "no prime skipped " + between.to_string())
      between = between + 1
    }
    prev = q
    count = count + 1
    q = next_prime(prev)
  }
  let mut brute = 0
  for n = 2; n <= 500; n = n + 1 {
    if prime_trial(n) {
      brute = brute + 1
    }
  }
  check(count == brute, "pi(500): chain " + count.to_string() + " vs trial " + brute.to_string())
  check(count == 95, "pi(500) == 95, got " + count.to_string())
}
