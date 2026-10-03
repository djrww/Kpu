// ===========================================================================
//  Primality, factorisation, quadratic residues.
// ===========================================================================

///| Modular exponentiation `b^e mod p` (binary method), `p ≥ 1`, `e ≥ 0`.
///| **A11 powmod**.
pub fn powmod(b : Int, e : Int, p : Int) -> Int {
  if p == 1 {
    return 0
  }
  let mut acc = 1
  let mut base = @kernel.modpos(b, p)
  let mut rest = e
  while rest > 0 {
    if rest % 2 == 1 {
      acc = @kernel.mulmod(acc, base, p)
    }
    base = @kernel.mulmod(base, base, p)
    rest = rest / 2
  }
  acc
}

///| Modular inverse via Bézout; returns `0` when `a` is not invertible mod `p`.
///| **A12 invmod**.
pub fn invmod(a : Int, p : Int) -> Int {
  let bz = egcd(a, p)
  if bz.g != 1 {
    return 0
  }
  @kernel.modpos(bz.x, p)
}

///| Modular inverse via Fermat's little theorem (requires `p` prime).
pub fn invmod_fermat(a : Int, p : Int) -> Int {
  powmod(a, p - 2, p)
}

///| Miller–Rabin witness test.  Deterministic for `n < 3_317_044_064_679_887_385_961_981`
///| when the first 13 prime bases are used; we use them all.
///| **A13 is_prime**.
pub fn is_prime(n : Int) -> Bool {
  if n < 2 {
    return false
  }
  let small = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]
  for i = 0; i < small.length(); i = i + 1 {
    if n == small[i] {
      return true
    }
    if n % small[i] == 0 {
      return false
    }
  }
  // n - 1 = d * 2^r with d odd
  let mut d = n - 1
  let mut r = 0
  while d % 2 == 0 {
    d = d / 2
    r = r + 1
  }
  let bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41]
  for k = 0; k < bases.length(); k = k + 1 {
    let a = bases[k]
    if a >= n {
      continue
    }
    let mut x = powmod(a, d, n)
    if x == 1 || x == n - 1 {
      continue
    }
    let mut composite = true
    let mut j = 0
    while j < r - 1 {
      x = @kernel.mulmod(x, x, n)
      if x == n - 1 {
        composite = false
        j = r
      } else {
        j = j + 1
      }
    }
    if composite {
      return false
    }
  }
  true
}

///| Smallest prime `> n`.
pub fn next_prime(n : Int) -> Int {
  let mut c = if n < 2 {
    2
  } else {
    n + 1
  }
  while !(is_prime(c)) {
    c = c + 1
  }
  c
}

///| Prime factorisation as sorted `(prime, exponent)` pairs.
///| Trial division plus **A14 Pollard rho** for the surviving cofactor.
pub fn factorize(n : Int) -> Array[(Int, Int)] {
  let out : Array[(Int, Int)] = []
  let mut m = n.abs()
  if m <= 1 {
    return out
  }
  let mut d = 2
  while d * d <= m && d < 100000 {
    while m % d == 0 {
      push_factor(out, d)
      m = m / d
    }
    d = d + 1
  }
  if m > 1 {
    split_semiprime(out, m)
  }
  out
}

fn push_factor(out : Array[(Int, Int)], p : Int) -> Unit {
  if out.length() > 0 && out[out.length() - 1].0 == p {
    let e = out[out.length() - 1]
    out[out.length() - 1] = (e.0, e.1 + 1)
  } else {
    out.push((p, 1))
  }
}

fn split_semiprime(out : Array[(Int, Int)], m : Int) -> Unit {
  if is_prime(m) {
    push_factor(out, m)
    return
  }
  let f = pollard_rho(m)
  split_semiprime(out, f)
  split_semiprime(out, m / f)
}

///| Pollard's rho with Brent's improvement.
pub fn pollard_rho(n : Int) -> Int {
  if n % 2 == 0 {
    return 2
  }
  let mut c = 1
  while c < n {
    let mut x = 2
    let mut y = 2
    let mut d = 1
    let mut steps = 0
    while d == 1 && steps < 100000 {
      x = @kernel.addmod(@kernel.mulmod(x, x, n), c, n)
      y = @kernel.addmod(@kernel.mulmod(y, y, n), c, n)
      y = @kernel.addmod(@kernel.mulmod(y, y, n), c, n)
      d = gcd((x - y).abs(), n)
      steps = steps + 1
    }
    if d != 1 && d != n {
      return d
    }
    c = c + 1
  }
  n
}

///| Jacobi symbol `(a/n)` for odd `n ≥ 1`.  **A15 jacobi**.
pub fn jacobi(a : Int, n : Int) -> Int {
  if n <= 0 || n % 2 == 0 {
    abort("jacobi: n must be a positive odd integer")
  }
  let mut x = @kernel.modpos(a, n)
  let mut m = n
  let mut result = 1
  while x != 0 {
    while x % 2 == 0 {
      x = x / 2
      if m % 8 == 3 || m % 8 == 5 {
        result = -result
      }
    }
    let t = x
    x = m
    m = t
    if x % 4 == 3 && m % 4 == 3 {
      result = -result
    }
    x = x % m
  }
  if m == 1 {
    result
  } else {
    0
  }
}

///| Legendre symbol `(a/p)` for prime `p`.
pub fn legendre(a : Int, p : Int) -> Int {
  if a % p == 0 {
    0
  } else {
    jacobi(a, p)
  }
}

///| Euler's criterion as an executable identity: `a^((p-1)/2) ≡ (a/p)`.
pub fn euler_criterion(a : Int, p : Int) -> Int {
  powmod(a, (p - 1) / 2, p)
}

///| Is `a` a quadratic residue mod prime `p`?
pub fn is_qr(a : Int, p : Int) -> Bool {
  @kernel.modpos(a, p) == 0 || legendre(a, p) == 1
}

///| **A16 Tonelli–Shanks**: a square root of `a` mod prime `p`, or `-1`.
pub fn sqrt_mod(a : Int, p : Int) -> Int {
  let aa = @kernel.modpos(a, p)
  if aa == 0 {
    return 0
  }
  if legendre(aa, p) != 1 {
    return -1
  }
  if p % 4 == 3 {
    return powmod(aa, (p + 1) / 4, p)
  }
  // write p - 1 = q * 2^s with q odd
  let mut q = p - 1
  let mut s = 0
  while q % 2 == 0 {
    q = q / 2
    s = s + 1
  }
  let mut z = 2
  while legendre(z, p) != -1 {
    z = z + 1
  }
  let mut m = s
  let mut c = powmod(z, q, p)
  let mut t = powmod(aa, q, p)
  let mut r = powmod(aa, (q + 1) / 2, p)
  while t != 1 {
    let mut i = 0
    let mut tt = t
    while tt != 1 {
      tt = @kernel.mulmod(tt, tt, p)
      i = i + 1
      if i == m {
        return -1
      }
    }
    let mut b = c
    for _k = 0; _k < m - i - 1; _k = _k + 1 {
      b = @kernel.mulmod(b, b, p)
    }
    m = i
    c = @kernel.mulmod(b, b, p)
    t = @kernel.mulmod(t, c, p)
    r = @kernel.mulmod(r, b, p)
  }
  r
}

///| Chinese remainder theorem for two coprime moduli.  **A17 crt**.
pub fn crt(r1 : Int, m1 : Int, r2 : Int, m2 : Int) -> Int {
  let g = gcd(m1, m2)
  if g != 1 {
    abort("crt: moduli must be coprime")
  }
  let bz = egcd(m1, m2)
  let m = m1 * m2
  @kernel.modpos(r1 * bz.y * m2 + r2 * bz.x * m1, m)
}
