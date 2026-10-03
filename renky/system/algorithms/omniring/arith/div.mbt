// ===========================================================================
//  Integer kernel: divisibility, primality, valuations.
//  (No third-party dependency: only `kernel`, which is itself proved.)
// ===========================================================================

///| Bézout certificate: `a*x + b*y == g == gcd(a,b)`.
pub struct Bezout {
  g : Int
  x : Int
  y : Int
}

pub fn Bezout::to_string(self : Bezout) -> String {
  "Bezout{g=" +
  self.g.to_string() +
  ",x=" +
  self.x.to_string() +
  ",y=" +
  self.y.to_string() +
  "}"
}

///| **A08 gcd** — Euclid's algorithm.
///| Definition. `gcd(a,b)` is the greatest common divisor, `gcd(a,0)=|a|`.
pub fn gcd(a : Int, b : Int) -> Int {
  let mut x = a.abs()
  let mut y = b.abs()
  while y != 0 {
    let r = x % y
    x = y
    y = r
  }
  x
}

///| **A09 extended gcd** — returns `g,x,y` with `a*x + b*y = g`.
pub fn egcd(a : Int, b : Int) -> Bezout {
  let sa = if a < 0 {
    -1
  } else {
    1
  }
  let sb = if b < 0 {
    -1
  } else {
    1
  }
  let mut r0 = a.abs()
  let mut r1 = b.abs()
  let mut x0 = 1
  let mut y0 = 0
  let mut x1 = 0
  let mut y1 = 1
  while r1 != 0 {
    let q = r0 / r1
    let r2 = r0 - q * r1
    let x2 = x0 - q * x1
    let y2 = y0 - q * y1
    r0 = r1
    r1 = r2
    x0 = x1
    y0 = y1
    x1 = x2
    y1 = y2
  }
  { g: r0, x: sa * x0, y: sb * y0 }
}

///| Least common multiple.
pub fn lcm(a : Int, b : Int) -> Int {
  let g = gcd(a, b)
  if g == 0 {
    0
  } else {
    (a.abs() / g) * b.abs()
  }
}

///| **A10 p-adic valuation** `v_p(n)`: the exponent of `p` in `n`.
///| `v_p(0) = -1` is our sentinel for "infinite valuation".
pub fn valuation(n : Int, p : Int) -> Int {
  if n == 0 {
    return -1
  }
  let mut m = n.abs()
  let mut v = 0
  while m % p == 0 {
    m = m / p
    v = v + 1
  }
  v
}

///| Integer power `b^e` (exact, no modular reduction), `e >= 0`.
pub fn pow_int(b : Int, e : Int) -> Int {
  let mut acc = 1
  for _i = 0; _i < e; _i = _i + 1 {
    acc = acc * b
  }
  acc
}

///| Base-`b` digits of `n ≥ 0`, least significant first.
pub fn digits(n : Int, b : Int) -> Array[Int] {
  let out : Array[Int] = []
  if n == 0 {
    out.push(0)
    return out
  }
  let mut m = n
  while m > 0 {
    out.push(m % b)
    m = m / b
  }
  out
}

///| Value of base-`b` digits (least significant first).
pub fn from_digits(ds : Array[Int], b : Int) -> Int {
  let mut acc = 0
  let mut i = ds.length() - 1
  while i >= 0 {
    acc = acc * b + ds[i]
    i = i - 1
  }
  acc
}

///| All positive divisors of `n > 0`, sorted.
pub fn divisors(n : Int) -> Array[Int] {
  let out : Array[Int] = []
  let mut d = 1
  while d * d <= n {
    if n % d == 0 {
      out.push(d)
      if d * d != n {
        out.push(n / d)
      }
    }
    d = d + 1
  }
  out.sort()
  out
}

///| Euler's totient `φ(n) = #{1 ≤ k ≤ n : gcd(k,n) = 1}`.
pub fn totient(n : Int) -> Int {
  let mut c = 0
  for k = 1; k <= n; k = k + 1 {
    if gcd(k, n) == 1 {
      c = c + 1
    }
  }
  c
}

///| Möbius function μ(n).
pub fn mobius(n : Int) -> Int {
  if n == 1 {
    return 1
  }
  let f = factorize(n)
  for i = 0; i < f.length(); i = i + 1 {
    if f[i].1 >= 2 {
      return 0
    }
  }
  if f.length() % 2 == 0 {
    1
  } else {
    -1
  }
}
