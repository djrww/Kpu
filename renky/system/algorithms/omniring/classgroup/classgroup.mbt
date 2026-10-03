// ===========================================================================
//  classgroup — the ideal class group of a quadratic order and its
//  correspondence ring (「理想類群對應環」, structure 8).
//
//  Fix a negative discriminant D = 0 or 1 (mod 4) and the quadratic order
//      O_D = Z[omega],   omega = (D + sqrt(D))/2,   omega^2 = D*omega - D(D-1)/4,
//  so an element of O_D is a pair (x, y) meaning x + y*omega and
//      (x + y w)(x' + y' w) = (x x' - y y' D(D-1)/4) + (x y' + y x' + y y' D) w.
//
//  Definition (binary quadratic form).  f = (a, b, c) is a*x^2 + b*x*y + c*y^2
//  of discriminant disc(f) = b^2 - 4ac.  f is *primitive* if gcd(a,b,c) = 1 and
//  *reduced* (D < 0, a > 0) if |b| <= a <= c and (b >= 0 whenever |b| = a or
//  a = c).
//
//  Definition (correspondence).  The primitive form (a,b,c) of discriminant D
//  corresponds to the invertible O_D-ideal
//      I(f) = Z*a + Z*((-b + sqrt(D))/2) = Z*(a, 0) + Z*((-b-D)/2, 1)
//  written in the (x, y) coordinates of omega.  Proper equivalence classes of
//  primitive forms correspond bijectively to ideal classes, and composition of
//  forms corresponds to ideal multiplication:
//      I(f * g) = I(f) * I(g),
//  the product being computed as the Z-span of the four products of the basis
//  elements, put back into a basis by the Hermite normal form (A29).
//
//  Definition (correspondence ring).  Cl(D) is a finite abelian group; its group
//  algebra F_p[Cl(D)] is a commutative ring graded by Cl(D): the product of the
//  homogeneous basis elements [g] and [h] is [g*h].  This is the
//  "ideal class group correspondence ring".
//
//  Theorems (verified exhaustively in `classgroup_wbtest`):
//   K1 Gauss reduction: every primitive form is equivalent to a unique reduced
//      form; reduction is idempotent and preserves the discriminant.
//   K2 Cl(D) is a finite abelian group under composition: closure,
//      commutativity, associativity, identity = the principal form, and
//      (a,b,c)^{-1} = (a,-b,c).
//   K3 the class numbers are the classical values
//      h(-3)=h(-4)=h(-7)=h(-8)=h(-11)=h(-163)=1, h(-15)=h(-20)=h(-24)=2,
//      h(-23)=h(-31)=3, h(-47)=5.
//   K4 Lagrange: g^{h(D)} = 1 for every class, and the order of every class
//      divides h(D).
//   K5 the norm is multiplicative and corresponds to the ideal norm; a prime q
//      is ramified iff q | D, split iff the Kronecker symbol (D/q) = 1, inert
//      iff (D/q) = -1; for split q the form (q, b, (b^2-D)/4q) is primitive of
//      discriminant D and its conjugate is its inverse in Cl(D).
//   K6 genus theory: the number of genera is 2^{t-1} with t the number of
//      distinct prime divisors of D; every genus character is a well defined
//      group homomorphism Cl(D) -> {+-1} and their product is trivial.
//   K7 F_p[Cl(D)] is a commutative associative unital ring graded by Cl(D), and
//      the augmentation epsilon(sum c_g [g]) = sum c_g is a ring homomorphism.
// ===========================================================================

///| A binary quadratic form a*x^2 + b*x*y + c*y^2.
pub struct BQF {
  a : Int
  b : Int
  c : Int
}

pub fn bqf(a : Int, b : Int, c : Int) -> BQF {
  { a, b, c }
}

pub fn BQF::disc(self : BQF) -> Int {
  self.b * self.b - 4 * self.a * self.c
}

pub fn BQF::to_string(self : BQF) -> String {
  "(" + self.a.to_string() + "," + self.b.to_string() + "," + self.c.to_string() + ")"
}

pub fn BQF::equal(self : BQF, other : BQF) -> Bool {
  self.a == other.a && self.b == other.b && self.c == other.c
}

pub fn BQF::is_primitive(self : BQF) -> Bool {
  @arith.gcd(@arith.gcd(self.a.abs(), self.b.abs()), self.c.abs()) == 1
}

///| Is the form reduced (D < 0, a > 0)?
pub fn BQF::is_reduced(self : BQF) -> Bool {
  if self.a <= 0 {
    return false
  }
  let ab = self.b.abs()
  if ab > self.a || self.a > self.c {
    return false
  }
  if (ab == self.a || self.a == self.c) && self.b < 0 {
    return false
  }
  true
}

///| **A82 Gauss reduction** for negative discriminants.  Every step preserves
///| the discriminant: after b is shifted by 2ak, c is recomputed as (b^2-D)/4a.
pub fn bqf_reduce(f : BQF) -> BQF {
  let d = f.disc()
  if d >= 0 {
    abort("[classgroup.bqf_reduce] only negative discriminants are reduced here")
  }
  let mut a = f.a
  let mut b = f.b
  let mut c = f.c
  if a < 0 {
    a = -a
    b = -b
    c = -c
  }
  let mut steps = 0
  while steps < 10000 {
    steps = steps + 1
    // shift b into (-a, a] and recompute c from the discriminant
    let nb = shift_b(a, b)
    if nb != b {
      b = nb
      c = (b * b - d) / (4 * a)
    }
    if a > c {
      let na = c
      let nc = a
      a = na
      c = nc
      b = -b
      continue
    }
    if b < 0 && (b.abs() == a || a == c) {
      b = -b
      continue
    }
    break
  }
  if b * b - 4 * a * c != d {
    abort("[classgroup.bqf_reduce] discriminant changed")
  }
  { a, b, c }
}

fn shift_b(a : Int, b : Int) -> Int {
  // move b into (-a, a]
  let mut x = b
  while x > a {
    x = x - 2 * a
  }
  while x <= -a {
    x = x + 2 * a
  }
  x
}

///| The principal form of discriminant D.
pub fn bqf_principal(d : Int) -> BQF {
  let b = @kernel.modpos(d, 2)
  { a: 1, b, c: (b * b - d) / 4 }
}

///| The inverse class: the conjugate form.
pub fn BQF::inverse(self : BQF) -> BQF {
  { a: self.a, b: -self.b, c: self.c }
}

///| All primitive reduced forms of discriminant D (< 0).
pub fn reduced_forms(d : Int) -> Array[BQF] {
  let r = @kernel.modpos(d, 4)
  if d >= 0 || (r != 0 && r != 1) {
    abort("[classgroup.reduced_forms] D must be a negative discriminant")
  }
  let out : Array[BQF] = []
  // a <= sqrt(|D|/3)
  let mut amax = 1
  while (amax + 1) * (amax + 1) * 3 <= -d {
    amax = amax + 1
  }
  for a = 1; a <= amax + 1; a = a + 1 {
    for b = -a; b <= a; b = b + 1 {
      let num = b * b - d
      if num % (4 * a) != 0 {
        continue
      }
      let c = num / (4 * a)
      let f : BQF = { a, b, c }
      if !(f.is_primitive()) || !(f.is_reduced()) {
        continue
      }
      out.push(f)
    }
  }
  out
}

///| **K3 class number** h(D).
pub fn class_number(d : Int) -> Int {
  reduced_forms(d).length()
}

// ---------------------------------------------------------------------------
//  A53 arithmetic in the quadratic order and ideal multiplication
// ---------------------------------------------------------------------------

///| **A83 quadratic-order arithmetic**: product in O_D = Z[omega], omega = (D + sqrt D)/2.
fn omega_mul(d : Int, x : Int, y : Int, x2 : Int, y2 : Int) -> (Int, Int) {
  let k = d * (d - 1) / 4
  (x * x2 - y * y2 * k, x * y2 + y * x2 + y * y2 * d)
}

///| **K5 form -> ideal**: the Z-basis (a, 0), ((-b-D)/2 mod a, 1) of I(f).
pub fn form_to_ideal(d : Int, f : BQF) -> @linalg.Mat {
  let y0 = @kernel.modpos((-f.b - d) / 2, f.a)
  [[f.a, 0], [y0, 1]]
}

///| N(x + y*omega) in O_D.
fn omega_norm(d : Int, x : Int, y : Int) -> Int {
  x * x + x * y * d + y * y * (d * d - d) / 4
}

///| Tr((x1 + y1 w) * conj(x2 + y2 w)).
fn omega_tr(d : Int, x1 : Int, y1 : Int, x2 : Int, y2 : Int) -> Int {
  2 * x1 * x2 + (x1 * y2 + x2 * y1) * d + 2 * y1 * y2 * (d * d - d) / 4
}

///| **K5 ideal -> form**: for a Z-basis (alpha, beta) of a proper O_D-ideal the
///| associated form is
///|     f(x,y) = N(alpha*x + beta*y) / N(I),   N(I) = |alpha*conj(beta) - conj(alpha)*beta|/sqrt(D)
///| which in omega-coordinates is just |det| of the basis matrix.  Changing the
///| basis by an element of SL_2(Z) gives an equivalent form, so the Hermite
///| normal form yields a canonical class, and reduction gives its representative.
pub fn ideal_to_form(d : Int, basis : @linalg.Mat) -> BQF {
  let rows : Array[Array[Int]] = []
  for i = 0; i < basis.length(); i = i + 1 {
    if basis[i][0] != 0 || basis[i][1] != 0 {
      rows.push(basis[i])
    }
  }
  if rows.length() != 2 {
    abort("[classgroup.ideal_to_form] the ideal must have rank 2")
  }
  let det = rows[0][0] * rows[1][1] - rows[0][1] * rows[1][0]
  let n = det.abs()
  if n == 0 {
    abort("[classgroup.ideal_to_form] degenerate basis")
  }
  let na = omega_norm(d, rows[0][0], rows[0][1])
  let nb = omega_norm(d, rows[1][0], rows[1][1])
  let tr = omega_tr(d, rows[0][0], rows[0][1], rows[1][0], rows[1][1])
  if na % n != 0 || nb % n != 0 || tr % n != 0 {
    abort("[classgroup.ideal_to_form] the lattice is not a proper O_D-ideal")
  }
  let a = na / n
  let b = -tr / n
  let c = nb / n
  if b * b - 4 * a * c != d {
    abort("[classgroup.ideal_to_form] wrong discriminant")
  }
  { a, b, c }
}

///| **A83 (cont.) ideal multiplication**: multiply two O_D-ideals given by Z-bases: the span of the four products,
///| returned in Hermite normal form (A29).
pub fn ideal_mul(d : Int, i : @linalg.Mat, j : @linalg.Mat) -> @linalg.Mat {
  let prods : @linalg.Mat = []
  for r = 0; r < i.length(); r = r + 1 {
    for s = 0; s < j.length(); s = s + 1 {
      let (x, y) = omega_mul(d, i[r][0], i[r][1], j[s][0], j[s][1])
      prods.push([x, y])
    }
  }
  @linalg.hnf(prods)
}

///| **A84 Gauss composition via ideal multiplication**, followed by reduction.
pub fn bqf_compose(d : Int, f : BQF, g : BQF) -> BQF {
  let prod = ideal_mul(d, form_to_ideal(d, f), form_to_ideal(d, g))
  bqf_reduce(ideal_to_form(d, prod))
}

///| Composition then reduction, restricted to the class group.
pub fn class_mul(d : Int, classes : Array[BQF], f : BQF, g : BQF) -> BQF {
  let r = bqf_compose(d, f, g)
  for i = 0; i < classes.length(); i = i + 1 {
    if classes[i].equal(r) {
      return r
    }
  }
  abort("[classgroup.class_mul] the product left the class group")
}

///| f^n by binary exponentiation in the class group.
pub fn bqf_pow(d : Int, f : BQF, n : Int) -> BQF {
  let mut acc = bqf_principal(d)
  let mut base = f
  let mut e = n.abs()
  while e > 0 {
    if e % 2 == 1 {
      acc = bqf_compose(d, acc, base)
    }
    base = bqf_compose(d, base, base)
    e = e / 2
  }
  if n < 0 {
    return bqf_compose(d, acc, bqf_principal(d))
  }
  acc
}

///| Is the form in the principal class?
pub fn is_principal(d : Int, f : BQF) -> Bool {
  bqf_reduce(f).equal(bqf_principal(d))
}

///| Proper equivalence of two primitive forms: f ~ g iff f * g^{-1} is principal.
pub fn bqf_equivalent(d : Int, f : BQF, g : BQF) -> Bool {
  is_principal(d, bqf_compose(d, bqf_reduce(f), bqf_reduce(g).inverse()))
}

///| The order of a class.
pub fn class_order(d : Int, f : BQF) -> Int {
  let mut k = 1
  let mut cur = f
  while k <= class_number(d) {
    if is_principal(d, cur) {
      return k
    }
    cur = bqf_compose(d, cur, f)
    k = k + 1
  }
  0
}

///| The norm of the ideal corresponding to f (= |a|).
pub fn form_norm(f : BQF) -> Int {
  f.a.abs()
}

// ---------------------------------------------------------------------------
//  K6  Kronecker symbol and genus characters
// ---------------------------------------------------------------------------

///| **A85 Kronecker symbol** (d/n) for a discriminant d and n >= 1.
pub fn kronecker(d : Int, n : Int) -> Int {
  if n == 0 {
    return if d == 1 || d == -1 {
      1
    } else {
      0
    }
  }
  if n < 0 {
    abort("[classgroup.kronecker] n must be positive")
  }
  let mut result = 1
  let mut m = n
  // the power of two
  while m % 2 == 0 {
    m = m / 2
    if d % 2 == 0 {
      return 0
    }
    let d2 = @kernel.modpos(d, 8)
    if d2 == 3 || d2 == 5 {
      result = -result
    }
  }
  if m == 1 {
    return result
  }
  let fac = @arith.factorize(m)
  for i = 0; i < fac.length(); i = i + 1 {
    let (q, e) = fac[i]
    let s = @arith.jacobi(d, q)
    for _k = 0; _k < e; _k = _k + 1 {
      result = result * s
    }
  }
  result
}

///| The discriminant divisors of D: those d | D with d and D/d both
///| discriminants (0 or 1 mod 4) and d != 1.  Together with d = 1 they give the
///| 2^{t-1} genus characters.
pub fn genus_discriminants(d : Int) -> Array[Int] {
  let out : Array[Int] = []
  let ds = @arith.divisors(d.abs())
  for i = 0; i < ds.length(); i = i + 1 {
    let e = ds[i]
    for sgn in [-1, 1] {
      let cand = sgn * e
      if cand == 1 {
        continue
      }
      let other = d / cand
      if is_disc(cand) && is_disc(other) && cand * other == d {
        out.push(cand)
      }
    }
  }
  out.sort_by(fn(x, y) { x - y })
  out
}

fn is_disc(x : Int) -> Bool {
  let r = @kernel.modpos(x, 4)
  r == 0 || r == 1
}

///| All representatives of the proper equivalence class of `f` reachable by the
///| SL_2(Z) moves (a,b,c) -> (a, b+2ak, (b+2ak)^2-D)/4a) and (a,b,c) -> (c,-b,a),
///| explored by breadth-first search up to `limit` forms.
pub fn equivalent_forms(f : BQF, limit : Int) -> Array[BQF] {
  let d = f.disc()
  let seen : Array[BQF] = []
  let queue : Array[BQF] = [f]
  let mut qi = 0
  while qi < queue.length() && seen.length() < limit {
    let cur = queue[qi]
    qi = qi + 1
    let mut dup = false
    for i = 0; i < seen.length(); i = i + 1 {
      if seen[i].equal(cur) {
        dup = true
      }
    }
    if dup {
      continue
    }
    seen.push(cur)
    // b -> b + 2ak
    for k = -3; k <= 3; k = k + 1 {
      let nb = cur.b + 2 * cur.a * k
      let num = nb * nb - d
      if num % (4 * cur.a) == 0 {
        queue.push({ a: cur.a, b: nb, c: num / (4 * cur.a) })
      }
    }
    // the swap (a,b,c) -> (c,-b,a)
    if cur.c != 0 {
      queue.push({ a: cur.c, b: -cur.b, c: cur.a })
    }
  }
  seen
}

///| The value of the genus character chi_e on the class of f: (e/a) for any
///| representative with gcd(a, e) = 1.  Genus theory guarantees the value does
///| not depend on the choice; every coprime representative found is checked to
///| agree.
pub fn genus_char(e : Int, f : BQF, d : Int) -> Int {
  ignore(d)
  let cands = equivalent_forms(bqf_reduce(f), 300)
  let mut result = 0
  let mut found = false
  for i = 0; i < cands.length(); i = i + 1 {
    let a = cands[i].a
    if a <= 0 || @arith.gcd(a.abs(), e.abs()) != 1 {
      continue
    }
    let v = kronecker(e, a.abs())
    if !(found) {
      result = v
      found = true
    } else if v != result {
      abort("[classgroup.genus_char] the genus character is not well defined")
    }
  }
  if !(found) {
    abort("[classgroup.genus_char] no representative coprime to the conductor")
  }
  result
}

///| **K6 number of genera**: 2^{t-1} with t the number of distinct prime
///| divisors of D.
pub fn genus_count(d : Int) -> Int {
  let t = @arith.factorize(d.abs()).length()
  @arith.pow_int(2, t - 1)
}

// ---------------------------------------------------------------------------
//  K7  the correspondence ring F_p[Cl(D)]
// ---------------------------------------------------------------------------

///| The group algebra F_p[Cl(D)], graded by the class group.
pub struct ClassRing {
  p : Int
  d : Int
  classes : Array[BQF]
  table : Array[Array[Int]]
  one : Int
}

///| **A86 correspondence ring**: build the group algebra F_p[Cl(D)].
pub fn class_ring(p : Int, d : Int) -> ClassRing {
  let classes = reduced_forms(d)
  let n = classes.length()
  let table : Array[Array[Int]] = []
  for i = 0; i < n; i = i + 1 {
    let row : Array[Int] = []
    for j = 0; j < n; j = j + 1 {
      let prod = bqf_compose(d, classes[i], classes[j])
      let mut idx = -1
      for k = 0; k < n; k = k + 1 {
        if classes[k].equal(prod) {
          idx = k
        }
      }
      if idx < 0 {
        abort("[classgroup.class_ring] composition left the class group")
      }
      row.push(idx)
    }
    table.push(row)
  }
  let mut one = 0
  for k = 0; k < n; k = k + 1 {
    if classes[k].equal(bqf_principal(d)) {
      one = k
    }
  }
  { p, d, classes, table, one }
}

pub fn ClassRing::dim(self : ClassRing) -> Int {
  self.classes.length()
}

///| Product in the group algebra.
pub fn ClassRing::mul_vec(self : ClassRing, u : Array[Int], v : Array[Int]) -> Array[Int] {
  let n = self.dim()
  let out : Array[Int] = Array::make(n, 0)
  for i = 0; i < n; i = i + 1 {
    for j = 0; j < n; j = j + 1 {
      let k = self.table[i][j]
      out[k] = @kernel.addmod(
        out[k],
        @kernel.mulmod(self.at(u, i), self.at(v, j), self.p),
        self.p,
      )
    }
  }
  out
}

pub fn ClassRing::add_vec(self : ClassRing, u : Array[Int], v : Array[Int]) -> Array[Int] {
  let n = self.dim()
  let out : Array[Int] = Array::make(n, 0)
  for i = 0; i < n; i = i + 1 {
    out[i] = @kernel.addmod(self.at(u, i), self.at(v, i), self.p)
  }
  out
}

fn ClassRing::at(self : ClassRing, v : Array[Int], i : Int) -> Int {
  if i < v.length() {
    @kernel.modpos(v[i], self.p)
  } else {
    0
  }
}

///| The identity element [1] of the group algebra.
pub fn ClassRing::unit_vec(self : ClassRing) -> Array[Int] {
  let v : Array[Int] = Array::make(self.dim(), 0)
  v[self.one] = 1
  v
}

///| The augmentation epsilon(sum c_g [g]) = sum c_g.
pub fn ClassRing::augmentation(self : ClassRing, v : Array[Int]) -> Int {
  let mut acc = 0
  for i = 0; i < self.dim(); i = i + 1 {
    acc = @kernel.addmod(acc, self.at(v, i), self.p)
  }
  acc
}

///| All elements of the group algebra (p^n of them; keep n small).
pub fn ClassRing::elements(self : ClassRing) -> Array[Array[Int]] {
  let n = self.dim()
  let total = @arith.pow_int(self.p, n)
  let out : Array[Array[Int]] = []
  for k = 0; k < total; k = k + 1 {
    let ds = @arith.digits(k, self.p)
    let v : Array[Int] = Array::make(n, 0)
    for i = 0; i < n; i = i + 1 {
      v[i] = if i < ds.length() {
        ds[i]
      } else {
        0
      }
    }
    out.push(v)
  }
  out
}
