// ===========================================================================
//  義務自證 · X1–X5 for perfect F-split rings.  The whole family of monomial
//  quotients F_p[x,y]/(x^a, y^b, M) is enumerated (every subset M of the box),
//  and for each member F-splitness, reducedness and perfection are computed by
//  independent routes and compared.
// ===========================================================================

fn check(cond : Bool, msg : String) -> Unit {
  if cond {
    ()
  } else {
    println("CHECK FAILED: " + msg)
    abort("[fsplit law failed] " + msg)
  }
}

///| The monomials strictly inside the box 0 <= e_i < box_i.
fn box_monomials(box : Array[Int]) -> Array[Array[Int]] {
  let mut total = 1
  for i = 0; i < box.length(); i = i + 1 {
    total = total * box[i]
  }
  let out : Array[Array[Int]] = []
  for k = 0; k < total; k = k + 1 {
    let e : Array[Int] = []
    let mut rest = k
    for i = 0; i < box.length(); i = i + 1 {
      e.push(rest % box[i])
      rest = rest / box[i]
    }
    out.push(e)
  }
  out
}

fn monom_ideal(
  p : Int,
  nvars : Int,
  box : Array[Int],
  mask : Int,
) -> @ideal.IdealState {
  let st = @ideal.ideal_new(p, nvars)
  for i = 0; i < nvars; i = i + 1 {
    let e : Array[Int] = Array::make(nvars, 0)
    e[i] = box[i]
    st.add(@poly.mp_monom(p, nvars, e, 1))
  }
  let mons = box_monomials(box)
  for i = 0; i < mons.length(); i = i + 1 {
    if ((mask >> i) & 1) == 1 {
      st.add(@poly.mp_monom(p, nvars, mons[i], 1))
    }
  }
  st
}

test "X2 for finite monomial quotients: F-split <=> reduced <=> perfect <=> rank(F) = dim" {
  let boxes = [[2, 2], [2, 3], [3, 2]]
  for p in [2, 3] {
    for bi = 0; bi < boxes.length(); bi = bi + 1 {
      let box = boxes[bi]
      let nm = box[0] * box[1]
      for mask = 0; mask < (1 << nm); mask = mask + 1 {
        let st = monom_ideal(p, 2, box, mask)
        if !(st.finite) {
          continue
        }
        let r = finite_ring(st)
        let reduced = r.is_reduced()
        check(
          reduced == r.is_reduced_bruteforce(),
          "X2 two notions of reducedness agree, p=" +
            p.to_string() +
            " mask=" +
            mask.to_string(),
        )
        check(
          r.is_f_split() == reduced,
          "X2 F-split <=> reduced, p=" + p.to_string() + " mask=" + mask.to_string(),
        )
        check(
          r.is_perfect() == reduced,
          "X2 perfect <=> reduced, p=" + p.to_string() + " mask=" + mask.to_string(),
        )
        check(
          (@linalg.rank(r.frob, p) == r.dim()) == reduced,
          "X2 rank(F) = dim <=> reduced",
        )
      }
    }
  }
}

test "X5 Frobenius is a ring homomorphism and equals the p-th power map" {
  let boxes = [[2, 2], [2, 3]]
  for p in [2, 3] {
    for bi = 0; bi < boxes.length(); bi = bi + 1 {
      let box = boxes[bi]
      let nm = box[0] * box[1]
      for mask = 0; mask < (1 << nm); mask = mask + 1 {
        let st = monom_ideal(p, 2, box, mask)
        let r = finite_ring(st)
        if @arith.pow_int(p, r.dim()) > 64 {
          continue
        }
        let els = r.elements()
        for i = 0; i < els.length(); i = i + 1 {
          let a = els[i]
          let mut pw = unit(r.dim())
          for _k = 0; _k < p; _k = _k + 1 {
            pw = r.mul_coords(pw, a)
          }
          check(vec_eq(pw, r.frob_coords(a)), "X5 F(a) = a^p")
          for j = 0; j < els.length(); j = j + 1 {
            let b = els[j]
            check(
              vec_eq(
                r.frob_coords(r.mul_coords(a, b)),
                r.mul_coords(r.frob_coords(a), r.frob_coords(b)),
              ),
              "X5 F is multiplicative",
            )
            check(
              vec_eq(
                r.frob_coords(vec_add(a, b, p)),
                vec_add(r.frob_coords(a), r.frob_coords(b), p),
              ),
              "X5 F is additive",
            )
          }
        }
      }
    }
  }
}

fn unit(n : Int) -> Array[Int] {
  let v : Array[Int] = Array::make(n, 0)
  if n > 0 {
    v[0] = 1
  }
  v
}

fn vec_add(a : Array[Int], b : Array[Int], p : Int) -> Array[Int] {
  let out : Array[Int] = []
  for i = 0; i < a.length(); i = i + 1 {
    out.push(@kernel.modpos(a[i] + (if i < b.length() {
      b[i]
    } else {
      0
    }), p))
  }
  out
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

fn mat_eq(a : @linalg.Mat, b : @linalg.Mat) -> Bool {
  if a.length() != b.length() {
    return false
  }
  for i = 0; i < a.length(); i = i + 1 {
    if !(vec_eq(a[i], b[i])) {
      return false
    }
  }
  true
}

test "X1/X3 splittings form a torsor over the R-linear maps annihilating F" {
  // small rings, so every matrix can be enumerated
  let cases : Array[(Int, Array[Int], Int)] = [
    (2, [2, 2], 0), // F_2[x,y]/(x^2,y^2) : not reduced, no splitting
    (2, [2, 2], 2), // F_2[x,y]/(x^2,y^2,x) = F_2[y]/(y^2) : not reduced
    (2, [2, 2], 6), // F_2[x,y]/(x^2,y^2,x,y) = F_2
    (2, [2, 2], 12), // F_2[x,y]/(x^2,y^2,xy) : reduced? (x*y = 0, no nilpotents)
    (3, [2], 0), // F_3[x]/(x^2) : not reduced
    (3, [2], 2), // F_3[x]/(x^2, x) = F_3
    (2, [2], 2), // F_2[x]/(x^2, x) = F_2
    (2, [3], 6), // F_2[x]/(x^3, x^2, x) = F_2
  ]
  for c = 0; c < cases.length(); c = c + 1 {
    let p = cases[c].0
    let box = cases[c].1
    let mask = cases[c].2
    let st = monom_ideal(p, box.length(), box, mask)
    let r = finite_ring(st)
    let dim = r.dim()
    let (sol, nulldim) = r.splitting_system()
    // brute force: every matrix, tested directly against the axioms
    let mut count = 0
    let total = @arith.pow_int(p, dim * dim)
    for k = 0; k < total; k = k + 1 {
      let ds = @arith.digits(k, p)
      let m = @linalg.mat_make(dim, dim)
      for i = 0; i < dim; i = i + 1 {
        for j = 0; j < dim; j = j + 1 {
          m[i][j] = if i * dim + j < ds.length() {
            ds[i * dim + j]
          } else {
            0
          }
        }
      }
      if is_splitting_matrix(r, m) {
        count = count + 1
        // X3: it differs from the particular solution by an R-linear map
        // annihilating F
        match sol {
          Some(xs) => {
            let phi0 = flat_to_mat(xs, dim)
            check(
              is_homogeneous(r, mat_diff(m, phi0, p)),
              "X3 every splitting is phi0 + psi with psi o F = 0",
            )
          }
          None => abort("[fsplit] no particular solution but a splitting exists")
        }
      }
    }
    match sol {
      Some(_) =>
        check(
          count == @arith.pow_int(p, nulldim),
          "X1/X3 #splittings = p^dim(psi space), case " +
            c.to_string() +
            " got " +
            count.to_string() +
            " want " +
            @arith.pow_int(p, nulldim).to_string(),
        )
      None =>
        check(count == 0, "X1 no splitting exists, case " + c.to_string())
    }
  }
}

fn flat_to_mat(xs : Array[Int], dim : Int) -> @linalg.Mat {
  let m = @linalg.mat_make(dim, dim)
  for i = 0; i < dim; i = i + 1 {
    for j = 0; j < dim; j = j + 1 {
      m[i][j] = xs[i * dim + j]
    }
  }
  m
}

fn mat_diff(a : @linalg.Mat, b : @linalg.Mat, p : Int) -> @linalg.Mat {
  let out = @linalg.mat_make(a.length(), a[0].length())
  for i = 0; i < a.length(); i = i + 1 {
    for j = 0; j < a[0].length(); j = j + 1 {
      out[i][j] = @kernel.modpos(a[i][j] - b[i][j], p)
    }
  }
  out
}

///| phi o F = id and phi(r^p s) = r phi(s) for all r, s.
fn is_splitting_matrix(r : FiniteRing, phi : @linalg.Mat) -> Bool {
  let p = r.p
  let dim = r.dim()
  let els = r.elements()
  for i = 0; i < els.length(); i = i + 1 {
    let a = els[i]
    // phi(F(x)) = x
    if !(vec_eq(apply_phi(phi, r.frob_coords(a), p), a)) {
      return false
    }
    for j = 0; j < els.length(); j = j + 1 {
      let b = els[j]
      let lhs = apply_phi(phi, r.mul_coords(r.frob_coords(a), b), p)
      let rhs = r.mul_coords(a, apply_phi(phi, b, p))
      if !(vec_eq(lhs, rhs)) {
        return false
      }
    }
  }
  ignore(dim)
  true
}

///| R-linear and annihilating the image of F.
fn is_homogeneous(r : FiniteRing, psi : @linalg.Mat) -> Bool {
  let p = r.p
  let els = r.elements()
  for i = 0; i < els.length(); i = i + 1 {
    let a = els[i]
    if !(is_zero_vec_local(apply_phi(psi, r.frob_coords(a), p))) {
      return false
    }
    for j = 0; j < els.length(); j = j + 1 {
      let b = els[j]
      let lhs = apply_phi(psi, r.mul_coords(r.frob_coords(a), b), p)
      let rhs = r.mul_coords(a, apply_phi(psi, b, p))
      if !(vec_eq(lhs, rhs)) {
        return false
      }
    }
  }
  true
}

fn is_zero_vec_local(v : Array[Int]) -> Bool {
  for i = 0; i < v.length(); i = i + 1 {
    if v[i] != 0 {
      return false
    }
  }
  true
}

test "X4 the ideals of F_p^n are exactly the coordinate subspaces, all compatibly split" {
  // R = F_p[x]/(x^n - x) = F_p^n for n <= p
  let cases = [(2, 2), (3, 2), (3, 3)]
  for c = 0; c < cases.length(); c = c + 1 {
    let p = cases[c].0
    let n = cases[c].1
    let st = @ideal.ideal_new(p, 1)
    st.add(@poly.mp_make(p, 1, [([n], 1), ([1], p - 1)]))
    let r = finite_ring(st)
    check(r.dim() == n, "X4 dim F_p^n = n, case " + c.to_string())
    // F is the identity on a product of copies of F_p
    let els = r.elements()
    for i = 0; i < els.length(); i = i + 1 {
      check(vec_eq(r.frob_coords(els[i]), els[i]), "X4 F = id on F_p^n")
    }
    let phi = match r.splitting() {
      Some(m) => m
      None => abort("[fsplit] F_p^n must be F-split")
    }
    let idm = @linalg.mat_identity(n)
    check(mat_eq(phi, idm), "X4 the unique splitting of F_p^n is the identity")
    let (_, nulldim) = r.splitting_system()
    check(nulldim == 0, "X4 the splitting is unique")
    // enumerate subspaces as spans of at most two vectors, plus the whole space
    let spans = all_spans(p, r, 2)
    let mut ideals = 0
    for i = 0; i < spans.length(); i = i + 1 {
      let s = spans[i]
      if is_ideal(r, s) {
        ideals = ideals + 1
        // every ideal is compatibly split by phi
        let gens : Array[Array[Int]] = []
        for k = 0; k < s.length(); k = k + 1 {
          gens.push(s[k])
        }
        check(compatibly_split(phi, gens, p), "X4 ideals are compatibly split")
        // and it is generated by an idempotent (basis independent form of
        // "coordinate subspace")
        check(
          generated_by_idempotent(r, s),
          "X4 every ideal of F_p^n is generated by an idempotent",
        )
      }
    }
    check(
      ideals == (1 << n),
      "X4 F_p^n has 2^n ideals, got " + ideals.to_string(),
    )
  }
}

fn all_spans(p : Int, r : FiniteRing, maxgen : Int) -> Array[Array[Array[Int]]] {
  let els = r.elements()
  let dim = r.dim()
  let seen : Array[Array[Array[Int]]] = []
  let keys : Array[String] = []
  let whole = els
  ignore(whole)
  // 0, singles, pairs
  let gens_list : Array[Array[Array[Int]]] = [[]]
  for i = 0; i < els.length(); i = i + 1 {
    gens_list.push([els[i]])
  }
  if maxgen >= 2 {
    for i = 0; i < els.length(); i = i + 1 {
      for j = i + 1; j < els.length(); j = j + 1 {
        gens_list.push([els[i], els[j]])
      }
    }
  }
  if dim >= 3 && maxgen >= 2 {
    // triples are needed to reach the whole space when dim = 3
    for i = 0; i < els.length(); i = i + 1 {
      for j = i + 1; j < els.length(); j = j + 1 {
        for k = j + 1; k < els.length(); k = k + 1 {
          gens_list.push([els[i], els[j], els[k]])
        }
      }
    }
  }
  for g = 0; g < gens_list.length(); g = g + 1 {
    let s = span_of(r, gens_list[g], p)
    let key = span_key(s)
    let mut dup = false
    for t = 0; t < keys.length(); t = t + 1 {
      if keys[t] == key {
        dup = true
      }
    }
    if !(dup) {
      keys.push(key)
      seen.push(s)
    }
  }
  seen
}

fn span_of(r : FiniteRing, gens : Array[Array[Int]], p : Int) -> Array[Array[Int]] {
  let acc : Array[Array[Int]] = []
  acc.push(Array::make(r.dim(), 0))
  let mut changed = true
  while changed {
    changed = false
    for g = 0; g < gens.length(); g = g + 1 {
      for c = 0; c < p; c = c + 1 {
        let cv = scale_vec(gens[g], c, p)
        for i = 0; i < acc.length(); i = i + 1 {
          let s = vec_add(acc[i], cv, p)
          let mut present = false
          for j = 0; j < acc.length(); j = j + 1 {
            if vec_eq(acc[j], s) {
              present = true
            }
          }
          if !(present) {
            acc.push(s)
            changed = true
          }
        }
      }
    }
  }
  acc
}

fn scale_vec(v : Array[Int], c : Int, p : Int) -> Array[Int] {
  let out : Array[Int] = []
  for i = 0; i < v.length(); i = i + 1 {
    out.push(@kernel.modpos(c * v[i], p))
  }
  out
}

fn span_key(s : Array[Array[Int]]) -> String {
  let parts : Array[String] = []
  for i = 0; i < s.length(); i = i + 1 {
    let digits : Array[String] = []
    for j = 0; j < s[i].length(); j = j + 1 {
      digits.push(s[i][j].to_string())
    }
    parts.push(digits.join(""))
  }
  parts.sort_by(fn(a, b) { a.compare(b) })
  parts.join("|")
}

fn is_ideal(r : FiniteRing, s : Array[Array[Int]]) -> Bool {
  let els = r.elements()
  for i = 0; i < s.length(); i = i + 1 {
    for j = 0; j < els.length(); j = j + 1 {
      let prod = r.mul_coords(els[j], s[i])
      let mut present = false
      for k = 0; k < s.length(); k = k + 1 {
        if vec_eq(s[k], prod) {
          present = true
        }
      }
      if !(present) {
        return false
      }
    }
  }
  true
}

///| Is the ideal S equal to R*e for some idempotent e in S?
fn generated_by_idempotent(r : FiniteRing, s : Array[Array[Int]]) -> Bool {
  for i = 0; i < s.length(); i = i + 1 {
    let e = s[i]
    if !(vec_eq(r.mul_coords(e, e), e)) {
      continue
    }
    // R*e must equal S
    let els = r.elements()
    let prod : Array[Array[Int]] = []
    for k = 0; k < els.length(); k = k + 1 {
      let v = r.mul_coords(els[k], e)
      let mut present = false
      for t = 0; t < prod.length(); t = t + 1 {
        if vec_eq(prod[t], v) {
          present = true
        }
      }
      if !(present) {
        prod.push(v)
      }
    }
    if prod.length() != s.length() {
      continue
    }
    let mut same = true
    for t = 0; t < prod.length(); t = t + 1 {
      let mut found = false
      for u = 0; u < s.length(); u = u + 1 {
        if vec_eq(s[u], prod[t]) {
          found = true
        }
      }
      if !(found) {
        same = false
      }
    }
    if same {
      return true
    }
  }
  false
}

test "X6 the splitting of a finite field is the inverse Frobenius" {
  // F_{p^2} = F_p[x]/(irreducible quadratic)
  for p in [2, 3, 5] {
    let f = @ff.fp_find_irreducible(p, 2)
    {
      {
        let st = @ideal.ideal_new(p, 1)
        let terms : Array[(Array[Int], Int)] = []
        for i = 0; i < f.c.length(); i = i + 1 {
          if f.c[i] != 0 {
            terms.push(([i], f.c[i]))
          }
        }
        st.add(@poly.mp_make(p, 1, terms))
        let r = finite_ring(st)
        check(r.dim() == 2, "X6 dim F_{p^2} = 2")
        check(r.is_reduced(), "X6 a field is reduced")
        check(r.is_perfect(), "X6 a finite field is perfect")
        let phi = match r.splitting() {
          Some(m) => m
          None => abort("[fsplit] a field must be F-split")
        }
        // phi o F = id, and on F_{p^2} we have F^2 = id, so phi = F
        check(mat_eq(phi, r.frob), "X6 the splitting of F_{p^2} is F itself")
        let (_, nulldim) = r.splitting_system()
        check(nulldim == 0, "X6 the splitting of a field is unique")
        // phi o F = id verified pointwise
        let els = r.elements()
        for i = 0; i < els.length(); i = i + 1 {
          check(
            vec_eq(apply_phi(phi, r.frob_coords(els[i]), p), els[i]),
            "X6 phi(F(x)) = x",
          )
        }
      }
    }
  }
}
