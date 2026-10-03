// ===========================================================================
//  homology — a generic chain-complex engine over Z and over F_p
//  (「同調」: the second half of the closed loop).
//
//  Definition (chain complex).  A chain complex of free modules is a sequence
//      ... -> C_{i+1} -d_{i+1}-> C_i -d_i-> C_{i-1} -> ...
//  with d_i * d_{i+1} = 0.  Here C_i = R^{ranks[i]} with R = Z (p = 0) or
//  R = F_p, and `diffs[i]` is the matrix of d_{i+1} : C_{i+1} -> C_i.
//
//  Definition (homology).  H_i = ker(d_i) / im(d_{i+1}).
//
//  Theorem (computability, used by the engine).
//   * Over a field F_p: dim ker d = n - rank d (A27/A28), so
//         dim H_i = (ranks[i] - rank d_i) - rank d_{i+1}.
//   * Over Z: ker(d_i) is a direct summand of C_i and, since
//     im(d_{i+1}) ⊆ ker(d_i), a Smith normal form (A30) of d_{i+1} gives
//         H_i = Z^{ranks[i] - rank d_i - rank d_{i+1}}  (+)  (+)_j Z/s_j
//     where the s_j > 1 are the invariant factors of d_{i+1}.
//
//  Theorems (verified exhaustively in `homology_wbtest`):
//   H1 Euler characteristic: chi(C) = sum (-1)^i ranks[i] = sum (-1)^i rank H_i.
//   H2 Koszul homology of x on R = F_p[x]/(x^n):  H_0 = R/(x) = F_p and
//      H_1 = ann_R(x) = (x^{n-1}) = F_p.
//   H3 simplicial homology of a graph: beta_0 = #components and
//      beta_1 = E - V + #components = the cycle rank (A51 dual-graph genus).
//   H4 Kunneth over a field: beta_k(A (x) B) = sum_{i+j=k} beta_i(A) beta_j(B).
//   H5 mapping cone (A49): d^2 = 0 on cone(f); cone(id) is acyclic;
//      cone(0) = B (+) A[-1] on homology; and for every chain map the
//      Betti numbers satisfy the long-exact-sequence alternating sum.
//   H6 torsion over Z is read off the Smith normal form, and chi is additive in
//      short exact sequences of complexes.
// ===========================================================================

pub type Mat = @linalg.Mat

///| Homology of one degree: free rank and torsion invariant factors.
pub struct HGroup {
  betti : Int
  torsion : Array[Int]
}

pub fn HGroup::to_string(self : HGroup) -> String {
  let ts : Array[String] = []
  for i = 0; i < self.torsion.length(); i = i + 1 {
    ts.push("Z/" + self.torsion[i].to_string())
  }
  "Z^" + self.betti.to_string() + (if ts.is_empty() {
    ""
  } else {
    " + " + ts.join(" + ")
  })
}

///| A chain complex: `p = 0` means coefficients in Z, otherwise F_p.
pub struct Complex {
  p : Int
  ranks : Array[Int]
  diffs : Array[Mat]
}

///| Build a complex, rejecting it unless d*d = 0.
pub fn complex_new(p : Int, ranks : Array[Int], diffs : Array[Mat]) -> Complex {
  if diffs.length() + 1 != ranks.length() {
    abort("[homology] a complex on n modules has n-1 differentials")
  }
  for i = 0; i < diffs.length(); i = i + 1 {
    if @linalg.mat_rows(diffs[i]) != ranks[i] {
      abort("[homology] differential shape mismatch (rows)")
    }
    if @linalg.mat_cols(diffs[i]) != ranks[i + 1] {
      abort("[homology] differential shape mismatch (cols)")
    }
  }
  let c : Complex = { p, ranks, diffs }
  if !(c.is_valid()) {
    abort("[homology] d*d != 0: not a chain complex")
  }
  c
}

///| Integer matrix product (coefficients in Z), used when p = 0.
fn mat_mul_int(a : Mat, b : Mat) -> Mat {
  let n = @linalg.mat_rows(a)
  let k = @linalg.mat_cols(a)
  let m = @linalg.mat_cols(b)
  let out = @linalg.mat_make(n, m)
  for i = 0; i < n; i = i + 1 {
    for j = 0; j < m; j = j + 1 {
      let mut acc = 0
      for t = 0; t < k; t = t + 1 {
        acc = acc + a[i][t] * b[t][j]
      }
      out[i][j] = acc
    }
  }
  out
}

fn mat_mul_any(a : Mat, b : Mat, p : Int) -> Mat {
  if p == 0 {
    mat_mul_int(a, b)
  } else {
    @linalg.mat_mul(a, b, p)
  }
}

///| Is d_i * d_{i+1} = 0 for every i?
pub fn Complex::is_valid(self : Complex) -> Bool {
  for i = 0; i + 1 < self.diffs.length(); i = i + 1 {
    let prod = mat_mul_any(self.diffs[i], self.diffs[i + 1], self.p)
    for r = 0; r < @linalg.mat_rows(prod); r = r + 1 {
      for c = 0; c < @linalg.mat_cols(prod); c = c + 1 {
        let v = prod[r][c]
        if self.p == 0 {
          if v != 0 {
            return false
          }
        } else if v % self.p != 0 {
          return false
        }
      }
    }
  }
  true
}

///| Rank of a differential (over F_p or over Z).
fn rank_of(m : Mat, p : Int) -> Int {
  if p == 0 {
    @linalg.rank_z(m)
  } else {
    @linalg.rank(m, p)
  }
}

///| **A72 chain-complex homology** in every degree (Gaussian elimination over F_p,
///| Smith normal form over Z).
pub fn Complex::homology(self : Complex) -> Array[HGroup] {
  let out : Array[HGroup] = []
  for i = 0; i < self.ranks.length(); i = i + 1 {
    let r_out = if i > 0 {
      rank_of(self.diffs[i - 1], self.p)
    } else {
      0
    }
    let r_in = if i < self.diffs.length() {
      rank_of(self.diffs[i], self.p)
    } else {
      0
    }
    let betti = self.ranks[i] - r_out - r_in
    let torsion : Array[Int] = []
    if self.p == 0 && i < self.diffs.length() {
      let inv = @linalg.snf(self.diffs[i])
      for k = 0; k < inv.length(); k = k + 1 {
        if inv[k] > 1 {
          torsion.push(inv[k])
        }
      }
    }
    out.push({ betti, torsion })
  }
  out
}

///| Betti numbers.
pub fn Complex::betti(self : Complex) -> Array[Int] {
  let h = self.homology()
  let out : Array[Int] = []
  for i = 0; i < h.length(); i = i + 1 {
    out.push(h[i].betti)
  }
  out
}

///| chi(C) = sum (-1)^i rank C_i.
pub fn Complex::euler(self : Complex) -> Int {
  let mut acc = 0
  for i = 0; i < self.ranks.length(); i = i + 1 {
    acc = if i % 2 == 0 {
      acc + self.ranks[i]
    } else {
      acc - self.ranks[i]
    }
  }
  acc
}

///| chi(H) = sum (-1)^i rank H_i (torsion does not contribute).
pub fn Complex::homology_euler(self : Complex) -> Int {
  let h = self.homology()
  let mut acc = 0
  for i = 0; i < h.length(); i = i + 1 {
    acc = if i % 2 == 0 {
      acc + h[i].betti
    } else {
      acc - h[i].betti
    }
  }
  acc
}

///| Total torsion size (product of invariant factors), for reporting.
pub fn Complex::torsion_order(self : Complex) -> Int {
  let h = self.homology()
  let mut acc = 1
  for i = 0; i < h.length(); i = i + 1 {
    for k = 0; k < h[i].torsion.length(); k = k + 1 {
      acc = acc * h[i].torsion[k]
    }
  }
  acc
}

// ---------------------------------------------------------------------------
//  Constructions: direct sum, tensor product, mapping cone, Koszul, graphs
// ---------------------------------------------------------------------------

///| Direct sum of complexes (used for split short exact sequences).
pub fn direct_sum(a : Complex, b : Complex) -> Complex {
  if a.p != b.p {
    abort("[homology.direct_sum] coefficient mismatch")
  }
  let ranks : Array[Int] = []
  for i = 0; i < a.ranks.length(); i = i + 1 {
    ranks.push(a.ranks[i] + b.ranks[i])
  }
  let diffs : Array[Mat] = []
  for i = 0; i < a.diffs.length(); i = i + 1 {
    diffs.push(block_diag(a.diffs[i], b.diffs[i]))
  }
  complex_new(a.p, ranks, diffs)
}

fn block_diag(a : Mat, b : Mat) -> Mat {
  let ra = @linalg.mat_rows(a)
  let ca = @linalg.mat_cols(a)
  let rb = @linalg.mat_rows(b)
  let cb = @linalg.mat_cols(b)
  let out = @linalg.mat_make(ra + rb, ca + cb)
  for i = 0; i < ra; i = i + 1 {
    for j = 0; j < ca; j = j + 1 {
      out[i][j] = a[i][j]
    }
  }
  for i = 0; i < rb; i = i + 1 {
    for j = 0; j < cb; j = j + 1 {
      out[ra + i][ca + j] = b[i][j]
    }
  }
  out
}

///| **A73/H4 tensor product of complexes** over a field: (A (x) B)_n = (+)_{i+j=n}
///| A_i (x) B_j with the sign-twisted differential.
pub fn tensor_product(a : Complex, b : Complex) -> Complex {
  if a.p == 0 || a.p != b.p {
    abort("[homology.tensor_product] field coefficients only")
  }
  let p = a.p
  let n = a.ranks.length() + b.ranks.length() - 1
  // index list of the summands A_i (x) B_j of total degree k
  let idx : Array[Array[(Int, Int)]] = []
  for k = 0; k < n; k = k + 1 {
    let row : Array[(Int, Int)] = []
    for i = 0; i < a.ranks.length(); i = i + 1 {
      let j = k - i
      if j >= 0 && j < b.ranks.length() {
        row.push((i, j))
      }
    }
    idx.push(row)
  }
  let ranks : Array[Int] = []
  for k = 0; k < n; k = k + 1 {
    let mut r = 0
    for t = 0; t < idx[k].length(); t = t + 1 {
      r = r + a.ranks[idx[k][t].0] * b.ranks[idx[k][t].1]
    }
    ranks.push(r)
  }
  // offsets of each summand inside degree k
  let offs : Array[Array[Int]] = []
  for k = 0; k < n; k = k + 1 {
    let row : Array[Int] = []
    let mut acc = 0
    for t = 0; t < idx[k].length(); t = t + 1 {
      row.push(acc)
      acc = acc + a.ranks[idx[k][t].0] * b.ranks[idx[k][t].1]
    }
    offs.push(row)
  }
  let diffs : Array[Mat] = []
  for k = 0; k + 1 < n; k = k + 1 {
    // d : degree k+1 -> degree k
    let m = @linalg.mat_make(ranks[k], ranks[k + 1])
    for t = 0; t < idx[k + 1].length(); t = t + 1 {
      let (i, j) = idx[k + 1][t]
      let src = offs[k + 1][t]
      // d_A (x) 1 : A_i (x) B_j -> A_{i-1} (x) B_j
      if i > 0 {
        let dst = offset_of(idx[k], offs[k], (i - 1, j))
        let da = a.diffs[i - 1]
        for r = 0; r < a.ranks[i - 1]; r = r + 1 {
          for c = 0; c < a.ranks[i]; c = c + 1 {
            for s = 0; s < b.ranks[j]; s = s + 1 {
              m[dst + r * b.ranks[j] + s][src + c * b.ranks[j] + s] = @kernel.modpos(
                da[r][c],
                p,
              )
            }
          }
        }
      }
      // (-1)^i * 1 (x) d_B : A_i (x) B_j -> A_i (x) B_{j-1}
      if j > 0 {
        let dst = offset_of(idx[k], offs[k], (i, j - 1))
        let db = b.diffs[j - 1]
        let sign = if i % 2 == 0 {
          1
        } else {
          -1
        }
        // rows of db index B_{j-1}, columns index B_j; the A index s is fixed
        for r = 0; r < b.ranks[j - 1]; r = r + 1 {
          for c = 0; c < b.ranks[j]; c = c + 1 {
            for s = 0; s < a.ranks[i]; s = s + 1 {
              m[dst + s * b.ranks[j - 1] + r][src + s * b.ranks[j] + c] = @kernel.modpos(
                sign * db[r][c],
                p,
              )
            }
          }
        }
      }
    }
    diffs.push(m)
  }
  complex_new(p, ranks, diffs)
}

fn offset_of(idx : Array[(Int, Int)], offs : Array[Int], key : (Int, Int)) -> Int {
  for t = 0; t < idx.length(); t = t + 1 {
    if idx[t].0 == key.0 && idx[t].1 == key.1 {
      return offs[t]
    }
  }
  abort("[homology] summand not found")
}

///| A chain map f : A -> B, given by one matrix per degree.
pub struct ChainMap {
  a : Complex
  b : Complex
  maps : Array[Mat]
}

///| Build a chain map, validating d_B f = f d_A.
pub fn chain_map_new(a : Complex, b : Complex, maps : Array[Mat]) -> ChainMap {
  let f : ChainMap = { a, b, maps }
  if !(f.is_valid()) {
    abort("[homology] not a chain map: d_B f != f d_A")
  }
  f
}

///| Validate the chain-map condition d_B f = f d_A.
pub fn ChainMap::is_valid(self : ChainMap) -> Bool {
  for i = 0; i + 1 < self.maps.length(); i = i + 1 {
    let lhs = mat_mul_any(self.b.diffs[i], self.maps[i + 1], self.a.p)
    let rhs = mat_mul_any(self.maps[i], self.a.diffs[i], self.a.p)
    for r = 0; r < @linalg.mat_rows(lhs); r = r + 1 {
      for c = 0; c < @linalg.mat_cols(lhs); c = c + 1 {
        let d = @kernel.modpos(lhs[r][c] - rhs[r][c], self.a.p)
        if self.a.p == 0 {
          if lhs[r][c] != rhs[r][c] {
            return false
          }
        } else if d != 0 {
          return false
        }
      }
    }
  }
  true
}

///| **A74 mapping cone**: cone(f)_n = B_n (+) A_{n-1} with
///| d(b, a) = (d_B b + f(a), -d_A a).
pub fn mapping_cone(f : ChainMap) -> Complex {
  let a = f.a
  let b = f.b
  if a.p != b.p {
    abort("[homology.mapping_cone] coefficient mismatch")
  }
  let n = Int::max(a.ranks.length(), b.ranks.length()) + 1
  let ranks : Array[Int] = []
  for k = 0; k < n; k = k + 1 {
    let rb = if k < b.ranks.length() {
      b.ranks[k]
    } else {
      0
    }
    let ra = if k - 1 >= 0 && k - 1 < a.ranks.length() {
      a.ranks[k - 1]
    } else {
      0
    }
    ranks.push(rb + ra)
  }
  let diffs : Array[Mat] = []
  for k = 0; k + 1 < n; k = k + 1 {
    let rows = ranks[k]
    let cols = ranks[k + 1]
    let m = @linalg.mat_make(rows, cols)
    let rb_k = if k < b.ranks.length() {
      b.ranks[k]
    } else {
      0
    }
    let rb_k1 = if k + 1 < b.ranks.length() {
      b.ranks[k + 1]
    } else {
      0
    }
    let ra_k = if k < a.ranks.length() && k - 1 >= 0 {
      a.ranks[k - 1]
    } else {
      0
    }
    let ra_k1 = if k < a.ranks.length() {
      a.ranks[k]
    } else {
      0
    }
    ignore(rb_k1)
    // d_B block
    if k < b.diffs.length() {
      let db = b.diffs[k]
      for r = 0; r < @linalg.mat_rows(db); r = r + 1 {
        for c = 0; c < @linalg.mat_cols(db); c = c + 1 {
          m[r][c] = db[r][c]
        }
      }
    }
    // -d_A block
    if k - 1 >= 0 && k - 1 < a.diffs.length() {
      let da = a.diffs[k - 1]
      for r = 0; r < @linalg.mat_rows(da); r = r + 1 {
        for c = 0; c < @linalg.mat_cols(da); c = c + 1 {
          m[rb_k + r][rb_k1 + c] = -da[r][c]
        }
      }
    }
    // f block: A_k -> B_k (the a-component of cone_{k+1} = B_{k+1} (+) A_k)
    if k < f.maps.length() {
      let fm = f.maps[k]
      for r = 0; r < @linalg.mat_rows(fm); r = r + 1 {
        for c = 0; c < @linalg.mat_cols(fm); c = c + 1 {
          if r < rb_k && rb_k1 + c < cols {
            m[r][rb_k1 + c] = fm[r][c]
          }
        }
      }
    }
    ignore(ra_k)
    ignore(ra_k1)
    diffs.push(m)
  }
  complex_new(a.p, ranks, diffs)
}

///| **A75 Koszul complex** K(x) of multiplication by `m` on R^n over F_p:
///|     0 -> R^n -m-> R^n -> 0
///| whose homology is H_1 = ker m and H_0 = coker m.
pub fn koszul_two_term(m : Mat, p : Int) -> Complex {
  complex_new(p, [m.length(), m[0].length()], [m])
}

///| Simplicial chain complex of a graph (0-simplices = vertices, 1-simplices =
///| edges) over F_p; H_0 = components, H_1 = cycle space.
pub fn graph_complex(nverts : Int, edges : Array[(Int, Int)], p : Int) -> Complex {
  let d = @linalg.mat_make(nverts, edges.length())
  for e = 0; e < edges.length(); e = e + 1 {
    d[edges[e].0][e] = @kernel.modpos(-1, p)
    d[edges[e].1][e] = 1
  }
  complex_new(p, [nverts, edges.length()], [d])
}

///| **A76 cycle rank / genus** of a graph: E - V + #components.
pub fn cycle_rank(nverts : Int, edges : Array[(Int, Int)]) -> Int {
  let comps = components(nverts, edges)
  edges.length() - nverts + comps
}

///| Number of connected components, by union-find.
pub fn components(nverts : Int, edges : Array[(Int, Int)]) -> Int {
  let parent : Array[Int] = []
  for i = 0; i < nverts; i = i + 1 {
    parent.push(i)
  }
  for e = 0; e < edges.length(); e = e + 1 {
    let a = find(parent, edges[e].0)
    let b = find(parent, edges[e].1)
    if a != b {
      parent[a] = b
    }
  }
  let mut cnt = 0
  for i = 0; i < nverts; i = i + 1 {
    if find(parent, i) == i {
      cnt = cnt + 1
    }
  }
  cnt
}

fn find(parent : Array[Int], x : Int) -> Int {
  let mut r = x
  while parent[r] != r {
    r = parent[r]
  }
  r
}
