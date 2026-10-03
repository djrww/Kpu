// ---------------------------------------------------------------------------
// A06  argmin      verified arg-minimum over a fixed array
// A07  bsearch     verified binary search
//
// `argmin` is the primitive on which the tropical engine builds its
// min-plus determinant (assignment problem): the tropical determinant is a
// minimum over permutations, and every greedy/shortest-path step is an
// `argmin`.  Having it machine-checked means the whole tropical half of the
// loop rests on a proved selector.
// ---------------------------------------------------------------------------

///| Index of a minimal entry.
///|
///| Theorem (proved).  If `xs.length() ≥ 1` then `argmin xs` is a valid index
///| and `xs[argmin xs] ≤ xs[j]` for every valid index `j`.
pub fn argmin(xs : FixedArray[Int]) -> Int where {
  proof_require: 0 < xs.length(),
  proof_ensure: result => argmin_ok(xs, result),
} {
  let mut best = 0
  for i = 1; i < xs.length(); i = i + 1 {
    if xs[i] < xs[best] {
      proof_assert xs[i] <= xs[best]
      best = i
    }
  } where {
    proof_invariant: 1 <= i,
    proof_invariant: 0 <= best,
    proof_invariant: best < xs.length(),
    proof_invariant: prefix_min(xs, best, i),
    proof_reasoning: (
      #| `best` always points inside the array and dominates the scanned
      #| prefix `[0, i)`.  Scanning position `i` either keeps `best` (then the
      #| prefix grows by an entry that is not smaller) or replaces it by `i`
      #| (then the new entry is ≤ every previous one).  At exit `i = length`,
      #| so the prefix is the whole array.
    ),
  }
  best
}

///| Binary search on a weakly increasing array.
///|
///| Theorem (proved).  `Some(i)` is returned iff `xs[i] == key`; `None` is
///| returned only when `key` does not occur in `xs`.
pub fn bsearch(xs : FixedArray[Int], key : Int) -> Int? where {
  proof_require: sorted(xs),
  proof_ensure: result => bsearch_ok(xs, key, result),
} {
  for lo = 0, hi = xs.length(); lo < hi; {
    let mid = lo + (hi - lo) / 2
    if xs[mid] < key {
      proof_assert ∀ j : Int, (0 <= j) && (j < mid + 1) → xs[j] < key
      continue mid + 1, hi
    } else if key < xs[mid] {
      proof_assert ∀ j : Int, (mid <= j) && (j < xs.length()) → key < xs[j]
      continue lo, mid
    } else {
      proof_assert xs[mid] == key
      break Some(mid)
    }
  } nobreak {
    None
  } where {
    proof_invariant: 0 <= lo,
    proof_invariant: lo <= hi,
    proof_invariant: hi <= xs.length(),
    proof_invariant: ∀ j : Int, (0 <= j) && (j < lo) → xs[j] < key,
    proof_invariant: ∀ j : Int, (hi <= j) && (j < xs.length()) → key < xs[j],
    proof_yield: res => bsearch_ok(xs, key, res),
    proof_reasoning: (
      #| The window `[lo, hi)` shrinks; everything left of it is `< key` and
      #| everything from `hi` on is `> key`.  An empty window therefore proves
      #| absence, and `xs[mid] == key` proves presence.
    ),
  }
}
