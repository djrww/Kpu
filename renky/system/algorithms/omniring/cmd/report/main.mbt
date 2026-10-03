// ===========================================================================
//  omniring · report — prints the public ledger of the axiom book:
//  the closed loop over F_p for p = 2, 3, 5 and the class-group table.
//
//      moon run ky678/renky/system/algorithms/omniring/cmd/report
// ===========================================================================

fn main {
  println("==================================================================")
  println(" omniring — 公理書 ledger (definitions / formulas / self-proof)")
  println(" no third-party dependencies; algorithms A01..A5x; 122 executable laws")
  println("==================================================================")
  for p in [2, 3] {
    let rep = @loop.run_all(p)
    println(rep.to_string())
    println("")
  }
  println("------------------------------------------------------------------")
  println(" ideal class groups (K3) and their correspondence rings (K7)")
  println("------------------------------------------------------------------")
  for d in [-3, -4, -7, -8, -11, -15, -20, -23, -24, -31, -40, -47, -71, -163] {
    let cl = @classgroup.reduced_forms(d)
    let parts : Array[String] = []
    for i = 0; i < cl.length(); i = i + 1 {
      parts.push(cl[i].to_string())
    }
    println(
      "  h(" +
        d.to_string() +
        ") = " +
        cl.length().to_string() +
        "   " +
        parts.join("  ") +
        "   genera = " +
        @classgroup.genus_count(d).to_string(),
    )
  }
  println("------------------------------------------------------------------")
  println(" F-split sample (X2): reduced <=> F-split <=> perfect")
  println("------------------------------------------------------------------")
  for p in [2, 3] {
    let cases : Array[(String, Int, Array[Int])] = [
      ("F_p[x]/(x^2)", p, [2]),
      ("F_p[x]/(x^2-x)", p, [2]),
      ("F_p[x,y]/(x^2,y^2)", p, [2, 2]),
      ("F_p[x,y]/(x^2,y^2,xy)", p, [2, 2]),
    ]
    for c = 0; c < cases.length(); c = c + 1 {
      let (name, pp, box) = cases[c]
      let st = build_case(pp, box, name)
      if !(st.finite) {
        println("  " + name + " : infinite")
        continue
      }
      let r = @fsplit.finite_ring(st)
      println(
        "  p=" +
          pp.to_string() +
          " " +
          name +
          " : dim=" +
          r.dim().to_string() +
          " reduced=" +
          r.is_reduced().to_string() +
          " F-split=" +
          r.is_f_split().to_string() +
          " perfect=" +
          r.is_perfect().to_string(),
      )
    }
  }
  println("==================================================================")
}

fn build_case(p : Int, box : Array[Int], name : String) -> @ideal.IdealState {
  let n = box.length()
  let st = @ideal.ideal_new(p, n)
  if name.has_prefix("F_p[x]/(x^2)") {
    st.add(@poly.mp_monom(p, 1, [2], 1))
  } else if name.has_prefix("F_p[x]/(x^2-x)") {
    st.add(@poly.mp_make(p, 1, [([2], 1), ([1], p - 1)]))
  } else if name.has_prefix("F_p[x,y]/(x^2,y^2,xy)") {
    st.add(@poly.mp_monom(p, 2, [2, 0], 1))
    st.add(@poly.mp_monom(p, 2, [0, 2], 1))
    st.add(@poly.mp_monom(p, 2, [1, 1], 1))
  } else {
    for i = 0; i < n; i = i + 1 {
      let e : Array[Int] = Array::make(n, 0)
      e[i] = box[i]
      st.add(@poly.mp_monom(p, n, e, 1))
    }
  }
  st
}
