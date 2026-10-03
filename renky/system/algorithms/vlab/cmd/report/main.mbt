///|
// 自證執行器：運行全部演算法義務驗證並輸出報告
fn main {
  println("========== 向量空間與內積空間演算法庫：義務自證 ==========")
  let entries = @vlab.registry()
  println("演算法總數：\{entries.length()} 條")
  println("----------------------------------------------------------")
  let (passed, total, fails) = @vlab.registry_run()
  for e in entries {
    let status = if (e.check)() { "PASS" } else { "FAIL" }
    println("  [\{status}] \{e.id}  \{e.name}")
  }
  println("----------------------------------------------------------")
  println("結果：\{passed} / \{total} 通過")
  if fails.length() > 0 {
    println("失敗項目：")
    for f in fails {
      println("  ✗ \{f}")
    }
  }
  if passed == total {
    println("全部演算法自證通過。")
  }
}
