///|
// 註冊表義務自證之總驗收

test "registry: 演算法數目多於 80" {
  assert_true(registry().length() > 80)
}

test "registry: 全部自證通過" {
  let (passed, total, fails) = registry_run()
  if fails.length() > 0 {
    fail("失敗：\{fails}")
  }
  assert_eq(passed, total)
}

test "registry: 每個條目有唯一編號與名稱" {
  let entries = registry()
  let ids : Array[String] = []
  for e in entries {
    let mut dup = false
    for x in ids {
      if x == e.id {
        dup = true
      }
    }
    assert_false(dup)
    ids.push(e.id)
    assert_true(e.name.length() > 0)
  }
}
