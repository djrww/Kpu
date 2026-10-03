#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
自研突變測試器（零第三方依賴）
=================================
對 *.mbt 源碼（排除 *_test.mbt / *_wbtest.mbt）做語法層突變：
  1. 算術運算子  + ↔ −,  * ↔ /
  2. 比較運算子  < → <=,  <= → <,  > → >=,  >= → >
  3. 邏輯運算子  && ↔ ||
  4. 數值字面量  0.0 ↔ 1.0,  2.0 → 1.0,  0.5 → 1.5,  整數邊界 ±1
每個突變體：
  moon check 失敗 → UNVIABLE（編譯器擊殺）
  moon test  失敗 → KILLED（測試擊殺）
  兩者皆過        → SURVIVED（倖存）
擊殺率 = KILLED / (KILLED + SURVIVED)
"""
import json
import os
import random
import re
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MOON = os.path.expanduser("~/.moon/bin/moon")
SEED = 20261003
MAX_MUTANTS = int(os.environ.get("MAX_MUTANTS", "400"))

# ---------------------------------------------------------------- tokenizer
TOKEN_RE = re.compile(
    r"""
    (?P<comment>//[^\n]*|/\*.*?\*/)
  | (?P<string>"(?:[^"\\]|\\.)*")
  | (?P<char>b?'(?:[^'\\]|\\.)*')
  | (?P<float>\d+\.\d+(?:[eE][+-]?\d+)?)
  | (?P<int>\b\d+\b)
  | (?P<op><=|>=|==|!=|&&|\|\||->|=>|[<>+\-*/])
  | (?P<other>.)
    """,
    re.S | re.X,
)


def tokenize(src):
    toks = []
    for m in TOKEN_RE.finditer(src):
        kind = m.lastgroup
        toks.append((kind, m.group(0), m.start(), m.end()))
    return toks


def is_code(tok):
    return tok[0] not in ("comment", "string", "char")


# ------------------------------------------------------------- mutation rules
def candidates_for(kind, text, prev_code, next_code):
    """Yield replacement strings for a token."""
    if kind == "op":
        if text == "+":
            yield "-"
        elif text == "-":
            # avoid unary minus / arrows: only mutate if preceded by
            # an identifier/closing token or a digit
            if prev_code and re.search(r"[\w)\]\d]$", prev_code):
                yield "+"
        elif text == "*":
            yield "/"
        elif text == "/":
            yield "*"
        elif text == "<":
            yield "<="
        elif text == "<=":
            yield "<"
        elif text == ">":
            yield ">="
        elif text == ">=":
            yield ">"
        elif text == "&&":
            yield "||"
        elif text == "||":
            yield "&&"
    elif kind == "float":
        if text == "0.0":
            yield "1.0"
        elif text == "1.0":
            yield "0.0"
        elif text == "2.0":
            yield "1.0"
        elif text == "0.5":
            yield "1.5"
    elif kind == "int":
        if text == "0":
            yield "1"
        elif text == "1":
            yield "0"
        elif text == "2":
            yield "3"
        elif text == "3":
            yield "2"


def collect_mutants(path):
    src = open(path, encoding="utf-8").read()
    toks = tokenize(src)
    out = []
    code_texts = [t[1] for t in toks if is_code(t)]
    ci = 0
    for t in toks:
        if is_code(t):
            prev_code = "".join(code_texts[max(0, ci - 3):ci])
            next_code = "".join(code_texts[ci + 1:ci + 3])
            for rep in candidates_for(t[0], t[1], prev_code, next_code):
                out.append((t[2], t[3], rep))
            ci += 1
    return out


def apply_mutant(path, mut):
    src = open(path, encoding="utf-8").read()
    start, end, rep = mut
    return src[:start] + rep + src[end:]


# ------------------------------------------------------------------ runner
def run(cmd, timeout=180):
    try:
        p = subprocess.run(
            cmd, cwd=ROOT, capture_output=True, timeout=timeout, text=True
        )
        return p.returncode, p.stdout + p.stderr
    except subprocess.TimeoutExpired:
        return 124, "TIMEOUT"


def main():
    random.seed(SEED)
    # registry.mbt 之 check 閉包屬「測試預言（oracle）」，按突變測試慣例排除
    targets = sorted(
        f
        for f in os.listdir(ROOT)
        if f.endswith(".mbt")
        and not f.endswith("_test.mbt")
        and not f.endswith("_wbtest.mbt")
        and f != "registry.mbt"
    )
    all_muts = []
    for f in targets:
        path = os.path.join(ROOT, f)
        for m in collect_mutants(path):
            all_muts.append((f, m))
    print(f"[i] {len(targets)} files, {len(all_muts)} raw mutants")
    if len(all_muts) > MAX_MUTANTS:
        all_muts = random.sample(all_muts, MAX_MUTANTS)
    print(f"[i] sampled {len(all_muts)} mutants (seed {SEED})")

    # baseline sanity
    rc, _ = run([MOON, "check"])
    assert rc == 0, "baseline moon check failed"
    rc, _ = run([MOON, "test"])
    assert rc == 0, "baseline moon test failed"
    print("[i] baseline green")

    results = []
    counts = {"KILLED": 0, "SURVIVED": 0, "UNVIABLE": 0, "TIMEOUT": 0}
    t0 = time.time()
    for idx, (f, mut) in enumerate(all_muts):
        path = os.path.join(ROOT, f)
        backup = open(path, encoding="utf-8").read()
        mutated = apply_mutant(path, mut)
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(mutated)
        try:
            rc, out = run([MOON, "check"], timeout=90)
            if rc != 0:
                status = "UNVIABLE"
            else:
                rc2, out2 = run([MOON, "test"], timeout=45)
                if rc2 == 124:
                    status = "TIMEOUT"
                elif rc2 != 0:
                    status = "KILLED"
                else:
                    status = "SURVIVED"
        finally:
            with open(path, "w", encoding="utf-8") as fh:
                fh.write(backup)
        counts[status] += 1
        results.append(
            {
                "file": f,
                "pos": mut[0],
                "repl": mut[2],
                "status": status,
            }
        )
        done = idx + 1
        viable = counts["KILLED"] + counts["SURVIVED"]
        rate = counts["KILLED"] / viable if viable else 0.0
        elapsed = time.time() - t0
        print(
            f"[{done}/{len(all_muts)}] {status:9s} {f}:{mut[0]} -> {mut[2]!r}"
            f"  (kill {rate:.1%}, {elapsed:.0f}s)",
            flush=True,
        )

    viable = counts["KILLED"] + counts["SURVIVED"]
    summary = {
        "total_mutants": len(all_muts),
        "killed": counts["KILLED"],
        "survived": counts["SURVIVED"],
        "unviable": counts["UNVIABLE"],
        "timeout": counts["TIMEOUT"],
        "kill_rate_viable": counts["KILLED"] / viable if viable else 0.0,
        "kill_rate_all": counts["KILLED"] / len(all_muts),
        "elapsed_s": time.time() - t0,
    }
    with open(os.path.join(ROOT, "tools", "mutation_report.json"), "w") as fh:
        json.dump({"summary": summary, "results": results}, fh, indent=2)
    print("\n===== SUMMARY =====")
    for k, v in summary.items():
        print(f"  {k}: {v:.4f}" if isinstance(v, float) else f"  {k}: {v}")
    surv = [r for r in results if r["status"] == "SURVIVED"]
    print(f"\n倖存突變體 ({len(surv)}):")
    for r in surv[:80]:
        print(f"  {r['file']}:{r['pos']} -> {r['repl']!r}")


if __name__ == "__main__":
    main()
