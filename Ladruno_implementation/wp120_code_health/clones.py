#!/usr/bin/env python3
"""WP-120 token-window clone finder (dependency-free stand-in for PMD CPD / jscpd).

Method (CPD's, re-implemented):
  * each file is cleaned by the quirk lint's `clean()` (comments, string contents,
    `#if 0` and preprocessor lines removed; line numbers kept) and tokenized;
  * a window of --min-tokens consecutive tokens (default 100, CPD's default) is
    rolling-hashed; windows seen in >= 2 places are clone seeds;
  * each seed pair is extended right to its maximal length (a pair whose left
    neighbours also match is interior to a longer pair and skipped); every pair
    is re-verified token by token, so hash collisions cannot count;
  * `--mode exact` compares raw tokens (CPD default, "type-1");
    `--mode norm` replaces identifiers and literals by placeholders, keeping
    keywords/punctuation ("type-2": renamed copies, CPD --ignore-identifiers
    --ignore-literals).

A line is "duplicated" if any token on it lies inside a clone.  Percent =
duplicated code lines / code lines (non-blank after cleaning).

    python clones.py                    # fork-stamped files, exact
    python clones.py --mode norm        # renamed-identifier clones
    python clones.py --scope cross      # fork files vs vanilla (copied-from-upstream)
    python clones.py --json out.json
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (KEYWORDS, clean, code_lines, functions, read, rel,  # noqa: E402
                    split_fork_vanilla, tokens)

MOD = (1 << 61) - 1
BASE = 1_000_003
MAX_BUCKET = 40          # a window repeated > 40 times is a micro-idiom, not a copy


class F:
    __slots__ = ("path", "rel", "cl", "toks", "lines", "ids", "fn_of_line", "code", "fork")


def load(paths, mode, fork_flag, vocab):
    out = []
    for p in paths:
        f = F()
        f.path, f.rel, f.fork = p, rel(p), fork_flag
        f.cl = clean(read(p))
        tk = tokens(f.cl)
        f.lines = [li for _, li in tk]
        if mode == "norm":
            seq = [t if (t in KEYWORDS or not (t[0].isalnum() or t[0] == "_")) else
                   ("#" if t[0].isdigit() else "$") for t, _ in tk]
        else:
            seq = [t for t, _ in tk]
        f.toks = [vocab.setdefault(t, len(vocab) + 1) for t in seq]
        f.code = code_lines(f.cl)
        fn = [None] * len(f.cl)
        for fu in functions(f.cl):
            for li in range(fu.start, fu.end + 1):
                if fn[li] is None or fn[li].count("::") < fu.name.count("::"):
                    fn[li] = fu.name
        f.fn_of_line = fn
        out.append(f)
    return out


def _hashes(t, W, powW):
    if len(t) < W:
        return
    h = 0
    for k in range(W):
        h = (h * BASE + t[k]) % MOD
    yield 0, h
    for k in range(W, len(t)):
        h = (h * BASE + t[k] - t[k - W] * powW) % MOD
        yield k - W + 1, h


def seeds(files, W, scope):
    """Groups [(fi, pos)] of windows of W tokens seen in >= 2 places.
    Two passes (a set of seen hashes, then positions of repeated hashes only) keep
    memory flat on the 1.3M-line vanilla tree.  scope == "cross": only windows
    present in a fork file are kept."""
    powW = pow(BASE, W, MOD)
    seen, rep = set(), set()
    for f in files:
        if scope == "cross" and not f.fork:
            continue
        for _, h in _hashes(f.toks, W, powW):
            (rep if h in seen else seen).add(h)
    if scope == "cross":
        rep = seen                                   # any fork window may match vanilla
    buckets = defaultdict(list)
    for fi, f in enumerate(files):
        for pos, h in _hashes(f.toks, W, powW):
            if h in rep:
                buckets[h].append((fi, pos))
    return [v for v in buckets.values() if 2 <= len(v) <= MAX_BUCKET
            and (scope != "cross" or len({files[fi].fork for fi, _ in v}) == 2)]


def maximal_pairs(files, W, groups, scope):
    pairs = []
    for g in groups:
        for a in range(len(g)):
            for b in range(a + 1, len(g)):
                (fi, i), (fj, j) = g[a], g[b]
                A, B = files[fi], files[fj]
                if scope == "cross" and A.fork == B.fork:
                    continue
                ta, tb = A.toks, B.toks
                if i > 0 and j > 0 and ta[i - 1] == tb[j - 1]:
                    continue                      # interior of a longer match
                if ta[i:i + W] != tb[j:j + W]:
                    continue                      # hash collision
                L = W
                lim = min(len(ta) - i, len(tb) - j)
                if fi == fj:
                    lo, hi = min(i, j), max(i, j)
                    if hi - lo < W:
                        continue                  # self-overlapping run (repetitive code)
                    lim = min(lim, hi - lo)
                while L < lim and ta[i + L] == tb[j + L]:
                    L += 1
                pairs.append((fi, i, fj, j, L))
    return pairs


def summarize(files, pairs, scope):
    dup = defaultdict(set)          # fi -> duplicated line indices
    per_pair = defaultdict(lambda: {"tokens": 0, "lines_a": set(), "lines_b": set(), "funcs": set()})
    for fi, i, fj, j, L in pairs:
        A, B = files[fi], files[fj]
        la = set(A.lines[i:i + L]); lb = set(B.lines[j:j + L])
        if scope != "cross" or A.fork:
            dup[fi] |= la
        if scope != "cross" or B.fork:
            dup[fj] |= lb
        key = (fi, fj) if fi <= fj else (fj, fi)
        rec = per_pair[key]
        rec["tokens"] += L
        (rec["lines_a"] if key[0] == fi else rec["lines_b"]).update(la)
        (rec["lines_b"] if key[0] == fi else rec["lines_a"]).update(lb)
        for li in la:
            if A.fn_of_line[li]:
                rec["funcs"].add(A.fn_of_line[li].split("::")[-1])
        for li in lb:
            if B.fn_of_line[li]:
                rec["funcs"].add(B.fn_of_line[li].split("::")[-1])
    return dup, per_pair


def families(files, per_pair, min_lines):
    """Connected components of files linked by >= min_lines shared lines."""
    parent = list(range(len(files)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    edges = []
    for (a, b), rec in per_pair.items():
        n = min(len(rec["lines_a"]), len(rec["lines_b"]))
        if a != b and n >= min_lines:
            parent[find(a)] = find(b)
            edges.append((a, b, n, rec))
    comp = defaultdict(list)
    for a, b, n, rec in edges:
        comp[find(a)].append((a, b, n, rec))
    fams = []
    for root, es in comp.items():
        members = sorted({x for a, b, *_ in es for x in (a, b)}, key=lambda k: files[k].rel)
        funcs = defaultdict(int)
        for a, b, n, rec in es:
            for fn in rec["funcs"]:
                funcs[fn] += 1
        fams.append({
            "files": [files[k].rel for k in members],
            "copies": len(members),
            "shared_lines_max_pair": max(n for *_, n, _ in es),
            "shared_lines_sum_pairs": sum(n for *_, n, _ in es),
            "pairs": [(files[a].rel, files[b].rel, n) for a, b, n, _ in sorted(es, key=lambda e: -e[2])],
            "functions": sorted(funcs, key=lambda k: -funcs[k]),
        })
    fams.sort(key=lambda d: -d["shared_lines_sum_pairs"])
    return fams


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=("exact", "norm"), default="exact")
    ap.add_argument("--scope", choices=("fork", "cross", "vanilla"), default="fork")
    ap.add_argument("--min-tokens", type=int, default=100)
    ap.add_argument("--family-lines", type=int, default=30)
    ap.add_argument("--vanilla-dirs", nargs="*", default=None,
                    help="restrict vanilla files to these SRC-relative dirs (e.g. element material)")
    ap.add_argument("--no-intra", action="store_true", help="ignore clones within one file")
    ap.add_argument("--add-fork", help="file listing extra repo-relative paths to treat as fork "
                                       "(e.g. fork-added files that lack the stamp)")
    ap.add_argument("--json")
    a = ap.parse_args()

    fork_p, van_p = split_fork_vanilla()
    if a.add_fork:
        extra = {ln.strip() for ln in Path(a.add_fork).read_text().splitlines() if ln.strip()}
        fork_p += [p for p in van_p if rel(p) in extra]
        van_p = [p for p in van_p if rel(p) not in extra]
    if a.vanilla_dirs:
        pre = tuple("SRC/" + d.strip("/") + "/" for d in a.vanilla_dirs)
        van_p = [p for p in van_p if rel(p).startswith(pre)]
    vocab = {}
    if a.scope == "fork":
        files = load(fork_p, a.mode, True, vocab)
    elif a.scope == "vanilla":
        files = load(van_p, a.mode, False, vocab)
    else:
        files = load(fork_p, a.mode, True, vocab) + load(van_p, a.mode, False, vocab)
    groups = seeds(files, a.min_tokens, a.scope)
    pairs = maximal_pairs(files, a.min_tokens, groups, a.scope)
    if a.no_intra:
        pairs = [q for q in pairs if q[0] != q[2]]
    dup, per_pair = summarize(files, pairs, a.scope)

    counted = [k for k, f in enumerate(files) if a.scope != "cross" or f.fork]
    total = sum(files[k].code for k in counted)
    duplines = sum(len(dup[k]) for k in counted)
    per_file = sorted(((files[k].rel, len(dup[k]), files[k].code) for k in counted if dup[k]),
                      key=lambda r: -r[1])
    fams = families(files, per_pair, a.family_lines)
    intra = sorted(((files[x].rel, len(rec["lines_a"] | rec["lines_b"]), sorted(rec["funcs"]))
                    for (x, y), rec in per_pair.items() if x == y), key=lambda r: -r[1])
    if a.scope == "cross":
        # per fork file: the vanilla file it shares most with
        best = {}
        for (x, y), rec in per_pair.items():
            fk, vk = (x, y) if files[x].fork else (y, x)
            n = len(rec["lines_a"] if fk == x else rec["lines_b"])
            if n > best.get(fk, (None, 0))[1]:
                best[fk] = (files[vk].rel, n)
        cross = sorted(((files[k].rel, v, n, files[k].code) for k, (v, n) in best.items()),
                       key=lambda r: -r[2])
    else:
        cross = None
    res = {
        "mode": a.mode, "scope": a.scope, "min_tokens": a.min_tokens,
        "files": len(counted), "code_lines": total, "dup_lines": duplines,
        "dup_pct": round(100.0 * duplines / max(total, 1), 2),
        "clone_pairs": len(pairs), "per_file": per_file, "families": fams, "cross": cross, "intra": intra,
    }
    print(f"[{a.mode}/{a.scope}/W={a.min_tokens}] files={res['files']} code_lines={total} "
          f"dup_lines={duplines} ({res['dup_pct']}%) clone_pairs={len(pairs)} families={len(fams)}")
    for fm in fams[:15]:
        print(f"  {fm['copies']} files, max-pair {fm['shared_lines_max_pair']} lines: "
              f"{', '.join(Path(x).name for x in fm['files'][:8])}"
              f"{' ...' if fm['copies'] > 8 else ''}")
    if cross:
        for r in cross[:15]:
            print(f"  {r[0]}  <-  {r[1]}  {r[2]}/{r[3]} lines")
    if a.json:
        Path(a.json).write_text(json.dumps(res, indent=1), encoding="utf-8")


if __name__ == "__main__":
    main()
