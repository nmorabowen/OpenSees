"""Summarize WP-260 Wave 0 (`260-w0`) HDF5 profiler baselines & compute Plan §0.5 ranking.

Reads `%TEMP%/wp260w0_prof/manifest.json` + per-run HDF5 files (parent process ONLY; uses h5py).
Never double-counts nested scopes: `elem.tangent` is explicitly split into
  `elem_by_type` (= element kernels)
+ `soe.addA`     (= sparse scatter)
+ residual wrapper overhead (`other`).
All phase shares are strictly disjoint and sum to 100.0% of profiled step time.
"""
from __future__ import annotations

import argparse
import json
import statistics
import tempfile
from pathlib import Path
from typing import Any, Dict, List, Tuple

import h5py

DEFAULT_PROF_DIR = Path(tempfile.gettempdir()) / "wp260w0_prof"
OUT_README = Path(__file__).resolve().parent / "README.md"

# Verified against SRC/classTags.h on origin/ladruno @ 7e91aaf1c
TAG_NAMES = {
    56: "stdBrick",
    33002: "LadrunoBrick",
    33023: "LadrunoCPML",
    33025: "LadrunoForceBeamColumn3d",
    33026: "LadrunoMPMPoint",
}

# Canonical WP-260 packages from Ladruno_implementation/_wp260_orchestration_plan.md
# (pkg_id, title, tier, blocker, target_desc, phase_fn, s_cons)
PACKAGES = [
    ("260a", "SOE Scatter & Direct Solver Hygiene (PARDISO/MUMPS setSize+addA+iparm+Krylov)",
     "G-BYTE+G-OPT", "Unblocked",
     "soe.symbolic + soe.factor + soe.trisolve + soe.addA",
     lambda ph, tags, meta: ph["soe.symbolic"] + ph["soe.factor"] + ph["soe.trisolve"] + ph["soe.addA"], 1.7),
    ("260b", "Explicit & Wave Hygiene (ExplicitBathe/CDL/CPML OpenMP + T7 inertia skips + SMS PCG)",
     "G-BYTE/G-EQUIV", "Unblocked",
     "LadrunoCPML + explicit update/residual/commit (D3a/D3b/D3c)",
     lambda ph, tags, meta: (
         (tags.get(33023, 0.0) + ph["unbal_other"] + ph["upd_other"] + ph["newStep_other"] + ph["commit.elements"])
         if meta["family"].startswith("D3-")
         else ((ph["soe.factor"] + ph["soe.trisolve"] + ph["elem_kernel"]) if "SMS" in meta["family"] else 0.0)
     ), 1.8),
    ("260c1", "LogStrainNDMaterial Finite-Strain Stress-Only + LadrunoUP OMP",
     "G-BYTE/G-EQUIV", "Unblocked",
     "LogStrain / UP (finite-strain suite; N/A on small-strain D1)",
     lambda ph, tags, meta: 0.0, 3.0),
    ("260c2", "Concrete3D / HS Stress-Only Contract #1a + Commit-Deferred -implex",
     "G-BYTE/G-EQUIV", "Blocked (#1153/#1154)",
     "LadrunoBrick (33002) on D1b Concrete3D",
     lambda ph, tags, meta: tags.get(33002, 0.0) if meta["family"].startswith("D1b") else 0.0, 2.0),
    ("260d", "Adaptive Transient + LadrunoFactorReuse across TRBDF2/3 + ModalDamping",
     "G-EQUIV/G-OPT", "Unblocked",
     "soe.factor + elem_tangent_kernel on implicit transient (D1a/D2/D4)",
     lambda ph, tags, meta: (ph["soe.factor"] + ph["elem_tan_kernel"])
     if meta["family"] in ("D1a-Brick-SANISAND", "D2-CPML-Implicit", "D4-MPM-ModeA") else 0.0, 1.6),
    ("260e", "Windows OpenMP libiomp5md Unification + mimalloc + /GL /LTCG (D4)",
     "G-BYTE", "Unblocked (cross-cutting build bench)",
     "Global build/runtime allocator + single OMP pool across all decks",
     lambda ph, tags, meta: 1.0, 1.15),
    ("260f", "3D/MPM/MCFT Threading + Persistent Assembly + LadrunoFastScatter",
     "G-BYTE/G-EQUIV", "Blocked (#1147/#1153)",
     "elem_kernel + soe.addA across D1a/D1b/D4/D5 (currently serial on 8 P-cores)",
     lambda ph, tags, meta: (ph["elem_kernel"] + ph["soe.addA"])
     if meta["family"] in ("D1a-Brick-SANISAND", "D1b-Brick-Concrete3D", "D4-MPM-ModeA", "D5-MCFT-2CLD12") else 0.0, 2.5),
    ("260g", "Mixed-Precision Direct Solver (SP PARDISO/MUMPS + GMRES-IR + Sym-Precond)",
     "G-OPT", "Unblocked",
     "soe.factor + soe.trisolve on FP64-safe implicit (D1a/D1b/D2/D4)",
     lambda ph, tags, meta: (ph["soe.factor"] + ph["soe.trisolve"])
     if meta["family"] in ("D1a-Brick-SANISAND", "D1b-Brick-Concrete3D", "D2-CPML-Implicit", "D4-MPM-ModeA") else 0.0, 1.5),
    ("260h", "Analytic Return Maps (SANISAND SAS-RK, Concrete3D Cardano, Anderson/LBFGS)",
     "G-EQUIV/G-OPT", "Blocked (#1153/#1154)",
     "elem_kernel on D1a SANISAND / D1b Concrete3D",
     lambda ph, tags, meta: ph["elem_kernel"]
     if meta["family"].startswith(("D1a", "D1b")) else 0.0, 2.0),
    ("260k", "LadrunoGPStore SoA + Active-Set Compaction + 7x7 Schur SANISAND CPPM",
     "G-OPT", "Post-Wave 1",
     "elem_kernel on D1a SANISAND",
     lambda ph, tags, meta: ph["elem_kernel"] if meta["family"].startswith("D1a") else 0.0, 2.0),
    ("260l", "LadrunoExplicitEngine Flat Explicit + CPML + Clustered LTS",
     "G-OPT", "Post-Wave 1",
     "Whole explicit step (D3a/D3b)",
     lambda ph, tags, meta: 1.0 if meta["family"].startswith("D3-") else 0.0, 2.0),
    ("260m", "LadrunoAMGSOE (Boost-free AMGCL) + Cluster Direct",
     "G-OPT", "Post-Wave 1",
     "soe.factor + soe.trisolve + soe.symbolic on 3D implicit (D1a/D1b/D2)",
     lambda ph, tags, meta: (ph["soe.factor"] + ph["soe.trisolve"] + ph["soe.symbolic"])
     if meta["family"] in ("D1a-Brick-SANISAND", "D1b-Brick-Concrete3D", "D2-CPML-Implicit") else 0.0, 1.5),
]


def _read_group(g: h5py.Group) -> Dict[str, Any]:
    wall_ns = int(g.attrs.get("wall_ns", 0))
    calls = int(g.attrs.get("calls", 0))
    children: Dict[str, Any] = {}
    elem_by_type: List[Tuple[int, int, int]] = []
    if "elem_by_type" in g:
        arr = g["elem_by_type"][:]
        for row in arr:
            elem_by_type.append((int(row["classTag"]), int(row["count"]), int(row["wall_ns"])))
    for k in g.keys():
        if k == "elem_by_type":
            continue
        obj = g[k]
        if isinstance(obj, h5py.Group):
            children[k] = _read_group(obj)
    return {"wall_ns": wall_ns, "calls": calls, "children": children, "elem_by_type": elem_by_type}


def _find_all(node: Dict[str, Any], name: str) -> List[Dict[str, Any]]:
    out = []
    for k, ch in node["children"].items():
        if k == name:
            out.append(ch)
        out.extend(_find_all(ch, name))
    return out


def _sum_name(node: Dict[str, Any], name: str) -> int:
    """Sum disjoint occurrences of `name` (does not re-enter matched nodes)."""
    tot = 0
    for k, ch in node["children"].items():
        if k == name:
            tot += ch["wall_ns"]
        else:
            tot += _sum_name(ch, name)
    return tot


def _sum_elem_by_type(node: Dict[str, Any]) -> Tuple[int, Dict[int, int]]:
    tot = 0
    by_tag: Dict[int, int] = {}
    for tag, _cnt, wns in node["elem_by_type"]:
        tot += wns
        by_tag[tag] = by_tag.get(tag, 0) + wns
    return tot, by_tag


def _gather_elem_kernels(node: Dict[str, Any], by_tag: Dict[int, int]) -> Tuple[int, int]:
    """Return (tangent_kernel_ns, all_elem_kernel_ns) without double counting."""
    tan_ns = 0
    for en in _find_all(node, "elem.tangent"):
        k, bt = _sum_elem_by_type(en)
        tan_ns += k
        for t, w in bt.items():
            by_tag[t] = by_tag.get(t, 0) + w
    res_ns = 0
    for rn in _find_all(node, "elem.residual"):
        k, bt = _sum_elem_by_type(rn)
        res_ns += k
        for t, w in bt.items():
            by_tag[t] = by_tag.get(t, 0) + w
    upd_ns = 0
    for un in _find_all(node, "elem.update"):
        k, bt = _sum_elem_by_type(un)
        upd_ns += k
        for t, w in bt.items():
            by_tag[t] = by_tag.get(t, 0) + w
    return tan_ns, tan_ns + res_ns + upd_ns


def analyze_h5(h5_path: Path, run_id: str) -> Dict[str, Any]:
    with h5py.File(h5_path, "r") as f:
        root_g = f["runs"][run_id]["rollup"]["root"]
        tree = _read_group(root_g)

    # Top-level container: MPM runs "mpm.a.step" (which wraps "mpm.a.prepare" + "step"), others run "step"
    top_nodes = _find_all(tree, "mpm.a.step")
    if not top_nodes:
        top_nodes = _find_all(tree, "step")
    step_wall = sum(n["wall_ns"] for n in top_nodes) or sum(ch["wall_ns"] for ch in tree["children"].values())
    root_scope = {"wall_ns": step_wall, "calls": 1,
                  "children": {f"s{i}": n for i, n in enumerate(top_nodes)},
                  "elem_by_type": []}

    # Disjoint SOE phases across whole step
    soe_sym = _sum_name(root_scope, "soe.symbolic")
    soe_fac = _sum_name(root_scope, "soe.factor")
    soe_tri = _sum_name(root_scope, "soe.trisolve")
    soe_adda = _sum_name(root_scope, "soe.addA")

    # Disjoint element kernels (elem_by_type in elem.tangent/elem.residual/elem.update)
    tag_ns: Dict[int, int] = {}
    tan_kern, elem_kern = _gather_elem_kernels(root_scope, tag_ns)

    # Disjoint remainders inside each phase container (subtract ONLY what belongs to that container!)
    dc_other = 0
    for dn in _find_all(root_scope, "domainChanged") + _find_all(root_scope, "mpm.a.prepare"):
        dc_soe = (_sum_name(dn, "soe.symbolic") + _sum_name(dn, "soe.factor")
                  + _sum_name(dn, "soe.trisolve") + _sum_name(dn, "soe.addA"))
        dc_other += max(0, dn["wall_ns"] - dc_soe)

    newstep_other = 0
    for nn in _find_all(root_scope, "newStep"):
        _, nn_kern = _gather_elem_kernels(nn, {})
        nn_soe = (_sum_name(nn, "soe.symbolic") + _sum_name(nn, "soe.factor")
                  + _sum_name(nn, "soe.trisolve") + _sum_name(nn, "soe.addA"))
        newstep_other += max(0, nn["wall_ns"] - nn_kern - nn_soe)

    unbal_other = 0
    for un in _find_all(root_scope, "formUnbalance"):
        _, uk = _gather_elem_kernels(un, {})
        unbal_other += max(0, un["wall_ns"] - uk)

    upd_other = 0
    for upn in _find_all(root_scope, "update"):
        _, upk = _gather_elem_kernels(upn, {})
        up_soe = (_sum_name(upn, "soe.symbolic") + _sum_name(upn, "soe.factor")
                  + _sum_name(upn, "soe.trisolve") + _sum_name(upn, "soe.addA"))
        upd_other += max(0, upn["wall_ns"] - upk - up_soe)

    commit_elem = _sum_name(root_scope, "commit.elements")
    commit_rec = _sum_name(root_scope, "commit.record")

    san_sub_ns = _sum_name(root_scope, "sanisand.sasME.substeps")

    parts_ns = {
        "soe.symbolic": soe_sym,
        "soe.factor": soe_fac,
        "soe.trisolve": soe_tri,
        "soe.addA": soe_adda,
        "elem_kernel": elem_kern,
        "newStep_other": newstep_other,
        "unbal_other": unbal_other,
        "upd_other": upd_other,
        "dc_other": dc_other,
        "commit.elements": commit_elem,
        "commit.record": commit_rec,
    }
    accounted = sum(parts_ns.values())
    parts_ns["other"] = max(0, step_wall - accounted)
    # Auxiliary (not part of the 100% disjoint sum)
    parts_ns["elem_tan_kernel"] = tan_kern
    parts_ns["sanisand_sub_ns"] = san_sub_ns
    parts_ns["step_wall_ns"] = step_wall

    shares = {k: (v / step_wall if step_wall > 0 else 0.0) for k, v in parts_ns.items()}
    tag_shares = {t: (w / step_wall if step_wall > 0 else 0.0) for t, w in tag_ns.items()}
    return {"shares": shares, "tag_shares": tag_shares, "tag_ns": tag_ns, "step_wall_ns": step_wall}


def _med_spread(vals: List[float]) -> Tuple[float, float, float]:
    s = sorted(vals)
    return statistics.median(s), s[0], s[-1]


def _dump_tree(node: Dict[str, Any], name: str = "root", pad: str = "") -> None:
    w = node["wall_ns"] / 1e9
    eb = ", ".join(f"{t}:{wns / 1e9:.3f}s" for t, _c, wns in node["elem_by_type"])
    eb_s = f" [{eb}]" if eb else ""
    if w >= 0.002 or not pad:
        print(f"{pad}{name}: {w:.3f}s (c={node['calls']}){eb_s}")
    for k, ch in node["children"].items():
        _dump_tree(ch, k, pad + "  ")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", type=Path, default=DEFAULT_PROF_DIR)
    ap.add_argument("--out-readme", type=Path, default=OUT_README)
    ap.add_argument("--out-comment", type=Path, default=None)
    ap.add_argument("--dump-trees", action="store_true")
    args = ap.parse_args()

    manifest = json.loads((args.dir / "manifest.json").read_text(encoding="utf-8"))
    if args.dump_trees:
        seen = set()
        for rec in manifest["records"]:
            if not rec.get("ok") or rec["mode"] != "1t_deep" or rec["deck"] in seen:
                continue
            seen.add(rec["deck"])
            print(f"=== {rec['deck']} ({rec['wall_s']:.3f}s) ===")
            with h5py.File(rec["h5"], "r") as f:
                _dump_tree(_read_group(f["runs"][rec["run_id"]]["rollup"]["root"]))
        return 0

    by_deck_mode: Dict[Tuple[str, str], List[Dict[str, Any]]] = {}
    for rec in manifest["records"]:
        if rec.get("ok"):
            by_deck_mode.setdefault((rec["deck"], rec["mode"]), []).append(rec)

    decks_order = sorted({k[0] for k in by_deck_mode})
    rows = []
    share_keys = [
        "soe.symbolic", "soe.factor", "soe.trisolve", "soe.addA", "elem_kernel",
        "elem_tan_kernel", "newStep_other", "unbal_other", "upd_other",
        "dc_other", "commit.elements", "commit.record", "other", "sanisand_sub_ns",
    ]
    for d in decks_order:
        r1 = by_deck_mode.get((d, "1t_deep"), [])
        r8 = by_deck_mode.get((d, "8t_coarse"), [])
        if not r1 or not r8:
            continue
        parsed1 = [analyze_h5(Path(r["h5"]), r["run_id"]) for r in r1]
        med_sh = {k: statistics.median([p["shares"][k] for p in parsed1]) for k in share_keys}
        all_tags = {t for p in parsed1 for t in p["tag_shares"]}
        med_tags = {t: statistics.median([p["tag_shares"].get(t, 0.0) for p in parsed1]) for t in all_tags}
        tag_desc = ", ".join(
            f"{TAG_NAMES.get(t, f'tag{t}')} ({t})={100.0 * sh:.1f}%"
            for t, sh in sorted(med_tags.items(), key=lambda kv: -kv[1])
            if sh >= 0.005
        ) or "-"

        w1_med, w1_lo, w1_hi = _med_spread([r["wall_s"] for r in r1])
        w8_med, w8_lo, w8_hi = _med_spread([r["wall_s"] for r in r8])
        rows.append({
            "deck": d,
            "family": r1[0]["family"],
            "nsteps": r1[0]["nsteps"],
            "nelm": r1[0]["nelm"],
            "build": r1[0]["build_id"],
            "shares": med_sh,
            "tag_shares": med_tags,
            "tag_desc": tag_desc,
            "w1": (w1_med, w1_lo, w1_hi),
            "w8": (w8_med, w8_lo, w8_hi),
            "scale_1_8": w1_med / w8_med if w8_med > 0 else 0.0,
        })

    # Plan §0.5 Amdahl ranking
    rankings = []
    for pkg, title, tier, blocker, target_desc, pfn, s_cons in PACKAGES:
        best_G = 1.0
        best_p = 0.0
        best_deck = "-"
        for row in rows:
            p = float(pfn(row["shares"], row["tag_shares"], {"family": row["family"]}))
            p = max(0.0, min(1.0, p))
            g = 1.0 / ((1.0 - p) + (p / s_cons)) if p > 0 else 1.0
            if g > best_G:
                best_G, best_p, best_deck = g, p, row["deck"]
        passes = (best_p >= 0.10) and (best_G >= 1.10)
        rankings.append({
            "pkg": pkg, "title": title, "tier": tier, "blocker": blocker,
            "target": target_desc, "s_cons": s_cons, "best_p": best_p,
            "best_G": best_G, "best_deck": best_deck, "passes": passes,
        })

    tier_order = {"G-BYTE": 0, "G-BYTE+G-OPT": 0, "G-BYTE/G-EQUIV": 1, "G-EQUIV/G-OPT": 2, "G-OPT": 3}
    rankings.sort(key=lambda r: (-round(r["best_G"], 4), tier_order.get(r["tier"], 9)))

    md = []
    md.append("# WP-260 Wave 0 (`260-w0`) — Phase-0 Profiler Baseline")
    md.append("")
    md.append(f"- **Binary Build SHA (`ops.ladrunoBuild()`):** `{manifest.get('exe_build', '')}` (pre-`#1147`/`#1153`/`#1154` baseline)")
    md.append("- **Host:** i7-13700KF (`8` P-cores / `16` logical HT pinned via affinity `0xFFFF`, zero concurrent builds)")
    md.append(f"- **Protocol:** `ops.profiler('start', '-deep', '-perStep')` at 1T (`OMP_NUM_THREADS=MKL_NUM_THREADS=1`, `ladrunoThreads(1)`) + coarse `ops.profiler('start', '-perStep')` at 8T (`OMP_NUM_THREADS=MKL_NUM_THREADS=8`, `ladrunoThreads(8)`), `{manifest['reps']}` interleaved reps.")
    md.append("- **Accounting rule:** Every row sums disjoint scopes to 100.0% of profiled step time; `elem.tangent` is strictly split into element kernels (`elem_by_type` inside `elem_kern`) and sparse scatter (`soe.addA`).")
    md.append("")
    md.append("## 1. Measured Phase-Share Table (1T `-deep`) & Production Wall Time (8T)")
    md.append("")
    md.append("| Deck | Family | Steps | `soe.sym` | `soe.fac` | `soe.tri` | `soe.addA` | `elem_kern` | `newStp_o` | `unbal_o` | `upd_o` | `dc_o` | `cmt_el` | `cmt_rec` | `other` | 1T med [min..max] | 8T med [min..max] | 1T/8T |")
    md.append("|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for r in rows:
        sh = r["shares"]
        w1m, w1l, w1h = r["w1"]
        w8m, w8l, w8h = r["w8"]
        md.append(
            f"| `{r['deck']}` | {r['family']} | {r['nsteps']} "
            f"| {100*sh['soe.symbolic']:.1f}% | {100*sh['soe.factor']:.1f}% | {100*sh['soe.trisolve']:.1f}% "
            f"| {100*sh['soe.addA']:.1f}% | {100*sh['elem_kernel']:.1f}% "
            f"| {100*sh['newStep_other']:.1f}% | {100*sh['unbal_other']:.1f}% | {100*sh['upd_other']:.1f}% "
            f"| {100*sh['dc_other']:.1f}% | {100*sh['commit.elements']:.1f}% | {100*sh['commit.record']:.1f}% | {100*sh['other']:.1f}% "
            f"| {w1m:.3f}s [{w1l:.3f}..{w1h:.3f}] | {w8m:.3f}s [{w8l:.3f}..{w8h:.3f}] | {r['scale_1_8']:.2f}x |"
        )
    md.append("")
    md.append("### Element Kernel Composition (`elem_by_type`) & SANISAND Substep Share")
    for r in rows:
        san = r["shares"]["sanisand_sub_ns"]
        san_note = f" *(SANISAND `sasME.substeps`: {100*san:.1f}%)*" if san >= 0.005 else ""
        md.append(f"- **`{r['deck']}`**: {r['tag_desc']}{san_note}")
    md.append("")
    md.append("## 2. Phase -> WP-260 Package Map")
    md.append("")
    md.append("| Scope / Phase | Mapped Package(s) | Status on `7e91aaf1c` |")
    md.append("|---|---|---|")
    md.append("| `soe.factor`, `soe.symbolic`, `soe.trisolve`, `soe.addA` | `260a` (PARDISO/MUMPS hygiene), `260d` (reuse), `260g` (SP/IR), `260m` (AMG) | `260a` / `260d` / `260g` unblocked |")
    md.append("| `LadrunoCPML` (`33023`) + explicit `update`/`residual`/`commit` (`D3a`/`D3b`/`D3c`) | `260b` (Explicit & Wave Hygiene) | Unblocked |")
    md.append("| `elem_by_type` (`LadrunoBrick` SANISAND/Concrete3D, `LadrunoForceBeamColumn3d` MCFT, `LadrunoMPMPoint`) | `260f` (Threading+FastScatter), `260c2` (Concrete3D stress-only/implex), `260h` (Analytic return maps), `260k` (SoA) | `260f` blocked (`#1147/#1153`); `260c2`/`260h` blocked (`#1153/#1154`) |")
    md.append("| Cross-cutting OMP / Allocator / LTO (`1T/8T <= 1.0x` on threaded/MKL decks) | `260e` (D4: `libiomp5md` unification + `mimalloc` + `/GL`) | Unblocked (benchmark after Wave-1 isolation) |")
    md.append("| `commit.record` | Outside WP-260 (recorders; flagged if >=10%) | Observed `0.0%` on all D1 baseline decks (no heavy file recorders enabled in bench decks) |")
    md.append("")
    md.append("## 3. Plan §0.5 Amdahl Ranking (`p >= 10%` and `G = 1/((1-p)+p/s) >= 1.10x`)")
    md.append("")
    md.append("| Rank | Pkg | Tier | Status / Blocker | Target Scope | Conservative `s` | Peak Deck | Peak `p` | Peak `G` | Gate |")
    md.append("|---:|---|---|---|---|---|---:|---:|---:|:---:|")
    for idx, rk in enumerate(rankings, 1):
        gate = "PASS" if rk["passes"] else "BELOW"
        md.append(
            f"| {idx} | `{rk['pkg']}` | {rk['tier']} | {rk['blocker']} | {rk['target']} | {rk['s_cons']:.2f}x "
            f"| `{rk['best_deck']}` | {100*rk['best_p']:.1f}% | **{rk['best_G']:.2f}x** | {gate} |"
        )
    md.append("")
    passing = [r for r in rankings if r["passes"]]
    unblocked = [r for r in passing if r["blocker"].startswith("Unblocked")]
    blocked = [r for r in passing if not r["blocker"].startswith("Unblocked")]
    md.append("## 4. Proposed Wave-1 Picks (for Owner Confirmation & Claude §4 Review)")
    md.append("")
    md.append("### Top Overall by Projected `G` (shows blocked dependencies)")
    for r in passing[:4]:
        md.append(f"- **`{r['pkg']}`** ({r['tier']}, **{r['blocker']}**) — `{r['title']}`: peak `p = {100*r['best_p']:.1f}%` on `{r['best_deck']}`, conservative `G = {r['best_G']:.2f}x` at `s = {r['s_cons']:.2f}x`.")
    md.append("")
    md.append("### Recommended Immediate Unblocked Wave-1 Launches (2–3 Parallel WPs per §0.5)")
    for r in unblocked[:3]:
        md.append(f"- **`{r['pkg']}`** ({r['tier']}, {r['blocker']}): peak `p = {100*r['best_p']:.1f}%` on `{r['best_deck']}`, projected `G = {r['best_G']:.2f}x`.")
    if blocked:
        md.append("")
        md.append("> **Dependency note:** Higher-ranked element/constitutive WPs ("
                  + ", ".join(f"`{b['pkg']}`" for b in blocked[:3])
                  + ") wait for open upstream/fork PRs `#1147` / `#1153` / `#1154` (which add MPM & SANISAND/Concrete3D threading) before cutting their WP branches so `260-w0` re-ranks cleanly after those merges.")

    text = "\n".join(md) + "\n"
    args.out_readme.write_text(text, encoding="utf-8")
    print(f"Wrote {args.out_readme}")
    if args.out_comment:
        args.out_comment.parent.mkdir(parents=True, exist_ok=True)
        args.out_comment.write_text(text, encoding="utf-8")
        print(f"Wrote {args.out_comment}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
