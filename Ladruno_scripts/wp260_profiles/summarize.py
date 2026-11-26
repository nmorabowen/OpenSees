"""Summarize WP-260 Wave 0 (`260-w0`) profiles & compute Plan §0.5 + brief §4 ranking.

Supports both raw HDF5 directories (`--dir`) and committed JSON snapshots (`--from-json`)
without needing `.h5` files or `opensees.pyd`. Always writes `RESULTS.md` (never overwrites `README.md`).

Accounting rule (`wp-briefs.md` §3):
  - `elem.tangent` = `elem_by_type` + `soe.addA` + unbucketed wrapper (`other`).
  - `mpm.a.prepare` gets its own column (never folded into `dc_other`).
  - `soe.cgs` gets its own column.
  - Prints `residual = step_wall - accounted`; NEVER clamps `other` at 0; fails if sum > 100.5%.
"""
from __future__ import annotations

import argparse
import json
import statistics
import tempfile
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

DEFAULT_PROF_DIR = Path(tempfile.gettempdir()) / "wp260w0_prof"
HERE = Path(__file__).resolve().parent
OUT_RESULTS = HERE / "RESULTS.md"

TAG_NAMES = {
    56: "stdBrick",
    33002: "LadrunoBrick",
    33023: "LadrunoCPML",
    33025: "LadrunoForceBeamColumn3d",
    33026: "LadrunoMPMPoint",
}


# --------------------------------------------------------------------------- #
#  Exact Brief §4 Package Table (Plan hypothesis, rule-2 s, source, p(deck))  #
# --------------------------------------------------------------------------- #
def _p_260a(sh: Dict[str, float], tags: Dict[int, float], meta: Dict[str, Any]) -> Optional[float]:
    # G-OPT PARDISO iparm levers target soe.factor on Pardiso decks; MUMPS items don't run on Pardiso D1 decks
    return 0.0 if meta["family"].startswith("D3-") else sh["soe.factor"]


def _p_260b(sh: Dict[str, float], tags: Dict[int, float], meta: Dict[str, Any]) -> Optional[float]:
    fam = meta["family"]
    if fam in ("D3-CentralDifferenceLadruno", "D3-ExplicitBathe"):
        p = tags.get(33023, 0.0) + sh["commit.elements"]
        if fam == "D3-ExplicitBathe":
            p += sh["formTangent_all"]
        return p
    if fam == "D3c-SMS-Consistent":
        return sh["soe.trisolve"]
    return 0.0


def _p_260c1(sh: Dict[str, float], tags: Dict[int, float], meta: Dict[str, Any]) -> Optional[float]:
    return 0.0  # No finite-strain LogStrain deck in D1 small-strain suite


def _p_260c2(sh: Dict[str, float], tags: Dict[int, float], meta: Dict[str, Any]) -> Optional[float]:
    if meta["family"] == "D1b-Brick-Concrete3D":
        # Concrete3D work in D1b: LadrunoBrick (33002) buckets + note unbucketed enforceSPs update when present
        return tags.get(33002, 0.0)
    return 0.0


def _p_260f_narrow(sh: Dict[str, float], tags: Dict[int, float], meta: Dict[str, Any]) -> Optional[float]:
    # Threading #1147/#1153/#1154 deliver is theirs; 260f owns persistent assembly + soe.addA scatter only
    return sh["soe.addA"]


def _p_260h_narrow(sh: Dict[str, float], tags: Dict[int, float], meta: Dict[str, Any]) -> Optional[float]:
    if meta["family"] == "D1a-Brick-SANISAND":
        return sh["sanisand_sub_ns"]
    if meta["family"] == "D1b-Brick-Concrete3D":
        return tags.get(33002, 0.0)
    return 0.0


PACKAGE_SPECS = [
    # (pkg, tier, rank_status, blocker, plan_hyp, s_rule2, source, p_desc, p_fn)
    ("260a", "G-BYTE+G-OPT", "eligible", "Unblocked",
     "1.3-2x PARDISO factor", 1.15, "plan §1 rule-2: 1+(1.3-1)/2 (MUMPS 3-20x N/A on Pardiso)",
     "soe.factor (Pardiso)", _p_260a),
    ("260b", "G-BYTE/G-EQUIV", "eligible", "Unblocked",
     "1.4-2.2x explicit", 1.20, "plan §1 rule-2: 1+(1.4-1)/2",
     "CPML + cmt_el (+Bathe tan / SMS tri)", _p_260b),
    ("260c1", "G-BYTE/G-EQUIV", "eligible", "Unblocked",
     "3-6x finite-strain GP", 2.00, "plan §1 rule-2: 1+(3.0-1)/2",
     "LogStrain GP (none in D1)", _p_260c1),
    ("260c2", "G-BYTE/G-EQUIV", "eligible", "Blocked (#1153/#1154)",
     "1.5-4x Concrete3D/HS", 1.25, "plan §1 rule-2: 1+(1.5-1)/2",
     "LadrunoBrick Concrete3D (D1b)", _p_260c2),
    ("260d", "G-EQUIV/G-OPT", "not phase-rankable", "Unblocked",
     "2-5x adaptive runs", None, "changes step count -> whole-run bench required",
     "not phase-rankable", None),
    ("260e", "G-BYTE", "not phase-rankable", "D4 + build bench",
     "to be measured", None, "global allocator/OMP -> whole-run bench required",
     "not phase-rankable", None),
    ("260f", "G-BYTE/G-EQUIV", "unranked", "Blocked (#1147/#1153)",
     "no number in plan", None, "unranked (threading belongs to #1147/#1153/#1154)",
     "soe.addA only", _p_260f_narrow),
    ("260g", "G-OPT", "unranked", "Unblocked",
     "no number in plan", None, "unranked (mixed-precision PARDISO requires bench)",
     "soe.factor + soe.trisolve", lambda sh, tags, m: sh["soe.factor"] + sh["soe.trisolve"]),
    ("260h", "G-EQUIV/G-OPT", "unranked", "Blocked (#1153/#1154)",
     "no number in plan", None, "unranked (SAS-RK / analytic Concrete3D micro-bench needed)",
     "sasME.substeps / Concrete3D", _p_260h_narrow),
    ("260i", "G-OPT", "unranked", "Post-Wave 1",
     "no number in plan", None, "unranked (global IMPL-EX / H5DRM v2)",
     "-", None),
    ("260j", "G-OPT", "unranked", "Post-Wave 1",
     "no number in plan", None, "unranked (contact BVH)",
     "-", None),
    ("260k", "G-OPT", "not eligible (Wave 3)", "Rule 6", "Wave 3", None, "not eligible for Wave 1 pick", "-", None),
    ("260l", "G-OPT", "not eligible (Wave 3)", "Rule 6", "Wave 3", None, "not eligible for Wave 1 pick", "-", None),
    ("260m", "G-OPT", "not eligible (Wave 3)", "Rule 6", "Wave 3", None, "not eligible for Wave 1 pick", "-", None),
    ("260n", "G-OPT", "not eligible (Wave 3)", "Rule 6", "Wave 3", None, "not eligible for Wave 1 pick", "-", None),
    ("260o", "G-OPT", "not eligible (Wave 3)", "Rule 6", "Wave 3", None, "not eligible for Wave 1 pick", "-", None),
]


# --------------------------------------------------------------------------- #
#  Tree helpers (compatible with HDF5 loader AND trees.json)                  #
# --------------------------------------------------------------------------- #
def _load_h5_tree(h5_path: Path, run_id: str) -> Dict[str, Any]:
    import h5py

    def rd(g) -> Dict[str, Any]:
        node = {
            "wall_ns": int(g.attrs.get("wall_ns", 0)),
            "calls": int(g.attrs.get("calls", 0)),
            "elem_by_type": [
                [int(r["classTag"]), int(r["count"]), int(r["wall_ns"])] for r in g["elem_by_type"][:]
            ] if "elem_by_type" in g else [],
            "children": {},
        }
        for k in g.keys():
            if k != "elem_by_type" and isinstance(g[k], h5py.Group):
                node["children"][k] = rd(g[k])
        return node

    with h5py.File(h5_path, "r") as f:
        return rd(f["runs"][run_id]["rollup"]["root"])


def _find_all(node: Dict[str, Any], name: str) -> List[Dict[str, Any]]:
    out = []
    for k, ch in node["children"].items():
        if k == name:
            out.append(ch)
        out.extend(_find_all(ch, name))
    return out


def _sum_name(node: Dict[str, Any], name: str) -> int:
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
    for tag, _cnt, wns in node.get("elem_by_type", []):
        tot += wns
        by_tag[int(tag)] = by_tag.get(int(tag), 0) + int(wns)
    return tot, by_tag


def _gather_elem_kernels(node: Dict[str, Any], by_tag: Dict[int, int]) -> int:
    tot = 0
    for scope_name in ("elem.tangent", "elem.residual", "elem.update"):
        for sn in _find_all(node, scope_name):
            k, bt = _sum_elem_by_type(sn)
            tot += k
            for t, w in bt.items():
                by_tag[t] = by_tag.get(t, 0) + w
    return tot


def analyze_tree(tree: Dict[str, Any]) -> Dict[str, Any]:
    # MPM runs wrap "mpm.a.prepare" + "step" inside "mpm.a.step"
    top_nodes = _find_all(tree, "mpm.a.step")
    if not top_nodes:
        top_nodes = _find_all(tree, "step")
    step_wall = sum(n["wall_ns"] for n in top_nodes) or sum(ch["wall_ns"] for ch in tree["children"].values())
    root_scope = {"wall_ns": step_wall, "calls": 1,
                  "children": {f"s{i}": n for i, n in enumerate(top_nodes)},
                  "elem_by_type": []}

    soe_sym = _sum_name(root_scope, "soe.symbolic")
    soe_fac = _sum_name(root_scope, "soe.factor")
    soe_tri = _sum_name(root_scope, "soe.trisolve")
    soe_cgs = _sum_name(root_scope, "soe.cgs")
    soe_adda = _sum_name(root_scope, "soe.addA")

    mpm_prep = _sum_name(root_scope, "mpm.a.prepare")

    tag_ns: Dict[int, int] = {}
    elem_kern = _gather_elem_kernels(root_scope, tag_ns)

    def sub_soe(n: Dict[str, Any]) -> int:
        return (_sum_name(n, "soe.symbolic") + _sum_name(n, "soe.factor")
                + _sum_name(n, "soe.trisolve") + _sum_name(n, "soe.cgs") + _sum_name(n, "soe.addA"))

    dc_other = sum(n["wall_ns"] - sub_soe(n) for n in _find_all(root_scope, "domainChanged"))
    newstep_other = sum(n["wall_ns"] - _gather_elem_kernels(n, {}) - sub_soe(n)
                        for n in _find_all(root_scope, "newStep"))
    unbal_other = sum(n["wall_ns"] - _gather_elem_kernels(n, {})
                      for n in _find_all(root_scope, "formUnbalance"))
    upd_other = sum(n["wall_ns"] - _gather_elem_kernels(n, {}) - sub_soe(n)
                    for n in _find_all(root_scope, "update"))

    commit_elem = _sum_name(root_scope, "commit.elements")
    commit_rec = _sum_name(root_scope, "commit.record")

    # Never clamp! Compute exact residual = step_wall - sum(disjoint accounted parts)
    disjoint_ns = {
        "soe.symbolic": soe_sym,
        "soe.factor": soe_fac,
        "soe.trisolve": soe_tri,
        "soe.cgs": soe_cgs,
        "soe.addA": soe_adda,
        "elem_kernel": elem_kern,
        "mpm.a.prepare": mpm_prep,
        "newStep_other": newstep_other,
        "unbal_other": unbal_other,
        "upd_other": upd_other,
        "dc_other": dc_other,
        "commit.elements": commit_elem,
        "commit.record": commit_rec,
    }
    accounted_ns = sum(disjoint_ns.values())
    residual_ns = step_wall - accounted_ns
    overcount_ratio = (accounted_ns / step_wall) if step_wall > 0 else 0.0
    if overcount_ratio > 1.005:
        raise RuntimeError(
            f"Scope accounting overcount detected: accounted/step_wall = {100 * overcount_ratio:.2f}% (> 100.5%); "
            f"residual_ns={residual_ns}, parts={disjoint_ns}"
        )

    shares = {k: (v / step_wall if step_wall > 0 else 0.0) for k, v in disjoint_ns.items()}
    shares["residual"] = residual_ns / step_wall if step_wall > 0 else 0.0
    shares["formTangent_all"] = _sum_name(root_scope, "formTangent") / step_wall if step_wall > 0 else 0.0
    shares["sanisand_sub_ns"] = _sum_name(root_scope, "sanisand.sasME.substeps") / step_wall if step_wall > 0 else 0.0
    tag_shares = {t: (w / step_wall if step_wall > 0 else 0.0) for t, w in tag_ns.items()}
    return {
        "shares": shares,
        "tag_shares": tag_shares,
        "step_wall_ns": step_wall,
        "overcount_ratio": overcount_ratio,
    }


def _med_spread(vals: List[float]) -> Tuple[float, float, float]:
    s = sorted(vals)
    return statistics.median(s), s[0], s[-1]


def _dump_tree_str(node: Dict[str, Any], name: str = "root", pad: str = "") -> List[str]:
    w = node["wall_ns"] / 1e9
    eb = ", ".join(f"{t}:n={c}:{wns / 1e9:.4f}s" for t, c, wns in node.get("elem_by_type", []))
    eb_s = f" [{eb}]" if eb else ""
    lines = []
    if w >= 0.001 or not pad:
        lines.append(f"{pad}{name}: {w:.4f}s (c={node['calls']}){eb_s}")
    for k, ch in node["children"].items():
        lines.extend(_dump_tree_str(ch, k, pad + "  "))
    return lines


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", type=Path, default=None, help="temp directory containing .h5 + manifest.json")
    ap.add_argument("--from-json", type=Path, default=None, help="committed results/<sha9> dir (manifest.json + trees.json)")
    ap.add_argument("--export-json", type=Path, default=None, help="write manifest.json, trees.json, dump_trees.txt here")
    ap.add_argument("--out-results", type=Path, default=OUT_RESULTS)
    ap.add_argument("--out-comment", type=Path, default=None)
    ap.add_argument("--dump-trees", action="store_true")
    args = ap.parse_args()

    if args.from_json:
        src_dir = args.from_json.resolve()
        manifest = json.loads((src_dir / "manifest.json").read_text(encoding="utf-8"))
        trees = json.loads((src_dir / "trees.json").read_text(encoding="utf-8"))
        get_tree = lambda r: trees[r["run_id"]]["tree"]
    else:
        src_dir = (args.dir or DEFAULT_PROF_DIR).resolve()
        manifest = json.loads((src_dir / "manifest.json").read_text(encoding="utf-8"))
        trees_file = src_dir / "trees.json"
        cached_trees = json.loads(trees_file.read_text(encoding="utf-8")) if trees_file.exists() else {}

        def get_tree(r):
            rid = r["run_id"]
            if rid in cached_trees:
                return cached_trees[rid]["tree"]
            t = _load_h5_tree(Path(r["h5"]), rid)
            cached_trees[rid] = {"tree": t, "series_wall_ms": r.get("series_wall_ms", [])}
            return t

    if args.dump_trees:
        seen = set()
        for rec in manifest["records"]:
            if not rec.get("ok") or rec["mode"] != "1t_deep":
                continue
            key = (rec["deck"], rec.get("size", "S"))
            if key in seen:
                continue
            seen.add(key)
            print(f"=== {rec['deck']}[{rec.get('size', 'S')}] ({rec['wall_s']:.3f}s) ===")
            for ln in _dump_tree_str(get_tree(rec)):
                print(ln)
        return 0

    # Assert SHA consistency between OpenSees.exe and every opensees.pyd record
    exe_build = manifest.get("exe_build", "")
    ok_recs = [r for r in manifest["records"] if r.get("ok")]
    if not ok_recs:
        raise RuntimeError("No successful records in manifest.json")
    for r in ok_recs:
        if r.get("build_id") != exe_build:
            raise RuntimeError(f"Mixed build SHA! exe={exe_build} vs record {r['run_id']}={r.get('build_id')}")

    if args.export_json:
        exp = args.export_json.resolve()
        exp.mkdir(parents=True, exist_ok=True)
        for r in ok_recs:
            get_tree(r)
        (exp / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
        if not args.from_json:
            (exp / "trees.json").write_text(json.dumps(cached_trees, indent=1), encoding="utf-8")
        dt_lines = []
        seen = set()
        for r in ok_recs:
            if r["mode"] == "1t_deep" and (r["deck"], r.get("size", "S")) not in seen:
                seen.add((r["deck"], r.get("size", "S")))
                dt_lines.append(f"=== {r['deck']}[{r.get('size', 'S')}] ===")
                dt_lines.extend(_dump_tree_str(get_tree(r)))
        (exp / "dump_trees.txt").write_text("\n".join(dt_lines) + "\n", encoding="utf-8")

    # Group by (deck, size, mode)
    by_dsm: Dict[Tuple[str, str, str], List[Dict[str, Any]]] = {}
    for r in ok_recs:
        by_dsm.setdefault((r["deck"], r.get("size", "S"), r["mode"]), []).append(r)

    combos = sorted({(k[0], k[1]) for k in by_dsm})
    share_keys = [
        "soe.symbolic", "soe.factor", "soe.trisolve", "soe.cgs", "soe.addA",
        "elem_kernel", "mpm.a.prepare", "newStep_other", "unbal_other",
        "upd_other", "dc_other", "commit.elements", "commit.record",
        "residual", "formTangent_all", "sanisand_sub_ns",
    ]
    rows: List[Dict[str, Any]] = []
    unexplained_flags = []
    for d, sz in combos:
        r1d = by_dsm.get((d, sz, "1t_deep"), [])
        r1c = by_dsm.get((d, sz, "1t_coarse"), [])
        rnc = by_dsm.get((d, sz, "nt_coarse"), [])
        if not r1d or not r1c or not rnc:
            continue
        parsed1 = [analyze_tree(get_tree(r)) for r in r1d]
        med_sh = {k: statistics.median([p["shares"][k] for p in parsed1]) for k in share_keys}
        all_tags = {t for p in parsed1 for t in p["tag_shares"]}
        med_tags = {t: statistics.median([p["tag_shares"].get(t, 0.0) for p in parsed1]) for t in all_tags}
        tag_desc = ", ".join(
            f"{TAG_NAMES.get(t, f'tag{t}')}({t})={100.0 * sh:.1f}%"
            for t, sh in sorted(med_tags.items(), key=lambda kv: -kv[1])
            if sh >= 0.005
        ) or "-"

        w1d = _med_spread([r.get("steady_wall_s", r["wall_s"]) for r in r1d])
        w1c = _med_spread([r.get("steady_wall_s", r["wall_s"]) for r in r1c])
        wnc = _med_spread([r.get("steady_wall_s", r["wall_s"]) for r in rnc])
        deep_ovh = w1d[0] / w1c[0] if w1c[0] > 0 else 0.0
        scale = w1c[0] / wnc[0] if wnc[0] > 0 else 0.0

        # Flag unexplained remainders > 5%
        for fn in ("residual", "upd_other", "newStep_other", "unbal_other", "dc_other"):
            if med_sh[fn] > 0.05:
                unexplained_flags.append(f"`{d}[{sz}]` `{fn}` = {100 * med_sh[fn]:.1f}%")

        rows.append({
            "deck": d, "size": sz, "label": f"{d}[{sz}]", "family": r1d[0]["family"],
            "nsteps": r1d[0]["nsteps"], "nelm": r1d[0]["nelm"], "ndof": r1d[0].get("ndof_est", 0),
            "shares": med_sh, "tag_shares": med_tags, "tag_desc": tag_desc,
            "w1d": w1d, "w1c": w1c, "wnc": wnc,
            "deep_ovh": deep_ovh, "scale": scale,
            "threaded": rnc[0].get("threaded", "unknown"),
        })

    # Per-deck p/s/G matrix (Rule 2)
    matrix_rows = []
    for pkg, tier, rstatus, blocker, hyp, s_r2, src, p_desc, pfn in PACKAGE_SPECS:
        deck_cells = {}
        best_G = 1.0
        best_p = 0.0
        best_lbl = "-"
        # Rule 1: rank on L if present for that deck family, else S/smoke
        has_L = {r["deck"] for r in rows if r["size"] == "L"}
        for r in rows:
            if pfn is None or rstatus in ("not phase-rankable", "not eligible (Wave 3)"):
                deck_cells[r["label"]] = "-"
                continue
            p = max(0.0, min(1.0, float(pfn(r["shares"], r["tag_shares"], {"family": r["family"]}))))
            if s_r2 is not None and p > 0.0:
                g = 1.0 / ((1.0 - p) + (p / s_r2))
                deck_cells[r["label"]] = f"{100 * p:.1f}% / {g:.2f}x"
                # Prefer L deck for ranking selection when available
                is_rank_size = (r["size"] == "L") or (r["deck"] not in has_L)
                if is_rank_size and g > best_G:
                    best_G, best_p, best_lbl = g, p, r["label"]
            elif p > 0.0:
                deck_cells[r["label"]] = f"{100 * p:.1f}% / unranked"
            else:
                deck_cells[r["label"]] = "0.0%"
        passes = (rstatus == "eligible") and (best_p >= 0.10) and (best_G >= 1.10)
        matrix_rows.append({
            "pkg": pkg, "tier": tier, "status": rstatus, "blocker": blocker,
            "hyp": hyp, "s": s_r2, "source": src, "p_desc": p_desc,
            "best_p": best_p, "best_G": best_G, "best_lbl": best_lbl,
            "passes": passes, "cells": deck_cells,
        })

    # Format 5-part deliverable (`wp-briefs.md` §5)
    mode_tag = "SMOKE / PRE-#1153/#1154 VERIFICATION" if manifest.get("smoke") else "FULL BASELINE"
    md: List[str] = []
    md.append(f"# WP-260 Wave 0 (`260-w0`) — {mode_tag}")
    md.append("")
    md.append(f"- **Binary SHA (verified `OpenSees.exe` == all `opensees.pyd` records):** `{exe_build}`")
    md.append(f"- **Environment (`OMP_*`/`KMP_*`/`MKL_*`/`LADRUNO_*`):** `{json.dumps(manifest.get('parent_env', {}))}`")
    md.append(f"- **Affinity arg:** `{manifest.get('affinity_arg', 'none')}` | **Production threads:** `{manifest.get('production_threads', 8)}`")
    md.append("- **Step-1 drop:** steady-state wall times computed from `-perStep` `series/wall_ms[1:]`. `other/residual` is NEVER clamped (`sum > 100.5%` fails).")
    md.append("")

    # Part 1: Phase-share table
    md.append("## 1. Phase-Share Table (`1t_deep`) & 3-Mode Wall Times")
    md.append("")
    md.append("| Deck[Size] | DOF~ | Steps | `soe.sym` | `soe.fac` | `soe.tri` | `soe.cgs` | `soe.addA` | `elem_kern` | `mpm.a.prep` | `newStp_o` | `unbal_o` | `upd_o` | `dc_o` | `cmt_el` | `cmt_rec` | `residual` | `1t_deep` | `1t_coarse` | `nt_coarse` | `deep_ovh` | `1tc/ntc` | `threaded?` |")
    md.append("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|")
    for r in rows:
        sh = r["shares"]
        md.append(
            f"| `{r['label']}` | {r['ndof']} | {r['nsteps']} "
            f"| {100 * sh['soe.symbolic']:.1f}% | {100 * sh['soe.factor']:.1f}% | {100 * sh['soe.trisolve']:.1f}% "
            f"| {100 * sh['soe.cgs']:.1f}% | {100 * sh['soe.addA']:.1f}% | {100 * sh['elem_kernel']:.1f}% "
            f"| {100 * sh['mpm.a.prepare']:.1f}% | {100 * sh['newStep_other']:.1f}% | {100 * sh['unbal_other']:.1f}% "
            f"| {100 * sh['upd_other']:.1f}% | {100 * sh['dc_other']:.1f}% | {100 * sh['commit.elements']:.1f}% "
            f"| {100 * sh['commit.record']:.1f}% | {100 * sh['residual']:+.1f}% "
            f"| {r['w1d'][0]:.3f}s [{r['w1d'][1]:.3f}..{r['w1d'][2]:.3f}] "
            f"| {r['w1c'][0]:.3f}s [{r['w1c'][1]:.3f}..{r['w1c'][2]:.3f}] "
            f"| {r['wnc'][0]:.3f}s [{r['wnc'][1]:.3f}..{r['wnc'][2]:.3f}] "
            f"| {r['deep_ovh']:.2f}x | {r['scale']:.2f}x | {r['threaded']} |"
        )
    md.append("")
    md.append("**`elem_by_type` breakdown:** " + " ; ".join(f"`{r['label']}`: {r['tag_desc']}" for r in rows))
    if unexplained_flags:
        md.append("")
        md.append("> **Remainders > 5 % (see `dump_trees.txt`):** " + ", ".join(unexplained_flags))
    md.append("")

    # Part 2: Size trend (S vs L)
    md.append("## 2. Size Trend (`p` at S vs L)")
    md.append("")
    by_deck = {}
    for r in rows:
        by_deck.setdefault(r["deck"], {})[r["size"]] = r
    trend_rows = [d for d, smap in sorted(by_deck.items()) if "S" in smap and "L" in smap]
    if trend_rows:
        md.append("| Deck | `p_solver` S (`sym+fac+tri+cgs`) | `p_solver` L | `p_elem` S (`elem_kern`) | `p_elem` L |")
        md.append("|---|---:|---:|---:|---:|")
        for d in trend_rows:
            s_r, l_r = by_deck[d]["S"], by_deck[d]["L"]
            ps_s = sum(s_r["shares"][k] for k in ("soe.symbolic", "soe.factor", "soe.trisolve", "soe.cgs"))
            ps_l = sum(l_r["shares"][k] for k in ("soe.symbolic", "soe.factor", "soe.trisolve", "soe.cgs"))
            md.append(f"| `{d}` | {100 * ps_s:.1f}% | {100 * ps_l:.1f}% | {100 * s_r['shares']['elem_kernel']:.1f}% | {100 * l_r['shares']['elem_kernel']:.1f}% |")
    else:
        md.append("_Only one size present in this run (smoke/partial run); full S vs L table prints once S and L both run._")
    md.append("")

    # Part 3: Phase -> package map + gaps
    md.append("## 3. Phase -> Package Map & Coverage Gaps")
    md.append("")
    md.append("- **`soe.{symbolic,factor,trisolve,cgs}`** -> `260a` (`soe.factor` on Pardiso via opt-in `iparm`; `soe.addA`/`dc.setSize` on MUMPS only), `260g` (mixed-precision IR, unranked), `260m` (AMG, Wave 3).")
    md.append("- **`soe.addA`** -> `260a` (MUMPS only), `260f` (persistent assembly/fast scatter, unranked).")
    md.append("- **`elem_by_type`** -> `260c1` (`LogStrain` finite-strain, N/A in small-strain D1), `260c2` (`LadrunoConcrete3D` on D1b, blocked `#1154`), `260h` (analytic return maps, unranked & blocked `#1153/#1154`), `260k` (Wave 3).")
    md.append("- **Explicit `LadrunoCPML` + `commit.elements` + `ExplicitBathe formTangent` + SMS `soe.trisolve`** -> `260b` (`Unblocked`).")
    md.append("- **`mpm.a.prepare`** -> MPM host work (program `#1075`, tracked cleanly in its own column).")
    md.append("- **`commit.record`** -> recorders (outside WP-260; none >= 10% here).")
    md.append("- **Coverage gaps >= 10 %:** D5 `LadrunoMCFTStrip` fibre/section kernel (`elem_kern` on D5 — no package in WP-260 targets it); `TransformationConstraintHandler::enforceSPs()` double-update (`upd_other` / `newStep_other` on D1b and explicit decks).")
    md.append("")

    # Part 4: Full per-deck p / s / G matrix (brief §4)
    md.append("## 4. Plan §0.5 + Brief §4 Per-Deck Ranking (`p >= 10 %`, `G >= 1.10x`, rule-2 `s`)")
    md.append("")
    deck_hdr = " | ".join(f"`{r['label']}`" for r in rows)
    deck_sep = " | ".join("---:" for _ in rows)
    md.append(f"| Pkg | Tier | Status / Blocker | Plan Hypothesis | Rule-2 `s` | Source | {deck_hdr} | Peak L `p / G` | Gate |")
    md.append(f"|---|---|---|---|---:|---|{deck_sep}|---:|:---:|")
    for mr in matrix_rows:
        s_str = f"{mr['s']:.2f}x" if mr["s"] is not None else "—"
        cells = " | ".join(mr["cells"].get(r["label"], "-") for r in rows)
        if mr["status"] in ("not phase-rankable", "unranked", "not eligible (Wave 3)"):
            gate = mr["status"]
            peak = "—"
        else:
            gate = "PASS" if mr["passes"] else "BELOW"
            peak = f"{100 * mr['best_p']:.1f}% / **{mr['best_G']:.2f}x** (`{mr['best_lbl']}`)" if mr["best_lbl"] != "-" else "—"
        md.append(
            f"| `{mr['pkg']}` | {mr['tier']} | {mr['blocker']} | {mr['hyp']} | {s_str} | {mr['source']} "
            f"| {cells} | {peak} | {gate} |"
        )
    md.append("")

    # Part 5: Proposed pick (ties: G-BYTE > G-EQUIV > G-OPT, then G descending)
    tier_rank = {"G-BYTE": 0, "G-BYTE+G-OPT": 1, "G-BYTE/G-EQUIV": 2, "G-EQUIV/G-OPT": 3, "G-OPT": 4}
    md.append("## 5. Proposed Pick")
    md.append("")
    passing = sorted(
        [m for m in matrix_rows if m["passes"]],
        key=lambda m: (-round(m["best_G"], 4), tier_rank.get(m["tier"], 9)),
    )
    if passing:
        unblocked = [m for m in passing if m["blocker"] == "Unblocked"]
        if unblocked:
            md.append("**Immediate unblocked picks (eligible now):**")
            for m in unblocked[:3]:
                md.append(f"- **`{m['pkg']}`** ({m['tier']}, blocker: {m['blocker']}): peak on `{m['best_lbl']}`, `p = {100 * m['best_p']:.1f}%`, rule-2 `s = {m['s']:.2f}x` -> **`G = {m['best_G']:.2f}x`**.")
        md.append("")
        md.append("**All passing eligible packages (with blockers):**")
        for m in passing[:3]:
            md.append(f"- **`{m['pkg']}`** ({m['tier']}, blocker: {m['blocker']}): peak on `{m['best_lbl']}`, `p = {100 * m['best_p']:.1f}%`, rule-2 `s = {m['s']:.2f}x` -> **`G = {m['best_G']:.2f}x`**.")
    else:
        md.append("- No eligible package cleared both `p >= 10%` and `G >= 1.10x` on this run.")

    text = "\n".join(md) + "\n"
    args.out_results.write_text(text, encoding="utf-8")
    print(f"Wrote {args.out_results}")
    if args.out_comment:
        args.out_comment.parent.mkdir(parents=True, exist_ok=True)
        args.out_comment.write_text(text, encoding="utf-8")
        print(f"Wrote {args.out_comment}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
