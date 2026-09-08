"""ADR-95 cross-check #3: does LadrunoSANISAND wall on the note-95 Prandtl-Reissner
strip-footing deck at the SAME trigger (first footing-edge GPs reaching mean stress
>= 0, i.e. p -> 0 in SANISAND's own compression-positive convention) that killed the
vanilla UW `DruckerPrager` (dead cutoff branch, fixed) and `ASDPlasticMaterial3D`
(unimplemented apex return, open ADR-94)?

Reuses `h20_prandtl.py`'s mesh, consistent Q8/Q4 surcharge, and ADR-63 D16 step
ladders UNCHANGED (imported, not copied), exactly as `asd_path_diag.py` does --
`build_model_sanisand()` mirrors `build_model_asd()` line for line except for the
material assignment (LadrunoSANISAND, not ASDPlasticMaterial3D) and ONE extra
fixity: the footing here is ROUGH (u_x = 0 on footing nodes, the ADR-92 CP1 deck's
own footing convention -- `h20_prandtl`'s footing is smooth and that is NOT reused).

Material command and staging are the ADR-92 CP1 deck's (`sanisand_tau0_band.py`),
reused verbatim: Gorini's calibrated `_PARAMS` (`tests/test_ladruno_sanisand.py`,
ADR-86 sec.5, e_init = 0.6944), the positional optionals `IntScheme=1 TanType=2
JacoType=1 TolF=TolR=1e-7`, the flags `-Presidual 0 -Pmin 1e-3*P_atm -honorTolR 0
-maxSubsteps 0` (deck default: UNCAPPED), `updateMaterialStage -stage 0` (elastic)
before the surcharge ramp and `-stage 1` (plastic) before the push, with the
eta_max/M_c < 1 guard at the flip.

Does NOT modify h20_prandtl.py, asd_path_diag.py, quad_path_diag.py, or
sanisand_tau0_band.py -- other runs depend on those files being byte-stable.

q_exact: NONE quoted. SANISAND under self weight has never produced a capacity in
the ADR-92 CP1 campaign (`_adr92_cp1_surcharge_results.md`: every leg WALL/FLOOR, no
plateau, no peak) and this deck is weightless with a different (smaller) domain, so
no number from CP1 is comparable. Only q (kPa), termination mode, and the tension /
refusal census are reported.

Run:
    py -3.12 sanisand_path_diag.py --elem h8bbar --sfrac 0.15 --tmax 5400 --suffix _sanisand
    py -3.12 sanisand_path_diag.py --elem h20uri --sfrac 0.15 --tmax 5400 --suffix _sanisand
    py -3.12 sanisand_path_diag.py --elem h20uri --sfrac 0.15 --tmax 5400 --implex --suffix _sanisand_implex
"""
import argparse
import csv
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import h20_prandtl as HP                                        # noqa: E402

ops = HP.ops

_REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
_TESTS = os.path.join(_REPO, "tests")
if _TESTS not in sys.path:
    sys.path.append(_TESTS)

from test_ladruno_sanisand import _PARAMS as _SANISAND_PARAMS   # noqa: E402

ELEMS = {"h8bbar": dict(order=1, form="bbar", label="LadrunoBrick -formulation bbar"),
         "h20uri": dict(order=2, form="uri",  label="LadrunoBrick20 -formulation uri")}
GP_COUNTS = {"h8bbar": 8, "h20uri": 8}
COND_TRIGGER = 0.25
TENSION_EVERY = 10           # converged steps between tension-GP census rows

# --- the ADR-92 CP1 deck's material staging, reused verbatim -----------------
E_INIT = _SANISAND_PARAMS[2]                  # Gorini's calibrated 0.6944
NU_SAN = _SANISAND_PARAMS[1]                  # 0.3129
M_C = _SANISAND_PARAMS[3]                     # 1.33090
P_ATM = _SANISAND_PARAMS[8]                   # 101.0 kPa

INT_SCHEME = 1                 # ModifiedEuler with substepping (the default)
TAN_TYPE = 2                   # consistent elastoplastic -- NOT the parser default 0
JACO_TYPE = 1
TOL_F, TOL_R = 1.0e-7, 1.0e-7   # the parser defaults, emitted explicitly
OPT_PRESIDUAL = 0.0
OPT_PMIN = 1.0e-3 * P_ATM       # ADR-86 default, explicitly emitted

N_GRAV = 10                     # elastic-stage ramp steps (mirrors the ADR-92 deck)
GRAV_TOL_REL = 1.0e-9


def log(*a):
    print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)


# ---------------------------------------------------------------------------
# the model: h20_prandtl.build_model(), verbatim geometry, material swapped for
# LadrunoSANISAND + a rough footing (u_x = 0 added on the footing nodes).
# ---------------------------------------------------------------------------
def build_model_sanisand(nodes, cells, sets, w, form, max_substeps, implex):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for i, (x, y, z) in enumerate(nodes, start=1):
        ops.node(i, float(x), float(y), float(z))

    par = list(_SANISAND_PARAMS)
    ops.nDMaterial(
        "LadrunoSANISAND", 1, *par,
        INT_SCHEME, TAN_TYPE, JACO_TYPE, TOL_F, TOL_R,
        "-Presidual", OPT_PRESIDUAL, "-Pmin", OPT_PMIN,
        "-honorTolR", 0, "-maxSubsteps", int(max_substeps),
        *(("-implex",) if implex else ()),
    )

    for e, conn in enumerate(cells, start=1):
        if cells.shape[1] == 20:
            ops.element("LadrunoBrick20", e, *[int(c) + 1 for c in conn], 1,
                        "-formulation", form, "-b", 0.0, 0.0, -HP.GAMMA)
        else:
            ops.element("LadrunoBrick", e, *[int(c) + 1 for c in conn], 1,
                        "-geom", "linear", "-b", 0.0, 0.0, -HP.GAMMA,
                        "-formulation", form)

    # ONE fix() per node. u_y = 0 everywhere (plane strain). The footing is
    # ROUGH here (u_x = 0 on footing nodes too, the ADR-92 CP1 convention) --
    # h20_prandtl's own footing is smooth; that fixity is NOT reused.
    bt = set(sets["bottom"].tolist())
    xf = set(sets["xface"].tolist())
    ft = set(sets["footing"].tolist())
    for n in range(len(nodes)):
        if n in bt:
            ops.fix(n + 1, 1, 1, 1)
        else:
            ops.fix(n + 1, 1 if (n in xf or n in ft) else 0, 1, 0)

    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for n in sets["top"]:
        if abs(w[n]) > 0:
            ops.load(int(n) + 1, 0.0, 0.0, -HP.Q0 * float(w[n]))
    return NU_SAN


def _dev_and_p(st):
    """(p, q, eta) in SANISAND's own COMPRESSION-positive convention, matching
    `sanisand_tau0_band.py::_dev_and_p` exactly (element/material stress `st`
    is tension-positive; p = -trace/3)."""
    p = -(st[0] + st[1] + st[2]) / 3.0
    s = (st[0] + p, st[1] + p, st[2] + p, st[3], st[4], st[5])
    j2 = 0.5 * (s[0] ** 2 + s[1] ** 2 + s[2] ** 2) + s[3] ** 2 + s[4] ** 2 + s[5] ** 2
    q = math.sqrt(3.0 * j2)
    return p, q, (q / p if p > 1.0e-9 else float("nan"))


def count_tension_gps(ncells, ngp):
    """Mean stress (TENSION-positive) >= 0 census, copied unmodified from
    `asd_path_diag.py`'s `count_tension_gps`, PLUS the per-GP substep/cap-hit
    counters LadrunoSANISAND exposes (`eleResponse(..., "substeps")` ->
    [substepsTakenInME, capHit]; capHit is 0 by construction at the deck's
    uncapped default and reported anyway so that is explicit, not assumed)."""
    n_seen = n_tension = n_nan = n_cap_hit = 0
    first_nan = None
    i1_min, i1_max = float("inf"), float("-inf")
    substeps_max = 0
    for e in range(1, ncells + 1):
        for gp in range(1, ngp + 1):
            try:
                resp = ops.eleResponse(e, "material", gp, "stress")
            except Exception:
                resp = None
            if resp is None or len(resp) < 6:
                continue
            n_seen += 1
            s = np.asarray(resp[:6], dtype=float)
            if not np.all(np.isfinite(s)):
                n_nan += 1
                if first_nan is None:
                    first_nan = (e, gp)
                continue
            mean = float(s[0] + s[1] + s[2]) / 3.0
            i1 = float(s[0] + s[1] + s[2])
            i1_min = min(i1_min, i1)
            i1_max = max(i1_max, i1)
            if mean >= 0.0:
                n_tension += 1
            try:
                sub = ops.eleResponse(e, "material", gp, "substeps")
            except Exception:
                sub = None
            if sub is not None and len(sub) >= 2:
                substeps_max = max(substeps_max, int(round(sub[0])))
                if sub[1] > 0.5:
                    n_cap_hit += 1
    return dict(n_seen=n_seen, n_tension=n_tension, n_nan=n_nan,
                first_nan=first_nan, i1_min=i1_min, i1_max=i1_max,
                substeps_max=substeps_max, n_cap_hit=n_cap_hit)


def implex_refusals():
    """Process-wide `LadrunoImplexGlobals` refusal ledger. Only meaningful on
    the `--implex` leg; returns zeros (harmlessly) if the response is not
    armed. Cumulative for the process, which is fine here since one leg runs
    per invocation (the ADR-92 deck's own caveat)."""
    try:
        v = ops.eleResponse(1, "material", 1, "implexRefusals")
        return dict(total=int(round(v[0])), d2=int(round(v[1])),
                    control=int(round(v[2])), companion=int(round(v[3])))
    except Exception:
        return dict(total=0, d2=0, control=0, companion=0)


def run(args):
    spec = ELEMS[args.elem]
    order, form = spec["order"], spec["form"]
    base, floor, dmax = HP.DS_LADDER[order]
    if args.ds_base is not None:
        base = args.ds_base
    if args.floor is not None:
        floor = args.floor
    if args.ds_max is not None:
        dmax = args.ds_max
    ladder_kind = "quad" if order == 2 else "linear"
    cond_at = COND_TRIGGER * base

    stamp = ops.ladrunoBuild()
    log(f"[control 0] engine {ops.__file__}")
    log(f"[control 0] ladrunoBuild() = {stamp}")

    nodes, cells, sets, xg, yg, zg = HP.strip_mesh(args.h0, order)
    w, nfaces = HP.consistent_surcharge(nodes, cells, order)
    HP.verify_surcharge(nodes, cells, w, nfaces, order)

    tag = f"{args.elem}_h{str(args.h0).replace('0.', '')}{args.suffix}"
    log(f"=== leg {tag}: {spec['label']} on LadrunoSANISAND "
        f"(Gorini e_init={E_INIT}, IntScheme={INT_SCHEME} TanType={TAN_TYPE} "
        f"maxSubsteps={args.maxsubsteps}{' -implex' if args.implex else ' (implicit)'}), "
        f"h0 = {args.h0} m, {len(nodes)} nodes / {3 * len(nodes)} DOF, "
        f"{len(cells)} elements, ROUGH footing")

    build_model_sanisand(nodes, cells, sets, w, form, args.maxsubsteps, args.implex)

    # ---- stage 0: elastic, ramp the surcharge under it ----------------------
    ops.updateMaterialStage("-material", 1, "-stage", 0)
    tol = 1.0e-5 * max(HP.Q0 * float(w.sum()), 1.0)
    HP.surcharge_step(nodes, sets, w, tol)

    # --- control 1: 1-D elastic patch test, SANISAND's OWN K0 (nu = 0.3129) --
    k0 = NU_SAN / (1.0 - NU_SAN)
    e_zz = e_xx = 0.0
    eta_max = 0.0
    for e in range(1, len(cells) + 1):
        for gp in range(1, GP_COUNTS[args.elem] + 1):
            st = ops.eleResponse(e, "material", gp, "stress")
            if not st or len(st) < 6:
                continue
            st = np.asarray(st[:6], dtype=float)
            e_zz = max(e_zz, abs(st[2] / HP.Q0 + 1.0))
            e_xx = max(e_xx, abs(st[0] / (k0 * HP.Q0) + 1.0))
            _, _, eta = _dev_and_p(st)
            if eta == eta:
                eta_max = max(eta_max, eta)
    patch_err = max(e_zz, e_xx)
    log(f"    [control 1] 1-D elastic stress patch (SANISAND K0={k0:.4f}): "
        f"max|szz/-q0-1| = {e_zz:.2e}, max|sxx/-K0q0-1| = {e_xx:.2e}  "
        f"-> {'PASS' if patch_err < 1e-6 else '*** FAIL ***'}")
    if patch_err >= 1e-6:
        raise SystemExit("patch test failed; every number downstream is void")
    log(f"    [control gravity-analogue] eta_max/M_c at the stage flip = "
        f"{eta_max / M_C:.4f}  -> {'PASS' if eta_max < M_C else '*** FAIL: '
        'Elastic2Plastic will inflate the calibrated friction constant ***'}")
    if eta_max >= M_C:
        raise SystemExit("eta_max >= M_c at the stage flip; not Gorini's soil")

    # ---- stage flip: plastic ------------------------------------------------
    ops.updateMaterialStage("-material", 1, "-stage", 1)

    # ---------------- the push ------------------------------------------------
    foot = [int(n) + 1 for n in sets["footing"]]
    r_corr = HP.Q0 * float(w[sets["footing"]].sum())
    area = HP.B_FOOT * HP.THICK
    uz0 = ops.nodeDisp(foot[0], 3)
    smax = args.sfrac * HP.B_FOOT
    ops.loadConst("-time", uz0)
    ops.timeSeries("Linear", 2)
    ops.pattern("Plain", 2, 2)
    for t in foot:
        ops.sp(t, 3, 1.0)

    def setup(sysname):
        ops.wipeAnalysis()
        ops.constraints("Transformation")
        ops.numberer("Plain" if sysname == "FullGeneral" else "RCM")
        if sysname == "FullGeneral":
            ops.system("FullGeneral")
        else:
            try:
                ops.system("Pardiso", "-matrixType", 0)   # unsymmetric flow rule
            except Exception:
                ops.system("UmfPack")
        ops.analysis("Static")

    setup("sparse")
    ladder = ([("KrylovNewton", tol, 25, 0),
               ("NewtonLineSearch", tol, 40, 0),
               ("KrylovNewton", 10.0 * tol, 60, 1)]
              if ladder_kind == "quad" else
              [("Newton", tol, 25, 0),
               ("NewtonLineSearch", tol, 40, 0),
               ("KrylovNewton", 10.0 * tol, 60, 1)])

    path = os.path.join(HERE, f"sanisand_{tag}.csv")
    fh = open(path, "w", newline="")
    wr = csv.writer(fh)
    wr.writerow(["s_m", "s_over_B", "q_kPa", "ds_mm", "relaxed", "wall_s"])

    tpath = os.path.join(HERE, f"sanisand_{tag}_tension.csv")
    tfh = open(tpath, "w", newline="")
    twr = csv.writer(tfh)
    twr.writerow(["s_over_B", "q_kPa", "n_seen", "n_tension", "n_nan",
                  "first_nan_ele", "first_nan_gp", "I1_min", "I1_max",
                  "substeps_max", "n_cap_hit",
                  "implex_total", "implex_d2", "implex_control", "implex_companion"])

    rows, ds, good, nfail, nrelax, nsub = [], base, 0, 0, 0, 0
    mode, verdict = "TARGET", "reached the target settlement"
    stall_window = HP.STALL_WINDOW
    stall_advance = HP.STALL_ADVANCE
    grow_after = HP.GROW_AFTER
    nan_event, tension_event = None, None

    t0 = time.time()
    while True:
        s_now = uz0 - ops.getTime()
        if s_now >= smax - 1e-12:
            break
        if time.time() - t0 > args.tmax:
            mode, verdict = "WALL", f"wall-clock cap at s/B = {s_now / HP.B_FOOT:.5f}"
            break
        ds = min(ds, smax - s_now)
        ops.integrator("LoadControl", -ds)
        ok, relaxed = False, 0
        for algo, tl, it, rl in ladder:
            ops.test("NormUnbalance", tl, it, 0)
            ops.algorithm(algo)
            rc = ops.analyze(1)
            if rc == 0:
                ok, relaxed = True, rl
                break
            nfail += 1
        if not ok:
            good, nsub = 0, nsub + 1
            ds *= 0.5
            if nsub > args.budget:
                mode = "BUDGET"
                verdict = f"subdivision budget of {args.budget} spent at s/B = {s_now / HP.B_FOOT:.5f}"
                break
            if ds < floor:
                mode = "FLOOR"
                verdict = (f"step collapsed to the {floor * 1e3:.4g} mm floor at "
                           f"s/B = {s_now / HP.B_FOOT:.5f} (every ladder rung failed)")
                break
            continue
        nrelax += relaxed
        good += 1
        if good >= grow_after and ds < dmax:
            ds, good = min(2 * ds, dmax), 0
        ops.reactions()
        q = (-sum(ops.nodeReaction(t, 3) for t in foot) + r_corr) / area
        s = uz0 - ops.getTime()
        rows.append((s, s / HP.B_FOOT, q, ds * 1e3, relaxed, time.time() - t0))
        wr.writerow([f"{v:.9g}" for v in rows[-1]])
        fh.flush()

        if len(rows) % TENSION_EVERY == 0:
            tc = count_tension_gps(len(cells), GP_COUNTS[args.elem])
            ir = implex_refusals() if args.implex else dict(total=0, d2=0, control=0, companion=0)
            twr.writerow([f"{s / HP.B_FOOT:.9g}", f"{q:.6g}", tc["n_seen"],
                          tc["n_tension"], tc["n_nan"],
                          tc["first_nan"][0] if tc["first_nan"] else "",
                          tc["first_nan"][1] if tc["first_nan"] else "",
                          f"{tc['i1_min']:.6g}", f"{tc['i1_max']:.6g}",
                          tc["substeps_max"], tc["n_cap_hit"],
                          ir["total"], ir["d2"], ir["control"], ir["companion"]])
            tfh.flush()
            log(f"    [tension] s/B {s / HP.B_FOOT:.5f}  q {q:8.2f}  "
                f"GP seen {tc['n_seen']}  tension(mean>=0) {tc['n_tension']}  "
                f"NaN {tc['n_nan']}  I1 [{tc['i1_min']:.3g},{tc['i1_max']:.3g}]  "
                f"substeps_max {tc['substeps_max']} capHit {tc['n_cap_hit']}"
                + (f"  implexRefusals total={ir['total']} companion={ir['companion']}"
                   if args.implex else ""))
            if tc["n_tension"] and tension_event is None:
                tension_event = dict(s_over_B=s / HP.B_FOOT, q=q, n=tc["n_tension"])
                log(f"    *** first p>=0 (tension) GP census *** at s/B="
                    f"{tension_event['s_over_B']:.5f} q={tension_event['q']:.3f} "
                    f"n_tension={tension_event['n']}")
            if tc["n_nan"] and nan_event is None:
                nan_event = dict(s_over_B=s / HP.B_FOOT, q=q, ele=tc["first_nan"][0],
                                  gp=tc["first_nan"][1])
                log(f"    *** NaN event *** first at s/B={nan_event['s_over_B']:.5f} "
                    f"q={nan_event['q']:.3f} ele={nan_event['ele']} gp={nan_event['gp']}")

        if len(rows) > stall_window and \
                rows[-1][0] - rows[-1 - stall_window][0] < stall_advance * smax:
            mode = "STALL"
            verdict = f"stalled at s/B = {s / HP.B_FOOT:.5f}"
            break
    fh.close()
    tfh.close()
    wall = time.time() - t0
    log(f"    MODE = {mode}   [{verdict}]   wall = {wall:.1f}s  "
        f"converged steps = {len(rows)}  nsub = {nsub}  nfail = {nfail}")
    if rows:
        q_arr = np.array([r[2] for r in rows])
        s_arr = np.array([r[1] for r in rows])
        log(f"    last s/B = {s_arr[-1]:.5f}  q_end = {q_arr[-1]:.4f} kPa  "
            f"q_max = {q_arr.max():.4f} kPa   (NO exact reference for SANISAND; "
            f"do not quote a ratio)")
    if tension_event:
        log(f"    SUMMARY: first p>=0 GP census at s/B={tension_event['s_over_B']:.5f}")
    if nan_event:
        log(f"    SUMMARY: NaN event recorded at s/B={nan_event['s_over_B']:.5f}, "
            f"ele={nan_event['ele']} gp={nan_event['gp']}")
    return dict(tag=tag, mode=mode, verdict=verdict, rows=rows,
                tension_event=tension_event, nan_event=nan_event,
                nsub=nsub, nfail=nfail, wall=wall)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--elem", default="h20uri", choices=sorted(ELEMS))
    ap.add_argument("--h0", type=float, default=1.0)
    ap.add_argument("--sfrac", type=float, default=0.15)
    ap.add_argument("--budget", type=int, default=200)
    ap.add_argument("--ds-base", type=float, default=None)
    ap.add_argument("--floor", type=float, default=None)
    ap.add_argument("--ds-max", type=float, default=None)
    ap.add_argument("--tmax", type=float, default=5400.0)
    ap.add_argument("--maxsubsteps", type=int, default=0,
                    help="LadrunoSANISAND -maxSubsteps; 0 = uncapped (the ADR-92 deck default)")
    ap.add_argument("--implex", action="store_true",
                    help="run with -implex (IMPL-EX) instead of the default implicit scheme")
    ap.add_argument("--suffix", default="")
    args = ap.parse_args()
    r = run(args)
    print()
    print(f"{r['tag']:>28} mode={r['mode']:<9} verdict={r['verdict']}")


if __name__ == "__main__":
    main()
