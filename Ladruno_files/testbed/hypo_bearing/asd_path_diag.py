"""ADR-95 cross-check: is the quadratic-wall finding specific to the vanilla UW
`DruckerPrager` material, or does the fork's OTHER Drucker-Prager --
`ASDPlasticMaterial3D` with `DruckerPrager_YF`/`DruckerPrager_PF` -- reproduce
it on the identical Prandtl-Reissner strip-footing deck?

Reuses h20_prandtl.py's mesh, consistent-surcharge, and push-loop machinery
UNCHANGED (imported, not copied) so the deck is bit-for-bit the one note 95
measured on the vanilla material.  Only the material assignment is new: a
`build_model_asd()` that mirrors `h20_prandtl.build_model()` line for line
except for the `nDMaterial` call, which is `ASDPlasticMaterial3D` instead of
`DruckerPrager`.

Does NOT modify quad_path_diag.py or h20_prandtl.py -- other runs depend on
those files being byte-stable.

PARAMETER MAPPING (stated once, used once)
-------------------------------------------
UW `DruckerPrager`'s cone (note 95's own form, `quad_path_diag.mobilisation`):
    f1_UW = ||dev sigma|| + rho * I1 - sqrt(2/3) * SY        (I1 = tr(sigma),
                                                               tension-positive)
ASD `DruckerPrager_YF`'s cone (`YieldFunctions/DruckerPrager_YF.h`):
    f_ASD = sqrt(J2) + eta * p - xi_c                        (p = tr(sigma)/3,
                                              ALSO tension-positive: R4's
                                              adjudication in adr94_verdict.md
                                              says the DP_YF header comment
                                              claiming compression-positive is
                                              WRONG, meanStress() = trace/3 is
                                              tension-positive throughout ASDP)
    sqrt(J2) = ||dev sigma|| / sqrt(2)   =>   f_ASD = 0  <=>
        ||dev sigma|| = sqrt(2)*xi_c - (sqrt(2)/3)*eta*I1
Matching the two cones' ||dev sigma|| vs I1 lines coefficient-by-coefficient:
    -rho = -(sqrt(2)/3)*eta        =>  eta   = (3/sqrt(2)) * rho = 3 * alpha
                                        (rho = sqrt(2)*alpha, so the sqrt(2)s
                                        cancel exactly)
    sqrt(2/3)*SY = sqrt(2)*xi_c     =>  xi_c  = SY / sqrt(3)
Apex check (||dev||=0): UW's T = sqrt(2/3)*SY/rho; ASD's I1 at f=0, sqrt(J2)=0
is 3*xi_c/eta = 3*(SY/sqrt(3)) / (3*alpha) = SY/(sqrt(3)*alpha) = SY*sqrt(2/3)/rho
(since rho = sqrt(2)*alpha) -- identical apex.  Non-associated psi=0 => the UW
leg uses rhoBar=0 (no volumetric term in the flow potential); the ASD analogue
is `DP_etabar = 0` (PF's pressure_part is etabar/3 * I, same role as rhoBar).
Cohesion/back-stress HARDENING is turned off (TensorLinearHardeningParameter =
ScalarLinearHardeningParameter = 0) to match perfect plasticity; the internal
variables (BackStress, DP_cohesion) still have to be declared and given an
initial value (0) because the parser requires them, but with H=0 their trial
value never moves and DP_cohesion's value is dead code in DruckerPrager_YF's
YIELD_FUNCTION (xi_c is read as a fixed model PARAMETER, not from the IV).

E, nu: E = 9*K*G/(3*K+G) with K = h20_prandtl.K_EL, G = h20_prandtl.g_el at
nu = h20_prandtl.NU -- the same E,nu <-> K,G identity h20_prandtl.py itself
uses to build its "elastic control" material, so ASD's Lame constants (built
from E, nu) reproduce the vanilla leg's K, G exactly.

Run:
    py -3.12 asd_path_diag.py --elem h8bbar --h0 1.0 --sfrac 0.15 --tmax 5400 --suffix _asd
    py -3.12 asd_path_diag.py --elem h20uri --h0 1.0 --sfrac 0.15 --tmax 5400 --suffix _asd
"""
import argparse
import csv
import math
import os
import sys
import time
from collections import deque

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import h20_prandtl as HP                                        # noqa: E402

ops = HP.ops

SQ23 = math.sqrt(2.0 / 3.0)

ELEMS = {"h8bbar": dict(order=1, form="bbar", label="LadrunoBrick -formulation bbar"),
         "h20uri": dict(order=2, form="uri",  label="LadrunoBrick20 -formulation uri")}
GP_COUNTS = {"h8bbar": 8, "h20uri": 8}
COND_TRIGGER = 0.25
TENSION_EVERY = 10           # converged steps between tension-GP census rows


def log(*a):
    print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)


# ---------------------------------------------------------------------------
# the model: h20_prandtl.build_model(), verbatim, with the material swapped
# ---------------------------------------------------------------------------
def build_model_asd(nodes, cells, sets, w, form):
    alpha = HP.alpha_from_phi_txc(HP.PHI_TXC)
    rho = math.sqrt(2.0) * alpha                      # UW's rho, for reporting
    g_el = 3.0 * HP.K_EL * (1.0 - 2.0 * HP.NU) / (2.0 * (1.0 + HP.NU))
    e_el = 9.0 * HP.K_EL * g_el / (3.0 * HP.K_EL + g_el)

    eta_asd = 3.0 * alpha
    xi_c_asd = HP.SY / math.sqrt(3.0)
    etabar_asd = 0.0                                   # psi = 0

    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for i, (x, y, z) in enumerate(nodes, start=1):
        ops.node(i, float(x), float(y), float(z))

    iv_type = "BackStress(TensorLinearHardeningFunction):DP_cohesion(ScalarLinearHardeningFunction):"
    ops.nDMaterial(
        "ASDPlasticMaterial3D", 1,
        "DruckerPrager_YF", "DruckerPrager_PF", "LinearIsotropic3D_EL", iv_type,
        "Begin_Internal_Variables",
        "BackStress", 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        "DP_cohesion", 0.0,
        "End_Internal_Variables",
        "Begin_Model_Parameters",
        "YoungsModulus", e_el,
        "PoissonsRatio", HP.NU,
        "DP_xi_c", xi_c_asd,
        "DP_eta", eta_asd,
        "DP_etabar", etabar_asd,
        "TensorLinearHardeningParameter", 0.0,
        "ScalarLinearHardeningParameter", 0.0,
        "End_Model_Parameters",
        "Begin_Integration_Options",
        "integration_method", "Backward_Euler",
        "tangent_type", "Continuum",
        "n_max_iterations", 100,
        "strict_convergence", 1,
        "End_Integration_Options",
    )

    for e, conn in enumerate(cells, start=1):
        if cells.shape[1] == 20:
            ops.element("LadrunoBrick20", e, *[int(c) + 1 for c in conn], 1,
                        "-formulation", form, "-b", 0.0, 0.0, -HP.GAMMA)
        else:
            ops.element("LadrunoBrick", e, *[int(c) + 1 for c in conn], 1,
                        "-geom", "linear", "-b", 0.0, 0.0, -HP.GAMMA,
                        "-formulation", form)

    bt, xf = set(sets["bottom"].tolist()), set(sets["xface"].tolist())
    for n in range(len(nodes)):
        if n in bt:
            ops.fix(n + 1, 1, 1, 1)
        else:
            ops.fix(n + 1, 1 if n in xf else 0, 1, 0)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for n in sets["top"]:
        if abs(w[n]) > 0:
            ops.load(int(n) + 1, 0.0, 0.0, -HP.Q0 * float(w[n]))
    return alpha, rho, eta_asd, xi_c_asd, e_el


def tangent_health():
    a = np.asarray(ops.printA("-ret"), dtype=float)
    n = int(round(math.sqrt(a.size)))
    if n * n != a.size:
        return None
    K = a.reshape(n, n)
    scale = float(np.trace(K)) / n
    sv = np.linalg.svd(K, compute_uv=False)
    S = 0.5 * (K + K.T)
    ev = np.linalg.eigvalsh(S)
    return dict(s_min=float(sv[-1]), s_max=float(sv[0]),
                cond=float(sv[0] / sv[-1]) if sv[-1] > 0 else float("inf"),
                s_min_rel=float(sv[-1] / scale) if scale else float("nan"),
                lam_min_sym_rel=float(ev[0] / scale) if scale else float("nan"))


def sample_tangent(setup):
    try:
        setup("FullGeneral")
        ops.test("NormUnbalance", 1.0e30, 1, 0)
        ops.algorithm("Linear")
        ops.integrator("LoadControl", 0.0)
        ops.analyze(1)
        th = tangent_health()
    except Exception as exc:                              # pragma: no cover
        log(f"    [diag 3] tangent sample failed: {exc}")
        th = None
    finally:
        setup("sparse")
    return th


def count_tension_gps(ncells, ngp):
    """Count Gauss points with mean stress (tension-positive convention) >= 0,
    plus any NaN encountered (ADR-94 B4's apex-NaN signature) and the ele/gp
    of the first NaN seen."""
    n_seen = n_tension = n_nan = 0
    first_nan = None
    i1_min, i1_max = float("inf"), float("-inf")
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
    return dict(n_seen=n_seen, n_tension=n_tension, n_nan=n_nan,
                first_nan=first_nan, i1_min=i1_min, i1_max=i1_max)


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
    log(f"=== leg {tag}: {spec['label']} on ASDPlasticMaterial3D/DruckerPrager_YF, "
        f"h0 = {args.h0} m, {len(nodes)} nodes / {3 * len(nodes)} DOF, "
        f"{len(cells)} elements")

    alpha, rho, eta_asd, xi_c_asd, e_el = build_model_asd(nodes, cells, sets, w, form)
    T_uw = SQ23 * HP.SY / rho
    T_asd = 3.0 * xi_c_asd / eta_asd
    log(f"    [mapping] alpha={alpha:.6f} rho={rho:.6f} -> ASD eta={eta_asd:.6f} "
        f"xi_c={xi_c_asd:.6f}; apex I1: UW T={T_uw:.6f}, ASD T={T_asd:.6f} "
        f"(should match to fp precision) ; E={e_el:.4f} nu={HP.NU}")
    assert abs(T_uw - T_asd) < 1e-9 * max(1.0, abs(T_uw)), "apex mapping mismatch"

    tol = 1.0e-5 * max(HP.Q0 * float(w.sum()), 1.0)
    HP.surcharge_step(nodes, sets, w, tol)

    # --- control 1: 1-D elastic patch test (material-agnostic elastic check) --
    k0 = HP.NU / (1.0 - HP.NU)
    e_zz = e_xx = 0.0
    for e in range(1, len(cells) + 1):
        st = np.asarray(ops.eleResponse(e, "stress"), dtype=float).reshape(-1, 6)
        e_zz = max(e_zz, float(np.abs(st[:, 2] / HP.Q0 + 1.0).max()))
        e_xx = max(e_xx, float(np.abs(st[:, 0] / (k0 * HP.Q0) + 1.0).max()))
    patch_err = max(e_zz, e_xx)
    log(f"    [control 1] 1-D elastic stress patch: max|szz/-q0-1| = {e_zz:.2e}, "
        f"max|sxx/-K0q0-1| = {e_xx:.2e}  -> {'PASS' if patch_err < 1e-6 else '*** FAIL ***'}")
    if patch_err >= 1e-6:
        raise SystemExit("patch test failed; every number downstream is void")

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
                ops.system("Pardiso")
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

    path = os.path.join(HERE, f"asd_{tag}.csv")
    fh = open(path, "w", newline="")
    wr = csv.writer(fh)
    wr.writerow(["s_m", "s_over_B", "q_kPa", "ds_mm", "relaxed", "wall_s"])

    tpath = os.path.join(HERE, f"asd_{tag}_tension.csv")
    tfh = open(tpath, "w", newline="")
    twr = csv.writer(tfh)
    twr.writerow(["s_over_B", "q_kPa", "n_seen", "n_tension", "n_nan",
                  "first_nan_ele", "first_nan_gp", "I1_min", "I1_max",
                  "s_min_rel", "cond", "lam_min_sym_rel"])

    rows, ds, good, nfail, nrelax, nsub = [], base, 0, 0, 0, 0
    mode, verdict = "TARGET", "reached the target settlement"
    cond_on = False
    stall_window = HP.STALL_WINDOW
    stall_advance = HP.STALL_ADVANCE
    grow_after = HP.GROW_AFTER
    nan_event = None

    t0 = time.time()
    while True:
        s_now = uz0 - ops.getTime()
        if s_now >= smax - 1e-12:
            break
        if time.time() - t0 > args.tmax:
            mode, verdict = "WALL", f"wall-clock cap at s/B = {s_now / HP.B_FOOT:.5f}"
            break
        ds = min(ds, smax - s_now)
        if not cond_on and ds < cond_at:
            cond_on = True
            log(f"    [diag 3] ds fell to {ds * 1e3:.6g} mm at s/B = "
                f"{s_now / HP.B_FOOT:.5f}: arming tangent + tension census")
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
                           f"s/B = {s_now / HP.B_FOOT:.5f}")
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
            th = sample_tangent(setup) if (args.cond and cond_on) else None
            tc = count_tension_gps(len(cells), GP_COUNTS[args.elem])
            twr.writerow([f"{s / HP.B_FOOT:.9g}", f"{q:.6g}", tc["n_seen"],
                          tc["n_tension"], tc["n_nan"],
                          tc["first_nan"][0] if tc["first_nan"] else "",
                          tc["first_nan"][1] if tc["first_nan"] else "",
                          f"{tc['i1_min']:.6g}", f"{tc['i1_max']:.6g}",
                          f"{th['s_min_rel']:.3e}" if th else "",
                          f"{th['cond']:.3e}" if th else "",
                          f"{th['lam_min_sym_rel']:.3e}" if th else ""])
            tfh.flush()
            log(f"    [tension] s/B {s / HP.B_FOOT:.5f}  q {q:8.2f}  "
                f"GP seen {tc['n_seen']}  tension(mean>=0) {tc['n_tension']}  "
                f"NaN {tc['n_nan']}  I1 [{tc['i1_min']:.3g},{tc['i1_max']:.3g}]"
                + (f"  |  sigma_min/scale {th['s_min_rel']:.3e} cond {th['cond']:.3e}"
                   if th else ""))
            if tc["n_nan"] and nan_event is None:
                nan_event = dict(s_over_B=s / HP.B_FOOT, q=q, ele=tc["first_nan"][0],
                                  gp=tc["first_nan"][1])
                log(f"    *** NaN/apex event *** first at s/B={nan_event['s_over_B']:.5f} "
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
        log(f"    last s/B = {s_arr[-1]:.5f}  q_end = {q_arr[-1]:.4f}  "
            f"q_max = {q_arr.max():.4f}")
    if nan_event:
        log(f"    SUMMARY: NaN/apex event recorded at s/B={nan_event['s_over_B']:.5f}, "
            f"ele={nan_event['ele']} gp={nan_event['gp']}")
    return dict(tag=tag, mode=mode, verdict=verdict, rows=rows, nan_event=nan_event,
                wall=wall)


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
    ap.add_argument("--cond", action="store_true")
    ap.add_argument("--suffix", default="")
    args = ap.parse_args()
    r = run(args)
    print()
    print(f"{r['tag']:>24} mode={r['mode']:<9} verdict={r['verdict']}")


if __name__ == "__main__":
    main()
