"""ADR-79 bearing-backbone campaign — the geometric-nonlinearity ladder on a
square footing on saturated PDMY sand.

THE QUESTION (ADR 78 §1 / ADR 79 P3 follow-through): the TIMs footing model
hardens monotonically to ~9.3x the Vesic capacity with no limit point. ADR 78
blamed (in part) the missing large-deformation support. With the full ladder
now shipped — linear ⊂ corot (rotation) ⊂ hypo (genuine large strain + n(J)
porosity) ⊂ hypo -kozenyCarman (+ k(J)) — this campaign measures what each
rung actually does to the bearing backbone of a clean benchmark:

  * a GRADED all-hex domain, 20 x 20 m plan x 8 m deep (4.5 B of side
    clearance, 4 B deep), 0.5 m hexes under and near the footing coarsening
    outward — built by build_mesh.py with apeGmsh, cached in bearing_mesh.npz.
    The clearance is what scoping finding 3 demands: every uniform-hex box
    small enough for a quick run was an OEDOMETER, not a bearing problem;
  * saturated PDMY sand (medium-loose, phi=33), drained (slow push), buoyant
    conditions via -body accel + a free-draining top;
  * 10 kPa surcharge over the top face BESIDE the footing — scoping finding 2:
    without it the zero-confinement surface DOFs are near-singular under
    stage-1 PDMY and every geometry method dies on the first push step. It is
    also the physical Vesic q0, so it enters the normalization;
  * square SMOOTH footing B x B = 2 x 2 m = the central 4 x 4-element patch,
    pushed DISPLACEMENT-controlled (scoping finding 1: force control diverges
    at the zero-confinement surface; -geom linear grinds and is excluded) to
    s = 15% B in ~2.5 mm steps (finding 4);
  * backbone R(s) per method, normalized by Vesic
    q_ult = q0*Nq*sq + 0.5*gamma'*B*Ngamma*s_gamma.

Output: raw CSV per leg (written incrementally — a killed leg keeps its data)
plus a truncation-honest summary. A truncated leg NEVER reads as "softening".

Run (worktree, built dist/bin):
    py -3.12 bearing_backbone.py [corot|hypo|hypo_kc|corot_bbar|all]

`all` runs the three original ladder legs. `corot_bbar` is the separate
LOCKING leg (see LEGS) — it is paired against `corot`, not part of the ladder.
"""
import os
import sys
import time
import math
import csv

import numpy as np

_WT = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))
_BIN = os.path.join(_WT, "dist", "bin")
os.environ.setdefault("LADRUNO_OPENSEES_QUIET", "1")
os.environ["PATH"] = _BIN + os.pathsep + os.environ.get("PATH", "")
if hasattr(os, "add_dll_directory"):
    os.add_dll_directory(_BIN)
sys.path.insert(0, _BIN)
import opensees as ops  # noqa: E402

# The INSTALLED Ladruno wires `C:\Program Files\Ladruno\...` into every venv it
# touches through a site .pth (`ladruno_opensees.pth` -> `_ladruno_opensees_boot`)
# that runs `import opensees` at INTERPRETER STARTUP. In such a venv `opensees`
# is already in sys.modules before line 1 of this script, so the sys.path
# insert above is silently ignored and the campaign would run on the installed
# build — which predates ADR 78/79 and refuses `-geom corot|hypo` (measured:
# "only 'linear' is accepted"). Use a clean interpreter (the base
# C:\...\Programs\Python\Python312\python.exe has no such .pth) and let this
# guard prove which engine is actually loaded.
if os.path.normcase(os.path.dirname(os.path.abspath(ops.__file__))) != \
        os.path.normcase(_BIN):
    raise SystemExit(
        f"WRONG ENGINE: imported {ops.__file__}\n"
        f"  expected the worktree build in {_BIN}\n"
        f"  a site .pth pre-imported another opensees — run this with an "
        f"interpreter that has no Ladruno .pth (base Python 3.12).")

HERE = os.path.dirname(os.path.abspath(__file__))
MESH = os.path.join(HERE, "bearing_mesh.npz")

G9 = 9.81
RHO_SAT, RHO_W = 2.0, 1.0        # Mg/m^3
KF, PORO = 2.2e6, 0.4            # kPa
# k = 1e-4 is already fully drained at these steps (cv = k*M ~ 20 m^2/s,
# t_push = 3000 s >> t99 of the box) AND is the value the P3 smoke proved
# convergent — k = 1e-2 made Newton fail from the very first 2.5 mm push
# step (measured; the p-row conditioning, not the physics, is the issue).
PERM = 1.0e-4
PDMY = dict(rho=RHO_SAT, G=5.5e4, B=1.5e5, phi=33.0, gammaPeak=0.1,
            refPress=80.0, d=0.5, PT=27.0, contract=0.05, dil1=0.0, dil2=0.0,
            liq1=0.0, liq2=0.0, liq3=0.0, NYS=20)

B_FOOT = 2.0
SFRAC = 0.15                     # final penetration / footing width B
Q_SURCH = 10.0                                      # kPa

# --- push increment (scoping finding 4, RE-MEASURED on the graded mesh) -----
# P3's 2.5 mm increment was calibrated on 1 m hexes. On this mesh (0.5 m hexes
# at the surface, 4.5 B of clearance) it does not converge for ANY geometry
# method — including `-geom linear`, so it is the model/mesh, not the geometry
# lane — and it is NOT a tolerance artifact: Newton diverges (‖Δu‖ ~ 0.2 m) and
# KrylovNewton stalls at ‖Δu‖ ~ 3e-5 for every test/tolerance tried down to
# 1e-4, including relative and energy tests. Measured threshold: 0.05 mm
# converges 12/12 at ~2.4 s/step, 0.10 mm and above fail. `dt` is irrelevant
# (2.5 / 25 / 250 s all fail identically at 2.5 mm) — the limit is the
# DISPLACEMENT increment.
#
# So the push is a linear displacement ramp in pseudo-time and the increment is
# carried by `dt`, which lets it ADAPT: halve on failure, grow back after a run
# of successes. T_PUSH is chosen so that even the smallest admissible increment
# still spans a fully DRAINED interval — c_v ≈ k·M/γ_w ≈ 2.3 m²/s, so
# T_v = c_v·dt/B² ≫ 1 needs dt ≳ 5 s, which DS_MIN maps to exactly.
DS_BASE, DS_MIN, DS_MAX = 5.0e-5, 6.25e-6, 4.0e-4   # m per step
GROW_AFTER = 5                   # consecutive good steps before doubling ds
T_PUSH = (SFRAC * B_FOOT / DS_MIN) * 5.0            # s, => dt(DS_MIN) = 5 s
# smoke/sizing knob only: cap the step count without touching the increment.
MAXSTEP = int(os.environ.get("ADR79_MAXSTEP", 200000))

GAMMA_B = (RHO_SAT - RHO_W) * G9                    # buoyant unit weight kN/m^3
PHI = math.radians(PDMY["phi"])
NQ = math.exp(math.pi * math.tan(PHI)) * math.tan(math.pi / 4 + PHI / 2) ** 2
NGAMMA = 2.0 * (NQ + 1.0) * math.tan(PHI)           # Vesic
SGAMMA = 0.6                                        # square shape factor (Ngamma)
SQ = 1.0 + math.tan(PHI)                            # square shape factor (Nq)
Q_VESIC = Q_SURCH * NQ * SQ + 0.5 * GAMMA_B * B_FOOT * NGAMMA * SGAMMA
A_FOOT = B_FOOT * B_FOOT

# leg -> (geom, kozenyCarman, formulation). `corot_bbar` is the LOCKING leg:
# every rung above runs `-formulation std`, so the whole geometry ladder was
# measured on a volumetrically LOCKED discretisation and the reported
# corot->hypo deltas are the geometric content of a locked mesh. `-geom hypo`
# refuses bbar (OPS_LadrunoUP.cpp, ADR 79 P2: "bbar in rate form is reserved"),
# but corot composes with it unchanged (ADR 78 §3.1 — bbar is a core-frame
# choice, and Q is barred through the same dNuBar_), so corot is where the
# locking lever and its cross-term with geometry ARE measurable today. Paired
# against the `corot` leg, same mesh / material / push.
# --- undrained pair (see UNDRAINED note below) ------------------------------
PERM_UND, TSCALE_UND = 1.0e-9, 100.0


def _leg(geom, kc=False, form="std", perm=PERM, tscale=1.0, grav_dt=500.0):
    return dict(geom=geom, kc=kc, form=form, perm=perm, tscale=tscale,
                grav_dt=grav_dt)


# 20 gravity steps must span T_v >> 1 so the initial state is CONSOLIDATED.
# At k=1e-9, B^2/c_v = 1.8e5 s, so 500 s/step (the drained legs' value) is
# T_v = 0.06 — undrained, and it breaks the plastic settle. 5e4 s/step gives
# T_v = 5.7. See gravity() and scoping finding 8.
GRAV_DT_UND = 5.0e4


LEGS = {"corot": _leg("corot"), "hypo": _leg("hypo"),
        "hypo_kc": _leg("hypo", kc=True),
        "corot_std": _leg("corot"),
        "corot_bbar": _leg("corot", form="bbar"),
        "corot_std_und": _leg("corot", perm=PERM_UND, tscale=TSCALE_UND,
                              grav_dt=GRAV_DT_UND),
        "corot_bbar_und": _leg("corot", form="bbar", perm=PERM_UND,
                               tscale=TSCALE_UND, grav_dt=GRAV_DT_UND)}
LADDER = ["corot", "hypo", "hypo_kc"]
# `corot_std` is `corot` re-run on WHATEVER build is in dist/bin right now, so
# the locking ratio is engine-matched. It exists because the committed
# backbone_corot.csv came from a different build than a fresh worktree carries,
# and a ratio that mixes formulation with engine is not a measurement. Run the
# pair (corot_std, corot_bbar) together; it also keeps backbone_corot.csv — the
# committed §8 evidence — from being overwritten.
PAIR = ["corot_std", "corot_bbar"]
PAIR_UND = ["corot_std_und", "corot_bbar_und"]
# UNDRAINED pair — why these two numbers and not others. bbar's lever is
# NEAR-INCOMPRESSIBILITY (the repo's own gate: bbar/std = 2.31 at nu=0.499 but
# 1.008 at nu=0.3), and undrained loading is what makes the pore fluid impose
# it. Undrained-ness is the consolidation time factor T_v = c_v*t/B^2 << 1,
# c_v = k*M/gamma_w, M = B + 4G/3 ~ 2.23e5 kPa => c_v = 2.28e4*k.
#   * PERMEABILITY ALONE IS NOT ENOUGH. The drained legs sit at T_v ~ 8.5e4;
#     even clay-grade k=1e-8 only reaches T_v=8.5 (still drained), because the
#     adaptive push spends ~1.5e5 s of pseudo-time getting to s/B=0.06.
#   * THE RATE HAS A FLOOR. Vs = sqrt(G/rho) = 166 m/s, so a shear wave crosses
#     B in 0.0121 s. T_PUSH/1000 would put dt at the SMALLEST increment at
#     0.005 s = 0.4x transit — that leg would measure inertial stiffening, not
#     undrained locking. T_PUSH/100 keeps dt in [0.05, 3.2] s = 4.1x..265x
#     transit, quasi-static throughout.
# So: tscale=100 (the most rate the wave floor allows) and k=1e-9 supplies the
# rest => T_v = 8.5e-3. Both knobs are needed; neither alone gets there.
# The runner logs T_v and records pore pressure under the footing so the
# undrained claim is MEASURED, not just designed.


def log(*a):
    print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)


def load_mesh():
    if not os.path.exists(MESH):
        raise SystemExit(f"{MESH} missing — run build_mesh.py first "
                         "(needs the apeGmsh env, see its docstring)")
    z = np.load(MESH)
    return (z["nodes"], z["hexes"], z["tributary"],
            {k[4:]: z[k] for k in z.files if k.startswith("set_")})


def pick_system():
    """PARDISO if this build has it (LadrunoUP's tangent is UNSYMMETRIC, so no
    -matrixType), else UmfPack. Both are exact direct solves — this is a
    wall-clock choice, not a physics one."""
    try:
        ops.system("Pardiso")
        return "Pardiso"
    except Exception:
        ops.system("UmfPack")
        return "UmfPack"


def build(cfg, nodes, hexes, trib, sets):
    geom, kc, form, perm = cfg["geom"], cfg["kc"], cfg["form"], cfg["perm"]
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 4)
    for i, (x, y, zz) in enumerate(nodes, start=1):
        ops.node(i, float(x), float(y), float(zz))
    p = PDMY
    ops.nDMaterial("PressureDependMultiYield", 1, 3, p["rho"], p["G"], p["B"],
                   p["phi"], p["gammaPeak"], p["refPress"], p["d"], p["PT"],
                   p["contract"], p["dil1"], p["dil2"], p["liq1"], p["liq2"],
                   p["liq3"], p["NYS"])
    extra = ["-kozenyCarman"] if kc else []
    # std is the parser default; pass nothing for it so the three original legs
    # stay byte-identical to the runs behind the committed CSVs.
    if form != "std":
        extra += ["-formulation", form]
    for e, conn in enumerate(hexes, start=1):
        ops.element("LadrunoUP", e, *[int(c) + 1 for c in conn], 1,
                    "-Kf", KF, "-poro", PORO, "-rhoF", RHO_W,
                    "-perm", perm, perm, perm, "-stab", "auto", 0.25,
                    "-body", 0.0, 0.0, -G9, "-geom", geom, *extra)

    # roller sides / fixed base / free-draining top (p = 0 at z = 0)
    fix = {}
    for n in sets["xface"]:
        fix.setdefault(int(n), [0, 0, 0, 0])[0] = 1
    for n in sets["yface"]:
        fix.setdefault(int(n), [0, 0, 0, 0])[1] = 1
    for n in sets["bottom"]:
        f = fix.setdefault(int(n), [0, 0, 0, 0])
        f[0] = f[1] = f[2] = 1
    for n in sets["top"]:
        fix.setdefault(int(n), [0, 0, 0, 0])[3] = 1
    for n, f in fix.items():
        ops.fix(n + 1, *f)

    # surcharge over the WHOLE top face (tributary-weighted nodal loads),
    # applied WITH gravity so the stage-1 surface tangent is regular — scoping
    # finding 2, and it must cover the footing patch too (see build_mesh.py).
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for n in sets["top"]:
        a = float(trib[n])
        if a > 0.0:
            ops.load(int(n) + 1, 0.0, 0.0, -Q_SURCH * a, 0.0)


def gravity(tag, grav_dt):
    """Establish the initial effective-stress state. `grav_dt` MUST be long
    enough that the gravity + surcharge state has CONSOLIDATED, whatever this
    leg's permeability is — the initial state is drained even when the push is
    not (scoping finding 8). At k = 1e-9 the stock 500 s leaves gravity
    undrained: pore pressure carries the overburden, effective stress stays
    ~0, and `updateMaterialStage 1` then hands PDMY a near-zero-confinement
    tangent that will not settle ("plastic settle failed") — the same
    near-singular surface DOFs as finding 2, reached through drainage instead
    of a missing surcharge."""
    ops.constraints("Transformation")
    ops.numberer("RCM")
    sysname = pick_system()
    ops.test("NormDispIncr", 1e-7, 100, 0)
    ops.algorithm("KrylovNewton")
    ops.integrator("Newmark", 0.6, 0.3025)
    ops.analysis("Transient")
    for _ in range(10):
        assert ops.analyze(1, grav_dt) == 0, f"{tag} elastic gravity failed"
    ops.updateMaterialStage("-material", 1, "-stage", 1)
    for _ in range(10):
        assert ops.analyze(1, grav_dt) == 0, f"{tag} plastic settle failed"
    return sysname


def run_leg(name):
    cfg = LEGS[name]
    geom, kc, form = cfg["geom"], cfg["kc"], cfg["form"]
    nodes, hexes, trib, sets = load_mesh()
    t_start = time.time()
    build(cfg, nodes, hexes, trib, sets)
    log(name, f"model built: {len(nodes)} nodes, {len(hexes)} hex "
        f"[-geom {geom}, -formulation {form}{', -kozenyCarman' if kc else ''}, "
        f"k={cfg['perm']:.1e}]")
    # the drainage regime this leg will actually be in, stated up front:
    # T_v = c_v*t/B^2 over the whole push, c_v = k*M/gamma_w, M = B + 4G/3.
    t_push = T_PUSH / cfg["tscale"]
    c_v = cfg["perm"] * (PDMY["B"] + 4.0 * PDMY["G"] / 3.0) / (RHO_W * G9)
    t_v = c_v * t_push / B_FOOT ** 2
    log(name, f"drainage: c_v={c_v:.3g} m2/s, t_push={t_push:.3g} s => "
        f"T_v={t_v:.3g} ({'UNDRAINED' if t_v < 0.05 else 'partially drained' if t_v < 1 else 'DRAINED'})")
    sysname = gravity(name, cfg["grav_dt"])
    log(name, f"gravity done [{sysname}, {time.time() - t_start:.0f}s]")

    foot = [int(n) + 1 for n in sets["footing"]]
    # nodeReaction at a constrained DOF is (internal - applied), so the raw
    # sum understates the load through the footing by exactly the surcharge
    # sitting on the footing patch. Correct it in closed form (= q0 * B^2).
    r_corr = Q_SURCH * float(trib[sets["footing"]].sum())
    # settlement datum: the footing nodes have already settled under gravity +
    # surcharge, and an `sp` value is an ABSOLUTE nodal displacement. Record
    # both the commanded and the MEASURED settlement so the abscissa is a fact,
    # not an assumption about sp semantics.
    uz0 = ops.nodeDisp(foot[0], 3)
    # Undrained-ness must be MEASURED, not inferred from T_v: pick the interior
    # node nearest (footing centre, B/2 below the surface) and record its pore
    # pressure. The footing nodes themselves are useless for this — they sit on
    # the free-draining top face, where dof 4 is fixed to p = 0. Drained => p
    # stays ~0; undrained => p tracks the applied pressure.
    fxy = nodes[[n - 1 for n in foot], :]
    z_top = float(nodes[:, 2].max())
    probe_xyz = np.array([fxy[:, 0].mean(), fxy[:, 1].mean(),
                          z_top - 0.5 * B_FOOT])
    p_node = int(np.argmin(((nodes - probe_xyz) ** 2).sum(axis=1))) + 1
    # Datum matters: the probe is 1 m down, where HYDROSTATIC p is already
    # 9.81 kPa, so absolute p is ~10 kPa before the push even starts and a raw
    # p/q ratio reads as undrained no matter the permeability. Record EXCESS
    # pore pressure over the consolidated gravity state — that is the part the
    # push generates, and the part that is ~0 in a drained leg.
    p0 = ops.nodeDisp(p_node, 4)
    log(name, f"pore-pressure probe: node {p_node} at "
        f"({nodes[p_node - 1, 0]:.2f}, {nodes[p_node - 1, 1]:.2f}, "
        f"{nodes[p_node - 1, 2]:.2f}), datum p0={p0:.3f} kPa "
        f"(hydrostatic would be {RHO_W * G9 * (z_top - nodes[p_node - 1, 2]):.3f})")
    ops.loadConst("-time", 0.0)
    smax = SFRAC * B_FOOT
    rate = smax / t_push                       # m of settlement per second
    # a 2-point linear ramp: settlement is exactly rate * pseudo-time, so the
    # per-step increment is carried by dt and can be adapted freely.
    ops.timeSeries("Path", 2, "-time", 0.0, t_push,
                   "-values", uz0, uz0 - smax, "-useLast")
    ops.pattern("Plain", 2, 2)
    for tt in foot:
        ops.sp(tt, 3, 1.0)

    ops.wipeAnalysis()
    ops.constraints("Transformation")
    ops.numberer("RCM")
    pick_system()
    # KrylovNewton is the PRIMARY here, not a fallback: on this mesh plain
    # Newton diverges at every increment tested down to 0.25 mm, while
    # KrylovNewton converges cleanly at 0.05 mm. The adaptive halving below
    # replaces P3's fixed dt/10 ladder.
    ops.test("NormDispIncr", 1e-7, 100, 0)
    ops.algorithm("KrylovNewton")
    ops.integrator("Newmark", 0.6, 0.3025)
    ops.analysis("Transient")

    # CSV-collision guard. These legs run for hours, in parallel, writing
    # incrementally to backbone_<leg>.csv. A short smoke of the SAME leg name
    # reopens that path with "w", truncates it under the live process, and the
    # live process then keeps writing at its old file offset — leaving the
    # smoke's rows, a block of NUL padding, and the tail of the real run, with
    # the real run's whole early backbone gone. (Measured, the hard way.) So:
    # a capped run is a smoke and gets its own filename, and any leg refuses to
    # clobber a file another process is actively appending to.
    smoke = MAXSTEP < 200000
    path = os.path.join(HERE,
                        f"backbone_{name}{'__smoke' if smoke else ''}.csv")
    if os.path.exists(path) and time.time() - os.path.getmtime(path) < 180 \
            and os.environ.get("ADR79_FORCE") != "1":
        raise SystemExit(
            f"REFUSING to write {os.path.basename(path)}: modified "
            f"{time.time() - os.path.getmtime(path):.0f}s ago, so another "
            f"process is probably still appending to it. Wait for it, or set "
            f"ADR79_FORCE=1 if you are certain it is dead.")
    f = open(path, "w", newline="")
    w = csv.writer(f)
    w.writerow(["s_m", "s_meas_m", "R_kN", "q_kPa", "q_over_qVesic",
                "ds_mm", "wall_s", "p_exc_kPa"])
    rows = []
    truncated_at = None
    ds, good, nfail = DS_BASE, 0, 0
    t0 = time.time()
    nstep = 0
    while nstep < MAXSTEP:
        s_now = ops.getTime() * rate
        if s_now >= smax - 1e-12:
            break
        ds = min(ds, smax - s_now)
        if ops.analyze(1, ds / rate) != 0:
            nfail += 1
            good = 0
            ds *= 0.5
            if ds < DS_MIN:
                truncated_at = s_now
                log(name, f"NO CONVERGENCE at s/B={s_now / B_FOOT:.4f} even at "
                    f"ds={ds * 1000:.4f} mm < DS_MIN — backbone TRUNCATED")
                break
            continue
        nstep += 1
        good += 1
        if good >= GROW_AFTER and ds < DS_MAX:
            ds, good = min(2.0 * ds, DS_MAX), 0
        ops.reactions()
        R = -sum(ops.nodeReaction(tt, 3) for tt in foot) + r_corr   # kN, compression +
        s = ops.getTime() * rate
        smeas = uz0 - ops.nodeDisp(foot[0], 3)
        row = (s, smeas, R, R / A_FOOT, R / A_FOOT / Q_VESIC, ds * 1000,
               time.time() - t0, ops.nodeDisp(p_node, 4) - p0)
        rows.append(row)
        w.writerow([f"{v:.9g}" for v in row])
        f.flush()
        if nstep % 100 == 0:
            log(name, f"s/B={s / B_FOOT:.4f} q={row[3]:.1f} kPa "
                f"q/qV={row[4]:.2f} dp={row[7]:.1f} kPa ds={ds * 1000:.3f}mm "
                f"[{nstep} steps, {nfail} retries, {row[6]:.0f}s]")
    f.close()
    if truncated_at is None and nstep >= MAXSTEP:
        truncated_at = rows[-1][0] if rows else 0.0
        log(name, f"hit MAXSTEP={MAXSTEP} — backbone TRUNCATED")
    log(name, f"leg done: {nstep} steps to s/B="
        f"{(rows[-1][0] / B_FOOT) if rows else 0:.4f}, {nfail} retries, "
        f"{time.time() - t_start:.0f}s -> {os.path.basename(path)}")
    return rows, truncated_at


def summarize(results):
    log("=== summary (q/q_Vesic at s/B checkpoints; '--' = not reached)")
    for frac in (0.01, 0.025, 0.05, 0.10, 0.15):
        line = f"s/B={frac:.3f}: "
        for name, (rows, _) in results.items():
            tgt = frac * B_FOOT
            if rows and rows[-1][0] >= tgt - 1e-9:
                best = min(rows, key=lambda r: abs(r[0] - tgt))
                line += f"{name}={best[4]:.2f}  "
            else:
                line += f"{name}=--  "
        log(line)
    for name, (rows, trunc) in results.items():
        if not rows:
            log(f"{name}: NO backbone (failed at step 1)")
            continue
        q = [r[3] for r in rows]
        complete = trunc is None and rows[-1][0] >= SFRAC * B_FOOT - 1e-9
        if not complete:
            log(f"{name}: TRUNCATED at s/B={rows[-1][0] / B_FOOT:.4f} "
                f"(q={q[-1]:.1f} kPa) — no softening verdict")
            continue
        tail = q[3 * len(q) // 4:]
        soft = tail[-1] < 0.995 * max(tail)
        log(f"{name}: q_final={q[-1]:.1f} kPa ({q[-1] / Q_VESIC:.2f}x Vesic), "
            f"{'SOFTENING (limit point) in the tail' if soft else 'still hardening'}")


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    # "all" = the geometry ladder only; corot_bbar is a paired locking probe and
    # must be asked for by name so `all` keeps reproducing the §8 campaign.
    names = (LADDER if which == "all" else PAIR if which == "pair" else
             PAIR_UND if which == "pair_und" else [which])
    for n in names:
        if n not in LEGS:
            raise SystemExit(f"unknown leg {n!r}; pick from {list(LEGS)} or 'all'")
    log(f"footing B={B_FOOT} m, q_Vesic={Q_VESIC:.1f} kPa "
        f"(q0*Nq*sq={Q_SURCH * NQ * SQ:.1f} + gamma term="
        f"{0.5 * GAMMA_B * B_FOOT * NGAMMA * SGAMMA:.1f}), "
        f"push to s/B={SFRAC} ({SFRAC * B_FOOT * 1000:.0f} mm), adaptive "
        f"increment base={DS_BASE * 1000:.3f} mm "
        f"in [{DS_MIN * 1000:.4f}, {DS_MAX * 1000:.3f}] mm")
    results = {}
    for n in names:
        log("=== leg:", n)
        results[n] = run_leg(n)
    summarize(results)


if __name__ == "__main__":
    main()
