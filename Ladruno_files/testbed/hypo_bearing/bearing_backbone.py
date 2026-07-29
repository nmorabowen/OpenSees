"""ADR-79 bearing-backbone campaign — the geometric-nonlinearity ladder on a
square footing on saturated PDMY sand.

THE QUESTION (ADR 78 §1 / ADR 79 P3 follow-through): the TIMs footing model
hardens monotonically to ~9.3x the Vesic capacity with no limit point. ADR 78
blamed (in part) the missing large-deformation support. With the full ladder
now shipped — linear ⊂ corot (rotation) ⊂ hypo (genuine large strain + n(J)
porosity) ⊂ hypo -kozenyCarman (+ k(J)) — this campaign measures what each
rung actually does to the bearing backbone of a clean benchmark:

  * box 6B x 6B x 3B of saturated PDMY sand (medium-loose, phi=33), drained
    (k large, slow push), buoyant conditions via -body accel + hydrostatic p;
  * square rough footing B x B = the central 2x2-element patch, pushed
    DISPLACEMENT-controlled (the P3 scoping finding: force control diverges
    at the zero-confinement surface; -geom linear grinds and is excluded) to
    s = 15% B;
  * backbone R(s) per method, normalized by the Vesic surface capacity
    q_ult = 0.5*gamma'*B*Ngamma*s_gamma (phi=33: Nq=26.1, Ngamma=35.2,
    s_gamma=0.6).

Output: bearing_backbone_results.md next to this script (a table of R(s) and
q/q_Vesic per method + the verdict), plus raw CSV per leg.

Run (worktree, built dist/bin):
    py -3.12 Ladruno_files/testbed/hypo_bearing/bearing_backbone.py [quick]
`quick` runs a 1/2-resolution sanity pass first.
"""
import os
import sys
import time
import math
import csv

_WT = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))
_BIN = os.path.join(_WT, "dist", "bin")
os.environ.setdefault("LADRUNO_OPENSEES_QUIET", "1")
os.environ["PATH"] = _BIN + os.pathsep + os.environ.get("PATH", "")
if hasattr(os, "add_dll_directory"):
    os.add_dll_directory(_BIN)
sys.path.insert(0, _BIN)
import opensees as ops  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))

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

QUICK = len(sys.argv) > 1 and sys.argv[1] == "quick"
# GEOMETRY LESSON (first quick run): a 3 m box under a 2 m footing is an
# OEDOMETER, not a bearing problem — roller sides + 0.5 m clearance confine
# the mechanism and q/q_Vesic blows past 2500 with no meaning. The campaign
# box is 8 m (1.5B clearance each side, depth 4B); still modest — the
# CROSS-METHOD comparison is the point, not absolute capacity.
N = 4 if QUICK else 8            # elements per edge (box N x N x N)
L = 4.0 if QUICK else 8.0        # box edge (m); element size 1 m either way
SFRAC = 0.15                     # final penetration / footing width B
# ~2.5 mm per step — the increment size the P3 smoke proved convergent; a
# 25 mm first increment hard-fails every geometry method (measured).
NPUSH = int(round(SFRAC * 2.0 / 0.0025))
DT = 25.0

# footing = central 2x2-element patch => B = 2 m (both resolutions)
B_FOOT = 2.0
GAMMA_B = (RHO_SAT - RHO_W) * G9                    # buoyant unit weight kN/m^3
# Surface SURCHARGE q0 (whole top face, gravity stage): the standard PDMY
# footing-model regularization. Without it the free LATERAL DOFs of the
# zero-confinement surface nodes are nearly singular under stage-1 PDMY
# (G ~ p'^d): Newton limit-cycles at 25 m increments and even KrylovNewton
# stalls at ~1e-6 — measured on the N=4 quick pass; the N=3 P3 smoke sat
# just on the convergent side of the same marginality.
Q_SURCH = 10.0                                      # kPa
PHI = math.radians(PDMY["phi"])
NQ = math.exp(math.pi * math.tan(PHI)) * math.tan(math.pi / 4 + PHI / 2) ** 2
NGAMMA = 2.0 * (NQ + 1.0) * math.tan(PHI)           # Vesic
SGAMMA = 0.6                                        # square shape factor (Ngamma)
SQ = 1.0 + math.tan(PHI)                            # square shape factor (Nq)
# Vesic with the surcharge term: q_ult = q0*Nq*sq + 0.5*gamma'*B*Ngamma*sgamma
Q_VESIC = Q_SURCH * NQ * SQ + 0.5 * GAMMA_B * B_FOOT * NGAMMA * SGAMMA
A_FOOT = B_FOOT * B_FOOT


def log(*a):
    print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)


def run_leg(geom, kc=False):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 4)
    tags = {}
    t = 1
    for k in range(N + 1):
        for j in range(N + 1):
            for i in range(N + 1):
                ops.node(t, L * i / N, L * j / N, L * k / N)
                tags[(i, j, k)] = t
                t += 1
    p = PDMY
    ops.nDMaterial("PressureDependMultiYield", 1, 3, p["rho"], p["G"], p["B"],
                   p["phi"], p["gammaPeak"], p["refPress"], p["d"], p["PT"],
                   p["contract"], p["dil1"], p["dil2"], p["liq1"], p["liq2"],
                   p["liq3"], p["NYS"])
    conns = []
    for k in range(N):
        for j in range(N):
            for i in range(N):
                conns.append([tags[(i, j, k)], tags[(i + 1, j, k)],
                              tags[(i + 1, j + 1, k)], tags[(i, j + 1, k)],
                              tags[(i, j, k + 1)], tags[(i + 1, j, k + 1)],
                              tags[(i + 1, j + 1, k + 1)],
                              tags[(i, j + 1, k + 1)]])
    extra = ["-kozenyCarman"] if kc else []
    for e, conn in enumerate(conns, start=1):
        ops.element("LadrunoUP", e, *conn, 1,
                    "-Kf", KF, "-poro", PORO, "-rhoF", RHO_W,
                    "-perm", PERM, PERM, PERM, "-stab", "auto", 0.25,
                    "-body", 0.0, 0.0, -G9, "-geom", geom, *extra)
    for (i, j, k), tt in tags.items():
        fx = 1 if (i == 0 or i == N) else 0
        fy = 1 if (j == 0 or j == N) else 0
        fz = 1 if k == 0 else 0
        if k == 0:
            fx = fy = 1
        fp = 1 if k == N else 0
        ops.fix(tt, fx, fy, fz, fp)

    # surface surcharge (whole top face, tributary-weighted nodal loads) —
    # applied WITH gravity so the stage-1 surface tangent is regular.
    h = L / N
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for (i, j, k), tt in tags.items():
        if k != N:
            continue
        wx = 0.5 if (i == 0 or i == N) else 1.0
        wy = 0.5 if (j == 0 or j == N) else 1.0
        ops.load(tt, 0.0, 0.0, -Q_SURCH * wx * wy * h * h, 0.0)

    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("UmfPack")
    ops.test("NormDispIncr", 1e-7, 100, 0)
    ops.algorithm("KrylovNewton")
    ops.integrator("Newmark", 0.6, 0.3025)
    ops.analysis("Transient")
    for _ in range(10):
        assert ops.analyze(1, 500.0) == 0, f"{geom} elastic gravity failed"
    ops.updateMaterialStage("-material", 1, "-stage", 1)
    for _ in range(10):
        assert ops.analyze(1, 500.0) == 0, f"{geom} plastic settle failed"
    log(geom, ("(kc) " if kc else ""), "gravity done")

    # displacement-controlled footing push: the central 2x2-element patch =
    # the 3x3 node square centred on c (N even => a symmetric patch). Smooth
    # footing (only uz driven) — matches the P3 smoke.
    ops.loadConst("-time", 0.0)
    c = N // 2
    foot = [tags[(i, j, N)]
            for i in range(c - 1, c + 2) for j in range(c - 1, c + 2)]
    smax = SFRAC * B_FOOT
    times = [DT * (s + 1) for s in range(NPUSH)] + [DT * (NPUSH + 5)]
    vals = [-smax * (s + 1) / NPUSH for s in range(NPUSH)] + [-smax]
    ops.timeSeries("Path", 2, "-time", *times, "-values", *vals, "-useLast")
    ops.pattern("Plain", 2, 2)
    for tt in foot:
        ops.sp(tt, 3, 1.0)
    ops.wipeAnalysis()
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("UmfPack")
    ops.test("NormDispIncr", 1e-8, 60, 0)
    ops.algorithm("Newton")
    ops.integrator("Newmark", 0.6, 0.3025)
    ops.analysis("Transient")
    rows = []
    t0 = time.time()
    for s in range(NPUSH):
        ok = ops.analyze(1, DT)
        if ok != 0:
            ops.test("NormDispIncr", 1e-7, 100, 0)
            ops.algorithm("KrylovNewton")
            ok = 0
            for _ in range(10):
                if ops.analyze(1, DT / 10) != 0:
                    ok = -1
                    break
            ops.test("NormDispIncr", 1e-8, 60, 0)
            ops.algorithm("Newton")
        if ok != 0:
            log(geom, ("(kc)" if kc else ""),
                f"step {s + 1} HARD FAIL at s/B="
                f"{smax * (s + 1) / NPUSH / B_FOOT:.3f} — backbone TRUNCATED")
            break
        ops.reactions()
        R = -sum(ops.nodeReaction(tt, 3) for tt in foot)      # kN, compression +
        sset = smax * (s + 1) / NPUSH
        rows.append((sset, R, R / A_FOOT, R / A_FOOT / Q_VESIC))
        if (s + 1) % 5 == 0:
            log(geom, ("(kc)" if kc else ""), f"s/B={sset / B_FOOT:.3f} "
                f"q/qV={rows[-1][3]:.2f} [{time.time() - t0:.0f}s]")
    name = geom + ("_kc" if kc else "")
    with open(os.path.join(HERE, f"backbone_{name}{'_quick' if QUICK else ''}.csv"),
              "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["settlement_m", "R_kN", "q_kPa", "q_over_qVesic"])
        w.writerows(rows)
    return rows


def main():
    log(f"box {N}x{N}x{N} (L={L} m), footing B={B_FOOT} m, "
        f"q_Vesic={Q_VESIC:.1f} kPa, push to s/B={SFRAC}")
    legs = [("corot", False), ("hypo", False), ("hypo", True)]
    results = {}
    for geom, kc in legs:
        name = geom + ("_kc" if kc else "")
        log("=== leg:", name)
        results[name] = run_leg(geom, kc)
    # summary — truncation-honest: checkpoints are only reported for legs
    # whose backbone actually REACHED them, and the softening verdict only
    # runs on a completed backbone (a truncated leg reads 'TRUNCATED', never
    # 'softening' — the first quick run mislabeled a hard-fail freeze).
    log("=== summary (q/q_Vesic at s/B checkpoints; '--' = not reached)")
    for frac in (0.025, 0.05, 0.10, 0.15):
        line = f"s/B={frac:.3f}: "
        for name, rows in results.items():
            tgt = frac * B_FOOT
            if rows and rows[-1][0] >= tgt - 1e-9:
                best = min(rows, key=lambda r: abs(r[0] - tgt))
                line += f"{name}={best[3]:.2f}  "
            else:
                line += f"{name}=--  "
        log(line)
    for name, rows in results.items():
        if not rows:
            log(f"{name}: NO backbone (failed at step 1)")
            continue
        complete = rows[-1][0] >= SFRAC * B_FOOT - 1e-9
        q = [r[2] for r in rows]
        if not complete:
            log(f"{name}: TRUNCATED at s/B={rows[-1][0] / B_FOOT:.3f} "
                f"(q={q[-1]:.1f} kPa) — no softening verdict")
            continue
        tail = q[3 * len(q) // 4:]
        soft = tail[-1] < 0.995 * max(tail)
        log(f"{name}: q_final={q[-1]:.1f} kPa "
            f"({q[-1] / Q_VESIC:.2f}x Vesic), "
            f"{'SOFTENING (limit point) in the tail' if soft else 'still hardening'}")


if __name__ == "__main__":
    main()
