"""ADR-77 T3 — implicit scheme shootout on COST PER SIMULATED SECOND at MATCHED ACCURACY.

This is the ADR-49a Rank-3 "assess TRBDF2 vs modern implicit Bathe before
committing" finally executed, and it feeds gate G3.

WHY THIS IS A CONVERGENCE STUDY AND NOT A RACE. "Cost per simulated second" is
meaningless on its own: any scheme can be made cheap by taking a bigger dt and
being wrong. The deliverable is cost AT MATCHED ACCURACY, so every arm is run
over a dt ladder to a FIXED end time, its error is measured against a common
fine-dt reference, and the arms are compared at equal error -- not at equal dt.

MATCHED-RHO-INFINITY IS THE OTHER FAIRNESS AXIS. Comparing "HHT" to
"GeneralizedAlpha" without saying how much numerical damping each carries
compares two arbitrary points, not two schemes. Conventions, derived from the
constructors in source (2026-07-26) rather than assumed:

  HHT(a_ops):    beta=(2-a)^2/4, gamma=1.5-a  => a_lit = a_ops-1 in [-1/3,0]
                 rho_inf = a_ops/(2-a_ops)
                 a_ops=1.0 -> rho=1.000 (trapezoidal, no damping)
                 a_ops=0.9 -> rho=0.818
                 a_ops=0.8 -> rho=0.667
  GeneralizedAlpha(aM,aF): beta=(1+aM-aF)^2/4, gamma=0.5+aM-aF
                 aM=(2-rho)/(1+rho), aF=1/(1+rho)
                 rho=1 -> (0.5,0.5) = trapezoidal  [verified against the ctor]

So HHT(0.9) and GeneralizedAlpha at rho=0.818 are the SAME damping level and
their difference is the scheme, which is the comparison worth making.

ALGORITHM IS HELD FIXED AT KrylovNewton -- T1's winner. Using Newton would
penalise every arm on the plastic deck where it diverges, and the point here is
to vary the INTEGRATOR, not to re-run T1.

DECK: SMOOTH RAMP LOADING, not T0/T1's step load, and the reason is a
methodological one measured on 2026-07-26. T0/T1 release a Constant load at
t=0. That is a SHOCK: the element wave-transit time is L/sqrt(E/rho) =
100/5.05e6 = 2.0e-5 s, so T1's dt=2e-4 is TEN transit times per step. A first
attempt at this stage inherited that deck and produced (i) a reference that
failed its own self-convergence check at 7.9e-3 and (ii) ladder errors of
24-62%. Both are correct readings of a study that was measuring how each scheme
survives a discontinuity, not how accurate it is. A convergence study needs a
SMOOTH solution or every scheme is order-limited by the forcing and the
comparison is meaningless. So T3 ramps the load linearly (timeSeries Linear)
to the same x4 plastic level by T_end, and sizes dt against the first-mode
period (~4*H/c) rather than against T1's step count.

LOAD LEVEL x1, NOT T1's x4, AND THAT TOO WAS MEASURED. Under the ramp the x4
deck yields 86.8% of sampled Gauss points by T_end (uz = -2.6 mm) -- it is
approaching a COLLAPSE MECHANISM, and a first pass there showed the symptoms:
~170 Newton iterations/step, and error ratios of ~1.84 per dt doubling (i.e.
FIRST-order convergence from schemes that are second-order by construction).
Both are properties of a near-singular tangent and of return-map error at yield
transitions, not of the integrators being compared. x1 under the ramp yields
12.5% of sampled Gauss points: solidly, genuinely plastic without being
degenerate. That is the regime where an integrator comparison means something.
x0.5 was also checked and is fully elastic (0% yielded), so it would have been
the mirror-image mistake.

G3's wording asks for "contact/sharp-nonlinearity"; this covers the PLASTICITY
axis only. Contact is a different deck and is NOT covered here -- G3 can
therefore only be decided on the plasticity axis by this run, which the results
doc must say.

Run:
  OPS_PYD=<worktree>\\dist\\bin MKL_NUM_THREADS=4 python3.12 -S t3_scheme_shootout.py
"""
from __future__ import annotations

import json
import os
import statistics
import sys
import time
from pathlib import Path

sys.stdout.reconfigure(line_buffering=True)
os.environ.setdefault("MKL_NUM_THREADS", "4")
os.environ.setdefault("OMP_NUM_THREADS", os.environ["MKL_NUM_THREADS"])
os.environ.setdefault("PMI_RANK", "0")
os.environ.setdefault("LADRUNO_OPENSEES_QUIET", "1")

DIST = os.environ.get(
    "OPS_PYD",
    r"C:\Users\nmb\Documents\Github\OpenSees"
    r"\.claude\worktrees\pardisio-profiling-0a03b1\dist\bin",
)
sys.path.insert(0, DIST)
if os.path.isdir(DIST):
    os.add_dll_directory(DIST)
import opensees as ops  # noqa: E402

HERE = Path(__file__).resolve().parent
N = int(os.environ.get("OPS_T3_N", "10"))
LOAD = float(os.environ.get("OPS_T3_LOAD", "1.0"))
SYSTEM = os.environ.get("OPS_T3_SYSTEM", "Pardiso")
MAXITER = int(os.environ.get("OPS_T3_MAXITER", "1000"))
REPEATS = int(os.environ.get("OPS_T3_REPEATS", "2"))
TAG = os.environ.get("OPS_T3_TAG", "")

E, nu = 200000.0, 0.3
K = E / (3.0 * (1.0 - 2.0 * nu))
G = E / (2.0 * (1.0 + nu))
s0, Qinf, b, Hiso = 250.0, 0.0, 0.0, 2000.0
RHO = 7.85e-9
L = 100.0
TOL = 1.0e-7

# Sized against the physics, not against T1: wave speed c = sqrt(E/rho) =
# 5.05e6 mm/s; element transit L/c = 2.0e-5 s; first-mode period ~ 4*H/c which
# for the n=10 block (H=1000mm) is ~7.9e-4 s. T_end covers ~1 period and the
# ladder spans T/79 .. T/10.
T_END = float(os.environ.get("OPS_T3_TEND", "8.0e-4"))
DT_BASE = float(os.environ.get("OPS_T3_DTBASE", "1.0e-5"))
DT_LADDER = [DT_BASE * m for m in (1, 2, 4, 8)]
DT_REF = DT_BASE / 8.0
DT_REF2 = DT_BASE / 16.0      # convergence check ON THE REFERENCE ITSELF


def rho_to_genalpha(rho):
    return (2.0 - rho) / (1.0 + rho), 1.0 / (1.0 + rho)


# (key, label, rho_inf_or_None, integrator argv)
SCHEMES = [
    ("newmark_trap", "Newmark(0.5,0.25) trapezoidal", 1.000, ["Newmark", 0.5, 0.25]),
    ("newmark_damp", "Newmark(0.6,0.3025) damped",    None,  ["Newmark", 0.6, 0.3025]),
    ("hht_090",      "HHT(0.90)",                     0.818, ["HHT", 0.9]),
    ("hht_080",      "HHT(0.80)",                     0.667, ["HHT", 0.8]),
    ("ga_0818",      "GeneralizedAlpha(rho=0.818)",   0.818,
     ["GeneralizedAlpha", *rho_to_genalpha(0.818)]),
    ("ga_0667",      "GeneralizedAlpha(rho=0.667)",   0.667,
     ["GeneralizedAlpha", *rho_to_genalpha(0.667)]),
    ("trbdf2",       "TRBDF2 (composite Trap+BDF2)",  None,  ["TRBDF2"]),
    ("backeuler",    "BackwardEuler (1st order)",     0.0,   ["BackwardEuler"]),
]


def build(n, integrator_argv):
    nx = ny = nz = n

    def nid(i, j, k):
        return 1 + i + (nx + 1) * (j + (ny + 1) * k)

    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("LadrunoJ2", 1, K, G, "-iso", "voce", s0, Qinf, b, Hiso,
                   "-rho", RHO)
    for k in range(nz + 1):
        for j in range(ny + 1):
            for i in range(nx + 1):
                ops.node(nid(i, j, k), i * L, j * L, k * L)
    for j in range(ny + 1):
        for i in range(nx + 1):
            ops.fix(nid(i, j, 0), 1, 1, 1)
    e = 1
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                ops.element(
                    "LadrunoBrick", e,
                    nid(i, j, k), nid(i+1, j, k), nid(i+1, j+1, k), nid(i, j+1, k),
                    nid(i, j, k+1), nid(i+1, j, k+1), nid(i+1, j+1, k+1), nid(i, j+1, k+1),
                    1)
                e += 1
    # Linear series => factor = t, so the load RAMPS smoothly from zero and
    # reaches the x4 plastic level at T_end. Scaled by 1/T_END for that.
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    P = 4.0e5 * LOAD / T_END
    for j in range(ny + 1):
        for i in range(nx + 1):
            ops.load(nid(i, j, nz), P, 0.0, -P)

    ops.constraints("Plain")
    ops.numberer("RCM")
    ops.system(SYSTEM)
    ops.test("NormDispIncr", TOL, MAXITER)
    ops.algorithm("KrylovNewton")           # T1's winner, held fixed
    ops.integrator(*integrator_argv)
    ops.analysis("Transient")
    return nid(nx // 2, ny // 2, nz)


SAMPLE_DT = float(os.environ.get("OPS_T3_SAMPLEDT", "8.0e-5"))


def run(dt, integrator_argv):
    """Integrate to T_END, sampling uz at SAMPLE_DT.

    Returns (wall, history[list], iters_total, nsteps). The history -- not the
    endpoint -- is the accuracy signal; see the module docstring.
    """
    tip = build(N, integrator_argv)
    nsteps = int(round(T_END / dt))
    if nsteps < 1:
        raise RuntimeError(f"dt={dt:g} exceeds T_end={T_END:g}: 0 steps")
    per_sample = int(round(SAMPLE_DT / dt))
    if per_sample < 1 or abs(per_sample * dt - SAMPLE_DT) > 1e-15:
        raise RuntimeError(f"dt={dt:g} does not divide sample interval {SAMPLE_DT:g}")
    nsamples = int(round(T_END / SAMPLE_DT))
    hist, iters = [], 0
    t0 = time.perf_counter()
    for _ in range(nsamples):
        for s in range(per_sample):
            if ops.analyze(1, dt) != 0:
                raise RuntimeError(f"analyze FAILED at dt={dt:g}")
            iters += int(ops.numIter())
        hist.append(ops.nodeDisp(tip, 3))
    wall = time.perf_counter() - t0
    return wall, hist, iters, nsteps


def l2_err(hist, ref):
    num = sum((a - b) ** 2 for a, b in zip(hist, ref))
    den = sum(b * b for b in ref)
    return (num / den) ** 0.5 if den else float("nan")


def main():
    print(f"ADR-77 T3 scheme shootout | n={N} load x{LOAD:g} system={SYSTEM} "
          f"T_end={T_END:g} algorithm=KrylovNewton threads={os.environ['MKL_NUM_THREADS']}")
    print(f"dt ladder: {[f'{d:g}' for d in DT_LADDER]}   reference dt={DT_REF:g}")

    # ---- reference, plus a convergence check ON THE REFERENCE -------------
    # Without this check the whole study rests on an unverified assumption.
    print("\nbuilding reference (Newmark trapezoidal, fine dt)...")
    t0 = time.perf_counter()
    _, h_ref, _, ns_ref = run(DT_REF, ["Newmark", 0.5, 0.25])
    print(f"  dt={DT_REF:g} ({ns_ref} steps)  [{time.perf_counter()-t0:.1f}s]")
    t0 = time.perf_counter()
    _, h_ref2, _, ns_ref2 = run(DT_REF2, ["Newmark", 0.5, 0.25])
    print(f"  dt={DT_REF2:g} ({ns_ref2} steps)  [{time.perf_counter()-t0:.1f}s]")
    ref_self = l2_err(h_ref, h_ref2)
    ok = ref_self < 1e-4
    print(f"  reference self-convergence (L2): {ref_self:.2e}  "
          + ("OK" if ok else "** MARGINAL -- errors near this floor are not resolved **"))
    uz_ref = h_ref2

    out = {"n": N, "load": LOAD, "T_end": T_END, "dt_ladder": DT_LADDER,
           "dt_ref": DT_REF2, "uz_ref": uz_ref, "ref_self_conv": ref_self,
           "algorithm": "KrylovNewton", "schemes": {}}

    for key, label, rho, argv in SCHEMES:
        rec = {"label": label, "rho_inf": rho, "points": {}}
        print(f"\n=== {label}" + (f"  [rho_inf={rho:.3f}]" if rho is not None else "") + " ===")
        for dt in DT_LADDER:
            try:
                samples = []
                for _ in range(REPEATS):
                    w, hist, iters, nsteps = run(dt, argv)
                    samples.append(w)
                wall = statistics.median(samples)
                err = l2_err(hist, uz_ref)
                rec["points"][f"{dt:g}"] = {
                    "dt": dt, "wall_s": wall, "uz_end": hist[-1], "rel_err": err,
                    "iters": iters, "nsteps": nsteps,
                    "iters_per_step": iters / nsteps,
                    "cost_per_sim_s": wall / T_END,
                }
                print(f"  dt={dt:<9g} wall={wall:7.3f}s  err={err:.3e}  "
                      f"it/step={iters/nsteps:6.2f}  cost/sim-s={wall/T_END:9.1f}")
            except RuntimeError as exc:
                rec["points"][f"{dt:g}"] = {"dt": dt, "failed": str(exc)}
                print(f"  dt={dt:<9g} FAILED: {exc}")
        out["schemes"][key] = rec

    (HERE / f"t3_scheme_shootout{TAG}.json").write_text(json.dumps(out, indent=2))
    print(f"\nwrote t3_scheme_shootout{TAG}.json")


if __name__ == "__main__":
    main()
