"""Can the soil model fail at all? A single-element deviatoric probe.

The bearing campaign (ADR 79 §8 / vault note 16) found no limit point out to
15 % of the footing width. Before spending anything on meshes, interfaces or
domains, one question has to be settled at the material level: does
`PressureDependMultiYield`, in the parameter set the campaign uses, reach a
DEVIATORIC PLATEAU? If a single Gauss point cannot become perfectly plastic in
shear, no boundary-value problem built from it can develop a bearing limit
point, and every other explanation is moot.

The probe is the standard PDMY simple-shear element test. One unit cube,
bottom face fully fixed, top face tied to a single control node so the four
top nodes translate together — which imposes eps_xx = eps_yy = 0 and makes the
kinematics simple shear exactly. The element is consolidated under a vertical
load in the elastic stage, switched to the plastic stage, then sheared under
DISPLACEMENT control at constant vertical load while the stress tensor is read
back from the element.

Reported per confinement: the shear stress path tau(gamma), whether it
plateaus, and the mobilised friction angle
    sin(phi_mob) = (sigma_1 - sigma_3) / (sigma_1 + sigma_3)
which is the quantity that says whether strength is fully mobilised. `-geom
linear` is used deliberately: the question is what the MATERIAL does, and the
engineering strain the element feeds it is then exactly the imposed strain.

Run (base Python 3.12, see the engine guard in bearing_backbone.py):
    python element_probe.py
"""
import os
import sys
import time

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

if os.path.normcase(os.path.dirname(os.path.abspath(ops.__file__))) != \
        os.path.normcase(_BIN):
    raise SystemExit(f"WRONG ENGINE: {ops.__file__}, expected {_BIN}")

# identical to the campaign's soil (bearing_backbone.py PDMY)
PDMY = dict(rho=2.0, G=5.5e4, B=1.5e5, phi=33.0, gammaPeak=0.1,
            refPress=80.0, d=0.5, PT=27.0, contract=0.05, dil1=0.0, dil2=0.0,
            liq1=0.0, liq2=0.0, liq3=0.0, NYS=20)
GAMMA_MAX = float(os.environ.get("ADR79_PROBE_GAMMA", 0.60))
NSTEP = int(os.environ.get("ADR79_PROBE_NSTEP", 600))
CONF = [float(v) for v in os.environ.get(
    "ADR79_PROBE_CONF", "25,50,100,200,400").split(",")]
PROGRESS = int(os.environ.get("ADR79_PROBE_PROGRESS", 0))   # 0 = quiet


def principal(s):
    """Principal stresses from the OpenSees Voigt order (xx yy zz xy yz zx)."""
    xx, yy, zz, xy, yz, zx = s[:6]
    T = np.array([[xx, xy, zx], [xy, yy, yz], [zx, yz, zz]])
    return np.sort(np.linalg.eigvalsh(T))[::-1]      # s1 >= s2 >= s3


def probe(sv):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for i, (x, y, z) in enumerate([(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
                                   (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)],
                                  start=1):
        ops.node(i, float(x), float(y), float(z))
    for i in (1, 2, 3, 4):
        ops.fix(i, 1, 1, 1)
    # The top face is driven directly on all four nodes, and x is RESTRAINED
    # during consolidation so that eps_xx = eps_yy = 0 and the strain field is
    # genuinely uniform. Two earlier versions of this were wrong and are worth
    # recording. Tying the top face onto a control node and putting the shear
    # `sp` there placed an SP on the RETAINED node of an MP constraint, which
    # the Transformation handler does not support cleanly and which diverged
    # immediately. Leaving x free during consolidation was worse but quieter:
    # the top face SPREAD while the fixed base could not, so the strain was not
    # uniform, sigma_xx != sigma_yy, and a spurious sigma_zx of -6.2 kPa
    # appeared under purely vertical load — i.e. it was a boundary-value
    # problem masquerading as a material point.
    for i in (5, 6, 7, 8):
        ops.fix(i, 0, 1, 0)                 # y restrained; x is prescribed below
    # Hold the top face HORIZONTAL by tying its vertical degrees of freedom.
    # Without this the four top z-DOFs are independent and admit a rocking mode
    # (u_z varying linearly with x) which shear excites and which a single
    # element barely resists: the tangent goes near-singular, and the analysis
    # then either fails outright under a force-residual test or - worse -
    # "converges" under NormDispIncr to a state that violates equilibrium.
    # The tie is on dof 3 only, while the shear sp acts on dof 1, so no degree
    # of freedom carries both an SP and an MP.
    for i in (6, 7, 8):
        ops.equalDOF(5, i, 3)

    p = PDMY
    ops.nDMaterial("PressureDependMultiYield", 1, 3, p["rho"], p["G"], p["B"],
                   p["phi"], p["gammaPeak"], p["refPress"], p["d"], p["PT"],
                   p["contract"], p["dil1"], p["dil2"], p["liq1"], p["liq2"],
                   p["liq3"], p["NYS"])
    ops.element("LadrunoBrick", 1, *range(1, 9), 1, "-geom", "linear")

    # One monotone pseudo-time drives both phases, which is what lets the
    # horizontal degree of freedom stay PRESCRIBED throughout — at zero while
    # the sample consolidates (so eps_xx = 0 and the strain is uniform) and then
    # ramping to impose the shear. Time 0 -> 1 consolidates, 1 -> 2 shears.
    #
    # The obvious alternative — `fix` x for consolidation, then `remove` the fix
    # and add the shear `sp` — does NOT work: `remove('sp', node, dof)` does not
    # match a constraint created by `fix`, which carries no load-pattern tag, so
    # the restraint survives and fights the new sp for the same degree of
    # freedom. That produces no error, just a degenerate system that grinds
    # silently (measured: no progress in 300 s on a four-unknown problem).
    # Avoiding `loadConst` as well keeps the two Path series on one clock.
    ops.timeSeries("Path", 1, "-time", 0.0, 1.0, 2.0,
                   "-values", 0.0, 1.0, 1.0, "-useLast")
    ops.pattern("Plain", 1, 1)
    for i in (5, 6, 7, 8):                  # unit top area, shared four ways
        ops.load(i, 0.0, 0.0, -sv * 0.25)
    ops.timeSeries("Path", 2, "-time", 0.0, 1.0, 2.0,
                   "-values", 0.0, 0.0, GAMMA_MAX, "-useLast")
    ops.pattern("Plain", 2, 2)
    for i in (5, 6, 7, 8):
        ops.sp(i, 1, 1.0)                   # u_x = series_2(t) = gamma
    ops.constraints("Transformation")
    ops.numberer("Plain")
    # FullGeneral, not UmfPack: at four free degrees of freedom UMFPACK's
    # symbolic analysis returns -8 and setSize fails outright, which surfaces as
    # a -1 domainChanged error rather than anything diagnostic. A dense solver
    # is the right instrument at this size and is exact.
    ops.system("FullGeneral")
    ops.test("NormUnbalance", 1e-6 * max(sv, 1.0), 100, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 0.1)
    ops.analysis("Static")
    if PROGRESS:
        print(f"    sv={sv}: consolidating...", flush=True)
    assert ops.analyze(10) == 0, f"consolidation failed at sv={sv}"
    if PROGRESS:
        print(f"    sv={sv}: consolidated, t={ops.getTime():.4g}, "
              f"szz={ops.eleResponse(1, 'stress')[2]:.2f}", flush=True)

    # Switching to the plastic stage redistributes stress (the elastic
    # consolidation state is not on the plastic model's yield surfaces), so
    # this step carries a real residual and needs a realistic tolerance and an
    # accelerated algorithm — 1e-10 under plain Newton does not converge.
    ops.updateMaterialStage("-material", 1, "-stage", 1)
    ops.test("NormUnbalance", 1e-6 * max(sv, 1.0), 400, 0)
    ops.algorithm("KrylovNewton")
    ops.integrator("LoadControl", 0.0)      # hold time, re-equilibrate
    assert ops.analyze(1) == 0, f"stage-1 re-equilibration failed at sv={sv}"
    if PROGRESS:
        print(f"    sv={sv}: stage 1 ok, "
              f"szz={ops.eleResponse(1, 'stress')[2]:.2f}", flush=True)

    # ---- shear: advance pseudo-time from 1 to 2, which ramps series 2 from 0
    # to GAMMA_MAX while series 1 holds the vertical load constant.
    #
    # `DisplacementControl` is the obvious instrument and is the wrong one here:
    # it scales a reference load pattern to hit a target displacement, and the
    # motion we want is imposed by a constraint rather than carried by a load,
    # so there is nothing for it to scale.
    # 1e-6, not 1e-8. The fork's own PDMY element batteries run at 1e-5 with a
    # 1e-4 fallback, and 1e-8 simply does not converge against this material —
    # it burns 100 iterations per step (~18 s on a four-unknown problem) and
    # never arrives. Little rides on it here in any case: gamma is PRESCRIBED,
    # so the strain path and hence the material's tau are fixed exactly, and
    # equilibrium only sets eps_zz on the four free vertical degrees of freedom.
    ops.integrator("LoadControl", 1.0 / NSTEP)
    ops.test("NormUnbalance", 1e-5 * max(sv, 1.0), 200, 0)
    ops.algorithm("KrylovNewton")
    rows = []
    t0 = time.time()
    for k in range(NSTEP):
        if ops.analyze(1) != 0:
            print(f"    sv={sv:6.1f}: diverged at gamma={ops.nodeDisp(5, 1):.4f}",
                  flush=True)
            break
        s = ops.eleResponse(1, "stress")[:6]
        s1, _, s3 = principal(s)
        den = -(s1 + s3)          # compression magnitudes
        rows.append((ops.nodeDisp(5, 1), s[5], s1, s3,
                     (s1 - s3) / den if den != 0 else np.nan))
        if PROGRESS and (k + 1) % PROGRESS == 0:
            print(f"    sv={sv:6.1f} step {k + 1}/{NSTEP} gamma={rows[-1][0]:.4f} "
                  f"tau={rows[-1][1]:8.2f} sin(phi_mob)={rows[-1][4]:.4f} "
                  f"[{time.time() - t0:.0f}s]", flush=True)
    return np.array(rows)


def main():
    phi = np.radians(PDMY["phi"])
    print(f"PDMY phi = {PDMY['phi']}deg  ->  sin(phi) = {np.sin(phi):.4f}, "
          f"gammaPeak = {PDMY['gammaPeak']} at refPress = {PDMY['refPress']} kPa")
    print(f"shearing to gamma = {GAMMA_MAX} in {NSTEP} steps\n")
    print(f"{'sigma_v':>8} {'tau_max':>9} {'gamma@max':>10} {'tau_end':>9} "
          f"{'tau_end/max':>12} {'sin(phi_mob)':>13} {'phi_mob':>8} {'verdict':>12}")
    out = {}
    for sv in CONF:
        r = probe(sv)
        if len(r) < 10:
            print(f"{sv:8.1f}   probe produced no usable path")
            continue
        g, tau = r[:, 0], r[:, 1]
        i = int(np.argmax(tau))
        smob = r[-1, 4]
        # plateau test: over the last 25 % of the path, does tau still climb
        # appreciably? A perfectly plastic response is flat to within a few 1e-3.
        tail = tau[int(0.75 * len(tau)):]
        climb = (tail[-1] - tail[0]) / tail[0]
        verdict = ("PLATEAU" if abs(climb) < 0.02 else
                   "softening" if climb < -0.02 else "still rising")
        print(f"{sv:8.1f} {tau.max():9.2f} {g[i]:10.4f} {tau[-1]:9.2f} "
              f"{tau[-1] / tau.max():12.4f} {smob:13.4f} "
              f"{np.degrees(np.arcsin(min(smob, 1.0))):8.2f} {verdict:>12}")
        out[sv] = r
    if out:
        np.savez_compressed(
            os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         "element_probe.npz"),
            **{f"sv_{int(k)}": v for k, v in out.items()})
        print("\nwrote element_probe.npz")
        print("columns: gamma, tau=-sigma_zx, sigma_1, sigma_3, sin(phi_mob)")


if __name__ == "__main__":
    main()
