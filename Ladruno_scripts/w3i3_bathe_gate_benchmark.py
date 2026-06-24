"""ADR-52 W3-I3 DECISION GATE — is a tunable implicit Bathe (Noh 2012) worth building?

TRBDF2 (SRC/analysis/integrator/TRBDF2.cpp) already IS the Bathe-2007 Trap+BDF2 composite,
with HARD-CODED constants (TRBDF2.cpp:114-144) -> it is L-stable, i.e. fixed at MAXIMAL
high-frequency dissipation (rho_inf = 0). The Bathe-Noh 2012 variant would make rho_inf
tunable in [0,1] (partial high-frequency dissipation) while keeping the sub-stepped composite
structure. The fork ALREADY has tunable MONOLITHIC dampers: HHT(alpha) and GeneralizedAlpha.

The strategic question this gate answers (build only if YES):
  Is there a regime that needs BOTH (a) tunable/partial high-frequency dissipation AND
  (b) the sub-stepped composite's robustness on sharp transients -- that neither fixed
  TRBDF2 (only rho_inf=0) nor the tunable monolithic HHT/GeneralizedAlpha already cover?

Test problem: a mass with initial velocity impacting an elastic-no-tension (ENT) stop -- the
canonical sharp-stiffness-onset / high-frequency-excitation benchmark from Bathe's own papers.
For an ENT (elastic) impact the EXACT solution conserves energy and the rebound speed equals
the impact speed; under a coarse dt the contact is under-resolved, so the integrators trade
ENERGY CONSERVATION against HIGH-FREQUENCY SMOOTHNESS -- exactly the axis a tunable rho_inf
dials. We map where each integrator sits on that tradeoff, and stress robustness by stiffening
the contact / coarsening dt until schemes fail.
"""
import math
import os
import sys

DIST = os.environ.get("LADRUNO_DIST", r"C:\Users\nmb\Documents\Github\OpenSees\dist\bin")
if os.path.isdir(DIST):
    os.add_dll_directory(DIST)
    sys.path.insert(0, DIST)
os.environ.setdefault("LADRUNO_OPENSEES_QUIET", "1")
import opensees as ops


K, M, V0 = 1.0e6, 1.0, -1.0                            # contact omega=1000, half-period ~3.14e-3 s
TC = math.pi * math.sqrt(M / K)                        # elastic-impact contact duration (half-sine)
F_ANALYTIC = abs(V0) * math.sqrt(K * M)                # exact peak contact force = v0*sqrt(km) = 1000


def run_impact(integ, dt, t_end=0.05, k=K):
    """Mass M at a zeroLength ENT spring (stiff k in compression, 0 in tension), launched at V0
    toward the stop. Returns metrics. integ = ('Name', *args) for ops.integrator."""
    f_analytic = abs(V0) * math.sqrt(k * M)
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(1, 0.0)
    ops.node(2, 0.0)
    ops.fix(1, 1)
    ops.mass(2, M)
    ops.uniaxialMaterial("ENT", 1, k)                  # elastic-no-tension: stiff push, no pull
    ops.element("zeroLength", 1, 1, 2, "-mat", 1, "-dir", 1)
    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("BandGen")
    ops.test("NormDispIncr", 1e-10, 100, 0)
    ops.algorithm("Newton")
    ops.integrator(*integ)
    ops.analysis("Transient")
    ops.setNodeVel(2, 1, V0, "-commit")                # launch toward the stop

    ke0 = 0.5 * M * V0 * V0
    n = int(round(t_end / dt))
    completed = True
    contact_force_peak = 0.0
    n_contact_intervals = 0                            # spurious re-impacts (chatter) > 1
    in_contact = False
    for _ in range(n):
        if ops.analyze(1, dt) != 0:
            completed = False
            break
        f = ops.eleResponse(1, "force")
        fmag = abs(f[0]) if f else 0.0
        contact_force_peak = max(contact_force_peak, fmag)
        active = fmag > 1e-6 * f_analytic
        if active and not in_contact:
            n_contact_intervals += 1
        in_contact = active
    v_reb = ops.nodeVel(2, 1)
    ke_end = 0.5 * M * v_reb * v_reb
    return {
        "completed": completed,
        "energy_retention": (ke_end / ke0) if ke0 > 0 else float("nan"),
        "force_overshoot": contact_force_peak / f_analytic,   # 1.0 = exact; >>1 = HF artifact
        "n_contacts": n_contact_intervals,                    # 1 = clean; >1 = chatter
    }


def run_bar_wave(integ, nel=50, Ebar=1.0e4, rho=1.0, Lbar=1.0, P=1.0, courant=1.0):
    """Fixed-free 1D elastic bar (truss chain), step axial load P at the free end -> stress wave.
    dt = courant * element transit time. Returns mid-bar axial-force overshoot past the analytic
    plateau P and the late-time ripple (peak-to-peak of the gauge force after the wave passed)."""
    Le = Lbar / nel
    c = math.sqrt(Ebar / rho)                          # wave speed
    dt = courant * Le / c
    A = 1.0
    mnode = rho * A * Le                               # lumped nodal mass
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    for i in range(nel + 1):
        ops.node(i + 1, i * Le)
    ops.fix(1, 1)
    for i in range(2, nel + 2):
        ops.mass(i, mnode)
    ops.uniaxialMaterial("Elastic", 1, Ebar)
    for e in range(1, nel + 1):
        ops.element("Truss", e, e, e + 1, A, 1)
    ops.timeSeries("Constant", 1)                      # step load (constant from t=0)
    ops.pattern("Plain", 1, 1)
    ops.load(nel + 1, P)
    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("BandGen")
    ops.test("NormDispIncr", 1e-10, 100, 0)
    ops.algorithm("Newton")
    ops.integrator(*integ)
    ops.analysis("Transient")

    gauge = nel // 2                                   # mid-bar element
    arrival = (gauge * Le) / c                         # time the wavefront reaches the gauge
    nsteps = int(round((1.8 * Lbar / c) / dt))         # ~1.8 transit times
    completed, fmax, after = True, 0.0, []
    for s in range(nsteps):
        if ops.analyze(1, dt) != 0:
            completed = False
            break
        f = ops.eleResponse(gauge, "axialForce")
        fv = f[0] if f else 0.0
        fmax = max(fmax, abs(fv))
        if s * dt > arrival + 0.3 * Lbar / c:          # well after the front passed
            after.append(fv)
    mean_after = sum(after) / len(after) if after else 0.0
    ripple = (max(after) - min(after)) if after else 0.0
    return {"completed": completed, "overshoot": fmax / P, "ripple": ripple}


INTEGRATORS = [
    ("Newmark-trap (rho=1)", ("Newmark", 0.5, 0.25)),
    ("HHT a=0.90", ("HHT", 0.90)),
    ("HHT a=0.80", ("HHT", 0.80)),
    ("HHT a=0.67 (max)", ("HHT", 0.67)),
    ("TRBDF2 (rho=0,comp)", ("TRBDF2",)),
]

# well-resolved -> severely under-resolved (the regime where HF dissipation matters)
DT_FACTORS = [("Tc/8", TC / 8), ("Tc/2", TC / 2), ("Tc", TC), ("2*Tc", 2 * TC)]


def main():
    print("=" * 96)
    print(f"ADR-52 W3-I3 gate — mass impacting an ENT stop  (k={K:g}, m={M:g}, v0={V0:g}; "
          f"contact Tc={TC:.2e}s, exact Fpeak={F_ANALYTIC:g})")
    print("EXACT elastic impact: energy_retention=1.000, force_overshoot=1.000, n_contacts=1")
    print("As dt coarsens past ~Tc the impact is UNDER-RESOLVED -> trapezoidal rings/overshoots;")
    print("high-frequency-dissipative schemes should stay smooth. Where does each sit + survive?")
    print("=" * 96)
    for flabel, dt in DT_FACTORS:
        print(f"\n--- dt = {flabel} = {dt:.3e}s  ({TC/dt:.1f} steps per contact) ---")
        print(f"{'integrator':<22}{'done':>5}{'E_retain':>11}{'F_overshoot':>13}{'n_contacts':>12}")
        for label, integ in INTEGRATORS:
            try:
                r = run_impact(integ, dt)
                print(f"{label:<22}{('Y' if r['completed'] else 'N'):>5}"
                      f"{r['energy_retention']:>11.4f}{r['force_overshoot']:>13.3f}"
                      f"{r['n_contacts']:>12d}")
            except Exception as e:
                print(f"{label:<22} ERROR: {str(e)[:50]}")
    # ROBUSTNESS STRESS: the composite's headline claim is robustness on stiff/sharp transients.
    # Fix a coarse dt and drive the contact stiffness up toward the rigid limit; if HHT FAILS to
    # converge where TRBDF2 SURVIVES, that is the unique-composite-robustness signal (-> GO).
    print("\n" + "=" * 96)
    print("ROBUSTNESS STRESS — fixed dt=1e-3s, stiffness k -> rigid limit (does any scheme uniquely survive?)")
    print("=" * 96)
    for k in (1.0e6, 1.0e8, 1.0e10, 1.0e12):
        steps = math.pi * math.sqrt(M / k) / 1.0e-3
        print(f"\n--- k = {k:.0e}  ({steps:.2e} steps per contact at dt=1e-3) ---")
        print(f"{'integrator':<22}{'done':>5}{'E_retain':>11}{'F_overshoot':>13}")
        for label, integ in INTEGRATORS:
            try:
                r = run_impact(integ, 1.0e-3, k=k)
                print(f"{label:<22}{('Y' if r['completed'] else 'N'):>5}"
                      f"{r['energy_retention']:>11.4f}{r['force_overshoot']:>13.3f}")
            except Exception as e:
                print(f"{label:<22} ERROR: {str(e)[:50]}")

    # COMPOSITE HOME TURF — 1D elastic bar, step end-load -> stress wave. Trapezoidal famously
    # rings (spurious high-frequency oscillation) at the wavefront; this is exactly where the
    # composite / dissipative schemes are advertised. Gauge the mid-bar axial-force overshoot
    # past the analytic plateau P (P/A with A=1). If HHT(tuned) damps the ringing as well as
    # TRBDF2, the composite has no unique home-turf edge either.
    print("\n" + "=" * 96)
    print("WAVE PROPAGATION — fixed-free elastic bar (50 elem), step end load, mid-bar force overshoot")
    print("(exact = clean step to P=1.0 at wave arrival; overshoot>1 = spurious HF ringing)")
    print("=" * 96)
    for cf in (1.0, 4.0):
        print(f"\n--- Courant = {cf:g} (dt = {cf:g} x element transit time) ---")
        print(f"{'integrator':<22}{'done':>5}{'mid_overshoot':>15}{'late_ripple':>14}")
        for label, integ in INTEGRATORS:
            try:
                r = run_bar_wave(integ, courant=cf)
                print(f"{label:<22}{('Y' if r['completed'] else 'N'):>5}"
                      f"{r['overshoot']:>15.3f}{r['ripple']:>14.4f}")
            except Exception as e:
                print(f"{label:<22} ERROR: {str(e)[:50]}")

    print("\n" + "=" * 96)
    print("Verdict logic: if HHT (tunable monolithic) spans clean->dissipative AND stays robust")
    print("at the under-resolved/stiff regime where it matters, a tunable COMPOSITE (Bathe-Noh)")
    print("adds little beyond {fixed TRBDF2, tunable HHT} -> NO-GO. If TRBDF2's composite is")
    print("uniquely robust/smooth there and partial dissipation is also wanted -> GO.")


if __name__ == "__main__":
    main()
