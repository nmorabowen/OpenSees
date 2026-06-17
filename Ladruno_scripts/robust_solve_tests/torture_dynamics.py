"""Torture fixture #4 - the rung-5 dynamics fall-through + R-HANDOFF contract.

Exercises the matrix-free DR floor and the atomic implicit<->DR handoff the
robust-solve driver uses when the implicit ladder (and stabilization) is
exhausted (ADR-31 rung 5; LadrunoDynamicRelaxation 33005 + the `ladrunoDR`
settling-signal command).

The model is the softening Concrete02 truss. With control=None and
stabilize=False, plain load control dies at the strength peak and the driver
falls through to DR: it freezes the load (loadConst), relaxes the matrix-free
dynamics to a quasi-static REST state (gated on the mass-free force residual
`ladrunoDR residualNorm`, NOT a physical-mass KE ratio), then restores the load
factor exactly on return. A DR rest state is `regularized` -- a relaxed rest, not
a traced equilibrium branch -- so it is NEVER reported truthy.

Run standalone:
    python Ladruno_scripts/robust_solve_tests/torture_dynamics.py
"""
import os
import sys

DIST = os.environ.get("LADRUNO_DIST",
                      r"C:\Users\nmb\Documents\Github\OpenSees\dist\bin")
if os.path.isdir(DIST):
    os.add_dll_directory(DIST)
    sys.path.insert(0, DIST)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import opensees as ops  # noqa: E402
import torture_softening as soft  # noqa: E402
from robust_drive import robust_drive  # noqa: E402


def dr_command_surface():
    """Bring up a DR integrator directly and read the ladrunoDR getters."""
    soft.build()
    ops.integrator("LoadControl", 0.1)
    ops.analyze(2)
    ops.loadConst("-time", ops.getTime())
    ops.wipeAnalysis()
    ops.constraints("Plain")
    ops.numberer("RCM")
    ops.system("BandGen")
    ops.test("NormUnbalance", 1.0e-3, 25, 0)
    ops.algorithm("Newton")
    ops.integrator("LadrunoDynamicRelaxation")
    ops.analysis("Transient")
    ops.analyze(1, 1.0)
    return ops.ladrunoDR("residualNorm"), ops.ladrunoDR("kineticEnergy")


def run_rung5(settle_tol=1.0e-3, dr_max_steps=3000):
    """Drive softening with control=None, stabilize=False -> rung-5 DR."""
    soft.build()
    ops.integrator("LoadControl", soft.PEAK / 200.0)
    res = robust_drive(
        ops, done=lambda: soft.strain() <= -0.010,
        load_increment=soft.PEAK / 200.0, control=None, stabilize=False,
        dynamics=True, dr_settle_tol=settle_tol, dr_max_steps=dr_max_steps,
        max_substeps=8000)
    return res, ops.getTime()


if __name__ == "__main__":
    rn, ke = dr_command_surface()
    print(f"ladrunoDR residualNorm={rn:.4e}  kineticEnergy={ke:.4e}")

    res, t_after = run_rung5()
    print(f"\n{res}")
    print(f"mode={res.mode}  switches={res.switches}  verdict={res.verdict}")
    print(f"dr_settled={res.dr_settled}  dr_lambda={res.dr_lambda}")
    print(f"getTime() after handoff = {t_after}  (R-HANDOFF: == dr_lambda)")

    ok = (res.mode == "Dynamics" and res.switches >= 1 and res.dr_settled
          and res.dr_lambda is not None and t_after == res.dr_lambda
          and not bool(res) and res.verdict == "regularized")
    print("\nSELF-TEST:", "PASS - rung-5 DR settled and the handoff restored "
          "lambda exactly" if ok else "FAIL - inspect above")
    raise SystemExit(0 if ok else 1)
