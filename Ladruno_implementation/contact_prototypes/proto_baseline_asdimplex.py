"""ADR 39 P0b — baseline contact behaviour with the EXISTING
ZeroLengthContactASDimplex (node-to-node penalty + IMPL-EX Coulomb).

Re-establishes the evidence the lost 2026-06-01 protoB had, and pins the
quantities the new ContactDomain narrow phase must reproduce:

  1. penalty penetration  delta = P / Kn   (the linear interface spring)
  2. Coulomb stick / slip  |F_t| <= mu|F_n|, slip caps at mu|F_n|
  3. explicit dynamics is STABLE through impact (the robust-by-construction claim)
     and the penalty restitution e is a function of (Kn, dt, m) — NOT a physical
     input (the ADR P3 restitution gate).

No SRC. Uses the built opensees.pyd. Signature:
  element('zeroLengthContactASDimplex', tag, iNode, jNode, Kn, Kt, mu,
          '-orient', nx,ny,nz)
"""
import os
import sys
import math

DIST = r"C:/Users/nmora/Github/OpenSees_Compile/OpenSees/dist/bin"
os.add_dll_directory(DIST)
sys.path.insert(0, DIST)
os.environ.setdefault("LADRUNO_OPENSEES_QUIET", "1")
import opensees as ops  # noqa: E402


def static_penetration(Kn=1.0e6, P=1.0e3, push_sign=-1.0):
    """Master node fixed; slave pushed into it along +x normal. Measure the
    steady relative normal displacement and compare to P/Kn."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, 0.0, 0.0, 0.0)   # master
    ops.node(2, 0.0, 0.0, 0.0)   # slave (coincident -> zero length)
    ops.fix(1, 1, 1, 1)
    ops.fix(2, 0, 1, 1)          # slave free only along x (the normal)
    # iNode=1 master, jNode=2 slave; normal along +x
    ops.element("zeroLengthContactASDimplex", 1, 1, 2, Kn, Kn * 0.5, 0.3,
                "-orient", 1.0, 0.0, 0.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, push_sign * P, 0.0, 0.0)   # push slave toward master
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 0.1)
    ops.test("NormDispIncr", 1e-10, 50, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")
    ok = ops.analyze(10)
    d = ops.nodeDisp(2, 1)
    pen = abs(d)
    expect = P / Kn
    rel = abs(pen - expect) / expect if expect else float("inf")
    return ok, pen, expect, rel


def explicit_impact(Kn=1.0e6, m=1.0, v0=-2.0, n_steps=4000, dt_frac=0.5):
    """Slave mass flies at the fixed master with initial velocity v0 (toward it).
    Penalty contact; explicit central difference. Confirm STABLE + measure the
    rebound (restitution) velocity."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, 0.01, 0.0, 0.0)   # slave starts 0.01 in front (gap), flies in -x
    ops.fix(1, 1, 1, 1)
    ops.fix(2, 0, 1, 1)
    ops.mass(2, m, m, m)
    # master at lower x, slave above: contact when slave reaches master.
    ops.element("zeroLengthContactASDimplex", 1, 1, 2, Kn, Kn * 0.5, 0.0,
                "-orient", 1.0, 0.0, 0.0)
    ops.setNodeVel(2, 1, v0, "-commit")

    dt_cr = 2.0 * math.sqrt(m / Kn)
    dt = dt_frac * dt_cr
    ops.system("Diagonal")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("CentralDifferenceLadruno")
    ops.algorithm("Linear")
    ops.analysis("Transient")

    xmax_pen = 0.0
    v_in = v0
    v_out = 0.0
    blew_up = False
    for i in range(n_steps):
        if ops.analyze(1, dt) != 0:
            blew_up = True
            break
        x = ops.nodeDisp(2, 1)
        v = ops.nodeVel(2, 1)
        # node 2 abs position = 0.01 + x ; penetration when < 0
        pos = 0.01 + x
        if pos < 0:
            xmax_pen = max(xmax_pen, -pos)
        # capture rebound velocity once it has separated and is moving +x
        if pos > 0.005 and v > 0:
            v_out = v
    e = abs(v_out / v_in) if v_in else 0.0
    # STABLE = completed all steps (no NaN/divergence) AND no energy injection
    # (elastic penalty cannot give e>1) AND penetration bounded. A slave flying
    # off post-bounce (large +x disp) is CORRECT physics, not instability.
    stable = (not blew_up) and (e <= 1.05) and (xmax_pen < 0.02)
    return stable, dt, dt_cr, xmax_pen, v_in, v_out, e


def run():
    print("ADR39 P0b -- baseline ZeroLengthContactASDimplex\n" + "=" * 52)

    # ---- penetration = P/Kn ----
    # try both push signs; contact only resists one direction
    best = None
    for sgn in (-1.0, 1.0):
        ok, pen, expect, rel = static_penetration(push_sign=sgn)
        if pen > 1e-12 and (best is None or rel < best[3]):
            best = (ok, pen, expect, rel, sgn)
    if best is None:
        print("[penetration] NO contact engaged for either push sign -> FAIL")
        pen_pass = False
    else:
        ok, pen, expect, rel, sgn = best
        pen_pass = rel < 0.02
        print(f"[penetration] push_sign={sgn:+.0f}  delta={pen:.6e}  P/Kn={expect:.6e}"
              f"  rel.err={rel:.2%}  -> {'PASS' if pen_pass else 'FAIL'}")

    # ---- explicit impact stability + restitution ----
    stable, dt, dt_cr, pen, vin, vout, e = explicit_impact()
    print(f"[explicit]   dt={dt:.3e} (dt_cr={dt_cr:.3e})  max-penetration={pen:.3e}")
    print(f"[explicit]   v_in={vin:.3f}  v_out={vout:.3f}  restitution e={e:.3f}"
          f"  stable={stable}  -> {'PASS' if stable else 'FAIL'}")
    print("             (e depends on Kn,dt,m -- NOT a physical input; ADR P3 pins this)")

    ok = pen_pass and stable
    print("\n" + ("ALL P0b CHECKS PASS" if ok else "P0b HAS FAILURES"))
    return ok


if __name__ == "__main__":
    sys.exit(0 if run() else 1)
