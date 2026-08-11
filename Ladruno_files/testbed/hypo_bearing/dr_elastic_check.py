"""dr_elastic_check -- is the DR path finding the SAME equilibrium as Newton?

The DR leg came out ~2x too SOFT against the static leg over the whole curve,
including the first, still-elastic increment.  A factor of two in the elastic
range cannot be plasticity, so this isolates the mechanics from the material:
the SAME mesh, the SAME consistent surcharge and the SAME imposed footing
displacement are driven once by static Newton and once by DR, on an
ELASTIC-ISOTROPIC material where the answer is unique and path-free.  Newton and
DR MUST agree.

It also prints the footing displacement actually reached, because "the imposed
displacement is only partly applied" and "the structure is half as stiff"
produce the same reaction and are otherwise indistinguishable.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import h20_prandtl as H                                          # noqa: E402
import quad_dr_probe as P                                        # noqa: E402

ops = H.ops
log = H.log
DS = 1.0e-3


class A:                       # a stand-in for the argparse namespace
    mass, mass_scale = "gershgorin", 1.0
    damping, zeta = "kinetic", 1.0
    no_auto_refresh, recompute = False, 0
    dt, dtfac = 1.0, 0.5
    chunk, maxit = 500, 40000
    tol_r, tol_res = 1.0e-5, 1.0e-6
    trace = False
    dt_march = 0.5


def setup(order, form):
    nodes, cells, sets, x, y, z = H.strip_mesh(0.5, order)
    w, nfaces = H.consistent_surcharge(nodes, cells, order)
    H.build_model(nodes, cells, sets, w, form, assoc=False, elastic=True)
    ftot = H.Q0 * float(w.sum())
    H.surcharge_step(nodes, sets, w, 1.0e-5 * ftot)
    foot = [int(n) + 1 for n in sets["footing"]]
    r_corr = H.Q0 * float(w[sets["footing"]].sum())
    return foot, r_corr, ftot, ops.nodeDisp(foot[0], 3)


def run(order, form, dtfac):
    log(f"===== ELASTIC {form} order {order} =====")

    # --- static Newton reference -----------------------------------------
    foot, r_corr, ftot, uz0 = setup(order, form)
    R0 = P.reaction(foot, r_corr)
    ops.loadConst("-time", uz0)
    ops.timeSeries("Linear", 2)
    ops.pattern("Plain", 2, 2)
    for t in foot:
        ops.sp(t, 3, 1.0)
    ops.wipeAnalysis()
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("UmfPack")
    ops.test("NormUnbalance", 1.0e-5 * ftot, 25, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", -DS)
    ops.analysis("Static")
    assert ops.analyze(1) == 0, "static push failed"
    Rs = P.reaction(foot, r_corr)
    us = ops.nodeDisp(foot[0], 3)
    log(f"  STATIC : uz {uz0:.9g} -> {us:.9g} (want {uz0 - DS:.9g}, "
        f"err {abs(us - (uz0 - DS)):.2e}); R {R0:.6f} -> {Rs:.6f} kN, "
        f"dR/ds = {(Rs - R0) / DS * 1e-3:.4f} kN/mm")

    # --- DR ---------------------------------------------------------------
    foot, r_corr, ftot, uz0 = setup(order, form)
    R0d = P.reaction(foot, r_corr)
    ops.loadConst("-time", 0.0)
    ops.timeSeries("Constant", 2)
    A.dtfac = dtfac
    A.dt_march = A.dt * dtfac
    P.set_dr_analysis(A)
    P.set_push(foot, uz0 - DS, first=True)
    ok, nit, Rd, ke, res, why, dR = P.relax(A, foot, r_corr, 1.0e-6 * ftot)
    ud = ops.nodeDisp(foot[0], 3)
    log(f"  DR     : uz {uz0:.9g} -> {ud:.9g} (want {uz0 - DS:.9g}, "
        f"err {abs(ud - (uz0 - DS)):.2e}); R {R0d:.6f} -> {Rd:.6f} kN, "
        f"dR/ds = {(Rd - R0d) / DS * 1e-3:.4f} kN/mm  "
        f"[{nit} steps, {why}, res {res:.2e}]")
    log(f"  VERDICT: DR/static secant stiffness = "
        f"{(Rd - R0d) / (Rs - R0):.6f}   "
        f"(1.000000 = the two paths found the same elastic equilibrium)")


if __name__ == "__main__":
    log(f"[control 0] ladrunoBuild() = {ops.ladrunoBuild()}")
    run(1, "bbar", 0.5)
    run(2, "uri", 0.25)
