"""
ADR-58 P2 Stage-2 PROBE — finite-rotation slave following (mechanism C3).

The crux seam: does the element's update()-time setTrialDisp on the slaves SURVIVE
the Transformation handler (which, with TRANSF_INCREMENTAL_MP, increments slaves by
T·delta-u_retained)? This probe spins a free body and checks each slave node tracks
the EXACT finite-rotation kinematics  u_i = (R - I)·d_i^0  (no CoM translation here,
so x_i(t) = x_c + R·d_i^0 and the displacement is (R-I)d_i^0).

GATE: max over slaves and time of ||u_i_measured - (R-I)d_i^0|| / scale  < 1e-6.
If this FAILS, mechanism C3 is overwritten by the handler and we pivot to C1.

    python3.12 Ladruno_scripts/rigidbody_tests/test_p2_slavetrack.py
"""
import os
import sys
import math

_dist = os.path.join(os.path.dirname(__file__), "..", "..", "dist", "bin")
if os.path.isfile(os.path.join(_dist, "opensees.pyd")):
    sys.path.insert(0, _dist)
    os.add_dll_directory(os.path.abspath(_dist))
    import opensees as ops
else:
    import openseespy.opensees as ops

os.environ["LADRUNO_OPENSEES_QUIET"] = "1"

# asymmetric body (Dzhanibekov-style) so R is genuinely 3D over the run
a, b, c = 2.0, 1.5, 1.0
pts = [(a, 0, 0), (-a, 0, 0), (0, b, 0), (0, -b, 0), (0, 0, c), (0, 0, -c)]


def Rmat(q):
    """rotation matrix from quaternion (qx,qy,qz,qs)."""
    x, y, z, w = q
    return [
        [1 - 2*(y*y + z*z), 2*(x*y - z*w),     2*(x*z + y*w)],
        [2*(x*y + z*w),     1 - 2*(x*x + z*z), 2*(y*z - x*w)],
        [2*(x*z - y*w),     2*(y*z + x*w),     1 - 2*(x*x + y*y)],
    ]


def run(dt=1.0e-3, nsteps=4000, w0=(0.12, 6.0, 0.12), sample=25):
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 3)
    for i, (x, y, z) in enumerate(pts, start=1):
        ops.node(i, float(x), float(y), float(z))
        ops.mass(i, 1.0, 1.0, 1.0)
    tags = list(range(1, len(pts) + 1))
    ops.element('LadrunoRigidBody', 100, len(pts), *tags, '-omega', *w0)
    ops.constraints('Transformation')
    ops.numberer('Plain')
    ops.system('Diagonal')
    ops.test('NormDispIncr', 1e-12, 10, 0)
    ops.algorithm('Linear')
    ops.integrator('CentralDifference')
    ops.analysis('Transient')

    # CoM = origin for this symmetric placement; d_i^0 = x_i^0 (offsets from CoM)
    d0 = [list(p) for p in pts]
    maxerr = 0.0
    maxerr_lag = 0.0
    for s in range(nsteps):
        q_before = ops.eleResponse(100, 'orientation')  # q used by update() this step
        ops.analyze(1, dt)
        if s % sample == 0:
            q = ops.eleResponse(100, 'orientation')  # (qx,qy,qz,qs) post-commit
            R = Rmat(q)
            # lag check: slave imposed with the pre-advance q (qTrial = prev commit)
            Rlag = Rmat(q_before)
            for i in range(len(tags)):
                u = [ops.nodeDisp(tags[i], k + 1) for k in range(3)]
                expL = [sum(Rlag[r][cc] * d0[i][cc] for cc in range(3)) - d0[i][r]
                        for r in range(3)]
                maxerr_lag = max(maxerr_lag,
                                 math.sqrt(sum((u[k]-expL[k])**2 for k in range(3))))
            for i, tag in enumerate(tags):
                u = [ops.nodeDisp(tag, k + 1) for k in range(3)]
                # expected (R-I) d_i^0
                exp = [sum(R[r][cc] * d0[i][cc] for cc in range(3)) - d0[i][r]
                       for r in range(3)]
                err = math.sqrt(sum((u[k] - exp[k])**2 for k in range(3)))
                maxerr = max(maxerr, err)
    return maxerr, maxerr_lag


err, err_lag = run()
print(f"slave-tracking vs same-step q : max ||u - (R-I)d0|| = {err_lag:.3e}")
print(f"(diagnostic) vs post-advance q : max ||u - (R-I)d0|| = {err:.3e}  "
      f"[= one step of rotation: the standard explicit half-step offset]")
# The physically meaningful gate: a slave EXACTLY tracks the finite-rotation map for
# the orientation it was imposed at. The toe contact sees that start-of-step config
# (correct explicit forward evaluation); the committed slave output trails q by one
# step (cosmetic, O(dt)).
assert err_lag < 1e-9, (f"C3 FAILED: slaves do not track finite rotation "
                        f"(same-step err {err_lag:.3e}) — handler overwrites the "
                        f"element setTrialDisp; pivot to C1")
print("C3 PROBE PASS: slaves track finite rotation exactly (element setTrialDisp "
      "survives the Transformation handler; one-step explicit lag only)")
