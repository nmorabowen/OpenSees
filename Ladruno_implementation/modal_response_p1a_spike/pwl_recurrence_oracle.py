"""
FINAL portable closed form in (omega, d) parametrization, where
    q'' + d q' + w^2 q = f(t),   d = c_alpha / m_alpha  (= 2 xi w for w>0).
This unifies all branches by the discriminant  Delta = w^2 - (d/2)^2  and
handles the rigid-body Rayleigh-damped case (w=0, d>0) that the pure
q''=f rigid branch misses.  Verified vs van-Loan expm FOH oracle for ANY (w,d).
This is the exact algebra to port to C++.
"""
import numpy as np
from scipy.linalg import expm

CRIT_TOL = 1e-8   # relative closeness to the critical discriminant


def _foh_blocks_expm(w, d, dt):
    F = np.array([[0.0, 1.0], [-w * w, -d]])
    G = np.array([[0.0], [1.0]])
    Z = np.zeros((4, 4)); Z[0:2, 0:2] = F; Z[0:2, 2:3] = G; Z[2, 3] = 1.0
    E = expm(Z * dt)
    Phi = E[0:2, 0:2]; g1 = E[0:2, 2]; g2 = E[0:2, 3]
    return Phi, np.column_stack([g1 - g2 / dt, g2 / dt])


def _homog_propagator(w, d, dt):
    """A = e^{-(d/2)dt} [[c+(d/2)sd, sd],[-w^2 sd, c-(d/2)sd]], with
    sd, c chosen by the sign of Delta = w^2-(d/2)^2.  Handles w=0 too."""
    half_d = 0.5 * d
    if w == 0.0:
        if d == 0.0:
            return 1.0, dt, 0.0, 1.0                     # double integrator
        e = np.exp(-d * dt)
        # roots 0, -d :  A = [[1, (1-e)/d],[0, e]]
        return 1.0, (1.0 - e) / d, 0.0, e
    e = np.exp(-half_d * dt)
    Delta = w * w - half_d * half_d                      # w^2 - (d/2)^2
    scale = w * w                                        # for relative crit test
    if abs(Delta) <= CRIT_TOL * scale:
        c = 1.0; sd = dt                                 # critical limit
    elif Delta > 0.0:
        wd = np.sqrt(Delta); c = np.cos(wd * dt); sd = np.sin(wd * dt) / wd
    else:
        wh = np.sqrt(-Delta); c = np.cosh(wh * dt); sd = np.sinh(wh * dt) / wh
    p11 = e * (c + half_d * sd)
    p12 = e * sd
    p21 = e * (-w * w * sd)
    p22 = e * (c - half_d * sd)
    return p11, p12, p21, p22


def closed_form_blocks(w, d, dt):
    p11, p12, p21, p22 = _homog_propagator(w, d, dt)
    A = np.array([[p11, p12], [p21, p22]])

    if w == 0.0:
        if d == 0.0:
            B = np.array([[dt * dt / 3.0, dt * dt / 6.0], [dt / 2.0, dt / 2.0]])
            return A, B
        # w=0, d>0 : particular for q'+... reduce via velocity v=q'
        #   v' + d v = f0 + r t ; q' = v
        #   v_p = (f0 + r t)/d - r/d^2 ; q_p = (f0 t + r t^2/2)/d - r t/d^2
        B = np.zeros((2, 2))
        for j, (f0, f1) in enumerate([(1.0, 0.0), (0.0, 1.0)]):
            r = (f1 - f0) / dt
            def vp(t): return (f0 + r * t) / d - r / (d * d)
            def qp(t): return (f0 * t + r * t * t / 2.0) / d - r * t / (d * d)
            qp0, vp0 = qp(0.0), vp(0.0)
            qh0, qhd0 = -qp0, -vp0
            qh_dt = p11 * qh0 + p12 * qhd0
            qhd_dt = p21 * qh0 + p22 * qhd0
            B[0, j] = qp(dt) + qh_dt
            B[1, j] = vp(dt) + qhd_dt
        return A, B

    # w>0 general (any d): particular  q_p = (f0+rt)/w^2 - d r / w^4
    w2 = w * w
    inv_w2 = 1.0 / w2
    inv_w4 = inv_w2 * inv_w2
    B = np.zeros((2, 2))
    for j, (f0, f1) in enumerate([(1.0, 0.0), (0.0, 1.0)]):
        r = (f1 - f0) / dt
        def qp(t): return (f0 + r * t) * inv_w2 - d * r * inv_w4
        def qpd(t): return r * inv_w2
        qp0, qpd0 = qp(0.0), qpd(0.0)
        qh0, qhd0 = -qp0, -qpd0
        qh_dt = p11 * qh0 + p12 * qhd0
        qhd_dt = p21 * qh0 + p22 * qhd0
        B[0, j] = qp(dt) + qh_dt
        B[1, j] = qpd(dt) + qhd_dt
    return A, B


# --- validation over (w, xi) AND rigid Rayleigh (w=0,d>0) ------------------
print("Closed form (w,d) vs expm-FOH oracle:")
worst = 0.0
tests = []
for f0 in [0.5, 1.0, 2.0, 3.0]:
    for xi in [0.005, 0.05, 0.2, 0.7, 0.999, 1.0, 1.0001, 1.5, 3.0]:
        w = 2 * np.pi * f0
        tests.append((f"w={f0}Hz xi={xi}", w, 2 * xi * w, 0.01))
# rigid modes
tests.append(("rigid d=0", 0.0, 0.0, 0.01))
tests.append(("rigid Rayleigh a0=0.5", 0.0, 0.5, 0.01))
tests.append(("rigid Rayleigh a0=3.0", 0.0, 3.0, 0.02))
# big dt
tests.append(("bigdt under", 2*np.pi*2, 2*0.05*2*np.pi*2, 0.2))

for name, w, d, dt in tests:
    A, B = closed_form_blocks(w, d, dt)
    Ae, Be = _foh_blocks_expm(w, d, dt)
    e = max(np.max(np.abs(A - Ae)), np.max(np.abs(B - Be)))
    worst = max(worst, e)
print(f"  worst over {len(tests)} cases = {worst:.2e}")
assert worst < 1e-9, "MISMATCH"

# rigid Rayleigh full-history sanity: q'' + d q' = f const -> steady q' = f/d
w, d, dt = 0.0, 2.0, 0.001
N = 20000
f = np.ones(N + 1) * 1.0
A, B = closed_form_blocks(w, d, dt)
q = 0.0; qd = 0.0
for n in range(N):
    q, qd = A @ np.array([q, qd]) + B @ np.array([f[n], f[n + 1]])
print(f"  rigid-Rayleigh steady q' -> {qd:.8f} (analytic f/d = {1.0/d:.8f})")
assert abs(qd - 1.0 / d) < 1e-6

print("ALL OK — portable (w,d) closed form verified.")
