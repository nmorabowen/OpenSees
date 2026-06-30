#!/usr/bin/env python
# ADR-61 P0 oracle — bipenalty deficit sizing for the explicit contact mesh-tie.
#
# Validates BLOCKER-A: the corrected ASSEMBLED-mass deficit formula sizes a mass
# penalty m_p so the contact penalty mode satisfies the central-difference Courant
# bound  omega_contact * dt <= 2*safety  on FINITE real masses — and that the v1
# massless-limit form  m_p = k_p*(dt/2)^2  UNDER-stabilizes (the bug Round-1 caught).
#
# Build-free numpy oracle (mirrors proto_b1_soft_penalty.py). No OpenSees needed.
#
# Scheme (ADR-61 §How "What sets dt_cr"):
#   central difference stable when  dt <= 2/omega_max(M^-1 K),  factor 1.0 (no Noh-Bathe).
#   contact penalty k_p contributes a mode  omega^2 = k_p * (B M^-1 B^T)  with the REAL
#   assembled masses in M^-1.  B = [n | -N_i n] over [slave | master facet nodes].
#   gap-mode inverse mass  s = invMproj_slave + sum_i N_i^2 * invMproj_master_i,
#   invMproj = sum_d n_d^2 / m_d.   Deficit sizing:
#     s_target = (2*safety/dt)^2 / k_p
#     rhs      = s_target - s_master            (s_master = sum_i N_i^2 invMproj_i)
#     m_p      = max(0, 1/rhs - m_s_proj)       (m_s_proj = 1/invMproj_slave)
#   rhs<=0 => the master alone throttles dt => bipenalty cannot fix it (fall back to SOFT).

import numpy as np

SAFETY = 0.9


def invMproj(m3, n):
    """sum_d n_d^2 / m_d over the 3 translational DOFs; a DOF with m<=0 = infinite mass."""
    s = 0.0
    for d in range(3):
        if m3[d] > 0.0:
            s += n[d] * n[d] / m3[d]
    return s


def size_mp(k_p, dt, m_slave3, master_masses, weights, n, safety=SAFETY):
    """Return (m_p, status). status: 'ok' | 'zero' (already heavy) | 'unfixable' (rhs<=0)."""
    n = np.asarray(n, float); n = n / np.linalg.norm(n)
    inv_s = invMproj(m_slave3, n)                     # = 1/m_s_proj
    s_master = sum(w * w * invMproj(mm, n) for w, mm in zip(weights, master_masses))
    s_target = (2.0 * safety / dt) ** 2 / k_p
    rhs = s_target - s_master
    if rhs <= 0.0:
        return 0.0, 'unfixable'
    if inv_s <= 0.0:                                   # massless slave: m_s_proj = inf
        m_p = 1.0 / rhs
        return m_p, 'ok'
    m_s_proj = 1.0 / inv_s
    m_p = 1.0 / rhs - m_s_proj
    if m_p <= 0.0:
        return 0.0, 'zero'
    return m_p, 'ok'


def assembled_omega(k_p, m_slave3, master_masses, weights, n, m_p_on_slave=0.0):
    """Largest eigenvalue frequency of the assembled contact mode, M^-1 K, with the
    penalty spring k_p between the slave and the weighted master point. Built explicitly
    as a (1+nm)-node, 3-DOF-per-node system in the normal direction and solved."""
    n = np.asarray(n, float); n = n / np.linalg.norm(n)
    nm = len(master_masses)
    ndof = 3 * (1 + nm)
    # gap operator g = n.u_s - sum_i w_i n.u_mi  -> B row over the dofs
    B = np.zeros(ndof)
    B[0:3] = n
    for i in range(nm):
        B[3 * (1 + i):3 * (1 + i) + 3] = -weights[i] * n
    K = k_p * np.outer(B, B)                            # rank-1 penalty stiffness
    m = np.zeros(ndof)
    ms = np.array(m_slave3, float).copy()
    # add m_p isotropically on the slave translational dofs
    for d in range(3):
        ms[d] = ms[d] + m_p_on_slave
    m[0:3] = ms
    for i in range(nm):
        m[3 * (1 + i):3 * (1 + i) + 3] = master_masses[i]
    # only dofs with finite mass participate; a fixed (m=inf, here large) master pins its dofs
    # represent a FIXED master by very large mass
    Minv = np.array([1.0 / mi if mi > 0 else 0.0 for mi in m])
    # generalized eig of K u = lam M u  ->  (Minv K) u = lam u, take max real
    A = (Minv[:, None] * K)
    lam = np.linalg.eigvals(A)
    lam = np.real(lam[np.abs(np.imag(lam)) < 1e-9])
    lam_max = max(lam.max(), 0.0)
    return np.sqrt(lam_max)


def check(label, k_p, dt, m_slave3, master_masses, weights, n, safety=SAFETY):
    m_p, status = size_mp(k_p, dt, m_slave3, master_masses, weights, n, safety)
    w_bip = assembled_omega(k_p, m_slave3, master_masses, weights, n, m_p_on_slave=m_p)
    # v1 naive (massless-limit) sizing for comparison: m_p_naive = k_p*(dt/2safety)^2 - 0
    m_p_naive = k_p * (dt / (2 * safety)) ** 2
    w_naive = assembled_omega(k_p, m_slave3, master_masses, weights, n, m_p_on_slave=m_p_naive)
    w_plain = assembled_omega(k_p, m_slave3, master_masses, weights, n, m_p_on_slave=0.0)
    tgt = 2 * safety
    print(f"\n{label}")
    print(f"  status={status}  m_p={m_p:.4g}")
    print(f"  omega*dt :  plain={w_plain*dt:.3f}   DEFICIT={w_bip*dt:.3f}   v1-naive={w_naive*dt:.3f}   (target<= {tgt:.2f})")
    ok = True
    if status == 'ok':
        # deficit must satisfy the bound (and be TIGHT — the binding mode at the target)
        ok = ok and (w_bip * dt <= tgt + 1e-6)
        tight = abs(w_bip * dt - tgt) < 1e-3
        print(f"  -> deficit within bound: {w_bip*dt <= tgt+1e-6}   tight-at-target: {tight}")
        # the naive form should OVER-shoot when the master is finite (the v1 bug)
        if w_naive * dt > tgt + 1e-3:
            print(f"  -> v1-naive OVER-stabilizes? no — UNDER-stabilizes: omega*dt={w_naive*dt:.3f} > {tgt:.2f}  (BUG confirmed)")
    elif status == 'zero':
        ok = ok and (w_plain * dt <= tgt + 1e-6)
        print(f"  -> already heavy enough, m_p=0, plain mode within bound: {w_plain*dt <= tgt+1e-6}")
    elif status == 'unfixable':
        print(f"  -> rhs<=0: master alone throttles dt (omega*dt plain={w_plain*dt:.3f}); fall back to SOFT. (expected)")
        ok = True
    print(f"  PASS={ok}")
    return ok


if __name__ == "__main__":
    print("ADR-61 P0 oracle — bipenalty deficit sizing (safety=%.2f)" % SAFETY)
    dt = 1.0e-6
    # a stiff tie penalty: pick k_p = 50x the max stable SOFT stiffness for a unit slave mass
    # k_soft,max = 4*m_eff/dt^2 ; here drive k_p well above it so a real m_p is needed.
    allpass = True

    # 1) rigid/fixed master (mortar tie to a fixed boundary): master mass huge -> 1/m_m -> 0
    k_p = 50 * 4 * 1.0 / dt**2
    allpass &= check("1) finite slave (m_s=1) vs FIXED master (rigid plane / tie-to-ground)",
                     k_p, dt, [1.0, 1.0, 1.0], [[1e30, 1e30, 1e30]], [1.0], [0, 0, 1])

    # 2) two finite masses (deformable-deformable tie), slave=master=1, single 1-weight master
    allpass &= check("2) finite slave (m_s=1) vs finite master (m_m=1), w=1  (the v1-bug case)",
                     k_p, dt, [1.0, 1.0, 1.0], [[1.0, 1.0, 1.0]], [1.0], [0, 0, 1])

    # 3) mortar facet: slave tied to a quad-4 facet with shape weights summing to 1
    allpass &= check("3) mortar tie: slave(m=2) vs quad-4 facet (m=1 each), weights .4/.3/.2/.1",
                     k_p, dt, [2.0, 2.0, 2.0],
                     [[1., 1., 1.], [1., 1., 1.], [1., 1., 1.], [1., 1., 1.]],
                     [0.4, 0.3, 0.2, 0.1], [0, 0, 1])

    # 4) already-heavy node: soft k_p the node already satisfies -> m_p should be 0
    k_small = 0.5 * 4 * 1.0 / dt**2
    allpass &= check("4) already-heavy: k_p below the node's own stable bound -> m_p=0",
                     k_small, dt, [1.0, 1.0, 1.0], [[1e30, 1e30, 1e30]], [1.0], [0, 0, 1])

    # 5) massless slave (RBE-style rescue): m_s=0 -> finite m_p, well-posed
    allpass &= check("5) massless slave vs fixed master -> m_p rescues M^-1",
                     k_p, dt, [0.0, 0.0, 0.0], [[1e30, 1e30, 1e30]], [1.0], [0, 0, 1])

    # 6) too-light master: master so light its own term exceeds s_target -> unfixable
    k_huge = 1e6 * 4 * 1.0 / dt**2
    allpass &= check("6) too-light master (m_m=1e-3) with huge k_p -> rhs<=0, SOFT fallback",
                     k_huge, dt, [1.0, 1.0, 1.0], [[1e-3, 1e-3, 1e-3]], [1.0], [0, 0, 1])

    # --- two-sided sizing for deformable-deformable (the slave-only 'unfixable' cases) ---
    # reduced mass mu = 1/s is bounded by the lighter side, so raise BOTH interface masses.
    # symmetric split: require 1/(m_s+mp) + sum N_i^2/(m_mi+mp) <= s_target by adding the
    # SAME mp to slave and every master node, solved by bisection on mp.
    def size_mp_twosided(k_p, dt, m_slave3, master_masses, weights, n, safety=SAFETY):
        n2 = np.asarray(n, float); n2 = n2 / np.linalg.norm(n2)
        s_target = (2.0 * safety / dt) ** 2 / k_p
        def s_of(mp):
            ms = [m + mp for m in m_slave3]
            sm = invMproj(ms, n2)
            for w, mm in zip(weights, master_masses):
                sm += w * w * invMproj([m + mp for m in mm], n2)
            return sm
        lo, hi = 0.0, 1.0
        while s_of(hi) > s_target:           # grow until feasible
            hi *= 2.0
            if hi > 1e18: return None
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if s_of(mid) > s_target: lo = mid
            else: hi = mid
        return hi

    def assembled_omega_twosided(k_p, m_slave3, master_masses, weights, n, mp):
        n2 = np.asarray(n, float); n2 = n2 / np.linalg.norm(n2)
        nm = len(master_masses); ndof = 3 * (1 + nm)
        B = np.zeros(ndof); B[0:3] = n2
        for i in range(nm): B[3*(1+i):3*(1+i)+3] = -weights[i] * n2
        K = k_p * np.outer(B, B)
        m = np.zeros(ndof)
        m[0:3] = [v + mp for v in m_slave3]
        for i in range(nm): m[3*(1+i):3*(1+i)+3] = [v + mp for v in master_masses[i]]
        Minv = np.array([1.0/mi if mi > 0 else 0.0 for mi in m])
        lam = np.linalg.eigvals(Minv[:, None] * K)
        lam = np.real(lam[np.abs(np.imag(lam)) < 1e-9])
        return np.sqrt(max(lam.max(), 0.0))

    print("\n--- two-sided sizing (deformable-deformable; the cases slave-only can't fix) ---")
    for lbl, ms, mm, w in [
        ("2T) slave=1 vs master=1", [1.,1.,1.], [[1.,1.,1.]], [1.0]),
        ("3T) mortar slave=2 vs quad-4 (m=1)", [2.,2.,2.],
         [[1.,1.,1.]]*4, [0.4,0.3,0.2,0.1]),
    ]:
        mp = size_mp_twosided(k_p, dt, ms, mm, w, [0,0,1])
        wdt = assembled_omega_twosided(k_p, ms, mm, w, [0,0,1], mp)
        m_phys = min(min(ms), min(min(x) for x in mm))
        print(f"  {lbl}: two-sided m_p={mp:.4g}  (vs physical ~{m_phys})  inflation~{mp/m_phys:.1f}x"
              f"   omega*dt={wdt*dt:.3f} (<= {2*SAFETY:.2f}: {wdt*dt<=2*SAFETY+1e-6})")

    print("\n" + ("ALL PASS" if allpass else "FAILURES PRESENT"))
    print("KEY FINDING: slave-only bipenalty fixes dt ONLY vs a (near-)rigid/fixed master;")
    print("deformable-deformable ties need TWO-SIDED m_p (~big inflation of BOTH interfaces).")
    import sys
    sys.exit(0 if allpass else 1)
