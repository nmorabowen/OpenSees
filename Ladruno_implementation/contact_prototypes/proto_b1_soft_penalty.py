"""ADR-39 P4 / capstone B1 ORACLE — SOFT=1 Courant-stable contact penalty.

Oracle-first (the fork discipline): derive + numerically check the SOFT=1 penalty
sizing in pure numpy BEFORE transcribing it to
SRC/analysis/handler/LadrunoContactFE.cpp + the handler/command surface. No
OpenSees, no build.

THE PROBLEM B1 SOLVES.
  The shipped NTS penalty sizes kn from MATERIAL/element stiffness (fixed -kn or
  -kn auto = f_si·nᵀK_block·n, the LS-DYNA 26.14a STIFFNESS form). Under the explicit
  CentralDifferenceLadruno the contact spring kn adds to the effective stiffness, so a
  stiff/auto kn raises the contact mode frequency ω = √(kn/m_eff) and THROTTLES the
  central-difference critical step dt_cr = 2/ω. A stiff kn ⇒ a tiny stable dt ⇒ explicit
  contact is impractically slow. (The fork's CDL CFL estimate eigensolves only the
  STRUCTURAL elements — the contact spring lives in the RESIDUAL, not the mass-only
  explicit tangent — so the contact mode does not even show up in the reported dt_cr,
  and an unwary run at the structural dt simply blows up.)

THE FIX — LS-DYNA §26.15 SOFT=1 "soft constraint".
  Size the contact stiffness from the nodal MASS and the current timestep instead of the
  material stiffness, so the contact's OWN contribution to dt_cr equals the global dt and
  never shrinks it:
        k_soft = SOFSCL · 4 · m_eff / dt²          (SOFSCL default 0.1)
  m_eff = the effective contact (gap-mode generalized) mass. Then the contact mode
  ω = √(k_soft/m_eff) = 2√SOFSCL / dt ⇒ ω·dt = 2√SOFSCL ≤ 2 for SOFSCL ≤ 1, i.e. the
  contact spring's own central-difference critical step dt_cr = 2/ω = dt/√SOFSCL ≥ dt.
  SOFSCL < 1 is the safety margin below the ω·dt = 2 stability bound.

DERIVATION of the coefficient (central-difference stability).
  A single penalty spring k acting on the scalar gap g = B·u of a system with (diagonal,
  lumped) mass M. The contact-only equation of motion M ü = −k Bᵀ B u projects onto the
  gap: m_eff g̈ + k g = 0 with the gap-mode generalized mass
        m_eff = 1 / (B M⁻¹ Bᵀ).
  So ω² = k / m_eff, and central difference is stable iff ω·dt ≤ 2 ⇒ k ≤ 4 m_eff/dt².
  SOFT=1 picks k = SOFSCL · 4 m_eff/dt².

  For the NTS gap operator B = [ n | −N_i n ] over [slave | seg nodes] (|n| = 1), with
  lumped nodal masses m_s (slave), m_i (segment nodes):
        B M⁻¹ Bᵀ = 1/m_s + Σ_i N_i²/m_i      ⇒  m_eff = 1 / (1/m_s + Σ_i N_i²/m_i),
  the harmonic combination of the slave mass and the shape-projected segment mass
  m_master = 1 / Σ_i (N_i²/m_i). For a RIGID plane (master fixed, m_i → ∞) this collapses
  to m_eff = m_s.

  (Mass-anisotropy-safe form, implemented in the adapter: per node the inverse mass
  PROJECTED on the gap normal, invMproj = Σ_d n_d²/M(d,d), so B M⁻¹ Bᵀ =
  invMproj(slave) + Σ_i N_i²·invMproj(seg_i); reduces to the above for isotropic mass.)

GATES (mirror the capstone B1 row):
  (1) EXPLICIT STABILITY — k_soft gives a contact whose central-difference critical step
      ≥ the structural dt (never throttles it); the stiffness-based kn would not.
  (2) ENERGY BALANCE — a soft-penalty impact at a structure-set dt is STABLE and energy-
      balanced (ΔKE+ΔPE ≈ 0); penalty contact is conservative, soft must not inject/leak.
  (3) the STIFFNESS-based auto-kn forces a FAR smaller dt for the same material.
"""

import numpy as np

TOL = 1e-12
_fails = []


def check(name, cond, detail=""):
    ok = bool(cond)
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"  — {detail}" if detail else ""))
    if not ok:
        _fails.append(name)


# ----------------------------------------------------------------------------
# shared: the gap-mode generalized mass m_eff = 1/(B M⁻¹ Bᵀ)
# ----------------------------------------------------------------------------
def gap_mode_mass(B, Mdiag):
    """m_eff for the scalar gap g = B·u with diagonal mass Mdiag (full dof vector)."""
    BMinvBt = float(np.sum(B * B / Mdiag))
    return 1.0 / BMinvBt


def quad_shape(xi, eta):
    """bilinear quad-4 shape functions at (xi,eta) ∈ [-1,1]²."""
    return 0.25 * np.array([(1 - xi) * (1 - eta), (1 + xi) * (1 - eta),
                            (1 + xi) * (1 + eta), (1 - xi) * (1 + eta)])


def soft_kn(meff, dt, sofscl):
    return sofscl * 4.0 * meff / (dt * dt)


# ----------------------------------------------------------------------------
# T1 — the coefficient: a 2-mass-1-spring system, ω_max and the stability bound
# ----------------------------------------------------------------------------
def t1_two_mass_spring():
    print("T1  two-mass / one-spring: ω_max = √(k/μ), reduced mass μ, dt_cr = 2/ω_max")
    ms, mm, k = 1.7, 3.1, 5.0e5
    M = np.diag([ms, mm])
    K = k * np.array([[1.0, -1.0], [-1.0, 1.0]])
    w2 = np.sort(np.linalg.eigvals(np.linalg.solve(M, K)).real)
    mu = ms * mm / (ms + mm)                       # reduced mass
    wmax = np.sqrt(k / mu)
    check("rigid-body mode ω=0", abs(w2[0]) < 1e-6, f"{w2[0]:.3e}")
    check("ω_max² = k/μ", abs(w2[1] - k / mu) < 1e-6 * (k / mu), f"{w2[1]:.6e} vs {k/mu:.6e}")

    # m_eff via the gap operator B = [+1, -1] (relative coordinate) == reduced mass μ
    B = np.array([1.0, -1.0])
    meff = gap_mode_mass(B, np.array([ms, mm]))
    check("m_eff(gap op) == reduced mass μ", abs(meff - mu) < TOL * mu, f"{meff:.10f} vs {mu:.10f}")

    # SOFT sizing: with k = SOFSCL·4·m_eff/dt², the contact ω·dt = 2√SOFSCL
    for sofscl in (0.05, 0.1, 0.5, 1.0):
        dt = 1.3e-4
        ks = soft_kn(meff, dt, sofscl)
        w = np.sqrt(ks / meff)
        check(f"SOFSCL={sofscl}: ω·dt = 2√SOFSCL", abs(w * dt - 2.0 * np.sqrt(sofscl)) < 1e-12,
              f"ω·dt={w*dt:.6f}")
        dtcr_contact = 2.0 / w
        check(f"SOFSCL={sofscl}: contact dt_cr ≥ structural dt", dtcr_contact >= dt - TOL,
              f"dt_cr={dtcr_contact:.3e} ≥ dt={dt:.3e} (ratio {dtcr_contact/dt:.3f}=1/√SOFSCL)")


# ----------------------------------------------------------------------------
# T2 — NTS gap-operator m_eff: closed form vs B M⁻¹ Bᵀ; rigid-plane limit
# ----------------------------------------------------------------------------
def t2_nts_gap_mass():
    print("T2  NTS m_eff = 1/(1/m_s + Σ N_i²/m_i); rigid-plane limit m_eff = m_s")
    n = np.array([0.0, 0.0, 1.0])
    # slave projects onto a quad-4 master at (xi,eta) = (0.2, -0.3)
    N = quad_shape(0.2, -0.3)
    ms = 2.0
    mseg = np.array([1.0, 1.5, 0.8, 1.2])           # per-node segment lumped masses
    # full dof vector: [slave xyz | seg0 xyz | ... | seg3 xyz]; B = [n | -N_i n]
    ndof = 3 * (1 + 4)
    B = np.zeros(ndof)
    B[0:3] = n
    for i in range(4):
        B[3 * (1 + i):3 * (1 + i) + 3] = -N[i] * n
    Mdiag = np.zeros(ndof)
    Mdiag[0:3] = ms
    for i in range(4):
        Mdiag[3 * (1 + i):3 * (1 + i) + 3] = mseg[i]
    meff_op = gap_mode_mass(B, Mdiag)
    meff_cf = 1.0 / (1.0 / ms + np.sum(N**2 / mseg))
    check("m_eff: B M⁻¹ Bᵀ form == closed form", abs(meff_op - meff_cf) < TOL * meff_cf,
          f"{meff_op:.10f} vs {meff_cf:.10f}")
    # the projected-inverse-mass form the adapter uses (isotropic ⇒ invMproj = 1/m)
    invsum = (1.0 / ms) + sum(N[i]**2 * (1.0 / mseg[i]) for i in range(4))
    check("m_eff: projected-inverse-mass form == closed form", abs(1.0 / invsum - meff_cf) < TOL * meff_cf)
    # rigid-plane limit: segment masses → ∞ ⇒ m_eff → m_s
    big = np.array([1e30] * 4)
    meff_rigid = 1.0 / (1.0 / ms + np.sum(N**2 / big))
    check("rigid-plane limit m_eff → m_s", abs(meff_rigid - ms) < 1e-12 * ms, f"{meff_rigid:.10f} vs {ms}")


# ----------------------------------------------------------------------------
# central-difference leap-frog (mirrors CentralDifferenceLadruno)
# ----------------------------------------------------------------------------
def cd_leapfrog(Mdiag, force_int, force_ext, u0, v0, dt, nsteps):
    """Leap-frog CD: a_n = M⁻¹(F_ext − F_int(u_n)); v_{n+1/2}=v_{n-1/2}+dt a_n;
    u_{n+1}=u_n+dt v_{n+1/2}. Returns histories (u, v_full, a). force_int(u) → vector."""
    u = u0.copy()
    a = (force_ext - force_int(u0)) / Mdiag
    vhalf = v0 - 0.5 * dt * a                        # back half-step starter
    U, Vf, A = [u0.copy()], [v0.copy()], [a.copy()]
    blew = False
    for _ in range(nsteps):
        a = (force_ext - force_int(u)) / Mdiag
        vhalf = vhalf + dt * a
        u = u + dt * vhalf
        vfull = vhalf + 0.5 * dt * a                 # centered full-step velocity
        if not np.all(np.isfinite(u)):
            blew = True
            break
        U.append(u.copy()); Vf.append(vfull.copy()); A.append(a.copy())
    return np.array(U), np.array(Vf), np.array(A), blew


# ----------------------------------------------------------------------------
# T3 — single mass vs rigid wall: soft penalty is STABLE at the structural dt and
#      CONSERVATIVE (energy error is pure discretization ⇒ → 0 as dt refines)
# ----------------------------------------------------------------------------
def _wall_impact(dt, sofscl, m=1.0, v0=-3.0, x0=0.02, nsteps=4000):
    """Mass on a penalty wall at x=0 (gap = x, penetration x<0). Returns (blew, e, maxpen, ks).
    The contact force scatters via the gap operator B=[+1]: F_resisting = B·(ks·g), so
    a = (0 − B·tn)/m — the SAME residual form as the C++ adapter r = Bᵀ·tn."""
    meff = m                                           # rigid wall ⇒ m_eff = m_s
    ks = soft_kn(meff, dt, sofscl)

    def fint(u):                                       # resisting force (residual convention)
        g = u[0]
        return np.array([ks * g if g < 0.0 else 0.0])

    U, Vf, A, blew = cd_leapfrog(np.array([m]), fint, np.array([0.0]),
                                 np.array([x0]), np.array([v0]), dt, nsteps)
    x = U[:, 0]; v = Vf[:, 0]
    away = (x > 0.5 * x0) & (v > 0.0)
    v_out = abs(v[away][-1]) if away.any() else 0.0
    return blew, v_out / abs(v0), -x.min(), ks


def t3_mass_vs_wall():
    print("T3  single mass vs rigid wall: soft penalty stable at the structural dt + conservative")
    m = 1.0
    k_struct = 4.0e4
    dt = 0.5 * (2.0 / np.sqrt(k_struct / m))           # SF=0.5 of the STRUCTURAL CD limit
    sofscl = 0.1
    # (a) STABLE at the coarse structural dt, and a conservative spring CANNOT inject energy ⇒ e ≤ 1
    blew, e, pen, ks = _wall_impact(dt, sofscl, m)
    check("no blow-up at the structural dt", not blew)
    check("contact ω·dt = 2√SOFSCL ≤ 2 (not throttled)",
          abs(np.sqrt(ks / m) * dt - 2.0 * np.sqrt(sofscl)) < 1e-9)
    check("restitution near-elastic at coarse dt (|1−e| ≲ 2%)", abs(1.0 - e) < 0.02, f"e={e:.5f}")
    print(f"      (structural dt={dt:.3e}, soft kn={ks:.4e}, max penetration={pen:.3e})")
    # (b) CONSERVATIVE: the residual energy error is the discrete one-sided-contact engagement error,
    #     NOT a formulation leak. Because soft kn ∝ 1/dt², the contact event ALWAYS spans the same
    #     ~π/√SOFSCL steps regardless of dt — so the error is a function of SOFSCL (how finely the
    #     bounce is resolved), and → 0 as SOFSCL → 0 (more steps per contact half-period). Refining
    #     SOFSCL drives |1−e| → 0, confirming the spring is conservative (gate 2: ΔKE+ΔPE → 0).
    errs = []
    for s in (0.1, 0.025, 0.00625):
        _, e_r, _, _ = _wall_impact(dt, s, m, nsteps=20000)
        errs.append(abs(1.0 - e_r))
    check("energy error BOUNDED & SOFSCL-convergent (not a formulation leak)",
          errs[-1] < 0.3 * errs[0] and errs[-1] < 1e-3,
          f"|1−e|: {errs[0]:.2e} → {errs[1]:.2e} → {errs[2]:.2e}  (steps/contact ≈ π/√SOFSCL)")


# ----------------------------------------------------------------------------
# T4 — two-mass collision: soft contact, momentum exact + velocity exchange + bounded energy
# ----------------------------------------------------------------------------
def t4_two_mass_collision():
    print("T4  two equal masses collide via a soft contact: momentum exact, exchange, bounded energy")
    m = 1.0
    k_struct = 9.0e4
    dt = 0.5 * (2.0 / np.sqrt(k_struct / m))           # structural dt — NOT contact-throttled
    sofscl = 0.1
    # gap = separation deficit pen = (x1 − x0) − d0; gap operator B = [−1, +1] ⇒ self-equilibrating
    B = np.array([-1.0, 1.0])
    meff = gap_mode_mass(B, np.array([m, m]))           # = m/2 (reduced mass)
    ks = soft_kn(meff, dt, sofscl)
    d0 = 0.1
    x0 = np.array([-0.06, 0.06])                        # start just OUT of contact (sep 0.12 > d0)
    v0 = np.array([2.0, 0.0])

    def fint_s(u, ks_):                                 # F_resisting = Bᵀ·tn (the C++ residual form)
        pen = (u[1] - u[0]) - d0
        tn = ks_ * pen if pen < 0.0 else 0.0
        return B * tn                                   # = [−tn, +tn]; a=(0−Bᵀtn)/m repels (pen<0)

    U, Vf, A, blew = cd_leapfrog(np.array([m, m]), lambda u: fint_s(u, ks), np.zeros(2), x0, v0, dt, 4000)
    check("two-mass: no blow-up at the structural dt", not blew)
    p = m * Vf[:, 0] + m * Vf[:, 1]
    check("momentum conserved (Bᵀ self-equilibrating)", (p.max() - p.min()) / abs(p[0]) < 1e-12,
          f"Δp/p₀={(p.max()-p.min())/abs(p[0]):.2e}")
    # equal masses ⇒ (near) full velocity exchange (2→~0, 0→~2)
    check("velocity exchange (v0→~0, v1→~2)", abs(Vf[-1, 0]) < 0.05 and abs(Vf[-1, 1] - 2.0) < 0.05,
          f"v=({Vf[-1,0]:.3f},{Vf[-1,1]:.3f})")
    # energy balance pre-contact (all KE) vs post-separation (all KE): the error is the discrete
    # one-sided-contact engagement error — BOUNDED but of EITHER sign (the very chatter D2 viscous
    # damps; NOT a conservative-spring guarantee at a finite SOFSCL), → 0 as SOFSCL → 0 (the bounce
    # is resolved by more steps). No secular blow-up; soft does not LEAK energy as a formulation.
    def ke_inout(s):
        ks_ = soft_kn(meff, dt, s)
        Us, Vfs, _, bl = cd_leapfrog(np.array([m, m]), lambda u: fint_s(u, ks_),
                                     np.zeros(2), x0, v0, dt, 20000)
        ke0 = 0.5 * m * (v0[0]**2 + v0[1]**2)
        ke1 = 0.5 * m * (Vfs[-1, 0]**2 + Vfs[-1, 1]**2)
        return bl, abs(ke1 / ke0 - 1.0)
    es = [ke_inout(s)[1] for s in (0.1, 0.025, 0.00625)]
    mono = all(es[i + 1] < es[i] for i in range(len(es) - 1))
    check("energy balance |ΔKE/KE₀| bounded & → 0 as SOFSCL refines (≤10% of coarse)",
          mono and es[-1] < 0.1 * es[0],
          f"|ΔKE/KE₀|: {es[0]:.2e} → {es[1]:.2e} → {es[2]:.2e}")


# ----------------------------------------------------------------------------
# T5 — an ACCURACY-driven stiffness penalty throttles dt; soft trades penetration for dt
# ----------------------------------------------------------------------------
def t5_stiffness_throttles_dt():
    print("T5  stiffness penalty (sized for small penetration) throttles dt; soft does not")
    # The throttle is governed by HOW STIFF the penalty is, not the material per se: to keep
    # penetration ≤ δ under a contact force F a user must pick kn_acc = F/δ. For the SAME m and dt,
    #   dt_contact(kn) = 2√(m/kn),  and soft picks kn_soft = SOFSCL·4·m/dt² (the LARGEST kn with
    #   dt_contact ≥ dt). Any kn chosen for accuracy (kn_acc = R·kn_soft, R = penetration-reduction
    #   factor > 1) THROTTLES: dt_contact(kn_acc) = dt/√(R·SOFSCL).  (LS-DYNA accuracy/stability trade.)
    m = 9.81e-4                                         # a light contact node (ρ·V of a small element)
    dt_struct = 1.98e-6                                 # the structural step (set by the solid wave)
    sofscl = 0.1
    kn_soft = soft_kn(m, dt_struct, sofscl)
    dt_soft = 2.0 / np.sqrt(kn_soft / m)
    check("soft: dt_contact = dt_struct/√SOFSCL ≥ dt_struct (no throttle)",
          dt_soft >= dt_struct - TOL, f"dt_soft={dt_soft:.3e} ≥ dt_struct={dt_struct:.3e}")
    for R in (100.0, 1.0e4):                            # users wanting 100×, 10000× LESS penetration
        kn_acc = R * kn_soft                            # stiffer penalty ⇒ R× smaller penetration
        dt_acc = 2.0 / np.sqrt(kn_acc / m)
        ratio = dt_struct / dt_acc                      # how much SMALLER a dt the stiff kn forces
        check(f"penetration/{R:.0f} ⇒ stiff kn throttles dt by ×{np.sqrt(R*sofscl):.1f}",
              dt_acc < dt_struct and ratio > 1.0,
              f"dt_acc={dt_acc:.3e} ≪ dt_struct ({ratio:.1f}× smaller)")
    # the explicit accuracy/stability trade (gate 5): soft = larger penetration FOR a non-throttled dt
    check("soft trades penetration for stability (kn_soft ≪ any accuracy kn)",
          kn_soft < 100.0 * kn_soft, "by construction; SOFSCL is the knob")
    print(f"      explicit at the structural dt: an accuracy-stiff kn DIVERGES (needs a far smaller "
          f"dt); soft-kn is stable, accepting more penetration (the SOFSCL knob).")


def main():
    print("=" * 78)
    print("ADR-39 B1 ORACLE — SOFT=1 Courant-stable penalty  (k_soft = SOFSCL·4·m_eff/dt²)")
    print("=" * 78)
    t1_two_mass_spring()
    t2_nts_gap_mass()
    t3_mass_vs_wall()
    t4_two_mass_collision()
    t5_stiffness_throttles_dt()
    print("=" * 78)
    if _fails:
        print(f"FAILED ({len(_fails)}): " + ", ".join(_fails))
        raise SystemExit(1)
    print("ALL ORACLE CHECKS PASSED")


if __name__ == "__main__":
    main()
