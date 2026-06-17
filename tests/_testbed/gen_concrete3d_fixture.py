"""gen_concrete3d_fixture.py — dump the numpy oracle's return-map + consistent-tangent
reference numbers to a plain-text fixture the standalone g++ kernel self-check
(concrete3d_kernel_check.cpp) reads and diffs. This is the "oracle-numeric-dump diff"
deliverable (ADR §5): it pins the C++ LadrunoConcrete3DKernel return map + analytic
tangent to the verified numpy oracle at the cross-platform tolerance floors
(stress 1e-7, kappa_p 1e-8, tangent 1e-6).

Run:  python tests/_testbed/gen_concrete3d_fixture.py [output_path]
Writes: tests/_testbed/concrete3d_oracle_fixture.txt   (committed; default), OR the path given as
        argv[1] (the pytest regenerates to a tmp dir so it never overwrites/validates-against-itself
        the committed artifact — PR #249 review).

Param block order (12 numbers): E nu fc ft e Df qh0 Hp Ah Bh Ch Dh
  (the C++ side sets m0 = m0Of(fc,ft,e); both sides derive K,G from E,nu identically.)

Deterministic (no Date/random) — regenerates byte-identically, so CI can assert the committed
fixture is up to date via a fresh regen + diff.
"""
import os
import sys
import numpy as np
import concrete3d_ref as ref

DEFAULT_OUT = os.path.join(os.path.dirname(__file__), "concrete3d_oracle_fixture.txt")


def _mp(E, nu, fc, ft, e, Df, qh0=0.3, Hp=0.5, Ah=0.08, Bh=0.003, Ch=2.0, Dh=1.0e-6):
    return ref.make_material(E, nu, fc, ft, Df=Df, e=e, qh0=qh0, Hp=Hp,
                             Ah=Ah, Bh=Bh, Ch=Ch, Dh=Dh)


def _pblock(mp):
    return [mp["E"], mp["nu"], mp["fc"], mp["ft"], mp["e"], mp["Df"],
            mp["qh0"], mp["Hp"], mp["Ah"], mp["Bh"], mp["Ch"], mp["Dh"]]


def _fmt(xs):
    return " ".join(repr(float(x)) for x in xs)


def run_path(mp, deps_list, hardening):
    """Replay a sequence of strain increments through the oracle's atomic tensor map,
    committing each step (sig_n, kp_n). Returns per-step (sig[6], kp)."""
    sig = np.zeros(6)
    kp = 0.0
    rows = []
    for deps in deps_list:
        sig, kp, _, _ = ref.return_map_tensor(sig, np.asarray(deps, float), mp, kp, hardening)
        rows.append((sig.copy(), kp))
    return rows


def driven_strain_path(mp, eps11_path, hardening, confine="free", sigma3=0.0):
    """Physical uniaxial-STRESS (or active-confined) driver: at each axial strain,
    Newton-solve the lateral strain el so the lateral stress residual vanishes, then
    record the REALIZED full total-strain increment deps[6] and the (sig, kp) the
    atomic tensor map returns. This keeps the path on the well-posed (compressive
    meridian, normal-return) regime — the pure-strain-control path drives into the
    low-kappa_p tiny-surface APEX regime, a documented known gap (handoff §6). Both
    the oracle and the C++ kernel must agree HERE; the apex regime is inherited, not
    gated.  Returns (deps_list, rows[(sig,kp)])."""
    eps = np.zeros(6)
    sig = np.zeros(6)
    kp = 0.0
    el = 0.0
    deps_list, rows = [], []
    for e11 in eps11_path:
        for _ in range(80):                      # inner lateral Newton (FD jacobian)
            deps = np.array([e11 - eps[0], el - eps[1], el - eps[2], 0, 0, 0])
            sn, _, _, _ = ref.return_map_tensor(sig, deps, mp, kp, hardening)
            res = ref.lateral_residual(0.5 * (sn[1] + sn[2]), el, confine, p=sigma3)
            if abs(res) < 1.0e-11 * (mp["fc"] + 1.0):
                break
            d = 1.0e-8 * (abs(el) + 1.0e-6)
            deps2 = np.array([e11 - eps[0], (el + d) - eps[1], (el + d) - eps[2], 0, 0, 0])
            sn2, _, _, _ = ref.return_map_tensor(sig, deps2, mp, kp, hardening)
            res2 = ref.lateral_residual(0.5 * (sn2[1] + sn2[2]), el + d, confine, p=sigma3)
            Jd = (res2 - res) / d
            if abs(Jd) < 1.0e-12:
                Jd = 1.0e-12 if Jd >= 0 else -1.0e-12
            el -= res / Jd
        deps = np.array([e11 - eps[0], el - eps[1], el - eps[2], 0, 0, 0])
        sig, kp, _, _ = ref.return_map_tensor(sig, deps, mp, kp, hardening)
        eps = np.array([e11, el, el, 0, 0, 0])
        deps_list.append(list(deps))
        rows.append((sig.copy(), kp))
    return deps_list, rows


def main(out=None):
    out = out or DEFAULT_OUT
    lines = []
    emitted = []   # (label, pblock, hardening, deps_list, rows)

    e30 = ref.eccentricity_from_kupfer(30.0, 3.0, 1.16)
    mp_pp = _mp(30000.0, 0.2, 30.0, 3.0, e30, 1.0)
    mp_na = _mp(30000.0, 0.2, 30.0, 3.0, e30, 0.3)
    mp_h = _mp(30000.0, 0.2, 30.0, 3.0, e30, 1.0)

    def add_fixed(mp, deps_list, hardening, label):
        emitted.append((label, _pblock(mp), hardening, deps_list, run_path(mp, deps_list, hardening)))

    def add_driven(mp, eps11_path, hardening, label, confine="free", sigma3=0.0):
        dl, rows = driven_strain_path(mp, eps11_path, hardening, confine, sigma3)
        emitted.append((label, _pblock(mp), hardening, dl, rows))

    # ---- perfect-plastic (failure surface) tensor paths ----
    # uniaxial-strain compression/tension + off-meridian shear (well-posed for perfect plasticity)
    add_fixed(mp_pp, [[-0.006 / 60, 0, 0, 0, 0, 0]] * 60, False, "pp_uniax_comp")
    add_fixed(mp_pp, [[0.0008 / 40, 0, 0, 0, 0, 0]] * 40, False, "pp_uniax_tens")
    add_fixed(mp_pp, [[-3.0e-3 / 50, 8.0e-4 / 50, 1.0e-4 / 50, 5.0e-4 / 50, 0, 0]] * 50,
              False, "pp_offmeridian_shear")
    add_fixed(mp_na, [[-2.0e-3 / 50, 6.0e-4 / 50, 2.0e-4 / 50, 1.0e-4 / 50, 0, 0]] * 50,
              False, "pp_nonassoc_shear")

    # ---- hardening (full CDPM2) tensor paths — PHYSICAL uniaxial/confined stress paths ----
    add_driven(mp_h, np.linspace(0, -0.006, 120), True, "hard_uniax_comp", confine="free")
    add_driven(mp_h, np.linspace(0, -0.012, 120), True, "hard_confined_comp",
               confine="active", sigma3=0.10 * 30.0)
    add_driven(mp_h, np.linspace(0, 0.0006, 80), True, "hard_uniax_tens", confine="free")

    lines.append(f"NPATH {len(emitted)}")
    for label, pblock, hardening, deps_list, rows in emitted:
        lines.append(f"PATH {label} {_fmt(pblock)} {1 if hardening else 0} {len(deps_list)}")
        for deps, (sig, kp) in zip(deps_list, rows):
            lines.append(f"{_fmt(deps)}  {_fmt(sig)}  {repr(float(kp))}")

    # ---- tangent reference cases (committed state + deps -> oracle numerical 6x6) ----
    tans = []   # (mp, sig_n[6], kp_n, deps[6], hardening, label)

    def add_tan(mp, deps, hardening, label, prestep=None):
        """Optionally pre-load a committed plastic state by replaying `prestep` deps."""
        sig = np.zeros(6); kp = 0.0
        if prestep is not None:
            for d in prestep:
                sig, kp, _, _ = ref.return_map_tensor(sig, np.asarray(d, float), mp, kp, hardening)
        tans.append((mp, sig.copy(), kp, np.asarray(deps, float), hardening, label))

    # elastic step (perfect-plastic): tangent == elastic C
    add_tan(mp_pp, [-1.0e-5, 0, 0, 0, 0, 0], False, "tan_elastic")
    # plastic non-associated off-meridian w/ shear (the headline asymmetric case)
    add_tan(mp_na, [-2.0e-3, 6.0e-4, 2.0e-4, 1.0e-4, 0, 0], False, "tan_pp_nonassoc")
    # plastic associated (e=1, Df=1) w/ shear
    mp_as = _mp(30000.0, 0.2, 30.0, 3.0, 1.0, 1.0)
    add_tan(mp_as, [-2.0e-3, 6.0e-4, 2.0e-4, 1.0e-4, 0, 0], False, "tan_pp_assoc")
    # plastic, diagonal (axisymmetric -> repeated eigenvalues path)
    add_tan(mp_pp, [-2.0e-3, 5.0e-4, 7.0e-4, 0, 0, 0], False, "tan_pp_diag")
    # hardening tangent from a PHYSICAL mid-plastic committed state (uniaxial-stress driven to
    # ~half the pre-peak ramp), then a small probe increment. Reuses driven_strain_path so the
    # committed (sig_n,kp_n) is well-posed (compressive meridian), not an apex-regime artifact.
    dl_h, rows_h = driven_strain_path(mp_h, np.linspace(0, -0.0025, 60), True, "free")
    sig_h = rows_h[-1][0]; kp_h = rows_h[-1][1]
    tans.append((mp_h, sig_h.copy(), kp_h, np.array([-5.0e-5, 1.0e-5, 1.0e-5, 0, 0, 0]),
                 True, "tan_hard_uniax"))
    # same physical committed state, probe WITH shear => exercises the spectral spin / eigenvector
    # recompose AND the 4x4 hardening principal Jacobian on a non-axisymmetric trial.
    tans.append((mp_h, sig_h.copy(), kp_h, np.array([-5.0e-5, 1.0e-5, 1.0e-5, 8.0e-6, 0, 0]),
                 True, "tan_hard_shear"))

    lines.append(f"NTAN {len(tans)}")
    for mp, sig_n, kp_n, deps, hardening, label in tans:
        C = ref.consistent_tangent(sig_n, deps, mp, kp_n, hardening=hardening)
        lines.append(f"TAN {label} {_fmt(_pblock(mp))} {1 if hardening else 0}")
        lines.append(_fmt(sig_n))
        lines.append(repr(float(kp_n)))
        lines.append(_fmt(deps))
        lines.append(_fmt(C.flatten()))

    with open(out, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"wrote {out}: {len(emitted)} paths, {len(tans)} tangent cases")
    return out


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else None)
