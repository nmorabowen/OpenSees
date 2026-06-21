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
    # COMPRESSION + CONFINEMENT only: these stay on the well-posed compressive meridian. A hardening
    # uniaxial-TENSION path is deliberately NOT included — it is apex-adjacent (the tensile vertex is
    # near), the documented KNOWN-GAP regime where the C++ kernel (post-PR-#249 admissibility fix)
    # correctly bails SAFE (honest non-convergence + elastic fallback) while the numpy oracle does
    # its buggy apex teleport, so they legitimately disagree there. Testing C++<->oracle AGREEMENT in
    # that regime is wrong; the SAFE behaviour is covered instead by the run_robustness fuzzer in
    # concrete3d_kernel_check.cpp (0 inadmissible over 180k random trials incl. tension). Tension is
    # damage-governed (P2) anyway.
    add_driven(mp_h, np.linspace(0, -0.006, 120), True, "hard_uniax_comp", confine="free")
    add_driven(mp_h, np.linspace(0, -0.012, 120), True, "hard_confined_comp",
               confine="active", sigma3=0.10 * 30.0)

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

    # ---- P2 DAMAGE cases (committed damage state + deps -> oracle NOMINAL stress) ----
    # Build a well-posed committed damage state with the path drivers, then probe one more increment;
    # the C++ kernel reproduces it via returnMap (effective return on in.sigEff + damagedUpdate).
    Gf, Gc, As = 0.1, 5.0, 2.0
    dmgs = []   # (label, mp, Gf, Gc, lch, As, state_dict, deps[6], sig_nom[6])

    def add_dmg(label, mp, lch, build_path, deps):
        st = ref.make_damage_state(mp)
        st, _, _, _ = ref._advance_damaged(st, build_path, mp, Gf, Gc, lch, As)
        sig_nom, _, _ = ref.damaged_step_tensor(st, np.asarray(deps, float), mp, Gf, Gc, lch, As)
        dmgs.append((label, mp, lch, st, np.asarray(deps, float), sig_nom))

    lch = 50.0
    # tension-damaged, probe a small further tension increment
    tpath = [np.array([e, 0, 0, 0, 0, 0]) for e in np.linspace(0, 4.5e-4, 300)]
    add_dmg("dmg_tension", mp_h, lch, tpath, [1.0e-6, 0, 0, 0, 0, 0])
    # tension-damaged, probe a COMPRESSION increment (reversal -> unilateral re-split routing)
    add_dmg("dmg_reversal", mp_h, lch, tpath, [-2.0e-6, 0, 0, 0, 0, 0])
    # CONFINED compression-damaged (stress-controlled, off the sigma_lat=0 kink), probe compression
    dconf = ref.drive_damaged_unified(mp_h, np.linspace(0, -0.05, 2000), Gf, Gc, lch, As, sigma3=0.05 * 30.0)
    ic = int(np.argmax(dconf["wc"] > 0.5))
    cpath = [np.array([dconf["eps11"][i], dconf["eps_lat"][i], dconf["eps_lat"][i], 0, 0, 0]) for i in range(ic)]
    add_dmg("dmg_compression", mp_h, lch, cpath,
            [dconf["eps11"][ic] - cpath[-1][0], dconf["eps_lat"][ic] - cpath[-1][1], dconf["eps_lat"][ic] - cpath[-1][2], 0, 0, 0])
    # sheared damaged state, probe with shear (exercises the spectral recompose of the nominal split)
    spath = [np.array([e, -0.2 * e, -0.2 * e, 0.3 * e, 0, 0]) for e in np.linspace(0, 5.0e-4, 300)]
    add_dmg("dmg_shear", mp_na, lch, spath, [3.0e-6, -0.6e-6, -0.6e-6, 0.9e-6, 0, 0])
    # P2g CYCLIC no-heal (discriminating): load tension into softening, THEN elastically unload so the
    # committed sigt_max EXCEEDS the live drive; probe a further unload step. The kernel must drive
    # omega_t with the MONOTONE sigt_max (not the live stress) to reproduce the secant nominal stress —
    # a kernel that re-solved omega against the live drive would HEAL and diverge here.
    cyc_path = ([np.array([e, 0, 0, 0, 0, 0]) for e in np.linspace(0, 9.0e-4, 300)]
                + [np.array([e, 0, 0, 0, 0, 0]) for e in np.linspace(9.0e-4, 5.0e-4, 120)[1:]])
    add_dmg("dmg_cyclic_unload", mp_h, lch, cyc_path, [-2.0e-5, 0, 0, 0, 0, 0])
    # P2h ctTemper (discriminating). _CT maps the mode to the kernel int.
    _CT = {"none": 0, "alphat": 1, "proj": 2}
    # proj: PURE TENSION (loading) — the tensile-stress-projected weight w_t<1 reduces the damage vs the
    # 'none' kernel (discriminating STRESS) AND exercises the d(w_t)/deps tangent term under loading (B4).
    mp_pr = dict(mp_h); mp_pr["ct_temper"] = "proj"
    tp_pr = [np.array([e, 0, 0, 0, 0, 0]) for e in np.linspace(0, 5.0e-4, 300)]
    add_dmg("dmg_cttemper_proj", mp_pr, lch, tp_pr, [2.0e-6, 0, 0, 0, 0, 0])
    # alphat: a BIAXIAL loading state (axial TENSION + lateral COMPRESSION) where alpha_c is partial, so
    # w_t=1-alpha_c < 1 ACTIVELY shields the tensile damage in a smooth loading step: alphat damages LESS
    # than 'none' (wt ~0.78 vs ~0.82 => higher tensile stress), discriminating the mode, while the
    # damaged tangent stays well-defined (a 'none' kernel diverges on the STRESS here). PURE compression
    # gives ac=1 (full shield) but no tensile stress to show it; pure tension gives ac=0 (w_t=1, no shield);
    # the biaxial state exercises a genuine partial w_t in both the stress and the tangent.
    mp_at = dict(mp_h); mp_at["ct_temper"] = "alphat"
    bpth = [np.array([e, -0.5 * e, -0.5 * e, 0, 0, 0]) for e in np.linspace(0, 6.0e-4, 300)]
    add_dmg("dmg_cttemper_alphat", mp_at, lch, bpth, [2.0e-6, -1.0e-6, -1.0e-6, 0, 0, 0])
    # P2i MULTIAXIAL apportioning (discriminating): an UNEQUAL biaxial-tension damaged state — two POSITIVE
    # principals, so E*eps_tilde (Eq.37) drives BOTH the omega-solve and (in the tangent) its gradient
    # E*det_deps ABOVE the extreme tensile principal. A pre-P2i kernel (extreme-principal drive) gives a
    # different nominal stress AND tangent here; uniaxial cases (dmg_tension/_cyclic_unload/_cttemper_proj)
    # stay byte-identical because E*eps_tilde == sig_bar_t there. Unequal (e, 0.6e) keeps the tensile
    # principals DISTINCT (an equal biaxial has a degenerate eigenpair where the FD rotates frozen
    # eigenvectors — the P2e/I4 frozen-eigenvector limitation, not a drive bug).
    bipath = [np.array([e, 0.6 * e, 0, 0, 0, 0]) for e in np.linspace(0, 6.0e-4, 300)]
    add_dmg("dmg_biaxial_tension", mp_h, lch, bipath, [1.0e-6, 0.6e-6, 0, 0, 0, 0])

    lines.append(f"NDMG {len(dmgs)}")
    for label, mp, lch, st, deps, sig_nom in dmgs:
        lines.append(f"DMG {label} {_fmt(_pblock(mp))} {repr(float(Gf))} {repr(float(Gc))} "
                     f"{repr(float(lch))} {repr(float(As))} {_CT[mp.get('ct_temper', 'none')]}")
        lines.append(_fmt(st["eps"]))
        lines.append(_fmt(st["sig_bar"]))
        lines.append(repr(float(st["kp"])))
        lines.append(_fmt([st["et_max"], st["kdt1"], st["kdt2"], st["kdc"], st["kdc1"], st["kdc2"],
                           st["sigt_max"], st["sigc_max"]]))   # P2g monotone-drive history (8 fields)
        lines.append(_fmt(deps))
        lines.append(_fmt(sig_nom))

    # ---- (B5) P3 Tier-2 IMPL-EX: a committed IMPLEX state (with the extrapolation increments) + a
    #      step -> the REPORTED explicit stress sig_rep. Pins the C++ returnMap(implex=true) reported
    #      stress to the oracle damaged_step_implex. Includes a dt-JUMP case (r clamped to implexRmax).
    rmax = float(ref._IMPLEX_RMAX)
    implexes = []   # (label, mp, lch, committed_state, deps[6], dt, sig_rep[6])

    def add_implex(label, mp, build_path, deps, dt):
        st = ref._advance_implex(mp, build_path, Gf, Gc, lch, As)          # uniform dt=1 => committed dt_n=1
        sig_rep, _C, _ns, _d = ref.damaged_step_implex(st, np.asarray(deps, float), dt, mp, Gf, Gc, lch, As)
        implexes.append((label, mp, st, np.asarray(deps, float), float(dt), sig_rep))

    itpath = [np.array([e, 0, 0, 0, 0, 0]) for e in np.linspace(4.5e-4 / 300, 4.5e-4, 300)]   # to softening
    add_implex("implex_tension", mp_h, itpath, [1.0e-6, 0, 0, 0, 0, 0], 1.0)                  # uniform dt (r=1)
    add_implex("implex_reversal", mp_h, itpath, [-2.0e-6, 0, 0, 0, 0, 0], 1.0)                # tension-damaged -> compression
    add_implex("implex_dtjump", mp_h, itpath, [1.0e-6, 0, 0, 0, 0, 0], 10.0)                  # r=10 -> CLAMPED to rmax

    lines.append(f"NIMPLEX {len(implexes)}")
    for label, mp, st, deps, dt, sig_rep in implexes:
        lines.append(f"IMPLEX {label} {_fmt(_pblock(mp))} {repr(float(Gf))} {repr(float(Gc))} "
                     f"{repr(float(lch))} {repr(float(As))} {repr(rmax)}")
        lines.append(_fmt(st["eps"]))
        lines.append(_fmt(st["sig_bar"]))
        lines.append(repr(float(st["kp"])))
        lines.append(_fmt([st["et_max"], st["kdt1"], st["kdt2"], st["kdc"], st["kdc1"], st["kdc2"],
                           st["sigt_max"], st["sigc_max"]]))   # P2g monotone-drive history (8 fields)
        lines.append(_fmt([st["wt"], st["wc"], st["dwt"], st["dwc"], st["dt_n"]]))
        lines.append(_fmt(st["depl"]))
        lines.append(repr(float(dt)))
        lines.append(_fmt(deps))
        lines.append(_fmt(sig_rep))

    # ---- (B6) P3 Duvaut-Lions -eta: a committed (INVISCID) damage state + a step -> the RELAXED nominal
    #      stress at eta>0/dt>0 (viscous) AND at dt<=0 (the inviscid fallback). Pins the C++
    #      returnMap(eta,dt) to the oracle damaged_step_tensor(...,dt) AND verifies the eta->0/dt<=0
    #      byte-identity (= the inviscid Tier-1 path). The committed history is eta-independent (built
    #      inviscid); only the current step relaxes. eta_plastic uses a clean uniaxial-STRESS pre-onset
    #      state (NOT uniaxial-strain, which hits the deep-compression apex chaos) so the viscous-inviscid
    #      gap is a genuine plastic overstress, not a near-elastic ~0 (the PV5b/PV6 oracle lesson).
    etas = []   # (label, mp_eta, lch, state, deps[6], eta, dt, sig_visc[6], sig_inv[6])

    def add_eta(label, mp, lch, build_path, deps, eta, dt):
        st = ref.make_damage_state(mp)
        st, _, _, _ = ref._advance_damaged(st, build_path, mp, Gf, Gc, lch, As)   # inviscid committed history
        me = dict(mp); me["eta"] = eta
        deps = np.asarray(deps, float)
        sig_visc, _, _ = ref.damaged_step_tensor(st, deps, me, Gf, Gc, lch, As, dt=dt)
        sig_inv, _, _ = ref.damaged_step_tensor(st, deps, me, Gf, Gc, lch, As, dt=0.0)   # dt<=0 => inviscid
        etas.append((label, me, lch, st, deps, float(eta), float(dt), sig_visc, sig_inv))

    add_eta("eta_tension", mp_h, lch, tpath, [1.0e-6, 0, 0, 0, 0, 0], 0.5, 1.0)
    add_eta("eta_compression", mp_h, lch, cpath,
            [dconf["eps11"][ic] - cpath[-1][0], dconf["eps_lat"][ic] - cpath[-1][1],
             dconf["eps_lat"][ic] - cpath[-1][2], 0, 0, 0], 0.5, 1.0)
    dpp = ref.drive_damaged_unified(mp_h, np.linspace(0.0, -3.0e-3, 120), Gf, Gc, lch, As)
    ipp = next(k for k in range(1, 120) if 0.5 < dpp["kp"][k] < 0.98 and dpp["wc"][k] < 1.0e-12)
    pppath = [np.array([dpp["eps11"][i], dpp["eps_lat"][i], dpp["eps_lat"][i], 0, 0, 0]) for i in range(ipp)]
    add_eta("eta_plastic", mp_h, lch, pppath,
            [dpp["eps11"][ipp] - pppath[-1][0], dpp["eps_lat"][ipp] - pppath[-1][1],
             dpp["eps_lat"][ipp] - pppath[-1][2], 0, 0, 0], 0.3, 1.0)

    lines.append(f"NETA {len(etas)}")
    for label, me, lch, st, deps, eta, dt, sig_visc, sig_inv in etas:
        lines.append(f"ETA {label} {_fmt(_pblock(me))} {repr(float(Gf))} {repr(float(Gc))} "
                     f"{repr(float(lch))} {repr(float(As))} {repr(eta)} {repr(dt)}")
        lines.append(_fmt(st["eps"]))
        lines.append(_fmt(st["sig_bar"]))
        lines.append(repr(float(st["kp"])))
        lines.append(_fmt([st["et_max"], st["kdt1"], st["kdt2"], st["kdc"], st["kdc1"], st["kdc2"],
                           st["sigt_max"], st["sigc_max"]]))   # P2g monotone-drive history (8 fields)
        lines.append(_fmt(deps))
        lines.append(_fmt(sig_visc))
        lines.append(_fmt(sig_inv))

    with open(out, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"wrote {out}: {len(emitted)} paths, {len(tans)} tangent cases, {len(dmgs)} damage cases, "
          f"{len(implexes)} implex cases, {len(etas)} eta cases")
    return out


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else None)
