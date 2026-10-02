"""WP-144 G2 (Zone B), part 2: LadrunoNorSand wrapped by the fork's LogStrain ND wrapper, at large strain.

DRIVER.  A unit LadrunoBrick `-geom finite` with `nDMaterial LogStrain 2 1` over `nDMaterial LadrunoNorSand 1`;
the eight nodes are fully prescribed (g2_common.build_finite) so every Gauss point sees the uniform deformation
gradient F_k = prod_(j<=k) f_j of the path.  Plan 144 section 2.8: "get finite strain by wrapping the small-strain
kernel in LogStrain"; sheet 144a section 1.4 / 9.5: the return map is identical in log stretches, v = v0 J (exponential update), the
spatial tangent is (S.34).  What the wrapper hands out: material response `stress` = CAUCHY sigma = tau / J
(Kirchhoff tau is the inner's stress), `state` / `elasticStrain` / `D` = the inner's, `strain` = the trial Hencky
strain; the element stiffness carries c + sigma-geometry, from which the AB06 spatial tangent
a^ep_ijkl = A_iJkL F_jJ F_lL (Kirchhoff based, tau(+)1 included) is recovered exactly (g2_common.a4_from_stiffness).

THREE ORACLES, three different questions.
  O2s  O2 in SMALL-strain mode fed the same principal LOG increments ln f (stress read as tau), with the (S.34)
       assembly on its last return (`tangent_finite` of the small-mode state).  That is exactly the algorithm the
       wrapper runs (small-strain kernel, v_{n+1} = v_n exp(tr d_eps), (S.34) assembly on the kernel's tangent), so
       the wrapper must agree with it to round-off: plumbing gate, 1e-10.
  O2f  O2 in FINITE mode (State.finite, v = v0 J): the sheet's reference for the finite-strain path.  Since the G2
       owner decision (2026-10-01, decision 1 = option b; sheet 1.2 / 1.4) the specific-volume update is
       EXPONENTIAL in both modes, v = v0 exp(tr eps), so O2f and O2s are the same algorithm on these coaxial paths
       and the wrapper must equal BOTH at round-off.  (Before the decision the small-strain kernel was linear,
       v = v0 (1 + tr eps), and the wrapper differed from O2f by v0 (1 + x - e^x) ~ -v0 x^2 / 2, x = ln J:
       2e-5 on K2, 1e-4 .. 1.9e-4 on a 20 % drained path, amplified 25-30 x into tau / pi_i once the path dilates:
       2.2e-3 TXC, 5.6e-3 TXE fork.  Those four wrapper-vs-O2f gates were strict xfails, "owner decision pending:
       finite-strain v-update"; the decision is taken and they are REAL gates below.)
  closed forms: v = v0 exp(tr eps) = v0 J with J = det F of the prescribed deformation gradient (sheet 1.2, 1.4,
  13.10; computed here from the nodal F, not from any oracle), objectivity of a rigid rotation.

WHAT IS MEASURED AND WHAT IS A FINDING (reported, never hidden).
  F1  v-update (CLOSED by the owner decision; kept as the record).  The wrapper feeds the kernel exactly ln f on the
      coaxial paths driven here (the recovered elastic strain cancels, see F2), so tr d_eps_feed = ln(J_{n+1}/J_n)
      and the kernel's own v_{n+1} = v_n exp(tr d_eps) is v0 J exactly (sheet 1.4).  Gated at 1e-12 against det F.
      The gate has power: the superseded linear update would read v0 (1 + ln J), off by v0 (1 + x - e^x) (>= 1e-5
      relative on every volume-changing path here, asserted from the closed form); the isochoric undrained path
      (ln J = 0) never discriminated and is the control.
  F2  elastic-strain recovery (CLOSED by the G2 owner decision 2, 2026-10-01, option c).  LogStrainNDMaterial used to
      recover eps^e_{n+1} = inv(D0) tau with D0 the inner's INITIAL tangent ("v1 assumes a linear-elastic inner law").
      NorSand's BA06 energy is nonlinear: K = -p / kappa_hat, mu = mu0 + alpha0 p~.  An isotropic error in eps^e cancels
      in b^e,tr = F_d b^e_n F_d^T -> eps~ - eps^e_n only through the log, not in b^e itself; with alpha0 != 0 the
      deviatoric recovery was wrong and the model was NOT OBJECTIVE: a rigid rotation (no deformation) changed sigma
      (relative to R sigma R^T) and DRIFTED pi_i.  The two alpha0 != 0 objectivity gates were xfail(strict) "owner
      decision 2 pending".  Decision: an inner that carries its own elastic strain implements the mixin
      LadrunoElasticStrainProvider (LadrunoNorSand does) and the wrapper builds b^e = exp(2 eps^e) from the PROVIDED
      strain; every other inner keeps the inv(D0) recovery unchanged.  The two gates are REAL gates below, with the
      provider identity (committed b^e = exp(2 eps^e)) and the non-provider fallback gate.

GATES (the sheet / closed forms / O2; written before the wrapper was run).
  plumbing  wrapper == O2s: tau = J sigma_cauchy, pi_i, v, eps_e, eps_p_v, eps_p_s, D to 1e-10 (kernel_parity
            scales; strain-like zero floor as in test_g2_shell_parity), the a^ep tensor to 1e-9 (the extraction
            contracts a 24 x 24 stiffness with the affine modes: round-off ~1e-12), 1e-7 inside the corner / vertex /
            coalescent bands of kernel_parity (same bands, same reason).
  v-update  v = v0 det F to 1e-12 relative in the wrapper, O2s and O2f; closed-form power check (above).
  finite    wrapper vs O2f (finite mode), every step, every path: tau, pi_i, v, eps_e, eps_p_v, eps_p_s, D at 1e-10,
            and the spatial tangent at 1e-9 / 1e-7 in the bands (sheet 1.4: "no v0 -> v replacement anywhere").
  K2        ordering n_0.7 < n_1.0, gap in [2, 6], nominal bands n_0.7 in [19, 25], n_1.0 in [23, 29] (sheet 14, first-step
            criterion, nominal pi_i0 = -60.4, chi = -3.5, v_c0 = 1.81); wrapper == O2s min-det curve to 1e-6 (a
            cancellation amplifies the 1e-9 tangent round-off by <= 1e2, G1 note); wrapper vs O2f: the same 1e-6 on
            the min-det curve, the same first step, n_interp within 1e-3 step (the 1e-6 det error over the ~0.1
            relative drop of the normalised det per step near the crossing is ~2e-5 step; 50 x margin).
  rotation  alpha0 = 0: Cauchy stress = R sigma R^T to 1e-12 relative, pi_i and v unchanged to 1e-12.
            alpha0 = 2, 50 (G2 decision 2): a 0.2 rad rotation after 10 plastic steps rotates the stress, sigma' = R
            sigma R^T to 1e-10 relative to max|sigma|, leaves pi_i unchanged to 1e-10 |pi_i| (v to 1e-12); a second
            oblique rotation (0.35 rad) and a held step repeat the same at 1e-10: the second rotation starts from the
            committed b^e, shear components included.
  provider  the committed b^e is exp(2 eps^e) of the inner's own elastic strain (read back through a held step after
            an oblique-axis rotation, 1e-12), at alpha0 = 0, 2, 50.
  fallback  ElasticIsotropic and LadrunoJ2 (plastic) under LogStrain keep the inv(D0) recovery: closed-form Hencky
            stress (1e-12), objectivity (1e-10), committed b^e = exp(2 C tau) with C the isotropic compliance (1e-12).

Runtime (measured at G2 round 1, before the decision, whole file, venv py312g2): 50-75 s, of which the two K2
localization runs are 4 s + 7 s and the drained-path generation (brentq on O2f) ~25 s; nothing exceeds 2 minutes,
so nothing is @slow.  Re-measured at G2 close with the provider / fallback gates added: both G2 files together
(113 tests) measured at 97-100 s; still nothing @slow.

MUTANTS killed: wrapper plumbing (tau vs Cauchy: sigma * J instead of / J; the Hencky feed ε_feed bookkeeping;
epsFeed / Be commit; getCopy history) by the plumbing tests on every path (ln J reaches 1.4e-2 on the TXC
path and 1.9e-2 on the TXE path; a Cauchy/Kirchhoff mix-up is a 1e-2 error); the v update (linear v0 (1 + tr),
or v_n + v0 tr: the v == v0 J gate and the O2f gate, 1e-5 relative signal on every volume-changing path) and the
v term of the tangent (vfac = v_n or v0 instead of v_{n+1}: 1.7e-6 / 1.8e-5 on the tangent per the sheet, seen
through the 1e-9 a^ep gates on the plastic dilative paths); the missing tau(+)1 element term (K2 tensor a^ep:
8.6e-3 per the sheet, gate 1e-9); the 1/2 on the spin sum (0.137).  The provider route (G2 decision 2) is gated by mutants MX1 (route disabled), MX2
(shear not doubled) in g2/mutation_gate.md.
"""
from __future__ import annotations

import functools
import math
import time

import numpy as np
import pytest
from scipy.linalg import expm
from scipy.optimize import brentq

import g2_common as G

if G.ops is None:
    pytest.skip(f"opensees.pyd not found/loadable in {G.DIST_BIN}", allow_module_level=True)
ops = G.ops

import o2_algo as O2                       # noqa: E402
import test_kernel_parity as KP            # noqa: E402
STRAIN_NATURAL = G.STRAIN_NATURAL      # the same strain-noise floor as test_g2_shell_parity

# no_gmsh: nothing here meshes; without it tests/conftest.py's session-wide gmsh gate skips this whole file when
# collected from the repo root on a box without gmsh (the G2 harness bug)
pytestmark = [pytest.mark.zone_b, pytest.mark.no_gmsh]

GATE = KP.GATE
BAND_TANGENT_GATE = KP.BAND_TANGENT_GATE
A4_GATE = 1.0e-9
SHEAR_ENG = np.array([1.0, 1.0, 1.0, 2.0, 2.0, 2.0])
I3 = np.eye(3)
SIG0 = -100.0 * I3

# ---- sheet 14 (S.43) protocol and gates
LAM1, LAM2 = 1.0e-3, 4.0e-4
N1 = 10
F1_STRETCH = np.array([1.0 + LAM2, 1.0 - LAM1, 1.0])
F2_STRETCH = np.array([1.0, 1.0 - LAM2, 1.0 + LAM1])
V0_K2 = 1.59
PI0_K2 = -60.4
RHO_07, RHO_10 = (0.7, 0.8), (1.0, 1.0)
GAP_LO, GAP_HI = 2, 6
NOM_N07_BAND, NOM_N10_BAND = (19, 25), (23, 29)
K2_N = 30                              # steps run: past the later localization (n_1.0 ~ 27)
V_EXACT_TOL = 1.0e-12                  # sheet 1.4 / 13.10: v = v0 J to round-off (exponential update, G2 decision)
V_POWER = 1.0e-5                       # the superseded linear update reads v0 (1 + ln J): off by >= this on the
                                       # volume-changing paths here (closed form v0 (1 + x - e^x) ~ v0 x^2 / 2)
MINDET_O2S_TOL = 1.0e-6
MINDET_O2F_TOL = 1.0e-6                # same math as O2s now (module docstring)
NINTERP_O2F_TOL = 1.0e-3               # step: 1e-6 det error / ~0.1 relative det drop per step ~ 2e-5, x50 margin


def k2_params(rho_pair):
    return O2.Params(**dict(KP.K2, rho=rho_pair[0], rho_bar=rho_pair[1])).validate()


# ----------------------------------------------------------------------------------------------
# paths: principal log increments per step
# ----------------------------------------------------------------------------------------------
def k2_logdeps(n):
    return [np.log(F1_STRETCH if k <= N1 else F2_STRETCH) for k in range(1, n + 1)]


def undrained_logdeps(ax_total, n):
    da = ax_total / n
    return [np.array([-0.5 * da, -0.5 * da, da])] * n


def drained_logdeps_finite(P, sigma0, v0, pi0, ax_total, n):
    """Drained triaxial in FINITE mode (O2f): constant CAUCHY lateral stress, tau_lat / J = sigma_lat0, found per
    increment by a scalar root solve on the lateral log increment.  The increments are replayed by every other
    oracle and by the wrapper; the physical drainage is what makes the path a 'drained TXC', nothing else uses it."""
    st = O2.initial_state(P, sigma0, v0, pi0, finite=True)
    lat0 = sigma0[0, 0]
    da = ax_total / n
    out = []
    for _ in range(n):
        def f(dl, st=st):
            s = O2.step(P, st, np.diag([dl, dl, da]))
            if s.flags["refused"]:
                return float("nan")
            return s.sigma[0, 0] / math.exp(2.0 * dl + da) - lat0
        lo, hi = sorted((0.0, -2.0 * da))        # compression expands laterally (dl > 0), extension contracts it
        flo, fhi = f(lo), f(hi)
        assert flo * fhi < 0.0, f"drained root not bracketed: f({lo}) = {flo}, f({hi}) = {fhi}"
        dl = brentq(f, lo, hi, xtol=1e-15, rtol=1e-14, maxiter=200)
        out.append(np.array([dl, dl, da]))
        st = O2.step(P, st, np.diag(out[-1]))
        assert not st.flags["refused"]
    return out


def _make_cases():
    out = {}
    for lab, rp in (("K2_rho0.7", RHO_07), ("K2_rho1.0", RHO_10)):
        out[lab] = dict(P=k2_params(rp), v0=V0_K2, pi0=PI0_K2, sigma0=SIG0, logdeps=lambda: k2_logdeps(K2_N))
    c = KP.case("paper", "drained", KP.DENSE, 40, ax=0.0)
    Pd = KP.o2_params(c)
    v0d = KP.v0_of(Pd, c)
    out["TXC_drained_20pct"] = dict(P=Pd, v0=v0d, pi0=KP.DENSE[0], sigma0=SIG0,
                                    logdeps=lambda P=Pd, v0=v0d: drained_logdeps_finite(
                                        P, SIG0, v0, KP.DENSE[0], math.log(0.8), 40))
    cu = KP.case("paper", "undrained", KP.LOOSE, 30, ax=0.0)
    Pu = KP.o2_params(cu)
    v0u = KP.v0_of(Pu, cu)
    out["TXC_undrained_15pct"] = dict(P=Pu, v0=v0u, pi0=KP.LOOSE[0], sigma0=SIG0,
                                      logdeps=lambda: undrained_logdeps(math.log(0.85), 30))
    ce = KP.case("fork", "drained", KP.DENSE, 40, ax=0.0)
    Pe = KP.o2_params(ce)
    v0e = KP.v0_of(Pe, ce)
    out["TXE_drained_fork_20pct"] = dict(P=Pe, v0=v0e, pi0=KP.DENSE[0], sigma0=SIG0,
                                         logdeps=lambda P=Pe, v0=v0e: drained_logdeps_finite(
                                             P, SIG0, v0, KP.DENSE[0], math.log(1.2), 40))
    return out


CASES = _make_cases()
PATHS = list(CASES)


# ----------------------------------------------------------------------------------------------
# the three oracle runs and the wrapper run (cached per path)
# ----------------------------------------------------------------------------------------------
@functools.lru_cache(maxsize=None)
def logdeps_of(name):
    return [np.asarray(d, float) for d in CASES[name]["logdeps"]()]


@functools.lru_cache(maxsize=None)
def oracles(name):
    c = CASES[name]
    P = c["P"]
    ld = logdeps_of(name)
    s0 = O2.initial_state(P, c["sigma0"], c["v0"], c["pi0"], finite=False)
    f0 = O2.initial_state(P, c["sigma0"], c["v0"], c["pi0"], finite=True)
    small, fin, subs_small = [], [], []
    ss, ff = s0, f0
    for d in ld:
        dd = np.diag(d)
        ssn = O2.step(P, ss, dd)
        ffn = O2.step(P, ff, dd)
        small.append(ssn)
        fin.append(ffn)
        subs_small.append((ss, dd))
        ss, ff = ssn, ffn
        if ssn.flags["refused"] or ffn.flags["refused"]:
            break
    return P, s0, small, fin, subs_small


def stretches(name):
    """cumulative F_k = prod f_j (diagonal) and J_k, from the log increments."""
    Fs, F = [], np.eye(3)
    for d in logdeps_of(name):
        F = np.diag(np.exp(d)) @ F
        Fs.append(F.copy())
    return Fs


@functools.lru_cache(maxsize=None)
def wrapper_run(name, with_a4=True, n_steps=None):
    t0 = time.time()
    c = CASES[name]
    P, s0, small, fin, _ = oracles(name)
    Fs = stretches(name)[:len(small)]
    if n_steps is not None:
        Fs = Fs[:n_steps]
    G.build_finite(G.norsand_args(P, c["v0"], c["pi0"], c["sigma0"]), Fs)
    recs = []
    for k, F in enumerate(Fs):
        rc = ops.analyze(1)
        assert rc == 0, f"{name}: analyze failed at step {k + 1} ({rc})"
        r = dict(F=F, J=float(np.linalg.det(F)), cauchy=G.mat_response("stress"), state=G.mat_response("state"),
                 eps_e=G.mat_response("elasticStrain"), D=G.mat_response("D")[0],
                 info=G.mat_response("stepInfo"), hencky=G.mat_response("strain"))
        if with_a4:
            r["a4"] = G.a4_from_stiffness(F)
        recs.append(r)
    return recs, time.time() - t0


def rel(a, b, natural):
    return KP._rel(a, b, natural)


def nat_scales(P):
    nat = KP._natural_scales(P)
    for q in ("eps_e", "eps_p_v", "eps_p_s"):
        nat[q] = STRAIN_NATURAL
    return nat


# ----------------------------------------------------------------------------------------------
# plumbing: wrapper == O2s (small-strain return on log increments + (S.34) assembly)
# ----------------------------------------------------------------------------------------------
@functools.lru_cache(maxsize=None)
def plumbing_errors(name):
    P, s0, small, fin, subs = oracles(name)
    recs, wall = wrapper_run(name)
    nat = nat_scales(P)
    errs = {q: 0.0 for q in ("tau", "pi_i", "v", "eps_e", "eps_p_v", "eps_p_s", "D", "a4", "a4_band")}
    n_sub, n_band, n_a4, flag_mismatch = 0, 0, 0, []
    D_w, D_o = [], []
    ost = s0
    for k, (o, r) in enumerate(zip(small, recs)):
        st = r["state"]
        tau = r["cauchy"] * r["J"]                               # Kirchhoff = J sigma_cauchy
        errs["tau"] = max(errs["tau"], rel(tau, G.t6(o.sigma), nat["sigma"]))
        errs["pi_i"] = max(errs["pi_i"], rel(st[0], o.pi_i, nat["pi_i"]))
        errs["v"] = max(errs["v"], rel(st[2], o.v, 1.0))
        errs["eps_p_v"] = max(errs["eps_p_v"], rel(st[4], o.eps_p_v, nat["eps_p_v"]))
        errs["eps_p_s"] = max(errs["eps_p_s"], rel(st[5], o.eps_p_s, nat["eps_p_s"]))
        errs["eps_e"] = max(errs["eps_e"], rel(r["eps_e"], G.t6(o.eps_e) * SHEAR_ENG, nat["eps_e"]))
        D_w.append(r["D"])
        D_o.append(o.D)
        for f, val in (("plastic", r["info"][1]), ("substeps", r["info"][6])):
            if int(val) != int(o.flags[f]):
                flag_mismatch.append((k, f, o.flags[f], val))
        # (S.34) assembly of the small-mode return, valid against the wrapper (chained D) only when not substepped
        m = o.flags["substeps"]
        if m == 1:
            a_o = O2.tangent_finite(P, o)
            e = float(np.abs(r["a4"] - a_o).max() / np.abs(a_o).max())
            cv, coal = KP._bands([o])
            if cv or coal:
                errs["a4_band"] = max(errs["a4_band"], e)
                n_band += 1
            else:
                errs["a4"] = max(errs["a4"], e)
            n_a4 += 1
        else:
            n_sub += 1
    if max(abs(x) for x in D_o) > 0.0:
        errs["D"] = max(abs(x - y) for x, y in zip(D_w, D_o)) / max(abs(x) for x in D_o)
    return dict(errs=errs, n=len(recs), n_substepped=n_sub, n_band=n_band, n_a4=n_a4, flag_mismatch=flag_mismatch,
                wall=wall)


@pytest.mark.parametrize("name", PATHS)
def test_wrapper_equals_o2_small_strain_on_log_increments(name):
    """Plumbing gate: LogStrain(LadrunoNorSand) == O2s on every step: Kirchhoff tau = J sigma_cauchy, pi_i, v,
    eps_e, eps_p, D at 1e-10; branch flags equal.  Kills a Cauchy / Kirchhoff mix-up (up to 2e-2 off here), a
    wrong Hencky feed bookkeeping, a commit that does not advance eps_feed / b^e, a shared-state copy."""
    rep = plumbing_errors(name)
    e = rep["errs"]
    print(f"\n[{name}] steps {rep['n']} (substepped {rep['n_substepped']}) wall {rep['wall']:.1f}s  "
          + "  ".join(f"{q}={e[q]:.2e}" for q in ("tau", "pi_i", "v", "eps_e", "eps_p_v", "eps_p_s", "D")))
    assert not rep["flag_mismatch"], rep["flag_mismatch"]
    bad = {q: e[q] for q in ("tau", "pi_i", "v", "eps_e", "eps_p_v", "eps_p_s", "D") if not e[q] <= GATE}
    assert not bad, f"wrapper vs O2s gate {GATE:.0e} exceeded: {bad}"


@pytest.mark.parametrize("name", PATHS)
def test_wrapper_spatial_tangent_equals_o2_s34_assembly(name):
    """The element's consistent tangent through the wrapper (c + geometric term, contracted to a^ep_ijkl) equals the
    sheet (S.34) assembly of the same small-strain return: 1e-9 off the bands, 1e-7 in them.  Kills a missing
    tau(+)1 / geometric term (8.6e-3 per the sheet), the 1/2 on the spin sum (0.137), a wrong J scaling (1e-3), a
    transposed (k, l) pair."""
    rep = plumbing_errors(name)
    e = rep["errs"]
    print(f"\n[{name}] a^ep vs O2s on {rep['n_a4']} unsubstepped steps ({rep['n_band']} in the band, "
          f"{rep['n_substepped']} substepped steps skipped): off-band {e['a4']:.2e}  band {e['a4_band']:.2e}")
    assert rep["n_a4"] >= 5, "too few comparable steps"
    assert e["a4"] <= A4_GATE, e["a4"]
    assert e["a4_band"] <= BAND_TANGENT_GATE, e["a4_band"]


# ----------------------------------------------------------------------------------------------
# F1: the finite-strain v-update
# ----------------------------------------------------------------------------------------------
@pytest.mark.parametrize("name", PATHS)
def test_specific_volume_is_v0_times_detF_in_the_wrapper_and_both_oracles(name):
    """Sheet 1.2 / 1.4 / 13.10 (G2 owner decision 2026-10-01): v = v0 exp(tr eps) and, under LogStrain, v = v0 J with
    J = det F of the prescribed deformation gradient.  The closed form is computed HERE from the nodal F (a product
    of the path's diagonal stretches; no oracle, no shell output) and the committed v of the wrapper, of O2 in
    small-strain mode (O2s) and of O2 in finite mode (O2f) must each equal v0 J to 1e-12 relative, at every step.
    Isochoric path (ln J = 0, undrained TXC): v stays v0 (control; it never discriminated the two updates).
    Kills: a linear v-update v0 (1 + tr eps) or v_n + v0 tr (signal v0 x^2 / 2), v_n + v_n tr (x^2 / 2 again),
    an update that uses tr of the TOTAL feed instead of the increment, a v0 := v slip."""
    P, s0, small, fin, _ = oracles(name)
    recs, _ = wrapper_run(name, with_a4=False)
    v0 = CASES[name]["v0"]
    worst = dict(wrapper=0.0, o2s=0.0, o2f=0.0)
    gap_end = 0.0
    for k, (F, r, o, f) in enumerate(zip(stretches(name), recs, small, fin)):
        J = float(np.linalg.det(F))
        worst["wrapper"] = max(worst["wrapper"], abs(r["state"][2] - v0 * J) / v0)
        worst["o2s"] = max(worst["o2s"], abs(o.v - v0 * J) / v0)
        worst["o2f"] = max(worst["o2f"], abs(f.v - v0 * J) / v0)
        x = math.log(J)
        gap_end = max(gap_end, abs(v0 * (1.0 + x - math.exp(x))) / (v0 * J))      # the superseded rule's miss
    print(f"\n[{name}] |v - v0 det F| / v0, worst over {len(recs)} steps: wrapper {worst['wrapper']:.2e}  "
          f"O2s {worst['o2s']:.2e}  O2f {worst['o2f']:.2e};  superseded linear update would miss by {gap_end:.2e}")
    bad = {k: v for k, v in worst.items() if not v <= V_EXACT_TOL}
    assert not bad, f"v != v0 det F beyond {V_EXACT_TOL:.0e}: {bad}"
    if name == "TXC_undrained_15pct":
        assert gap_end <= 1e-12, f"isochoric control moved: {gap_end:.2e}"
    else:
        assert gap_end >= V_POWER, f"gate has no power on this path: the linear rule would only miss by {gap_end:.2e}"


@functools.lru_cache(maxsize=None)
def finite_errors(name):
    """Wrapper vs O2f (finite mode), max over the steps: tau (Kirchhoff = J sigma_cauchy), pi_i, v, eps_e, eps_p_v,
    eps_p_s at the kernel_parity scales, and D path-wide (as in the plumbing gate); by_step keeps (n, ln J, errors)."""
    P, s0, small, fin, _ = oracles(name)
    recs, _ = wrapper_run(name, with_a4=False)
    nat = nat_scales(P)
    keys = ("tau", "pi_i", "v", "eps_e", "eps_p_v", "eps_p_s")
    errs = {q: 0.0 for q in keys}
    by_step = []
    D_w, D_f = [], []
    for k, (f, r) in enumerate(zip(fin, recs)):
        tau = r["cauchy"] * r["J"]
        st = r["state"]
        e = dict(tau=rel(tau, G.t6(f.sigma), nat["sigma"]), pi_i=rel(st[0], f.pi_i, nat["pi_i"]),
                 v=rel(st[2], f.v, 1.0), eps_e=rel(r["eps_e"], G.t6(f.eps_e) * SHEAR_ENG, nat["eps_e"]),
                 eps_p_v=rel(st[4], f.eps_p_v, nat["eps_p_v"]), eps_p_s=rel(st[5], f.eps_p_s, nat["eps_p_s"]))
        by_step.append((k + 1, float(np.log(r["J"])), e))
        for q in keys:
            errs[q] = max(errs[q], e[q])
        D_w.append(r["D"])
        D_f.append(f.D)
    errs["D"] = (max(abs(x - y) for x, y in zip(D_w, D_f)) / max(abs(x) for x in D_f)
                 if D_f and max(abs(x) for x in D_f) > 0.0 else max([abs(x) for x in D_w] + [0.0]))
    return errs, by_step


@pytest.mark.parametrize("name", PATHS)
def test_wrapper_equals_o2_finite_mode_at_1e_10(name):
    """THE GATE the sheet's finite-strain mapping (1.4) asks for, now a real gate on every path (G2 owner decision
    2026-10-01, exponential v-update): tau, pi_i, v, eps_e, eps_p, D of the LogStrain-wrapped shell equal O2's FINITE
    mode at 1e-10 on every step (kernel_parity scales).  Before the decision the four volume-changing paths were
    strict xfails (2e-5 on K2, ~1e-3 at 20 % strain: the linear v-update); the isochoric undrained path was the
    control that passed.
    Kills: any wrapper / kernel error larger than round-off against the sheet's finite-strain reference: a linear or
    otherwise wrong v-update, a Cauchy/Kirchhoff slip, a wrong v term in the return map's Pi_v, a Hencky-feed
    bookkeeping error."""
    errs, _ = finite_errors(name)
    print(f"\n[{name}] wrapper vs O2 finite: " + "  ".join(f"{q}={v:.2e}" for q, v in errs.items()))
    bad = {q: v for q, v in errs.items() if not v <= GATE}
    assert not bad, f"{name}: wrapper vs O2 finite beyond {GATE:.0e}: {bad}"


@functools.lru_cache(maxsize=None)
def finite_tangent_errors(name):
    """The wrapper's spatial tangent (a^ep_ijkl from the element stiffness) vs (S.34) on O2's FINITE-mode state, on
    the unsubstepped steps (the chained substep tangent has no (S.34) closed form from one return; skipped as in the
    plumbing gate).  Off-band / in-band split exactly as kernel_parity."""
    P, s0, small, fin, _ = oracles(name)
    recs, _ = wrapper_run(name)
    off, band, n_cmp, n_band = 0.0, 0.0, 0, 0
    for f, r in zip(fin, recs):
        if f.flags["substeps"] != 1:
            continue
        a_o = O2.tangent_finite(P, f)
        e = float(np.abs(r["a4"] - a_o).max() / np.abs(a_o).max())
        cv, coal = KP._bands([f])
        if cv or coal:
            band, n_band = max(band, e), n_band + 1
        else:
            off = max(off, e)
        n_cmp += 1
    return dict(off=off, band=band, n=n_cmp, n_band=n_band)


@pytest.mark.parametrize("name", PATHS)
def test_wrapper_spatial_tangent_equals_o2_finite_mode_s34(name):
    """The element's consistent tangent through the wrapper equals (S.34) assembled on O2's FINITE-mode return: 1e-9
    off the bands, 1e-7 in them (the kernel_parity bands).  Sheet 1.4: the local residual, Jacobian and a~^ep are
    identical in finite and small strain with no v0 -> v replacement, so this is the same number as the O2s gate
    above, reached independently through the finite-mode state.
    Kills: a vfac = v_n / v0 term in the tangent (1.7e-6 / 1.8e-5 on a^ep per the sheet), a missing tau(+)1, the 1/2
    on the spin sum, a transposed (k, l) pair."""
    e = finite_tangent_errors(name)
    print(f"\n[{name}] a^ep vs O2 finite (S.34) on {e['n']} unsubstepped steps ({e['n_band']} in the band): "
          f"off-band {e['off']:.2e}  band {e['band']:.2e}")
    assert e["n"] >= 5, "too few comparable steps"
    assert e["off"] <= A4_GATE, e["off"]
    assert e["band"] <= BAND_TANGENT_GATE, e["band"]


def test_v_update_record_superseded_linear_rule_signal():
    """RECORD of the superseded G0/G1 rule (sheet 16.6 item 14), asserted from the closed form only: along each
    volume-changing path the linear update v0 (1 + ln J) misses v0 J by v0 (1 + x - e^x) ~ -v0 x^2 / 2, relative
    size >= V_POWER = 1e-5 at the path end (K2: x = ln J = 6e-3 -> 2e-5; TXC 20 %: 1e-4; TXE fork 20 %: 1.9e-4).  The
    discrimination is therefore ~1e7 x the 1e-12 gate of test_specific_volume_is_v0_times_detF....
    Kills: nothing by itself; it is the proof that those gates have power on every path they claim to cover."""
    print("\nsuperseded linear update, relative miss of v0 J at the path end (closed form):")
    for name in PATHS:
        x = math.log(float(np.linalg.det(stretches(name)[-1])))
        miss = abs(1.0 + x - math.exp(x)) / math.exp(x)
        print(f"  {name:>26s}  ln J {x:+.3e}  miss {miss:.2e}")
        if name != "TXC_undrained_15pct":
            assert miss >= V_POWER, (name, x, miss)
        else:
            assert abs(x) <= 1e-14


# ----------------------------------------------------------------------------------------------
# K2: the localization step through the wrapper
# ----------------------------------------------------------------------------------------------
@functools.lru_cache(maxsize=None)
def k2_localization(name):
    """Wrapper min-det curve (from the element stiffness, n >= N1), plus the O2s curve ((S.34) assembly of the small
    return) and the O2f k2_path result. Returns dict with arrays indexed by step n."""
    t0 = time.time()
    P, s0, small, fin, _ = oracles(name)
    recs, wall = wrapper_run(name)
    dets_w, dets_s = {}, {}
    n_first_w = None
    for n in range(N1, len(recs) + 1):
        a_w = recs[n - 1]["a4"]
        dets_w[n] = O2.acoustic_min_det(P, small[n - 1], a_w)[0]
        dets_s[n] = O2.acoustic_min_det(P, small[n - 1], O2.tangent_finite(P, small[n - 1]))[0]
        if n_first_w is None and dets_w[n] <= 0.0:
            n_first_w = n
            break
    st0 = O2.initial_state(P, SIG0, V0_K2, PI0_K2, finite=True)
    o2f = O2.k2_path(P, st0, K2_N)

    def interp(dets, n_first):
        if n_first is None or n_first < 2:
            return None
        d0, d1 = dets[n_first - 1], dets[n_first]
        return (n_first - 1) + d0 / (d0 - d1)

    n_first_s = next((n for n in sorted(dets_s) if dets_s[n] <= 0.0), None)
    dets_f = {n: float(o2f["min_det"][n - 1]) for n in range(N1, len(o2f["min_det"]) + 1)}
    return dict(dets_w=dets_w, dets_s=dets_s, dets_f=dets_f, n_first_w=n_first_w, n_first_s=n_first_s,
                n_interp_w=interp(dets_w, n_first_w), n_interp_s=interp(dets_s, n_first_s),
                n_first_f=o2f["n_first"], n_interp_f=o2f["n_interp"], wall=time.time() - t0)


def test_k2_localization_ordering_gap_and_bands_through_the_wrapper():
    """Sheet 14 K2 gate, nominal combination, first-step criterion, through OpenSees: rho = 0.7 / rho_bar = 0.8
    localizes strictly before rho = rho_bar = 1, the gap is in [2, 6], and the nominal bands hold.  (The paper has
    22 / 26: sanity, printed, not asserted.)  Kills: a tangent error through the wrapper large enough to move the
    localization (missing tau(+)1, the 1/2 spin factor, a transposed (k, l) pair) or to lose the rho ordering."""
    a, b = k2_localization("K2_rho0.7"), k2_localization("K2_rho1.0")
    n07, n10 = a["n_first_w"], b["n_first_w"]
    print(f"\nK2 through LogStrain(LadrunoNorSand): n_first = ({n07}, {n10}), n_interp = "
          f"({a['n_interp_w']:.2f}, {b['n_interp_w']:.2f}); O2 small+(S.34): ({a['n_first_s']}, {b['n_first_s']}); "
          f"O2 finite k2_path: ({a['n_first_f']}, {b['n_first_f']}) interp ({a['n_interp_f']:.2f}, "
          f"{b['n_interp_f']:.2f}); paper (22, 26)  wall {a['wall']:.0f}s + {b['wall']:.0f}s")
    assert n07 is not None and n10 is not None, (n07, n10)
    assert n07 < n10, f"ordering violated: n_0.7 = {n07} !< n_1.0 = {n10}"
    assert GAP_LO <= n10 - n07 <= GAP_HI, f"gap {n10 - n07} outside [{GAP_LO}, {GAP_HI}]"
    assert NOM_N07_BAND[0] <= n07 <= NOM_N07_BAND[1], n07
    assert NOM_N10_BAND[0] <= n10 <= NOM_N10_BAND[1], n10


@pytest.mark.parametrize("name", ("K2_rho0.7", "K2_rho1.0"))
def test_k2_min_det_curve_equals_o2_small_assembly(name):
    """The wrapper's min acoustic determinant equals the O2s one (same math, (S.34) on the small return) at every
    step n >= 10, normalised by the step-10 value, to 1e-6; and the localization step is the same.
    Kills: the same tangent errors as above, at the resolution of the determinant (a cancellation)."""
    k = k2_localization(name)
    d10 = k["dets_s"][N1]
    worst = max(abs(k["dets_w"][n] - k["dets_s"][n]) / abs(d10) for n in k["dets_w"])
    print(f"\n[{name}] max |det_wrapper - det_O2s| / |det_10| = {worst:.2e} over n = {N1}..{max(k['dets_w'])}")
    assert worst <= MINDET_O2S_TOL, worst
    assert k["n_first_w"] == k["n_first_s"]


@pytest.mark.parametrize("name", ("K2_rho0.7", "K2_rho1.0"))
def test_k2_localization_step_vs_o2_finite_mode(name):
    """Wrapper vs O2's FINITE k2_path, now the same algorithm (exponential v-update in both modes, G2 decision): the
    min-acoustic-determinant curve equal to 1e-6 (normalised by the step-10 value, as against O2s: the 1e-9 tangent
    round-off through the det cancellation), the SAME first localization step, and the interpolated crossing within
    1e-3 step.  Before the decision this was 'within 0.1 step, first step within 1' (the linear v-update shifted the
    curve by 1e-3 relative).
    Kills: a v-update / tangent regression that moves the localization step against the sheet's finite-strain mode."""
    k = k2_localization(name)
    d10 = k["dets_f"][N1]
    common = [n for n in k["dets_w"] if n in k["dets_f"]]
    worst = max(abs(k["dets_w"][n] - k["dets_f"][n]) / abs(d10) for n in common)
    dn = abs(k["n_interp_w"] - k["n_interp_f"])
    print(f"\n[{name}] n_interp wrapper {k['n_interp_w']:.4f}  O2 finite {k['n_interp_f']:.4f}  (diff {dn:.2e});"
          f" n_first wrapper {k['n_first_w']}  O2 finite {k['n_first_f']};  max |det_w - det_f| / |det_10| = "
          f"{worst:.2e} over {len(common)} steps")
    assert len(common) >= 10, "too few common steps to compare the min-det curves"
    assert worst <= MINDET_O2F_TOL, worst
    assert k["n_first_w"] == k["n_first_f"]
    assert dn <= NINTERP_O2F_TOL, dn


# ----------------------------------------------------------------------------------------------
# F2: objectivity.  CLOSED by the G2 owner decision 2 (2026-10-01, option c): LogStrain takes the elastic strain from
# an inner that PROVIDES it (LadrunoElasticStrainProvider, implemented by LadrunoNorSand) instead of recovering it as
# inv(D0) : tau.  Every gate below was written from the closed forms stated in its docstring before it was run.
# ----------------------------------------------------------------------------------------------
ROT_TOL = 1.0e-10          # task gate: stress (rotated back) relative to max|sigma|, pi_i relative to |pi_i|
V_ROT_TOL = 1.0e-12
B_IDENT_TOL = 1.0e-12      # exp(2 eps_e) vs the wrapper's committed b^e: round trip of one eigen-decomposition


def rot_z(theta):
    c, s = math.cos(theta), math.sin(theta)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def rot_axis(axis, theta):
    """Rodrigues: R = I + sin K + (1 - cos) K^2, K = [axis]_x of the unit axis."""
    a = np.asarray(axis, float)
    a = a / np.linalg.norm(a)
    K = np.array([[0.0, -a[2], a[1]], [a[2], 0.0, -a[0]], [-a[1], a[0], 0.0]])
    return I3 + math.sin(theta) * K + (1.0 - math.cos(theta)) * (K @ K)


def drive_rotation(build, F_stretch, R, nsteps=10, nhold=1, R2=None):
    """nsteps of the stretch increment F_stretch (plastic history), then ONE step that rotates the body by R with no
    deformation (F -> R F), then (if R2 is given) a SECOND rigid rotation (F -> R2 R F), then nhold steps that HOLD
    the last F (zero increment).  Returns snapshots {before, rot, [rot2], hold} of the material response at the last
    step of each phase: Cauchy stress (3x3), state vector, elasticStrain and the wrapper's Hencky strain (engineering
    Voigt).  A held step cannot see an error in the committed b^e by itself (it feeds the inner eps_tr - eps_n = 0
    whatever b^e is); only a LATER NON-ZERO increment does, hence the second rotation."""
    Fs, F = [], np.eye(3)
    for _ in range(nsteps):
        F = np.diag(F_stretch) @ F
        Fs.append(F.copy())
    Fs.append(R @ F)
    last = R @ F
    if R2 is not None:
        last = R2 @ R @ F
        Fs.append(last)
    Fs += [last] * nhold
    build(Fs)

    def snap():
        return dict(sig=G.m3(G.mat_response("stress")), state=G.mat_response("state"),
                    eps_e=G.mat_response("elasticStrain"), hencky=G.mat_response("strain"))
    for _ in range(nsteps):
        assert ops.analyze(1) == 0
    out = dict(before=snap())
    assert ops.analyze(1) == 0
    out["rot"] = snap()
    if R2 is not None:
        assert ops.analyze(1) == 0
        out["rot2"] = snap()
    for _ in range(nhold):
        assert ops.analyze(1) == 0
    out["hold"] = snap()
    return out


def norsand_rotation(alpha0, R, nsteps=10, R2=None):
    P = O2.Params(**dict(KP.K2, rho=0.7, rho_bar=0.8, alpha0=alpha0)).validate()
    return P, drive_rotation(lambda Fs: G.build_finite(G.norsand_args(P, V0_K2, PI0_K2, SIG0), Fs),
                             F1_STRETCH, R, nsteps, R2=R2)


def rigid_rotation_violation(alpha0, nsteps=10, theta=0.2):
    """After nsteps of the K2 f1 increment (plastic state) a rigid rotation R_z(theta) is applied, then a second rigid
    rotation R2 (oblique axis (0.3, -1, 0.5), 0.35 rad) and one held step.  Closed form: no deformation => Cauchy
    sigma_new = R sigma R^T at each rotation, and the state (pi_i, v) is unchanged throughout.  The second rotation
    starts from the b^e COMMITTED after the first (shear components included), the held step from the second.
    Returns a dict of the measured violations."""
    R = rot_z(theta)
    R2 = rot_axis((0.3, -1.0, 0.5), 0.35)
    P, r = norsand_rotation(alpha0, R, nsteps, R2=R2)
    b, a, a2, h = r["before"], r["rot"], r["rot2"], r["hold"]
    exp = R @ b["sig"] @ R.T
    exp2 = R2 @ a["sig"] @ R2.T
    scale = float(np.abs(b["sig"]).max())
    return dict(stress=float(np.abs(a["sig"] - exp).max() / scale),
                pi_i=abs(a["state"][0] - b["state"][0]) / abs(b["state"][0]),
                v=abs(a["state"][2] - b["state"][2]),
                size=float(np.abs(exp - b["sig"]).max()),
                stress2=float(np.abs(a2["sig"] - exp2).max() / scale),
                pi_i2=abs(a2["state"][0] - b["state"][0]) / abs(b["state"][0]),
                size2=float(np.abs(exp2 - a["sig"]).max()),
                hold_stress=float(np.abs(h["sig"] - exp2).max() / scale),
                hold_pi_i=abs(h["state"][0] - b["state"][0]) / abs(b["state"][0]),
                eps_p_s=float(b["state"][5]), pi_i_moved=abs(b["state"][0] - PI0_K2) / abs(PI0_K2))


def _print_rot(tag, e):
    print(f"\n{tag}: rotation 0.2 rad after 10 plastic steps (|R s R^T - s| = {e['size']:.2f} kPa): stress "
          f"{e['stress']:.2e}  pi_i {e['pi_i']:.2e}  v {e['v']:.2e}  | second rotation 0.35 rad (|dsigma| = "
          f"{e['size2']:.2f} kPa): stress {e['stress2']:.2e}  pi_i {e['pi_i2']:.2e}  | one held step: stress "
          f"{e['hold_stress']:.2e}  pi_i {e['hold_pi_i']:.2e}")


def _assert_nonvacuous(e):
    assert e["size"] > 1.0 and e["size2"] > 1.0, "the rotations must actually move the stress"
    assert e["eps_p_s"] > 0.0 and e["pi_i_moved"] > 1.0e-3, "the history must be plastic (pi_i evolved) before the rotation"


def test_rigid_rotation_is_objective_for_constant_shear_modulus():
    """alpha0 = 0 (BA06 energy: mu constant): a rigid rotation after a plastic K2 history changes nothing, nor does a
    further held step.  (Before the provider route this held through an exact deviatoric recovery; it is the control
    that must not regress.)
    Kills: a wrapper that commits b^e from a wrong elastic strain / mis-rotates the committed b^e (objectivity)."""
    e = rigid_rotation_violation(0.0)
    _print_rot("alpha0 = 0", e)
    _assert_nonvacuous(e)
    assert e["stress"] <= 1e-12 and e["pi_i"] <= 1e-12 and e["v"] <= 1e-12
    assert e["stress2"] <= ROT_TOL and e["pi_i2"] <= ROT_TOL
    assert e["hold_stress"] <= ROT_TOL and e["hold_pi_i"] <= ROT_TOL


@pytest.mark.parametrize("alpha0", (2.0, 50.0))
def test_rigid_rotation_is_objective_with_pressure_dependent_shear_modulus(alpha0):
    """REAL GATE (G2 owner decision 2, 2026-10-01, option c; was xfail(strict) 'owner decision 2 pending').  With
    alpha0 != 0 the BA06 shear modulus mu = mu0 + alpha0 p~ is pressure dependent, so tau != D0 : eps^e and the old
    inv(D0) recovery was not objective.  Closed form (objectivity of a rigid rotation, no oracle): a 0.2 rad rotation
    applied after 10 plastic steps leaves the state unchanged and rotates the stress, sigma' = R sigma R^T:
      stress  |sigma' - R sigma R^T| / max|sigma|  <= 1e-10,
      pi_i    |pi_i' - pi_i| / |pi_i|              <= 1e-10,
      v       |v' - v|                             <= 1e-12 (J does not change),
    then a SECOND rigid rotation (oblique axis (0.3, -1, 0.5), 0.35 rad) must again give sigma_new = R2 sigma R2^T and
    an unchanged pi_i, v at 1e-10, and ONE held step (zero increment) leaves that state unchanged at 1e-10.  The
    second rotation is the sensitive one for the shear convention of the provided eps^e: it starts from the b^e
    COMMITTED at the first rotation, whose shear components are the only place a tensor-vs-engineering slip can show
    (the coaxial steps before it have none, and a held step feeds the inner eps_tr - eps_n = 0 whatever b^e is).
    Kills: MX1 (route disabled: inv(D0) fallback), MX2 (provider returns tensor instead of engineering shear), a
    provider that returns the committed instead of the trial eps^e, a wrapper that ignores the provided strain."""
    e = rigid_rotation_violation(alpha0)
    _print_rot(f"alpha0 = {alpha0}", e)
    _assert_nonvacuous(e)
    assert e["stress"] <= ROT_TOL, e["stress"]
    assert e["pi_i"] <= ROT_TOL, e["pi_i"]
    assert e["v"] <= V_ROT_TOL, e["v"]
    assert e["stress2"] <= ROT_TOL, e["stress2"]
    assert e["pi_i2"] <= ROT_TOL, e["pi_i2"]
    assert e["hold_stress"] <= ROT_TOL, e["hold_stress"]
    assert e["hold_pi_i"] <= ROT_TOL, e["hold_pi_i"]


def _tensor_from_eng(v):
    return G.m3(np.asarray(v, float) * np.array([1.0, 1.0, 1.0, 0.5, 0.5, 0.5]))


@pytest.mark.parametrize("alpha0", (0.0, 2.0, 50.0))
def test_provider_identity_committed_be_is_exp_two_eps_e(alpha0):
    """The committed elastic left Cauchy-Green tensor of LogStrain(LadrunoNorSand) is exp(2 eps^e) of the material's
    OWN elastic strain (owner decision 2, option c).  The committed b^e is not an output, but it is exactly what the
    NEXT step starts from: with F held (F_d = I) the trial Hencky strain the wrapper reports ("strain") is
    1/2 ln b^e_committed.  So after 10 plastic K2 steps and an OBLIQUE-axis rigid rotation (0.2 rad about (1, 2, 3):
    all three shear components of eps^e are non-zero) one held step gives
        exp(2 * strain_hold)  ==  expm(2 * eps^e_rot)
    where eps^e_rot is the inner's own `elasticStrain` (engineering Voigt) at the rotation step, compared as 3x3
    tensors at 1e-12 relative to max|b|; and, as the direct sensitivity statement, the held Hencky strain equals
    eps^e_rot component by component (engineering Voigt, shear included) to 1e-11 of max|eps|.
    alpha0 = 0 included: there the old inv(D0) recovery already had the right deviator but a wrong isotropic part
    (K = -p / kappa_hat is nonlinear), which does not cancel in b^e itself.
    Kills: MX1 (route disabled), MX2 (shear not doubled: 50 % of every shear component), the committed instead of the
    trial eps^e, an eps^e taken before the plastic correction, a transposed Voigt order."""
    R = rot_axis((1.0, 2.0, 3.0), 0.2)
    P, r = norsand_rotation(alpha0, R)
    eps_rot, hencky_hold = r["rot"]["eps_e"], r["hold"]["hencky"]
    assert np.abs(eps_rot[3:]).min() > 1.0e-7, f"all three shear components must be non-trivial: {eps_rot}"
    b_exp = expm(2.0 * _tensor_from_eng(eps_rot))
    b_obs = expm(2.0 * _tensor_from_eng(hencky_hold))
    e_b = float(np.abs(b_obs - b_exp).max() / np.abs(b_exp).max())
    e_eps = float(np.abs(hencky_hold - eps_rot).max() / np.abs(eps_rot).max())
    print(f"\nalpha0 = {alpha0}: |exp(2 hencky_hold) - exp(2 eps_e)| / max|b| = {e_b:.2e};  max|hencky_hold - eps_e| / "
          f"max|eps_e| = {e_eps:.2e};  shear eps_e = {eps_rot[3:]}")
    assert e_b <= B_IDENT_TOL, e_b
    assert e_eps <= 1.0e-11, e_eps


# ----------------------------------------------------------------------------------------------
# the fallback: a NON-provider inner under LogStrain keeps the v1 inv(D0) recovery
# ----------------------------------------------------------------------------------------------
K_FB, G_FB = 1500.0, 700.0                                  # the K, G of tests/test_finite_strain_L1_analytical.py
E_FB = 9.0 * K_FB * G_FB / (3.0 * K_FB + G_FB)
NU_FB = (3.0 * K_FB - 2.0 * G_FB) / (2.0 * (3.0 * K_FB + G_FB))
FB_STRETCH = F1_STRETCH
FB_TOL = 1.0e-12


def _elastic_inner(tag):
    ops.nDMaterial("ElasticIsotropic", tag, E_FB, NU_FB)


def _j2_inner(tag):
    ops.nDMaterial("LadrunoJ2", tag, K_FB, G_FB, "-iso", "voce", 10.0, 0.0, 1.0, 60.0, "-kin", 0)


def _hencky_tau(lnF):
    """Closed-form Kirchhoff stress of the elastic Hencky law for principal log stretches (the diagonal path)."""
    tr = float(np.sum(lnF))
    return np.diag(K_FB * tr + 2.0 * G_FB * (lnF - tr / 3.0))


def _compliance_eps(tau):
    """eps^e = C : tau of the isotropic linear law, as a tensor: tr(tau)/(9K) I + dev(tau)/(2G)."""
    tr = float(np.trace(tau))
    return tr / (9.0 * K_FB) * I3 + (tau - tr / 3.0 * I3) / (2.0 * G_FB)


def _eng(t):
    return np.array([t[0, 0], t[1, 1], t[2, 2], 2.0 * t[0, 1], 2.0 * t[1, 2], 2.0 * t[0, 2]])


@pytest.mark.parametrize("inner", ("ElasticIsotropic", "LadrunoJ2"))
def test_non_provider_inner_keeps_the_d0_inversion_fallback(inner):
    """ElasticIsotropic and LadrunoJ2 (plastic) are NOT providers: the wrapper's elastic strain is the unchanged v1
    recovery eps^e = inv(D0) : tau, exact for a linear-elastic inner.  Closed forms (no oracle, no wrapper output):
      (a) elastic inner, coaxial path: Cauchy sigma = tau / J with the Hencky law tau = K tr(ln F) I + 2G dev(ln F)
          after the 10 steps (1e-12 relative); the J2 inner is checked to be PLASTIC instead (its Mises stress is
          below the elastic line by > 10 %), i.e. the recovery is exercised where tau != D : feed;
      (b) both inners, after a 0.2 rad rigid rotation: sigma' = R sigma R^T (1e-10, relative to max|sigma|), and a
          held step leaves the stress unchanged (1e-10);
      (c) both inners: the Hencky strain of the held step, which is 1/2 ln b^e_committed, equals C : tau of the
          rotation step with C the isotropic compliance of (K, G) and tau = J sigma (1e-12 of max|eps|): the
          committed b^e is exp(2 C tau).
    Bit-identity with the pre-G2 LogStrain on these inners is measured separately (mutation_gate.md, MX0): there is
    no recorded reference file because a bit pattern of libm results is not portable across CPUs.
    Kills: an edit that routes a non-provider through the provider branch or changes the fallback arithmetic; a
    fallback that is skipped when the provider cast fails."""
    define = _elastic_inner if inner == "ElasticIsotropic" else _j2_inner
    R = rot_z(0.2)
    r = drive_rotation(lambda Fs: G.build_finite_inner(define, Fs), FB_STRETCH, R, nsteps=10, nhold=1)
    b, a, h = r["before"], r["rot"], r["hold"]
    lnF = 10.0 * np.log(FB_STRETCH)
    J = math.exp(float(lnF.sum()))
    scale = float(np.abs(b["sig"]).max())
    tau_el = _hencky_tau(lnF)
    if inner == "ElasticIsotropic":
        e_a = float(np.abs(b["sig"] - tau_el / J).max() / np.abs(tau_el / J).max())
        print(f"\n[{inner}] (a) Cauchy vs tau/J of the Hencky closed form: {e_a:.2e}")
        assert e_a <= FB_TOL, e_a
    else:
        def dev(t):
            return t - np.trace(t) / 3.0 * I3
        q_el = math.sqrt(1.5) * np.linalg.norm(dev(tau_el))
        q_j2 = math.sqrt(1.5) * np.linalg.norm(dev(b["sig"] * J))
        print(f"\n[{inner}] (a) Mises(tau) {q_j2:.3f} vs elastic line {q_el:.3f}")
        assert q_j2 < 0.9 * q_el, "the J2 path must be plastic"
    exp = R @ b["sig"] @ R.T
    e_rot = float(np.abs(a["sig"] - exp).max() / scale)
    e_hold = float(np.abs(h["sig"] - exp).max() / scale)
    eps_cl = _eng(_compliance_eps(J * a["sig"]))
    e_c = float(np.abs(h["hencky"] - eps_cl).max() / np.abs(eps_cl).max())
    print(f"[{inner}] (b) rotation {e_rot:.2e}  hold {e_hold:.2e};  (c) hencky_hold vs C:tau {e_c:.2e}")
    assert float(np.abs(exp - b["sig"]).max()) > 1.0, "the rotation must move the stress"
    assert e_rot <= ROT_TOL and e_hold <= ROT_TOL, (e_rot, e_hold)
    assert e_c <= FB_TOL, e_c
