"""WP-144 round 3b, gate G2 (Zone A): the HAR energy option, the p' floor and the unified pi_i0 rule of `nDMaterial LadrunoNorSand`,
THROUGH THE SHELL (openseespy only; no oracle, scipy or sympy is imported -- every expected value is a closed form of the equation
sheet 144a written in this file, or a number PRINTED in the sheet).

Scope.  The material-level shell contract is `tests/test_ladruno_norsand.py`, the element-level contract `tests/test_ladruno_norsand_element.py`;
the Zone B parity of the shell with the O2 oracle on HAR / floor / pi_i0 paths is `norsand_oracle/g2/test_g2_shell_parity.py` and
`test_g2_logstrain.py` (HAR objectivity).  This file is what is NEW in round 3b at the shell:

  (1) THE ECHO (the material says which energy, which p_ref, which floor it runs), and EVERY PARSER REFUSAL of sheet 2.4 with its code:
      201 (a BA06 constant with -energy HAR), 202 (a HAR constant with BA06), 203 (both stiffness forms), 204 (one of a pair), 206 (DM04
      range), 207 (-pi0 with -pi0_auto), the kernel codes 21 (k <= 0), 22 (g <= 0), 23 (n outside [0, 1)), 24 (p_a <= 0), 25 (p_min < 0),
      26 (the gated (S.56) smooth-cap refusal; planar / none never), a bad -energy name, HAR without -p_a.  Positive controls beside each.
  (2) HAR CLOSED FORMS: K1.1h (p(eps_v) and K at the printed TIMs values), K1.11 (constant-volume shear: eta = 3 g eps_s exactly, printed p and q),
      the DM04 mapping (-G0 -nu -e_ref -> the printed g, k), the gate-table row eigenvalues of the 6 x 6 tangent (lambda_min / K_iso at
      eta = 0, 1.331, 2.1: 0.8551, 0.9750, 1.0980), the symmetric elastic tangent.
  (3) THE FLOOR: K1.12 (BA06: p = -0.5, d eps^f_v = kappa ln 2, W_f = 3.465735903e-3, E_f = 2.5e-3, tangent 2 mu0 (I - delta delta/3), counters, the
      idempotent second increment, leaving the floor), K1.13 (HAR in and out of the domain: no refusal, the same eps_v,f, 2G(p_min) tangent;
      -pmin 0 refuses; the BA06 control), K1.14 (HAR under shear: eps_s unchanged, q_f = 14.8362051, delta:C = 0, the (S.51a) tangent), the
      `floor`, `floorEnergy`, `floorInit` responses and stepInfo[9:11], the projected initial state, the unified pi_i0 (K1.15 -50.995881 /
      -46.4758 / -71.554175) behind -pi0_auto.
  (4) ELEMENTS: a zero-free-DOF cube driven past the HAR domain edge floors and is never refused (and is cut by the forwarding element with
      -pmin 0); the assembled free-DOF stiffness under HAR against the central FD of the resisting force (the element file's own gate on the HAR
      law: D12 != 0 and q/eps_s != D22 live in a real Newton iteration); a database round trip carries the floor counters.

Every test names, in a `KILLS:` line, the mutant (a silent regression of the shell) that makes it fail.  The numbers are the sheet's; none is
harvested from the shell.  THESE TESTS NEED THE ROUND-3B BUILD (`-energy`, `-pmin`, `-pi0_auto`, the `floor*` responses): against a pre-round-3b
binary they fail at construction.
"""
import math
import os
import re
import tempfile

import numpy as np
import pytest

import test_ladruno_norsand as T0               # engine binding + the cube / prescription helpers of the shell file
import test_ladruno_norsand_element as TE        # the FD machinery of the element file

ops = T0.ops
_OpsErr = T0._OpsErr
pytestmark = [pytest.mark.zone_a]

# ---- the TIMs HAR constants of the sheet (2.3 gate table, 13.1h) -------------------------------------------------------
KH, GH, NH, PA = 1889.48104361, 807.80387674, 0.5, 101.0
KN = KH * (1.0 - NH)
EDGE = 1.0 / KN                                  # 1.058491699e-3: p = 0 (13.1h)
PMIN_H = 5.0e-3 * PA                             # 0.505 kPa (1.3, 9.7)
SQ23, SQ32 = math.sqrt(2.0 / 3.0), math.sqrt(1.5)

_HAR_BASE = dict(energy="HAR", k=KH, g=GH, n=NH, p_a=PA, M=1.3309, N=0.4, N_bar=0.2, rho=0.71, rho_bar=0.71, chi=-3.5, h=280.0,
                 csl="fork", e0=0.83, lambda_c=0.027, xi=0.45, v0=1.70, pi0=-5000.0)
_K2 = dict(p0=-100.0, kappa_hat=0.01, mu0=5400.0, M=1.2, N=0.4, N_bar=0.2, rho=0.7, rho_bar=0.8, chi=-3.5, h=280.0, csl="paper",
           lambda_tilde=0.0135, v_c0=1.81, v0=1.70, pi0=-5000.0)
KAPPA, MU0 = 0.01, 5400.0


def flags(d):
    """{name: value} -> ['-name', value, ...]; a list value is splatted (-sigma0 s11 s22 s33 s12 s23 s13); None drops the flag."""
    out = []
    for k, v in d.items():
        if v is None:
            continue
        out.append("-" + k)
        if isinstance(v, (list, tuple)):
            out += list(v)
        elif v is not True:
            out.append(v)
    return out


def har(**over):
    d = dict(_HAR_BASE)
    d.update(over)
    return d


def k2(**over):
    d = dict(_K2)
    d.update(over)
    return d


def make(d, tag=1):
    ops.nDMaterial("LadrunoNorSand", tag, *flags(d))


def fresh():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)


def nd(e, tag=1):
    ops.NDTest("SetStrain", tag, *[float(x) for x in e])


def resp(name, tag=1):
    return list(ops.NDTest("GetResponse", tag, name))


def stress(tag=1):
    return list(ops.NDTest("GetStress", tag))


def tangent(tag=1):
    return np.array(ops.NDTest("GetTangentStiffness", tag), dtype=float).reshape(6, 6)


def pq_of(sig):
    """(p, q) of a Voigt stress (compression negative)."""
    s = np.array([[sig[0], sig[3], sig[5]], [sig[3], sig[1], sig[4]], [sig[5], sig[4], sig[2]]])
    w = np.linalg.eigvalsh(s)
    p = float(w.mean())
    return p, SQ32 * float(np.linalg.norm(w - p))


# ---- closed forms of the sheet (S.5h), (S.5h'), (S.5h''), (S.50), (S.4h): plain python -----------------------------------
def har_pq(ev, es):
    est = EDGE - ev
    u = math.sqrt(est * est + 3.0 * GH * es * es / KN)
    w = (KN * u) ** (NH / (1.0 - NH))
    return -PA * KN * est * w, 3.0 * GH * PA * es * w


def har_D(p, q):
    vp2 = p * p + KN * q * q / (3.0 * GH)
    fac = PA * (math.sqrt(vp2) / PA) ** NH
    Z = vp2 / (p * p)
    return (KH * fac * (1 - NH + NH / Z), NH * KH * p * q * fac / vp2, (3.0 * GH / (1 - NH)) * fac * (1 - NH / Z), 3.0 * GH * fac)


def har_inverse(p, q):
    vp = math.sqrt(p * p + KN * q * q / (3.0 * GH))
    return (1.0 / KN) * (1.0 - (abs(p) / PA) ** (1 - NH) * (abs(p) / vp) ** NH), q / (3.0 * GH * PA * (vp / PA) ** NH), vp


def har_psi(ev, es):
    est = EDGE - ev
    u = math.sqrt(est * est + 3.0 * GH * es * es / KN)
    return PA / (KH * (2 - NH)) * (KN * u) ** ((2 - NH) / (1 - NH))


def har_floor(es, pmin=PMIN_H):
    """(S.50), n = 1/2 closed: x, eps_v,f, q_f, eps'_f."""
    a, b = 3.0 * KN * GH * es * es, (pmin / PA) ** 2
    x = 0.5 * (a + math.sqrt(a * a + 4.0 * b))
    q_f = 3.0 * GH * PA * es * x ** NH
    return dict(x=x, ev_f=EDGE - (pmin / PA) / (KN * x ** NH), q_f=q_f,
                epsp=NH * pmin * q_f / ((1 - NH) * PA * PA * x * x + NH * pmin * pmin))


def ae_har(p, q, nh):
    """(S.3) a^e in the principal basis (HAR D)."""
    one = np.ones(3)
    D11, D12, D22, ratio = har_D(p, q)
    return (D11 * np.outer(one, one) + SQ23 * D12 * (np.outer(one, nh) + np.outer(nh, one)) + (2.0 / 3.0) * D22 * np.outer(nh, nh)
            + (2.0 * ratio / 3.0) * (np.eye(3) - np.outer(one, one) / 3.0 - np.outer(nh, nh)))


def shear_free_block(c_diag, c_off, c_shear):
    """6 x 6 Voigt tangent (engineering shear columns halved) of C = c_diag delta_ij + c_off (i != j) with shear entry c_shear."""
    T = np.full((6, 6), 0.0)
    for i in range(3):
        for j in range(3):
            T[i, j] = c_diag if i == j else c_off
    for i in range(3, 6):
        T[i, i] = c_shear
    return T


def text(capfd):
    c = capfd.readouterr()
    return c.out + c.err


def expect_refused(capfd, code, d, tag=1):
    """The command raises, prints 'REFUSED (code <code>)' (or the given message), and registers nothing."""
    capfd.readouterr()
    fresh()
    with pytest.raises(_OpsErr):
        make(d, tag)
    out = text(capfd)
    if code is not None:
        assert f"REFUSED (code {code})" in out, (code, out)
    return out


# ======================================================================================================================
# 1. THE ECHO
# ======================================================================================================================
def test_echo_har_says_the_energy_p_ref_floor_and_the_axis_poisson_ratio(capfd):
    """The echo of a HAR material (sheet 2.4 table: 'parameters ... Print echo'): the energy line with k, g, n, p_a and the axis Poisson ratio
    nu = (3k - 2g)/(6k + 2g) = 0.3129 (13.1h / 2.3), 'energy option: HAR, p_ref=101 (= p_a)', the floor ON at p_min = 0.505 = 0.005 p_ref (the default),
    and the unified-rule note when -pi0_auto is given.  Printed with 6 digits: compared to 5e-6.
    KILLS: an echo that still prints the BA06 constants for a HAR deck (the user thinks alpha0 is live), a p_ref of |p0| = 101 taken from the wrong source,
    a default floor of 0.5 (the BA06 number) under HAR (M-F8), the floor silently off."""
    capfd.readouterr()
    fresh()
    make(har(pi0=None, **{"pi0_auto": True}))
    out = text(capfd)
    m = re.search(r"HAR energy.*?k=(\S+) g=(\S+) n=(\S+) p_a=(\S+)\s", out)
    assert m, out
    k, g, n, pa = (float(x) for x in m.groups())
    assert abs(k - KH) <= 5e-6 * KH and abs(g - GH) <= 5e-6 * GH and n == NH and pa == PA, (k, g, n, pa)
    assert "p0=" not in out.split("energy option")[0].split("elastic")[-1], out                    # no BA06 constants in the HAR echo
    mnu = re.search(r"nu = (\S+?)[;,]", out)
    assert mnu and abs(float(mnu.group(1)) - 0.3129) <= 5e-5, out
    assert re.search(r"energy option: HAR, p_ref=101 \(= p_a\)", out), out
    m = re.search(r"p' floor: ON, p_min=(\S+) \((\S+) p_ref\)", out)
    assert m and abs(float(m.group(1)) - 0.505) <= 5e-6 and abs(float(m.group(2)) - 0.005) <= 5e-8, out
    assert "unified rule" in out and "pi_i0 from the unified rule" in out, out


def test_echo_ba06_default_and_floor_off(capfd):
    """BA06 (the default, paper mode): 'energy option: BA06 (default; paper mode), p_ref=100 (= |p0|)', default p_min = 0.5 = 0.005 p_ref; with -pmin 0 the
    echo says 'p' floor: OFF (-pmin 0)' and that a trial with p >= 0 is a REFUSAL as before round 3; -pmin 2 echoes p_min=2.
    KILLS: a p_ref taken as p_a under BA06, a default floor that follows the wrong reference (M-F8), an echo that hides the floor being off."""
    capfd.readouterr()
    fresh()
    make(k2())
    out = text(capfd)
    assert "energy option: BA06 (default; paper mode), p_ref=100 (= |p0|)" in out, out
    m = re.search(r"p' floor: ON, p_min=(\S+) \((\S+) p_ref\)", out)
    assert m and abs(float(m.group(1)) - 0.5) <= 5e-6 and abs(float(m.group(2)) - 0.005) <= 5e-8, out
    make(k2(pmin=0.0), tag=2)
    assert "p' floor: OFF (-pmin 0)" in text(capfd)
    make(k2(pmin=2.0), tag=3)
    m = re.search(r"p' floor: ON, p_min=(\S+)", text(capfd))
    assert m and float(m.group(1)) == 2.0


# ======================================================================================================================
# 2. EVERY PARSER REFUSAL (sheet 2.4) -- codes, with positive controls
# ======================================================================================================================
@pytest.mark.parametrize("flag,val", [("p0", -100.0), ("kappa_hat", 0.01), ("eps_v0", 0.0), ("mu0", 5400.0), ("alpha0", 0.0)])
def test_har_refuses_each_ba06_constant_code_201(capfd, flag, val):
    """Sheet 2.4: any of -p0 -kappa_hat -eps_v0 -mu0 -alpha0 GIVEN with -energy HAR is REFUSED (code 201), never ignored -- also at its default value
    (alpha0 = 0, eps_v0 = 0): HAR replaces the alpha0 coupling entirely.  The refusal names the flag.  Positive control: the same deck without it.
    KILLS: a parser that silently ignores the BA06 constants under HAR (the user believes alpha0 = 5 is active)."""
    out = expect_refused(capfd, 201, har(**{flag: val}))
    assert "-" + flag in out, out
    fresh()
    make(har())
    text(capfd)


@pytest.mark.parametrize("flag,val", [("k", 1000.0), ("g", 500.0), ("n", 0.5), ("G0", 264.32), ("nu", 0.3129), ("e_ref", 0.6944)])
def test_ba06_refuses_each_har_constant_code_202(capfd, flag, val):
    """Sheet 2.4: -k -g -n (and the DM04 inputs -G0 -nu -e_ref) given under the BA06 energy (the default) are REFUSED (code 202).  `-p_a` is NOT in this
    list (the fork-CSL parameter too): a BA06 deck with -csl fork -p_a 101 and a BA06 paper deck with -p_a (only a NOTE) are accepted.
    KILLS: HAR constants silently dropped on a BA06 deck (the user believes the energy is HAR); -p_a wrongly in the refusal list."""
    expect_refused(capfd, 202, k2(**{flag: val}))
    capfd.readouterr()
    fresh()
    make(k2(csl="fork", lambda_tilde=None, v_c0=None, e0=0.83, lambda_c=0.027, xi=0.45, p_a=101.0))
    make(k2(p_a=101.0), tag=2)
    assert "-p_a is not used (BA06 energy with the paper CSL)" in text(capfd)


def test_har_stiffness_forms_codes_203_204_206_207(capfd):
    """Sheet 2.3 / 2.4 (shell): both stiffness forms (-k -g AND -G0 -nu) -> 203; one of a pair (-k alone, -g alone, -G0 alone, -nu alone, -e_ref alone) -> 204;
    an out-of-range DM04 input (G0 <= 0, nu >= 1/2, nu <= -1, e_ref <= -1) -> 206; -pi0 together with -pi0_auto -> 207.
    KILLS: a parser that lets the two forms silently fight (one wins), accepts half a pair (k without g -> g garbage), lets the DM04 mapping divide by
    (1 - 2 nu) = 0, or lets two pi_i0 rules both pass."""
    dm = dict(G0=264.32, nu=0.3129, e_ref=0.6944)
    expect_refused(capfd, 203, har(**dm))
    expect_refused(capfd, 204, har(g=None))
    expect_refused(capfd, 204, har(k=None))
    base = har(k=None, g=None)
    expect_refused(capfd, 204, dict(base, G0=264.32))
    expect_refused(capfd, 204, dict(base, nu=0.3129))
    expect_refused(capfd, 204, dict(base, e_ref=0.6944))
    for bad in (dict(G0=-1.0), dict(G0=0.0), dict(nu=0.5), dict(nu=-1.0), dict(e_ref=-1.0), dict(e_ref=-1.5)):
        expect_refused(capfd, 206, dict(base, **dict(dm, **bad)))
    expect_refused(capfd, 207, har(pi0=-5000.0, **{"pi0_auto": True}))


@pytest.mark.parametrize("label,over,code", [
    ("k = 0", dict(k=0.0), 21), ("k < 0", dict(k=-5.0), 21), ("g = 0", dict(g=0.0), 22), ("g < 0", dict(g=-5.0), 22),
    ("n = 1", dict(n=1.0), 23), ("n > 1", dict(n=1.5), 23), ("n < 0", dict(n=-0.1), 23),
    ("p_a = 0", dict(p_a=0.0), 24), ("p_a < 0", dict(p_a=-101.0), 24), ("pmin < 0", dict(pmin=-0.1), 25)])
def test_har_kernel_range_refusals_codes_21_to_25(capfd, label, over, code):
    """Sheet 2.4 hard refusals through the shell with the kernel's codes: k <= 0 (21), g <= 0 (22), n < 0 or n >= 1 (23: n = 1 is HAR05 eq 47-48, another
    closed form, not shipped), p_a <= 0 (24), p_min < 0 (25).  n = 0 is accepted (HAR05 eq 22) and so is -pmin 0 (positive controls).
    KILLS: a missing range check (n = 1 divides by zero at k(1-n); p_a <= 0 inverts p), a wrong code."""
    expect_refused(capfd, code, har(**over))
    fresh()
    make(har(n=0.0))
    make(har(pmin=0.0), tag=2)
    text(capfd)


def test_pmin_negative_refused_code_25_under_ba06_too(capfd):
    """-pmin < 0 is refused under BA06 as well (sheet 9.7: 'p_min < 0 refused'; 0 switches the floor off).  KILLS: the check living only in the HAR branch."""
    expect_refused(capfd, 25, k2(pmin=-1.0))


def test_bad_energy_name_and_har_without_p_a_are_errors_with_their_message(capfd):
    """`-energy foo` (anything but BA06 | HAR, case-insensitive) is an error that names the choices; `-energy HAR` without -p_a says that -p_a is needed
    (the reference pressure, one flag shared with the fork CSL); lower-case `har` is accepted.  KILLS: a typo silently read as BA06, HAR with a hidden
    default p_a (a unit-blind 101.325)."""
    capfd.readouterr()
    fresh()
    with pytest.raises(_OpsErr):
        make(har(energy="HAR05"))
    assert "-energy wants BA06|HAR" in text(capfd)
    with pytest.raises(_OpsErr):
        make(har(p_a=None))
    assert "-energy HAR needs -p_a" in text(capfd)
    make(har(energy="har"), tag=2)
    text(capfd)


def test_smooth_cap_scan_gate_code_26_only_for_smooth(capfd):
    """(S.56) round 3b A3 through the shell: a smooth cap with c1 = 0.05, c2 = 0.06 (W_ramp = 0.0061 < 10 PI_SCAN_REL) is REFUSED (code 26); c2 = 0.07
    (W_ramp = 0.0122) and the K2 default (0.05, 0.15) are accepted; the planar cap (c1 = c2) and no cap are NEVER refused (W_ramp = 0 there: an ungated
    check would reject every planar and no-cap model).  The message carries the W_ramp number.
    KILLS: no refusal, an ungated refusal (planar / none refused), the round-3 inverted W_ramp (-0.0645: the K2 default refused), the factor 10 changed."""
    out = expect_refused(capfd, 26, k2(cap="smooth", c1=0.05, c2=0.06))
    assert "W_ramp=0.0061" in out, out
    fresh()
    make(k2(cap="smooth", c1=0.05, c2=0.07))
    make(k2(cap="smooth", c1=0.05, c2=0.15), tag=2)
    make(k2(cap="planar", c1=0.10), tag=3)
    make(k2(cap="planar", c1=0.01, c2=0.01), tag=4)
    make(k2(cap="none"), tag=5)
    text(capfd)


# ======================================================================================================================
# 3. HAR CLOSED FORMS THROUGH THE SHELL
# ======================================================================================================================
def tims_material(**over):
    fresh()
    make(har(**over))


def test_k1_1h_isotropic_compression_printed_values_and_bulk_modulus():
    """13.1h: p(eps_v) = -p_a [1 - k(1-n) eps_v]^(1/(1-n)), K = k p_a (|p|/p_a)^n.  Printed (TIMs): eps_v = -1e-3 -> p = -381.983585 kPa, K = 371129.585 kPa;
    eps_v = +5e-4 -> p = -28.117707 kPa (the shell compared to the printed digits: 5e-9 relative).  The initial state is the HAR origin (sigma0 = -p_a, eps^e = 0,
    p0 := -p_a).  Closed form at 6 more strains to 1e-10.  The floor is inert (|p| >= 28 >> 0.505): its counters stay zero.
    KILLS: HAR replaced by BA06 (p = p0 exp(-eps_v/kappa): hundreds of kPa off), n hard-coded, p_a left at 101.325, K from the wrong formula."""
    tims_material()
    assert abs(stress()[0] + PA) <= 1e-12 * PA
    for ev, p_printed in ((-1e-3, -381.983585), (5e-4, -28.117707)):
        nd([ev / 3.0] * 3 + [0, 0, 0])
        p = sum(stress()[:3]) / 3.0
        assert abs(p - p_printed) <= 5e-9 * abs(p_printed) + 5e-7, (ev, p)   # + half a unit of the last printed digit (6 decimals)
    nd([-1e-3 / 3.0] * 3 + [0, 0, 0])
    T = tangent()
    K = float(T[:3, :3].sum()) / 9.0
    assert abs(K - 371129.585) <= 5e-9 * 371129.585, K
    for ev in (-2e-3, -5e-4, 0.0, 1e-4, 3e-4, 7e-4, 9e-4):
        nd([ev / 3.0] * 3 + [0, 0, 0])
        p = sum(stress()[:3]) / 3.0
        assert abs(p - har_pq(ev, 0.0)[0]) <= 1e-10 * abs(p), (ev, p)
    assert resp("floor") == [0.0, 0.0, 0.0, 0.0, 0.0]


def test_k1_11_constant_volume_shear_p_grows_and_eta_is_3_g_eps_s():
    """13.11: at eps_v = 0 from eps^e = 0, p(eps_s) = -p_a [1 + 3 g k(1-n) eps_s^2]^(n/(2(1-n))), q = 3 g p_a eps_s [.]^same, eta = q/|p| = 3 g eps_s EXACTLY (the
    stress-induced |p| growth).  Printed: eps_s = 8.66546968e-4 -> eta = 2.1, p = -166.548671, q = 349.752210; eps_s = 1e-3 -> p = -183.183351, q = 443.928664, eta =
    2.42341163.  Strain: tr = 0, e = eps_s sqrt(3/2) (1,1,-2)/sqrt(6) (TXC meridian); stress sigma_11 - sigma_33 = q.  M = 50 keeps the path elastic
    (eta = 2.4 is above the usual surface).  The BA06 control (K2, alpha0 = 0): p = -100 constant on the same strain.
    KILLS: HAR replaced by BA06 (p constant), D12 / the |p| growth dropped, q computed without the w factor."""
    tims_material(M=50.0)
    for es, eta, p_pr, q_pr in ((8.66546968e-4, 2.1, -166.548671, 349.752210), (1e-3, 2.42341163, -183.183351, 443.928664)):
        nd([0.5 * es, 0.5 * es, -es, 0, 0, 0])
        p, q = pq_of(stress())
        assert abs(p - p_pr) <= 3e-9 * abs(p_pr) + 5e-7 and abs(q - q_pr) <= 5e-9 * q_pr + 5e-7, (p, q)
        assert abs(q / abs(p) - 3.0 * GH * es) <= 1e-10 * 3.0 * GH * es and abs(q / abs(p) - eta) <= 3e-9 * eta
        pc, qc = har_pq(0.0, es)
        assert abs(p - pc) <= 1e-10 * abs(pc) and abs(q - qc) <= 1e-10 * qc
    fresh()
    make(k2(M=50.0))
    nd([0.5e-3, 0.5e-3, -1e-3, 0, 0, 0])
    p, q = pq_of(stress())
    assert abs(p + 100.0) <= 1e-9 and abs(q - 3.0 * MU0 * 1e-3) <= 1e-9 * q                          # the discriminating control


def test_dm04_mapping_gives_the_printed_g_and_k_and_the_same_response(capfd):
    """Sheet 2.3 / 15: -G0 264.32 -nu 0.3129 -e_ref 0.6944 (TIMs) gives g = G0 (2.97 - e)^2/(1 + e) = 807.80387674 and k = g 2(1+nu)/(3(1-2nu)) = 1889.48104361
    (the printed constants, 1e-9), announced in the echo; e_ref defaults to v0 - 1 (-v0 1.6944); the response at eps_v = -1e-3 equals the explicit -k -g material
    (381.983585 printed) to 1e-8.  KILLS: a DM04 mapping with the wrong f(e), (1 + nu) <-> (1 - nu), a lost factor 3, e_ref ignored."""
    capfd.readouterr()
    fresh()
    make(har(k=None, g=None, G0=264.32, nu=0.3129, e_ref=0.6944))
    out = text(capfd)
    m = re.search(r"-> g=(\S+) k=(\S+)", out)
    assert m, out
    assert abs(float(m.group(1)) - GH) <= 5e-6 * GH and abs(float(m.group(2)) - KH) <= 5e-6 * KH
    nd([-1e-3 / 3.0] * 3 + [0, 0, 0])
    p = sum(stress()[:3]) / 3.0
    assert abs(p + 381.983585) <= 1e-8 * 381.983585, p
    fresh()
    make(har(k=None, g=None, G0=264.32, nu=0.3129, v0=1.6944))
    assert "(= v0 - 1)" in text(capfd)
    nd([-1e-3 / 3.0] * 3 + [0, 0, 0])
    assert abs(sum(stress()[:3]) / 3.0 + 381.983585) <= 1e-8 * 381.983585


GATE_ROWS = [(0.0, "0.8551"), (1.331, "0.9750"), (2.1, "1.0980")]       # sheet 2.3 gate table, lambda_min(6-D)/K_iso


@pytest.mark.parametrize("p", [-10.0, -3.5])
@pytest.mark.parametrize("eta,printed", GATE_ROWS)
def test_har_elastic_tangent_gate_table_eigenvalues_and_symmetry(p, eta, printed):
    """Sheet 2.3 'Gate table' through the shell: lambda_min of the 6-D Mandel Hessian (D T D, D = diag(1,1,1,sqrt2,sqrt2,sqrt2) from the engineering-shear
    Voigt tangent) over K_iso(p) = k p_a (|p|/p_a)^n at eta = 0, 1.331, 2.1: 0.8551, 0.9750, 1.0980 (printed, half a unit in the last digit); the same for
    p = -10 and -3.5 (a function of eta alone); the six eigenvalues are {eig of [[3 D11, sqrt2 D12],[sqrt2 D12, 2 D22/3]]} u {2 G x 4}, the tangent is symmetric
    (an elastic state of a conservative law: D12 != 0 enters symmetrically).  The strain is the inverse map (S.5h'') of (p, eta), M = 50 keeps it elastic.
    KILLS: AB06 eq 64 (the (2/3)(D22 - q/eps_s) n n term dropped: the rows at eta >= 1 move by 1e-2), HAR replaced by BA06, D12 = 0, a Voigt shear-factor slip
    (the Mandel conversion would not give 4 equal shear eigenvalues)."""
    tims_material(M=50.0)
    q = eta * abs(p)
    ev, es, vp = har_inverse(p, q)
    nh = np.array([-1.0, 0.2, 0.8])
    nh = nh / np.linalg.norm(nh)
    e = ev / 3.0 + SQ32 * es * nh
    nd(list(e) + [0, 0, 0])
    p_s, q_s = pq_of(stress())
    assert abs(p_s - p) <= 1e-9 * abs(p) and abs(q_s - q) <= 1e-9 * max(q, 1e-3 * abs(p))
    T = tangent()
    D = np.diag([1, 1, 1, math.sqrt(2), math.sqrt(2), math.sqrt(2)])
    M6 = D @ T @ D
    assert np.abs(M6 - M6.T).max() <= 1e-9 * np.abs(M6).max()
    eig = np.linalg.eigvalsh(0.5 * (M6 + M6.T))
    Kiso = KH * PA * (abs(p) / PA) ** NH
    dec = len(printed.split(".")[1])
    assert abs(eig[0] / Kiso - float(printed)) <= 0.5 * 10.0 ** (-dec) * 1.0001 + 1e-9, (eig[0] / Kiso, printed)
    D11, D12, D22, _ = har_D(p, q)
    two_G = 2.0 * GH * PA * (vp / PA) ** NH
    blk = np.array([[3 * D11, math.sqrt(2) * D12], [math.sqrt(2) * D12, 2 * D22 / 3]])
    want = np.sort(np.concatenate([np.linalg.eigvalsh(blk), [two_G] * 4]))
    assert np.abs(eig - want).max() <= 1e-9 * want.max(), (eig, want)


# ======================================================================================================================
# 4. THE FLOOR THROUGH THE SHELL
# ======================================================================================================================
DEV_K112 = KAPPA * math.log(4.0)             # p: -1 -> -0.25 = -p_min/2 (13.12)


def ba06_floor_material(**over):
    fresh()
    make(k2(sigma0=[-1.0, -1.0, -1.0, 0.0, 0.0, 0.0], **over))


def test_k1_12_ba06_floor_stress_counters_energy_and_tangent():
    """13.12 (K2 set, alpha0 = 0, default p_min = 0.5): from p = -1 a trial at p^tr = -p_min/2 = -0.25 (eps_v,tr = 0.0599146455) is projected: committed p = -0.5,
    eps_v,f = -kappa ln(p_min/|p0|) = 0.0529831737 (printed), d eps^f_v = kappa ln 2 = 6.931471806e-3, W_f = p_min d eps^f_v = 3.465735903e-3 (printed), E_f =
    kappa (p_min - |p^tr|) = 2.5e-3.  Responses: floor = [1, 1, 0, d eps^f_v, W_f]; floorEnergy = [E_f, E_f, W_f]; floorInit = 0; stepInfo[9:11] = (1, 0); refusal
    all zero (counted, NEVER refused).  pi_i = -5000 and v = v0 exp(tr d eps) (the RAW trace) unchanged by the projection (M-F7).  Tangent: 2 mu0 (I - delta
    delta/3) = 7200 / -3600 normal block, shear entry mu0 = 5400 (engineering), delta:C = 0 (no bulk stiffness faked).  A second increment with the same strain is
    NOT a floor event (idempotent: counters stay, stepInfo (0, 0), at_floor stays 1).
    KILLS: M-F1 (no floor: p = -0.25), M-F2 (not counted), M-F7 (pi_i or v moved), a bulk regularisation in the tangent, a re-projection chatter, slots of
    floor / floorEnergy swapped."""
    ba06_floor_material()
    nd([DEV_K112 / 3.0] * 3 + [0, 0, 0])
    assert max(abs(x + 0.5) for x in stress()[:3]) <= 1e-12 and max(abs(x) for x in stress()[3:]) <= 1e-12
    dfv, Wf = KAPPA * math.log(2.0), 0.5 * KAPPA * math.log(2.0)
    assert abs(dfv - 6.931471806e-3) <= 6e-13 and abs(Wf - 3.465735903e-3) <= 6e-13                   # the printed digits
    fl = resp("floor")
    assert fl[0] == 1.0 and fl[1] == 1.0 and fl[2] == 0.0 and abs(fl[3] - dfv) <= 1e-13 and abs(fl[4] - Wf) <= 1e-14, fl
    fe = resp("floorEnergy")
    assert abs(fe[0] - 2.5e-3) <= 1e-13 and abs(fe[1] - 2.5e-3) <= 1e-13 and abs(fe[2] - Wf) <= 1e-14 and 0.0 <= fe[0] <= fe[2], fe
    assert resp("floorInit") == [0.0]
    info = resp("stepInfo")
    assert info[0] == 0.0 and info[9] == 1.0 and info[10] == 0.0, info
    assert resp("refusal") == [0.0, 0.0, 0.0, 0.0, 0.0]
    ee = resp("elasticStrain")
    assert abs(sum(ee[:3]) - 0.0529831737) <= 6e-11, ee                                                # printed eps_v,f
    st = resp("state")
    assert st[0] == -5000.0 and abs(st[2] - 1.70 * math.exp(DEV_K112)) <= 1e-13 and st[3] == 1.70, st
    T = tangent()
    want = shear_free_block(2.0 * MU0 * 2.0 / 3.0, -2.0 * MU0 / 3.0, MU0)
    assert np.abs(T - want).max() <= 1e-9 * 2.0 * MU0, np.abs(T - want).max()
    assert np.abs(T[:3, :].sum(axis=0)).max() <= 1e-9 * 2.0 * MU0, "delta:C != 0 at a floored state"
    # idempotence through the shell: commit, same total strain again
    ops.NDTest("CommitState", 1)
    nd([DEV_K112 / 3.0] * 3 + [0, 0, 0])
    assert resp("floor")[:3] == [1.0, 1.0, 0.0] and resp("stepInfo")[9:11] == [0.0, 0.0]
    assert max(abs(x + 0.5) for x in stress()[:3]) <= 1e-12
    # unloading by half the expansion leaves the floor with the full bulk stiffness: p = -1 (the start), at_floor 0, K = -p/kappa
    nd([DEV_K112 / 6.0] * 3 + [0, 0, 0])
    p = sum(stress()[:3]) / 3.0
    assert abs(p + 1.0) <= 1e-10 and resp("floor")[0] == 0.0 and resp("stepInfo")[9:11] == [0.0, 0.0]
    K = float(tangent()[:3, :3].sum()) / 9.0
    assert abs(K - 1.0 / KAPPA) <= 1e-8 / KAPPA


def test_floor_off_pmin_zero_is_an_ordinary_elastic_step_and_counts_nothing():
    """-pmin 0 (sheet 9.7, 1.3): the very increment of K1.12 is an elastic step to p = p0 exp(-eps_v/kappa) = -0.25; floor = zeros; no event.
    KILLS: a hard-coded floor that ignores -pmin 0, a counter moved with the floor off."""
    ba06_floor_material(pmin=0.0)
    nd([DEV_K112 / 3.0] * 3 + [0, 0, 0])
    assert abs(sum(stress()[:3]) / 3.0 + 0.25) <= 1e-10
    assert resp("floor") == [0.0, 0.0, 0.0, 0.0, 0.0] and resp("stepInfo")[9:11] == [0.0, 0.0] and resp("floorEnergy") == [0.0, 0.0, 0.0]


@pytest.mark.parametrize("dev,dfv_printed,out_of_domain", [(1.0e-4, 6.9522805e-5, False), (1.1e-4, 7.9522805e-5, True)])
def test_k1_13_har_floor_in_and_out_of_the_domain_is_counted_not_refused(dev, dfv_printed, out_of_domain):
    """13.13, TIMs HAR, p_min = 0.505: an isotropic state at p = -1 (eps_v = 9.53167838e-4, the inverse map) with d eps_v = +1e-4 is in the domain (p^tr =
    -2.56e-3) and floors to eps_v,f = 9.83645033e-4 with d eps^f_v = 6.9522805e-5, W_f = 3.5109017e-5; with +1.1e-4 the trial (eps_v = 1.0632e-3 > 1.0585e-3 =
    the domain edge) is OUT of the domain and floors to the SAME eps_v,f with d eps^f_v = 7.9522805e-5: NO REFUSAL (refusal response zero).  Committed p =
    -0.505, q = 0; floor = [1, 1, 0, d eps^f_v, W_f]; floorEnergy[0] = Psi(eps_f) - Psi(eps_pre) (S.4h) in the domain, the bound W_f out of it (S.52);
    tangent 2 G(p_min) (I - delta delta/3), G(p_min) = g p_a (p_min/p_a)^(1/2) = 5769.156 (normal 2/3 2G = 7692.2, off -3846.1, shear G).
    KILLS: M-F1 (refusal or pass-through), M-F5 (trial floor skipped: the out-of-domain trial refuses), M-F4 (the BA06 inverse: eps_v,f = 0.0530), a bulk
    stiffness added, the shear stiffness regularised."""
    tims_material(sigma0=[-1.0, -1.0, -1.0, 0.0, 0.0, 0.0])
    assert abs(sum(resp("elasticStrain")[:3]) - 9.53167838e-4) <= 1e-10
    assert (EDGE < 9.53167838e-4 + dev) == out_of_domain
    nd([dev / 3.0] * 3 + [0, 0, 0])
    assert resp("refusal") == [0.0, 0.0, 0.0, 0.0, 0.0], "a floor event was REFUSED"
    s = stress()
    assert max(abs(x + PMIN_H) for x in s[:3]) <= 1e-12 and max(abs(x) for x in s[3:]) <= 1e-12
    assert abs(sum(resp("elasticStrain")[:3]) - 9.83645033e-4) <= 6e-12
    fl = resp("floor")
    assert fl[:3] == [1.0, 1.0, 0.0] and abs(fl[3] - dfv_printed) <= 6e-13 and abs(fl[4] - PMIN_H * fl[3]) <= 1e-15, fl
    if not out_of_domain:
        assert abs(fl[4] - 3.5109017e-5) <= 6e-13
    fe = resp("floorEnergy")
    if out_of_domain:
        assert abs(fe[0] - fl[4]) <= 1e-16, "out-of-domain trial: only the bound W_f is stated (S.52)"
    else:
        want = har_psi(har_floor(0.0)["ev_f"], 0.0) - har_psi(9.53167838e-4 + dev, 0.0)
        assert abs(fe[0] - want) <= 1e-8 * abs(want) and 0.0 <= fe[0] <= fl[4], (fe, want)
    G = GH * PA * math.sqrt(PMIN_H / PA)
    assert abs(G - 5769.156) <= 6e-4
    T = tangent()
    assert np.abs(T - shear_free_block(2.0 * G * 2.0 / 3.0, -2.0 * G / 3.0, G)).max() <= 1e-8 * 2.0 * G
    assert np.abs(T[:3, :].sum(axis=0)).max() <= 1e-8 * 2.0 * G


def test_k1_13_controls_ba06_has_no_event_and_pmin_zero_refuses_the_out_of_domain_trial():
    """13.13 controls: (a) under BA06 (K2, p_min 0.5) the same two isotropic increments are ordinary elastic steps (M-F4: no floor); (b) with -pmin 0 the HAR
    out-of-domain trial (+1.1e-4) is the pre-round-3 REFUSAL: the refusal response is non-zero, finest 'p or pi_i not negative' (6), the committed stress
    is the initial one, the in-domain +1e-4 is accepted; the same deck with the default floor is not refused.
    KILLS: a floor computed with the wrong energy (events where there are none), -pmin 0 not restoring the refusal, a refused trial that moves the state."""
    for dev in (1.0e-4, 1.1e-4):
        ba06_floor_material()
        nd([dev / 3.0] * 3 + [0, 0, 0])
        assert abs(sum(stress()[:3]) / 3.0 + math.exp(-dev / KAPPA)) <= 1e-9 and resp("floor")[:3] == [0.0, 0.0, 0.0]
    tims_material(sigma0=[-1.0, -1.0, -1.0, 0.0, 0.0, 0.0], pmin=0.0)
    s0 = stress()
    nd([1.0e-4 / 3.0] * 3 + [0, 0, 0])
    assert resp("refusal")[0] == 0.0
    ops.NDTest("CommitState", 1)
    tims_material(sigma0=[-1.0, -1.0, -1.0, 0.0, 0.0, 0.0], pmin=0.0)
    nd([1.1e-4 / 3.0] * 3 + [0, 0, 0])
    r = resp("refusal")
    assert r[0] != 0.0 and r[1] >= 1.0 and int(r[3]) == 6, r
    assert stress() == s0, "a refused trial moved the stress"


def test_k1_14_har_floor_under_shear_printed_values_and_tangent():
    """13.14, TIMs HAR (M = 50 keeps the path elastic: eta_f = 29 lies above any physical surface), p_min = 0.505: from the origin (sigma0 = -p_a) a single trial to
    (eps_v, eps_s) = (1.050e-3, 2e-4) (in the domain, p^tr = -0.245 > -p_min) floors to x = varpi_f/p_a = 0.09185198375, eps_v,f = 1.04102893e-3, q_f = 14.8362051 kPa
    (eta_f = 29.38); the elastic shear strain is UNCHANGED by Pi_f (M-F6: the projection keeps the deviatoric strain, q is not kept); committed p = -0.505.  Tangent:
    C_f = a^e(eps_f) Phi with Phi = I - 1/3 + (1/3) eps'_f sqrt(2/3) 1 n_hat^T, eps'_f = 0.08679793 (printed), a^e from (S.3) with the HAR D, compared on the normal
    block to 5e-8 (the printed digits of eps'_f); delta:C_f = 0 to 1e-12.
    KILLS: M-F3b (eps'_f dropped: delta:C = sqrt(6) D12 n_hat != 0), M-F6 (q kept), M-F4 (BA06 floor value), a^e built without the t2 / t4 terms of (S.3)."""
    tims_material(M=50.0)
    nh = np.array([-1.0, 0.2, 0.8])
    nh = nh / np.linalg.norm(nh)
    es = 2e-4
    e = 1.050e-3 / 3.0 + SQ32 * es * nh
    nd(list(e) + [0, 0, 0])
    assert resp("refusal")[0] == 0.0
    p, q = pq_of(stress())
    assert abs(p + PMIN_H) <= 1e-12 and abs(q - 14.8362051) <= 6e-8 and abs(q / PMIN_H - 29.38) <= 5e-3, (p, q)
    ee = np.array(resp("elasticStrain")[:3])
    assert abs(ee.sum() - 1.04102893e-3) <= 6e-12
    es_out = SQ23 * float(np.linalg.norm(ee - ee.mean()))
    assert abs(es_out - es) <= 1e-12 * es + 1e-15, "eps_s changed by the projection (M-F6)"
    vp = math.sqrt(p * p + KN * q * q / (3.0 * GH))
    assert abs(vp / PA - 0.09185198375) <= 6e-12
    fl = resp("floor")
    assert fl[:3] == [1.0, 1.0, 0.0]
    one = np.ones(3)
    Phi = np.eye(3) - 1.0 / 3.0 + (1.0 / 3.0) * 0.08679793 * SQ23 * np.outer(one, nh)
    Cexp = ae_har(p, q, nh) @ Phi
    T = tangent()
    assert np.abs(T[:3, :3] - Cexp).max() <= 5e-8 * np.abs(Cexp).max(), np.abs(T[:3, :3] - Cexp).max() / np.abs(Cexp).max()
    assert np.abs(T[:3, :].sum(axis=0)).max() <= 1e-12 * np.abs(T).max(), "delta:C_f != 0 at a floored HAR state"
    Cmut = ae_har(p, q, nh) @ (np.eye(3) - 1.0 / 3.0)
    dm = Cmut.sum(axis=0)
    assert np.abs(dm).max() / np.abs(Cmut).max() >= 1e-2, "M-F3b would go unseen"


def test_initial_state_is_projected_counted_and_the_stress_replaced_ba06_and_har(capfd):
    """9.7 `initialState`: sigma0 below the floor is replaced by sigma(eps^e_f) and COUNTED: BA06 K2 sigma0 = isotropic -0.2 -> p = -0.5, floorInit = 1, floor = [1,
    0, 0, kappa ln(0.5/0.2) = 9.162907319e-3, W_f = 0.5 x that], the echo says 'PROJECTED by the p' floor (n_f_init = 1'; HAR TIMs sigma0 (p, q) = (-0.3, 0.5)
    -> p = -0.505, eps_s unchanged, q_f = q (varpi_f/varpi)^n.  An initial state above the floor (-5) is not projected (floorInit 0).
    KILLS: an initial state left below p_min, an uncounted projection, q kept (the BA06 rule) under HAR, a refusal of a deck whose first Gauss points sit near the surface."""
    capfd.readouterr()
    fresh()
    make(k2(sigma0=[-0.2, -0.2, -0.2, 0.0, 0.0, 0.0]))
    assert "PROJECTED by the p' floor (n_f_init = 1" in text(capfd)
    assert max(abs(x + 0.5) for x in stress()[:3]) <= 1e-12
    assert resp("floorInit") == [1.0]
    want = KAPPA * math.log(0.5 / 0.2)
    fl = resp("floor")
    assert abs(fl[3] - want) <= 1e-12 and abs(fl[4] - 0.5 * want) <= 1e-12 and fl[0] == 1.0
    nh = np.array([-1.0, 0.2, 0.8])
    nh = nh / np.linalg.norm(nh)
    p0_, q0_ = -0.3, 0.5
    sig = list(p0_ + SQ23 * q0_ * nh) + [0.0, 0.0, 0.0]
    fresh()
    make(har(sigma0=sig, M=50.0))
    ev_i, es_i, vp_i = har_inverse(p0_, q0_)
    f = har_floor(es_i)
    q_f = q0_ * (f["x"] * PA / vp_i) ** NH
    p, q = pq_of(stress())
    assert abs(p + PMIN_H) <= 1e-12 and abs(q - q_f) <= 1e-10 * q_f, (q, q_f)
    ee = np.array(resp("elasticStrain")[:3])
    assert abs(SQ23 * np.linalg.norm(ee - ee.mean()) - es_i) <= 1e-13 and resp("floorInit") == [1.0]
    assert abs(resp("floor")[3] - (ev_i - f["ev_f"])) <= 1e-12
    fresh()
    make(har(sigma0=[-5.0, -5.0, -5.0, 0.0, 0.0, 0.0]))
    assert resp("floorInit") == [0.0] and abs(stress()[0] + 5.0) <= 1e-12


@pytest.mark.parametrize("cap,sigma_eta,printed", [("smooth", 0.0, -50.995881), ("none", 0.0, -46.475800), ("planar", 0.0, -50.995881),
                                                    ("smooth", 0.75, -71.554175), ("none", 1.2, -100.0)])
def test_pi0_auto_is_the_unified_rule_k1_15(capfd, cap, sigma_eta, printed):
    """13.15, K2 set, p_init = -100, -pi0_auto: eta* = max(eta_init, c2 M): isotropic -> -50.995881 (smooth cap c2 = 0.15 and the planar cap c1 = c2 = 0.15),
    -46.475800 without a cap (the apex through p_init); a compression-meridian start with eta_init = 0.75 -> -71.554175; eta_init = M -> -100.  Read from the
    state response pi_i (before any step), 1e-8 relative (8 printed digits); the echo says it came from the unified rule.  An explicit -pi0 overrides.
    KILLS: the pre-round-3 apex behind -pi0_auto under a capped model (-46.4758), c2 mapped to c1, eta* = min / sum, the planar cap's c2 taken as 0."""
    capfd.readouterr()
    fresh()
    sig = [-100.0] * 3 + [0.0] * 3
    if sigma_eta > 0.0:
        q = sigma_eta * 100.0
        nh = np.array([1.0, 1.0, -2.0]) / math.sqrt(6.0)                   # theta = pi/3 (zeta = 1): eta_init = q / |p|
        sig = list(-100.0 + SQ23 * q * nh) + [0.0] * 3
    over = dict(sigma0=sig, pi0=None, **{"pi0_auto": True})
    if cap == "smooth":
        over.update(cap="smooth", c1=0.05, c2=0.15)
    elif cap == "planar":
        over.update(cap="planar", c1=0.15, c2=0.15)
    make(k2(**over))
    out = text(capfd)
    st = resp("state")
    assert abs(st[0] - printed) <= 1e-8 * abs(printed) + 6e-7, (st[0], printed)
    assert "pi_i0 from the unified rule" in out, out
    if sigma_eta == 0.0:                                       # an explicit -pi0 overrides the rule (inside the surface for the isotropic start)
        fresh()
        make(k2(**dict(over, pi0=-80.0, pi0_auto=None)))
        assert resp("state")[0] == -80.0


# ======================================================================================================================
# 5. ELEMENTS
# ======================================================================================================================
def har_cube(elem, e_final, n, **over):
    """A zero-free-DOF cube (every node DOF prescribed, u = f(t) E x) of the HAR material, LoadControl 1/n; returns nothing."""
    fresh()
    make(har(**over))
    ids = T0._cube(elem, 1, 1)
    ops.timeSeries("Linear", 1)
    T0._prescribe(ids, e_final, 1, 1)
    T0._analysis(1.0 / n)


def test_cube_driven_past_the_har_domain_edge_floors_and_is_never_refused_but_pmin_zero_is_cut():
    """A zero-free-DOF LadrunoBrick (it FORWARDS a material refusal) of the HAR TIMs material, sigma0 = isotropic -2 kPa (eps_v = 9.095e-4, 1.49e-4 below the domain
    edge), driven by an isotropic expansion of 6e-4 in 10 steps: with the default floor every step converges (analyze == 0), the committed stress ends at p =
    -0.505 (1e-10), and the Gauss point's counters are exactly the closed form: d eps^f_v total = 6e-4 - (eps_v,f - eps_v0) (every step after the first floor
    event is absorbed whole), W_f = p_min eps^f_v, floorEnergy cumulative in [0, W_f]; refusal zero.  With -pmin 0 the same deck is CUT by the forwarding element
    (analyze != 0) once the trial leaves the domain, and the point counts refusals.  The stdBrick (which discards the code) floors identically.
    KILLS: M-F1 (the trial leaves the domain and the step is refused), a floor that does not reach the element (the stress at the Gauss points is the raw
    HAR extrapolation), counters lost between the shell and the Gauss-point clones."""
    ev0 = har_inverse(-2.0, 0.0)[0]
    total = 6e-4
    want_efv = total - (har_floor(0.0)["ev_f"] - ev0)
    for elem in ("LadrunoBrick", "stdBrick"):
        har_cube(elem, [total / 3.0] * 3 + [0, 0, 0], 10, sigma0=[-2.0, -2.0, -2.0, 0.0, 0.0, 0.0])
        for k in range(10):
            assert ops.analyze(1) == 0, f"{elem}: step {k + 1} refused / failed"
        s = T0._gp("stress")
        assert max(abs(x + PMIN_H) for x in s[:3]) <= 1e-10 and max(abs(x) for x in s[3:]) <= 1e-10, s
        fl = T0._gp("floor")
        assert fl[0] == 1.0 and fl[1] >= 1.0 and abs(fl[3] - want_efv) <= 1e-12, (fl, want_efv)
        assert abs(fl[4] - PMIN_H * fl[3]) <= 1e-15 + 1e-12 * fl[4]
        fe = T0._gp("floorEnergy")
        assert -1e-18 <= fe[1] <= fl[4] * (1.0 + 1e-12) and fe[2] == fl[4], (fe, fl)
        assert T0._gp("refusal") == [0.0, 0.0, 0.0, 0.0, 0.0]
    har_cube("LadrunoBrick", [total / 3.0] * 3 + [0, 0, 0], 10, sigma0=[-2.0, -2.0, -2.0, 0.0, 0.0, 0.0], pmin=0.0)
    rcs = [ops.analyze(1) for _ in range(10)]
    assert any(rc != 0 for rc in rcs), "pmin = 0: the trial outside the domain must be REFUSED (the step cut)"
    assert T0._gp("refusal")[1] >= 1.0


def _deck_brick_har(Q, T):
    """The element file's drained-triaxial load-control deck (confining -100 kPa, axial ramp Q, shear couple T) on the TIMs HAR material:
    sigma0 = isotropic -100, pi_i0 = -80, psi_i0 = -0.05 (fork CSL: v0 = 1.7556886184)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    make(har(v0=1.0 + 0.83 - 0.05 - 0.027 * (80.0 / PA) ** 0.45, pi0=-80.0, sigma0=[-100.0, -100.0, -100.0, 0.0, 0.0, 0.0]))
    for i, c in enumerate(TE._BRICK_XYZ, 1):
        ops.node(i, float(c[0]), float(c[1]), float(c[2]))
    ops.fix(1, 1, 1, 1)
    ops.fix(2, 0, 1, 1)
    ops.fix(4, 0, 0, 1)
    ops.element("stdBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
    ops.timeSeries("Constant", 2)
    ops.pattern("Plain", 2, 2)
    for n in TE._XF1:
        ops.load(n, -25.0, 0.0, 0.0)
    for n in TE._XF0:
        ops.load(n, +25.0, 0.0, 0.0)
    for n in TE._YF1:
        ops.load(n, 0.0, -25.0, 0.0)
    for n in TE._YF0:
        ops.load(n, 0.0, +25.0, 0.0)
    for n in TE._ZF1:
        ops.load(n, 0.0, 0.0, -25.0)
    for n in TE._ZF0:
        ops.load(n, 0.0, 0.0, +25.0)
    ops.timeSeries("Linear", 3)
    ops.pattern("Plain", 3, 3)
    for n in TE._ZF1:
        ops.load(n, 0.0, 0.0, -Q / 4.0)
    for n in TE._ZF0:
        ops.load(n, 0.0, 0.0, +Q / 4.0)
    if T:
        for n in TE._ZF1:
            ops.load(n, T / 4.0, 0.0, 0.0)
        for n in TE._ZF0:
            ops.load(n, -T / 4.0, 0.0, 0.0)
        for n in TE._XF1:
            ops.load(n, 0.0, 0.0, T / 4.0)
        for n in TE._XF0:
            ops.load(n, 0.0, 0.0, -T / 4.0)
    TE._solver_stack()
    return {"nodes": list(range(1, 9)), "ndf": 3, "eles": [1], "ngp": 8}


@pytest.mark.parametrize("case", ["txc_1", "shear_1"])
def test_har_brick_assembled_tangent_vs_fd(case):
    """The element file's gate (assembled stdBrick stiffness, free equations, FullGeneral == central FD of the resisting force, best-h per-column <= 1e-6, plastic at
    every Gauss point, not refused, the SAME substep count at every FD point; the symmetrised and the committed-state tangents must FAIL the same gate) on the HAR
    law: D12 != 0 and q/eps_s != D22 are live in a real Newton iteration.  Drained TXC to q = 160 kPa in one 100 kPa step from q = 60 (yield at q ~ 80), and the
    non-coaxial family with tau_xz = 30 kPa (the same load history as the BA06 cases txc_1 / shear_1, m = 1 substep).  The expected tangent IS the FD of the
    element's own force (the definition of a consistent tangent); the tolerance is the plan sec.5.1 bound.
    KILLS: the HAR D12 / t4 terms of (S.3) missing from the shell tangent, an energy-dependent term dropped in the assembled K, a symmetrised tangent."""
    c = TE._BRICK_CASES[case]
    res = TE._fd_case(_deck_brick_har, c["Q"], c["T"], c["hist"], c["big"], c["m"])
    assert res["N"] == 18
    TE._assert_fd(res, f"stdBrick-HAR/{case}")


def test_database_roundtrip_carries_the_floor_counters_and_energy():
    """sendSelf / recvSelf through `database File`: a HAR material whose Gauss point has floored (6 of 10 expansion steps, the first floor event inside) is restored
    into a fresh skeleton with the SAME `floor` response (at_floor, n_f_tr, n_f_post, eps^f_v, W_f) and the same cumulative E_f (`floorEnergy[1:3]`), bit-exactly,
    and the continuation (4 more steps) reaches the never-saved reference (stress 1e-12, floor counters equal).
    KILLS: a wire that drops or reorders the floor block (counters reset to 0 on restore), E_f lost, a committed floor state not restored."""
    def deck():
        har_cube("stdBrick", [6e-4 / 3.0] * 3 + [0, 0, 0], 10, sigma0=[-2.0, -2.0, -2.0, 0.0, 0.0, 0.0])

    def run_to(n):
        deck()
        for _ in range(n):
            assert ops.analyze(1) == 0

    run_to(6)
    for _ in range(4):
        assert ops.analyze(1) == 0
    ref = dict(stress=T0._gp("stress"), floor=T0._gp("floor"), fe=T0._gp("floorEnergy"))
    run_to(6)
    before = dict(floor=T0._gp("floor"), fe=T0._gp("floorEnergy"))
    assert before["floor"][1] >= 1.0 and before["floor"][3] > 0.0, "the point must have floored before the save"
    with tempfile.TemporaryDirectory(prefix="ladruno_norsand_", ignore_cleanup_errors=True) as td:
        db = os.path.join(td, "norsand_floor_rt")
        ops.database("File", db)
        ops.save(1)
        deck()
        assert T0._gp("floor") != before["floor"]
        ops.database("File", db)
        ops.restore(1)
        after = dict(floor=T0._gp("floor"), fe=T0._gp("floorEnergy"))
        assert after["floor"] == before["floor"], (before["floor"], after["floor"])
        assert after["fe"][1:] == before["fe"][1:], (before["fe"], after["fe"])
        ops.wipeAnalysis()
        T0._analysis(1.0 / 10.0)
        for _ in range(4):
            assert ops.analyze(1) == 0
        cont = dict(stress=T0._gp("stress"), floor=T0._gp("floor"), fe=T0._gp("floorEnergy"))
        ops.wipe()
    assert T0._maxabs(cont["stress"], ref["stress"]) <= 1e-12 * max(abs(x) for x in ref["stress"])
    assert cont["floor"][:3] == ref["floor"][:3] and abs(cont["floor"][3] - ref["floor"][3]) <= 1e-14
