"""WP-144 G2 (Zone B), part 1: the LadrunoNorSand NDMaterial SHELL through OpenSees against O2, step by step.

WHAT IS COMPARED.  O2 (o2_algo) is the contract (backward-Euler spectral return map + chained substep tangent,
sheet 144a sections 9.1-9.6).  The G1 paths of kernel_parity (drained / undrained TXC and TXE, NONCOAXIAL, paper and
fork CSL, Willam-Warnke and Gudehus-Argyris, caps none / planar / smooth including the AMP_STOP n = 40 smooth-cap
path where substepping is normal, the planar / no-cap AMP_STOP paths that O2 refuses, the vertex / isotropic paths,
N = 0, alpha0 != 0, the neutral-shear path) are replayed through the real `nDMaterial LadrunoNorSand` in a unit
stdBrick (and a LadrunoBrick -geom linear) whose eight nodes are fully prescribed (no free DOF; g2_common.py): the
shell receives exactly O2's strain increments as engineering-shear Voigt vectors, and every quantity the shell
exposes is compared with O2 per step.  Three additional G2 paths exercise what the G1 paths never do:
all three shear components (the 13 slot, vmap index 5, is zero on every G1 path) and an off-axis start.

  quantity                      shell                                      O2
  stress (Voigt 11 22 33 12 23 13)  getStress                              sigma tensor, shear = same number
  elastic strain                elasticStrain response (gamma = 2 eps)     eps_e tensor * (1,1,1,2,2,2)
  pi_i, v, v0, eps_p_v, eps_p_s state response [pi_i, psi_i, v, v0, ..]    State fields
  psi_i, psi                    state[1], psi response                     sheet (S.22) closed form on the shell's
                                                                           own (v, pi_i, p): a formula check
  D                             D response                                 State.D (sum over sub-increments)
  stepInfo                      [refusal plastic vertex cap local pi substeps finest finest_sub]   State.flags
  substeps response             [last, number of COMMITTED steps substepped]   flags['substeps'] > 1 count
  tangent                       tangent response, 36 values row-major      c4_to_c6(O2.tangent) with the shear
                                (T[a][b] = C[a][b] w_b, w_b = 1/2 for      COLUMNS halved (the engineering-shear
                                the shear columns, rows unchanged)         convention of LadrunoNorSand.h)

GATES (expected values from O2 / the sheet, written before the shell was run; none harvested from the shell).
  * 1e-10 relative on every quantity and every step: the shell is the same algorithm as the kernel that
    kernel_parity gates at 1e-10; what is new here is the shell (Voigt <-> tensor mapping, commit / revert, response
    packing) and the OpenSees element that feeds it.  The relative measure, the natural scales and the zero floor
    are those of kernel_parity (max|a - b| / max(max|b|, 1e-6 x natural scale)).
  * Tangent: 1e-10 off the three bands of kernel_parity/test_kernel_parity.py (corner |sin 3 theta| < 1e-7 or
    near-vertex R < 1e-3 |p|; near-coalescent trial / final eigenvalues; both documented in that module's docstring),
    1e-7 inside them.  The band decision is taken on O2's states (identical to kernel_parity's `_bands`).
  * Flags: plastic / vertex / cap_active / substeps are equal on every step, local and nested iteration counts
    equal (same algorithm, same constants), the committed-substepped census equals O2's.
  * Refusal: at O2's refusal step analyze() fails (stdBrick discards the material code: the commit aborts and the
    point LATCHES; LadrunoBrick forwards it: the trial is cut, no latch); the `refusal` response then reads
    SUBSTEPS_EXHAUSTED with O2's finest reason and sub-reason; the stress is the frozen committed one.

THE OPENSEES NOISE (measured, not assumed).  The brick forms the strain from eight nodal displacements, so the
material sees O2's strains to ~1e-16 absolute (measured 1.1e-16 on the cumulative strain, test_shell_noise_floor)
rather than bitwise: an exactly axisymmetric strain (eps_11 = eps_22) reaches the shell with a 1-ulp gap.  Strain-like
quantities that pass through zero (eps_e at the apex, eps_p_v on undrained paths) therefore get an absolute
tolerance of 1e-15 (STRAIN_NATURAL) instead of kernel_parity's 1e-18; D is normalised path-wide as in
kernel_parity.  The 1-ulp gap is the "near-coalescent" regime of kernel_parity's tangent band, entered at an
exactly-zero gap too; no step needed it (measured: the off-band tangent worst case is 2e-11 on the neutral-shear
path, `tangent_coalescent` is 0 on every path), but `test_shell_noise_floor` states and bounds the input noise.

Runtime: measured in the module summary test (see the printed wall time); every case is well below 2 minutes, so
nothing here is @slow.

MUTANTS each test kills (shell mutants; the C++ is mutated by the Mutate step, this file names what it must catch):
  parity paths           shear strain not halved in setTrialStrain (NONCOAXIAL, ALLSHEAR); tangent shear columns not
                         halved or ROWS halved instead (NONCOAXIAL, ALLSHEAR); row-major / column-major swap
                         (the asymmetry witness makes this visible: ~1e-2 on the NONCOAXIAL tangent); the
                         12 <-> 23 <-> 13 index map (ALLSHEAR); commit not copying epsC / sC (every plastic
                         path: wrong deps from step 2); the v term of the tangent (vfac = v_n or v0 instead of
                         v_{n+1}, any plastic path with a volume change: 1.7e-6 / 1.8e-5 on a^ep per sheet 1.2);
                         the v update itself (linear v0 (1 + tr eps) or v_n + v0 tr instead of v_n exp(tr d_eps):
                         test_specific_volume_is_v0_exp_of_the_total_trace_strain); state-response slot order (pi_i, psi_i, v, v0, eps_p_v, eps_p_s);
                         stepInfo / substeps response order; `D` slot; elasticStrain shear factor.
  refusal paths          refusal not propagated / latch not set on the stdBrick path / wrong finest reason.
  test_reset_replay      revertToStart not restoring sC / v0 / the strain history.
  test_gauss_points      per-GP state shared through a static (history isolation).
"""
from __future__ import annotations

import functools
import math
import os
import subprocess
import time

import numpy as np
import pytest

import g2_common as G

if G.ops is None:
    pytest.skip(f"opensees.pyd not found/loadable in {G.DIST_BIN}", allow_module_level=True)
ops = G.ops

import o2_algo as O2                       # noqa: E402
from o2_algo import api as O2api           # noqa: E402
from o2_algo import kernel as K            # noqa: E402
import ns_kernel as NK                     # noqa: E402
import test_kernel_parity as KP            # noqa: E402  (cases, O2 drivers, band helpers: one source of truth)

# no_gmsh: nothing here meshes; without it tests/conftest.py's session-wide gmsh gate skips this whole file when
# collected from the repo root on a box without gmsh (the G2 harness bug)
pytestmark = [pytest.mark.zone_b, pytest.mark.no_gmsh]

GATE = KP.GATE                             # 1e-10
BAND_TANGENT_GATE = KP.BAND_TANGENT_GATE   # 1e-7
I3 = np.eye(3)
SHEAR_ENG = np.array([1.0, 1.0, 1.0, 2.0, 2.0, 2.0])

# ----------------------------------------------------------------------------------------------
# cases: the G1 set of kernel_parity + the G2-only paths
# ----------------------------------------------------------------------------------------------
E2_ALL = np.array([[0.0005, 0.006, 0.002], [0.006, 0.0005, 0.004], [0.002, 0.004, -0.002]])


def allshear_deps(n):
    """The NONCOAXIAL path with a third shear (13 = 0.002 tensor): E1 over the first 60 %, then E2_ALL."""
    def eps(s):
        return KP.E1_NC * min(s, KP.KNOT) / KP.KNOT + E2_ALL * max(s - KP.KNOT, 0.0) / (1.0 - KP.KNOT)
    t = np.arange(n + 1) / n
    return [eps(t[k + 1]) - eps(t[k]) for k in range(n)]


SIG_OFFAXIS = np.array([[-90.0, 6.0, 2.0], [6.0, -100.0, 3.0], [2.0, 3.0, -125.0]])
CASES = dict(KP.CASES)
CASES["ALLSHEAR_paper"] = KP.case("paper", "custom", KP.LOOSE, 25, deps=allshear_deps(25))
CASES["ALLSHEAR_fork_GA"] = KP.case("fork", "custom", KP.LOOSE, 25, over=dict(KP.GA, rho=0.85, rho_bar=0.9),
                                    deps=allshear_deps(25))
CASES["OFFAXIS_onsurface_start_ALLSHEAR"] = KP.case("paper", "custom", (None, None), 25, v0=1.75,
                                                    over=dict(alpha0=2.0), sigma0=SIG_OFFAXIS,
                                                    deps=allshear_deps(25))
REFUSING = ("CAP_planar_AMPSTOP_n40", "CAP_none_AMPSTOP_n40")
ELEMENTS = ("stdBrick", "LadrunoBrick")
CASE_NAMES = list(CASES)
# the cases that are the substance of the substep / chained-tangent claim
AMP_SMOOTH = ("CAP_smooth_AMPSTOP_n40", "CAP_smooth_AMPSTOP_n40_fork")


# ----------------------------------------------------------------------------------------------
# O2 side (cached; the same drivers as kernel_parity)
# ----------------------------------------------------------------------------------------------
@functools.lru_cache(maxsize=None)
def o2_run(name):
    c = CASES[name]
    P = KP.o2_params(c)
    v0 = KP.v0_of(P, c)
    pi0 = c["start"][0]
    st0 = O2.initial_state(P, c["sigma0"], v0, pi0)
    deps = KP.increments(P, c, st0)
    sts, tans = [], []
    st = st0
    for d in deps:
        st = O2.step(P, st, d)
        sts.append(st)
        tans.append(O2.tangent(P, st))
        if st.flags["refused"]:
            break
    return P, v0, st0, deps, sts, tans


def csl_psi(P, v, pressure):
    """Sheet (S.22), written out here from the sheet (not taken from O2 or the shell): psi = v - v_c(pressure)."""
    if P.csl_mode == "paper":
        return v - P.v_c0 + P.lam_tilde * math.log(-pressure)
    return (v - 1.0) - P.e0 + P.lam_c * (-pressure / P.p_a) ** P.xi


# ----------------------------------------------------------------------------------------------
# the OpenSees side
# ----------------------------------------------------------------------------------------------
def read_step():
    """Everything the shell exposes after one analyze(), at Gauss point 1."""
    return dict(stress=G.mat_response("stress"), state=G.mat_response("state"),
                eps_e=G.mat_response("elasticStrain"), D=G.mat_response("D")[0],
                T=G.mat_response("tangent").reshape(6, 6), info=G.mat_response("stepInfo"),
                subs=G.mat_response("substeps"), psi=G.mat_response("psi")[0],
                refusal=G.mat_response("refusal"), strain=G.mat_response("strain"))


def drive(name, ele):
    """Build the model and analyze one step per O2 increment. Returns (O2 bundle, list of per-step dicts, rcs)."""
    P, v0, st0, deps, sts, tans = o2_run(name)
    c = CASES[name]
    G.build_small(G.norsand_args(P, v0, st0.pi_i, c["sigma0"]), deps, ele)
    out, rcs = [], []
    for k in range(len(sts)):
        rc = ops.analyze(1)
        rcs.append(rc)
        out.append(read_step() if rc == 0 else None)
        if rc != 0:
            break
    return (P, v0, st0, deps, sts, tans), out, rcs


# ----------------------------------------------------------------------------------------------
# comparison
# ----------------------------------------------------------------------------------------------
QUANT = KP.QUANTITIES
EXTRA = ("psi_i", "psi", "v0", "eps_e_eng", "info_dec")
STRAIN_NATURAL = G.STRAIN_NATURAL


@functools.lru_cache(maxsize=None)
def compare(name, ele):
    t0 = time.time()
    (P, v0, st0, deps, sts, tans), out, rcs = drive(name, ele)
    nat = KP._natural_scales(P)
    # STRAIN-valued quantities reach the shell through the brick with ~1e-18 absolute noise (module docstring), so
    # their zero floor is 1e-6 x 10 = 1e-5 (abs tolerance 1e-15 at GATE, 10x the 1.1e-16 input noise measured by
    # test_shell_noise_floor) instead of kernel_parity's 1e-6 x kappa_hat.
    for q in ("eps_e", "eps_p_v", "eps_p_s"):
        nat[q] = STRAIN_NATURAL
    D_sh, D_o2 = [], []
    rep = dict(name=name, ele=ele, n_o2=len(sts), errs={q: 0.0 for q in KP.REPORTED + EXTRA},
               flag_mismatch=[], iter_mismatch=0, refusal=None, substepped=0, plastic=0, band_steps=0,
               coalescent_steps=0, lastsub_gap=0.0, asym=0.0, rcs=rcs, n_committed_substepped=0,
               census_mismatch=[], fail_step=None)
    errs = rep["errs"]
    ost = st0
    nsub_o2 = 0
    for k, (d, o) in enumerate(zip(deps, sts)):
        if o.flags["refused"]:
            rep["refusal"] = (k, o.flags["reason"])
            rep["fail_step"] = k if rcs[k] != 0 else None
            rep["o2_state_k"] = o
            rep["frozen_o2"] = ost
            break
        assert rcs[k] == 0, f"{name}/{ele} step {k}: analyze failed ({rcs[k]}) but O2 accepted the step"
        s = out[k]
        inf = s["info"]
        # ---- flags / counters
        for f, val in (("plastic", inf[1]), ("vertex", inf[2]), ("cap_active", inf[3]), ("substeps", inf[6])):
            if int(val) != int(o.flags[f]):
                rep["flag_mismatch"].append((k, f, o.flags[f], val))
        if int(inf[0]) != 0 or int(inf[7]) != 0 or int(inf[8]) != 0:
            rep["flag_mismatch"].append((k, "stepInfo refusal/finest", 0, tuple(inf[[0, 7, 8]])))
        if (int(inf[4]), int(inf[5])) != (o.flags["local_iters"], o.flags["pi_iters"]):
            rep["iter_mismatch"] += 1
        rep["substepped"] += o.flags["substeps"] > 1
        rep["plastic"] += bool(o.flags["plastic"])
        nsub_o2 += o.flags["substeps"] > 1
        if int(s["subs"][1]) != nsub_o2 or int(s["subs"][0]) != o.flags["substeps"]:
            rep["census_mismatch"].append((k, nsub_o2, tuple(s["subs"])))
        # ---- state
        st = s["state"]
        assert st[3] == v0, f"{name}/{ele} step {k}: v0 not carried unchanged ({st[3]!r} vs {v0!r})"
        errs["v0"] = max(errs["v0"], abs(st[3] - v0))
        for q, x, y in (("sigma", s["stress"], G.t6(o.sigma)), ("pi_i", st[0], o.pi_i), ("v", st[2], o.v),
                        ("eps_p_v", st[4], o.eps_p_v), ("eps_p_s", st[5], o.eps_p_s)):
            errs[q] = max(errs[q], KP._rel(x, y, nat[q]))
        D_sh.append(s["D"])
        D_o2.append(o.D)
        errs["eps_e"] = max(errs["eps_e"], KP._rel(s["eps_e"], G.t6(o.eps_e) * SHEAR_ENG, nat["eps_e"]))
        # ---- sheet (S.22) on the shell's own v, pi_i, p
        p_mean = float(np.sum(s["stress"][:3])) / 3.0
        errs["psi"] = max(errs["psi"], abs(s["psi"] - csl_psi(P, st[2], p_mean)))
        errs["psi_i"] = max(errs["psi_i"], abs(st[1] - csl_psi(P, st[2], st[0])))
        # ---- tangent (row-major, shear columns halved)
        Co = NK.c4_to_c6(tans[k])
        Te = Co * G.W_SHEAR
        et = float(np.abs(s["T"] - Te).max() / np.abs(Te).max())
        rep["asym"] = max(rep["asym"], float(np.abs(Te - Te.T).max() / np.abs(Te).max()))
        m = o.flags["substeps"]
        subs = [o] if m == 1 else KP._o2_subincrements(P, ost, d, [1.0 / m] * m)
        if m > 1:
            Cl = NK.c4_to_c6(O2.tangent_last_substep(P, o)) * G.W_SHEAR
            rep["lastsub_gap"] = max(rep["lastsub_gap"], float(np.abs(Cl - Te).max() / np.abs(Te).max()))
        cv, coal = KP._bands(subs)
        if m > 1 and not (cv or coal):
            errs["tangent_substepped"] = max(errs["tangent_substepped"], et)
        KP._gate_tangent(errs, rep, et, cv, coal,
                         lambda: KP._o2_1ulp(P, ost, d, Co) )
        ost = o
    # D as kernel_parity: max|D_shell - D_O2| over the path / max|D_O2| (a per-step ratio is meaningless on the
    # neutral-loading steps where D itself is ~1e-10)
    if D_o2 and max(abs(x) for x in D_o2) > 0.0:
        errs["D"] = max(abs(x - y) for x, y in zip(D_sh, D_o2)) / max(abs(x) for x in D_o2)
    else:
        errs["D"] = max([abs(x) for x in D_sh] + [0.0])
    # ---- refusal step
    if rep["refusal"] is not None:
        k = rep["refusal"][0]
        rep["rc_at_refusal"] = rcs[k] if k < len(rcs) else None
        rep["refusal_resp"] = G.mat_response("refusal")
        rep["stress_after"] = G.mat_response("stress")
        rep["state_after"] = G.mat_response("state")
    rep["wall"] = time.time() - t0
    return rep


# ----------------------------------------------------------------------------------------------
# tests
# ----------------------------------------------------------------------------------------------
def test_binary_is_not_older_than_the_material_source():
    """BUILD_GOTCHAS 4b: a stale binary silently tests the old code.  The build stamp (ladrunoBuild) must contain
    the last commit that touched the shell / kernel sources.  Skips only when git is unavailable.
    Kills: every test in this file silently running a stale dist/bin build."""
    stamp = ops.ladrunoBuild()
    assert isinstance(stamp, str) and len(stamp) >= 7, f"no build stamp: {stamp!r}"
    try:
        last = subprocess.run(["git", "log", "-1", "--format=%H", "--", "SRC/material/nD/LadrunoNorSand.cpp",
                               "SRC/material/nD/LadrunoNorSand.h", "SRC/material/nD/LadrunoNorSandKernel.h",
                               "SRC/material/nD/LadrunoNorSand3D.cpp", "SRC/material/nD/LadrunoNorSandPlaneStrain.cpp"],
                              cwd=G.REPO, capture_output=True, text=True, check=True).stdout.strip()
        anc = subprocess.run(["git", "merge-base", "--is-ancestor", last, stamp], cwd=G.REPO,
                             capture_output=True, text=True)
    except (OSError, subprocess.CalledProcessError):
        pytest.skip("git not available")
    assert anc.returncode == 0, (f"dist/bin build stamp {stamp[:12]} does not contain the last material source commit "
                                 f"{last[:12]}: rebuild with Ladruno_scripts\\build.bat before trusting any result")


@pytest.mark.parametrize("ele", ELEMENTS)
@pytest.mark.parametrize("name", CASE_NAMES)
def test_shell_matches_o2_on_path(name, ele):
    rep = compare(name, ele)
    errs = rep["errs"]
    line = "  ".join(f"{q}={errs[q]:.2e}" for q in KP.REPORTED)
    print(f"\n[{name}/{ele}] steps {rep['n_o2']} plastic {rep['plastic']} substepped {rep['substepped']} "
          f"band {rep['band_steps']} (coalescent only {rep['coalescent_steps']}) refusal {rep['refusal']} "
          f"iter-count mismatches {rep['iter_mismatch']} wall {rep['wall']:.2f}s"
          f"\n  {line}\n  psi_i={errs['psi_i']:.1e} psi={errs['psi']:.1e} v0={errs['v0']:.1e}")
    assert not rep["flag_mismatch"], f"flags differ (step, flag, O2, shell): {rep['flag_mismatch']}"
    assert rep["iter_mismatch"] == 0, f"{rep['iter_mismatch']} steps with different local/nested iteration counts"
    assert not rep["census_mismatch"], f"substeps response differs (step, O2 census, shell): {rep['census_mismatch']}"
    bad = {q: errs[q] for q in QUANT if not errs[q] <= GATE}
    assert not bad, f"shell vs O2 gate {GATE:.0e} exceeded: {bad}"
    assert errs["psi_i"] <= 1e-12 and errs["psi"] <= 1e-12, (errs["psi_i"], errs["psi"])
    for key in ("tangent_band", "tangent_coalescent"):
        assert errs[key] <= BAND_TANGENT_GATE, f"{key} {errs[key]:.3e} > {BAND_TANGENT_GATE:.0e}"


@pytest.mark.parametrize("ele", ELEMENTS)
@pytest.mark.parametrize("name", REFUSING)
def test_refusing_paths_refuse_through_the_shell(name, ele):
    """O2 refuses these (sheet 3.2 / 10.1) after 2^8 substeps.  Through OpenSees: analyze fails at O2's step, the
    `refusal` response carries SUBSTEPS_EXHAUSTED and O2's finest reason / sub-reason, the committed state is
    untouched.  stdBrick discards the material code (WP-99: the commit aborts and the point latches);
    LadrunoBrick forwards it (the trial is cut, no latch)."""
    rep = compare(name, ele)
    assert rep["refusal"] is not None, f"{name}: O2 did not refuse (the path is no longer the refusing one)"
    k, reason = rep["refusal"]
    assert rep["rc_at_refusal"] != 0, f"{name}/{ele}: O2 refuses at step {k} but analyze() succeeded"
    assert all(r == 0 for r in rep["rcs"][:k]) and len(rep["rcs"]) == k + 1
    last, nref, latched, fin, finsub = rep["refusal_resp"]
    f2, s2 = NK.parse_o2_reason(reason)
    assert int(last) == NK.REFUSAL.index("SUBSTEPS_EXHAUSTED"), (last, reason)
    assert (int(fin), int(finsub)) == (NK.REFUSAL.index(f2), NK.EVALERR.index(s2)), \
        f"finest reason differs: O2 {reason!r} -> {(f2, s2)}, shell {(fin, finsub)}"
    assert nref >= 1
    assert int(latched) == (1 if ele == "stdBrick" else 0), \
        f"{ele}: latch {latched} (stdBrick discards the code and latches; LadrunoBrick cuts at the trial)"
    # committed state untouched: the stress is that of O2's last accepted state (the frozen one)
    frozen = rep["frozen_o2"]
    P = o2_run(name)[0]
    err = KP._rel(rep["stress_after"], G.t6(frozen.sigma), KP._natural_scales(P)["sigma"])
    assert err <= GATE, f"committed stress changed by the refused step: {err:.2e}"
    assert not rep["flag_mismatch"] and rep["errs"]["sigma"] <= GATE   # every step before the refusal is parity


def test_smooth_cap_amp_path_is_substepped_and_chained():
    """O2 README: the AMP_STOP smooth-cap path at n = 40 substeps 30 of its 34 plastic increments, where the
    consistent tangent is the CHAINED one (sheet 9.6, owner decision 2026-10-01).  Through OpenSees: the same
    substep counts (compared per step above), the committed-substepped census equals O2's, and the chained tangent
    is gated off-band at 1e-10 while O2's last-sub-increment CTO is >= 0.1 away from it (a shell returning the
    last-sub CTO cannot pass)."""
    for name in AMP_SMOOTH:
        rep = compare(name, "stdBrick")
        assert rep["refusal"] is None
        assert rep["substepped"] >= 20 and rep["plastic"] >= 1, (name, rep["substepped"])
        assert rep["errs"]["tangent_substepped"] > 0.0, "no off-band substepped tangent was compared"
        assert rep["errs"]["tangent_substepped"] <= GATE, rep["errs"]["tangent_substepped"]
        assert rep["lastsub_gap"] > KP.LASTSUB_MIN_GAP, (name, rep["lastsub_gap"])
        print(f"\n[{name}] substepped {rep['substepped']}: chained tangent vs O2 {rep['errs']['tangent_substepped']:.2e} "
              f"(off-band); O2 last-sub CTO vs chain, max {rep['lastsub_gap']:.2f}")


def test_tangent_layout_is_discriminated():
    """Non-vacuity of the tangent comparison: on the NONCOAXIAL path the expected tangent is asymmetric by at
    least 1e-2 (the model's non-associativity, sheet 9.4: 1-5 %), so a row-major / column-major (transpose) slip in
    the response or the getTangent fill would exceed the 1e-10 gate by more than 1e8."""
    rep = compare("NONCOAXIAL_paper", "stdBrick")
    assert rep["asym"] > 1.0e-2, rep["asym"]


def test_gate_discriminates_mapping_mutants():
    """The 1e-10 gate must separate the shell from its plausible mapping mutants (a test of the test: the shell
    output is held fixed and the EXPECTED side is mutated the way a shell bug would mutate the actual side):
    tangent shear columns not halved, rows halved instead, tangent transposed, 12 <-> 23 shear stress swap,
    elastic shear strain not doubled.  Measured separation >= 5e-2 on both shear paths (gate 1e-10)."""
    for name in ("NONCOAXIAL_paper", "ALLSHEAR_paper"):
        (P, v0, st0, deps, sts, tans), out, rcs = drive(name, "stdBrick")
        worst = {k: 0.0 for k in ("cols_not_halved", "rows_halved", "transposed", "shear_stress_swap", "gamma_not_2")}
        for k, o in enumerate(sts):
            Co = NK.c4_to_c6(tans[k])
            Te = Co * G.W_SHEAR
            Ta = out[k]["T"]
            n = np.abs(Te).max()
            worst["cols_not_halved"] = max(worst["cols_not_halved"], np.abs(Ta - Co).max() / n)
            worst["rows_halved"] = max(worst["rows_halved"], np.abs(Ta - Co * G.W_SHEAR[:, None]).max() / n)
            worst["transposed"] = max(worst["transposed"], np.abs(Ta.T - Te).max() / n)
            so = G.t6(o.sigma)
            so[[4, 5]] = so[[5, 4]]
            worst["shear_stress_swap"] = max(worst["shear_stress_swap"],
                                             np.abs(out[k]["stress"] - so).max() / abs(P.p0))
            worst["gamma_not_2"] = max(worst["gamma_not_2"],
                                       np.abs(out[k]["eps_e"] - G.t6(o.eps_e)).max() / P.kappa_hat)
        print(f"\n[{name}] mutant separation: " + "  ".join(f"{k}={v:.1e}" for k, v in worst.items()))
        assert all(v > 1.0e-2 for v in worst.values()), worst


def test_shell_noise_floor():
    """The brick feeds the shell strains that are O2's to ~1e-18 absolute (not bitwise).  Stated and bounded so the
    1e-10 gates above cannot be satisfied by luck: (a) the strain the shell reports back (getStrain, engineering
    shear) equals the cumulative O2 strain to <= 1e-15 absolute on a path with all three shears (measured 1.1e-16),
    and (b) the stress of that whole path agrees to 1e-12 relative.
    Kills: a brick / driver change that degrades the input (and with it the meaning of every 1e-10 gate here)."""
    name = "ALLSHEAR_paper"
    (P, v0, st0, deps, sts, tans), out, rcs = drive(name, "stdBrick")
    E = np.zeros((3, 3))
    worst = 0.0
    for k, d in enumerate(deps):
        E = E + d
        want = G.t6(E) * SHEAR_ENG
        worst = max(worst, float(np.abs(out[k]["strain"] - want).max()))
    assert worst <= 1e-15, worst
    rep = compare(name, "stdBrick")
    assert rep["errs"]["sigma"] <= 1e-12, rep["errs"]["sigma"]
    print(f"\nnoise floor: reported strain vs cumulative O2 strain, max abs {worst:.2e}")


def test_gauss_points_are_isolated_and_identical():
    """All eight Gauss points of the brick carry independent copies of the material (history isolation) and, the
    deformation being uniform, identical states: stress / state / tangent of GP 1..8 agree to 1e-13 relative after a
    plastic path with substepping.  Kills: a static / shared history between Gauss points.  Stress / state agree to 1e-13; the TANGENT to 1e-11, because the eight Gauss
    points receive strains that differ at the ulp level (different B-matrix rounding) and the tangent amplifies that
    (measured 1.2e-13 between GPs).  A shared static would differ at O(1)."""
    name = "CAP_smooth_AMPSTOP_n40"
    drive(name, "stdBrick")
    ref = {r: G.mat_response(r, 1) for r in ("stress", "state", "tangent", "stepInfo")}
    for gp in range(2, 9):
        for r, v1 in ref.items():
            v = G.mat_response(r, gp)
            err = float(np.abs(v - v1).max() / max(np.abs(v1).max(), 1e-12))
            assert err <= (1e-11 if r == "tangent" else 1e-13), f"GP{gp} {r}: {err:.2e}"


def test_reset_replay_is_bitwise():
    """revertToStart (ops.reset) restores the INITIAL state (sigma0, v0, pi0, zero strain, no latch) and the replay
    of the same history reproduces the first pass.  Mutant: revertToStart leaving sC / v0 / epsC stale or v0 := v."""
    name = "CAP_smooth_AMPSTOP_n40"
    (P, v0, st0, deps, sts, tans), out, rcs = drive(name, "stdBrick")
    first = [(o["stress"].copy(), o["state"].copy(), o["T"].copy()) for o in out]
    ops.reset()
    s0 = G.mat_response("state")
    assert s0[0] == st0.pi_i and s0[2] == v0 and s0[3] == v0 and s0[4] == 0.0 and s0[5] == 0.0, s0
    assert np.abs(G.mat_response("strain")).max() == 0.0
    for k in range(len(sts)):
        assert ops.analyze(1) == 0
        s = (G.mat_response("stress"), G.mat_response("state"), G.mat_response("tangent").reshape(6, 6))
        for a, b in zip(s, first[k]):
            assert np.array_equal(a, b), f"replay differs at step {k}"


V_PATHS = ("TXC_drained_paper", "TXE_drained_fork", "NONCOAXIAL_paper", "CAP_smooth_AMPSTOP_n40",
           "ALLSHEAR_fork_GA", "OFFAXIS_onsurface_start_ALLSHEAR", "TXC_undrained_paper")
V_SIGNAL = 1.0e-6          # the superseded linear update misses v0 exp(x) by ~ v0 x^2 / 2: asserted on the paths below


@pytest.mark.parametrize("name", V_PATHS)
def test_specific_volume_is_v0_exp_of_the_total_trace_strain(name):
    """Sheet 1.2 / 13.10 (G2 owner decision 2026-10-01, exponential update): at every committed step the shell's v
    equals v0 exp(tr eps) with eps the TOTAL strain, to 1e-12 relative, from the closed form alone: tr eps is the
    cumulative sum of the driver's own increments' traces (no O2 state, no shell output on the right-hand side),
    for any increment sequence and any substepping (exp of a sum: the composed sub-increment updates equal the
    single-increment one; CAP_smooth_AMPSTOP_n40 is substepped).  Also checks O2's v against the same closed form
    (the contract) and the v0 slot (carried unchanged).  Undrained path: tr eps = 0, v = v0 (the control).
    The superseded linear update v0 (1 + x) would miss by v0 (1 + x - e^x) ~ -v0 x^2 / 2: asserted >= V_SIGNAL on a
    volume-changing path so the gate cannot pass vacuously.
    Kills: v += v0 tr d_eps (linear), v_n (1 + tr d_eps), v += v tr d_eps, v updated from the stress-point's
    elastic instead of total strain, v lost across a substep boundary, v0 := v."""
    (P, v0, st0, deps, sts, tans), out, rcs = drive(name, "stdBrick")
    x, worst_sh, worst_o2, signal = 0.0, 0.0, 0.0, 0.0
    for k, (d, o) in enumerate(zip(deps, sts)):
        if o.flags["refused"] or out[k] is None:
            break
        x += float(np.trace(d))
        want = v0 * math.exp(x)
        worst_sh = max(worst_sh, abs(out[k]["state"][2] - want) / want)
        worst_o2 = max(worst_o2, abs(o.v - want) / want)
        signal = max(signal, abs(v0 * (1.0 + x) - want) / want)
        assert out[k]["state"][3] == v0
    print(f"\n[{name}] steps {len(deps)}  tr eps_end {x:+.3e}  |v - v0 exp(tr eps)|/v: shell {worst_sh:.2e} O2 {worst_o2:.2e};"
          f"  superseded linear rule would miss by {signal:.2e}")
    assert worst_sh <= 1e-12 and worst_o2 <= 1e-12, (worst_sh, worst_o2)
    if name == "TXC_undrained_paper":
        assert abs(x) <= 1e-12 and signal <= 1e-12, (x, signal)
    else:
        assert signal >= V_SIGNAL, f"gate has no power on {name}: tr eps_end = {x:.2e}"


def test_psi_closed_forms_are_not_vacuous():
    """The (S.22) check in compare() is a closed form on the shell's own numbers: confirm it is sensitive, i.e.
    psi_i differs from psi (they use pi_i vs p) by > 1e-3 somewhere on a plastic path; otherwise the slot swap
    state[1] <-> psi response would pass.  Kills: that slot swap (and psi evaluated at pi_i / p mixed up)."""
    (P, v0, st0, deps, sts, tans), out, rcs = drive("TXC_drained_paper", "stdBrick")
    gaps = [abs(o["state"][1] - o["psi"]) for o in out]
    assert max(gaps) > 1e-3, max(gaps)


def test_g2_shell_summary():
    """Max error per quantity over every path and both elements (the G2 shell parity numbers)."""
    reps = {(n, e): compare(n, e) for n in CASE_NAMES for e in ELEMENTS}
    worst = {q: max(r["errs"][q] for r in reps.values()) for q in KP.REPORTED + EXTRA}
    print("\nshell (OpenSees) vs O2, max relative error per quantity over %d paths x %d elements (gate %.0e; band %.0e):"
          % (len(CASE_NAMES), len(ELEMENTS), GATE, BAND_TANGENT_GATE))
    for q in KP.REPORTED + ("psi_i", "psi"):
        arg = max(reps, key=lambda kk: reps[kk]["errs"][q])
        print(f"  {q:>19s}: {worst[q]:.3e}   (worst {arg[0]}/{arg[1]})")
    tot = {k: sum(r[k] for r in reps.values()) for k in ("n_o2", "plastic", "substepped", "band_steps", "iter_mismatch")}
    print("  steps compared %d, plastic %d, substepped %d, band %d, iter-count mismatches %d, total wall %.1fs"
          % (tot["n_o2"], tot["plastic"], tot["substepped"], tot["band_steps"], tot["iter_mismatch"],
             sum(r["wall"] for r in reps.values())))
    assert all(worst[q] <= GATE for q in KP.QUANTITIES)
    assert worst["tangent_band"] <= BAND_TANGENT_GATE and worst["tangent_coalescent"] <= BAND_TANGENT_GATE
