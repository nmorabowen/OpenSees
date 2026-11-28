"""WP-263: dual-number tangent oracle for the LadrunoRCConcrete return map (header-only, g++).

Compiles tests/_testbed/rc_dual_gpp.cpp (it includes only SRC/material/nD/LadrunoRCKernel.h,
which pulls in LadrunoDual.h) and runs 18 strain paths, committed history threaded step to step.
At every trial state it compares three tangents: the shipped mode 0 (the secant W_B C0 plus the
beta / tension-stiffening / interlock cross-terms), the new mode 3 (-dualTangent: one Dual<6> pass
of returnMap3DT) and a central finite difference of the map itself at nine steps h.

Asserts:

  * mode 3's stress and every committed-history field equal mode 0's BIT FOR BIT (the Dual run's
    value part is the production arithmetic);
  * the dual tangent equals the FD of the map: <= 1e-8 relative (Frobenius) on every smooth state
    class -- cracking, crushing, beta strut, tension-stiffened, interlock capped and cyclic,
    IMPL-EX (the extrapolated map), spatial add-on plane, plane stress (the sigma_33 condensation
    with its own Newton, FD of the CONDENSED map) -- with an FD error that falls as h^2 (measured
    worst 2.8e-10); and on equal principal strains (the degenerate class, a kink of max-eigenvalue
    functions) an FD error that falls monotonically, at least first order, to <= 1e-6 (measured
    3.4e-7 at h = 1e-10): the dual is the h -> 0 limit of the central FD there;
  * the break-the-gate legs: mode 0 differs from the dual by >= 1e-2 on every loading step of the
    projector-rotation states (it omits dP/deps, the damage rate, dp1/deps, dv_ci,max/de_n), while
    on the elastic path it equals the dual to round-off; and a SYMMETRIZED dual tangent fails the
    FD by >= 1e-3 wherever damage grows (the damaged tangent is unsymmetric).

It prints the omitted-term table (||C_dual - C_mode0|| / ||C_dual|| per state class) and the
per-call cost of mode 3 vs mode 0.

The second test is the WP-263 bit-identity gate for the double path: the same file compiled with
-DRC_DUAL_DOUBLE_ONLY touches only the pre-WP double API (modes 0/1/2, the sigma_33 condensation,
driveBeamFiber) and prints an FNV-1a digest of every stress, tangent and history field; built
against the current header and against the WP-263 base header (`git show PRE_SHA:...`), the two
outputs must be byte-identical. It skips when that git object is unavailable (a shallow clone).
A later WP that changes the RC double path ON PURPOSE re-pins PRE_SHA to its own base.
"""
import collections
import os
import shutil
import subprocess
import tempfile

import pytest

from _cxx import cxx_runtime_flags  # WP-199: static GCC runtime on Windows

pytestmark = [pytest.mark.zone_a, pytest.mark.t1]

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.normpath(os.path.join(_HERE, ".."))
_SRC_ND = os.path.join(_REPO, "SRC", "material", "nD")
_SRC = os.path.join(_HERE, "_testbed", "rc_dual_gpp.cpp")
PRE_SHA = "f183b9cd1f9e63c139997aa191b6bafa9e4660c8"   # ladruno at the WP-263 cut (#1157 merge)

SCENARIOS = {
    "elastic": "elastic",
    "crack_uniaxial": "cracking", "crack_shear": "cracking", "crackednu_beta": "cracking",
    "crush": "crushing", "beta_tc": "beta",
    "tensstiff": "tension-stiffened", "tensstiff_rot": "tension-stiffened",
    "interlock_capped": "interlock-capped",
    "interlock_cyclic": "interlock-cyclic", "xcrack_dsfm": "interlock-cyclic",
    "degenerate_tension": "degenerate", "degenerate_compression": "degenerate",
    "implex_crack_shear": "implex", "implex_ts_beta": "implex",
    "spatial_ts_beta": "spatial",
    "ps_crack_shear": "plane-stress", "ps_interlock_cyclic": "plane-stress",
}
FD_TOL = 1e-8          # smooth classes (measured worst 2.8e-10)
FD_TOL_KINK = 1e-6     # degenerate class at h = 1e-10 (measured 3.4e-7)
OMIT_FLOOR = 1e-2      # mode 0 must miss the dual by at least this on the rotation states
SYM_FLOOR = 1e-3       # the symmetrized dual must miss the FD by at least this where damage grows
ROTATION = ("crack_uniaxial", "crack_shear", "tensstiff_rot", "implex_crack_shear",
            "spatial_ts_beta", "ps_crack_shear")


def _gpp():
    gpp = shutil.which("g++") or shutil.which("c++") or shutil.which("clang++")
    if gpp is None:
        pytest.skip("no C++ compiler available")
    return gpp


def _build_run(gpp, td, name, extra, includes):
    exe = os.path.join(td, name + ".exe")
    inc = [a for d in includes for a in ("-I", d)]
    cc = subprocess.run([gpp, "-O2", "-std=c++17", *extra, *inc, _SRC, "-o", exe, *cxx_runtime_flags()],
                        stdin=subprocess.DEVNULL, capture_output=True, text=True)
    assert cc.returncode == 0, f"g++ failed:\n{cc.stderr}"
    rr = subprocess.run([exe], stdin=subprocess.DEVNULL, capture_output=True, text=True)
    assert rr.returncode == 0, f"run failed:\n{rr.stderr}"
    return rr.stdout


@pytest.fixture(scope="module")
def oracle():
    gpp = _gpp()
    with tempfile.TemporaryDirectory() as td:
        out = _build_run(gpp, td, "rc_dual", [], [_SRC_ND])
    rows, cls, timing, cur = collections.OrderedDict(), {}, {}, None
    for line in out.splitlines():
        p = line.split()
        if not p:
            continue
        if p[0] == "SCENARIO":
            cur = p[1]
            cls[cur] = p[2]
            rows[cur] = []
        elif p[0] == "STEP":
            kv = {k: float(v) for k, v in zip(p[2::2], p[3::2])}
            kv["step"] = int(p[1])
            rows[cur].append(kv)
        elif p[0] == "FDH":
            rows[cur][-1]["fdh"] = [float(x) for x in p[1:]]
        elif p[0] == "TIMING" and p[1] == "mode":
            timing[int(p[2])] = float(p[3])
    return rows, cls, timing


def test_every_path_ran_and_reached_its_state(oracle):
    rows, cls, _ = oracle
    assert cls == SCENARIOS
    assert all(len(r) >= 3 for r in rows.values())
    mx = lambda c, k: max(r[k] for r in rows[c])  # noqa: E731
    mn = lambda c, k: min(r[k] for r in rows[c])  # noqa: E731
    for c in ("crack_uniaxial", "crack_shear", "tensstiff", "degenerate_tension", "implex_crack_shear"):
        assert mx(c, "DT") > 0.5, f"{c}: tension damage never grew past 0.5"
    assert mx("crush", "DC") > 0.3 and mx("degenerate_compression", "DC") > 0.3
    assert mn("interlock_capped", "BSR") < 0.5, "interlock cap never bound"
    assert mn("interlock_cyclic", "BSR") == -1.0, "cyclic interlock never re-capped after the reversal"
    assert any(r["LOAD"] == 0 for r in rows["interlock_cyclic"]), "no unloading steps on the cyclic path"
    assert all(mx(c, "CRK") == 1 for c in ("interlock_capped", "interlock_cyclic", "ps_interlock_cyclic"))


def test_dual_value_part_is_bit_identical_to_mode0(oracle):
    rows, cls, _ = oracle
    for c, rr in rows.items():
        for r in rr:
            if cls[c] == "plane-stress":
                assert r["BITS"] == -1
                continue
            assert r["BITS"] == 1, f"{c} step {r['step']}: mode 3 changed the stress or the history"


def test_plane_stress_newtons_converged(oracle):
    rows, _, _ = oracle
    for c, rr in rows.items():
        for r in rr:
            assert r["NFAIL"] == 0, f"{c} step {r['step']}: a sigma_33 Newton did not converge"


def test_dual_tangent_equals_fd_of_the_map(oracle):
    rows, cls, _ = oracle
    for c, rr in rows.items():
        for r in rr:
            fdh = r["fdh"]
            if cls[c] == "degenerate" and r["FD"] > FD_TOL:
                # a kink (the max-eigenvalue of a tensile cluster): central FD converges to the
                # dual at first order; an equal COMPRESSIVE cluster is smooth and takes the else
                assert r["FD"] <= FD_TOL_KINK, f"{c} step {r['step']}: FD {r['FD']:.3e}"
                head = fdh[:7]   # h = 1e-6 .. 1e-9, above round-off
                assert all(a > b for a, b in zip(head, head[1:])), f"{c} step {r['step']}: FD not converging {fdh}"
                assert head[0] / head[-1] >= 100.0, f"{c} step {r['step']}: FD slower than first order {fdh}"
            else:
                assert r["FD"] <= FD_TOL, (
                    f"{c} step {r['step']}: dual vs FD {r['FD']:.3e} > {FD_TOL:.0e} (FD by h: {fdh})")
                if r["FD"] > 1e-13:   # an exactly linear (elastic) state has no truncation to show
                    assert fdh[0] > 1e2 * r["FD"], f"{c} step {r['step']}: no h-convergence {fdh}"


def test_mode0_misses_the_dual_where_it_omits_terms(oracle):
    """Break-the-gate leg: the oracle must SEE what mode 0 leaves out."""
    rows, _, _ = oracle
    for r in rows["elastic"]:
        assert r["OMIT"] <= 1e-9, "on the elastic path mode 0 IS the exact tangent"
    for c in ROTATION:
        for r in rows[c]:
            if r["LOAD"] and r["DT"] > 0.05:
                assert r["OMIT"] >= OMIT_FLOOR, f"{c} step {r['step']}: mode 0 within {r['OMIT']:.2e} of the dual"
                assert r["OMIT"] > 1e5 * FD_TOL


def test_symmetrized_dual_is_flagged(oracle):
    rows, cls, _ = oracle
    for c, rr in rows.items():
        if cls[c] in ("elastic",):
            continue
        dmg = [r for r in rr if r["LOAD"] and (r["DT"] > 0.1 or r["DC"] > 0.1)]
        if dmg:
            worst = max(r["SYM"] for r in dmg)
            assert worst >= SYM_FLOOR, f"{c}: symmetrized dual passes the FD (worst {worst:.2e})"


def test_print_omitted_term_table(oracle):
    rows, cls, timing = oracle
    by = collections.OrderedDict()
    for c, rr in rows.items():
        by.setdefault(cls[c], []).extend(rr)
    print("\nWP-263 omitted-term table: ||C_dual - C_mode0||_F / ||C_dual||_F, per state class")
    print(f"{'class':18s} {'steps':>5s} {'max FD':>9s} {'median':>9s} {'max':>9s}")
    for k, rr in by.items():
        om = sorted(r["OMIT"] for r in rr)
        print(f"{k:18s} {len(rr):5d} {max(r['FD'] for r in rr):9.2e} {om[len(om) // 2]:9.3e} {om[-1]:9.3e}")
    if timing:
        print(f"per-call cost: mode 0 {timing.get(0, float('nan')):.2f} us, mode 3 {timing.get(3, float('nan')):.2f} us")
    assert set(by) == set(SCENARIOS.values())


def test_double_path_bit_identical_to_the_wp263_base():
    gpp = _gpp()
    git = shutil.which("git")
    if git is None:
        pytest.skip("git not available")
    show = subprocess.run([git, "-C", _REPO, "show", f"{PRE_SHA}:SRC/material/nD/LadrunoRCKernel.h"],
                          stdin=subprocess.DEVNULL, capture_output=True)
    if show.returncode != 0:
        pytest.skip(f"base object {PRE_SHA[:9]} not in this clone")
    with tempfile.TemporaryDirectory() as td:
        pre = os.path.join(td, "pre")
        os.makedirs(pre)
        with open(os.path.join(pre, "LadrunoRCKernel.h"), "wb") as fh:
            fh.write(show.stdout)
        now = _build_run(gpp, td, "dig_now", ["-DRC_DUAL_DOUBLE_ONLY"], [_SRC_ND])
        old = _build_run(gpp, td, "dig_pre", ["-DRC_DUAL_DOUBLE_ONLY"], [pre, _SRC_ND])
    assert now.count("DIGEST") >= 56, now
    assert now == old, "the double path of LadrunoRCKernel.h changed:\n" + "\n".join(
        f"  now {a}\n  pre {b}" for a, b in zip(now.splitlines(), old.splitlines()) if a != b)
