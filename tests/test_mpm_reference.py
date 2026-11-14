"""WP-232b: the numpy MPM oracle checked against closed forms and its own frozen cases.

Every `expect` and `compare` line of tests/mpm_cases.txt, and every case of
tests/hypo_incremental_cases.txt, becomes one test, so each case is reproduced from the case files
alone (interface: tests/mpm_interface.md). The closed forms of the 232b brief are cases there:
the 1D vibrating bar (period; PIC vs FLIP dissipation), the self-weight column and the
cell-crossing meter, rigid translation (engine) and rigid rotation (kernel), the two-disc
impact, the frictional incline, the critical time step (refused above, stable just below).
The direct tests at the end check properties that need no case file.

numpy only (zone_a). Wall time about 45 s, dominated by the run cases, each run once per
session (cached).
"""
from __future__ import annotations

import functools
import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import hypo_incremental_reference as H  # noqa: E402
import mpm_reference as M  # noqa: E402

pytestmark = pytest.mark.zone_a

CASES = M.parse_cases()
HYPO = H.parse_cases()

# The frozen tolerance uGIMP cannot meet (interface Amendment 4; ledger quirk
# WP-232b-ugimp-fixed-domains-do-not-converge-under-large-compression). strict=True: when phase 3's
# domain-tracking basis makes it pass, the xfail must be removed. The closed-form gate for slice 1
# is column_selfweight_small_2d_gimp.
KNOWN_SHORTFALLS = {
    ("column_selfweight_2d_gimp", "selfweight_error"):
        "uGIMP with fixed domains does not converge at 20% compression (0.072 at h=0.5, "
        "0.067 at 0.25, 0.106 at 0.125 > 0.03); acceptance test for the phase-3 basis",
}


@functools.lru_cache(maxsize=None)
def result(name):
    return M.run_case(CASES[name])


def _close(got, exp, tol):
    got = np.asarray(got, dtype=float).ravel()
    exp = np.asarray(exp, dtype=float).ravel()
    assert got.shape == exp.shape, (got, exp)
    a, r = tol
    bad = np.abs(got - exp) > a + r * np.abs(exp)
    assert not np.any(bad), f"got {got[bad]} expected {exp[bad]} (tol abs {a}, rel {r})"


# --------------------------------------------------------------------------- #
#  Case-driven tests                                                          #
# --------------------------------------------------------------------------- #

def _expect_items():
    for name, lines in CASES.items():
        k = 0
        for t in lines:
            if t[0] in ("expect", "compare"):
                yield pytest.param(name, k, id=f"{name}:{k}:{t[0]}-{t[1]}")
                k += 1


def _line(name, k):
    return [t for t in CASES[name] if t[0] in ("expect", "compare")][k]


def _case_tol(lines):
    t = M._get(lines, "tol")
    return (float(t[0]), float(t[1])) if t else (0.0, 0.0)


@pytest.mark.parametrize("name,k", list(_expect_items()))
def test_case_expectation(name, k, request):
    lines = CASES[name]
    tok = _line(name, k)
    kind = M._get(lines, "kind")[0]

    if tok[0] == "compare":
        obs_args, op, other = tok[1:-2], tok[-2], tok[-1]
        obs, args = obs_args[0], obs_args[1:]
        a = M.observable(result(name), obs, args)[0]
        b = M.observable(result(other), obs, args)[0]
        assert (a < b) if op == "lt" else (a > b), f"{name} {obs}={a} not {op} {other}'s {b}"
        return

    body = tok[1:]
    if body[0] == "refused":
        assert result(name)["refused"], f"{name} was not refused"
        return
    if body[0] == "accepted":
        assert not result(name)["refused"], result(name).get("refusal")
        return
    if kind == "shape1d":
        assert body[0] == "s1"
        k_eq = body.index("eq")
        _close(M.shape1d_values(body), M._floats(body[k_eq + 1:k_eq + 3]), _case_tol(lines))
        return
    if body[0] == "node":
        nodes = result(name)["node"]
        nid = int(body[1])
        assert nid in nodes, f"node {nid} not in the support {sorted(nodes)}"
        _close(nodes[nid], M._floats(body[3:]), _case_tol(lines))
        return

    obs, args, op, vals, tol = M.split_expect(body)
    if (name, obs) in KNOWN_SHORTFALLS:
        request.applymarker(pytest.mark.xfail(reason=KNOWN_SHORTFALLS[(name, obs)], strict=True))
    res = result(name)
    assert not res.get("refused"), res.get("refusal")

    if obs == "particle_fields":
        assert vals == ["fields"], "particle_fields still pending"
        ndm = res["engine"].ndm
        stored = {int(t[1]): t[2:] for t in lines if t[0] == "field"}
        assert set(stored) == set(res["fields"]), "particle ids differ"
        for pid, toks in stored.items():
            nums = [float(x) for x in toks if x not in ("x", "v", "F", "V", "sig", "sighat", "epshat")]
            _close(res["fields"][pid], nums, tol or _case_tol(lines))
        assert ndm in (2, 3)
        return

    assert M.PENDING not in vals, f"{name}: {obs} still pending"
    got = M.observable(res, obs, args)
    if op == "eq":
        _close(got, M._floats(vals), tol or _case_tol(lines))
    elif op == "le":
        assert np.all(got <= float(vals[0])), f"{name}: {obs}={got} > {vals[0]}"
    elif op == "ge":
        assert np.all(got >= float(vals[0])), f"{name}: {obs}={got} < {vals[0]}"


def test_shape_cases_list_the_whole_support():
    """A 'kind shape' case lists every node with non-zero S or gradient; no other node has any."""
    for name, lines in CASES.items():
        if M._get(lines, "kind")[0] != "shape":
            continue
        listed = {int(t[2]) for t in lines if t[0] == "expect" and t[1] == "node"}
        assert set(result(name)["node"]) == listed, name


def test_no_pending_values_remain():
    for path in (M.CASES_PATH, H.CASES_PATH):
        with open(path, encoding="ascii") as fh:
            body = [ln for ln in fh if not ln.lstrip().startswith("#")]
        assert not any(M.PENDING in ln for ln in body), path


@pytest.mark.parametrize("idx", range(len(HYPO)), ids=[c["name"][0] for c in HYPO])
def test_hypo_incremental_case(idx):
    case = HYPO[idx]
    got = H.evaluate_case(case)
    assert int(got["status"][0]) == int(case["status"][0])
    if int(case["status"][0]) != 0:
        return
    tol = (float(case["tol"][0]), float(case["tol"][1]))
    for key in ("deps", "R1", "J1", "scauchy"):
        if key in case:
            assert case[key] != [H.PENDING], f"{key} pending"
            _close(got[key], H._arr(case[key]), tol)


# --------------------------------------------------------------------------- #
#  Direct property tests (no case file)                                       #
# --------------------------------------------------------------------------- #

def _rot(rng, max_angle):
    axis = rng.normal(size=3)
    axis /= np.linalg.norm(axis)
    th = rng.uniform(-max_angle, max_angle)
    K = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
    return np.eye(3) + np.sin(th) * K + (1 - np.cos(th)) * (K @ K)


def test_kernel_objective_under_rotation_increments():
    """F1 = Q Fn: deps = 0 and R1 = Q polar(Fn) (Hughes & Winget 1980; Cayley transform)."""
    rng = np.random.default_rng(2321)
    for _ in range(50):
        S = 0.2 * (lambda A: 0.5 * (A + A.T))(rng.normal(size=(3, 3)))
        Fn = _rot(rng, 1.0) @ (np.eye(3) + S)
        Q = _rot(rng, 3.0)
        st, deps, R1, J1 = H.hypo_increment_from_F(Fn, Q @ Fn)
        assert st == 0
        assert np.max(np.abs(deps)) < 1e-12
        assert np.allclose(R1, Q @ H.polar_rotation(Fn), atol=1e-12)
        assert np.isclose(J1, np.linalg.det(Fn), rtol=1e-13)


def test_kernel_plane_strain_is_the_exact_embedding():
    rng = np.random.default_rng(2322)
    for _ in range(20):
        th = rng.uniform(-1, 1)
        R = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
        Fn = R @ (np.eye(2) + 0.1 * rng.normal(size=(2, 2)))
        F1 = (np.eye(2) + 0.05 * rng.normal(size=(2, 2))) @ Fn
        st, d6, R3, J = H.hypo_increment_from_F(H.embed_plane_strain(Fn), H.embed_plane_strain(F1))
        assert st == 0
        assert d6[2] == 0.0 and d6[4] == 0.0 and d6[5] == 0.0   # exactly, by structure
        assert R3[2, 2] == 1.0 and not np.any(R3[:2, 2]) and not np.any(R3[2, :2])
        st2, d3, R2, J2 = H.hypo_increment_from_F_plane_strain(Fn, F1)
        assert np.array_equal(d3, d6[[0, 1, 3]]) and np.array_equal(R2, R3[:2, :2]) and J2 == J


def test_kernel_small_increment_is_the_symmetric_velocity_gradient():
    """For F1 = (I + dt L) Fn with small dt: deps -> R^T sym(L) R dt to O(dt^2)."""
    rng = np.random.default_rng(2323)
    L = rng.normal(size=(3, 3))
    Fn = _rot(rng, 0.7) @ (np.eye(3) + 0.1 * np.diag(rng.normal(size=3)))
    R = H.polar_rotation(Fn)
    D = 0.5 * (L + L.T)
    for dt in (1e-3, 1e-4):
        _, deps, _, _ = H.hypo_increment_from_F(Fn, (np.eye(3) + dt * L) @ Fn)
        Dm = R.T @ D @ R * dt
        ref = np.array([Dm[0, 0], Dm[1, 1], Dm[2, 2], 2 * Dm[0, 1], 2 * Dm[1, 2], 2 * Dm[2, 0]])
        assert np.max(np.abs(deps - ref)) < 50 * dt ** 2


def test_basis_partition_of_unity():
    """Sum S = 1 and sum grad S = 0 for every particle whose domain lies inside the grid."""
    rng = np.random.default_rng(2324)
    for ndm, basis in ((2, "gimp"), (3, "gimp"), (2, "linear"), (3, "linear")):
        g = M.Grid(np.zeros(ndm), 0.5, [6] * ndm)
        lp = rng.uniform(0.02, 0.25, size=(200, ndm))
        x = rng.uniform(lp, 3.0 - lp)
        nid, S, G = g.shape(x, lp, basis)
        assert np.allclose(S.sum(axis=1), 1.0, atol=1e-14)
        assert np.allclose(G.sum(axis=1), 0.0, atol=1e-12)


def test_basis_gradient_is_the_derivative():
    """dS/dr equals a central difference of S, in all four regions (C1 uGIMP)."""
    h, lp = 1.0, 0.3
    r = np.linspace(-1.35, 1.35, 541)
    r = r[np.min(np.abs(r[:, None] - np.array([-1.3, -0.7, -0.3, 0.3, 0.7, 1.3])[None]), axis=1) > 1e-3]
    e = 1e-6
    S_p, _ = M.gimp_1d(r + e, h, np.full(r.shape, lp))
    S_m, _ = M.gimp_1d(r - e, h, np.full(r.shape, lp))
    _, dS = M.gimp_1d(r, h, np.full(r.shape, lp))
    assert np.max(np.abs((S_p - S_m) / (2 * e) - dS)) < 1e-6


def test_momentum_is_conserved_without_conditions():
    """No faces, no gravity, no damping, no cutoff: total momentum is exact for FLIP and PIC."""
    for flip in (1.0, 0.0, 0.6):
        lines = [t[:] for t in CASES["impact_discs_2d"]]
        lines = [t for t in lines if t[0] not in ("steps", "flip", "expect")]
        lines += [["steps", "200"], ["flip", repr(flip)]]
        eng = M.build_engine(lines)
        p0 = np.sum(eng.m[:, None] * eng.v, axis=0)
        for _ in range(200):
            eng.step(0.001)
        assert np.allclose(np.sum(eng.m[:, None] * eng.v, axis=0), p0, atol=1e-11)


def test_failed_step_changes_nothing():
    """A refused step (dt above the estimate) leaves every particle array untouched."""
    eng = M.build_engine(CASES["critdt_refuse_above_2d"])
    before = {k: getattr(eng, k).copy() for k in ("x", "v", "F", "V", "sig", "eps_hat")}
    with pytest.raises(M.StepRefused):
        eng.step(1.01 * eng.critical_dt() * eng.safety)
    for k, a in before.items():
        assert np.array_equal(getattr(eng, k), a), k
