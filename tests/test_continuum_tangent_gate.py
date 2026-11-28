"""WP-265: CONTINUUM tangent gate -- K == d(resisting force)/d(u) at a plastic state.

The mutation gate's IDENT mode (SRC/Ladruno_mutation.h) replaces the element
tangent with the identity and leaves the resisting force exact. A suite can only
kill it if something CHECKS the tangent. Before WP-265 the only CONTINUUM tests
that compared K with a finite difference of the force were the `-geom finite`
and `-geom corot` ones (test_ladrunoBrick_finite.py, test_ladrunoBrick_corot.py),
on elastic or path-independent states. Nothing pinned the small-strain
LadrunoBrick formulations (std, bbar, eas, ssp, uri) or LadrunoBrick20 against an
FD at a genuinely NONLINEAR, history-dependent state.

The rig, per formulation:

  * one distorted element, LadrunoJ2 (Voce + two Armstrong-Frederick terms), every
    DOF imposed by `sp` under the Lagrange handler: the constraint rows land u
    exactly whatever the tangent, so the rig itself cannot be fooled by a wrong K
    (the WP-158 quirk's trap of `setNodeDisp` trial wiping is avoided by never
    probing with it);
  * two COMMITTED steps of an affine + per-node wiggle deformation, then the probe
    step to the target u, left as the trial (Constant-series pattern per step,
    the test_ladrunoBrick_hypo.py `_step_to` idiom);
  * K = eleResponse('stiff') at the probe; the FD re-runs the same committed path
    and perturbs only the probe target, so the central difference is the
    derivative of the trial map from the same committed history -- exactly what
    the consistent tangent is.

Asserted: the probe state is plastic (the FD departs from 'stiffInitial'), and
max|K - K_fd| <= FD_TOL * max|K|. Under IDENT K is the identity and the second
assertion fails by ~1 (relative); under a correct tangent it passes at the FD's
truncation level.

Measured 2026-11-28 on the WP-262 build (ladrunoBuild ac449a03c; no C++ in this WP
changes the binary), max|K - K_fd| / max|K|: brick std 7.2e-11, bbar 1.4e-10,
eas 1.6e-8, ssp 1.1e-10, uri 6.9e-11, brick20 std 3.1e-11. An identity tangent
misses by 0.95-0.99 and the elastic one by 0.28-0.63 (FD_TOL = 1e-5).

`-formulation eas` is probed as a trial state, not a committed step: see
_probe_trial and the quirk WP-265-eas-stiff-after-commit-is-a-zero-increment-tangent.
"""
import numpy as np
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

_E, _NU = 200.0, 0.3
_K = _E / (3.0 * (1.0 - 2.0 * _NU))
_G = _E / (2.0 * (1.0 + _NU))

_CORNERS = [
    (0.00, 0.00, 0.00), (1.00, 0.10, 0.05), (1.10, 1.00, 0.00), (0.05, 0.95, 0.10),
    (0.00, 0.05, 1.00), (1.00, 0.00, 1.05), (1.05, 1.00, 1.10), (0.00, 1.00, 0.95),
]
_EDGES = [(1, 2), (2, 3), (3, 4), (4, 1), (5, 6), (6, 7), (7, 8), (8, 5),
          (1, 5), (2, 6), (3, 7), (4, 8)]

# a moderate 3D deformation: stretch + shear (yield strain ~ 0.0015, so ~10x past yield)
_FBAR = np.array([[1.020, 0.008, -0.004],
                  [0.003, 0.990, 0.006],
                  [-0.002, 0.005, 1.012]])
_WIGGLE_AMP = 1.5e-3
_PATH = (0.4, 0.7)          # committed steps (fractions of the target)
FD_H = 1.0e-7
FD_TOL = 1.0e-5             # relative to max|K|
PLASTIC_MIN = 0.05          # max|K_fd - K_initial| / max|K_initial| at the probe

# probe mode: "commit" = the probe is one more imposed, COMMITTED step and K is read
# after the commit (valid where getTangentStiff re-reads the material's last trial
# tangent); "trial" = the probe displacements are set with setNodeDisp on the
# committed path and never committed (needed by eas: see _probe_trial).
CASES = {
    "brick_std":  ("LadrunoBrick", ["-formulation", "std"], "commit"),
    "brick_bbar": ("LadrunoBrick", ["-formulation", "bbar"], "commit"),
    "brick_eas":  ("LadrunoBrick", ["-formulation", "eas"], "trial"),
    "brick_ssp":  ("LadrunoBrick", ["-formulation", "ssp"], "commit"),
    "brick_uri":  ("LadrunoBrick", ["-formulation", "uri", "-hourglass", "stiffness"], "commit"),
    "brick20_std": ("LadrunoBrick20", [], "commit"),
}


def _nodes(ele):
    nodes = {i + 1: c for i, c in enumerate(_CORNERS)}
    if ele == "LadrunoBrick20":
        for k, (a, b) in enumerate(_EDGES):
            ca, cb = np.array(_CORNERS[a - 1]), np.array(_CORNERS[b - 1])
            nodes[9 + k] = tuple(0.5 * (ca + cb))
    return nodes


def _target(nodes):
    """Nodal displacement vector of the probe state (deterministic, no RNG)."""
    u = np.zeros(3 * len(nodes))
    I = np.eye(3)
    for tag, X in nodes.items():
        X = np.asarray(X)
        w = _WIGGLE_AMP * np.array([np.sin(3.1 * tag), np.cos(1.7 * tag), np.sin(2.3 * tag + 0.5)])
        u[3 * (tag - 1):3 * tag] = (_FBAR - I) @ X + w
    return u


def _build(ele, args, nodes):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for tag, (x, y, z) in nodes.items():
        ops.node(tag, x, y, z)
    ops.nDMaterial("LadrunoJ2", 1, _K, _G, "-iso", "voce", 0.30, 0.10, 20.0, 5.0,
                   "-kin", 2, 8.0, 50.0, 4.0, 20.0)
    ops.element(ele, 1, *sorted(nodes), 1, *args)
    ops.constraints("Lagrange")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-12, 50, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    ops.timeSeries("Constant", 1)


def _step_to(nodes, u, k):
    pat = 100 + k
    ops.pattern("Plain", pat, 1)
    for tag in nodes:
        for d in range(3):
            ops.sp(tag, d + 1, float(u[3 * (tag - 1) + d]))
    assert ops.analyze(1) == 0, f"imposed step {k} failed"
    ops.remove("loadPattern", pat)


def _committed_path(ele, args, nodes):
    _build(ele, args, nodes)
    u_t = _target(nodes)
    for k, a in enumerate(_PATH):
        _step_to(nodes, a * u_t, k)


def _read(n, want_k):
    f = np.array(ops.eleForce(1), dtype=float)
    if not want_k:
        return f
    K = np.array(ops.eleResponse(1, "stiff"), dtype=float)
    K0 = np.array(ops.eleResponse(1, "stiffInitial"), dtype=float)
    assert K.size == n * n and K0.size == n * n
    return f, K.reshape(n, n), K0.reshape(n, n)


def _probe_commit(ele, args, nodes, u_probe, want_k=False):
    """Committed path, then the probe as one more imposed step (committed by
    analyze). Rebuilt per call, so every FD probe starts from a bit-identical
    committed history."""
    _committed_path(ele, args, nodes)
    _step_to(nodes, u_probe, len(_PATH))
    return _read(len(u_probe), want_k)


def _probe_trial(nodes, u_probe, want_k=False):
    """The probe as a TRIAL state on the already-committed path (eas).

    `-formulation eas` re-solves its enhanced parameters and re-strains every GP on
    EVERY assembly (formEAStrue, called from getResistingForce and getTangentStiff
    alike), so after a commit 'stiff' is the tangent of a ZERO increment from the
    NEW committed state -- elastic here -- not the converged step's consistent
    tangent (measured on the commit rig: K == K_initial exactly at a plastic
    state). Because eas strains its GPs from the nodal trial displacements on
    every call, setNodeDisp alone drives it. WP-158 quirk: pass '-commit' on each
    component, or each call wipes the ones set before it; it commits the NODE,
    never the materials, which keep the path's committed history.
    """
    for tag in nodes:
        for d in range(3):
            ops.setNodeDisp(tag, d + 1, float(u_probe[3 * (tag - 1) + d]), "-commit")
    return _read(len(u_probe), want_k)


def _probe(case, nodes, u_probe, want_k=False):
    ele, args, mode = CASES[case]
    if mode == "commit":
        return _probe_commit(ele, args, nodes, u_probe, want_k)
    return _probe_trial(nodes, u_probe, want_k)


def _diag(case):
    """The measurements the gate asserts on (also printed by hand when re-measuring)."""
    ele, args, mode = CASES[case]
    nodes = _nodes(ele)
    u = _target(nodes)
    if mode == "trial":
        _committed_path(ele, args, nodes)          # once: trial probes never commit
    f0, K, K0 = _probe(case, nodes, u, want_k=True)
    Kfd = np.zeros_like(K)
    for d in range(len(u)):
        up = u.copy(); up[d] += FD_H
        um = u.copy(); um[d] -= FD_H
        Kfd[:, d] = (_probe(case, nodes, up) - _probe(case, nodes, um)) / (2.0 * FD_H)
    scale = np.abs(K).max()
    return {"finite": float(np.isfinite(K).all() and np.isfinite(f0).all()),
            "plastic": np.abs(Kfd - K0).max() / np.abs(K0).max(),
            "err": np.abs(K - Kfd).max() / scale,
            "ident_err": np.abs(np.eye(len(u)) - Kfd).max() / np.abs(Kfd).max()}


@pytest.mark.parametrize("case", sorted(CASES))
def test_continuum_tangent_matches_fd_at_a_plastic_state(case):
    r = _diag(case)
    assert r["finite"] == 1.0
    # the FD (not K, which IDENT replaces) proves the state is nonlinear
    assert r["plastic"] >= PLASTIC_MIN, (
        f"{case}: the probe state is not plastic enough (|K_fd-K0|/|K0| = {r['plastic']:.3e}); "
        f"the gate would only check the elastic tangent")
    assert r["err"] <= FD_TOL, (
        f"{case}: element tangent != FD of the resisting force at a plastic state: "
        f"max rel err {r['err']:.3e} (> {FD_TOL:.0e})")
    # break-the-gate legs, in process: the two wrong tangents this gate exists to
    # reject would miss by orders of magnitude more than FD_TOL -- the ELASTIC
    # tangent (the "plastic" margin above) and the IDENTITY (the IDENT mutant).
    assert r["plastic"] >= 1.0e3 * FD_TOL and r["ident_err"] >= 1.0e3 * FD_TOL
