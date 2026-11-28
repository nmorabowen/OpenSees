"""WP-265: SHELLMOD tangent gate -- a shell on LadrunoShellModifier, K == d(force)/d(u) plastic.

The SHELLMOD mutation gate sabotages the LadrunoShellModifier section accessors
(SRC/material/section/LadrunoShellModifierSection.cpp: getSectionTangent and
getInitialTangent under IDENT return the 8x8 identity; the stress resultant stays
exact). Before WP-265 no SHELLMOD test compared a tangent with a finite difference:
the section-tangent checks compared against the closed form S.D.S or a bare
section, on elastic states.

The rig: one ShellMITC4 on a LadrunoShellModifier with eight distinct, non-unit
modifiers, wrapping a 4-layer LayeredShell of LadrunoJ2 plate fibres (Voce + two
Armstrong-Frederick terms). Every DOF is driven: two COMMITTED imposed steps
(Lagrange) of a membrane + bending + twist field, then the probe as a TRIAL state
set with setNodeDisp. ShellMITC4 strains its section from the nodal trial
displacements on every force/tangent call, so a tangent read after a commit would
be the tangent of a zero increment from the new committed state (the WP-265 eas
quirk, same mechanism); setNodeDisp drives it without committing any material.
WP-158 quirk: '-commit' on every component, or each call wipes the ones set before
it (it commits the node, never the materials). Shells answer no 'stiff' response,
so K is read from printA (see _element_K).

Asserted: the probe is plastic (the FD departs from the elastic K) and
max|K - K_fd| <= FD_TOL * max|K|; the identity and the elastic tangents miss by
more than 1e3 * FD_TOL (break-the-gate legs).

Measured 2026-11-28 on the WP-262 build (ladrunoBuild ac449a03c): max|K - K_fd| /
max|K| = 2.0e-7; the identity misses by 0.77 and the elastic tangent by 0.55.
"""
import numpy as np
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

_E, _NU = 200.0, 0.3
_K = _E / (3.0 * (1.0 - 2.0 * _NU))
_G = _E / (2.0 * (1.0 + _NU))
_H, _NLAYERS = 0.1, 4
_MODS = {"-f11": 0.7, "-f22": 0.9, "-f12": 0.6, "-m11": 0.5, "-m22": 0.8, "-m12": 0.4,
         "-v13": 0.9, "-v23": 0.7}
_NODES = {1: (0.0, 0.0, 0.0), 2: (1.0, 0.05, 0.0), 3: (1.05, 1.0, 0.0), 4: (0.0, 0.95, 0.0)}
_PATH = (0.4, 0.7)
# h = 1e-6 is the measured optimum: 2.0e-7 (1e-5: 1.3e-6 truncation; 1e-7: 1.1e-6
# and 1e-8: 3.9e-6, the noise floor of the plate-fibre sigma_zz condensation).
FD_H = 1.0e-6
FD_TOL = 1.0e-5
PLASTIC_MIN = 0.05

# Only ShellMITC4: it re-strains its section on every force/tangent call, so the
# trial rig drives it. ASDShellQ4 does not (setNodeDisp leaves its force unchanged:
# K_fd == 0), and on the committed-step rig its tangent misses the FD by 1.0e-3,
# independent of h (1e-6..1e-8) and of the modifiers (all-unit: 1.3e-3), in the
# u, v, w, theta_x blocks: an ASDShellQ4 + LayeredShell(LadrunoJ2) property, not a
# SHELLMOD one. Reported by WP-265 as an out-of-scope finding; not gated here.
CASES = {"shellmitc4": "ShellMITC4"}


def _target():
    """24-vector (6 DOF x 4 nodes): membrane stretch/shear + bending + twist + a
    deterministic wiggle (no RNG). Every DOF is imposed, so any field is admissible."""
    exx, eyy, gxy = 0.012, -0.004, 0.006
    kx, ky, kxy = 0.15, -0.08, 0.10
    u = np.zeros(24)
    for tag, (x, y, _z) in _NODES.items():
        w = -(0.5 * kx * x * x + 0.5 * ky * y * y + kxy * x * y)
        rx = -(ky * y + kxy * x)          # dw/dy
        ry = (kx * x + kxy * y)           # -dw/dx
        wig = 2.0e-4 * np.sin(np.arange(6) * 1.3 + 0.7 * tag)
        u[6 * (tag - 1):6 * tag] = np.array([exx * x + 0.5 * gxy * y, eyy * y + 0.5 * gxy * x,
                                             w, rx, ry, 0.002 * tag]) + wig
    return u


def _build(ele):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    for tag, (x, y, z) in _NODES.items():
        ops.node(tag, x, y, z)
    ops.nDMaterial("LadrunoJ2", 1, _K, _G, "-iso", "voce", 0.30, 0.10, 20.0, 5.0,
                   "-kin", 2, 8.0, 50.0, 4.0, 20.0)
    layers = []
    for _ in range(_NLAYERS):
        layers += [1, _H / _NLAYERS]
    ops.section("LayeredShell", 10, _NLAYERS, *layers)
    mods = []
    for k, v in _MODS.items():
        mods += [k, v]
    ops.section("LadrunoShellModifier", 11, 10, *mods)
    ops.element(ele, 1, 1, 2, 3, 4, 11)
    ops.constraints("Lagrange")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-12, 50, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    ops.timeSeries("Constant", 1)


def _committed_path(ele, path=_PATH):
    """Imposed committed steps. The LAST step's pattern is kept, so its Lagrange
    constraint objects stay alive for printA's formTangent (removing it would leave
    the analysis model pointing at deleted SP constraints until the next analyze)."""
    _build(ele)
    u_t = _target()
    for k, a in enumerate(path):
        pat = 100 + k
        ops.pattern("Plain", pat, 1)
        for tag in _NODES:
            for d in range(6):
                ops.sp(tag, d + 1, float(a * u_t[6 * (tag - 1) + d]))
        assert ops.analyze(1) == 0, f"imposed step {k} failed"
        if k < len(path) - 1:
            ops.remove("loadPattern", pat)


def _element_K():
    """The shell's 24x24 tangent at the CURRENT trial state. Shells answer no
    'stiff' response, so it is read from the assembled system: printA re-forms the
    tangent (formTangent) before printing; the node-DOF block holds the element alone
    (the Lagrange SP rows only couple node and multiplier DOFs). WP-158 quirk:
    printA -ret is the column-major buffer, so the reshape is transposed."""
    A = np.array(ops.printA("-ret"), dtype=float)
    n = int(round(np.sqrt(A.size)))
    A = A.reshape(n, n).T
    eq = [e for tag in sorted(_NODES) for e in ops.nodeDOFs(tag)]
    assert len(eq) == 24 and min(eq) >= 0
    return A[np.ix_(eq, eq)]


def _probe_trial(u, want_k=False):
    for tag in _NODES:
        for d in range(6):
            ops.setNodeDisp(tag, d + 1, float(u[6 * (tag - 1) + d]), "-commit")
    f = np.array(ops.eleForce(1), dtype=float)
    if not want_k:
        return f
    return f, _element_K()


def _diag(case):
    ele = CASES[case]
    # the elastic reference: the same rig after one tiny imposed step (strain ~1e-6)
    _committed_path(ele, path=(1.0e-4,))
    K0 = _element_K()
    _committed_path(ele)
    probe = _probe_trial
    u = _target()
    Kfd = np.zeros((24, 24))
    for d in range(24):
        up = u.copy(); up[d] += FD_H
        um = u.copy(); um[d] -= FD_H
        Kfd[:, d] = (probe(up) - probe(um)) / (2.0 * FD_H)
    f0, K = probe(u, want_k=True)
    scale = np.abs(K).max()
    return {"finite": float(np.isfinite(K).all() and np.isfinite(f0).all()),
            "plastic": np.abs(Kfd - K0).max() / np.abs(K0).max(),
            "err": np.abs(K - Kfd).max() / scale,
            "ident_err": np.abs(np.eye(24) - Kfd).max() / np.abs(Kfd).max()}


@pytest.mark.parametrize("case", sorted(CASES))
def test_shellmod_tangent_matches_fd_at_a_plastic_state(case):
    r = _diag(case)
    assert r["finite"] == 1.0
    assert r["plastic"] >= PLASTIC_MIN, (
        f"{case}: probe not plastic enough (|K_fd-K0|/|K0| = {r['plastic']:.3e})")
    assert r["err"] <= FD_TOL, (
        f"{case}: shell tangent on LadrunoShellModifier != FD of the resisting force: "
        f"max rel err {r['err']:.3e} (> {FD_TOL:.0e})")
    assert r["plastic"] >= 1.0e3 * FD_TOL and r["ident_err"] >= 1.0e3 * FD_TOL
