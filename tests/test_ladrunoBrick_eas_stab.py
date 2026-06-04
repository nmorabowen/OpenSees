"""LadrunoBrick EAS — `-stab` tangent regularization (ADR 20) validation.

ADR 20 adds an optional `-stab <beta>` knob to `-formulation eas`: a Tikhonov
floor on the condensed tangent, `K* = Kdd - Kda (Kaa + beta*Kaa0)^-1 Kad`, applied
ONLY at the condensation `.Solve()` (never to the accumulated `Kaa`, and NOT to the
inner-Newton direction). Because `Kaa` is a Jacobian and never enters the element
RESIDUAL, the converged state is a *modified-Newton* result: in the convex regime
it is EXACTLY beta-independent. These gates are that theorem's teeth, run first as a
correctness firewall before the (heavier, deferred) DEN-bar cure + plateau sweep.

  1. parser: `-stab` is accepted with `eas`, clamps negative beta, and is warned-and-
     ignored (still builds) for non-eas formulations;
  2. beta=0 is BIT-IDENTICAL to bare `eas` (the `easStabBeta==0` fast path);
  3. beta-independence, convex ELASTIC, to ~Newton-tol: a distorted constant-strain
     patch and a coarse bending cantilever give the SAME converged response for
     beta in {0, 1e-3, 1e-1, 1} — the direct check that beta never leaks into the
     residual (a failure here = beta was added to the accumulated Kaa, a bug);
  4. beta-independence, convex PLASTIC (monotonic hardening J2), same assertion in
     the inelastic-but-still-convex regime that exercises the nonlinear tangent.

NOTE the theorem's caveat (ADR 20 §How): beta-independence holds only for a
CONVERGED, iterating algorithm. Every gate below uses `algorithm Newton` with a
tight residual test — under `algorithm Linear` a regularized tangent WOULD bias the
result (that is the modified-Newton nature of the knob, not a bug).

Deferred (ADR 20 §Testing #1/#4, separate campaign): the notched DEN-bar robustness
cure and the beta-sweep plateau (the genuinely non-convex regime).
"""
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

_CUBE = {1: (0, 0, 0), 2: (1, 0, 0), 3: (1, 1, 0), 4: (0, 1, 0),
         5: (0, 0, 1), 6: (1, 0, 1), 7: (1, 1, 1), 8: (0, 1, 1)}
_CONN = [1, 2, 3, 4, 5, 6, 7, 8]

_BETAS = [0.0, 1.0e-3, 1.0e-1, 1.0]


# ===========================================================================
# 1. parser: -stab accepted (eas), clamped (negative), warned+ignored (non-eas)
# ===========================================================================
def _build_one(form_args):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, c in _CUBE.items():
        ops.node(t, *[float(x) for x in c])
    ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.3)
    ops.element("LadrunoBrick", 1, *_CONN, 1, *form_args)
    tags = ops.getEleTags()
    tags = [tags] if isinstance(tags, int) else (tags or [])
    return 1 in tags


def test_stab_parser_accepts_and_guards():
    assert _build_one(["-formulation", "eas", "-stab", 0.5]), \
        "-formulation eas -stab 0.5 must build"
    assert _build_one(["-formulation", "eas", "-stab", 0.0]), \
        "-stab 0 must build (bit-identical fast path)"
    # negative beta is clamped to 0 (warns) but still builds
    assert _build_one(["-formulation", "eas", "-stab", -1.0]), \
        "-stab with negative beta must clamp and still build"
    # -stab on a non-eas formulation is warned-and-ignored, not an error
    assert _build_one(["-formulation", "std", "-stab", 1.0]), \
        "-stab on std must be ignored (warned) and still build"
    assert _build_one(["-formulation", "bbar", "-stab", 1.0]), \
        "-stab on bbar must be ignored (warned) and still build"


# ===========================================================================
# Newton-based coarse bending cantilever (the convex-regime workhorse). Returns
# the mean tip transverse deflection. Uses algorithm Newton + a tight residual
# test so the beta-independence theorem applies (NOT algorithm Linear).
# ===========================================================================
def _cantilever_tip_newton(form_args, nx=4, L=10.0, b=1.0, h=1.0,
                           E=1000.0, nu=0.0, P=1.0, nsteps=1):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    dx = L / nx

    def nt(i, j, k):
        return 1 + i + (nx + 1) * j + (nx + 1) * 2 * k

    for k in range(2):
        for j in range(2):
            for i in range(nx + 1):
                ops.node(nt(i, j, k), i * dx, j * b, k * h)
    ops.nDMaterial("ElasticIsotropic", 1, E, nu)
    e = 1
    for i in range(nx):
        ops.element("LadrunoBrick", e,
                    nt(i, 0, 0), nt(i + 1, 0, 0), nt(i + 1, 1, 0), nt(i, 1, 0),
                    nt(i, 0, 1), nt(i + 1, 0, 1), nt(i + 1, 1, 1), nt(i, 1, 1),
                    1, *form_args)
        e += 1
    for k in range(2):
        for j in range(2):
            ops.fix(nt(0, j, k), 1, 1, 1)
    ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1)
    tip = [nt(nx, j, k) for k in range(2) for j in range(2)]
    for n in tip:
        ops.load(n, 0.0, 0.0, P / len(tip))
    ops.system("FullGeneral"); ops.numberer("RCM"); ops.constraints("Plain")
    ops.test("NormDispIncr", 1e-13, 200)
    ops.integrator("LoadControl", 1.0 / nsteps)
    ops.algorithm("Newton"); ops.analysis("Static")
    assert ops.analyze(nsteps) == 0, f"bending solve diverged for {form_args}"
    return sum(ops.nodeDisp(n)[2] for n in tip) / len(tip)


# ===========================================================================
# 2. beta=0 is BIT-IDENTICAL to bare eas (the easStabBeta==0 fast path)
# ===========================================================================
def test_stab_zero_is_bit_identical_to_bare_eas():
    """`-stab 0` must take the unregularized code path and reproduce bare `eas` to
    the last bit (the cheapest regression firewall: a non-zero diff means the
    beta==0 branch is not actually bypassing the regularization)."""
    bare = _cantilever_tip_newton(["-formulation", "eas"])
    stab0 = _cantilever_tip_newton(["-formulation", "eas", "-stab", 0.0])
    assert bare == stab0, (
        f"-stab 0 ({stab0!r}) must be bit-identical to bare eas ({bare!r})"
    )


# ===========================================================================
# 3a. beta-independence (convex ELASTIC) — coarse bending sweep to ~Newton-tol
# ===========================================================================
def test_stab_beta_independence_elastic_bending():
    """The theorem's headline: in the convex regime the converged tip deflection is
    EXACTLY beta-independent (Kaa never enters the residual). Sweep beta and require
    the tips to agree to ~Newton tolerance — NOT merely engineering tolerance. A
    failure means beta leaked into the residual (regularization wrongly applied to
    the accumulated Kaa instead of only the .Solve() operator)."""
    ref = _cantilever_tip_newton(["-formulation", "eas", "-stab", 0.0])
    assert abs(ref) > 1e-3, "bending tip must be non-trivial (non-vacuous test)"
    for beta in _BETAS:
        tip = _cantilever_tip_newton(["-formulation", "eas", "-stab", beta])
        assert abs(tip - ref) <= 1e-8 * (1.0 + abs(ref)), (
            f"beta={beta}: tip {tip:.12e} != beta=0 {ref:.12e} "
            f"(converged convex-regime state must be beta-independent)"
        )


# ===========================================================================
# 3b. beta-independence (convex ELASTIC) — distorted constant-strain patch
# ===========================================================================
def _G_field(X):
    g = ((2.0e-3, 1.0e-3, 0.5e-3), (0.0, -1.0e-3, 0.0), (0.0, 0.0, 1.5e-3))
    x, y, z = X
    return (g[0][0] * x + g[0][1] * y + g[0][2] * z,
            g[1][0] * x + g[1][1] * y + g[1][2] * z,
            g[2][0] * x + g[2][1] * y + g[2][2] * z)


def _patch_interior_disp(beta, E=1000.0, nu=0.3):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)

    def nid(i, j, k):
        return 1 + i + 3 * j + 9 * k

    base = {}
    for k in range(3):
        for j in range(3):
            for i in range(3):
                base[nid(i, j, k)] = [float(i), float(j), float(k)]
    base[nid(1, 1, 1)] = [1.0 + 0.18, 1.0 - 0.13, 1.0 + 0.21]   # distort interior
    for t, c in base.items():
        ops.node(t, *c)
    ops.nDMaterial("ElasticIsotropic", 1, E, nu)
    form = ["-formulation", "eas"] + (["-stab", beta] if beta else [])
    e = 1
    for k in range(2):
        for j in range(2):
            for i in range(2):
                conn = [nid(i, j, k), nid(i + 1, j, k), nid(i + 1, j + 1, k), nid(i, j + 1, k),
                        nid(i, j, k + 1), nid(i + 1, j, k + 1), nid(i + 1, j + 1, k + 1), nid(i, j + 1, k + 1)]
                ops.element("LadrunoBrick", e, *conn, 1, *form)
                e += 1
    interior = nid(1, 1, 1)
    ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1)
    for t, c in base.items():
        if t == interior:
            continue
        for dof, uu in enumerate(_G_field(c), start=1):
            ops.sp(t, dof, uu)
    ops.system("FullGeneral"); ops.numberer("Plain")
    ops.constraints("Penalty", 1e14, 1e14)
    ops.test("NormDispIncr", 1e-13, 200)
    ops.integrator("LoadControl", 1.0); ops.algorithm("Newton"); ops.analysis("Static")
    assert ops.analyze(1) == 0, f"patch solve diverged for beta={beta}"
    return list(ops.nodeDisp(interior))


def test_stab_patch_test_preserved_under_beta():
    """The constant-strain patch test stays EXACT for beta=1: the free interior node
    still lands on the linear field u=G.X. (Per the theorem this is automatic — the
    converged state is beta-independent — so it is a regression assertion, no longer
    an open question.)"""
    interior_X = [1.0 + 0.18, 1.0 - 0.13, 1.0 + 0.21]
    exp = _G_field(interior_X)
    ref = _patch_interior_disp(0.0)
    for beta in _BETAS:
        got = _patch_interior_disp(beta)
        for d in range(3):
            # exact vs the analytic linear field ...
            assert abs(got[d] - exp[d]) <= 1e-8 + 1e-6 * abs(exp[d]), (
                f"beta={beta}: patch interior u[{d}]={got[d]:.3e} != field {exp[d]:.3e}"
            )
            # ... and beta-independent vs the beta=0 solution
            assert abs(got[d] - ref[d]) <= 1e-8 * (1.0 + abs(ref[d])), (
                f"beta={beta}: patch interior u[{d}] drifted from beta=0"
            )


# ===========================================================================
# 4. beta-independence (convex PLASTIC) — monotonic hardening J2 on a distorted hex
# ===========================================================================
def _push_j2_distorted_tip(beta, nsteps=8, push=6.0):
    """Distorted hex, J2Plasticity (linear hardening => convex), monotonic push past
    yield under Newton. Returns (tip_ux, tip_uz) at the end."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    coords = {1: (0, 0, 0), 2: (1, 0, 0), 3: (1.1, 1.05, 0.0), 4: (0.05, 1, 0),
              5: (0, 0, 1), 6: (1, 0.05, 1.1), 7: (1.05, 1, 1), 8: (0, 1, 1.05)}
    for t, c in coords.items():
        ops.node(t, *[float(x) for x in c])
    K = 1000.0 / (3.0 * (1.0 - 2.0 * 0.3)); G = 1000.0 / (2.0 * (1.0 + 0.3))
    ops.nDMaterial("J2Plasticity", 1, K, G, 2.0, 2.0, 0.0, 50.0)   # linear hardening H=50
    for n in (1, 4, 5, 8):
        ops.fix(n, 1, 1, 1)
    form = ["-formulation", "eas"] + (["-stab", beta] if beta else [])
    ops.element("LadrunoBrick", 1, *_CONN, 1, *form)
    ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1)
    for n in (2, 3, 6, 7):
        ops.load(n, push, 0.0, 0.0)   # push well past yield
    ops.system("FullGeneral"); ops.numberer("Plain"); ops.constraints("Plain")
    ops.test("NormDispIncr", 1e-12, 200)
    ops.integrator("LoadControl", 1.0 / nsteps)
    ops.algorithm("Newton"); ops.analysis("Static")
    assert ops.analyze(nsteps) == 0, f"J2 push diverged for beta={beta}"
    return [ops.nodeDisp(2)[0], ops.nodeDisp(7)[2]]


def test_stab_beta_independence_plastic_hardening():
    """In the convex INELASTIC regime (monotonic, linearly-hardening J2) the
    converged displacement is still beta-independent — the -stab regularization
    touches only the tangent, and the inner Newton (on the true Kaa) drives the same
    enhanced state. Confirms the nonlinear tangent path is wired correctly."""
    ref = _push_j2_distorted_tip(0.0)
    assert abs(ref[0]) > 1e-3, "the J2 push must be non-trivial (non-vacuous test)"
    for beta in _BETAS:
        got = _push_j2_distorted_tip(beta)
        for a, b in zip(got, ref):
            assert abs(a - b) <= 1e-7 * (1.0 + abs(b)), (
                f"beta={beta}: plastic tip {a:.10e} != beta=0 {b:.10e} "
                f"(convex inelastic converged state must be beta-independent)"
            )
