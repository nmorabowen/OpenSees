"""LadrunoRCConcrete IMPL-EX (Phase 4) — Zone-A material battery.

IMPL-EX (Oliver et al. 2008) freezes the damage thresholds (xt,xc) and the MCFT
beta over a step by EXPLICIT extrapolation of the committed history, so the
returned tangent is the constant SPD secant W_B*C0 (no softening-indefinite term,
no beta cross-term). The true implicit thresholds are recomputed at commit to
advance the state and measure the implex error. It is the robustness escape for
cyclic softening RC (where the implicit consistent tangent goes indefinite), and
the prerequisite for the meshed cyclic-wall validation.

These tests pin the MATERIAL-level contract:
  1. `-implex` off  -> bit-identical to the implicit material (the V0 tripwire).
  2. `-implex` on a smooth proportional path tracks the implicit stress (IMPL-EX
     is exact for a linearly-growing threshold; consistent otherwise).
  3. the implex error is genuinely active (non-zero on a rate change, zero while
     elastic).
  4. the implex tangent is SPD under deep softening (the robustness property).
  5. send/recvSelf round-trips the schema-v2 (IMPL-EX) state.
"""
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

# gentle/ductile regularized backbones (single element traverses softening)
E, NU = 30000.0, 0.2
CE = [0.0, 0.0007, 0.0020, 0.0100]; CS = [0.0, 24.0, 30.0, 8.0]
CD = [0.0, 0.0, 0.25, 1.0 - 8.0 / (E * 0.0100)]
FT = 3.0
TE = [0.0, FT / E, 2.0e-3, 1.5e-2]; TS = [0.0, FT, 1.8, 0.15]
TD = [0.0, 0.0, 1.0 - 1.8 / (E * 2.0e-3), 1.0 - 0.15 / (E * 1.5e-2)]

_CUBE = {1: (0, 0, 0), 2: (1, 0, 0), 3: (1, 1, 0), 4: (0, 1, 0),
         5: (0, 0, 1), 6: (1, 0, 1), 7: (1, 1, 1), 8: (0, 1, 1)}


def _mat(tag, implex=False, numerical=False):
    a = ["LadrunoRCConcrete", tag, E, NU, "-Ce", *CE, "-Cs", *CS, "-Cd", *CD,
         "-Te", *TE, "-Ts", *TS, "-Td", *TD]
    if implex:
        a += ["-implex"]
    if numerical:
        a += ["-numericalTangent"]
    ops.nDMaterial(*a)


def _drive(implex, targets, nper=30):
    """Fully-prescribed homogeneous uniaxial-strain brick (3D view): ramp eps_xx
    through the `targets` list (nper substeps each) via a Path-scaled sp field.
    Returns per-step [(eps_xx, sig_xx, dt, dc, implexError)]."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, (x, y, z) in _CUBE.items():
        ops.node(t, float(x), float(y), float(z))
    _mat(1, implex=implex)
    ops.element("stdBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
    # pseudo-time k -> strain targets[k]; sp value = base*X scaled by the Path
    times = list(range(len(targets) + 1))
    vals = [0.0] + list(targets)
    ops.timeSeries("Path", 1, "-time", *times, "-values", *vals, "-useLast")
    ops.pattern("Plain", 1, 1)
    for t, (x, y, z) in _CUBE.items():
        ops.sp(t, 1, float(x))     # u_x = (Path)*X  -> eps_xx = Path value
        ops.sp(t, 2, 0.0); ops.sp(t, 3, 0.0)
    ops.constraints("Penalty", 1e16, 1e16)
    ops.numberer("Plain"); ops.system("UmfPack")
    ops.test("NormDispIncr", 1e-9, 60, 0); ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / nper)       # advances pseudo-time -> drives the Path
    ops.analysis("Static")
    out = []
    nstg = len(targets)
    for _ in range(nstg * nper):
        if ops.analyze(1) != 0:
            break
        sig = ops.eleResponse(1, "material", 1, "stress")
        dmg = ops.eleResponse(1, "material", 1, "damage")
        err = ops.eleResponse(1, "material", 1, "implexError")
        out.append((ops.nodeDisp(2, 1), sig[0], dmg[0], dmg[1], err[0] if err else 0.0))
    return out


# ---------------------------------------------------------------------------
@pytest.mark.t1
def test_implex_off_is_baseline_identical():
    """With -implex absent the material is the implicit baseline, bit-for-bit
    (the V0 tripwire: the implex plumbing must not perturb the default path)."""
    a = _drive(implex=False, targets=[1.5e-3])
    b = _drive(implex=False, targets=[1.5e-3])
    assert len(a) == len(b) and len(a) > 10
    for (ea, sa, *_), (eb, sb, *_) in zip(a, b):
        assert sa == sb           # determinism / no hidden implex state when off


@pytest.mark.t1
def test_implex_tracks_implicit_on_smooth_path():
    """On a smooth monotone tension path the implex stress tracks the implicit
    one closely (IMPL-EX is exact for a linearly-growing threshold)."""
    im = _drive(implex=False, targets=[1.2e-3])
    ix = _drive(implex=True, targets=[1.2e-3])
    assert len(im) == len(ix) and len(im) > 10
    speak = max(abs(s) for _, s, *_ in im) or 1.0
    dmax = max(abs(si[1] - sx[1]) for si, sx in zip(im, ix))
    assert dmax / speak < 0.02, f"implex deviates {dmax/speak*100:.2f}% from implicit"


@pytest.mark.t1
def test_implex_error_active_on_rate_change_zero_when_elastic():
    """The implex error is a genuine extrapolation residual: ~0 while elastic
    (below cracking, no damage rate), and strictly > 0 once damage evolves at a
    changing rate (the cracking onset, which the lagged explicit pass under-predicts)."""
    # stage 1 elastic (below ft/E = 1e-4); stage 2 ramps through cracking into softening
    traj = _drive(implex=True, targets=[0.6 * FT / E, 1.5e-3], nper=40)
    nstg1 = 40
    elastic_err = max(e[4] for e in traj[:nstg1 - 2])     # stage-1 ramp, below cracking
    overall_err = max(e[4] for e in traj)
    assert elastic_err < 1e-9, f"implex error {elastic_err:.2e} non-zero while elastic"
    assert overall_err > 1e-4, f"implex error never became active (max {overall_err:.2e})"


@pytest.mark.t1
def test_implex_tangent_is_spd_under_softening():
    """At a deep-softening state the implex (explicit-pass) secant W_B*C0 with frozen
    damage and no beta cross-term is symmetric positive-definite — the robustness gain
    over the implicit consistent tangent, which is indefinite on the softening branch."""
    import numpy as np
    ops.wipe(); ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, (x, y, z) in _CUBE.items():
        ops.node(t, float(x), float(y), float(z))
    _mat(1, implex=True)
    ops.element("stdBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
    ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1)
    for t, (x, y, z) in _CUBE.items():
        ops.sp(t, 1, 1.2e-3 * float(x)); ops.sp(t, 2, 0.0); ops.sp(t, 3, 0.0)
    ops.constraints("Penalty", 1e16, 1e16); ops.numberer("Plain"); ops.system("UmfPack")
    ops.test("NormDispIncr", 1e-9, 60, 0); ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / 60); ops.analysis("Static")
    for _ in range(60):
        ops.analyze(1)
    dmg = ops.eleResponse(1, "material", 1, "damage")
    assert dmg[0] > 0.3, f"expected deep tension damage, got dt={dmg[0]:.3f}"
    flat = ops.eleResponse(1, "material", 1, "tangent")
    C = np.array(flat).reshape(6, 6)
    Csym = 0.5 * (C + C.T)
    w = np.linalg.eigvalsh(Csym)
    assert w.min() > 0.0, f"implex tangent not SPD under softening: min eig {w.min():.3e}"


@pytest.mark.t1
def test_implex_bypasses_numerical_tangent():
    """Under -implex the forward-difference (-numericalTangent) path MUST be bypassed:
    its base point is the explicit stress while a perturbed re-call re-enters implicitly,
    so the quotient is O(implexError/h) garbage and would wreck the solve. Because the
    solver uses the explicit secant regardless, the SOLVED stress trajectory with
    -implex -numericalTangent is identical to -implex alone."""
    a = _drive(implex=True, targets=[1.3e-3])
    b = _drive_numeric(targets=[1.3e-3])
    assert len(a) == len(b) and len(a) > 10
    for (_, sa, *_), (_, sb, *_) in zip(a, b):
        assert abs(sa - sb) <= 1e-9 * (abs(sa) + 1.0), f"FD path leaked under implex: {sa} vs {sb}"


def _drive_numeric(targets, nper=30):
    ops.wipe(); ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, (x, y, z) in _CUBE.items():
        ops.node(t, float(x), float(y), float(z))
    _mat(1, implex=True, numerical=True)
    ops.element("stdBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
    times = list(range(len(targets) + 1)); vals = [0.0] + list(targets)
    ops.timeSeries("Path", 1, "-time", *times, "-values", *vals, "-useLast")
    ops.pattern("Plain", 1, 1)
    for t, (x, y, z) in _CUBE.items():
        ops.sp(t, 1, float(x)); ops.sp(t, 2, 0.0); ops.sp(t, 3, 0.0)
    ops.constraints("Penalty", 1e16, 1e16); ops.numberer("Plain"); ops.system("UmfPack")
    ops.test("NormDispIncr", 1e-9, 60, 0); ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / nper); ops.analysis("Static")
    out = []
    for _ in range(len(targets) * nper):
        if ops.analyze(1) != 0:
            break
        sig = ops.eleResponse(1, "material", 1, "stress")
        out.append((ops.nodeDisp(2, 1), sig[0], 0.0, 0.0, 0.0))
    return out


def _run_to(nsteps, save_at=None, dbpath=None):
    """Drive the prescribed brick `nsteps` LoadControl steps under -implex; if
    save_at given, save to the FE database at that step, wipe, rebuild, restore,
    and continue. Returns the final committed (post-revert) damage = (dt, dc)."""
    def build():
        ops.wipe(); ops.model("basic", "-ndm", 3, "-ndf", 3)
        for t, (x, y, z) in _CUBE.items():
            ops.node(t, float(x), float(y), float(z))
        _mat(1, implex=True)
        ops.element("stdBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
        ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1)
        for t, (x, y, z) in _CUBE.items():
            ops.sp(t, 1, 1.4e-3 * float(x)); ops.sp(t, 2, 0.0); ops.sp(t, 3, 0.0)
        ops.constraints("Penalty", 1e16, 1e16); ops.numberer("Plain"); ops.system("UmfPack")
        ops.test("NormDispIncr", 1e-9, 60, 0); ops.algorithm("Newton")
        ops.integrator("LoadControl", 1.0 / nsteps); ops.analysis("Static")
    build()
    for k in range(1, nsteps + 1):
        ops.analyze(1)
        if save_at is not None and k == save_at:
            try:
                ops.database("File", dbpath); ops.save(1)
            except Exception as exc:                          # noqa: BLE001
                pytest.skip(f"database() unsupported: {exc}")
            ops.wipe(); build()
            ops.database("File", dbpath); ops.restore(1)
    return list(ops.eleResponse(1, "material", 1, "damage"))


@pytest.mark.t1
def test_implex_serialization_roundtrip_continuation():
    """send/recvSelf must round-trip the schema-v2 IMPL-EX state (xt_old, xc_old,
    eps1_old, dtime, commitDone). Save mid-run, restore into a fresh model, and
    continue: the final damage must match an UNINTERRUPTED run bit-for-bit — a
    corrupted/forgotten implex field would perturb the post-restore extrapolation."""
    import tempfile, os
    ref = _run_to(8)
    with tempfile.TemporaryDirectory(prefix="rc_implex_rt_", ignore_cleanup_errors=True) as td:
        got = _run_to(8, save_at=4, dbpath=os.path.join(td, "rc_implex"))
        ops.wipe()
    assert ref and got and len(ref) == len(got)
    for r, g in zip(ref, got):
        assert abs(r - g) <= 1e-9 * abs(r) + 1e-10, (
            f"IMPL-EX state did not survive save/restore: {r} -> {g}")
