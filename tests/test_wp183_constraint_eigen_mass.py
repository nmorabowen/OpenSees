"""WP-183: constraint FEs (Penalty / Lagrange MP, SP, EQ) must not enter the mass.

Two upstream defects (identical code on OpenSees/OpenSees master):

  A  `PenaltyMP_FE::getM_Force` (and PenaltyEQ, LagrangeMP/SP/EQ) returned
     `*resid`, the last constraint residual -alpha C^T C (u - u0). ARPACK's
     M*x product (`ArpackSolver::myMv`, used whenever M is not diagonal) and
     `IncrementalIntegrator::doMv` (modalDamping setup) summed it in, so after
     any step that evaluated the residual at a deformed state (Newton, and
     Linear after its first step) M*x became M*x + b: eigenvalues ~ +-1e-4.
  B  The eigen mass loops (`EigenIntegrator::formM`, `StaticAnalysis::eigen`,
     `DirectIntegrationAnalysis::eigen`, the two DD analyses) added every FE's
     `getTangent`; a constraint FE ignores the integrator and returns its
     penalty stiffness, which went into M. That made M non-diagonal (sending
     ARPACK down path A) and, with explicit-M solvers (-fullGenLapack,
     -symmBandLapack, -feast), gave lambda ~ 1 at every constrained dof.

Every handler must give the Transformation eigenvalues, at K0, after a Newton
gravity, after a multi-step Linear analysis, for every eigen solver; and a
modalDamping transient built on an eigen after gravity must match too.

Classic Tcl (OpenSees.exe): the defect sits in `StaticAnalysis::eigen`, which
is the path classic `eigen` takes once a static analysis exists.
"""
import os
import re
import subprocess
from pathlib import Path

import pytest

pytestmark = [pytest.mark.zone_a]

REPO = Path(__file__).resolve().parents[1]


def _find_exe():
    exe = os.environ.get("LADRUNO_TCL_EXE")
    if exe:
        return exe if Path(exe).exists() else None
    cand = REPO / "dist" / "bin" / ("OpenSees.exe" if os.name == "nt" else "OpenSees")
    return str(cand) if cand.exists() else None


TCL_EXE = _find_exe()
needs_exe = pytest.mark.skipif(TCL_EXE is None, reason="classic OpenSees exe not found (build it, or set LADRUNO_TCL_EXE)")

REL_TOL = 1.0e-6
NMODES = 3
TRANSFORMATION = "Transformation"
# handler -> eigenvalue tolerance vs Transformation. A penalty handler carries
# its own approximation error ~ k/alpha (measured on this portal after the fix:
# Penalty 1e13 <= 2.6e-7, Auto -autoPenalty 5 <= 1.7e-8, Auto at its default
# penalty 1.5e-6). The first two are held to the 1e-6 gate; default Auto gets
# 1e-5, still five orders below the defect (O(1) before WP-183). Lagrange is
# exact (no penalty error): its multiplier dofs now carry zero mass, which
# ARPACK, dense LAPACK (infinite eigenvalues dropped) and FEAST all handle.
HANDLERS = {
    "Penalty 1.0e13 1.0e13": REL_TOL,
    "Auto -autoPenalty 5": REL_TOL,
    "Auto": 1.0e-5,
    "Lagrange": REL_TOL,
}
# ARPACK (default = -genBandArpack), the explicit-M dense LAPACK solver, FEAST.
# -symmBandLapack is absent on purpose: SymBandEigenSolver refuses the
# generalized problem ("only does standard problem"), with or without M.
SOLVERS = ["default", "-fullGenLapack", "-feast"]
STAGES = ["K0", "NEWTON", "LINEAR"]

# 2D portal (2 columns, 2 beam segments) in kN, m, t. Lumped mass on every
# free dof. SI units keep the stiffness range narrow (lateral 12EI/L^3 ~ 4e4,
# axial EA/L ~ 2e6, rotational 4EI/L ~ 2e5 kN m), so a penalty of 1e13 is
# both stiff enough (approximation error ~ k/alpha ~ 1e-7) and well
# conditioned. (In N-mm the base rotational stiffness is 1.7e11, and a 1e12
# penalty SP is a semi-rigid base: a 10 % eigenvalue shift that is the
# penalty method itself, not this defect.)
MODEL = """
wipe
model basic -ndm 2 -ndf 3
node 1 0 0; node 2 5 0; node 3 0 3.5; node 4 5 3.5; node 5 2.5 3.5
fix 1 1 1 1; fix 2 1 1 1
mass 3 10 10 0.5; mass 4 10 10 0.5; mass 5 20 20 1.0
geomTransf Linear 1
element elasticBeamColumn 1 1 3 0.25 3.0e7 0.005 1
element elasticBeamColumn 2 2 4 0.25 3.0e7 0.005 1
element elasticBeamColumn 3 3 5 0.20 3.0e7 0.004 1
element elasticBeamColumn 4 5 4 0.20 3.0e7 0.004 1
"""

# mp=False drops the equalDOFs: the handler then sees only the `fix` SPs
# (the SP-only Penalty case, PenaltySP_FE alone).
GRAVITY = "pattern Plain 1 Linear { load 3 10 -200 0; load 4 0 -200 0; load 5 0 -400 0 }"

EIG_PROC = """
proc eig {tag} {
  global ES
  if {$ES eq "default"} {
    set l [eigen %(n)d]
  } elseif {$ES eq "-feast"} {
    set l [eigen -feast 0.0 60.0]
  } else {
    set l [eigen $ES %(n)d]
  }
  puts "EIG $tag [join [lmap x $l {format %%.15e $x}] { }]"
}
""" % {"n": NMODES}


def _deck(handler, solver, mp=True):
    lines = [MODEL]
    if mp:
        lines.append("equalDOF 3 4 1; equalDOF 3 5 1")
    lines += [
        f"set ES {solver}",
        f"constraints {handler}",
        "numberer RCM",
        "system UmfPack",
        EIG_PROC,
        "eig K0",
        GRAVITY,
        "test NormDispIncr 1e-12 20",
        "algorithm Newton",
        "integrator LoadControl 0.1",
        "analysis Static",
        "analyze 10",
        "eig NEWTON",
        "wipeAnalysis",
        "loadConst -time 0.0",
        "pattern Plain 2 Linear { load 3 20 0 0 }",
        f"constraints {handler}",
        "numberer RCM",
        "system UmfPack",
        "test NormDispIncr 1e-12 20",
        "algorithm Linear",
        "integrator LoadControl 0.1",
        "analysis Static",
        "analyze 10",
        "eig LINEAR",
    ]
    return "\n".join(lines) + "\n"


def _run(tmp_path, name, deck):
    d = tmp_path / name
    d.mkdir()
    (d / "m.tcl").write_text(deck)
    p = subprocess.run([TCL_EXE, "m.tcl"], cwd=d, capture_output=True, text=True, timeout=300)
    out = p.stdout + p.stderr
    assert p.returncode == 0, out
    return d, out


def _eigs(out):
    res = {}
    for m in re.finditer(r"^EIG (\w+)((?: \S+)*)\s*$", out, re.M):
        res[m.group(1)] = [float(x) for x in m.group(2).split()]
    return res


def _assert_same(got, ref, what, tol=REL_TOL):
    assert len(got) == len(ref) and len(ref) > 0, (what, got, ref)
    for g, r in zip(got, ref):
        assert r > 0.0, (what, ref)
        assert abs(g - r) <= tol * abs(r), (what, got, ref)


_ref_cache = {}


def _reference(tmp_path, solver, mp):
    key = (solver, mp)
    if key not in _ref_cache:
        _, out = _run(tmp_path, f"ref_{solver.strip('-')}_{int(mp)}", _deck(TRANSFORMATION, solver, mp))
        e = _eigs(out)
        assert set(e) == set(STAGES), out
        _ref_cache[key] = e
    return _ref_cache[key]


@needs_exe
@pytest.mark.parametrize("solver", SOLVERS)
@pytest.mark.parametrize("handler", list(HANDLERS))
def test_mp_eigen_matches_transformation(tmp_path, handler, solver):
    ref = _reference(tmp_path, solver, True)
    _, out = _run(tmp_path, "h", _deck(handler, solver, True))
    got = _eigs(out)
    for stage in STAGES:
        _assert_same(got.get(stage, []), ref[stage], f"{handler} {solver} {stage}", HANDLERS[handler])


@needs_exe
@pytest.mark.parametrize("solver", SOLVERS)
def test_sp_only_penalty_eigen_matches_transformation(tmp_path, solver):
    ref = _reference(tmp_path, solver, False)
    _, out = _run(tmp_path, "h", _deck("Penalty 1.0e13 1.0e13", solver, False))
    got = _eigs(out)
    for stage in STAGES:
        _assert_same(got.get(stage, []), ref[stage], f"SP-only Penalty {solver} {stage}")


TH = MODEL + """
equalDOF 3 4 1; equalDOF 3 5 1
constraints %(h)s
numberer RCM
system UmfPack
set l [eigen 3]
pattern Plain 1 Linear { load 3 0 -200 0; load 4 0 -200 0; load 5 0 -400 0 }
test NormDispIncr 1e-12 20; algorithm Newton; integrator LoadControl 0.1; analysis Static; analyze 10
loadConst -time 0.0
set l [eigen 3]
puts "EIG AFTER [join [lmap x $l {format %%.15e $x}] { }]"
%(damp)s
timeSeries Trig 2 0.0 1.0 0.5 -factor 2.0
pattern UniformExcitation 3 1 -accel 2
test NormDispIncr 1e-10 20; algorithm Newton; integrator Newmark 0.5 0.25; analysis Transient
set fo [open th.txt w]
for {set i 0} {$i < 400} {incr i} { analyze 1 0.005; puts $fo "[getTime] [format %%.15e [nodeDisp 3 1]] [format %%.15e [nodeDisp 4 2]]" }
close $fo
"""

DAMPING = {
    "modal": "modalDamping 0.05",
    # fixed coefficient (0.05 * omega_1 of the reference) so the control does
    # not depend on the eigen and can be compared bitwise across binaries
    "rayleigh": "rayleigh 1.7587 0 0 0",
}


def _history(d):
    rows = [list(map(float, ln.split())) for ln in (d / "th.txt").read_text().split("\n") if ln.strip()]
    assert len(rows) == 400
    return rows


@needs_exe
@pytest.mark.parametrize("damp", ["modal", "rayleigh"])
@pytest.mark.parametrize("handler", list(HANDLERS))
def test_transient_after_gravity_matches_transformation(tmp_path, handler, damp):
    """modalDamping is built on the eigen after gravity: it read the polluted
    M*phi through doMv (48 %% drift error before WP-183). Rayleigh is the
    control: it never went through getM_Force."""
    if handler == "Lagrange" and damp == "modal":
        pytest.skip("pre-existing, not WP-183: IncrementalIntegrator::setupModal asks every "
                    "DOF_Group for eigenvectors, and a LagrangeDOF_Group has no node "
                    "(DOF_Group::getEigenvectors exit(-1)s)")
    dref, out_ref = _run(tmp_path, "ref", TH % {"h": TRANSFORMATION, "damp": DAMPING[damp]})
    dh, out_h = _run(tmp_path, "h", TH % {"h": handler, "damp": DAMPING[damp]})
    tol = HANDLERS[handler]
    _assert_same(_eigs(out_h)["AFTER"], _eigs(out_ref)["AFTER"], f"{handler} eigen after gravity", tol)
    ref, got = _history(dref), _history(dh)
    for col in (1, 2):
        peak = max(abs(r[col]) for r in ref)
        err = max(abs(g[col] - r[col]) for g, r in zip(got, ref))
        assert err <= tol * peak, (handler, damp, col, err, peak)
