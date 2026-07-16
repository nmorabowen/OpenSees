"""LadrunoPorousOverlay (PATTERN 33022) - ADR-73 P1 FRAMEWORK battery (WP1.E-ii).

The *framework* half of the P1 battery (the physics half lives in
``tests/test_ladruno_overlay_physics.py``). Everything here exercises the
overlay's plumbing rather than its consolidation accuracy:

  1. pattern-timing bit-constancy <A-5>   -- the fluid/force refresh advances at
     COMMIT granularity, never per Newton iteration; holds under Static
     (LoadControl), Newmark, AND CentralDifference; multi-iteration Newton steps
     still produce exactly one record row/commit; runs are bit-deterministic.
  2. factor-ignored semantics             -- the overlay carries no TimeSeries and
     ignores the domain load-factor machinery; no spurious "ignored" notice.
  3. revertToLastCommit consistency       -- a failed+reverted step does NOT
     advance the fluid; after retry the p-history (incl. the dpCommitted_
     rate-form reference used at the NEXT step) is bit-identical to a clean run.
  4. region-validation fatals <A-13>      -- each in a fresh subprocess: unsupported
     cell, empty -drained, overlapping overlays, LadrunoUP in region, unknown
     flag (-substep / -bogus), missing -moduli, -fluidUpdate explicit NYI.
  5. pInit steady/hydrostatic closed form -- hydrostatic gamma_w*(z_w - z) clipped
     at 0 is node-exact; steady seepage under -fluidBody is a node-exact linear
     profile.
  6. serial DB round-trip (the 1.D gate)  -- committed p AND dpCommitted_ round-trip
     bit-exact (isolated with a rigid skeleton); a rich config restores without
     crashing and the overlay stays active; deformable transient restart
     reproduces the continuous run (hard gate since #577).
  7. -record everyN semantics             -- everyN=2 writes every other commit;
     the header row is the region node-tag set.

RUNNER (single command, exits nonzero on any unexpected failure):

    C:\\Users\\nmora\\AppData\\Local\\Python\\pythoncore-3.12-64\\python.exe -S \\
        tests/test_ladruno_overlay_framework.py

`-S` is mandatory (project memory): the installer's boot .pth otherwise PRELOADS
a stale opensees.pyd. `-S` also drops site-packages, so this file is numpy-free
by design (pure Python floats/CSV). The bootstrap below asserts the imported
engine is THIS worktree's build.

MEASURED FINDINGS (2026-07-14, reported upward - see the WP1.E-ii handoff):
  * FileDatastore does NOT restore a sibling Plain LoadPattern's nodal loads
    (base reaction collapses post-restore) - a framework limitation, so the DB
    gate is driven purely by the overlay (self-contained), never an external load.
  * DB round-trip of a DEFORMABLE model was NON-DETERMINISTIC - RESOLVED:
    root cause was UPSTREAM (#577, landed 2026-07-15): quad/tri elements never
    serialized their ELEMENT-LEVEL rho and the broker ctor left it
    uninitialized, so every DB/MPI-restored element rebuilt its mass matrix
    from heap garbage (overriding the correctly-restored material density via
    the `if (rho==0)` gate). The P1 adjudication that exonerated the overlay
    stands verified: its restored state (p/dp/uSnapshot, per-node d/v/a, SPs)
    was in-situ exact in diverging runs; the overlay only amplified the
    wrong-M response through its +Q*p forces. (The P1-era "any LoadPattern in
    the DB stream triggers it" attribution was an artifact - the no-pattern
    control was nearly motionless, so the corrupt M went unexercised; the true
    trigger is mass-driven response after restore.) Test 6c is now a HARD
    GATE (promoted per its XPASS contract): 6/6 fresh-process deformable
    restarts must reproduce the continuous run.
  * A -pInit list (explicit $nd $val ...) overlay additionally crashes the
    process on FileDatastore restore (observed with a frozen rig); hydrostatic
    pInit is crash-free. Reported, not separately gated.
"""
import os
import subprocess
import sys
import tempfile
from pathlib import Path

# --------------------------------------------------------------------------
# bootstrap: THIS worktree's engine (evict any boot-.pth preloaded stale pyd)
# --------------------------------------------------------------------------
_DIST = str(Path(__file__).resolve().parents[1] / "dist" / "bin")
if not os.path.isfile(os.path.join(_DIST, "opensees.pyd")):
    # Collection-safe skip. Under pytest (e.g. Zone-A `-m zone_a` still IMPORTS
    # every module before filtering) a module-level sys.exit raises SystemExit
    # during collection and crashes the WHOLE session — use pytest's module-
    # level skip there; only a direct `python -S` run exits.
    _msg = f"worktree engine not built: {_DIST}"
    if "pytest" in sys.modules:
        import pytest
        pytest.skip(_msg, allow_module_level=True)
    print(f"SKIP: {_msg}")
    raise SystemExit(0)

os.environ["PATH"] = _DIST + os.pathsep + os.environ.get("PATH", "")
try:
    os.add_dll_directory(_DIST)
except (FileNotFoundError, OSError):
    pass
if _DIST not in sys.path:
    sys.path.insert(0, _DIST)
for _m in ("opensees", "openseespy", "openseespy.opensees"):
    sys.modules.pop(_m, None)
import opensees as ops  # noqa: E402

assert os.path.normcase(os.path.dirname(ops.__file__)) == os.path.normcase(_DIST), (
    f"wrong opensees.pyd imported: {ops.__file__} (want {_DIST})"
)

# pytest support (optional; the primary runner is `python -S <file>`)
try:
    import pytest  # noqa: F401
    pytestmark = [pytest.mark.zone_b]
except Exception:  # pragma: no cover
    pytest = None

_SCRATCH = tempfile.mkdtemp(prefix="ov_fw_")

# common overlay knobs (k-bar = k/gamma_w; near-incompressible fluid)
_OV = ("-Kf", 2.2e7, "-rhoF", 1000.0, "-perm", 1e-6, 1e-6,
       "-poro", 0.4, "-moduli", 1.0e4, 0.2)


# ==========================================================================
# shared model builders
# ==========================================================================
def _column(nely=6, L=6.0, E=1.0e4, mat="elastic"):
    """2-wide plane-strain Q4 column, base uy fixed, ux fixed everywhere.
    Returns (nid map, element tags, top-node list, L, nely)."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    nid = {}
    k = 1
    for j in range(nely + 1):
        for i in range(2):
            ops.node(k, float(i), L * j / nely)
            nid[(i, j)] = k
            k += 1
    if mat == "j2":
        # yields under the test load -> genuine multi-iteration Newton
        ops.nDMaterial("J2Plasticity", 1, 5.0e3, 3.0e3, 8.0, 12.0, 0.1, 20.0)
    else:
        ops.nDMaterial("ElasticIsotropic", 1, E, 0.2, 2.0)
    eles = []
    for j in range(nely):
        et = j + 1
        ops.element("quad", et, nid[(0, j)], nid[(1, j)], nid[(1, j + 1)],
                    nid[(0, j + 1)], 1.0, "PlaneStrain", 1)
        eles.append(et)
    for j in range(nely + 1):
        for i in range(2):
            ops.fix(nid[(i, j)], 1, (1 if j == 0 else 0))
    top = [nid[(0, nely)], nid[(1, nely)]]
    return nid, eles, top, L, nely


def _newmark(gamma=0.6, beta=0.3025, tol=1e-12, maxit=25, algo="Newton"):
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("ProfileSPD")
    ops.test("NormDispIncr", tol, maxit)
    ops.algorithm(algo)
    ops.integrator("Newmark", gamma, beta)
    ops.analysis("Transient")


def _read_csv(path):
    """Return (header-tags[list[int]], rows[list[list[float]]])."""
    with open(path) as f:
        lines = f.read().strip().splitlines()
    hdr = [int(x) for x in lines[0].split(",")[1:]]
    rows = [[float(x) for x in r.split(",")] for r in lines[1:]]
    return hdr, rows


def _tmp(name):
    return os.path.join(_SCRATCH, name)


# ==========================================================================
# 1. pattern-timing bit-constancy <A-5>
# ==========================================================================
def test_bit_constancy_once_per_commit_all_integrators():
    """<A-5>: the overlay refreshes its fluid state (and thus the +Q*p force it
    injects) exactly once per COMMIT, never per Newton iteration.

    Strongest observables the Python surface allows, and what each proves:
      (a) ONCE-PER-COMMIT across three integrators (Newmark transient,
          LoadControl static with -staticMode steady, CentralDifference
          explicit): N committed analyze-steps produce exactly N record rows
          (everyN=1). If the fluid advanced per iteration the row count would
          exceed the commit count.
      (b) MULTI-ITERATION INVARIANCE: a J2-plastic column solved with Newton
          taking several iterations/step still yields exactly one row/step, and
          two byte-for-byte identical runs produce a BIT-IDENTICAL p record.
          A per-iteration-varying injected force would make the converged
          committed state (and therefore the next-step p) sensitive to the
          iteration trajectory, breaking determinism / the clean row count.
    """
    # (a) Newmark
    nid, eles, top, L, nely = _column()
    rec = _tmp("bc_newmark.csv")
    ops.pattern("LadrunoPorousOverlay", 77, "-region", *eles, *_OV,
                "-drained", *top, "-record", rec, 1)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for t in top:
        ops.load(t, 0.0, -5.0)
    _newmark()
    for _ in range(7):
        assert ops.analyze(1, 0.1) == 0
    hdr, rows = _read_csv(rec)
    assert len(rows) == 7, f"Newmark: {len(rows)} rows for 7 commits"

    # (a) LoadControl static + -staticMode steady
    nid, eles, top, L, nely = _column()
    rec = _tmp("bc_static.csv")
    ops.pattern("LadrunoPorousOverlay", 77, "-region", *eles, *_OV,
                "-drained", *top, "-staticMode", "steady",
                "-fluidBody", 0.0, -9.81, "-record", rec, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t in top:
        ops.load(t, 0.0, -5.0)
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("ProfileSPD")
    ops.test("NormDispIncr", 1e-10, 25)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    for _ in range(4):
        assert ops.analyze(1) == 0
    hdr, rows = _read_csv(rec)
    assert len(rows) == 4, f"static-steady: {len(rows)} rows for 4 commits"

    # (a) CentralDifference explicit
    nid, eles, top, L, nely = _column()
    rec = _tmp("bc_cd.csv")
    ops.pattern("LadrunoPorousOverlay", 77, "-region", *eles, *_OV,
                "-drained", *top, "-record", rec, 1)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for t in top:
        ops.load(t, 0.0, -5.0)
    ops.constraints("Plain")
    ops.numberer("RCM")
    ops.system("ProfileSPD")
    ops.integrator("CentralDifference")
    ops.algorithm("Linear")
    ops.analysis("Transient")
    for _ in range(6):
        assert ops.analyze(1, 1e-4) == 0
    hdr, rows = _read_csv(rec)
    assert len(rows) == 6, f"CentralDifference: {len(rows)} rows for 6 commits"

    # (b) multi-iteration Newton (J2) -> still one row/commit + determinism
    def run_j2(rec):
        nid, eles, top, L, nely = _column(mat="j2")
        ops.pattern("LadrunoPorousOverlay", 77, "-region", *eles, *_OV,
                    "-drained", *top, "-record", rec, 1)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        for t in top:
            ops.load(t, 0.0, -200.0)
        _newmark(tol=1e-11, maxit=60)
        iters = []
        for _ in range(5):
            assert ops.analyze(1, 0.05) == 0
            iters.append(ops.testIter())
        return iters

    rec1 = _tmp("bc_j2_a.csv")
    rec2 = _tmp("bc_j2_b.csv")
    it1 = run_j2(rec1)
    it2 = run_j2(rec2)
    h1, r1 = _read_csv(rec1)
    h2, r2 = _read_csv(rec2)
    assert len(r1) == 5, f"J2 multi-iter: {len(r1)} rows for 5 commits"
    assert max(it1) >= 2, f"J2 model did not iterate (iters={it1})"
    # bit-identical across the two runs
    dmax = max(abs(a - b) for ra, rb in zip(r1, r2) for a, b in zip(ra, rb))
    assert dmax == 0.0, f"non-deterministic p record across identical runs: {dmax}"
    assert it1 == it2, f"iteration counts differ across identical runs: {it1} {it2}"
    return f"once/commit Newmark=7 static=4 CD=6; J2 iters={it1} rows=5 det=bit-exact"


# ==========================================================================
# 2. factor-ignored semantics
# ==========================================================================
def test_factor_ignored_semantics():
    """The overlay owns its force amplitudes and ignores the LoadPattern factor
    machinery. On the Python surface a TimeSeries can never be ATTACHED to the
    overlay (the `pattern LadrunoPorousOverlay` command takes no $tsTag - unlike
    `pattern Plain`), so the design guarantee is structural. Checks:
      (a) the overlay runs to completion WITHOUT any timeSeries argument;
      (b) adding an unrelated Plain pattern carrying a wildly ramping Linear
          timeSeries (nonzero load on a disconnected, fully-fixed dummy node)
          leaves the overlay p record BIT-IDENTICAL - the overlay does not pick
          up the evolving domain load factor;
      (c) a normal run emits NO 'TimeSeries ... IGNORED' notice (checked in a
          subprocess that captures C++ stderr). The one-time notice fires only
          if setTimeSeries() is ever called, which the parser never does, so
          'at most once' is satisfied vacuously - documented here.
    """
    def run(rec, with_dummy_series):
        nid, eles, top, L, nely = _column()
        ops.pattern("LadrunoPorousOverlay", 77, "-region", *eles, *_OV,
                    "-drained", *top, "-record", rec, 1)
        ops.timeSeries("Constant", 1)
        ops.pattern("Plain", 1, 1)
        for t in top:
            ops.load(t, 0.0, -5.0)
        if with_dummy_series:
            # a disconnected, fully-fixed dummy node driven by a ramping factor
            dn = 9001
            ops.node(dn, 50.0, 50.0)
            ops.fix(dn, 1, 1)
            ops.timeSeries("Linear", 2, "-factor", 1.0e6)
            ops.pattern("Plain", 2, 2)
            ops.load(dn, 123.0, 456.0)
        _newmark()
        for _ in range(6):
            assert ops.analyze(1, 0.1) == 0

    rec_a = _tmp("factor_a.csv")
    rec_b = _tmp("factor_b.csv")
    run(rec_a, False)          # (a) no series at all -> must complete
    run(rec_b, True)           # (b) ramping factor present elsewhere
    ha, ra = _read_csv(rec_a)
    hb, rb = _read_csv(rec_b)
    assert ra and any(abs(x) > 1.0 for x in ra[-1][1:]), "overlay produced ~no p"
    dmax = max(abs(a - b) for x, y in zip(ra, rb) for a, b in zip(x, y))
    assert dmax == 0.0, (
        f"overlay p changed when a ramping load factor was present: {dmax} "
        "(overlay must ignore the domain load factor)")

    # (c) no spurious notice on a clean run (subprocess captures C++ stderr)
    body = f"""
{_boot_src()}
nid,eles,top,L,nely = column()
ops.pattern("LadrunoPorousOverlay",77,"-region",*eles,{_ov_src()},
            "-drained",*top)
ops.timeSeries("Constant",1); ops.pattern("Plain",1,1)
for t in top: ops.load(t,0.0,-5.0)
newmark()
for _ in range(3):
    assert ops.analyze(1,0.1)==0
print("RUN_OK")
"""
    out, rc = _run_child(body)
    assert rc == 0, f"clean factor run failed in child (rc={rc}):\n{out}"
    assert "RUN_OK" in out
    assert "IGNORED" not in out.upper() or "TIMESERIES" not in out.upper(), (
        f"a TimeSeries-ignored notice was emitted on a series-free run:\n{out}")
    return "no-series run OK; ramping factor changes nothing (bit-exact); no notice"


# ==========================================================================
# 3. revertToLastCommit consistency
# ==========================================================================
def test_revert_to_last_commit_consistency():
    """A step that fails convergence is reverted by the integrator
    (Domain::revertToLastCommit). Because the overlay advances its fluid ONLY in
    the Domain::commit hook, a failed step never reaches the hook -> the fluid is
    not advanced and there is nothing to double-count. After the retry succeeds
    the p-history - including dpCommitted_, the rate-form reference consumed at
    the NEXT step - must equal a clean run.

    Scenario (linear model, so the retry reproduces the clean state to machine
    precision): march 5 clean steps; on step 6 force a NormDispIncr(maxIter=1)
    failure (linear Newton needs 2 iters), let the domain auto-revert, retry with
    maxIter=25; then take a clean step 7. Compare the record at t6 AND t7 to a
    clean 7-step run.
    """
    def run_clean(rec):
        nid, eles, top, L, nely = _column()
        ops.pattern("LadrunoPorousOverlay", 77, "-region", *eles, *_OV,
                    "-drained", *top, "-record", rec, 1)
        ops.timeSeries("Constant", 1)
        ops.pattern("Plain", 1, 1)
        for t in top:
            ops.load(t, 0.0, -5.0)
        _newmark()
        for _ in range(7):
            assert ops.analyze(1, 0.1) == 0

    def run_revert(rec):
        nid, eles, top, L, nely = _column()
        ops.pattern("LadrunoPorousOverlay", 77, "-region", *eles, *_OV,
                    "-drained", *top, "-record", rec, 1)
        ops.timeSeries("Constant", 1)
        ops.pattern("Plain", 1, 1)
        for t in top:
            ops.load(t, 0.0, -5.0)
        _newmark()
        for _ in range(5):
            assert ops.analyze(1, 0.1) == 0
        ops.test("NormDispIncr", 1e-12, 1)      # starve iterations -> fail
        rc_fail = ops.analyze(1, 0.1)
        assert rc_fail != 0, "step 6 unexpectedly converged at maxIter=1"
        ops.test("NormDispIncr", 1e-12, 25)     # retry with budget
        assert ops.analyze(1, 0.1) == 0, "retry after revert failed"
        assert ops.analyze(1, 0.1) == 0         # clean step 7
        return rc_fail

    recC = _tmp("rev_clean.csv")
    recR = _tmp("rev_revert.csv")
    run_clean(recC)
    rc_fail = run_revert(recR)
    hC, rC = _read_csv(recC)
    hR, rR = _read_csv(recR)
    assert len(rR) == len(rC) == 7, (
        f"revert run wrote {len(rR)} rows vs clean {len(rC)} - a failed step "
        "leaked a fluid advance (extra commit row)")
    # t6 (reverted+retried) and t7 (uses dpCommitted_ from t6) must match clean
    for k in (-2, -1):
        d = max(abs(a - b) for a, b in zip(rC[k], rR[k]))
        assert d == 0.0, (
            f"row t={rR[k][0]:.3f} differs from clean by {d:.3e} after "
            "revert+retry (rate-form / dpCommitted_ inconsistency)")
    return f"forced-fail rc={rc_fail}; rows=7=7; t6 & t7 bit-identical to clean"


# ==========================================================================
# 4. region-validation fatals <A-13> (each in a fresh subprocess)
# ==========================================================================
def _boot_src():
    """Source of the child bootstrap + shared helpers (numpy-free)."""
    return f"""import os, sys
_D = r{_DIST!r}
os.environ["PATH"] = _D + os.pathsep + os.environ.get("PATH","")
try: os.add_dll_directory(_D)
except Exception: pass
sys.path.insert(0, _D)
for _m in ("opensees","openseespy","openseespy.opensees"): sys.modules.pop(_m, None)
import opensees as ops
assert os.path.normcase(os.path.dirname(ops.__file__))==os.path.normcase(_D)
def column(nely=6, L=6.0, ndf=2):
    ops.wipe(); ops.model("basic","-ndm",2,"-ndf",ndf)
    nid={{}}; k=1
    for j in range(nely+1):
        for i in range(2):
            ops.node(k,float(i),L*j/nely); nid[(i,j)]=k; k+=1
    ops.nDMaterial("ElasticIsotropic",1,1e4,0.2,2.0)
    eles=[]
    for j in range(nely):
        et=j+1
        ops.element("quad",et,nid[(0,j)],nid[(1,j)],nid[(1,j+1)],nid[(0,j+1)],1.0,"PlaneStrain",1)
        eles.append(et)
    for j in range(nely+1):
        for i in range(2):
            ops.fix(nid[(i,j)],1,(1 if j==0 else 0))
    top=[nid[(0,nely)],nid[(1,nely)]]
    return nid,eles,top,L,nely
def newmark():
    ops.constraints("Transformation"); ops.numberer("RCM"); ops.system("ProfileSPD")
    ops.test("NormDispIncr",1e-8,25); ops.algorithm("Newton")
    ops.integrator("Newmark",0.6,0.3025); ops.analysis("Transient")
"""


def _ov_src():
    return ('"-Kf",2.2e7,"-rhoF",1000.0,"-perm",1e-6,1e-6,'
            '"-poro",0.4,"-moduli",1e4,0.2')


def _run_child(body):
    """Run `python -S <tmpfile>` with the given source. Returns (output, rc)."""
    fd, path = tempfile.mkstemp(suffix=".py", dir=_SCRATCH)
    os.close(fd)
    with open(path, "w") as f:
        f.write(body)
    proc = subprocess.run([sys.executable, "-S", "-u", path],
                          capture_output=True, text=True, timeout=180)
    return proc.stdout + proc.stderr, proc.returncode


# child scenarios: (name, body, expected-substring). Each must exit NONZERO.
def _fatal_cases():
    ov = _ov_src()
    boot = _boot_src()
    cases = []

    # (a) unsupported cell type in the region (a 2-node Truss)
    cases.append(("unsupported_cell", f"""
{boot}
nid,eles,top,L,nely = column()
ops.uniaxialMaterial("Elastic",50,1e6)
ops.element("Truss",999,nid[(0,0)],nid[(0,1)],1.0,50)
ops.pattern("LadrunoPorousOverlay",1,"-region",*(eles+[999]),{ov},"-drained",*top)
newmark()
rc = ops.analyze(1,0.1)
print("ANALYZE_RC", rc)
sys.exit(0 if rc==0 else 7)
""", "unsupported cell"))

    # (b) empty / omitted -drained
    cases.append(("empty_drained", f"""
{boot}
nid,eles,top,L,nely = column()
ops.pattern("LadrunoPorousOverlay",1,"-region",*eles,{ov})
print("SHOULD_NOT_REACH")
sys.exit(0)
""", "-drained"))

    # (c) two overlays claiming the same element
    # overlay 2 goes inert (region element already claimed by overlay 1); the
    # inert overlay's onDomainCommit returns -1, so the first commit fails.
    cases.append(("overlapping_overlays", f"""
{boot}
nid,eles,top,L,nely = column()
ops.pattern("LadrunoPorousOverlay",1,"-region",*eles,{ov},"-drained",*top)
ops.pattern("LadrunoPorousOverlay",2,"-region",eles[0],{ov},"-drained",nid[(0,1)])
newmark()
rc = ops.analyze(1,0.1)
print("ANALYZE_RC", rc)
sys.exit(0 if rc==0 else 7)
""", "already claimed"))

    # (d) a monolithic LadrunoUP element in the region
    cases.append(("up_in_region", f"""
{boot}
ops.wipe(); ops.model("basic","-ndm",2,"-ndf",3)
for i,(x,y) in enumerate([(0,0),(1,0),(1,1),(0,1)],1): ops.node(i,float(x),float(y))
ops.nDMaterial("ElasticIsotropic",1,1e4,0.2,2.0)
ops.element("LadrunoUP",1,1,2,3,4,1,"-Kf",5000.0,"-poro",0.4,"-rhoF",1.0,"-perm",1e-4,1e-4)
ops.fix(1,1,1,0); ops.fix(2,1,1,0)
ops.pattern("LadrunoPorousOverlay",7,"-region",1,{ov},"-drained",3,4)
ops.timeSeries("Constant",1); ops.pattern("Plain",1,1)
ops.load(3,0.0,-5.0,0.0)
ops.constraints("Transformation"); ops.numberer("RCM"); ops.system("UmfPack")
ops.test("NormDispIncr",1e-8,25); ops.algorithm("Newton")
ops.integrator("Newmark",0.6,0.3025); ops.analysis("Transient")
rc = ops.analyze(1,0.1)
print("ANALYZE_RC", rc)
sys.exit(0 if rc==0 else 7)
""", "monolithic u-p"))

    # (e1) unknown flag -substep (E7.3b removal, named)
    cases.append(("substep_removed", f"""
{boot}
nid,eles,top,L,nely = column()
ops.pattern("LadrunoPorousOverlay",1,"-region",*eles,{ov},"-drained",*top,"-substep",5)
print("SHOULD_NOT_REACH")
sys.exit(0)
""", "-substep was REMOVED"))

    # (e2) a generic unknown flag
    cases.append(("unknown_flag", f"""
{boot}
nid,eles,top,L,nely = column()
ops.pattern("LadrunoPorousOverlay",1,"-region",*eles,{ov},"-drained",*top,"-bogus",5)
print("SHOULD_NOT_REACH")
sys.exit(0)
""", "unknown flag '-bogus'"))

    # (f) missing required -moduli
    cases.append(("missing_moduli", f"""
{boot}
nid,eles,top,L,nely = column()
ops.pattern("LadrunoPorousOverlay",1,"-region",*eles,
            "-Kf",2.2e7,"-rhoF",1000.0,"-perm",1e-6,1e-6,"-poro",0.4,"-drained",*top)
print("SHOULD_NOT_REACH")
sys.exit(0)
""", "-moduli"))

    # (bonus) -fluidUpdate explicit is fatal-NYI until P3b
    cases.append(("fluidupdate_explicit_nyi", f"""
{boot}
nid,eles,top,L,nely = column()
ops.pattern("LadrunoPorousOverlay",1,"-region",*eles,{ov},"-drained",*top,
            "-fluidUpdate","explicit")
print("SHOULD_NOT_REACH")
sys.exit(0)
""", "-fluidUpdate explicit is NYI"))
    return cases


def test_region_validation_fatals():
    """<A-13> + parser hygiene. Each bad configuration runs in its OWN
    subprocess (a snapshot/parse fatal must not poison the rest of the suite),
    and must (i) exit NONZERO and (ii) print the specific diagnostic substring.
    Parse-time fatals raise OpenSeesError (nonzero exit); snapshot-time fatals
    make the overlay inert -> the first analyze fails because the inert overlay's
    onDomainCommit returns -1.
    """
    results = []
    for name, body, substr in _fatal_cases():
        out, rc = _run_child(body)
        assert rc != 0, (
            f"[{name}] expected NONZERO exit but got {rc}; output:\n{out}")
        assert substr in out, (
            f"[{name}] expected substring {substr!r} not found; output:\n{out}")
        assert "SHOULD_NOT_REACH" not in out, (
            f"[{name}] parser accepted an invalid pattern; output:\n{out}")
        results.append(name)
    return "fatal (nonzero+msg): " + ", ".join(results)


# ==========================================================================
# 5. pInit steady / hydrostatic vs closed form
# ==========================================================================
def test_pinit_hydrostatic_and_steady_closed_form():
    """Node-exact closed-form checks of the two analytic pInit modes.

      (a) HYDROSTATIC: -pInit hydrostatic $gw $zw sets p = gw*(zw - z) clipped at
          0. Read it with a single tiny-dt (1e-9) zero-external-load step so the
          record row == the initial committed field (diffusion frozen). Tested
          with zw at the column top (no clipping) AND zw at mid-height (upper
          half clipped to 0). Node-exact.
      (b) STEADY: -staticMode steady + -pInit steady solves H*p = f_seep with the
          drained set held at 0. Under a gravity -fluidBody this is the node-exact
          linear seepage profile gw*(z_drain - z) (drained top -> p(top)=0,
          p(base)=gw*L). Node-exact for the P1 basis.
    """
    gw, g, L, nely = 9810.0, 9.81, 6.0, 6

    # (a1) hydrostatic, zw = top (no clip)
    nid, eles, top, _, _ = _column(L=L, nely=nely)
    rec = _tmp("pi_hydro_full.csv")
    ops.pattern("LadrunoPorousOverlay", 77, "-region", *eles,
                "-Kf", 2.2e7, "-rhoF", 1000.0, "-perm", 1e-9, 1e-9,
                "-poro", 0.4, "-moduli", 1e4, 0.2, "-drained", *top,
                "-pInit", "hydrostatic", gw, L, "-record", rec, 1)
    _newmark()
    assert ops.analyze(1, 1e-9) == 0
    hdr, rows = _read_csv(rec)
    p = dict(zip(hdr, rows[-1][1:]))
    for (i, j), nt in nid.items():
        z = L * j / nely
        exp = max(gw * (L - z), 0.0)
        assert abs(p[nt] - exp) <= 1e-6 * gw + 1e-6, (
            f"hydro(full) node z={z}: p={p[nt]} exp={exp}")

    # (a2) hydrostatic, zw = mid (upper half clipped to 0)
    zw = 3.0
    nid, eles, top, _, _ = _column(L=L, nely=nely)
    rec = _tmp("pi_hydro_clip.csv")
    ops.pattern("LadrunoPorousOverlay", 77, "-region", *eles,
                "-Kf", 2.2e7, "-rhoF", 1000.0, "-perm", 1e-9, 1e-9,
                "-poro", 0.4, "-moduli", 1e4, 0.2, "-drained", *top,
                "-pInit", "hydrostatic", gw, zw, "-record", rec, 1)
    _newmark()
    assert ops.analyze(1, 1e-9) == 0
    hdr, rows = _read_csv(rec)
    p = dict(zip(hdr, rows[-1][1:]))
    nclip = 0
    for (i, j), nt in nid.items():
        z = L * j / nely
        exp = max(gw * (zw - z), 0.0)
        if exp == 0.0 and z > zw:
            nclip += 1
        assert abs(p[nt] - exp) <= 1e-6 * gw + 1e-6, (
            f"hydro(clip) node z={z}: p={p[nt]} exp={exp}")
    assert nclip > 0, "clip test exercised no clipped nodes"

    # (b) steady seepage -> node-exact linear profile
    nid, eles, top, _, _ = _column(L=L, nely=nely)
    rec = _tmp("pi_steady.csv")
    ops.pattern("LadrunoPorousOverlay", 77, "-region", *eles,
                "-Kf", 2.2e7, "-rhoF", 1000.0, "-perm", 1e-7, 1e-7,
                "-poro", 0.4, "-moduli", 1e4, 0.2, "-drained", *top,
                "-staticMode", "steady", "-pInit", "steady",
                "-fluidBody", 0.0, -g, "-record", rec, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("ProfileSPD")
    ops.test("NormDispIncr", 1e-10, 25)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    hdr, rows = _read_csv(rec)
    assert len(rows) == 1
    p = dict(zip(hdr, rows[-1][1:]))
    for j in range(nely + 1):
        z = L * j / nely
        exp = gw * (L - z)
        nt = nid[(0, j)]
        assert abs(p[nt] - exp) <= 1e-6 * gw + 1e-6, (
            f"steady node z={z}: p={p[nt]} exp={exp}")
    return "hydro(full+clip) node-exact; steady linear node-exact"


# ==========================================================================
# 6. serial DB round-trip (the 1.D gate)
# ==========================================================================
def _db_child_src(E, mod, presave, postsave, rich=False, frozen=False,
                  gw=9810.0, L=6.0, nely=6):
    """Child that runs ONE overlay-driven (self-contained) save/restore cycle in
    a fresh process and prints 'DMAX <value> REL <rel>'. Fresh process + single
    cycle removes FileDatastore cross-invocation state. NO external solid load
    (FileDatastore does not restore sibling Plain loads); hydrostatic pInit is the
    drive.

    `frozen`: FREEZE every solid DOF (Penalty) so u==0 exactly -> uSnapshot_ and
    Qt*Du vanish and the round-trip depends ONLY on pCommitted_/dpCommitted_
    (reliably bit-exact -> the p/dp serialization gate).
    `rich`: attach the maximal serialized surface (-layer/-alpha/-Ks/-fsL/
    -onRemoval)."""
    extra = ('"-alpha",0.9,"-Ks",3.6e8,"-layer",eles[0],eles[1],'
             '"-perm",2e-6,2e-6,"-poro",0.35,"-fsL","oedometric",'
             '"-onRemoval","drain",10.0,') if rich else ""
    fix_line = ("for i in range(2): ops.fix(nid[(i,j)],1,1)" if frozen else
                "for i in range(2): ops.fix(nid[(i,j)],1,(1 if j==0 else 0))")
    if frozen:
        anal = ('def run():\n'
                '    ops.constraints("Penalty",1e14,1e14); ops.numberer("RCM")\n'
                '    ops.system("ProfileSPD"); ops.test("NormDispIncr",1e-10,25)\n'
                '    ops.algorithm("Newton"); ops.integrator("Newmark",0.6,0.3025)\n'
                '    ops.analysis("Transient")')
    else:
        anal = 'def run(): newmark()'
    return f"""
{_boot_src()}
import tempfile, os as _os
FR = tempfile.mkdtemp(prefix="ovdb_")
def col2(E):
    ops.wipe(); ops.model("basic","-ndm",2,"-ndf",2)
    nid={{}}; k=1
    for j in range({nely}+1):
        for i in range(2):
            ops.node(k,float(i),{L}*j/{nely}); nid[(i,j)]=k; k+=1
    ops.nDMaterial("ElasticIsotropic",1,E,0.2,2.0)
    eles=[]
    for j in range({nely}):
        et=j+1
        ops.element("quad",et,nid[(0,j)],nid[(1,j)],nid[(1,j+1)],nid[(0,j+1)],1.0,"PlaneStrain",1)
        eles.append(et)
    for j in range({nely}+1):
        {fix_line}
    return nid,eles,[nid[(0,{nely})],nid[(1,{nely})]]
{anal}
def mk(rec):
    nid,eles,top = col2({E})
    ops.pattern("LadrunoPorousOverlay",77,"-region",*eles,"-Kf",2.2e7,"-rhoF",1000.0,
        "-perm",1e-6,1e-6,"-poro",0.4,"-moduli",{mod},0.2,"-drained",*top,
        {extra}"-pInit","hydrostatic",{gw},{L},"-record",rec,1)
def rows(p):
    with open(p) as f: L=f.read().strip().splitlines()
    return [[float(x) for x in r.split(",")] for r in L[1:]]
recC=_os.path.join(FR,"c.csv"); mk(recC); run()
for _ in range({presave}+{postsave}): assert ops.analyze(1,0.1)==0
recR=_os.path.join(FR,"r.csv"); mk(recR); run()
for _ in range({presave}): assert ops.analyze(1,0.1)==0
dbp=_os.path.join(FR,"db")
ops.database("File",dbp); ops.save(1); ops.wipe()
ops.database("File",dbp); ops.restore(1); run()
for _ in range({postsave}): assert ops.analyze(1,0.1)==0
rC=rows(recC); rR=rows(recR); tC={{round(r[0],6):r for r in rC}}
dmax=max(max(abs(a-b) for a,b in zip(r,tC[round(r[0],6)])) for r in rR)
scale=max(abs(x) for r in rC for x in r[1:]) or 1.0
print("DMAX", dmax, "REL", dmax/scale)
"""


def _child_dmax(out):
    for line in out.splitlines():
        if line.startswith("DMAX"):
            parts = line.split()   # "DMAX <dmax> REL <rel>"
            return float(parts[1]), float(parts[3])
    raise AssertionError(f"child produced no DMAX line:\n{out}")


def test_db_roundtrip_committed_p_and_dp():
    """1.D CORE gate (subprocess, isolated fresh process): committed p AND
    dpCommitted_ (the rate-form reference; dpCommitted_ is NOT re-derivable)
    round-trip BIT-EXACT through ops.database('File')/save/wipe/restore.

    Isolation on TWO axes: (i) FREEZE every solid DOF (Penalty) so u==0 exactly ->
    uSnapshot_ / Qt*Du vanish and the round-trip depends ONLY on the overlay's
    serialized fluid state (this cleanly separates the p/dp gate from any
    mass/solid-path restore defect (the #577 rho class that test 6c now gates)); (ii) a fresh child process runs a single
    save/restore cycle, removing FileDatastore cross-invocation state. The fluid
    diffuses from a hydrostatic pInit under a rich serialized surface (-layer,
    per-layer perm/poro, -alpha, -Ks, -fsL oedometric, -onRemoval drain).
    Continuous 6 steps vs restart 3+3; restart steps 4-6 must be byte-identical -
    verified reproducible across repeated fresh processes (always 0.0)."""
    out, rc = _run_child(_db_child_src(1e4, 1e4, 3, 3, rich=True, frozen=True))
    assert rc == 0, f"frozen-solid DB child failed (rc={rc}):\n{out}"
    dmax, rel = _child_dmax(out)
    assert dmax == 0.0, (
        f"frozen-solid DB restart differs from continuous by {dmax:.3e} "
        f"(rel {rel:.2e}) - committed p / dpCommitted_ did not round-trip")
    return "frozen-solid restart 4-6 bit-identical (p & dpCommitted_ serialize)"


def test_db_roundtrip_rich_config_survives_and_active():
    """1.D robustness: a rich DEFORMABLE overlay config restores WITHOUT crashing
    and stays ACTIVE afterwards (in-process). Build (layer, oedometric, onRemoval,
    hydrostatic pInit, alpha/Ks), march 3, save, wipe, restore, rebuild the
    analysis, march 3 more; assert crash-free restore, the overlay keeps recording
    (3 rows), elements survive, and p stays finite (no NaN / no 1e88 blow-up).
    A gentle pInit (gw=100) keeps this robustness gate mild; exact deformable
    continuity is the HARD GATE below (6c; upstream rho bug fixed in #577)."""
    gw, L, nely = 100.0, 6.0, 6
    rec = _tmp("db_rich.csv")
    nid, eles, top, _, _ = _column(L=L, nely=nely)
    ops.pattern("LadrunoPorousOverlay", 77, "-region", *eles,
                "-Kf", 2.2e7, "-rhoF", 1000.0, "-perm", 1e-6, 1e-6,
                "-poro", 0.4, "-moduli", 1e4, 0.2, "-drained", *top,
                "-alpha", 0.9, "-Ks", 3.6e8,
                "-layer", eles[0], eles[1], "-perm", 2e-6, 2e-6, "-poro", 0.35,
                "-fsL", "oedometric", "-onRemoval", "drain", 10.0,
                "-pInit", "hydrostatic", gw, L, "-record", rec, 1)
    _newmark()
    for _ in range(3):
        assert ops.analyze(1, 0.1) == 0
    dbp = _tmp("db_rich_ds")
    ops.database("File", dbp)
    ops.save(1)
    ops.wipe()
    ops.database("File", dbp)
    ops.restore(1)
    assert sorted(ops.getEleTags()) == sorted(eles), "elements lost on restore"
    _newmark()
    for _ in range(3):
        assert ops.analyze(1, 0.1) == 0, "overlay inert / analyze failed post-restore"
    hdr, rows = _read_csv(rec)
    assert len(rows) == 3, f"overlay stopped recording post-restore ({len(rows)} rows)"
    scale = 1000.0 * gw
    for r in rows:
        for v in r[1:]:
            assert v == v and abs(v) < scale, f"post-restore p out of regime: {v}"
    return "rich config restore crash-free; overlay active; p finite"


def test_db_roundtrip_deformable():
    """HARD GATE (promoted from xfail 2026-07-15, per its own XPASS contract):
    a DEFORMABLE overlay model's transient DB restart reproduces the
    continuous run across SEVERAL fresh child processes.

    History: shipped at P1 (#576) as a KNOWN-FAIL pinning what looked like an
    upstream FileDatastore/Domain restore corruption (non-deterministic across
    fresh processes, restored-operator eigenvalues going negative, all
    python-visible restored state exact). ROOT CAUSE found + fixed in #577:
    quad/tri elements' ELEMENT-LEVEL rho was never serialized and the broker
    ctor left it uninitialized -> garbage heap value overriding the
    (correctly restored) material density via the `if (rho==0)` gate ->
    corrupt mass matrix on every DB/MPI-restored element. The overlay was
    exonerated (its p/dp/uSnapshot round-trip was verified exact in diverging
    runs); it merely amplified the wrong-M response through its +Q*p forces.

    The multi-fresh-process structure is retained deliberately - it is the
    regression net for exactly this class of uninitialized-serialization bug
    (a single in-process run can land on benign heap garbage and pass)."""
    rels = []
    for _ in range(6):
        out, rc = _run_child(_db_child_src(1e4, 1e4, 3, 3, rich=False))
        if rc != 0:
            rels.append(float("inf"))     # restart-induced instability
            continue
        try:
            _, rel = _child_dmax(out)
        except AssertionError:
            rels.append(float("inf"))
            continue
        rels.append(rel)
    worst = max(rels)
    assert worst <= 0.01, (
        f"deformable DB restart diverges from continuous - worst rel={worst:.2e} "
        f"over 6 fresh processes (rels={['%.2e' % r for r in rels]}). "
        "The #577 class of bug (unserialized element state corrupting the "
        "restored operator) has REGRESSED - check element sendSelf/recvSelf "
        "completeness before suspecting the overlay.")
    return (f"deformable restart reproduces continuous in 6/6 fresh processes "
            f"(worst rel={worst:.2e}; upstream rho bug fixed in #577)")


# ==========================================================================
# 7. -record everyN semantics
# ==========================================================================
def test_record_everyN_and_header():
    """everyN=2 writes every OTHER commit; the header row is the region node-tag
    set. Six committed steps -> three rows; header (sans the leading 'time') is a
    permutation of the region node tags."""
    nid, eles, top, L, nely = _column()
    rec = _tmp("everyN2.csv")
    ops.pattern("LadrunoPorousOverlay", 77, "-region", *eles, *_OV,
                "-drained", *top, "-record", rec, 2)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for t in top:
        ops.load(t, 0.0, -5.0)
    _newmark()
    for _ in range(6):
        assert ops.analyze(1, 0.1) == 0
    hdr, rows = _read_csv(rec)
    assert len(rows) == 3, f"everyN=2 over 6 commits -> {len(rows)} rows (want 3)"
    region_nodes = set()
    for (i, j) in nid:
        region_nodes.add(nid[(i, j)])
    assert set(hdr) == region_nodes, (
        f"header {sorted(hdr)} != region node set {sorted(region_nodes)}")
    # everyN=1 baseline: 6 commits -> 6 rows
    nid, eles, top, L, nely = _column()
    rec1 = _tmp("everyN1.csv")
    ops.pattern("LadrunoPorousOverlay", 77, "-region", *eles, *_OV,
                "-drained", *top, "-record", rec1, 1)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for t in top:
        ops.load(t, 0.0, -5.0)
    _newmark()
    for _ in range(6):
        assert ops.analyze(1, 0.1) == 0
    _, rows1 = _read_csv(rec1)
    assert len(rows1) == 6, f"everyN=1 over 6 commits -> {len(rows1)} rows (want 6)"
    return "everyN=2 -> 3 rows, everyN=1 -> 6 rows, header == region node set"


# ==========================================================================
# standalone runner (python -S tests/test_ladruno_overlay_framework.py)
# ==========================================================================
# (name, callable, is_xfail)
_TESTS = [
    ("1  bit-constancy / once-per-commit", test_bit_constancy_once_per_commit_all_integrators, False),
    ("2  factor-ignored semantics",         test_factor_ignored_semantics, False),
    ("3  revertToLastCommit consistency",   test_revert_to_last_commit_consistency, False),
    ("4  region-validation fatals",         test_region_validation_fatals, False),
    ("5  pInit hydrostatic/steady",         test_pinit_hydrostatic_and_steady_closed_form, False),
    ("6a DB p/dp round-trip (frozen)",      test_db_roundtrip_committed_p_and_dp, False),
    ("6b DB rich config active",            test_db_roundtrip_rich_config_survives_and_active, False),
    ("6c DB deformable restart reproduces", test_db_roundtrip_deformable, False),
    ("7  -record everyN semantics",         test_record_everyN_and_header, False),
]


def _main():
    import traceback
    n_pass = n_fail = n_xfail = n_xpass = 0
    failures = []
    for label, fn, is_xfail in _TESTS:
        try:
            info = fn()
            if is_xfail:
                n_xpass += 1
                print(f"[XPASS] {label}  (known-fail now PASSES - promote it!)")
            else:
                n_pass += 1
                print(f"[PASS ] {label}" + (f"  -- {info}" if info else ""))
        except BaseException as e:  # SystemExit from a child etc.
            if is_xfail:
                n_xfail += 1
                print(f"[xfail] {label}  -- {type(e).__name__}: {str(e)[:90]}")
            else:
                n_fail += 1
                failures.append(label)
                print(f"[FAIL ] {label}\n{traceback.format_exc()}")
    print("-" * 70)
    print(f"PASS={n_pass} FAIL={n_fail} xfail={n_xfail} XPASS={n_xpass} "
          f"(of {len(_TESTS)})")
    if failures:
        print("FAILED:", ", ".join(failures))
    # suite is green iff no unexpected failures (XPASS is a loud note, not a fail)
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(_main())
