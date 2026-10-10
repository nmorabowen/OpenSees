"""ADR-159 -- the smoothed (C1) mortar contact law: `contact ... -mortar -smoothN <g0> [-smoothT <r>]`.

R0.8 (pile-contact): Newton could not carry the 3D mortar interface through active-set changes
(lift-off behind a laterally loaded pile, a moving slip front), because the shipped normal law
p = min(0, epsN*gbar) has a kink at gbar = 0 and the cohesive friction cone jumps from c to 0 at
lift-off. -smoothN replaces the kink by a C1 quadratic onset over |gbar| < g0 and fades the friction
traction in over the same band (a C1 smoothstep); -smoothT rounds the stick/slip corner of the return
map over |rho - cap| < r*cap. Off => byte-identical (no flag => setMortarSmoothing never called).

The law (x = -epsN*gbar, S = epsN*g0):  P = 0 (x <= -S); (x+S)^2/(4S) (|x| < S); x (x >= S).
Oracle: Ladruno_implementation/contact_prototypes/proto_adr159_smooth_normal.py (law + FD tangent).

Gates:
  (law)   a stiff block on springs, force control, pressed to a pressure inside the band: the
          equilibrium matches the smoothed law exactly; pressed beyond the band it equals the shipped
          penalty law; pulled off past the band the contact carries NOTHING (no tension).
  (a)     block lift-off with a cohesive bond under a lateral load, pressed then pulled off past
          separation under force control: Newton converges in <= 8 iterations per step through the
          separation and the master reaction is exactly zero once open.
  (b)     the ADR-158 solid creased roof (mu = 0.2, shared ridge, -consistanttan) still matches the
          analytic descent with a small band, in <= 6 iterations per step.
  (c)     band -> 0 on a closed case: the state converges to the shipped answer (monotone, and
          identical once the band no longer reaches the penetration).
  (db)    database save -> wipe -> restore reproduces a smoothed contact exactly (wire v5), and a
          committed v4 stream written by the pre-ADR-159 binary restores on this one (data/adr159_v4_db).
  (ref)   every refusal is named and the test matches the NAME: no -mortar, -tie, an augmenting
          contact, -smoothT without friction, out-of-range values, -soft/-visc, and a 2D pair (a
          handle-time FATAL). A pre-ADR-159 binary also raises on these inputs ("unexpected token"),
          so matching the message is what makes the gate discriminate (the ADR-155 pattern).
"""
import math
import os
import shutil
import tempfile

import pytest

from _testbed import ops

import test_adr158_mortar_consistanttan_multipair as R158

pytestmark = [pytest.mark.zone_a]

KN = 1.0e6       # epsN
KD = 1.0e3       # spring per top node (4 top nodes)
A = 1.0          # block footprint (unit square)


def _block(flags=(), coh=0.0, mu=0.0, delta=0.0):
    """A stiff unit brick resting on a fixed, larger master quad (mortar, -augment never), exactly
    (delta = 0) or overlapping it by delta, each top node on (x, y, z) springs KD. Returns (top nodes, bottom nodes, master nodes)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    m = []
    for t, (x, y) in enumerate([(-0.5, -0.5), (1.5, -0.5), (1.5, 1.5), (-0.5, 1.5)], start=1):
        ops.node(t, x, y, 0.0)
        ops.fix(t, 1, 1, 1)
        m.append(t)
    bot, top = [11, 12, 13, 14], [15, 16, 17, 18]
    for t, (x, y) in zip(bot, [(0, 0), (1, 0), (1, 1), (0, 1)]):
        ops.node(t, float(x), float(y), -delta)
        ops.node(t + 4, float(x), float(y), 1.0 - delta)
    ops.nDMaterial("ElasticIsotropic", 1, 1.0e11, 0.2)
    ops.element("stdBrick", 1, *bot, *top, 1)
    ops.contactSurface(1, "-master", 4, *m)
    ops.contactSurface(2, "-slave-segments", 4, *bot)
    opts = ["-mortar", "-epsN", KN, "-epsT", KN, "-outward", 0.0, 0.0, 1.0, "-augment", "never"]
    if mu > 0.0:
        opts += ["-mu", mu]
    if coh > 0.0:
        opts += ["-cohesion", coh]
    ops.contact(1, 1, 2, *opts, *flags)
    ops.uniaxialMaterial("Elastic", 1, KD)
    for t in top:
        ops.node(100 + t, *ops.nodeCoord(t))
        ops.fix(100 + t, 1, 1, 1)
        ops.element("zeroLength", 200 + t, 100 + t, t, "-mat", 1, 1, 1, "-dir", 1, 2, 3)
    _solver()
    return top, bot, m


def _solver():
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormUnbalance", 1.0e-7, 40, 0)
    ops.algorithm("Newton")


def _law(x, S):
    if S <= 0.0:
        return max(x, 0.0)
    if x <= -S:
        return 0.0
    if x >= S:
        return x
    return (x + S) ** 2 / (4.0 * S)


def _w_analytic(F, g0):
    """Downward displacement w of the rigid block: F = A*P(KN*w) + 4*KD*w (bisection)."""
    S = KN * g0
    lo, hi = -10.0, 10.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if A * _law(KN * mid, S) + 4 * KD * mid - F > 0.0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def _apply(top, fz, fx=0.0, steps=1):
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t in top:
        ops.load(t, fx / 4.0, 0.0, -fz / 4.0)
    ops.integrator("LoadControl", 1.0 / steps)
    ops.analysis("Static")
    its = []
    for _ in range(steps):
        if ops.analyze(1) != 0:
            return its, False
        its.append(ops.testIter())
    return its, True


def _w(bot):
    return -sum(ops.nodeDisp(n, 3) for n in bot) / len(bot)


# ------------------------------------------------------------------------------------- (law)
@pytest.mark.parametrize("F", [0.5, 2.0, 50.0, -0.5])
def test_adr159_block_equilibrium_matches_smoothed_law(F):
    """In the band (F = 0.5, 2: x/S = -0.1..0.8), beyond it (F = 50: the shipped ramp), and pulled
    off past it (F = -0.5: w = F/(4 KD), the contact carries nothing)."""
    g0 = 1.0e-6                                    # S = 1.0 (pressure units)
    top, bot, m = _block(("-smoothN", g0))
    _, ok = _apply(top, F)
    assert ok
    w = _w(bot)
    wa = _w_analytic(F, g0)
    assert w == pytest.approx(wa, rel=1e-6, abs=1e-14), f"w {w:.6e} != smoothed law {wa:.6e}"
    if F >= 50.0:
        assert wa == pytest.approx(_w_analytic(F, 0.0), rel=1e-12)   # outside the band: shipped
    ops.reactions()
    rz = sum(ops.nodeReaction(t, 3) for t in m)
    if F < 0.0:
        assert KN * w < -KN * g0                   # open beyond the band
        assert rz == 0.0, f"tension through an open contact: master Rz = {rz}"


# --------------------------------------------------------------------------------------- (a)
def test_adr159_block_liftoff_newton_bounded():
    """Cohesive bond (c = 2) under a lateral pull, pressed then pulled off past separation, force
    control, 12 steps: every step converges in <= 8 iterations (measured 3, 3, 3, 3, 4, 4, 6, 2, 1,
    1, 1, 1), the master carries exactly nothing once the block is open beyond the band, and the
    lateral pull then rides the springs alone. A rigid block flips all four nodes at once, so the
    shipped law converges here too (1-2 iterations): this gate pins that the smoothed law and its
    tangent get THROUGH separation, not the pile-deck stall (that is the R3 repro, ADR-159 section 4)."""
    g0, coh, H = 1.0e-4, 2.0, 0.5
    top, bot, m = _block(("-consistanttan", "-smoothN", g0, "-smoothT", 0.1), coh=coh)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t in top:
        ops.load(t, H / 4.0, 0.0, -1.0 / 4.0)      # pattern 1: lateral H + press 1 (x lambda)
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    ops.loadConst("-time", 0.0)
    ops.timeSeries("Linear", 2)
    ops.pattern("Plain", 2, 2)
    for t in top:
        ops.load(t, 0.0, 0.0, 1.0 / 4.0)           # pull up: net vertical goes -1 -> +1.4
    ops.integrator("LoadControl", 2.4 / 12)
    its = []
    for _ in range(12):
        assert ops.analyze(1) == 0, f"Newton failed through lift-off after {its}"
        its.append(ops.testIter())
    assert max(its) <= 8, f"iterations per step {its} (want <= 8)"
    w = _w(bot)
    assert KN * w < -KN * g0, f"block not open beyond the band (w = {w})"
    ops.reactions()
    rx = sum(ops.nodeReaction(t, 1) for t in m)
    rz = sum(ops.nodeReaction(t, 3) for t in m)
    assert rx == 0.0 and rz == 0.0, f"open contact carries Rx {rx}, Rz {rz}"
    ux = sum(ops.nodeDisp(t, 1) for t in top) / 4.0
    assert ux == pytest.approx(H / (4 * KD), rel=1e-6)


# --------------------------------------------------------------------------------------- (b)
def test_adr159_creased_roof_analytic_with_band():
    """The ADR-158 solid roof (mu = 0.2, shared ridge, full slip) with -consistanttan and a band far
    inside the penetration: the analytic descent still holds (2e-4) in <= 6 iterations per step."""
    its, ok, surf, springs = R158._roof(flags=("-consistanttan", "-smoothN", 1.0e-8), **R158.SOLID)
    assert ok, f"Newton failed after {its}"
    assert max(its) <= 6, f"iterations per step {its}"
    kd = 1.0e3
    _, fz = R158._descent_analytic(len(surf), len(springs), R158.SOLID["pz"], kd, 1.0e6,
                                   R158.SOLID["mu"], R158.SOLID["delta"], 1.0)
    got = len(surf) * R158.SOLID["pz"] + sum(kd * ops.nodeDisp(n, 3) for n in springs)
    assert got == pytest.approx(fz, rel=2e-4)


# --------------------------------------------------------------------------------------- (c)
def test_adr159_band_to_zero_converges_to_shipped():
    """Closed case (the block starts 3e-6 into the foundation -- engaged from the first evaluation
    under both laws -- and is pressed by F = 2 with a lateral pull and a cohesive bond; penetration
    ~5e-6): the smoothed state -> the shipped state as g0 -> 0, monotone, and identical to the
    tolerance once the band no longer reaches the penetration (g0 = 1e-6)."""
    def state(flags):
        top, _, _ = _block(("-consistanttan",) + flags, coh=2.0, delta=3.0e-6)
        _, ok = _apply(top, 2.0, fx=0.5)
        assert ok, f"{flags}: Newton failed"
        return [v for t in top for v in ops.nodeDisp(t)]
    ref = state(())
    scale = max(abs(v) for v in ref)
    errs = []
    for g0 in (3.0e-5, 1.0e-5, 1.0e-6):
        got = state(("-smoothN", g0))
        errs.append(max(abs(a - b) for a, b in zip(got, ref)) / scale)
    assert errs[0] > errs[1] > errs[2], f"not monotone: {errs}"
    assert errs[0] > 1e-2 and errs[2] < 1e-6, f"errors {errs}"


# -------------------------------------------------------------------------------------- (db)
def test_adr159_database_roundtrip():
    """save -> wipe -> restore -> analyze reproduces a smoothed cohesive contact exactly (the v5
    mortar record carries smoothN/smoothT; a record that dropped them would restore the shipped
    law and land elsewhere)."""
    flags = ("-consistanttan", "-smoothN", 1.0e-6, "-smoothT", 0.2)
    top, bot, _ = _block(flags, coh=2.0)
    _, ok = _apply(top, 2.0, fx=0.5)
    assert ok
    ref = [ops.nodeDisp(t) for t in top]
    d = tempfile.mkdtemp(prefix="adr159_db_")
    top, bot, _ = _block(flags, coh=2.0)
    ops.database("File", os.path.join(d, "db"))
    ops.save(1)
    ops.wipe()
    ops.database("File", os.path.join(d, "db"))
    ops.restore(1)
    _solver()                                      # analysis objects are not part of the database
    _, ok = _apply(top, 2.0, fx=0.5)
    assert ok
    got = [ops.nodeDisp(t) for t in top]
    assert got == ref
    # and the shipped law lands elsewhere on the same model (the smoothing is live)
    top, bot, _ = _block(("-consistanttan",), coh=2.0)
    _, ok = _apply(top, 2.0, fx=0.5)
    assert ok and [ops.nodeDisp(t) for t in top] != ref


# ----------------------------------------------------------------------------------- (db v4)
V4_DB = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "adr159_v4_db")


def _v4_model():
    """The model inside the committed v4 fixture: the cohesive block, -augment never -consistanttan,
    no smoothing (the pre-ADR-159 binary cannot write it)."""
    return _block(("-consistanttan",), coh=2.0)


def write_v4_fixture(dirpath):
    """Regenerate data/adr159_v4_db on the PRE-ADR-159 binary (ladruno 3144e19ba), never on this one
    (which writes v5). From tests/:  python -c "import test_adr159_mortar_smooth_contact as t;
    t.write_v4_fixture('data/adr159_v4_db')"   (the README in that folder says the same)."""
    _v4_model()
    ops.database("File", os.path.join(dirpath, "db"))
    ops.save(1)
    ops.wipe()


def test_adr159_restores_a_v4_stream(tmp_path):
    """A database written by the pre-ADR-159 binary (definitions stream v4: 64 slots, no smoothN /
    smoothT) restores on this binary as an unsmoothed contact: the restored analysis equals a freshly
    built twin bit for bit. A v4 read that left the two new slots undefined, or refused v4, fails."""
    files = sorted(f for f in os.listdir(V4_DB) if f.startswith("db."))
    assert any(f.startswith("db.VECs.64.") for f in files), files   # v4 defs record (v5 is 66)
    for f in files:
        shutil.copy(os.path.join(V4_DB, f), str(tmp_path / f))
    top, _, _ = _v4_model()
    _, ok = _apply(top, 2.0, fx=0.5)
    assert ok
    ref = [ops.nodeDisp(t) for t in top]
    ops.wipe()
    ops.database("File", str(tmp_path / "db"))
    ops.restore(1)
    _solver()                                      # analysis objects are not part of the database
    _, ok = _apply(top, 2.0, fx=0.5)
    assert ok
    assert [ops.nodeDisp(t) for t in top] == ref


# ------------------------------------------------------------------------------------- (ref)
def _facets():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, (x, y) in enumerate([(0, 0), (1, 0), (1, 1), (0, 1)], start=1):
        ops.node(t, float(x), float(y), 0.0)
        ops.node(t + 10, float(x), float(y), 0.0)
    ops.contactSurface(1, "-master", 4, 1, 2, 3, 4)
    ops.contactSurface(2, "-slave-segments", 4, 11, 12, 13, 14)


def _contact_fails(opts):
    _facets()
    try:
        ops.contact(1, 1, 2, *opts)
    except Exception:
        return True
    return False


def _refused(capfd, why, opts):
    """Refused FOR THE NAMED REASON (the ADR-155 pattern): the base binary raises on these inputs too,
    but with "unexpected token '-smoothN'", which carries neither "ADR-159" nor `why`."""
    _facets()
    capfd.readouterr()
    with pytest.raises(Exception):
        ops.contact(1, 1, 2, *opts)
    err = capfd.readouterr()
    text = err.out + err.err
    assert "ADR-159" in text and why in text, text


@pytest.mark.parametrize("opts,why", [
    (("auto", "-smoothN", 1e-4), "-smoothN/-smoothT are -mortar options"),                  # NTS
    (("-mortar", "-epsN", 1e6, "-tie", "-augment", "never", "-smoothN", 1e-4),
     "-smoothN/-smoothT do not apply to -tie"),
    (("-mortar", "-epsN", 1e6, "-smoothN", 1e-4),                                           # commit
     "pure-penalty laws; add -augment never"),
    (("-mortar", "-epsN", 1e6, "-augment", "request", "-smoothN", 1e-4),
     "pure-penalty laws; add -augment never"),
    (("-mortar", "-epsN", 1e6, "-augment", "never", "-smoothT", 0.1),
     "this contact has no friction"),
    (("-mortar", "-epsN", 1e6, "-augment", "never", "-smoothN", 0.0), "-smoothN - need a band g0 > 0"),
    (("-mortar", "-epsN", 1e6, "-augment", "never", "-smoothN", -1e-4), "-smoothN - need a band g0 > 0"),
    (("-mortar", "-epsN", 1e6, "-mu", 0.3, "-augment", "never", "-smoothT", 1.0),
     "-smoothT - need a ratio 0 < r < 1"),
    (("-mortar", "-epsN", 1e6, "-mu", 0.3, "-augment", "never", "-smoothT", 0.0),
     "-smoothT - need a ratio 0 < r < 1"),
    (("-mortar", "-epsN", 1e6, "-augment", "never", "-visc", 1.0, "-smoothN", 1e-4),
     "not defined with -soft or -visc"),
    (("-mortar", "-epsN", 1e6, "-augment", "never", "-soft", 0.1, "-smoothN", 1e-4),
     "not defined with -soft or -visc"),
])
def test_adr159_refusals(capfd, opts, why):
    _refused(capfd, why, opts)


@pytest.mark.parametrize("flag", [("-smoothN", 1e-4), ("-mu", 0.3, "-smoothT", 0.1)])
def test_adr159_2d_pair_refused_at_handle(capfd, flag):
    """The smoothed law is wired to the 3D mortar lane only: a 2D pair draws the named handle-time
    FATAL (analyze < 0) instead of silently running the shipped law. The same 2D block with
    -augment never alone converges (the control)."""
    def block(extra):
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 2)
        ops.node(101, 0.0, 0.0)
        ops.node(102, 1.0, 0.0)
        ops.fix(101, 1, 1)
        ops.fix(102, 1, 1)
        ops.node(1, 0.0, -1e-4)
        ops.node(2, 1.0, -1e-4)
        ops.fix(1, 1, 0)
        ops.fix(2, 1, 0)
        ops.contactSurface(10, "-master", 2, 101, 102)
        ops.contactSurface(20, "-slave-segments", 2, 1, 2)
        ops.contact(1, 10, 20, "-mortar", "-epsN", 1e6, "-outward", 0.0, 1.0, "-augment", "never",
                    *extra)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(1, 0.0, -1.0)
        ops.load(2, 0.0, -1.0)
        _solver()
        ops.integrator("LoadControl", 1.0)
        ops.analysis("Static")
        try:
            return ops.analyze(1)
        except Exception:
            return -1
    assert block(()) == 0
    capfd.readouterr()
    assert block(flag) < 0
    err = capfd.readouterr()
    text = err.out + err.err
    assert "-smoothN/-smoothT are 3D-mortar options (ADR-159); this pair is 2D" in text, text


def test_adr159_accepts_the_recipe():
    assert not _contact_fails(("-mortar", "-epsN", 1e6, "-cohesion", 1.0, "-augment", "never",
                               "-smoothN", 1e-4, "-smoothT", 0.1))
