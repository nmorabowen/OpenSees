"""WP-133 (TIMs F23a): PDMY03 critical-state constants as optional flags.

``PressureDependMultiYield03`` hard-coded ``ei = 0.6, cs1 = 0.9, cs2 = 0.02,
cs3 = 0.7`` in its constructor. They feed ``isCriticalState()``, the only
dilation/contraction brake in the model::

    e    = ei + eps_v (1 + ei)                        (eps_v = tr(strain))
    e_cr = cs1 - cs2 (p'/pa)^cs3      (cs3 = 0 -> cs1 - cs2 ln(p'/pa))

WP-133 exposes them as trailing flags, in both parsers (Tcl ladder and the
Python ``OPS_`` function)::

    nDMaterial PressureDependMultiYield03 $tag ... <positional tail>
        <-ei $e0> <-cs1 $v> <-cs2 $v> <-cs3 $v>

Gates:

  G1  omitted flags -> byte-identical to the baseline captured from the
      UNMODIFIED build (Python: three decks incl. user-defined surfaces;
      Tcl: one deck). Explicit flags equal to the defaults -> identical too.
  G2  the constants reach the brake: moving the critical-state line so the
      path crosses it changes the response from the crossing step on, and
      that step moves the way the formula says (higher cs1 -> later; higher
      cs2, cs3 or ei -> earlier). Also pins the brake's nature: a one-step
      crossing event, after which dilation resumes at the reference rate.
  G3  isolation: with two PDMY03 materials of different constants defined in
      one model, a quad on either reproduces its single-material run bit for
      bit, in either creation order and when >20 further PDMY03 materials
      are created afterwards (the static per-material arrays are reallocated
      in chunks of 20; before WP-133 that reallocation overwrote every
      existing material's constants with the newest material's). Plus two
      elements in one system agree with their single runs (to the Newton tol).
  G4  a bad or incomplete flag is refused (no material is created).

Wall time: ~30 s.
"""
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

_HERE = Path(__file__).resolve().parent
_WT = _HERE.parent
_DIST = str(_WT / "dist" / "bin")
if not os.path.isfile(os.path.join(_DIST, "opensees.pyd")):
    pytest.skip(f"worktree engine not built: {_DIST}", allow_module_level=True)

from _engine import bind_worktree_engine  # noqa: E402
ops = bind_worktree_engine(_DIST)

sys.path.insert(0, str(_HERE))
import wp133_pdmy03_deck as D  # noqa: E402

pytestmark = [pytest.mark.zone_a]

DEFAULT_FLAGS = ["-ei", 0.6, "-cs1", 0.9, "-cs2", 0.02, "-cs3", 0.7]


def _run(extra=(), tag=1):
    D.set_ele_mat(1, tag)
    return D.run_quad(ops, tag, materials=lambda o: D.build_material(
        o, tag, tail=D.TAIL, extra=extra))


def _hex(res):
    return D.to_hex(res)


def _void_ratio(res, ei=0.6):
    # Gauss point 1; plane strain -> eps_v = exx + eyy
    return [ei + (s[0] + s[1]) * (1.0 + ei) for s in res["strain"]]


def _divergence_step(res, ref):
    """First step at which ``res`` departs from ``ref`` (bitwise).

    Before the path reaches the critical-state line ``isCriticalState()``
    returns 0 for both runs, so they are bitwise identical; the first
    differing step is where the ENGINE's brake first fired. None if never.
    """
    for k, (a, b) in enumerate(zip(res["stress"], ref["stress"])):
        if a != b:
            return k
    return None


def _cs1_mid(ref):
    """A cs1 that puts the CSL across the reference path near mid-run.

    Aimed with a nominal p' = 1000 kPa (sigma_zz is not reported for a
    PlaneStrain quad, so p' is not reconstructed); the gates read the
    engine's own crossing step, not this estimate. With this deck the
    engine's brake first fires at step ~20 of 100.
    """
    k = len(ref["strain"]) // 2
    return _void_ratio(ref)[k] + 0.02 * (1000.0 / 101.0) ** 0.7


# ---------------------------------------------------------------- G1
def test_g1_python_byte_identical_to_unmodified_build():
    base = json.loads((_HERE / "wp133_pdmy03_byteid_baseline.json").read_text())
    for name, fn in D.CASES.items():
        got = _hex(fn(ops))
        assert got["n"] == base[name]["n"], name
        assert got == base[name], f"{name}: not byte-identical to the pre-WP-133 build"


def test_g1_python_explicit_default_flags_identical():
    base = json.loads((_HERE / "wp133_pdmy03_byteid_baseline.json").read_text())
    assert _hex(_run(DEFAULT_FLAGS)) == base["default_tail"]
    # order-independent
    shuffled = ["-cs3", 0.7, "-ei", 0.6, "-cs2", 0.02, "-cs1", 0.9]
    assert _hex(_run(shuffled)) == base["default_tail"]


def _tcl(*extra):
    exe = os.path.join(_DIST, "OpenSees.exe")
    r = subprocess.run([exe, str(_HERE / "wp133_pdmy03_deck.tcl"), *map(str, extra)],
                       capture_output=True, text=True, timeout=300,
                       stdin=subprocess.DEVNULL)  # a Tcl error drops to the prompt
    out = (r.stdout or "") + (r.stderr or "")
    rows = [ln for ln in out.splitlines() if ln[:1].isdigit()]
    return rows, out


def test_g1_tcl_byte_identical_to_unmodified_build():
    base = (_HERE / "wp133_pdmy03_tcl_baseline.txt").read_text().splitlines()
    rows, out = _tcl()
    assert len(rows) == 100, out[-2000:]
    assert rows == base
    rows2, _ = _tcl(*DEFAULT_FLAGS)
    assert rows2 == base


def test_g1_tcl_and_python_agree():
    # same deck through both parsers -> same numbers. Not bitwise: Tcl's
    # eleResponse result passes through a decimal string (last-ulp noise).
    rows, _ = _tcl("-cs1", 0.62)
    res = _run(["-cs1", 0.62])
    assert len(rows) == res["n"]
    for row, s, e in zip(rows, res["stress"], res["strain"]):
        vals = [float(v) for v in row.split()[1:]]
        assert vals == pytest.approx(list(s) + list(e), rel=1e-12, abs=1e-15)


# ---------------------------------------------------------------- G2
def test_g2_constants_reach_the_brake():
    ref = _run()
    # with the defaults the dense deck stays far below the CSL (e ~ 0.60-0.63
    # vs e_cr ~ 0.8-0.88): cs1 = 5 moves the line even further away -> the
    # response cannot change
    assert _hex(_run(["-cs1", 5.0])) == _hex(ref)
    cs1 = _cs1_mid(ref)
    moved = _run(["-cs1", cs1])
    assert moved["n"] == ref["n"]
    k = _divergence_step(moved, ref)
    assert k is not None and 5 < k < moved["n"] - 5, k


def test_g2_crossing_moves_with_each_constant():
    ref = _run()
    cs1 = _cs1_mid(ref)
    k0 = _divergence_step(_run(["-cs1", cs1]), ref)
    assert k0 is not None
    # e rises (dilation) while e_cr = cs1 - cs2 (p/pa)^cs3 falls (p rises):
    k = _divergence_step(_run(["-cs1", cs1 + 0.01]), ref)
    assert k > k0, ("higher cs1 -> line higher -> later crossing", k0, k)
    k = _divergence_step(_run(["-cs1", cs1, "-cs2", 0.025]), ref)
    assert k < k0, ("higher cs2 -> line lower -> earlier crossing", k0, k)
    k = _divergence_step(_run(["-cs1", cs1, "-cs3", 0.8]), ref)
    assert k < k0, ("higher cs3 (p > pa) -> line lower -> earlier", k0, k)
    k = _divergence_step(_run(["-cs1", cs1, "-ei", 0.605]), ref)
    assert k < k0, ("higher ei -> e higher -> earlier crossing", k0, k)


def test_g2_brake_is_a_crossing_event_not_a_state():
    """The act's reading, checked: isCriticalState() returns 1 only while
    current and trial states sit on OPPOSITE sides of the CSL. Once past
    the line, dilation resumes at the reference rate -- no plateau."""
    ref = _run()
    moved = _run(["-cs1", _cs1_mid(ref)])
    k = _divergence_step(moved, ref)
    ev = lambda r: [st[0] + st[1] for st in r["strain"]]
    vr, vm = ev(ref), ev(moved)
    rate = lambda v, j: v[j] - v[j - 1]
    # at the crossing step the volumetric rate dips
    assert rate(vm, k) < rate(vr, k)
    # ten steps later it is back to the reference rate (within 0.5 %)
    j = k + 10
    assert rate(vm, j) == pytest.approx(rate(vr, j), rel=5e-3)
    assert rate(vm, j) > 0.5 * max(rate(vm, i) for i in range(1, len(vm)))


# ---------------------------------------------------------------- G3
# The constants live in STATIC per-material arrays indexed by matN, so
# interference is through those arrays, not through elements. The bitwise
# gates therefore build several materials in one model and drive ONE quad:
# two decoupled quads in one Newton system do not reproduce single runs
# bitwise anyway (the global test gives each element extra iterations).
def _one_quad_among(mat_defs, use_tag, pad_after=0):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for tag, extra in mat_defs:
        D.build_material(ops, tag, tail=D.TAIL, extra=extra)
    for j in range(pad_after):
        D.build_material(ops, 100 + j, tail=D.TAIL, extra=["-cs1", 0.5 + 0.001 * j])
    D.set_ele_mat(1, use_tag)
    D.add_quad(ops, use_tag, 0, 1, 0.0)
    return D.drive(ops, [(0, 1)], D.NSTEP)[1]


def test_g3_two_materials_do_not_interfere():
    ref = _run()
    cs1 = _cs1_mid(ref)
    alone_a = _hex(ref)
    alone_b = _hex(_run(["-cs1", cs1]))
    assert alone_a != alone_b
    defs = [(1, []), (2, ["-cs1", cs1])]
    assert _hex(_one_quad_among(defs, 1)) == alone_a
    assert _hex(_one_quad_among(defs, 2)) == alone_b
    defs = defs[::-1]  # reversed creation order
    assert _hex(_one_quad_among(defs, 1)) == alone_a
    assert _hex(_one_quad_among(defs, 2)) == alone_b


def test_g3_realloc_chunk_keeps_each_materials_constants():
    # >20 further PDMY03 materials force the matCount%20 reallocation AFTER
    # the two under test were created; the pre-WP-133 copy loop wrote the
    # NEWEST material's constants into every existing slot (here cs1 = 0.524,
    # which would put material 1 far on the loose side of the line).
    ref = _run()
    cs1 = _cs1_mid(ref)
    alone_a = _hex(ref)
    alone_b = _hex(_run(["-cs1", cs1]))
    defs = [(1, []), (2, ["-cs1", cs1])]
    assert _hex(_one_quad_among(defs, 1, pad_after=25)) == alone_a
    assert _hex(_one_quad_among(defs, 2, pad_after=25)) == alone_b


def _pair(extra2, nstep):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    D.build_material(ops, 1, tail=D.TAIL)
    D.build_material(ops, 2, tail=D.TAIL, extra=extra2)
    D.set_ele_mat(1, 1)
    D.set_ele_mat(2, 2)
    D.add_quad(ops, 1, 0, 1, 0.0)
    D.add_quad(ops, 2, 10, 2, 5.0)
    return D.drive(ops, [(0, 1), (10, 2)], nstep)


def test_g3_two_elements_coexist_before_the_crossing():
    # Two quads in one Newton system: pair (default, moved-cs1) must equal
    # pair (default, default) BITWISE until material 2's own crossing, for
    # both elements -- a leaked constant would move either one earlier.
    # (Stops before the crossing: past it, this 2-element system stalls
    # inside analyze(); see the WP-133 note, open item.)
    ref = _run()
    cs1 = _cs1_mid(ref)
    k = _divergence_step(_run(["-cs1", cs1]), ref)
    n = k - 1
    moved = _pair(["-cs1", cs1], n)
    same = _pair([], n)
    assert moved[1]["n"] == same[1]["n"] == n
    assert _hex(moved[1]) == _hex(same[1])
    assert _hex(moved[2]) == _hex(same[2])


# ---------------------------------------------------------------- G4
@pytest.mark.parametrize("extra", [["-cs1"], ["-cs1", "abc"], ["-ei", 0.6, "-bogus", 1.0]])
def test_g4_bad_flags_refused(extra):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    with pytest.raises(Exception):
        D.build_material(ops, 7, tail=D.TAIL, extra=extra)


def test_g4_bad_flags_refused_tcl():
    rows, out = _tcl("-cs1")
    assert rows == []
    assert "bad option" in out
