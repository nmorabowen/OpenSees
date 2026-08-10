"""The `.ladruno` recorder must not drop an `-E` request in silence (TIMs 5/6).

Two independent silent failures were found while chasing the TIMs bearing decks:

  * the ORCHESTRATOR's drop path. `LadrunoRecorder::initElementSources` probes
    every element with `elem->setResponse(argv, argc, stream)`. A null response
    (the element does not know that token) used to fall out of the if/else-if
    chain with no final `else`, an empty bucket was skipped, and a request that
    matched nothing produced no HDF5 group and no message at all. A typo in `-E`
    was therefore indistinguishable from a working recorder — the analysis ran to
    completion and the result simply was not in the file.

  * the ELEMENT's token vocabulary. The vanilla tetrahedra parsed their response
    tokens with a bare `strcmp` that accepted only the PLURAL spellings
    ("stresses"/"strains"), while their neighbours (LadrunoBrick, FourNodeQuad)
    accept both. `-E stress` therefore recorded data for some elements of a mixed
    model and nothing for the tets — again silently, because of the drop path
    above. Both tets now route through `LadrunoResp::is` (SRC/element/
    LadrunoResponseTokens.h), which matches both spellings and does NOT change
    what is emitted (same responseID, same ResponseType labels, same width).

The gates here are deliberately observable WITHOUT h5py: stderr (the new
warnings) plus `eleResponse` (the element-side alias fingerprint). The recorder
is exercised end-to-end all the same — the warnings are emitted from the real
`initElementSources` pass on a real domain.
"""

import pytest

from _testbed import ops


pytestmark = [pytest.mark.zone_a]

_E = 1.0e7
_NU = 0.0
_EPS = 1.0e-3
_PENALTY = 1.0e18

# The two new orchestrator messages, matched on a stable fragment each.
_MATCHED_NOTHING = "matched no element in the model"
_NO_RESPONSE = "returned no response"

# --- tet10 geometry (same layout as tests/test_tet10_response_size.py) -------
_NGP10 = 4          # TenNodeTetrahedron::NumGaussPoints
_NCOMP = 6          # Voigt {xx, yy, zz, xy, yz, zx}
_VERTS = [
    (0.0, 0.0, 0.0),
    (1.0, 0.0, 0.0),
    (0.0, 1.0, 0.0),
    (0.0, 0.0, 1.0),
]
_EDGES = [(1, 2), (2, 3), (1, 3), (1, 4), (3, 4), (2, 4)]


def _tet10_coords():
    """Node tag -> (x, y, z) for the 10 nodes, in element connectivity order."""
    pts = {i + 1: _VERTS[i] for i in range(4)}
    for k, (a, b) in enumerate(_EDGES, start=5):
        pa, pb = pts[a], pts[b]
        pts[k] = tuple(0.5 * (pa[j] + pb[j]) for j in range(3))
    return pts


def _static_rig():
    ops.constraints("Penalty", _PENALTY, _PENALTY)
    ops.numberer("RCM")
    ops.system("BandGen")
    ops.test("NormDispIncr", 1.0e-12, 20, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")


def _build_bricks():
    """Two stdBrick cubes stacked in x, bottom face fixed, top face pushed down.

    Small on purpose: the point is the recorder's setup pass, not the mechanics.
    """
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("ElasticIsotropic", 1, _E, 0.25)

    tag = 0
    for ix in range(3):            # 3 layers of 4 nodes -> 2 stacked bricks
        for (x, y) in [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]:
            tag += 1
            ops.node(tag, x, y, float(ix))
            if ix == 0:
                ops.fix(tag, 1, 1, 1)

    ops.element("stdBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
    ops.element("stdBrick", 2, 5, 6, 7, 8, 9, 10, 11, 12, 1)

    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for t in (9, 10, 11, 12):
        ops.sp(t, 3, -1.0e-3)


def _build_tet10():
    """One tet10 under the uniform-strain field u = (eps*x, 0, 0)."""
    pts = _tet10_coords()
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("ElasticIsotropic", 1, _E, _NU)
    for tag, (x, y, z) in pts.items():
        ops.node(tag, x, y, z)
    ops.element("TenNodeTetrahedron", 1, *range(1, 11), 1)

    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for tag, (x, _y, _z) in pts.items():
        ops.sp(tag, 1, _EPS * x)
        ops.sp(tag, 2, 0.0)
        ops.sp(tag, 3, 0.0)


def test_bogus_token_warns_instead_of_dropping_silently(tmp_path, capfd):
    """`-E banana`: the run still succeeds, but it must SAY that nothing matched.

    Pre-change this produced no HDF5 group, no message, and exit code 0 — the
    exact failure mode that made a mistyped request look like a working recorder.
    """
    _build_bricks()
    ops.recorder("ladruno", str(tmp_path / "out.ladruno"), "-E", "banana")
    _static_rig()
    rc = ops.analyze(1)
    ops.remove("recorders")
    err = capfd.readouterr().err

    assert rc == 0, f"the bogus -E request must not break the analysis (rc={rc})"
    assert _MATCHED_NOTHING in err, (
        "`-E banana` matched nothing and the recorder said nothing — the "
        f"silent-drop path is back. stderr tail={err[-800:]!r}")
    assert "banana" in err, (
        "the warning must name the offending request so the user can find the "
        f"typo. stderr tail={err[-800:]!r}")
    # the per-class aggregate must also fire, and exactly once per class (never
    # once per element — there are two stdBricks in this model)
    assert err.count(_NO_RESPONSE) == 1, (
        "the per-class 'returned no response' line must be emitted ONCE for the "
        f"stdBrick class, not once per element. stderr tail={err[-800:]!r}")


@pytest.mark.parametrize("token,plural", [("stress", "stresses"),
                                          ("strain", "strains")])
def test_tet10_accepts_the_singular_token(token, plural):
    """The alias fingerprint: singular now answers, with the plural's payload.

    Pre-change `setResponse` returned null for the singular spelling, so
    `eleResponse` came back empty and the recorder dropped the element.
    """
    _build_tet10()
    _static_rig()
    assert ops.analyze(1) == 0, "tet10 uniform-strain patch step failed"

    singular = ops.eleResponse(1, token)
    assert len(singular) == _NCOMP * _NGP10, (
        f"eleResponse(1, {token!r}) returned {len(singular)} values; expected "
        f"{_NCOMP * _NGP10} (6 components x {_NGP10} Gauss points). The singular "
        f"spelling is not reaching the {plural!r} branch.")

    # the alias must NOT change what is emitted: identical to the plural
    assert singular == ops.eleResponse(1, plural), (
        f"{token!r} and {plural!r} must return the SAME payload — an alias may "
        f"never change the response the element emits.")


def test_valid_request_emits_no_drop_warning(tmp_path, capfd):
    """Guard against warning spam: a fully answered request must stay quiet.

    Single-class model, single valid token — nothing is dropped, so neither of
    the two new lines may appear.
    """
    _build_tet10()
    ops.recorder("ladruno", str(tmp_path / "ok.ladruno"), "-E", "stress")
    _static_rig()
    rc = ops.analyze(1)
    ops.remove("recorders")
    err = capfd.readouterr().err

    assert rc == 0, f"the valid -E request must not break the analysis (rc={rc})"
    assert _MATCHED_NOTHING not in err, (
        f"'stress' IS answered by tet10 now; the matched-nothing warning must "
        f"not fire. stderr tail={err[-800:]!r}")
    assert _NO_RESPONSE not in err, (
        f"no element was dropped, so the per-class warning must not fire. "
        f"stderr tail={err[-800:]!r}")
