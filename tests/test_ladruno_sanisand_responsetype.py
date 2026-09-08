"""ResponseType component names on the seven fork material responses of
`LadrunoSANISAND` (WP-86d).

Before this change, `LadrunoSANISAND::setResponse` created all seven fork
responses (`substeps`, `implexError`, `avgImplexError`, `implexDetail`,
`implexRefusals`, `psi`/`stateParameter`, `yieldDistance`/`yieldFunction`) as
bare `new MaterialResponse(this, id, probe)` -- none of them called
`output.tag("ResponseType", ...)`, so any recorder (the classic XML element
recorder, and the fork's own `recorder ladruno`) fell back to the generic
`C1..Cn` column names (see `Ladruno_ElementResults.h`'s `"C" << i + 1`
fallback, and the sibling battery `test_asdplastic_response_tags.py` for the
same papercut on `ASDPlasticMaterial3D`). apeGmsh's reader has agreed on a
fixed set of column-name strings; this file pins that those exact strings,
in that exact order, are what both recorders now emit.

Mirrors the idiom used elsewhere in the tree (`FSAM::setResponse`,
`ASDConcrete3DMaterial::setResponse`):
`output.tag("NdMaterialOutput"); output.attr("matType", ...);
output.attr("matTag", ...); output.tag("ResponseType", "<name>"); ...;
output.endTag();` around each `new MaterialResponse(...)`.

Deck: the confine-first zero-free-DOF cube of `test_ladruno_sanisand.py`
(`_build_confined` / `_confine_leg`), `stdBrick`, same as
`test_ladruno_sanisand_responses.py`. The seven responses are read-only
diagnostics that exist on `LadrunoSANISAND` regardless of `-implex` (they are
plain `MaterialResponse` registrations in `setResponse`, gated on nothing but
the response name); one test below builds WITHOUT `-implex` and confirms the
IMPL-EX-flavoured responses still carry their ResponseType tags.
"""
import glob
import os
import re

import pytest

from _testbed import ops
from test_ladruno_sanisand import _build_confined, _confine_leg

pytestmark = [pytest.mark.zone_a]

# name -> expected ResponseType strings, in the exact order apeGmsh's reader
# expects (== the order LadrunoSANISAND::setResponse emits them == the order
# LadrunoSANISAND::getResponse fills the underlying Vector).
_EXPECTED = {
    "substeps": ["substeps_me", "substeps_capHit"],
    "implexError": ["implexError"],
    "avgImplexError": ["avgImplexError"],
    "implexDetail": [
        "implexDetail_total", "implexDetail_dev", "implexDetail_vol",
        "implexDetail_clampFired", "implexDetail_clampCount", "implexDetail_f",
    ],
    "implexRefusals": [
        "implexRefusals_total", "implexRefusals_signChange",
        "implexRefusals_control", "implexRefusals_companion",
    ],
    "psi": ["psi"],
    "yieldDistance": ["yieldDistance"],
}

_RESPONSETYPE_RE = re.compile(r"<ResponseType>([^<]*)</ResponseType>")


def _xml_response_types(path):
    """Every `<ResponseType>` value in the XML header, in document order,
    with the recorder's own `time` column (always first when `-time` isn't
    suppressed) filtered out -- it is not part of what we are pinning here."""
    with open(path, "r") as f:
        text = f.read()
    return [m for m in _RESPONSETYPE_RE.findall(text) if m != "time"]


@pytest.mark.parametrize("name", sorted(_EXPECTED))
def test_xml_recorder_header_names_the_components(name, tmp_path):
    """The classic `recorder Element -xml` header carries the exact,
    correctly-ordered ResponseType strings for each of the seven responses --
    not the generic C1..Cn fallback."""
    opts = ("-implex", "-maxSubsteps", 5000)
    _build_confined("LadrunoSANISAND", 1, opts)
    out = str(tmp_path / f"rt_{name}.xml")
    ops.recorder("Element", "-xml", out, "-ele", 1, "material", 1, name)
    _confine_leg(1)
    ops.wipe()   # flush/close the recorder before reading the file back

    files = glob.glob(out) or glob.glob(out + "*")
    assert files, f"{name}: no XML file produced at {out}"
    got = _xml_response_types(files[0])
    assert got == _EXPECTED[name], (name, got, _EXPECTED[name])


def test_implex_named_responses_exist_without_implex():
    """The IMPL-EX-flavoured responses are plain registrations in
    `setResponse` -- they must carry their ResponseType tags (and therefore
    be usable by a recorder) even on a deck that never passes `-implex`."""
    _build_confined("LadrunoSANISAND", 1, ())   # no -implex
    for name in ("implexError", "avgImplexError", "implexDetail",
                 "implexRefusals"):
        r = ops.eleResponse(1, "material", 1, name)
        assert len(r) == len(_EXPECTED[name]), (name, r)

    out_dir_case = "implexDetail"
    import tempfile
    tmp = tempfile.mkdtemp()
    out = os.path.join(tmp, "rt_no_implex.xml")
    _build_confined("LadrunoSANISAND", 1, ())   # fresh deck, still no -implex
    ops.recorder("Element", "-xml", out, "-ele", 1, "material", 1,
                 out_dir_case)
    _confine_leg(1)
    ops.wipe()
    files = glob.glob(out) or glob.glob(out + "*")
    assert files, "no XML file produced without -implex"
    got = _xml_response_types(files[0])
    assert got == _EXPECTED[out_dir_case], (got, _EXPECTED[out_dir_case])


def test_ladruno_recorder_column_map_names_the_components():
    """The fork's own `recorder ladruno` (HDF5) is the one apeGmsh reads.
    Same mechanism as `test_asdplastic_response_tags.py`: `ResponseType`
    tags land in each Gauss point's `COLUMN_MAP` `COMP_NAMES` attribute
    instead of the generic `C1..Cn` fallback (`Ladruno_ElementResults.h`,
    `OutputDescriptorStream::tag(name, value)` pushing onto `components`
    when `name == "ResponseType"`)."""
    h5py = pytest.importorskip("h5py")
    import numpy as np
    os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

    opts = ("-implex", "-maxSubsteps", 5000)
    # `-E material.<name>` (no explicit Gauss-point index) expands over EVERY
    # integration point of the element (LadrunoRecorder.cpp's `do_all_materials`
    # loop) -- stdBrick carries 8 (`NDMaterial *materialPointers[8]`,
    # SRC/element/brick/Brick.h), so COMP_NAMES holds one newline-separated
    # row PER Gauss point, each row the (possibly comma-joined) names for
    # that response -- same layout `test_asdplastic_response_tags.py` pins
    # for ASDPlasticMaterial3D.
    _NGP = 8

    checks = [
        ("material.psi", "psi", ["psi"] * _NGP),
        ("material.implexDetail", "implexDetail",
         [",".join(_EXPECTED["implexDetail"])] * _NGP),
    ]
    for tok, name, expect in checks:
        _build_confined("LadrunoSANISAND", 1, opts)
        import tempfile
        tmp = tempfile.mkdtemp()
        out = os.path.join(tmp, f"rt_{name}.ladruno")
        ops.recorder("ladruno", out, "-E", tok)
        _confine_leg(1)
        ops.wipe()

        files = glob.glob(out) or glob.glob(out + "*")
        assert files, f"{tok}: no .ladruno file"

        got_comp = []
        with h5py.File(files[0], "r") as h:
            def visit(hname, obj):
                if (isinstance(obj, h5py.Group)
                        and f"ON_ELEMENTS/{tok}/" in hname
                        and hname.count("/") == 4
                        and "COLUMN_MAP" in obj
                        and "COMP_NAMES" in obj["COLUMN_MAP"].attrs):
                    v = obj["COLUMN_MAP"].attrs["COMP_NAMES"]
                    comp = v[0].decode() if isinstance(v, np.ndarray) else str(v)
                    got_comp.append(comp)
            h.visititems(visit)

        assert got_comp, f"{tok}: no COLUMN_MAP/COMP_NAMES bucket found"
        rows = got_comp[0].splitlines()
        assert rows == expect, (tok, rows, expect)
