"""ASDPlasticMaterial3D emits real ResponseType column names for its per-GP
material responses recorded via `-E material.<token>`.

Context: recording the equivalent plastic strain of ASDPlasticMaterial3D through
a solid element was reported as producing an empty group. Root cause was the
recorder TOKEN, not the element: only a `material.`-prefixed token makes the
recorder iterate Gauss points and hit the element's material-forwarding branch
(both SixNodeTri/tri6n and TenNodeTetrahedron already have it). A bare
`-E eqpstrain` never reaches the material -> empty group.

A secondary papercut: ASDPlasticMaterial3D::setResponse emitted no
NdMaterialOutput/ResponseType XML, so the recorded columns landed under the
recorder's generic `C1..CN` fallback. This battery pins the fix:

  * material.eqpstrain -> a single column named "eqpstrain" per Gauss point,
  * material.pstrain   -> six columns "epsP11..epsP13" per Gauss point,
  * an unknown token (e.g. "plasticStrain", which is spelled "pstrain" on this
    material) records NOTHING (guarded) rather than a malformed bucket,
  * a bare (non `material.`) token records an empty group,
  * on BOTH the tri6n (plane strain wrapper) and TenNodeTetrahedron (3D) hosts.
"""
import glob
import os

import numpy as np
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

IV = "BackStress(TensorLinearHardeningFunction):YieldStress(ScalarLinearHardeningFunction):"
E, NU, SY = 70000.0, 0.3, 30.0


def _asdp(tag):
    """A VonMises ASDPlasticMaterial3D (perfect plasticity: zero hardening)."""
    ops.nDMaterial(
        "ASDPlasticMaterial3D", tag,
        "VonMises_YF", "VonMises_PF", "LinearIsotropic3D_EL", IV,
        "Begin_Model_Parameters",
        "YoungsModulus", E, "PoissonsRatio", NU,
        "ScalarLinearHardeningParameter", 0.0,
        "TensorLinearHardeningParameter", 0.0,
        "MassDensity", 0.0,
        "End_Model_Parameters",
        "Begin_Internal_Variables",
        "YieldStress", SY,
        "BackStress", 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        "End_Internal_Variables",
    )


def _has_asdp():
    """Skip cleanly if ASDPlasticMaterial3D isn't in this build."""
    try:
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 3)
        _asdp(1)
        return ops.nDMaterialType(1) is not None if hasattr(ops, "nDMaterialType") else True
    except Exception:
        return False


def _tet():
    ops.wipe(); ops.model("basic", "-ndm", 3, "-ndf", 3)
    nds = {1: (0, 0, 0), 2: (1, 0, 0), 3: (0, 1, 0), 4: (0, 0, 1),
           5: (.5, 0, 0), 6: (.5, .5, 0), 7: (0, .5, 0), 8: (0, 0, .5), 9: (.5, 0, .5), 10: (0, .5, .5)}
    for t, c in nds.items():
        ops.node(t, *map(float, c))
    for t in (1, 2, 3, 5, 6, 7):
        ops.fix(t, 1, 1, 1)
    top = (4, 8, 9, 10)
    for t in top:
        ops.fix(t, 1, 1, 0)
    _asdp(1)
    ops.element("TenNodeTetrahedron", 1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1)
    return "TenNodeTetrahedron", 4, [(t, 3, -0.02) for t in top]


def _tri6():
    ops.wipe(); ops.model("basic", "-ndm", 2, "-ndf", 2)
    nds = {1: (0, 0), 2: (1, 0), 3: (0, 1), 4: (.5, 0), 5: (.5, .5), 6: (0, .5)}
    for t, c in nds.items():
        ops.node(t, *map(float, c))
    for t in (1, 2, 4):
        ops.fix(t, 1, 1)
    for t in (3, 5, 6):
        ops.fix(t, 1, 0)
    _asdp(1)
    ops.nDMaterial("PlaneStrain", 2, 1)
    ops.element("tri6n", 1, 1, 2, 3, 4, 5, 6, 1.0, "PlaneStrain", 2)
    return "SixNodeTri", 3, [(t, 2, -0.02) for t in (3, 5, 6)]


def _drive(sps):
    ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1)
    for (nd, dof, val) in sps:
        ops.sp(nd, dof, val)
    ops.constraints("Penalty", 1e14, 1e14)
    ops.numberer("Plain"); ops.system("UmfPack")
    ops.test("NormDispIncr", 1e-8, 100, 0); ops.algorithm("Newton")
    ops.integrator("LoadControl", 0.05); ops.analysis("Static")
    assert ops.analyze(20) == 0


def _bucket(fn, result, h5py):
    with h5py.File(fn, "r") as h:
        got = {}

        def visit(name, obj):
            if isinstance(obj, h5py.Group) and f"ON_ELEMENTS/{result}/" in name and name.count("/") == 4:
                has = "DATA" in obj
                comp = None
                if "COLUMN_MAP" in obj and "COMP_NAMES" in obj["COLUMN_MAP"].attrs:
                    v = obj["COLUMN_MAP"].attrs["COMP_NAMES"]
                    comp = v[0].decode() if isinstance(v, np.ndarray) else str(v)
                data_max = float(np.max(np.abs(obj["DATA"][-1]))) if has else None
                got[name.rsplit("/", 1)[-1]] = (has, comp, data_max)

        h.visititems(visit)
    return got


@pytest.fixture(scope="module")
def _asdp_available():
    if not _has_asdp():
        pytest.skip("ASDPlasticMaterial3D not available in this build")


@pytest.mark.parametrize("builder", [_tet, _tri6], ids=["TenNodeTetrahedron", "tri6n"])
def test_asdp_material_response_tags(builder, tmp_path, _asdp_available):
    h5py = pytest.importorskip("h5py")
    os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

    # yield it, then confirm the element forwards the per-GP material request
    name, ngp, sps = builder()
    _drive(sps)
    eqp = ops.eleResponse(1, "material", 1, "eqpstrain")
    assert eqp and abs(eqp[0]) > 1e-4, f"{name}: expected plastic strain, got {eqp}"

    cases = {
        "material.eqpstrain": dict(present=True, ncol=1 * ngp,
                                   names=["eqpstrain"] * ngp),
        "material.pstrain": dict(present=True, ncol=6 * ngp,
                                 names=["epsP11,epsP22,epsP33,epsP12,epsP23,epsP13"] * ngp),
        "material.plasticStrain": dict(present=False),   # typo -> guarded empty
        "eqpstrain": dict(present=False),                # bare -> empty group
    }

    for tok, exp in cases.items():
        name2, ngp2, sps2 = builder()
        out = str(tmp_path / f"{name}_{tok.replace('.', '_')}.ladruno")
        ops.recorder("ladruno", out, "-E", tok)
        _drive(sps2)
        ops.wipe()
        files = glob.glob(out) or glob.glob(out + "*")
        assert files, f"{tok}: no .ladruno file"
        buckets = _bucket(files[0], tok, h5py)

        if not exp["present"]:
            assert not any(v[0] for v in buckets.values()), \
                f"{tok}: expected NO data bucket, got {buckets}"
            continue

        assert buckets, f"{tok}: expected a data bucket, got none"
        (has, comp, dmax) = next(iter(buckets.values()))
        assert has, f"{tok}: bucket has no DATA"
        rows = [r for r in (comp or "").splitlines() if r.strip()]
        assert rows == exp["names"], f"{tok}: comp names {rows} != {exp['names']}"
        assert dmax and dmax > 1e-4, f"{tok}: expected nonzero plastic data, got {dmax}"
