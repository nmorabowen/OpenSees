"""Self-tests for the MPCO_Ladruno test harness (runs today, no OpenSees build).

These validate the *harness itself* against a synthetic spec-conformant file:
  - L2 validator accepts a conformant file and rejects targeted corruptions
  - the reader + canonical normalizers extract the values we put in
  - L3 basis reconstruction maps quad Gauss points to the right global coords

When the C++ recorder lands, the same validator/normalizers run on real output and the
parity gate (test_parity.py, added in Phase 3) diffs against legacy .mpco.

Run:  & C:/Users/nmora/venv/opensees_venv/Scripts/python.exe -m pytest -q
"""

from __future__ import annotations

import numpy as np
import pytest

import ladruno_basis as lb
import ladruno_format as lf
import make_synthetic


@pytest.fixture()
def synth(tmp_path):
    p = str(tmp_path / "synthetic.ladruno")
    make_synthetic.build(p)
    return p


# ---- L2: validator ----------------------------------------------------------


def test_validator_accepts_conformant_file(synth):
    problems = lf.validate(synth)
    assert problems == [], f"synthetic file should be conformant, got: {problems}"


def test_validator_catches_bad_generator(synth, tmp_path):
    import h5py

    with h5py.File(synth, "a") as f:
        f["INFO"].attrs["GENERATOR"] = np.bytes_(b"STKO")
    problems = lf.validate(synth)
    assert any("GENERATOR" in p for p in problems)


def test_validator_catches_column_count_mismatch(synth):
    import h5py

    with h5py.File(synth, "a") as f:
        # break the NUM_COLUMNS == sum(MULT*NUM_COMP) invariant
        f["MODEL_STAGE[0]/RESULTS/ON_ELEMENTS/stress/156-FourNodeQuad[0:0]"].attrs[
            "NUM_COLUMNS"
        ] = 99
    problems = lf.validate(synth)
    assert any("NUM_COLUMNS" in p for p in problems)


def test_validator_catches_gp_param_shape(synth):
    import h5py

    with h5py.File(synth, "a") as f:
        g = f["MODEL_STAGE[0]/MODEL/ELEMENTS/156-FourNodeQuad[0]"]
        del g["QUADRATURE/GP_PARAM"]
        g["QUADRATURE"].create_dataset("GP_PARAM", data=np.zeros((3, 2)))  # wrong nGP
    problems = lf.validate(synth)
    assert any("GP_PARAM" in p for p in problems)


# ---- reader + normalizers ---------------------------------------------------


def test_reader_info_and_stage(synth):
    with lf.LadrunoReader(synth) as r:
        assert r.info()["GENERATOR"] == "MPCO_Ladruno"
        assert r.stages() == ["MODEL_STAGE[0]"]
        assert r.stage_attrs("MODEL_STAGE[0]")["KIND"] == "static"


def test_normalize_nodal(synth):
    with lf.LadrunoReader(synth) as r:
        nd = lf.normalize_nodal(r, "MODEL_STAGE[0]")
    # node 1, Ux, step 0 -> first entry of (k+1)*1e-3*arange(1..10) = 1e-3*1
    assert nd[(1, "DISPLACEMENT:Ux", 0)] == pytest.approx(1e-3)
    assert nd[(5, "DISPLACEMENT:Uy", 1)] == pytest.approx(2e-3 * 10)


def test_normalize_element_quad_stress(synth):
    with lf.LadrunoReader(synth) as r:
        ed = lf.normalize_element(r, "MODEL_STAGE[0]")
    # quad 10, gp 2, no section/fiber, sigma_yy, step 0.
    # row data = arange(12); gp2 block starts at col 6 -> [6,7,8]=[xx,yy,xy]; yy=7
    assert ed[(10, 2, -1, -1, "stress:sigma_yy", 0)] == pytest.approx(7.0)


def test_normalize_element_fiber_multiplicity(synth):
    with lf.LadrunoReader(synth) as r:
        ed = lf.normalize_element(r, "MODEL_STAGE[0]")
    # beam 20, gp0 has 3 fibers (multiplicity 3) -> fiber ids 0,1,2 at cols 0,1,2
    # data = arange(6)*10 -> fiber1 of gp0 = col1 = 10
    assert ed[(20, 0, 1, 1, "section.fiber.stress:sigma", 0)] == pytest.approx(10.0)
    # gp1 fiber2 -> col 5 = 50
    assert ed[(20, 1, 1, 2, "section.fiber.stress:sigma", 0)] == pytest.approx(50.0)


# ---- L3: basis reconstruction ----------------------------------------------


def test_reconstruct_quad_gauss_points(synth):
    """Quad on the unit square: GP at +/-1/sqrt3 -> global 0.5 +/- 0.5/sqrt3."""
    with lf.LadrunoReader(synth) as r:
        g = r.element_groups("MODEL_STAGE[0]")["156-FourNodeQuad[0]"]
        nid, coords = r.nodes("MODEL_STAGE[0]")
    tag2row = {int(t): i for i, t in enumerate(nid)}
    conn = g["connectivity"][0, 1:]  # node tags of element 10
    Xc = coords[[tag2row[int(t)] for t in conn]]

    half = 0.5
    d = 0.5 / np.sqrt(3.0)
    expected = np.array(
        [[half - d, half - d], [half + d, half - d],
         [half + d, half + d], [half - d, half + d]]
    )
    for j in range(g["num_gp"]):
        xg = lb.reconstruct_global(
            Xc, g["family"], g["topology"], g["order"], g["param_domain"], g["gp_param"][j]
        )
        assert np.allclose(xg, expected[j], atol=1e-12), f"gp {j}: {xg} != {expected[j]}"


def test_reconstruct_bernstein_line_endpoints():
    """Sanity: a linear Bernstein line still interpolates its endpoints."""
    Xc = np.array([[0.0, 0.0], [2.0, 0.0]])
    x0 = lb.reconstruct_global(Xc, "bernstein", "line", (1,), "[0,1]", np.array([0.0]))
    x1 = lb.reconstruct_global(Xc, "bernstein", "line", (1,), "[0,1]", np.array([1.0]))
    assert np.allclose(x0, [0.0, 0.0])
    assert np.allclose(x1, [2.0, 0.0])
