"""WP-152 (tension cutoff) test helper: a single-brick material-point driver with a
prescribed NORMAL-strain path, and the three element paths the oracle fixture uses.

Unit stdBrick, rollers on the three negative faces, every positive-face DOF prescribed
from three Path series (eps_xx, eps_yy, eps_zz; OpenSees convention, tension positive):
zero free equations, a homogeneous strain, eight identical Gauss points.

The path: N0 stage-0 (elastic) steps of isotropic compression to about p0 = 2 kPa, the
stage flip, then the plastic-stage increments (compression-positive, the oracle's
convention, negated here).  The oracle expectations are in
tests/data/wp152_oracle_paths.json (Ladruno_files/testbed/sanisand_tension_cutoff/
tc_oracle.py --fixture, CPython 3.11 + scipy), computed from the SAME post-flip state
this driver records (record_flip_state())."""
from __future__ import annotations

import json
import os

from _testbed import ops
import sanisand_replay as sr

HERE = os.path.dirname(os.path.abspath(__file__))
FIXTURE = os.path.join(HERE, "data", "wp152_oracle_paths.json")
FLIP_STATE = os.path.join(HERE, "data", "wp152_flip_state.json")

PARAMS = list(sr.CAMPAIGN_PARAMS)          # the TIMs campaign set (the oracle's CAMPAIGN)
BASE = (129, 0, 1, 1.0e-7, 1.0e-4, "-flipAlphaIn", "init", "-Pmin", 0.0101,
        "-maxSubsteps", 2000, "-Presidual", 0.0, "-honorTolR", 0)
R1 = ("-sasHFloor", 1.0, "-sasReseatHyst", 1.0, "-sasSoftCap", 0.5)
CUTOFF = ("-sasTensionCutoff", 0.5, 1.0)
XY = [(0., 0.), (1., 0.), (1., 1.), (0., 1.)]
N0 = 5                    # stage-0 steps
EV0 = 1.05e-5             # stage-0 volumetric compression (about p0 = 2 kPa)
D = 1.0e-5                # plastic-stage strain increment


def paths():
    """name -> list of compression-positive increments (6, Voigt; normal strains only)."""
    iso = [[-D, -D, -D, 0, 0, 0]] * 50 + [[D, D, D, 0, 0, 0]] * 60
    te = [[-D, 0, 0, 0, 0, 0]] * 60 + [[D, 0, 0, 0, 0, 0]] * 90
    cyc = []
    for _ in range(4):
        cyc += [[-D, -D, -D, 0, 0, 0]] * 60 + [[D, D, D, 0, 0, 0]] * 60
    return {"iso": iso, "te": te, "cyc": cyc}


def _series(tag0, incs, ev0=EV0):
    """Cumulative element-convention strains per step (stage 0, then the plastic
    increments negated), one Path series per normal component."""
    e0 = -ev0 / 3.0
    hist = [[e0 * (k + 1) / N0] * 3 for k in range(N0)]
    cur = list(hist[-1])
    for de in incs:
        cur = [cur[i] - de[i] for i in range(3)]
        hist.append(list(cur))
    for i in range(3):
        # padded with the last value once more: PathSeries returns 0.0 AT its last point
        # (getFactor needs index + 1 < size), so the final step would see zero strain
        ops.timeSeries("Path", tag0 + i, "-dt", 1.0, "-values", 0.0, *[h[i] for h in hist], hist[-1][i])
    return hist


def build(incs, opts, mat_tag=1, ev0=EV0):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for k in range(2):
        for j, (x, y) in enumerate(XY):
            ops.node(4 * k + j + 1, x, y, float(k))
    ops.nDMaterial("LadrunoSANISAND", mat_tag, *PARAMS, *opts)
    ops.element("stdBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, mat_tag)
    for k in range(2):
        for j, (x, y) in enumerate(XY):
            ops.fix(4 * k + j + 1, 1 if x == 0. else 0, 1 if y == 0. else 0, 1 if k == 0 else 0)
    hist = _series(1, incs, ev0)
    for i, ts in enumerate((1, 2, 3)):
        ops.pattern("Plain", 10 + i, ts)
        for k in range(2):
            for j, (x, y) in enumerate(XY):
                n = 4 * k + j + 1
                if i == 0 and x == 1.:
                    ops.sp(n, 1, 1.0)
                if i == 1 and y == 1.:
                    ops.sp(n, 2, 1.0)
                if i == 2 and k == 1:
                    ops.sp(n, 3, 1.0)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-12, 10, 0)
    ops.algorithm("Linear")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    return hist


def mresp(name):
    return list(ops.eleResponse(1, "material", 1, name))


def sig_comp():
    """Gauss point 1 stress, compression positive (the oracle's convention)."""
    return [-x for x in mresp("stress")]


def stage0(mat_tag=1):
    ops.updateMaterialStage("-material", mat_tag, "-stage", 0)
    for _ in range(N0):
        assert ops.analyze(1) == 0
    ops.updateMaterialStage("-material", mat_tag, "-stage", 1)


def run(incs, opts, per_step=None, ev0=EV0):
    """Drive the path; returns (flip stress, records).  Each record: rc, stress
    (compression positive), sasStats."""
    build(incs, opts, ev0=ev0)
    stage0()
    flip = sig_comp()
    out = []
    for k in range(len(incs)):
        rc = ops.analyze(1)
        rec = dict(k=k, rc=rc, sigma=sig_comp(), sas=mresp("sasStats"))
        if per_step:
            per_step(k, rec)
        out.append(rec)
        if rc != 0:
            break
    return flip, out


def e_flip():
    """The void ratio after stage 0: e_init - (1 + e_init) tr(eps), compression positive."""
    return PARAMS[2] - (1.0 + PARAMS[2]) * EV0


def record_flip_state(path=FLIP_STATE):
    """The post-flip state the oracle fixture must start from (run on any build: stage 0
    is elastic and untouched by WP-152)."""
    build([[0, 0, 0, 0, 0, 0]], BASE)
    stage0()
    st = dict(sigma=sig_comp(), e=e_flip(), note="stdBrick, stage 0 isotropic compression "
              f"eps_v = {EV0} in {N0} steps, then the flip (alpha = alpha_in = 0)")
    ops.wipe()
    json.dump(st, open(path, "w"), indent=1)
    return st


def with_max_substeps(opts, n):
    """`opts` with its -maxSubsteps value replaced (a cheap, deterministic code 9)."""
    o = list(opts)
    o[o.index("-maxSubsteps") + 1] = n
    return tuple(o)


# ------------------------------------------------------------------ review item 9
# A two-brick oedometric column solved by NEWTON (free equations, a force test):
# bottom face on rollers, every lateral face on rollers (u_x = 0 on x = 0, 1; u_y = 0
# on y = 0, 1), the four MIDDLE nodes (z = 1) free in z and loaded by the weight W
# (downward), the four TOP nodes (z = 2) prescribed in z. Opening the top separates the
# upper brick while the lower one carries W: the middle nodes are an equilibrium the
# material state machine must let Newton find, through the entry and the re-contact.
COL_W = 4.0            # total weight on the middle nodes (kPa x 1 m^2)


def column_series(top_u):
    """Top-face z displacement history (element convention: up positive)."""
    ops.timeSeries("Path", 1, "-dt", 1.0, "-values", 0.0, *top_u, top_u[-1])


def build_column(top_u, opts, mat_tag=1):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    tag = 0
    for k in range(3):
        for j, (x, y) in enumerate(XY):
            tag = 4 * k + j + 1
            ops.node(tag, x, y, float(k))
            ops.fix(tag, 1, 1, 1 if k == 0 else 0)     # lateral rollers everywhere, base fixed in z
    ops.nDMaterial("LadrunoSANISAND", mat_tag, *PARAMS, *opts)
    ops.element("stdBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, mat_tag)
    ops.element("stdBrick", 2, 5, 6, 7, 8, 9, 10, 11, 12, mat_tag)
    column_series(top_u)
    ops.pattern("Plain", 10, 1)
    for n in (9, 10, 11, 12):
        ops.sp(n, 3, 1.0)
    ops.timeSeries("Constant", 2)
    ops.pattern("Plain", 11, 2)
    for n in (5, 6, 7, 8):
        ops.load(n, 0.0, 0.0, -COL_W / 4.0)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormUnbalance", 1.0e-9, 30, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")


def eresp(ele, name):
    return list(ops.eleResponse(ele, "material", 1, name))


# ------------------------------------------------------------------ plane strain
def build_plane(hist_xy, opts, mat_tag=1):
    """A unit `quad` (PlaneStrain) with rollers on x = 0 and y = 0 and the right / top
    edges prescribed from two Path series of cumulative strains (tension positive)."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for j, (x, y) in enumerate(XY):
        ops.node(j + 1, x, y)
        ops.fix(j + 1, 1 if x == 0. else 0, 1 if y == 0. else 0)
    ops.nDMaterial("LadrunoSANISAND", mat_tag, *PARAMS, *opts)
    ops.element("quad", 1, 1, 2, 3, 4, 1.0, "PlaneStrain", mat_tag)
    for i in range(2):
        ops.timeSeries("Path", 1 + i, "-dt", 1.0, "-values", 0.0, *[h[i] for h in hist_xy], hist_xy[-1][i])
        ops.pattern("Plain", 10 + i, 1 + i)
        for j, (x, y) in enumerate(XY):
            if i == 0 and x == 1.:
                ops.sp(j + 1, 1, 1.0)
            if i == 1 and y == 1.:
                ops.sp(j + 1, 2, 1.0)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-12, 10, 0)
    ops.algorithm("Linear")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
