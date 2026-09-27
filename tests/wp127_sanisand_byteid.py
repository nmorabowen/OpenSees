"""WP-127 byte-identity decks for the SANISAND substep counters.

WP-127 adds per-instance ModifiedEuler counters and an optional substep trace
to `ManzariDafalias` (vanilla, `// Ladruno WP-127`) and `LadrunoSANISAND`.
They must not move a single computed bit.  This module runs a fixed set of
zero-free-DOF material-point decks and returns every committed stress /
back-stress / state vector as `float.hex` strings, so two builds can be
compared EXACTLY (not to a tolerance).

The pinned reference `wp127_sanisand_byteid_baseline.json` beside this file
was captured with the PRE-WP-127 binary (tree c03a1bd4b, built unmodified in
the WP-127 worktree before any source edit) by

    python -S <bootstrap> wp127_sanisand_byteid.py --write

and `test_ladruno_sanisand_replay_counters.py::test_counters_are_byte_identical`
compares the current build against it.  Regenerate it ONLY for a deliberate
numerical change, and say so in that PR.

Not collected by pytest (no `test_` prefix); imported by the test.

Decks (all `LoadControl 1.0` over Path series, zero free equations, so the
Gauss-point strain is exactly the prescribed one):
  md3d       vanilla `ManzariDafalias` on the 3D confine-first cube -- the
             vanilla class itself was edited, so it is pinned directly.
  ls3d       `LadrunoSANISAND` defaults, same cube.
  ls_ps      `LadrunoSANISAND` plane-strain quad with the TIMs campaign
             options (IntScheme 1, TanType 0, -maxSubsteps 20000, -Pmin
             0.0101, -Presidual 0, nu 0.312885).
  ls_ps_big  the same with 4 huge deviatoric steps (5 % strain each) at low
             confinement -- the leg that exercises error-test REJECTIONS and
             the dT_min branches the counters classify.
  ls_ps_cyc  the campaign material through a strain reversal (loading
             reversal -> alpha_in reset) at low confinement.
"""
import json
import os
import sys

from _testbed import ops
import test_ladruno_sanisand as sani

_HERE = os.path.dirname(os.path.abspath(__file__))
BASELINE = os.path.join(_HERE, "wp127_sanisand_byteid_baseline.json")

_XY = sani._XY
_PARAMS = list(sani._PARAMS)
# The TIMs campaign set (_tims_2d_model_requests_2026-09-25/README.md): the
# test module's set with nu = 0.312885 (the Jaky K0 substitution).
_CAMPAIGN = list(_PARAMS)
_CAMPAIGN[1] = 0.312885
_CAMPAIGN_OPTS = (1, 0, 1, 1.0e-7, 1.0e-7,
                  "-Presidual", 0.0, "-Pmin", 0.0101, "-maxSubsteps", 20000,
                  "-flipAlphaIn", "init")


def _analysis():
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-13, 25, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")


def _series(lat_path, ax_path):
    lat = list(lat_path) + [lat_path[-1]]      # HOLD point: PathSeries -> 0 past its end
    ax = list(ax_path) + [ax_path[-1]]
    ops.timeSeries("Path", 1, "-dt", 1.0, "-values", *lat)
    ops.timeSeries("Path", 2, "-dt", 1.0, "-values", *ax)


def _paths(n_conf, e_conf, dev_increments):
    """Unit-magnitude shapes: confine 0->1, then per-step (d_lat, d_ax) in
    units of e_conf (sign: + = further compression)."""
    lat = [i / n_conf for i in range(n_conf + 1)]
    ax = list(lat)
    for d_lat, d_ax in dev_increments:
        lat.append(lat[-1] + d_lat / e_conf)
        ax.append(ax[-1] + d_ax / e_conf)
    return lat, ax


def _build_3d(matcmd, params, opts, n_conf, e_conf, incs):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for k in range(2):
        for j, (x, y) in enumerate(_XY):
            ops.node(4 * k + j + 1, x, y, float(k))
    ops.nDMaterial(matcmd, 1, *params, *opts)
    ops.element("stdBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
    for k in range(2):
        for j, (x, y) in enumerate(_XY):
            ops.fix(4 * k + j + 1, 1 if x == 0. else 0, 1 if y == 0. else 0,
                    1 if k == 0 else 0)
    _series(*_paths(n_conf, e_conf, incs))
    ops.pattern("Plain", 1, 1)
    for k in range(2):
        for j, (x, y) in enumerate(_XY):
            n = 4 * k + j + 1
            if x == 1.:
                ops.sp(n, 1, -e_conf)
            if y == 1.:
                ops.sp(n, 2, -e_conf)
    ops.pattern("Plain", 2, 2)
    for j in range(4):
        ops.sp(4 + j + 1, 3, -e_conf)
    _analysis()


def _build_ps(params, opts, n_conf, e_conf, incs):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for j, (x, y) in enumerate(_XY):
        ops.node(j + 1, x, y)
    ops.nDMaterial("LadrunoSANISAND", 1, *params, *opts)
    ops.element("quad", 1, 1, 2, 3, 4, 1.0, "PlaneStrain", 1)
    for j, (x, y) in enumerate(_XY):
        ops.fix(j + 1, 1 if x == 0. else 0, 1 if y == 0. else 0)
    _series(*_paths(n_conf, e_conf, incs))
    ops.pattern("Plain", 1, 1)
    for j, (x, y) in enumerate(_XY):
        if x == 1.:
            ops.sp(j + 1, 1, -e_conf)
    ops.pattern("Plain", 2, 2)
    for j, (x, y) in enumerate(_XY):
        if y == 1.:
            ops.sp(j + 1, 2, -e_conf)
    _analysis()


def _hex(v):
    return [float(x).hex() for x in v]


def _snapshot():
    out = []
    for gp in (1,):
        out.extend(_hex(ops.eleResponse(1, "material", gp, "stress")))
        out.extend(_hex(ops.eleResponse(1, "material", gp, "strain")))
        out.extend(_hex(ops.eleResponse(1, "material", gp, "state")))
    return out


def _run(n_conf, n_steps):
    rows = []
    ops.updateMaterialStage("-material", 1, "-stage", 0)
    for _ in range(n_conf):
        rc = ops.analyze(1)
        rows.append([rc] + _snapshot())
    ops.updateMaterialStage("-material", 1, "-stage", 1)
    for _ in range(n_steps):
        rc = ops.analyze(1)
        rows.append([rc] + _snapshot())
    return rows


def _iso_dev(n, e_ax, lat):
    return [(-lat * e_ax / n, e_ax / n)] * n


def deck_md3d():
    incs = _iso_dev(40, 5.0e-3, 0.5)
    _build_3d("ManzariDafalias", _PARAMS, (), 10, 3.0e-6, incs)
    return _run(10, len(incs))


def deck_ls3d():
    incs = _iso_dev(40, 5.0e-3, 0.5)
    _build_3d("LadrunoSANISAND", _PARAMS, (), 10, 3.0e-6, incs)
    return _run(10, len(incs))


def deck_ls_ps():
    incs = _iso_dev(40, 5.0e-3, 1.0)
    _build_ps(_CAMPAIGN, _CAMPAIGN_OPTS, 10, 3.0e-6, incs)
    return _run(10, len(incs))


def deck_ls_ps_big():
    incs = _iso_dev(4, 0.2, 1.0)
    _build_ps(_CAMPAIGN, _CAMPAIGN_OPTS, 10, 3.0e-6, incs)
    return _run(10, len(incs))


def deck_ls_ps_cyc():
    e = 2.0e-3
    incs = ([(-e / 10, e / 10)] * 10 + [(e / 10, -e / 10)] * 20
            + [(-e / 10, e / 10)] * 20)
    _build_ps(_CAMPAIGN, _CAMPAIGN_OPTS, 10, 3.0e-6, incs)
    return _run(10, len(incs))


DECKS = {
    "md3d": deck_md3d,
    "ls3d": deck_ls3d,
    "ls_ps": deck_ls_ps,
    "ls_ps_big": deck_ls_ps_big,
    "ls_ps_cyc": deck_ls_ps_cyc,
}


def run_all():
    return {name: fn() for name, fn in DECKS.items()}


if __name__ == "__main__":
    res = run_all()
    if "--write" in sys.argv:
        with open(BASELINE, "w") as fh:
            json.dump({"build": ops.ladrunoBuild() if hasattr(ops, "ladrunoBuild") else "?",
                       "decks": res}, fh, indent=0)
        print("wrote", BASELINE)
    for name, rows in res.items():
        print(name, len(rows), "rows; rc set", sorted({r[0] for r in rows}))
