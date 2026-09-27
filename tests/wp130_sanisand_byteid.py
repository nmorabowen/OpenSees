"""WP-130 byte-identity decks for IntScheme 2 (BackwardEuler_CPPM).

WP-130 edits the vanilla `ManzariDafalias` CPPM path (`// Ladruno WP-130`):
the three `Matrix::Invert` calls in `NewtonSol` go to a stack-local LU, the
dead `NewtonIter`'s function statics become locals, the halving cap becomes a
seam, and new opt-in branches (refuse on failure, local line search, the
ModifiedEuler -> CPPM fallback) sit behind flags that default OFF.  At the
defaults none of that may move a single computed bit.  WP-127's decks
(`wp127_sanisand_byteid.py`) are all IntScheme 1 and never enter the CPPM, so
this module adds IntScheme-2 decks and returns every committed stress /
strain / state / tangent as `float.hex` strings.

The pinned reference `wp130_sanisand_byteid_baseline.json` beside this file
was captured with the PRE-WP-130 binary (source tree 234a75751 = WP-127's tip
minus a tests-only commit, i.e. no SRC difference from e8fb51cdb) by

    python -S <bootstrap> wp130_sanisand_byteid.py --write

and `test_ladruno_sanisand_cppm_newton.py::test_scheme2_defaults_are_byte_identical`
compares the current build against it.  Regenerate it ONLY for a deliberate
numerical change, and say so in that PR.

Not collected by pytest (no `test_` prefix).

Decks:
  md3d_s2      vanilla `ManzariDafalias`, IntScheme 2, TanType 2, 3D cube,
               all faces prescribed (zero free DOF).
  ls3d_s2      `LadrunoSANISAND` IntScheme 2 on the same cube.
  ls3d_s2_big  the same through 4 huge deviatoric steps (20 % strain each):
               the CPPM's halving ladder and its explicit fallback fire.
  ls_ps_s2     the TIMs campaign set in plane strain, IntScheme 2, TanType 2.
  ls_ps_s2_cyc the same through a strain reversal.
  ls_ps_s2_free  a quad with LOADED edges (genuine free DOF, global Newton):
               the algorithmic tangent (the three inverted 6x6s) steers the
               iterations, so the per-step Newton iteration count is pinned
               too.
"""
import json
import os
import sys

from _testbed import ops
import test_ladruno_sanisand as sani
import wp127_sanisand_byteid as b127

_HERE = os.path.dirname(os.path.abspath(__file__))
BASELINE = os.path.join(_HERE, "wp130_sanisand_byteid_baseline.json")

_PARAMS = list(sani._PARAMS)
_S2 = (2, 2, 1, 1.0e-7, 1.0e-7)
_CAMPAIGN_S2 = (2, 2, 1, 1.0e-7, 1.0e-7,
                "-Presidual", 0.0, "-Pmin", 0.0101, "-maxSubsteps", 20000,
                "-flipAlphaIn", "init")


def _hex(v):
    return [float(x).hex() for x in v]


def _snapshot():
    out = []
    out.extend(_hex(ops.eleResponse(1, "material", 1, "stress")))
    out.extend(_hex(ops.eleResponse(1, "material", 1, "strain")))
    out.extend(_hex(ops.eleResponse(1, "material", 1, "state")))
    try:
        out.extend(_hex(ops.eleResponse(1, "material", 1, "tangent")))
    except Exception:  # pragma: no cover -- vanilla class has no "tangent"
        pass
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
        rows.append([rc, ops.testIter()] + _snapshot())
    return rows


def deck_md3d_s2():
    incs = b127._iso_dev(40, 5.0e-3, 0.5)
    b127._build_3d("ManzariDafalias", _PARAMS, _S2, 10, 3.0e-6, incs)
    return _run(10, len(incs))


def deck_ls3d_s2():
    incs = b127._iso_dev(40, 5.0e-3, 0.5)
    b127._build_3d("LadrunoSANISAND", _PARAMS, _S2, 10, 3.0e-6, incs)
    return _run(10, len(incs))


def deck_ls3d_s2_big():
    incs = b127._iso_dev(4, 0.2, 0.5)
    b127._build_3d("LadrunoSANISAND", _PARAMS, _S2, 10, 3.0e-6, incs)
    return _run(10, len(incs))


def deck_ls_ps_s2():
    incs = b127._iso_dev(40, 5.0e-3, 1.0)
    b127._build_ps(b127._CAMPAIGN, _CAMPAIGN_S2, 10, 3.0e-6, incs)
    return _run(10, len(incs))


def deck_ls_ps_s2_cyc():
    e = 2.0e-3
    incs = ([(-e / 10, e / 10)] * 10 + [(e / 10, -e / 10)] * 20
            + [(-e / 10, e / 10)] * 20)
    b127._build_ps(b127._CAMPAIGN, _CAMPAIGN_S2, 10, 3.0e-6, incs)
    return _run(10, len(incs))


def build_free_quad(opts, lateral=5.0):
    """Plane-strain quad, LOADED right and top edges (free DOF), stage 0
    confinement to `2*lateral` kPa, then the stage flip.  Returns with the
    analysis set up for a deviatoric load push (pattern 2 not yet added)."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for j, (x, y) in enumerate(sani._XY):
        ops.node(j + 1, x, y)
    ops.nDMaterial("LadrunoSANISAND", 1, *b127._CAMPAIGN, *opts)
    ops.element("quad", 1, 1, 2, 3, 4, 1.0, "PlaneStrain", 1)
    for j, (x, y) in enumerate(sani._XY):
        ops.fix(j + 1, 1 if x == 0. else 0, 1 if y == 0. else 0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for j, (x, y) in enumerate(sani._XY):
        ops.load(j + 1, -lateral if x == 1. else 0.0, -lateral if y == 1. else 0.0)
    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-10, 30, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 0.1)
    ops.analysis("Static")
    ops.updateMaterialStage("-material", 1, "-stage", 0)
    for _ in range(10):
        assert ops.analyze(1) == 0
    ops.updateMaterialStage("-material", 1, "-stage", 1)
    ops.loadConst("-time", 0.0)


def deck_ls_ps_s2_free():
    build_free_quad(_CAMPAIGN_S2)
    ops.timeSeries("Linear", 2)
    ops.pattern("Plain", 2, 2)
    for j, (x, y) in enumerate(sani._XY):
        if y == 1.:
            ops.load(j + 1, 0.0, -1.0)
    ops.integrator("LoadControl", 0.5)
    rows = []
    for _ in range(12):
        rc = ops.analyze(1)
        rows.append([rc, ops.testIter()] + _snapshot())
        if rc != 0:
            break
    return rows


DECKS = {
    "md3d_s2": deck_md3d_s2,
    "ls3d_s2": deck_ls3d_s2,
    "ls3d_s2_big": deck_ls3d_s2_big,
    "ls_ps_s2": deck_ls_ps_s2,
    "ls_ps_s2_cyc": deck_ls_ps_s2_cyc,
    "ls_ps_s2_free": deck_ls_ps_s2_free,
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
        iters = [r[1] for r in rows if len(r) > 1 and isinstance(r[1], int)]
        print(name, len(rows), "rows; rc set", sorted({r[0] for r in rows}),
              "width", len(rows[-1]))
