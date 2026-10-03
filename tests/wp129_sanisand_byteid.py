"""WP-129 byte-identity decks: every EXISTING SANISAND IntScheme, before/after SAS-ME.

WP-129 adds a NEW integration scheme (IntScheme 129, "SAS-ME") to
`ManzariDafalias` (vanilla, `// Ladruno WP-129`) and `LadrunoSANISAND`.  The
existing schemes must not move a single computed bit.  This module runs one
zero-free-DOF material-point deck per existing scheme (0-9, 45 on
`LadrunoSANISAND`; 1, 2, 45 on vanilla `ManzariDafalias`; TanType 0/1/2 on the
IntScheme-1 plane-strain campaign material) and returns every committed
stress / strain / state vector, the delivered tangent and the WP-127
`substepStats` census as `float.hex` strings, so two builds can be compared
EXACTLY (not to a tolerance) -- the WP-127 method (wp127_sanisand_byteid.py).

The pinned reference `wp129_sanisand_byteid_baseline.json` beside this file was
captured with the UNMODIFIED WP-127 binary (ladrunoBuild 234a7575; `git diff
234a7575 e8fb51cdb -- SRC` is empty, so it is the binary of the WP-127 tree this
branch was cut from) by

    python -S <bootstrap> wp129_sanisand_byteid.py --write

and `test_ladruno_sanisand_sasme.py::test_existing_schemes_byte_identical`
compares the current build against it.  Regenerate it ONLY for a deliberate
numerical change to an existing scheme, and say so in that PR.

Re-pinned since (one deck at a time, every other deck left byte-for-byte):
  ls3d_s5  WP-158 -- ForwardEuler's shadowed `r` and its two tangent defects
           fixed; plastic rows 10-29 move, elastic rows 0-9 do not. Captured
           with the WP-158 build (origin/ladruno d63f49750 + the WP-158
           ManzariDafalias.cpp edit), 2026-10-01; the 17 other decks matched
           the pinned values exactly on that build.

Not collected by pytest (no `test_` prefix); imported by the test.
"""
import json
import os
import sys

from _testbed import ops
import test_ladruno_sanisand as sani
import wp127_sanisand_byteid as w127

_HERE = os.path.dirname(os.path.abspath(__file__))
BASELINE = os.path.join(_HERE, "wp129_sanisand_byteid_baseline.json")

_PARAMS = list(sani._PARAMS)
_CAMPAIGN = list(w127._CAMPAIGN)


def _hex(v):
    return [float(x).hex() for x in v]


def _snapshot(ladruno=True):
    out = []
    for gp in (1,):
        out.extend(_hex(ops.eleResponse(1, "material", gp, "stress")))
        out.extend(_hex(ops.eleResponse(1, "material", gp, "strain")))
        out.extend(_hex(ops.eleResponse(1, "material", gp, "state")))
        out.extend(_hex(ops.eleResponse(1, "material", gp, "tangent")))
        if ladruno:
            out.extend(_hex(ops.eleResponse(1, "material", gp, "substepStats")))
    return out


def _run(n_conf, n_steps, ladruno=True):
    rows = []
    ops.updateMaterialStage("-material", 1, "-stage", 0)
    for _ in range(n_conf):
        rc = ops.analyze(1)
        rows.append([rc] + _snapshot(ladruno))
    ops.updateMaterialStage("-material", 1, "-stage", 1)
    for _ in range(n_steps):
        rc = ops.analyze(1)
        rows.append([rc] + _snapshot(ladruno))
    return rows


def _deck_3d(matcmd, scheme, tantype, incs):
    opts = (scheme, tantype, 1, 1.0e-7, 1.0e-7)
    w127._build_3d(matcmd, _PARAMS, opts, 10, 3.0e-6, incs)
    return _run(10, len(incs), ladruno=(matcmd == "LadrunoSANISAND"))


def _deck_ps(tantype, incs):
    opts = (1, tantype) + tuple(w127._CAMPAIGN_OPTS[2:])
    w127._build_ps(_CAMPAIGN, opts, 10, 3.0e-6, incs)
    return _run(10, len(incs))


# moderate triaxial-type deviatoric path: plastic, error control binding
_INCS_3D = w127._iso_dev(20, 2.0e-3, 0.5)
# plane-strain campaign cycle (loading reversal -> alpha_in reset)
_E = 2.0e-3
_INCS_CYC = ([(-_E / 10, _E / 10)] * 8 + [(_E / 10, -_E / 10)] * 12
             + [(-_E / 10, _E / 10)] * 8)

LS_SCHEMES = (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 45)
# decks whose PLASTIC rows are not reproducible even on the unmodified binary
# (measured: two runs of ls3d_s4 in one process differ) -> number of leading
# rows (the elastic stage) that ARE pinned. Cause: MaxEnergyInc's uninitialised
# `double nG, nK` handed to ForwardEuler when it sub-steps (LEDGER_quirks).
NONDETERMINISTIC = {"ls3d_s4": 10}
MD_SCHEMES = (1, 2, 45)


def decks():
    d = {}
    for s in LS_SCHEMES:
        d[f"ls3d_s{s}"] = (lambda s=s: _deck_3d("LadrunoSANISAND", s, 2, _INCS_3D))
    d["ls3d_s1_tan1"] = lambda: _deck_3d("LadrunoSANISAND", 1, 1, _INCS_3D)
    for s in MD_SCHEMES:
        d[f"md3d_s{s}"] = (lambda s=s: _deck_3d("ManzariDafalias", s, 2, _INCS_3D))
    for t in (0, 1, 2):
        d[f"ls_ps_cyc_tan{t}"] = (lambda t=t: _deck_ps(t, _INCS_CYC))
    return d


def run_all():
    return {name: fn() for name, fn in decks().items()}


if __name__ == "__main__":
    res = run_all()
    if "--write" in sys.argv:
        with open(BASELINE, "w") as fh:
            json.dump({"build": ops.ladrunoBuild() if hasattr(ops, "ladrunoBuild") else "?",
                       "decks": res}, fh, indent=0)
        print("wrote", BASELINE)
    for name, rows in res.items():
        print(name, len(rows), "rows; rc set", sorted({r[0] for r in rows}))
