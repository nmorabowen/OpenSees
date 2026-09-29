"""WP-130 merge: WHAT changed in WP-129's byte-identity decks, before re-baselining.

WP-129's pinned rows carry the 17-column WP-127 `substepStats`; after the
WP-130 merge the census has LMS_COUNT = 29 columns (12 CPPM columns appended),
and `-cppmTangent fixed` is the LadrunoSANISAND default, so the IntScheme 2 +
TanType 2 deck `ls3d_s2` hands out the CPPM tangent with the corrected sign.
This script runs WP-129's decks on the current build and classifies every
difference against the OLD baseline:
  * census   -- the snapshot now has 12 more trailing census entries; the
                first 17 must be unchanged
  * tangent  -- only in ls3d_s2: the 36 tangent entries must be the old ones
                NEGATED (to 5 %: the corrected low-p D_factor derivative also
                changes the tangent itself), and the stress / strain / state
                entries may move by at most 1e-5 relative (that correction
                moves the local Newton's iterates; this deck is at p ~ 2 kPa)
  * anything else is a FAILURE and the baseline must NOT be regenerated.

usage: python -S <bootstrap> wp129_rebaseline_check.py [--write]
(--write regenerates wp129_sanisand_byteid_baseline.json only if the check passes)
"""
import json
import sys

import wp129_sanisand_byteid as B

N_OLD = 17
N_NEW = 29
N_TAN_3D = 36
N_TAN_PS = 9


def main():
    with open(B.BASELINE) as fh:
        old = json.load(fh)
    got = B.run_all()
    bad = []
    tangent_rows = 0
    state_dev, flip_err = [], []
    for name, rows in old["decks"].items():
        new_rows = got[name]
        assert len(new_rows) == len(rows), name
        ladruno = not name.startswith("md3d")
        ntan = N_TAN_PS if name.startswith("ls_ps") else N_TAN_3D
        for k, (a, b) in enumerate(zip(new_rows, rows)):
            if k >= B.NONDETERMINISTIC.get(name, len(rows)):
                break   # WP-129's own rule: rows past this are not reproducible (ls3d_s4)
            if not ladruno:
                if a != b:
                    bad.append((name, k, "vanilla row moved"))
                continue
            # layout: [rc] + stress + strain + state + tangent + census
            if len(a) != len(b) + (N_NEW - N_OLD):
                bad.append((name, k, "unexpected width", len(a), len(b)))
                continue
            head_a, cen_a = a[:len(a) - N_NEW], a[len(a) - N_NEW:]
            head_b, cen_b = b[:len(b) - N_OLD], b[len(b) - N_OLD:]
            if cen_a[:N_OLD] != cen_b:
                bad.append((name, k, "WP-127 census columns moved"))
            if head_a == head_b:
                continue
            if name != "ls3d_s2":
                bad.append((name, k, "non-census entry moved"))
                continue
            # stress/strain/state: the `fixed` default also corrects the low-p
            # D_factor derivative in the LOCAL Jacobian (this deck sits at
            # p ~ 2 kPa < 0.05 P_atm), which moves the local Newton's iterates,
            # so the converged root moves within the local tolerance
            xa = [float.fromhex(x) for x in head_a[1:-ntan]]
            xb = [float.fromhex(x) for x in head_b[1:-ntan]]
            sc = max(1.0, max(abs(x) for x in xb))
            d = max(abs(x - y) for x, y in zip(xa, xb)) / sc
            state_dev.append(d)
            if head_a[0] != head_b[0] or d > 1e-5:
                bad.append((name, k, "state moved beyond the local tolerance", d))
                continue
            ta = [float.fromhex(x) for x in head_a[-ntan:]]
            tb = [float.fromhex(x) for x in head_b[-ntan:]]
            nb = sum(y * y for y in tb) ** 0.5
            e = sum((x + y) ** 2 for x, y in zip(ta, tb)) ** 0.5 / nb
            flip_err.append(e)
            if e <= 5e-2:
                tangent_rows += 1
            else:
                bad.append((name, k, "tangent change is not a sign flip", e))
    print("tangent rows sign-flipped (ls3d_s2):", tangent_rows,
          "| max state deviation (rel):", max(state_dev) if state_dev else 0.0,
          "| max ||T_new + T_old||/||T_old||:", max(flip_err) if flip_err else 0.0)
    print("unexplained differences:", len(bad))
    for x in bad[:20]:
        print("  ", x)
    if not bad and "--write" in sys.argv:
        with open(B.BASELINE, "w") as fh:
            json.dump({"build": got and (B.ops.ladrunoBuild() if hasattr(B.ops, "ladrunoBuild") else "?"),
                       "decks": got}, fh, indent=0)
        print("rewrote", B.BASELINE)
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
