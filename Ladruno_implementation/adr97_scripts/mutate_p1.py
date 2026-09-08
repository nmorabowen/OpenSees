"""ADR-97 P1 mutation gate (ADR-87 D2 / ADR-97 gate 5).

Drops ONE term of the consistent tangent -- the `dl * dm/dsigma` term of the
algorithmic elastic modulus `Xi = (E^-1 + dl*dm/ds)^-1` -- from the closest-point
Jacobian's J_ss block, on a SCRATCH build.  With it gone, `Xi` collapses back to
`E` and `C_alg` degenerates into exactly the shipped `Continuum` operator, whose
error against a central difference of the material's own committed response
ADR-94 M3 measured at 57 %.  Gate 2 must therefore go RED.

Note this mutation changes ONLY the Jacobian, never the residual: the converged
stress stays exact (inexact Newton), so gate 1 and gate 4 should stay GREEN and
only the tangent gates should die.  That is the point -- it isolates the tangent.

usage:  python3.12 mutate_p1.py {apply|revert} [<worktree root>]
"""
import io
import os
import sys

MODE = sys.argv[1]
ROOT = sys.argv[2] if len(sys.argv) > 2 else \
    r"C:\Users\nmb\Documents\Github\OpenSees\.claude\worktrees\asdplastic-review-plan-62585c"
PATH = os.path.join(ROOT, "SRC", "material", "nD", "ASDPlasticMaterial3D",
                    "ASDPlasticMaterial3D.h")

ORIG = """                Jm(i, j) = ((i == j) ? 1.0 : 0.0) + dl * Edmds(i, j);"""
MUT = """                Jm(i, j) = ((i == j) ? 1.0 : 0.0); // ADR-97 MUTATION: dl*dm/ds dropped"""

src = io.open(PATH, encoding="utf-8", errors="surrogateescape", newline="").read()
crlf = "\r\n" in src
o = ORIG.replace("\n", "\r\n") if crlf else ORIG
m = MUT.replace("\n", "\r\n") if crlf else MUT

if MODE == "apply":
    assert src.count(o) == 1, "anchor not found/unique"
    src = src.replace(o, m)
elif MODE == "revert":
    assert src.count(m) == 1, "mutation not present"
    src = src.replace(m, o)
else:
    raise SystemExit("usage: mutate_p1.py {apply|revert}")

io.open(PATH, "w", encoding="utf-8", errors="surrogateescape",
        newline="").write(src)
print("mutate_p1: %s ok" % MODE)
