"""ADR-97 P2 mutation gate.

Mutation: in the 6D back-transform of the principal-space Koiter tangent,
DROP the eigenprojection ROTATION term -- the shear-slot coefficients

    T[3+s, 3+s] = (y_i - y_j) / (x_i - x_j)

that carry the fact that the return moves the principal VALUES while the
principal DIRECTIONS follow the trial stress.  They are replaced by 1.0, which
is the value they take for an elastic step (y == x): a plausible, subtle bug --
the committed stress is untouched, `Backward_Euler` is untouched, and every
region still returns the exact closest point.

Only the NON-degenerate branch is mutated.  The l'Hopital limit used when two
trial eigenvalues coincide is left alone, so the mutation is expected to kill
the FACE tests of gate 2 (three separated principal stresses, rotation term
live) and to LEAVE the degenerate edge test alive -- which is what isolates the
term rather than merely proving that the tangent matters at all.

    python3.12 Ladruno_implementation/adr97_scripts/mutate_p2.py apply
    python3.12 Ladruno_implementation/adr97_scripts/mutate_p2.py revert
"""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
TARGET = os.path.join(ROOT, "SRC", "material", "nD", "ASDPlasticMaterial3D",
                      "ASDPlasticMaterial3D.h")

ORIG = """                Tp(3 + s, 3 + s) = (adx > eps_deg)
                                   ? (yv(i) - yv(j)) / dx
                                   : (dydx(i, i) - dydx(i, j));"""

MUT = """                Tp(3 + s, 3 + s) = (adx > eps_deg)
                                   ? 1.0   /* ADR-97 P2 MUTATION */
                                   : (dydx(i, i) - dydx(i, j));"""


def main():
    if len(sys.argv) != 2 or sys.argv[1] not in ("apply", "revert"):
        print(__doc__)
        return 1
    with open(TARGET, "r", encoding="utf-8", newline="") as fh:
        txt = fh.read()
    crlf = "\r\n" in txt
    a = ORIG.replace("\n", "\r\n") if crlf else ORIG
    b = MUT.replace("\n", "\r\n") if crlf else MUT
    frm, to = (a, b) if sys.argv[1] == "apply" else (b, a)
    if to in txt:
        print("already %sed" % sys.argv[1])
        return 0
    if txt.count(frm) != 1:
        print("ANCHOR MISS (%d matches)" % txt.count(frm))
        return 3
    with open(TARGET, "w", encoding="utf-8", newline="") as fh:
        fh.write(txt.replace(frm, to, 1))
    print("%sed the eigenprojection rotation-term mutation" % sys.argv[1])
    return 0


if __name__ == "__main__":
    sys.exit(main())
