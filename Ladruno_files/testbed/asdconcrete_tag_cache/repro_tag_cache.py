"""Minimal user-level reproducer: ASDConcrete3D silently reuses the FIRST
hardening law ever registered under a material tag.

No tests, no gmsh, no elements. Define tag 1 with a ductile tension backbone,
wipe, redefine tag 1 with a brittle one, and pull a single cube each time.
If the cache were sound the two runs would differ; they do not.

Control: doing the brittle run FIRST under a fresh tag proves the two backbones
really are different, so a match in the first experiment can only be the cache.
"""
import opensees as ops

E, NU, FT = 30000.0, 0.2, 3.0
ET0 = FT / E

DUCTILE = dict(Te=[ET0, 2.0e-3, 1.5e-2], Ts=[FT, 1.8, 0.15])
BRITTLE = dict(Te=[ET0, 1.0e-3, 8.0e-3], Ts=[FT, 1.0, 0.05])
CE, CS, CD = [-ET0, -1.0e-3], [-FT, -1.0], [0.0, 0.0]

CONN = [1, 2, 3, 4, 5, 6, 7, 8]
CUBE = {1: (0, 0, 0), 2: (1, 0, 0), 3: (1, 1, 0), 4: (0, 1, 0),
        5: (0, 0, 1), 6: (1, 0, 1), 7: (1, 1, 1), 8: (0, 1, 1)}


def td(law):
    return [0.0] + [1.0 - y / (E * x) for x, y in zip(law["Te"][1:], law["Ts"][1:])]


def pull(tag, law, nsteps=400, dmax=5.0e-2, side=2.0):
    """Uniaxial-stress pull of one cube; returns the DISSIPATED WORK.

    NOT the peak reaction: peak is FT*area = 3.0 for every backbone here, so it
    reads identical whatever the cache does and cannot discriminate. The
    softening branch is the thing that differs, so integrate under the curve.

    side=2.0 with -autoRegularization 1.0 is REQUIRED, not cosmetic: regularize()
    early-returns when lch == lch_ref, and deRegularize() -- the call that reads
    the tag-keyed global cache -- sits just past that guard. A unit cube with
    lch_ref=1.0 therefore never touches the cache and the bug stays invisible.
    """
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, c in CUBE.items():
        ops.node(t, side*float(c[0]), side*float(c[1]), side*float(c[2]))
    ops.nDMaterial("ASDConcrete3D", tag, E, NU, "-rho", 0.0,
                   "-Te", *law["Te"], "-Ts", *law["Ts"], "-Td", *td(law),
                   "-Ce", *CE, "-Cs", *CS, "-Cd", *CD,
                   "-autoRegularization", 1.0)
    ops.element("LadrunoBrick", 1, *CONN, tag)
    for n in (1, 4, 5, 8):
        ops.fix(n, 1, 0, 0)
    for n in (1, 2, 5, 6):
        ops.fix(n, 0, 1, 0)
    for n in (1, 2, 3, 4):
        ops.fix(n, 0, 0, 1)
    ops.equalDOF(2, 3, 1); ops.equalDOF(2, 6, 1); ops.equalDOF(2, 7, 1)
    ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1)
    ops.load(2, 1.0, 0.0, 0.0)
    ops.constraints("Transformation"); ops.numberer("Plain")
    ops.system("FullGeneral"); ops.test("NormDispIncr", 1e-10, 100, 0)
    ops.algorithm("Newton")
    ops.integrator("DisplacementControl", 2, 1, dmax / nsteps)
    ops.analysis("Static")
    W, uprev, Rprev, du = 0.0, 0.0, 0.0, dmax / nsteps
    for i in range(nsteps):
        if ops.analyze(1) != 0:
            break
        ops.reactions()
        R = abs(sum(ops.nodeReaction(n)[0] for n in (1, 4, 5, 8)))
        u = (i + 1) * du
        W += 0.5 * (R + Rprev) * (u - uprev)
        uprev, Rprev = u, R
    return W


print("EXPERIMENT  same tag 1, ductile first then brittle")
a = pull(1, DUCTILE)
b = pull(1, BRITTLE)          # <-- asks for brittle, silently gets ductile
print("  tag1 ductile W = %.6f" % a)
print("  tag1 brittle W = %.6f" % b)
print("  identical         = %s   <-- BUG if True" % (abs(a - b) < 1e-6 * max(a, 1e-30)))

print("CONTROL     the two backbones really are different (fresh tags)")
c = pull(21, DUCTILE)
d = pull(22, BRITTLE)
print("  tag21 ductile W = %.6f" % c)
print("  tag22 brittle W = %.6f" % d)
print("  differ             = %s   <-- must be True" % (abs(c - d) > 1e-6 * max(c, 1e-30)))
