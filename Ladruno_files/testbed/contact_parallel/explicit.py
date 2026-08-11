"""P0 explicit lane — the half that actually tests ADR-78 Q-EXPLICIT.

Implicit (Mumps) assembles globally, so a ghost DOF is unremarkable there. The
explicit lane is the interesting one: the integrator inverts a DIAGONAL mass and
a ghost node has ZERO mass on the owner rank. The claim under test is that
DistributedDiagonalSolver::solve() sums BOTH A (mass) and B (residual) over
shared DOFs, so the ghost's contact force reaches the native rank's mass.

Same two-brick model, now transient under a steady load.
"""
NSTEP, DT = 200, 5.0e-5
M_NODE = 1.0


def add_mass(ops, tags):
    for t in tags:
        ops.mass(t, M_NODE, M_NODE, M_NODE)


def analyse_explicit(ops, numberer, system):
    ops.constraints('LadrunoContact')
    ops.numberer(numberer)
    ops.system(system)
    ops.integrator('CentralDifferenceLadruno')
    ops.algorithm('Linear')
    ops.analysis('Transient')
    return ops.analyze(NSTEP, DT)
