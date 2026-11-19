"""WP-239 test model: LadrunoMembrane floors on frames whose beams are
re-created on hinge-twin nodes AFTER the membranes exist.

The shape of the YOO RevB NLRHA hinge step, at toy size: a 2 x 1 bay floor at
z = H on four columns per side, two membranes (100, 101) on seven beams, and a
"replay" that removes each beam and re-adds it with the SAME tag on twin nodes
coincident with the original ends. The twins are tied to the originals by
`equalDOF` (all six DOFs) or by a `zeroLength` (rigid springs, or with the
vertical-bending DOF released -- a real moment hinge).

Used in-process by tests/test_wp239_membrane_frame_replace.py and as the body
of the save/restore child processes, so it must import only `ops`.
"""
from ._ops import ops

LX, LY, H = 6.0, 4.0, 3.0
Q = 7.5                    # slab pressure magnitude; the membrane normal is +z, so -pressure -Q pushes down
E, G = 30.0e6, 12.5e6
A, J, IY, IZ = 0.18, 0.002, 0.003, 0.004
K_RIGID = 1.0e12

# slab nodes 1..6 at z = H, column bases 11..16 at z = 0
SLAB = {1: (0.0, 0.0), 2: (LX, 0.0), 3: (2 * LX, 0.0),
        4: (0.0, LY), 5: (LX, LY), 6: (2 * LX, LY)}
BEAMS = {1: (1, 2), 2: (2, 3), 3: (4, 5), 4: (5, 6),   # along x
         5: (1, 4), 6: (2, 5), 7: (3, 6)}             # along y
COLUMNS = {20 + n: (10 + n, n) for n in SLAB}
MEMBRANES = {100: ((1, 2, 5, 4), (1, 6, 3, 5)),
             101: ((2, 3, 6, 5), (2, 7, 4, 6))}
DECOY_NODES = (901, 902)


def _beam(tag, ni, nj):
    # all beams are horizontal: transf 1 has vecxz = global z
    ops.element("elasticBeamColumn", tag, ni, nj, A, E, G, J, IY, IZ, 1)


def _membranes():
    for m, (nodes, frames) in MEMBRANES.items():
        ops.element("LadrunoMembrane", m, "-nodes", *nodes, "-frames", *frames,
                    "-method", "tributary", "-pressure", -Q,
                    "-stiff", 25.0e6, 0.2, 0.15)


def build(hinge=None, decoy=False, membranes_after_hinge=False):
    """Build the model. hinge in {None, 'equalDOF', 'zeroLength',
    'zeroLengthPinned'}; decoy=True allocates another beam between each
    remove and re-add (an allocator that would reuse the freed block for it
    sends a stale-pointer slab share to the decoy). membranes_after_hinge
    creates the membranes after the replay, so they resolve against frames
    already on the twins -- what a database restore does."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    for n, (x, y) in SLAB.items():
        ops.node(n, x, y, H)
        ops.node(10 + n, x, y, 0.0)
        ops.fix(10 + n, 1, 1, 1, 1, 1, 1)
    ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
    ops.geomTransf("Linear", 2, 1.0, 0.0, 0.0)
    for c, (nb, nt) in COLUMNS.items():
        ops.element("elasticBeamColumn", c, nb, nt, A, E, G, J, IY, IZ, 2)
    for b, (ni, nj) in BEAMS.items():
        _beam(b, ni, nj)
    if not membranes_after_hinge:
        _membranes()

    if decoy:
        for k, n in enumerate(DECOY_NODES):
            ops.node(n, 50.0 + 3.0 * k, 0.0, 0.0)
            ops.fix(n, 1, 1, 1, 1, 1, 1)

    if hinge is not None:
        ops.uniaxialMaterial("Elastic", 1, K_RIGID)
        ops.uniaxialMaterial("Elastic", 2, 1.0e-6)
        for b, (ni, nj) in BEAMS.items():
            twins = []
            for end, n in enumerate((ni, nj)):
                t = 1000 * b + end + 1
                x, y = SLAB[n]
                ops.node(t, x, y, H)
                if hinge == "equalDOF":
                    ops.equalDOF(n, t, 1, 2, 3, 4, 5, 6)
                else:
                    mats = [1] * 6
                    if hinge == "zeroLengthPinned":
                        # release vertical bending: about global y for a beam
                        # along x, about global x for a beam along y
                        along_x = abs(SLAB[nj][0] - SLAB[ni][0]) > 0.0
                        mats[4 if along_x else 3] = 2
                    ops.element("zeroLength", 5000 + 10 * b + end, n, t,
                                "-mat", *mats, "-dir", 1, 2, 3, 4, 5, 6)
                twins.append(t)
            ops.remove("element", b)
            if decoy:
                ops.element("elasticBeamColumn", 900 + b, *DECOY_NODES,
                            A, E, G, J, IY, IZ, 1)
            _beam(b, *twins)
    if membranes_after_hinge:
        _membranes()


def apply_loads(lateral=True):
    """Pattern 1: slab pressure on both membranes. Patterns 2-4: many
    -beamUniform element loads (one eleLoad per beam). Pattern 5: lateral."""
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.eleLoad("-ele", *MEMBRANES, "-type", "-surfaceLoad", 1.0)
    for p, w in ((2, -1.25), (3, -0.5), (4, -2.0)):
        ops.pattern("Plain", p, 1)
        for b in BEAMS:
            ops.eleLoad("-ele", b, "-type", "-beamUniform", 0.0, w)
    if lateral:
        ops.pattern("Plain", 5, 1)
        ops.load(3, 40.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        ops.load(6, 40.0, 15.0, 0.0, 0.0, 0.0, 0.0)


def analysis():
    ops.system("BandGeneral")
    ops.numberer("RCM")
    ops.constraints("Transformation")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Linear")
    ops.analysis("Static")


def total_applied_vertical():
    """Downward total: slab pressure on both bays + the -beamUniform loads."""
    slab = Q * (2 * LX * LY)
    beam_len = sum(LX if abs(SLAB[j][0] - SLAB[i][0]) > 0 else LY
                   for i, j in BEAMS.values())
    return slab + (1.25 + 0.5 + 2.0) * beam_len


def snapshot():
    """Everything the tests compare: slab displacements, beam local forces,
    support reactions (after ops.reactions())."""
    ops.reactions()
    disp = {n: list(ops.nodeDisp(n)) for n in SLAB}
    forces = {b: list(ops.eleResponse(b, "localForce")) for b in BEAMS}
    reac = {n: list(ops.nodeReaction(n)) for n in [10 + s for s in SLAB] + [*DECOY_NODES]
            if n in ops.getNodeTags()}
    return {"disp": disp, "forces": forces, "reac": reac}
