"""P0 model builder, shared by the serial reference and the 2-rank run.

Two stacked unit bricks with NTS penalty contact between them.
  bottom block = MASTER (top face -> master segment), base fixed
  top block    = SLAVE  (bottom nodes -> slave node set), loaded down
Lateral DOFs fixed everywhere => a clean 1-D column; only w (z) is free.

rank=None  -> build the whole model in one domain (serial reference)
rank=0     -> master-owning rank: bottom block + contact + GHOSTS of 11..14
rank=1     -> slave rank: top block + the load
"""

BOT_NODES = {
    1: (0.0, 0.0, 0.0), 2: (1.0, 0.0, 0.0), 3: (1.0, 1.0, 0.0), 4: (0.0, 1.0, 0.0),
    5: (0.0, 0.0, 1.0), 6: (1.0, 0.0, 1.0), 7: (1.0, 1.0, 1.0), 8: (0.0, 1.0, 1.0),
}
# The top block starts slightly INTERPENETRATED (gap0 = -1e-3). With an exactly
# zero gap the penalty is inactive at iteration 1, so the top block has a free
# rigid-body mode in z and Newton oscillates (measured: no convergence in 50
# iterations, residual stuck at 1.7e4). Starting inside contact makes the
# tangent non-singular from the first iteration. This is a property of the test
# model, not of the parallel question under study.
GAP0 = 1.0e-3
TOP_NODES = {
    11: (0.0, 0.0, 1.0 - GAP0), 12: (1.0, 0.0, 1.0 - GAP0),
    13: (1.0, 1.0, 1.0 - GAP0), 14: (0.0, 1.0, 1.0 - GAP0),
    15: (0.0, 0.0, 2.0 - GAP0), 16: (1.0, 0.0, 2.0 - GAP0),
    17: (1.0, 1.0, 2.0 - GAP0), 18: (0.0, 1.0, 2.0 - GAP0),
}
BASE = (1, 2, 3, 4)              # fully fixed
INTERFACE_SLAVE = (11, 12, 13, 14)   # ghosted onto rank 0
LOADED = (15, 16, 17, 18)
E, NU = 2.0e7, 0.0
P_NODE = -2500.0


def build(ops, rank=None, mutate=None):
    """mutate=None      -> the design under test
       mutate='noghost' -> rank 0 does NOT declare the slave interface nodes
                           (tests ADR-78 D5: does a missing node degrade silently?)
       mutate='double'  -> the contact verb is emitted TWICE (tests INV-1:
                           does duplicate assembly double-count the force?)"""
    return _build(ops, rank, mutate)


def _build(ops, rank=None, mutate=None):
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 3)
    ops.nDMaterial('ElasticIsotropic', 1, E, NU)
    ops.timeSeries('Linear', 1)

    def put(tags, table):
        for t in tags:
            ops.node(t, *table[t])
            ops.fix(t, 1, 1, 1 if t in BASE else 0)

    if rank is None or rank == 0:
        put(BOT_NODES, BOT_NODES)
        ops.element('stdBrick', 101, 1, 2, 3, 4, 5, 6, 7, 8, 1)

    if rank == 0 and mutate != 'noghost':
        # GHOSTS: the slave interface nodes live on rank 1. Declare them here
        # with the owner's SP state replayed (apeGmsh ADR 0027 INV-2).
        for t in INTERFACE_SLAVE:
            ops.node(t, *TOP_NODES[t])
            ops.fix(t, 1, 1, 0)

    if rank is None or rank == 0:
        # contact is emitted on the MASTER-owning rank only (ADR-78 D1/INV-1)
        ops.contactSurface(1, '-master', 4, 5, 6, 7, 8)
        ops.contactSurface(2, '-slave', *INTERFACE_SLAVE)
        ops.contact(1, 1, 2, 'auto', '-outward', 0.0, 0.0, 1.0)
        if mutate == 'double':
            ops.contactSurface(3, '-master', 4, 5, 6, 7, 8)
            ops.contactSurface(4, '-slave', *INTERFACE_SLAVE)
            ops.contact(2, 3, 4, 'auto', '-outward', 0.0, 0.0, 1.0)

    if rank is None or rank == 1:
        put(TOP_NODES, TOP_NODES)
        ops.element('stdBrick', 102, 11, 12, 13, 14, 15, 16, 17, 18, 1)

    ops.pattern('Plain', 1, 1)
    if rank is None or rank == 1:
        for t in LOADED:
            ops.load(t, 0.0, 0.0, P_NODE)


def analyse(ops, numberer, system):
    ops.constraints('LadrunoContact')
    ops.numberer(numberer)
    ops.system(system)
    ops.test('NormDispIncr', 1.0e-10, 50, 0)
    ops.algorithm('Newton')
    ops.integrator('LoadControl', 1.0)
    ops.analysis('Static')
    return ops.analyze(1)
