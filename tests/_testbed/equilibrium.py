"""assert_equilibrium — the universal guard (G7).

The apeGmsh bridge does NOT auto-emit loads/masses/SP (but DOES auto-emit MP
constraints), so a mis-wired model can converge with ZERO applied load and
"pass" any assertion that happens to match ~0. One global force balance per
case kills that entire class:

    Sum(reactions over all nodes, per global direction) + Sum(applied) == 0

After a static solve and ops.reactions(), reactions are nonzero only at
constrained DOFs and sum to -Sum(applied). So the residual must vanish.
"""
from ._ops import ops


def assert_equilibrium(node_tags, n_force_dofs, applied_total=None, atol=1e-6):
    """Assert global static FORCE equilibrium: sum of reactions + sum of applied
    loads == 0 over the translational DOFs.

    Only the leading `n_force_dofs` (translational) DOFs are checked. Moment DOFs
    are deliberately NOT summed: a force applied at a lever arm contributes to the
    moment balance through r x F, so a naive sum of reaction moments does not
    vanish (e.g. a tip-loaded cantilever has base moment -P*L with zero applied
    moment). Force balance is what actually catches the zero-load false-green
    class (G7); full moment balance needs geometry and is out of scope here.

    Because ndf alone is ambiguous (ndf=3 is a 2D beam = 2 forces + 1 moment OR a
    3D solid = 3 forces), the caller must state how many leading DOFs are forces:
      2D solid -> 2 | 3D solid -> 3 | 2D beam -> 2 | 3D beam -> 3

    node_tags     : iterable of every node tag in the model
    n_force_dofs  : number of leading translational DOFs to balance
    applied_total : length-n_force_dofs total applied force per global direction
                    (default zeros — e.g. a self-equilibrated or gravity case
                    where you pass the gravity resultant explicitly)
    """
    if applied_total is None:
        applied_total = [0.0] * n_force_dofs
    ops.reactions()
    react = [0.0] * n_force_dofs
    for n in node_tags:
        r = ops.nodeReaction(int(n))
        for d in range(min(n_force_dofs, len(r))):
            react[d] += r[d]
    residual = [react[d] + applied_total[d] for d in range(n_force_dofs)]
    bad = [(d, residual[d]) for d in range(n_force_dofs) if abs(residual[d]) > atol]
    assert not bad, (
        f"force equilibrium violated (reactions+applied != 0): {bad}; "
        f"reactions={react}, applied={list(applied_total)}. "
        "Most likely the load/mass/SP was never emitted onto the bridge."
    )
    return residual
