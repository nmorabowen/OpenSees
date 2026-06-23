"""analyze_augmented — held-load Uzawa augmentation loop for mortar ALM contact (ADR-41 C2.2).

Penalty enforcement alone leaves a residual interface penetration g~ = O(1/epsN): the
contact is only as accurate as epsN is large, and a large epsN ill-conditions the system.
Augmented-Lagrange (Uzawa) drives g~ -> 0 at a FINITE epsN by accumulating a per-slave-node
normal multiplier lambda_N, updated once per Domain::commit():

    lambda_I  <-  min(0, lambda_I + epsN * gbar_I)        (gbar_I = g~_I^global / a_I^global)

Across load steps this augments for free on stock Newton. To converge the penetration at a
HELD load (the headline accuracy gate), this proc re-commits at a ZERO load increment: each
analyze(1) re-solves the same external load with the latest lambda and commits another Uzawa
step. The SOE equation count is CONSTANT across augmentations (lambda_N is Domain-side state,
NOT a new DOF; the active set is frozen within the sweep), so no domainChanged()/re-number.

This is the ADR-41 Q-DRIVER resolution: NO custom EquiSolnAlgo — a tested driver recipe.

Recipe (OpenSeesPy)::

    ops.constraints('LadrunoContact'); ops.numberer('Plain'); ops.system('FullGeneral')
    ops.test('NormDispIncr', 1e-11, 40); ops.algorithm('Newton'); ops.analysis('Static')
    ops.integrator('LoadControl', 1.0)
    ops.analyze(1)                                  # apply load -> penalty equilibrium (1 aug)
    from analyze_augmented import analyze_augmented
    status, n_aug, pen = analyze_augmented(ops, maxAug=20, augTol=1e-9)

Returns ``(status, n_aug, penetration)``:
    status 0  = converged (penetration < augTol),
    status 1  = hit maxAug without converging,
    status <0 = a re-solve returned that nonzero analyze() code.
``penetration`` is the engine's max KKT-active ||gbar||_inf (ops.ladrunoMortarPenetration()).
"""


def analyze_augmented(ops, maxAug=20, augTol=1.0e-9, verbose=False):
    """Drive mortar contact penetration to ``augTol`` by held-load Uzawa augmentation.

    Must be called AFTER the load has been applied and a first equilibrium reached (so the
    penalty pressure / active set are established). Switches the integrator to a zero-increment
    LoadControl and re-solves+commits until ``ops.ladrunoMortarPenetration() < augTol`` or
    ``maxAug`` augmentations have run.
    """
    # hold the external load: a zero-increment LoadControl re-solves at the SAME load level,
    # so the only thing that changes between augmentations is the committed multiplier lambda.
    ops.integrator("LoadControl", 0.0)
    pen = ops.ladrunoMortarPenetration()
    for k in range(1, int(maxAug) + 1):
        rc = ops.analyze(1)
        if rc != 0:
            return rc, k, pen
        pen = ops.ladrunoMortarPenetration()
        if verbose:
            print(f"  [augment {k}] penetration = {pen:.3e}")
        if pen < augTol:
            return 0, k, pen
    return 1, int(maxAug), pen
