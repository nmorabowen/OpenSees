"""analyze_augmented — held-load within-step Uzawa augmentation for mortar ALM contact (ADR-41 D1).

Penalty enforcement alone leaves a residual interface penetration g~ = O(1/epsN): the contact is
only as accurate as epsN is large, and a large epsN ill-conditions the system. Augmented-Lagrange
(Uzawa) drives g~ -> 0 at a FINITE epsN by accumulating a per-slave-node normal multiplier lambda_N,
updated once per Domain::commit():

    lambda_I  <-  min(0, lambda_I + epsN * gbar_I)        (gbar_I = g~_I^global / a_I^global)

Across load steps this augments for free on stock Newton (ONE augmentation per physical step). To
converge the penetration at a HELD load WITHIN a single step (the headline accuracy gate), this proc
re-commits at a ZERO load increment: each analyze(1) re-solves the SAME external load with the latest
lambda and commits another Uzawa step, so multiple augmentations fire inside one physical step. The
SOE equation count is CONSTANT across augmentations (lambda_N is Domain-side state, NOT a new DOF;
the active set is frozen within the sweep), so no domainChanged()/re-number.

D1 — NO RECORDER / LOAD / TIME CORRUPTION (the capstone contract #3 "without recorder/load
corruption" clause). The held-load re-commits route through Domain::commit(), which normally fires
recorders + advances committedTime/commitTag. That would inject SPURIOUS duplicate recorder samples
at the same pseudo-time and jump commitTag. So this proc brackets the sweep with
``ops.ladrunoBeginAugment()`` / ``ops.ladrunoEndAugment()``: between them Domain::commit() performs
the contact Uzawa lambda update but fires NO recorders and bumps NO commitTag. Therefore a model with
a recorder produces IDENTICAL recorder output (same sample count + same committedTime sequence) with
vs without within-step augmentation. The bracket is restored in a ``finally`` so an exception cannot
leave recorders muted. (committedTime is independently a no-op: a zero-increment LoadControl holds
the load factor constant, so the held re-solves also do not double-count the load — gate "no load
corruption".)

This is the ADR-41 Q-DRIVER resolution: NO custom EquiSolnAlgo — a tested driver recipe.

Recipe (OpenSeesPy)::

    ops.constraints('LadrunoContact'); ops.numberer('Plain'); ops.system('FullGeneral')
    ops.test('NormDispIncr', 1e-11, 40); ops.algorithm('Newton'); ops.analysis('Static')
    ops.integrator('LoadControl', 1.0)
    ops.analyze(1)                                  # apply load -> penalty equilibrium (1 aug)
    from analyze_augmented import analyze_augmented
    status, n_aug, pen = analyze_augmented(ops, maxAug=20, augTol=1e-9)
    # -> the recorder (if any) has exactly ONE sample for this step (the penalty equilibrium);
    #    the AUGMENTED state is available via ops.nodeDisp(...) and propagates to the next step.

For ADR-41 C4 mesh-tying pass the tie residual query::

    status, n_aug, res = analyze_augmented(ops, query=ops.ladrunoMortarTieResidual, augTol=1e-12)

Returns ``(status, n_aug, measure)``:
    status 0  = converged (measure < augTol),
    status 1  = hit maxAug without converging,
    status <0 = a re-solve returned that nonzero analyze() code.
``measure`` is the engine convergence measure read by ``query`` (default the normal contact
penetration ``ops.ladrunoMortarPenetration()`` — the max KKT-active ||gbar||_inf).
"""


def analyze_augmented(ops, maxAug=20, augTol=1.0e-9, verbose=False, query=None):
    """Drive a mortar ALM measure to ``augTol`` by held-load Uzawa augmentation, WITHOUT recorder /
    load / time corruption (ADR-41 D1).

    Must be called AFTER the load has been applied and a first equilibrium reached (so the penalty
    pressure / active set are established). Switches the integrator to a zero-increment LoadControl
    and re-solves+commits until ``query() < augTol`` or ``maxAug`` augmentations run. The whole sweep
    is bracketed by ``ladrunoBeginAugment``/``ladrunoEndAugment`` so the held re-commits fire no
    recorders and advance no commitTag.

    ``query`` is the convergence measure (a zero-arg callable). Default = the NORMAL contact
    penetration ``ops.ladrunoMortarPenetration()`` (ADR-41 C2.2). For ADR-41 C4 mesh-tying pass
    ``query=ops.ladrunoMortarTieResidual`` (the weighted relative-displacement bond ‖r̄‖_∞) — the
    tie Uzawa (λ_tie ← λ_tie + epsTie·r/a, NO clamp) drives the SAME held-load loop.
    """
    if query is None:
        query = ops.ladrunoMortarPenetration
    # hold the external load: a zero-increment LoadControl re-solves at the SAME load level,
    # so the only thing that changes between augmentations is the committed multiplier lambda.
    ops.integrator("LoadControl", 0.0)
    # D1: mute recorders / freeze commitTag across the held-load sweep so the augmentation re-commits
    # produce NO spurious recorder samples and NO time advance. The begin() AND the first query() are
    # INSIDE the try so that ANY exception (including a throwing query) still hits `finally` — the
    # flag can never stay stuck ON / recorders muted for the rest of the session.
    try:
        ops.ladrunoBeginAugment()
        pen = query()
        for k in range(1, int(maxAug) + 1):
            rc = ops.analyze(1)
            if rc != 0:
                return rc, k, pen
            pen = query()
            if verbose:
                print(f"  [augment {k}] measure = {pen:.3e}")
            if pen < augTol:
                return 0, k, pen
        return 1, int(maxAug), pen
    finally:
        ops.ladrunoEndAugment()
