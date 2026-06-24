"""ADR-41 D1 — within-step held-load augmentation sign-off (the no-corruption first-class recipe).

C2.2 (#375) and C4 (#381) ship the commit-cycle ALM: ONE Uzawa augmentation per physical step, so
the penetration / tie residual reaches its epsN-INDEPENDENT value only across MANY steps. D1 turns
the held-load `analyze_augmented` proc into a FIRST-CLASS tested recipe: a within-step augmentation
loop that drives ‖ḡ‖_∞ → augTol (contact) / ‖r̄‖_∞ → augTol (tie) INSIDE a single physical step,
WITHOUT recorder / load / time corruption. The proc brackets its held-load re-commits with
ladrunoBeginAugment/ladrunoEndAugment, so Domain::commit() performs the Uzawa λ update but fires NO
recorders and bumps NO commitTag.

Gates (capstone D1 row):
  1. WITHIN-STEP CONVERGENCE — at a HELD load the loop drives the penetration (contact) and the tie
     residual (tie) below augTol inside ONE step, epsN-independent at finite epsN.
  2. NO RECORDER CORRUPTION — a model with a recorder produces IDENTICAL recorder output (same sample
     count + same time sequence) with vs without within-step augmentation; the unbracketed loop is
     shown to corrupt it (the NEGATIVE control), proving the bracket is what fixes it.
  3. NO LOAD CORRUPTION — the held re-solves do not re-apply / double-count the load (the converged
     contact reaction == the single-application load; the load factor / time is not advanced).
  4. OPT-IN BYTE-IDENTITY — a begin/end no-op pair leaves the recorder output identical to never
     touching augmentation (the flag, once ended, restores the stock commit path).

Mechanics validated oracle-first in contact_prototypes/proto_d1_augmentation.py (11/11).
"""
import os
import sys
import tempfile

import pytest

from _testbed import ops

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "Ladruno_scripts"))
from analyze_augmented import analyze_augmented  # noqa: E402

pytestmark = [pytest.mark.zone_a]

PTOT = 400.0


def _unit_square(tags, z):
    c = [(tags[0], 0.0, 0.0, z), (tags[1], 1.0, 0.0, z),
         (tags[2], 1.0, 1.0, z), (tags[3], 0.0, 1.0, z)]
    for t, x, y, zz in c:
        ops.node(t, x, y, zz)


def _static_setup():
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-11, 50, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")


def _build_contact(epsN, series="Linear"):
    """One 4-node slave facet (seeded penetrating) on a coincident FIXED master quad, under a uniform
    downward pressure. Returns the slave node tags (5..8)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    _unit_square([1, 2, 3, 4], 0.0)
    for t in (1, 2, 3, 4):
        ops.fix(t, 1, 1, 1)
    _unit_square([5, 6, 7, 8], -1.0e-4)
    for t in (5, 6, 7, 8):
        ops.fix(t, 1, 1, 0)
    ops.contactSurface(1, "-master", 4, 1, 2, 3, 4)
    ops.contactSurface(2, "-slave-segments", 4, 5, 6, 7, 8)
    ops.contact(1, 1, 2, "-mortar", "-epsN", epsN, "-outward", 0.0, 0.0, 1.0)
    ops.timeSeries(series, 1)
    ops.pattern("Plain", 1, 1)
    for t in (5, 6, 7, 8):
        ops.load(t, 0.0, 0.0, -PTOT / 4.0)
    _static_setup()
    return [5, 6, 7, 8]


def _build_tie(epsTie):
    """A NON-MATCHING tied interface: a fixed master quad at z=0 + a free slave quad seeded just off
    it, pulled UP by a tension load. The tie binds the slave to the master (r → 0). Returns slave
    tags (5..8)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    _unit_square([1, 2, 3, 4], 0.0)
    for t in (1, 2, 3, 4):
        ops.fix(t, 1, 1, 1)
    _unit_square([5, 6, 7, 8], 1.0e-4)         # slave just ABOVE the master (a gap a contact would open)
    for t in (5, 6, 7, 8):
        ops.fix(t, 1, 1, 0)
    ops.contactSurface(1, "-master", 4, 1, 2, 3, 4)
    ops.contactSurface(2, "-slave-segments", 4, 5, 6, 7, 8)
    ops.contact(1, 1, 2, "-mortar", "-tie", "-epsTie", epsTie, "-outward", 0.0, 0.0, 1.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t in (5, 6, 7, 8):
        ops.load(t, 0.0, 0.0, +PTOT / 4.0)     # tension AWAY from master — only the tie holds it
    _static_setup()
    return [5, 6, 7, 8]


# ============================================================ Gate 1: within-step convergence
def test_d1_within_step_contact_convergence():
    """At a HELD load the augmentation loop drives the contact penetration below augTol INSIDE one
    physical step, and the converged value is epsN-INDEPENDENT (penalty alone is O(1/epsN))."""
    augTol = 1.0e-9
    for epsN in (1.0e4, 1.0e6):
        _build_contact(epsN)
        ops.integrator("LoadControl", 1.0)
        assert ops.analyze(1) == 0, "penalty load step did not converge"
        penalty_pen = ops.ladrunoMortarPenetration()
        status, n_aug, alm_pen = analyze_augmented(ops, maxAug=60, augTol=augTol)
        assert status == 0, f"epsN={epsN:.0e}: within-step augmentation failed (pen {alm_pen:.2e})"
        assert alm_pen < augTol, f"epsN={epsN:.0e}: penetration {alm_pen:.2e} !< augTol within one step"
        # the penalty leaves a real penetration that ALM removes (the within-step win is non-trivial)
        assert penalty_pen > alm_pen, "no penetration to converge (test is vacuous)"


def test_d1_within_step_tie_convergence():
    """At a HELD load the augmentation loop drives the TIE residual ‖r̄‖_∞ below augTol inside one
    physical step (the tie Uzawa, no-clamp), epsTie-INDEPENDENT."""
    augTol = 1.0e-11
    for epsTie in (1.0e6, 1.0e8):
        _build_tie(epsTie)
        ops.integrator("LoadControl", 1.0)
        assert ops.analyze(1) == 0, "penalty tie step did not converge"
        penalty_res = ops.ladrunoMortarTieResidual()
        status, n_aug, alm_res = analyze_augmented(
            ops, maxAug=80, augTol=augTol, query=ops.ladrunoMortarTieResidual)
        assert status == 0, f"epsTie={epsTie:.0e}: tie augmentation failed (‖r‖ {alm_res:.2e})"
        assert alm_res < augTol, f"epsTie={epsTie:.0e}: ‖r‖ {alm_res:.2e} !< augTol within one step"
        assert penalty_res > alm_res, "no tie residual to converge (test is vacuous)"


# ============================================================ Gate 2: no recorder corruption
def _record_times(with_aug, bracket=True, nsteps=4, epsN=1.0e4):
    """Run an nsteps LoadControl analysis on the contact model with a Node disp recorder. Optionally
    run a within-step augmentation each step (bracketed by the proc, or — for the negative control —
    UNbracketed raw held-load re-commits). Returns the list of recorded time-column values."""
    fd, path = tempfile.mkstemp(suffix=".out", prefix="d1_rec_")
    os.close(fd)
    try:
        _build_contact(epsN)
        ops.recorder("Node", "-file", path, "-time", "-node", 5, "-dof", 3, "disp")
        for _ in range(nsteps):
            ops.integrator("LoadControl", 1.0)
            assert ops.analyze(1) == 0
            if with_aug and bracket:
                analyze_augmented(ops, maxAug=40, augTol=1.0e-10)
            elif with_aug and not bracket:
                # NEGATIVE control: the pre-D1 raw held-load loop — NO begin/end bracket, so every
                # held re-commit fires the recorder (the corruption D1 fixes).
                ops.integrator("LoadControl", 0.0)
                for _ in range(3):
                    assert ops.analyze(1) == 0
        ops.wipe()                       # flush + close the recorder file
        with open(path) as f:
            return [float(line.split()[0]) for line in f if line.strip()]
    finally:
        if os.path.exists(path):
            os.remove(path)


def test_d1_no_recorder_corruption():
    """A model with a recorder produces IDENTICAL recorder output (same sample count + same time
    sequence) WITH vs WITHOUT within-step augmentation — the augmentation passes fire no recorders."""
    base = _record_times(with_aug=False)
    augmented = _record_times(with_aug=True, bracket=True)
    assert len(base) == 4, f"baseline recorded {len(base)} samples (expected 4 = one per step)"
    assert len(augmented) == len(base), (
        f"within-step augmentation corrupted the sample COUNT: {len(augmented)} vs {len(base)}")
    assert augmented == base, (
        f"within-step augmentation corrupted the TIME sequence: {augmented} vs {base}")


def test_d1_recorder_corruption_without_bracket():
    """NEGATIVE control: the SAME held-load re-commits WITHOUT the begin/end bracket DO inject
    spurious recorder samples (3 extra per step) — proving the bracket is what prevents corruption,
    not the zero-load increment alone."""
    base = _record_times(with_aug=False)
    raw = _record_times(with_aug=True, bracket=False)
    assert len(raw) > len(base), (
        f"unbracketed held-load loop did NOT corrupt the recorder ({len(raw)} vs {len(base)}) — "
        "the negative control is vacuous, so the positive gate proves nothing")
    # 4 real steps + 3 held re-commits each = 4 + 12 = 16 spurious-inclusive samples
    assert len(raw) == 4 + 4 * 3, f"unexpected corrupted sample count {len(raw)}"


# ============================================================ Gate 3: no load corruption
def test_d1_no_load_corruption():
    """The held re-solves re-APPLY the same load factor (a zero-increment LoadControl SETS, not
    accumulates), so the converged contact reaction == the single-application load and the load
    factor / pseudo-time is NOT advanced by the augmentation passes."""
    _build_contact(epsN=1.0e4)
    ops.integrator("LoadControl", 1.0)
    assert ops.analyze(1) == 0
    assert ops.getTime() == 1.0
    status, n_aug, _ = analyze_augmented(ops, maxAug=60, augTol=1.0e-10)
    assert status == 0
    # the load factor / pseudo-time is unchanged by the whole sweep (no time advance).
    assert ops.getTime() == 1.0, f"augmentation advanced the load factor to {ops.getTime()}"
    assert n_aug >= 1, "augmentation loop did not run (gate is vacuous)"
    # the converged contact reaction balances exactly the SINGLE applied load (not 2×): the slave
    # free-DOF reaction carries the contact traction; its sum == |applied total|.
    ops.reactions()
    contact_reaction = sum(ops.nodeReaction(t, 3) for t in (5, 6, 7, 8))
    assert abs(contact_reaction - PTOT) < 1.0e-6 * PTOT, (
        f"contact reaction {contact_reaction:.4f} != single applied load {PTOT} "
        "(load double-counted by the held re-solves?)")

    # ADVERSARIAL discriminator: drive MANY more held re-solves at the same load. If a re-solve
    # accumulated (re-applied) the load pattern, the load factor would advance and/or the reaction
    # would grow with each pass. Both must stay PINNED — proving the zero-increment LoadControl SETS,
    # never accumulates (the actual "no load corruption" guarantee, not just the converged snapshot).
    ops.ladrunoBeginAugment()
    try:
        ops.integrator("LoadControl", 0.0)
        for _ in range(10):
            assert ops.analyze(1) == 0
            assert ops.getTime() == 1.0, f"a held re-solve advanced the load factor to {ops.getTime()}"
    finally:
        ops.ladrunoEndAugment()
    ops.reactions()
    reaction_after_many = sum(ops.nodeReaction(t, 3) for t in (5, 6, 7, 8))
    assert abs(reaction_after_many - PTOT) < 1.0e-6 * PTOT, (
        f"reaction grew to {reaction_after_many:.4f} over 10 extra held re-solves — the load is "
        "being accumulated, not held")


# ============================================================ Gate 4: opt-in byte-identity
def _record_times_noop_bracket(nsteps=4, epsN=1.0e4):
    """Run the plain analysis but wrap each step in a begin/end no-op pair (no held re-commits in
    between). The flag is toggled on then immediately off, so the recorded commit must be IDENTICAL
    to never touching augmentation — proving ladrunoEndAugment restores the stock commit path."""
    fd, path = tempfile.mkstemp(suffix=".out", prefix="d1_noop_")
    os.close(fd)
    try:
        _build_contact(epsN)
        ops.recorder("Node", "-file", path, "-time", "-node", 5, "-dof", 3, "disp")
        for _ in range(nsteps):
            ops.integrator("LoadControl", 1.0)
            ops.ladrunoBeginAugment()
            ops.ladrunoEndAugment()          # ended BEFORE the real commit ⇒ stock path records it
            assert ops.analyze(1) == 0
        ops.wipe()
        with open(path) as f:
            return [line.strip() for line in f if line.strip()]
    finally:
        if os.path.exists(path):
            os.remove(path)


def test_d1_opt_in_byte_identity():
    """A begin/end no-op pair (flag toggled on then off before the real commit) leaves the recorder
    output BYTE-IDENTICAL to a run that never touches augmentation — the opt-in is clean: when the
    flag is OFF the commit path is the stock path."""
    fd, path = tempfile.mkstemp(suffix=".out", prefix="d1_plain_")
    os.close(fd)
    try:
        _build_contact(epsN=1.0e4)
        ops.recorder("Node", "-file", path, "-time", "-node", 5, "-dof", 3, "disp")
        for _ in range(4):
            ops.integrator("LoadControl", 1.0)
            assert ops.analyze(1) == 0
        ops.wipe()
        with open(path) as f:
            plain = [line.strip() for line in f if line.strip()]
    finally:
        if os.path.exists(path):
            os.remove(path)
    noop = _record_times_noop_bracket()
    assert noop == plain, (
        "begin/end no-op pair changed the recorder output — the flag does not cleanly restore the "
        f"stock commit path:\n  plain={plain}\n  noop ={noop}")
