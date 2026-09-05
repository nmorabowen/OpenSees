"""Broker/database round-trip characterization for the vanilla
`ManzariDafaliasPlaneStrain` null-constructor classTag defect.

See `Ladruno_implementation/LEDGER_quirks.md` -- "`ManzariDafaliasPlaneStrain`'s
null constructor sets the WRONG classTag, and a broker restore never repairs
it" -- for the full writeup.  Short version:

    ManzariDafalias3D::ManzariDafalias3D()                     -- CORRECT
      : ManzariDafalias(ND_TAG_ManzariDafalias3D)
    ManzariDafaliasPlaneStrain::ManzariDafaliasPlaneStrain()   -- WRONG
      : ManzariDafalias()   // -> NDMaterial(0, ND_TAG_ManzariDafalias)

`ManzariDafalias`'s bare null constructor hardcodes the DIMENSIONLESS base's
own classTag, so a null-constructed `ManzariDafaliasPlaneStrain` -- which is
exactly what `FEM_ObjectBrokerAllClasses.cpp:2548` builds on every database or
parallel restore -- reports `getClassTag() == ND_TAG_ManzariDafalias` for the
rest of its life. `recvSelf` never repairs it (there is no `setClassTag`
anywhere in `SRC/` -- `MovableObject::classTag` is private).

WHY THIS IS A CHARACTERIZATION TEST, NOT A FIX GATE.  The owner explicitly
decided (2026-08-28, recorded in the ledger entry above) NOT to fix this:
the only base-class constructor that takes a classTag also resets the
process-wide, static `mElastFlag` stage flag as a side effect
(`ManzariDafalias.cpp:348`) -- a live behaviour change on every broker/MP
restore of ANY Manzari-family material, not just this one -- and the only
side-effect-free fix would add a new constructor overload to the vanilla
`ManzariDafalias` API. Both were weighed and declined; the defect is narrow
(it only bites a restored object that is later RE-serialised). ADR-90 WP-B
(2026-09-04) reconfirmed this reasoning still holds -- see the ledger entry's
updated status line -- so this test documents and reproduces the defect
end-to-end rather than gating a fix that was deliberately not made.

WHY IT TAKES TWO ROUND TRIPS.  After the FIRST restore, the defect is
invisible: `recvSelf` restores the full 97-slot state correctly, so the
element behaves normally even though the restored object's internal classTag
is wrong. The failure only appears once that mistagged object is asked to
`sendSelf` ITS OWN (wrong) classTag -- i.e. on a SECOND save -- and the
broker then builds a bare `ManzariDafalias` (not `-PlaneStrain`) on the
SECOND restore.

WHY THIS RUNS IN A SUBPROCESS, WITH `-u` AND MERGED STREAMS.
`ManzariDafalias::getCopy(void)` / `getType()` / `getOrder()` (the base
class's own "subclass responsibility" overrides) all call `exit(-1)` if
actually reached; an in-process assertion could never observe that failure,
so the interpreter must be a child we can inspect after the fact. What is
actually reached today is `NDMaterial`'s own generic
`setTrialStrain`/`getStress` defaults, which print the identical "subclass
responsibility" phrase and return failure codes rather than aborting.  The
assertions below need to know WHERE in the sequence that phrase appears (it
must be absent through the first restore, present only after the second), so
the child is run unbuffered (`-u`) with stderr merged into the SAME stream as
stdout (`stderr=subprocess.STDOUT`, via `run_python_script(...,
merge_stderr=True)`) -- MEASURED to matter: without `-u`, Python's stdout is
fully buffered once redirected to a pipe while C++'s `opserr` is not, so
"subclass responsibility" lines appear to precede `print()` calls that
actually executed first (a pure buffering artifact, reproduced and diagnosed
while writing this test). Without `merge_stderr=True`,
`subprocess.run(capture_output=True)` hands back stdout and stderr as two
SEPARATE strings; concatenating them (`stdout + stderr`, the convention this
file used before) puts EVERY stdout marker before EVERY stderr line
regardless of when either was actually written, which would make an
ordering assertion here vacuously true no matter what really happened.

THE LOAD-BEARING SIGNAL IS THE STRESS VECTOR LENGTH, NOT JUST THE WARNING
TEXT. `ops.eleResponse(1, 'material', 1, 'stress')` on the plane-strain
element must be a 3-component vector (`{sigma_xx, sigma_yy, sigma_xy}`) as
long as the Gauss-point material is really a `ManzariDafaliasPlaneStrain`.
Once the SECOND restore broker-builds a bare `ManzariDafalias` instead, that
response comes back some OTHER size (measured: length 1, `[0.0]`) because
`NDMaterial`'s generic default has no notion of the plane-strain order. That
length mismatch is a structural, message-independent witness of the defect
-- it does not depend on "subclass responsibility" appearing verbatim (a
future OpenSees version could reword the message) or on a crash (rc == 0 is
exactly what was measured; nothing here forces or requires a process exit).
"""
import ast
import os

import pytest

from _testbed.subprocess_run import run_python_script

_HERE = os.path.dirname(os.path.abspath(__file__))

# Gorini's calibrated set (test_ladruno_sanisand.py / test_fspm_over_manzari_family.py
# carry the same list); values are irrelevant to the classTag defect itself.
_PARAMS = [264.32, 0.3129, 0.6944, 1.33090, 0.71, 0.027, 0.83, 0.45, 101.0,
           0.005, 1.3, 0.968, 3.5, 0.05, 5.75, 12.5, 1100.0, 2.0]
_TAIL = [1, 0, 1, 1.0e-10, 1.0e-10]          # ModifiedEuler, elastic tangent, analytic Jacobian

pytestmark = [pytest.mark.zone_a, pytest.mark.t0m]

# Ordered markers the child MUST print, in this order, before the double
# round trip is considered to have run at all.  `_split_segments` uses these
# to slice the (chronologically faithful, thanks to `-u` + merge_stderr) child
# output into the two windows the two halves of the characterization need.
_MARKERS = ('BUILT', 'SAVE1_OK', 'RESTORE1_OK', 'STRESS_AFTER_RESTORE1',
            'SAVE2_OK', 'RESTORE2_OK', 'STRESS_AFTER_RESTORE2', 'ALL_DONE')

_SCRIPT = r'''
import sys
sys.path.insert(0, {here!r})
from _testbed import ops
P = {params!r}; T = {tail!r}
ops.wipe(); ops.model('basic', '-ndm', 2, '-ndf', 2)
for n, (x, y) in enumerate([(0,0),(1,0),(1,1),(0,1)], start=1):
    ops.node(n, x, y)
ops.fix(1, 1, 1); ops.fix(2, 0, 1); ops.fix(4, 1, 0)
ops.nDMaterial('ManzariDafalias', 1, *P, *T)
ops.element('quad', 1, 1,2,3,4, 1.0, 'PlaneStrain', 1)
print('BUILT')

dbpath = {dbpath!r}
try:
    ops.database('File', dbpath)
except Exception as exc:  # noqa: BLE001 - build without FE_Datastore
    print('DB_UNSUPPORTED', str(exc))
    sys.exit(0)

saved1 = ops.save(1)
if saved1 is not None and saved1 < 0:
    print('DB_SAVE_FAILED')
    sys.exit(0)
print('SAVE1_OK')

ops.wipe()
ops.database('File', dbpath)
ops.restore(1)
print('RESTORE1_OK')
s1 = ops.eleResponse(1, 'material', 1, 'stress')
print('STRESS_AFTER_RESTORE1', repr(list(s1)))

# The mistagged (but data-correct) object re-serialises ITS OWN classTag now.
ops.save(2)
print('SAVE2_OK')

ops.wipe()
ops.database('File', dbpath)
ops.restore(2)
print('RESTORE2_OK')

s2 = ops.eleResponse(1, 'material', 1, 'stress')
print('STRESS_AFTER_RESTORE2', repr(list(s2)))
ops.wipe()  # release FE_Datastore handles (Windows) before tmp cleanup
print('ALL_DONE')
'''


def _run(script):
    return run_python_script(script, cwd=_HERE, timeout=300, merge_stderr=True)


def _segment(out, start_marker, end_marker):
    """The substring of `out` strictly between two markers, in that order.

    Requires BOTH markers to be present, in that order -- callers use this
    only after already asserting the marker sequence, so a `ValueError` here
    means the marker check above is out of sync with this helper, not that
    the child misbehaved.
    """
    after_start = out.split(start_marker, 1)[1]
    return after_start.split(end_marker, 1)[0]


def test_manzari_planestrain_classtag_survives_one_roundtrip_but_not_two(tmp_path):
    """One round trip: invisible (data restores fine, response shape intact).
    Two round trips: the broker builds a bare `ManzariDafalias`, and its
    'stress' response comes back the wrong size."""
    dbpath = str(tmp_path / 'manzari_ps_classtag_rt')
    rc, out = _run(_SCRIPT.format(here=_HERE, params=_PARAMS, tail=_TAIL, dbpath=dbpath))

    if 'DB_UNSUPPORTED' in out:
        pytest.skip('database() unsupported in this build:\n' + out)
    if 'DB_SAVE_FAILED' in out:
        pytest.skip('database save returned failure on this build:\n' + out)

    missing = [m for m in _MARKERS if m not in out]
    if rc != 0 and missing:
        pytest.fail(
            'the double round trip did not complete: child exited with code '
            f'{rc} before printing {missing}. This is a crash somewhere OTHER '
            'than where the classTag defect is expected to bite (RESTORE2 '
            'onward) -- diagnose from the tail of the output below rather '
            f'than assuming it confirms the defect.\n--- output tail ---\n{out[-4000:]}'
        )
    assert not missing, (
        f'child output is missing marker(s) {missing} with rc={rc} -- the '
        f'deck itself no longer builds/saves/restores cleanly.\n{out}'
    )
    # From here every marker is present, IN ORDER (each `_segment`/`.index`
    # call below would raise if not, which is itself a hard failure).
    assert [out.index(m) for m in _MARKERS] == sorted(out.index(m) for m in _MARKERS), (
        'markers are present but out of order -- the child printed something '
        'unexpected between them', out)

    before_restore2 = _segment(out, 'RESTORE1_OK', 'SAVE2_OK')
    assert 'subclass responsibility' not in before_restore2, (
        'the FIRST restore already shows "subclass responsibility" -- the '
        'defect this test characterizes is supposed to be invisible after '
        'exactly one round trip (recvSelf restores the 97-slot state '
        'correctly regardless of the internal classTag). Either the wire '
        'format changed or the broker dispatch changed.', out)

    stress1_line = [l for l in out.splitlines() if l.startswith('STRESS_AFTER_RESTORE1')][0]
    stress1 = ast.literal_eval(stress1_line.split(' ', 1)[1])
    assert len(stress1) == 3, (
        'STRESS_AFTER_RESTORE1 is not a 3-component plane-strain vector -- '
        'the FIRST restore already produced the wrong material shape, which '
        'is not the defect this test characterizes (see module docstring)',
        stress1, out)

    after_restore2 = _segment(out, 'RESTORE2_OK', 'ALL_DONE')
    stress2_line = [l for l in out.splitlines() if l.startswith('STRESS_AFTER_RESTORE2')][0]
    stress2 = ast.literal_eval(stress2_line.split(' ', 1)[1])

    shape_broke = len(stress2) != 3
    warned = 'subclass responsibility' in after_restore2

    # The load-bearing assertion. Either discriminator alone is sufficient
    # evidence (see the module docstring for why the shape check is the one
    # that cannot be talked out of by a future reworded message); if NEITHER
    # fires, and every marker above was present with rc == 0, the defect has
    # plausibly been fixed (or the broker/element path changed) -- fail
    # loudly with the tail rather than passing silently on a stale claim.
    assert shape_broke or warned, (
        'after the SECOND restore, the "stress" response is still a 3-vector '
        'and no "subclass responsibility" warning appeared anywhere after '
        'RESTORE2_OK. Either the classTag defect has been fixed (update '
        'LEDGER_quirks.md and this test to say so) or the broker/element '
        f'path changed and no longer reproduces it.\nrc={rc}\n'
        f'STRESS_AFTER_RESTORE1={stress1} STRESS_AFTER_RESTORE2={stress2}\n'
        f'--- output tail ---\n{out[-4000:]}')
