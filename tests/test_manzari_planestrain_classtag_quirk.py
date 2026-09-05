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

WHY THIS RUNS IN A SUBPROCESS.  `ManzariDafalias::getCopy(void)` / `getType()`
/ `getOrder()` (the base class's own "subclass responsibility" overrides) all
call `exit(-1)` if actually reached. If some future OpenSees change starts
calling one of those during element `recvSelf` or `Response` construction (it
does not today -- see below), an in-process assertion could never observe the
failure; the interpreter would simply be gone. What IS reached today is
`NDMaterial`'s own generic `setTrialStrain`/`getStress`/`getTangent`
defaults, which print the identical "subclass responsibility" phrase and
return failure codes rather than aborting -- still visible, but only in a
process we can inspect after the fact regardless of which of the two
"subclass responsibility" families ends up firing.
"""
import os
import subprocess
import sys

import pytest

_HERE = os.path.dirname(os.path.abspath(__file__))

# Gorini's calibrated set (test_ladruno_sanisand.py / test_fspm_over_manzari_family.py
# carry the same list); values are irrelevant to the classTag defect itself.
_PARAMS = [264.32, 0.3129, 0.6944, 1.33090, 0.71, 0.027, 0.83, 0.45, 101.0,
           0.005, 1.3, 0.968, 3.5, 0.05, 5.75, 12.5, 1100.0, 2.0]
_TAIL = [1, 0, 1, 1.0e-10, 1.0e-10]          # ModifiedEuler, elastic tangent, analytic Jacobian

pytestmark = [pytest.mark.zone_a, pytest.mark.t0m]

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
'''


def _run(script):
    p = subprocess.run([sys.executable, '-c', script], cwd=_HERE,
                       capture_output=True, text=True, timeout=300)
    return p.returncode, p.stdout + p.stderr


def test_manzari_planestrain_classtag_survives_one_roundtrip_but_not_two(tmp_path):
    """One round trip: invisible (data restores fine). Two round trips: the
    broker builds a bare `ManzariDafalias` and something downstream notices.
    """
    dbpath = str(tmp_path / 'manzari_ps_classtag_rt')
    rc, out = _run(_SCRIPT.format(here=_HERE, params=_PARAMS, tail=_TAIL, dbpath=dbpath))

    if 'DB_UNSUPPORTED' in out:
        pytest.skip('database() unsupported in this build:\n' + out)
    if 'DB_SAVE_FAILED' in out:
        pytest.skip('database save returned failure on this build:\n' + out)

    assert 'BUILT' in out, out
    assert 'SAVE1_OK' in out, out
    assert 'RESTORE1_OK' in out, out

    before_restore2 = out.split('RESTORE1_OK', 1)[1].split('RESTORE2_OK', 1)[0]
    assert 'subclass responsibility' not in before_restore2, (
        'the FIRST restore already shows "subclass responsibility" -- the '
        'defect this test characterizes is supposed to be invisible after '
        'exactly one round trip (recvSelf restores the 97-slot state '
        'correctly regardless of the internal classTag). Either the wire '
        'format changed or the broker dispatch changed.', out)

    # The load-bearing assertion: the SECOND round trip must show the defect
    # firing, one way or another -- a process exit (rc != 0) if some future
    # code path reaches ManzariDafalias::getCopy/getType/getOrder (the
    # exit(-1) family), or the "subclass responsibility" phrase from
    # NDMaterial's own generic setTrialStrain/getStress/getTangent defaults
    # (what actually fires today, via `eleResponse ... stress`).
    assert rc != 0 or 'subclass responsibility' in out, (
        'the double round trip produced neither a crash nor a "subclass '
        'responsibility" warning after RESTORE2. Either the classTag defect '
        'has been fixed (update LEDGER_quirks.md and this test to say so) or '
        'the broker/element path changed and no longer reproduces it.',
        rc, out)
