"""ExplicitBathe correctness pair (review 2026-07-01, PR 2).

1. update() misuse guard: the Noh-Bathe protocol is exactly ONE external update()
   per step (the second solve is internal). The old guard fired only at
   updateCount > 2, so `algorithm Newton` iterating twice silently advanced the
   domain time by an extra (1-p)*dt and overwrote A_tpdt with a stale re-solve.
   Post-fix the guard fires at > 1: the second update() fails LOUDLY.

2. -sms default lumping: the unified command defaulted the sizing lumping to
   RowSum while the deprecated alias names (ExplicitBatheSMS, ...) default to
   Diagonal ("matches the system Diagonal run"). On rotational-DOF (consistent-
   mass beam) models rowsum != diagonal, so the SAME model injected DIFFERENT
   scaling mass depending on which spelling the user typed. Post-fix, `-sms`
   without an explicit `-lump` defaults to Diagonal, matching the aliases (the
   bare non-SMS dt_cr estimate keeps its historical RowSum default).
"""
import math

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]


# --------------------------------------------------------------------------
# 1. update() called twice per step must fail loudly (guard at > 1, not > 2)
# --------------------------------------------------------------------------
def test_second_update_per_step_fails_loudly(capfd):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(1, 0.0, 0.0)
    ops.fix(1, 1, 1)
    ops.node(2, 0.0, 0.0)
    ops.fix(2, 0, 1)
    ops.mass(2, 1.0, 0.0)
    ops.uniaxialMaterial("Elastic", 1, 100.0)
    ops.element("zeroLength", 1, 1, 2, "-mat", 1, "-dir", 1)
    ops.setNodeVel(2, 1, 1.0, "-commit")
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("FixedNumIter", 2)                 # forces TWO update() calls/step
    ops.algorithm("Newton")
    ops.integrator("ExplicitBathe", 0.54)
    ops.analysis("Transient")
    capfd.readouterr()
    rc = ops.analyze(1, 1.0e-3)
    err = capfd.readouterr().err
    assert rc != 0, (
        "a second update() in one Noh-Bathe step must FAIL (got rc=0 — the old "
        ">2 guard let a 2-iteration Newton silently double-advance domain time)"
    )
    assert "update()" in err and "once" in err.lower(), (
        "the misuse warning must explain the one-update()-per-step protocol; "
        "stderr:\n%s" % err
    )


# --------------------------------------------------------------------------
# 2. unified -sms defaults to Diagonal sizing lumping, matching the alias
# --------------------------------------------------------------------------
def _beam_sms(integrator_args):
    """Consistent-mass beam (rowsum != diagonal on rotational DOFs) + SMS prime.
    Returns the injected nodal translational mass on the free node."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    ops.node(2, 1.0, 0.0, 0.0)
    ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
    ops.element("elasticBeamColumn", 1, 1, 2, 1.0, 200.0, 80.0, 1.0, 1.0, 1.0, 1,
                "-mass", 1.0, "-cMass")
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    ops.integrator(*integrator_args)
    ops.analysis("Transient")
    ops.analyze(1, 1.0e-6)                      # domainChanged -> injection
    return ops.nodeMass(2, 2)                   # injected uy mass (baseline 0? no
                                                # nodal mass commands -> injection only)


def test_unified_sms_lumping_default_matches_alias():
    dtT = 0.05                                  # above the beam's dt_e -> scaled
    dm_alias = _beam_sms(("ExplicitBatheSMS", 0.54, dtT))
    dm_unified = _beam_sms(("ExplicitBathe", 0.54, "-sms", dtT))
    assert dm_alias > 0.0, "fixture: the alias must inject scaling mass"
    assert dm_unified == pytest.approx(dm_alias, rel=1e-12), (
        "unified `-sms` (no -lump) injected %.6g but the deprecated alias "
        "injected %.6g — the sizing lumping defaults diverge (RowSum vs "
        "Diagonal); byte-identity per combo is broken on rotational-mass models"
        % (dm_unified, dm_alias)
    )
    # explicit -lump still wins (rowsum chosen on purpose must differ here)
    dm_rowsum = _beam_sms(("ExplicitBathe", 0.54, "-sms", dtT, "-lump", "rowsum"))
    assert dm_rowsum != pytest.approx(dm_alias, rel=1e-9), (
        "fixture: an EXPLICIT -lump rowsum should size differently on a "
        "consistent-mass beam (else this model cannot detect the default flip)"
    )
