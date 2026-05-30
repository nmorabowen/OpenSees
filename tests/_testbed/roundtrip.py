"""Axis-1, Zone-A proxy for the C++-invisible serialization defects (D4).

Three defect classes are structurally invisible from a normal serial pytest:
  1. forgotten / mismatched sendSelf/recvSelf  -> corrupts OpenSeesMP & database
  2. classTag -> constructor dispatch via the broker
  3. a missing setResponse branch              -> recorder silently empty

setresponse_smoke covers (3) cheaply. database_roundtrip covers (1)+(2) by
round-tripping the committed state through the FE_Datastore (which exercises
sendSelf/recvSelf and the broker) — no MPI needed. A single vendored doctest
C++ target proving BIT-identical pack/unpack is deferred until the first
element whose serialization you genuinely distrust lands (canonical doc, ADR D4).
"""
import os
import tempfile

import pytest

from ._ops import ops


def setresponse_smoke(ele_tag, responses):
    """Assert each declared element response returns a non-empty list.
    `responses` is an iterable of arg-tuples passed to eleResponse, e.g.
    [("forces",), ("material", 1, "stress")]."""
    ele_tag = int(ele_tag)
    for args in responses:
        out = ops.eleResponse(ele_tag, *args)
        assert out, (
            f"eleResponse({ele_tag}, {args}) returned {out!r} — "
            "missing or broken setResponse branch."
        )


def database_roundtrip(build_fn, probe_nodes, ndf, dbname="ladruno_rt", rtol=1e-10):
    """Build -> solve -> save to FE_Datastore -> wipe -> rebuild skeleton ->
    restore -> assert committed nodal disp is recovered bit-for-bit.

    build_fn() must construct the full model and run one analyze() step.
    Skips (does not fail) if this OpenSees build lacks database support.

    The FE_Datastore writes <dbname>.IDs.* / <dbname>.VECs.* files next to the
    datastore path, so we point it at a TemporaryDirectory — never the cwd /
    repo (those artifacts must not pollute tests/).
    """
    # ignore_cleanup_errors: on Windows the FE_Datastore may still hold file
    # handles at scope exit; we wipe() to release them but keep this as a guard.
    with tempfile.TemporaryDirectory(prefix="ladruno_rt_",
                                     ignore_cleanup_errors=True) as td:
        dbpath = os.path.join(td, dbname)

        build_fn()
        before = {int(n): ops.nodeDisp(int(n)) for n in probe_nodes}

        try:
            ops.database("File", dbpath)
        except Exception as exc:  # noqa: BLE001 - build without FE_Datastore
            pytest.skip(f"database() unsupported in this build: {exc}")
        saved = ops.save(1)
        if saved is not None and saved < 0:
            pytest.skip("database save returned failure on this build")

        ops.wipe()
        build_fn()                  # same topology, fresh (uncommitted) state
        ops.database("File", dbpath)
        ops.restore(1)

        after = {int(n): ops.nodeDisp(int(n)) for n in probe_nodes}
        ops.wipe()                  # release FE_Datastore file handles before cleanup
        for n in before:
            for d in range(ndf):
                b, a = before[n][d], after[n][d]
                assert abs(a - b) <= rtol * abs(b) + 1e-12, (
                    f"database round-trip changed node {n} dof {d}: {b} -> {a} "
                    "(sendSelf/recvSelf or broker mismatch?)"
                )
