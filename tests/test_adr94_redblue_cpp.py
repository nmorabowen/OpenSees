"""ADR-94 R2 red lane 1 (C++ state machine / API contract) -- cheap,
source-level reproducers for findings not covered by R1-A/R1-B.

These mirror the structural-regex style already used for H2/H14 in
``test_adr94_hlist_mechanical.py`` (no build/runtime needed for a
class-static-sharing or dead-code claim -- a source read settles it, and
the regex flips red the moment someone fixes the underlying code). Zone-A,
source reads only, sub-second.

See ``Ladruno_implementation/_adr94_redblue/red1_cpp.md`` for the full
analysis (F1-F6).
"""
import os
import re

import pytest

pytestmark = [pytest.mark.zone_a]

_SRC_DIR = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), os.pardir,
    "SRC", "material", "nD", "ASDPlasticMaterial3D")

_MAIN_H = os.path.join(_SRC_DIR, "ASDPlasticMaterial3D.h")
_DP_YF = os.path.join(_SRC_DIR, "YieldFunctions", "DruckerPrager_YF.h")
_DOMAIN_CPP = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), os.pardir,
    "SRC", "domain", "domain", "Domain.cpp")


def _read(path):
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        return f.read()


def test_F1_getInitialTangent_writes_shared_static_stiffness():
    """FIXED by wp/94b.  ``getInitialTangent()`` used to be a "getter" that
    assigned the class-static ``Stiffness`` (the same static H1 showed is
    shared across every GP/element/tag of one YF x PF x EL combo) before
    copying to its own return buffer, so any ``-initial`` algorithm or
    diagnostic call silently overwrote the tangent every other instance's
    later ``getTangent()`` would read.  It now computes into a local and
    copies THAT; no member state is touched."""
    src = _read(_MAIN_H)
    m = re.search(r"getInitialTangent\(\s*\)\s*\{(.*?)\n    \}", src, re.S)
    assert m, "could not locate getInitialTangent() body"
    body = m.group(1)
    live = "\n".join(ln for ln in body.splitlines()
                     if ln.strip() and not ln.strip().startswith("//"))
    assert "Stiffness" not in live, (
        "getInitialTangent() mutates Stiffness again -- F1 "
        "(getter-with-side-effect) is back:\n" + live)
    assert "copyToMatrixReference(Eelastic, return_matrix);" in live, (
        "getInitialTangent() no longer copies its LOCAL elastic tangent to the "
        "return buffer; re-check the F1 fix")


def test_F3_sendSelf_recvSelf_report_success_while_doing_nothing():
    """sendSelf/recvSelf print 'not implemented' but return 0 (success) --
    unlike revertToStart(), which at least returns -1. A parallel/database
    caller has no signal that state was not transferred."""
    src = _read(_MAIN_H)
    for name in ("sendSelf", "recvSelf"):
        m = re.search(name + r"\([^)]*\)\s*\{(.*?)\n    \}", src, re.S)
        assert m, f"could not locate {name}() body"
        body = m.group(1)
        assert "not implemented" in body, f"{name}() no longer self-reports as a stub"
        # the only return statement in the (short) stub body must be `return 0`
        returns = re.findall(r"return\s+(-?\d+)\s*;", body)
        assert returns == ["0"], (
            f"{name}() returns {returns}, expected exactly ['0'] (false-success "
            "contract) -- if this is now -1, F3 is fixed; update red1_cpp.md")


def test_F4_domain_reverttostart_discards_element_return_code():
    """Domain::revertToStart() calls elePtr->revertToStart() in a bare
    statement, discarding its return -- the H4 swallow originates here, one
    layer below OPS_resetModel (which R1-B cited)."""
    src = _read(_DOMAIN_CPP)
    m = re.search(r"Domain::revertToStart\(void\)\s*\{(.*?)\n\}", src, re.S)
    assert m, "could not locate Domain::revertToStart() body"
    body = m.group(1)
    # the element loop must call revertToStart() without capturing/checking it
    assert re.search(r"elePtr->revertToStart\(\)\s*;", body), (
        "Domain::revertToStart() element loop changed shape; re-check whether "
        "the return code is now inspected (F4 would then be fixed)")
    assert not re.search(r"(if|while)\s*\([^)]*elePtr->revertToStart\(\)",
                          body), (
        "Domain::revertToStart() now appears to branch on "
        "elePtr->revertToStart()'s return value -- F4 may be fixed")


def test_F2_yf_and_pf_functors_return_class_static_buffers():
    """FIXED by wp/94b.  Every YF/PF functor used to return its VoigtVector
    results through a private ``static VoigtVector vv_out`` scoped to the
    FUNCTOR's own template parameters -- a strictly wider sharing key than
    ``Stiffness`` (H1), since two unrelated ``ASDPlasticMaterial3D<E,Y,P,tag>``
    combos that merely reuse the same YF or PF type shared one buffer.  They
    are per-instance ``mutable`` members now.  Checked across the whole
    YieldFunctions/ + PlasticFlowDirections/ tree, not just DruckerPrager."""
    for sub in ("YieldFunctions", "PlasticFlowDirections"):
        d = os.path.join(_SRC_DIR, sub)
        for fn in sorted(os.listdir(d)):
            if not fn.endswith(".h"):
                continue
            src = _read(os.path.join(d, fn))
            assert not re.search(r"static\s+VoigtVector\s+(vv_out|result)\s*;",
                                 src), (
                f"{sub}/{fn} declares a static return buffer again -- F2 is back")
    # the DruckerPrager buffer specifically is still written by both the
    # deviatoric-direction getter and apex_stress, i.e. the fix moved the
    # buffer, it did not delete the code paths.
    src = _read(_DP_YF)
    assert re.search(r"mutable\s+VoigtVector\s+vv_out\s*=", src), (
        "DruckerPrager_YF's vv_out is no longer a per-instance member")
    assignments = re.findall(r"vv_out\s*=", src)
    assert len(assignments) >= 3, (
        f"expected the declaration plus >=2 writes into vv_out (df_dsigma_ij + "
        f"apex_stress), found {len(assignments)}")
