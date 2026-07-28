"""InitialStateAnalysis on|off lifetime regression (backport of upstream 191c67c2d).

Pre-fix, `ops.InitialStateAnalysis(flag)` did

    Parameter *theP = new InitialStateParameter(...);
    theDomain->addParameter(theP);   // domain STORES the pointer (owns it)
    delete theP;                     // command frees it anyway

(`OPS_InitialStateAnalysis`, SRC/interpreter/OpenSeesMiscCommands.cpp) — the
domain's parameter container kept the freed pointer, and `wipe()` →
`Domain::clearAll` → `MapOfTaggedObjects::clearAll(true)` re-deleted it:
heap corruption (0xC0000374 / access-violation 0xC0000005) at wipe, plus a
second crash surface at interpreter teardown when `wipe()` was never called.
A second latent defect rode the same lines: `InitialStateParameter` has FIXED
tag 0, so a REPEAT call hit `Domain::addParameter`'s already-exists early
return *before* `setDomain` — the `ops_InitialStateAnalysis` global froze at
its first value ("off"-after-"on" never flipped it back; PM4Sand-family
`revertToStart` then silently kept state across `ops.reset()`).

The fix (verbatim backport of upstream 191c67c2d) uses a stack
`InitialStateParameter` + `setDomain()` — nothing persistent is registered:
  * `wipe()` and process exit are clean, and
  * repeat calls no longer collide, so the tag-0 "already exists" line must
    NOT appear and the flag genuinely toggles every call.

The child runs out-of-process (`-S` driver idiom, same as the other worktree
subprocess tests) so a pre-fix crash is contained in the subprocess and this
test just fails red. It only runs where a built engine exists
(dist/bin/opensees.pyd). Intentionally NO env-var gate: red until the fixed
pyd is built, green thereafter, red again on any regression.

Run:  py -3.12 -m pytest tests/test_initial_state_analysis_lifetime.py -x -q
"""
import os
import subprocess
import sys
from pathlib import Path

import pytest

# --------------------------------------------------------------------------
# bootstrap: THIS worktree's engine
# --------------------------------------------------------------------------
_DIST = str(Path(__file__).resolve().parents[1] / "dist" / "bin")
if not os.path.isfile(os.path.join(_DIST, "opensees.pyd")):
    pytest.skip(f"worktree engine not built: {_DIST}", allow_module_level=True)

pytestmark = [pytest.mark.zone_b]

# --------------------------------------------------------------------------
# out-of-process child: exercises both crash surfaces + the repeat-toggle path
# --------------------------------------------------------------------------
_CHILD = r'''
import os, sys
DIST = sys.argv[1]
os.add_dll_directory(DIST)
sys.path.insert(0, DIST)
os.environ["PATH"] = DIST + os.pathsep + os.environ.get("PATH", "")
import opensees as ops                 # -S => no boot-.pth stale-pyd preload
assert os.path.normcase(DIST) in os.path.normcase(ops.__file__), ops.__file__

def build_quad():
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(1, 0.0, 0.0); ops.node(2, 1.0, 0.0)
    ops.node(3, 1.0, 1.0); ops.node(4, 0.0, 1.0)
    ops.fix(1, 1, 1); ops.fix(2, 1, 1)
    ops.nDMaterial("ElasticIsotropic", 1, 25000.0, 0.25)
    ops.element("quad", 1, 1, 2, 3, 4, 1.0, "PlaneStrain", 1)

# both branches; "off" (revertToStart + parameter) first = the higher-risk path.
# The repeat call per branch locks in the repeat-toggle fix (no tag-0 collision).
for flag in ("off", "on"):
    build_quad()
    ops.InitialStateAnalysis(flag)
    ops.InitialStateAnalysis(flag)     # pre-fix: 'parameter with tag 0 already exists'
    # structural repeat-toggle guard, independent of opserr message wording:
    # post-backport NOTHING may stay registered in the parameter container.
    tags = ops.getParamTags()
    assert not tags, f"ISA left a Parameter registered in the domain: {tags}"
    ops.wipe()                          # pre-fix: heap corruption HERE

# at-exit surface: establish ISA state and do NOT wipe — the fixed build must
# also tear the interpreter down cleanly (nothing dangling in the domain).
build_quad()
ops.InitialStateAnalysis("on")
print("OK")
'''


def test_isa_survives_wipe_and_repeat(tmp_path):
    driver = tmp_path / "isa_lifetime_driver.py"
    driver.write_text(_CHILD, encoding="utf-8")  # child has non-ASCII comments
    env = dict(os.environ)
    env["PATH"] = _DIST + os.pathsep + env.get("PATH", "")
    try:
        proc = subprocess.run(
            [sys.executable, "-S", str(driver), _DIST],
            capture_output=True, text=True, timeout=120, env=env,
            encoding="utf-8", errors="replace")
    except subprocess.TimeoutExpired as exc:  # heap corruption can hang, not crash
        pytest.fail(
            f"ISA lifetime child hung (>120s) — likely heap corruption:\n{exc}\n"
            f"--- partial stdout ---\n{exc.stdout}\n--- partial stderr ---\n{exc.stderr}")

    combined = (proc.stdout or "") + "\n" + (proc.stderr or "")
    rc = proc.returncode

    # crash is contained in the child: pre-fix it dies at wipe() (or at exit)
    # with a Windows access-violation / heap code and never exits 0.
    assert rc == 0 and "OK" in proc.stdout, (
        f"ISA lifetime child died rc=0x{rc & 0xFFFFFFFF:08X} "
        f"(expect the UAF backport of upstream 191c67c2d in this pyd)\n"
        f"--- stdout ---\n{proc.stdout}\n--- stderr ---\n{proc.stderr}")

    # positive proof both command branches actually executed (not a silent no-op)
    assert ("InitialStateAnalysis ON" in combined
            and "InitialStateAnalysis OFF" in combined), (
        "ISA on/off banner missing — command did not run\n" + combined)

    # repeat-toggle fix: post-backport nothing persistent is registered, so the
    # duplicate InitialStateParameter (fixed tag 0) collision must NOT appear.
    assert "parameter with tag 0 already exists" not in combined, (
        "repeat ISA still registers a persistent InitialStateParameter "
        "(Domain::addParameter tag-0 collision) — repeat-toggle fix missing\n"
        + combined)
