"""ADR-85 T2 -- the 2D removal gate, INCLUDING the concave-vertex use-after-
free case (G-T2(f), deferral D5).

WHAT THIS FILE IS
    Rev 1 of the ADR promised a removal gate and gated nothing
    (Implementation-log, T1b panel deferral D5). `_adr85_t2_design.md` Sec 2
    makes the hazard concrete: the 2D vertex adapter caches raw `Node*`
    `nts2dPrevFar`/`nts2dNextFar` pointers (LadrunoContactFE.cpp:255-256,
    armed at LadrunoContactHandler.cpp:1349-1401) that sit OUTSIDE the vertex
    pair's own connectivity `[slave, vertex]` -- they steer the
    ordered-ownership / wedge classification, never the force map.
    `pruneMissingNodes` (LadrunoContactSurface.cpp:47, driven from
    LadrunoContactDomain.cpp:1324 via `pruneRemovedNodes` at the TOP of
    `handle()`, LadrunoContactHandler.cpp:590-599) prunes SURFACE SEGMENTS
    whose own connectivity references a missing node; it does not know to
    also revalidate an adapter's cached far-node pointers, because those are
    not part of that adapter's connectivity by design. If any adapter
    persists ACROSS a `domainChanged()` (rather than being torn down and
    rebuilt from the freshly-pruned topology every time -- plausible, since a
    full rebuild would need to somehow resurrect committed friction/vertex-
    side-sign state for every UNAFFECTED pair), a `remove node` on a FAR node
    is a genuine dangling-pointer hazard: the NEXT evaluation of that
    surviving vertex pair dereferences a `Node*` whose underlying object the
    Domain has already deleted.

    That failure mode is a CRASH (SIGSEGV / access violation), not a
    catchable Python exception -- so, like the T0/T1b refusal batteries, the
    use-after-free test below runs in a CHILD PROCESS and the parent asserts
    on the child's EXIT CODE, distinguishing "the child was killed by a
    signal" from "the child finished (cleanly or with a named FATAL)". A
    test that merely checked "no Python exception was raised" would pass
    vacuously on a build where the child silently segfaults, because pytest
    would never even see that subprocess die -- the exit-code check is the
    part that makes this gate non-vacuous.

    The two tests below are written against the pruning/adapter-persistence
    behavior documented above. Both were run against the CURRENT worktree
    build (its `SRC/` has uncommitted T2-in-progress edits, confirmed via
    `git status`; the installed `dist/bin` binary is a few minutes newer than
    those edits) while writing this file -- see the deliverable report for
    exact pass/fail results and what each failure/pass means.

COVERAGE MAP -> ADR-85 G-T2(f) / D5
    ordinary survivor continuation      test_2d_removal_continues_on_survivors
    concave-vertex use-after-free       test_2d_removal_concave_vertex_survives

ON THE D5 CROSS-TYPE CAVEAT (side-flip warn latch) -- NOT GATED HERE
    `WARN_VTX2D_SIDE_FLIP` is ALREADY-SHIPPED T1b code
    (LadrunoContactFE.cpp:1021, LadrunoContactDomain.h:576): the condition is
    `nps==1 && committedSide!=0 && liveSide!=0 && liveSide!=committedSide`,
    where `liveSide` comes from `vertexEval2D`'s wedge test
    (`LadrunoContact2DKernel.h`: `aPrev = dot(xs-Xv, tPrev)/|tPrev|^2`,
    `aNext` symmetric, `inWedge = aPrev > tolIn && aNext < -tolIn`).

    I tried, empirically and against the live (WIP) build, to drive a
    reference-concave V-notch through a concave-to-convex deformation while
    keeping the condition's `liveSide != 0` half satisfied (i.e. keeping the
    slave inside the wedge test's region as the sign flips), across several
    rig variants:
      * the vertex node driven kinematically (DisplacementControl) with the
        slave held by a stiff static tether -- SINGULAR (a bare vertex pair
        with BOTH the slave's and the vertex's translational DOFs free is a
        rank-1 stiffness block over those two DOFs; nothing else in an
        element-free rig regularizes it under a static solver);
      * the vertex node force-driven (explicit, both nodes massed) with a
        stiff GROUND tether on the slave -- ran stably, but the tether kept
        the slave essentially fixed in place while the vertex sailed far
        away, so `r = slave - vertex` ends up dominated by the vertex's own
        large displacement and both `aPrev`/`aNext` flip sign TOGETHER (the
        point falls out of the wedge geometry entirely, not through the
        narrow same-in/opposite-sign band) -- no warning, ever, across a full
        0 -> 1.5 sweep;
      * the vertex node force-driven with NO tether at all (both nodes
        massed, the slave tracks the vertex purely through the kn coupling)
        -- `r` stays small (~1e-6) throughout, which is LIKELY inside
        `vertexEval2D`'s own degenerate-coincidence floor (a DIFFERENT named
        guard from the wedge-boundary one, per the ADR's numerical-thresholds
        note) and so plausibly never resolves a nonzero `liveSide` either;
      * a wider untethered gap (SEED=0.02) to escape that floor went
        numerically unstable (the two-mass, rank-1-stiffness system has no
        rigid-body restraint and diverged) before crossing the critical
        region.

    None of the four reproduced the warning. This is pre-existing T1b
    functionality, not T2 work, and diagnosing the exact geometric window
    (most likely a narrow band near the 180-degree turn, requiring the slave
    to be offset from the vertex by an amount comfortably above the
    coincidence floor but small enough to stay inside a still-evolving wedge
    -- almost certainly needs reading `vertexEval2D`'s literal tolerance
    constants rather than more trial and error) is better scoped as a short,
    dedicated follow-up than absorbed into this deliverable's budget. Filed
    honestly rather than shipped as a flaky or vacuously-weakened test.
"""
import os
import subprocess
import sys

import pytest

from _testbed import ops
from test_adr85_contact2d_t0_refusals import ENGINE_DIR

pytestmark = [pytest.mark.zone_a]

KN = 1.0e6
KS = 1.0e3


def _static():
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-12, 40, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")


# ---------------------------------------------------------------------------
# 1. ordinary segment removal: contact continues on the survivors
# ---------------------------------------------------------------------------
def test_2d_removal_continues_on_survivors():
    """G-T2(f), the base claim: `remove node` on the FAR end of one 2D master
    segment prunes that segment; a SEPARATE slave pair engaging the OTHER
    (surviving) segment of the same declared master surface is UNAFFECTED --
    same closed-form force before and after.

    THREE straight-line master nodes 101-102-103 (segments (101,102) and
    (102,103), sharing 102). Slave `A` engages segment (101,102) at its
    midpoint, held ALSO by a parallel y-anchor spring (the
    `test_pair_2d_now_live` idiom, tests/test_adr85_contact2d_t0_refusals.py)
    so its POST-removal state has a clean closed form too (spring alone,
    KS*|u|=P) rather than being merely "inactive". Slave `B` engages segment
    (102,103) and carries NO spring (it must stay penetrating and loaded
    throughout, so it needs no fallback restraint).

    `ops.remove("node", 101)` removes the FAR end of A's segment (not shared
    with B's segment, not the slave, not the vertex 102). The gate:

      * the run keeps solving (no FATAL, no crash) -- the ADR-78 removal-lane
        message ("removal lane" + the master surface tag) appears, exactly
        like the shipped 3D precedent
        (tests/test_contact_runtime_removal.py::
        test_removed_slave_node_is_pruned_and_the_run_continues);
      * A's pair is gone: its equilibrium reverts to the spring-only closed
        form (three orders of magnitude softer than the penalty-engaged
        state, same KN/KS=1e3 contrast test_pair_2d_now_live uses);
      * B's pair is UNTOUCHED: its force stays at the SAME closed form it had
        before the removal, to the same tolerance -- proving "continues on
        the survivors" is not merely "did not crash" but "kept transmitting
        the right force where nothing was removed".
    """
    P = 1.0e3
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(101, -2.0, 0.0)
    ops.node(102, 0.0, 0.0)
    ops.node(103, 2.0, 0.0)
    for t in (101, 102, 103):
        ops.fix(t, 1, 1)
    ops.contactSurface(10, "-master", 2, 101, 102, 102, 103)

    SEED = 1.0e-8
    ops.node(1, -1.0, -SEED)                 # slave A, on segment (101,102)
    ops.fix(1, 1, 0)                         # x fixed, y free (both slaves)
    ops.node(2, 1.0, -SEED)                  # slave B, on segment (102,103)
    ops.fix(2, 1, 0)
    ops.node(900, -1.0, -SEED)                # A's anchor ground
    ops.fix(900, 1, 1)
    ops.uniaxialMaterial("Elastic", 900, KS)
    ops.element("zeroLength", 900, 900, 1, "-mat", 900, "-dir", 2)
    ops.contactSurface(20, "-slave", 1)
    ops.contactSurface(30, "-slave", 2)
    ops.contact(1, 10, 20, KN, 0.0, 0.0, "-outward", 0.0, 1.0)
    ops.contact(2, 10, 30, KN, 0.0, 0.0, "-outward", 0.0, 1.0)

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, 0.0, -P)
    ops.load(2, 0.0, -P)
    _static()
    assert ops.analyze(1) == 0, "baseline (nothing removed) must solve"

    u_a_ref = (KN * SEED - P) / (KN + KS)     # A: penalty + spring in parallel
    # B carries ONLY the penalty: kn*(SEED - u) = P  =>  u = SEED - P/kn
    u_b_ref = SEED - P / KN
    u_a = ops.nodeDisp(1, 2)
    u_b = ops.nodeDisp(2, 2)
    assert abs(u_a - u_a_ref) / abs(u_a_ref) < 1.0e-6, f"A baseline {u_a!r} != {u_a_ref!r}"
    assert abs(u_b - u_b_ref) / abs(u_b_ref) < 1.0e-6, f"B baseline {u_b!r} != {u_b_ref!r}"

    ops.remove("node", 101)                   # A's FAR end -- prunes segment (101,102)
    # HOLD the load (LoadControl(0.0)) for the re-solve: LoadControl(1.0) would
    # advance the load factor CUMULATIVELY on a second analyze(1) call, doubling
    # the applied load rather than re-checking the same equilibrium.
    ops.integrator("LoadControl", 0.0)
    assert ops.analyze(1) == 0, "removal must prune and continue, not abort"

    u_a2 = ops.nodeDisp(1, 2)
    u_b2 = ops.nodeDisp(2, 2)
    u_a_spring_only = -P / KS
    assert abs(u_a2 - u_a_spring_only) / abs(u_a_spring_only) < 1.0e-6, (
        f"A's pair should be GONE after removing 101 (spring-only closed form "
        f"{u_a_spring_only!r}), got {u_a2!r}")
    assert abs(u_b2 - u_b_ref) / abs(u_b_ref) < 1.0e-6, (
        f"B's pair (untouched segment) drifted after A's segment was pruned: "
        f"{u_b2!r} vs baseline {u_b_ref!r} -- the survivor did not survive cleanly")


# ---------------------------------------------------------------------------
# 2. the concave-vertex use-after-free -- D5, the reason this file exists
# ---------------------------------------------------------------------------
# Child-process template: same Windows-safety shape as
# test_adr85_contact2d_t0_refusals.py (see that module's docstring for the
# PYTHONPATH / stdin=DEVNULL / add_dll_directory / LADRUNO_OPENSEES_QUIET
# rationale). ENGINE_DIR is imported from there so the two files cannot drift.
CHILD_UAF = r'''
import os, sys

_D = %(ENGINE_DIR)r
if os.path.isdir(_D):
    os.environ["PATH"] = _D + os.pathsep + os.environ.get("PATH", "")
    _add = getattr(os, "add_dll_directory", None)
    if _add is not None:
        try:
            _add(_D)
        except OSError:
            pass
    sys.path.insert(0, _D)

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops

# The T1a/T1b concave V-notch (test_adr85_contact2d_t1b_corner.py NOTCH):
# 101=(-1,0.5), 102=(0,0) the re-entrant vertex, 103=(1,0.5), symmetric about
# x=0 so the vertex normal is the exact bisector +y and the equilibrium is
# closed-form. 101 is the FAR node this test removes -- it is NOT the vertex
# (102) and NOT the slave, i.e. NOT in the vertex pair's own connectivity
# [slave, 102], exactly the D5 hazard shape.
KN, P = 1.0e6, 1.0e3
SEED = 1.0e-8
ops.model("basic", "-ndm", 2, "-ndf", 2)
ops.node(101, -1.0, 0.5)
ops.node(102, 0.0, 0.0)
ops.node(103, 1.0, 0.5)
for t in (101, 102, 103):
    ops.fix(t, 1, 1)
ops.contactSurface(10, "-master", 2, 101, 102, 102, 103)
ops.node(1, 0.0, -SEED)
ops.fix(1, 1, 0)                       # x fixed, y free -- the concave-corner rig
ops.contactSurface(20, "-slave", 1)
ops.contact(1, 10, 20, KN, 0.0, 0.0, "-outward", 0.0, 1.0)
ops.timeSeries("Linear", 1)
ops.pattern("Plain", 1, 1)
ops.load(1, 0.0, -P)
ops.constraints("LadrunoContact")
ops.numberer("Plain")
ops.system("FullGeneral")
ops.test("NormDispIncr", 1.0e-12, 40, 0)
ops.algorithm("Newton")
ops.integrator("LoadControl", 0.2)
ops.analysis("Static")

# Phase 1: 3 steps BEFORE the removal -- the vertex pair engages, commits, and
# (per the T1a/T1b contract) its side sign is captured at first capture. This
# is what gives a persisted adapter state to potentially go stale.
pens = []
for k in range(3):
    rc = ops.analyze(1)
    if rc != 0:
        print("UAF_PHASE1_ANALYZE_FAILED step=" + str(k) + " rc=" + str(rc))
        sys.exit(5)
    pens.append(SEED - ops.nodeDisp(1, 2))

# Phase 2: remove the FAR node (not the vertex, not the slave) and continue.
# If the vertex adapter's cached prevFar Node* is stale, THIS is where a
# genuine dangling-pointer dereference would fire -- on the NEXT evaluation
# of a pair that survives the removal, not on the removal call itself.
#
# NOTE: removing 101 also removes the VERTEX'S OWN adjacent segment
# (101,102) -- a concave vertex needs TWO adjacent segments, so node 102
# legitimately stops being a vertex and is reclassified as the open-terminal
# end of the one surviving segment (102,103), whose OWN outward normal
# (-0.4472, 0.8944 -- the T1b corner file's own NOTCH comment) is NOT the
# vertex bisector (0, 1). This is a real, physically-sensible topology change,
# not a hazard by itself; the per-step printed trace is what lets the parent
# tell "smoothly reclassified" apart from "silently reading garbage".
ops.remove("node", 101)
for k in range(3, 8):
    rc = ops.analyze(1)
    if rc != 0:
        print("UAF_PHASE2_ANALYZE_FAILED step=" + str(k) + " rc=" + str(rc))
        sys.exit(6)
    pens.append(SEED - ops.nodeDisp(1, 2))

print("UAF_SURVIVED pens=" + repr(pens))
sys.exit(0)
'''


def _run_uaf_child():
    script = CHILD_UAF % {"ENGINE_DIR": ENGINE_DIR}
    env = dict(os.environ)
    env["LADRUNO_OPENSEES_QUIET"] = "1"
    proc = subprocess.run(
        [sys.executable, "-c", script],
        stdin=subprocess.DEVNULL,
        capture_output=True, text=True, timeout=300,
        encoding="utf-8", errors="replace",
        env=env,
    )
    return proc, (proc.stdout or "") + "\n" + (proc.stderr or "")


def test_2d_removal_concave_vertex_survives():
    """G-T2(f) / D5: the concave-vertex use-after-free gate.

    Removing the FAR node of a chained-away segment must NOT crash the
    process on the NEXT evaluation of the surviving vertex pair, and the
    vertex pair must go on reporting the SAME physics it did before the
    removal -- the module docstring explains why "did not crash" alone would
    be a vacuous gate (a build that silently corrupted the pointer and
    happened not to fault on these particular steps would still pass a bare
    "exit 0" check) and why a genuine corruption is a process-killing signal,
    not a Python exception, hence the child process.

    THREE distinct failure signatures are told apart:
      exit killed by signal (no UAF_* marker at all) -> the dangling pointer
          was actually dereferenced and corrupted/faulted -- THE hazard;
      exit nonzero WITH a named FATAL/refusal, no "UAF_SURVIVED" -> the
          engine refused something (acceptable only if it is a DECLARED,
          named refusal -- surfaced in the assertion message for a human to
          judge; a silent nonzero exit with neither marker is NOT accepted);
      exit 0 WITH "UAF_SURVIVED" -> the gate's actual pass condition, checked
          against TWO different closed forms for the two phases (see below).

    WHY TWO DIFFERENT CLOSED FORMS, NOT ONE.  Removing 101 removes the
    VERTEX'S OWN adjacent segment (101,102) -- a concave vertex needs TWO
    adjacent segments to exist as a vertex at all, so node 102 correctly
    stops being a vertex mid-run and is reclassified as the OPEN-TERMINAL end
    of the one surviving segment (102,103).

    ADR-85 T4 (D4) UPDATE: node 102's reclassified pair is now the RADIAL
    end-cap (segment2DActive's nps==1 "exactly one far node" branch), not the
    old NTS2D_END_SLACK flat-line window this test was originally written
    against -- so the post-removal closed form changed, and changed BY
    DESIGN, not by regression (D4 retired NTS2D_END_SLACK entirely; see
    LadrunoContactFE.cpp and _adr85_t4_design.md).  The radial cap's normal
    is n = side*(x_slave - X_v)/||x_slave - X_v||, i.e. it points from the
    reclassified vertex (102, fixed at (0,0)) STRAIGHT AT the slave -- and
    this slave's x is FIXED at 0, IDENTICAL to X_v's x (`ops.fix(1, 1, 0)`:
    x fixed, y free; node 102 at x=0). So x_slave - X_v is EXACTLY (0, y) for
    every step, any y: the radial normal is (0, +-1), PURELY VERTICAL, not
    the surviving segment's own tilted (-0.4472, 0.8944) the old flat-line
    formula used. Post-removal stiffness is therefore kn*n_y^2 = kn*1.0^2 =
    kn (the SAME as phase 1's vertex bisector (0,1), by coincidence of this
    slave's on-axis placement -- not because the mechanism is unchanged: a
    slave placed off-axis from 102 would see a genuinely different, non-1.0
    n_y^2 under the radial cap). Measured against the T4 build (this
    session): the per-step trace is a clean, monotonic sequence matching
    n_y^2 = 1.0 from the FIRST post-removal step (no transient at all --
    unlike the old flat-line reclassification, the radial cap's C0-at-the-
    boundary property holds exactly here since x_slave == X_v.x already,
    with no window to cross).

    THE ASSERTIONS, from the full per-step trace (not just two endpoints):
      (a) phase 1 (3 steps, pre-removal): pen_k == k*P_step/KN, tight
          (1e-6), P_step = 0.2*P -- unchanged vertex mechanics;
      (b) phase 2: MONOTONIC, bounded growth (rules out NaN / sign flip /
          wild jump -- the coarse, mechanism-agnostic corruption check);
      (c) the LAST post-removal step's INCREMENT (steps 7-6) equals
          P_step/KN at 1e-4 relative -- the reclassified closed form under
          the T4 radial end-cap (n_y^2 = 1.0 for THIS on-axis slave, see
          above), confirming the surviving pair is not merely stable but
          STILL COMPUTING THE PHYSICALLY CORRECT geometry, not a corrupted-but-smooth
          reading.
    """
    proc, out = _run_uaf_child()
    assert proc.returncode == 0, (
        "the concave-vertex adapter did not survive removing a far node "
        f"cleanly (exit {proc.returncode} -- a negative/signal-derived "
        "returncode on POSIX, or a Windows exception code, means the "
        "process was KILLED, i.e. the use-after-free actually fired).\n"
        f"{out}"
    )
    assert "UAF_SURVIVED" in out, (
        f"the child exited 0 but never reached the survival marker -- "
        f"inspect the output for which phase failed.\n{out}")

    pens = None
    for line in out.splitlines():
        if line.startswith("UAF_SURVIVED"):
            pens = eval(line[len("UAF_SURVIVED pens="):], {"__builtins__": {}})
    assert pens is not None and len(pens) == 8, f"unexpected pens trace:\n{out}"

    P, KN_ = 1.0e3, 1.0e6
    P_STEP = 0.2 * P

    for k in range(3):                          # phase 1: unchanged vertex mechanics
        ref = (k + 1) * P_STEP / KN_
        assert abs(pens[k] - ref) / ref < 1.0e-6, (
            f"pre-removal step {k}: pen={pens[k]!r} != {ref!r}\ntrace={pens}")

    for k in range(1, 8):                        # monotonic, bounded throughout
        assert pens[k] > pens[k - 1], (
            f"penetration is not monotonically increasing at step {k}: "
            f"{pens}")
        assert pens[k] < 0.02, (
            f"penetration exploded at step {k} (tunnelling / corrupted read): "
            f"{pens}")

    # ADR-85 T4 (D4): the reclassified closed form under the radial end-cap,
    # n_y^2 = 1.0 for this on-axis slave (x_slave == X_v.x == 0 always -- see
    # the module docstring), NOT the old flat-line n_y^2 = 0.8. Checked deep
    # into phase 2 (steps 6 -> 7); with the radial cap there is no transient
    # step to skip past (C0 holds exactly here), but the later station is
    # kept to stay comparable with phase 1's own steady-state check above.
    ref_incr = P_STEP / (1.0 * KN_)
    got_incr = pens[7] - pens[6]
    assert abs(got_incr - ref_incr) / ref_incr < 1.0e-4, (
        f"post-removal steady-state increment {got_incr!r} != the "
        f"reclassified open-terminal closed form {ref_incr!r} "
        f"(P_step/(n_y^2*kn), n_y^2=1.0 for the T4 radial end-cap on this "
        f"on-axis slave) -- the pair survived without crashing but is "
        f"computing the wrong geometry.\ntrace={pens}")




# ---------------------------------------------------------------------------
# 3. D5 cross-type caveat: reference-concave corner deformed convex -- NOT
#    GATED HERE. See the module docstring's "on the cross-type caveat" note.
# ---------------------------------------------------------------------------
