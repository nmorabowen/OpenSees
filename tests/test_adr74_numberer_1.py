"""ADR-74 N0 — the G1 numbering oracle, proven BEFORE any numberer code exists.

The `ladrunoNumbering <file>` verb (SRC/tcl/commands.cpp, ADR-74 N0) dumps each
rank's view of the assigned equation numbering: one sorted line per node,
`<nodeTag> <ndf> <id0> <id1> ...`. This file is the instrument the later
ADR-74 gates lean on (G1 bit-identity, G1b bijection, G1c fallback), so N0's
job is to prove the instrument itself.

DISCOVERED EN ROUTE (2026-07-22, the oracle's first real catch): under
`system MPIDiagonal`, the globally-consistent numbering the ParallelNumberer
computes exists only TRANSIENTLY. `MPIDiagonalSOE::setSize` builds its
shared-DOF exchange map from the global ids, then deliberately REWRITES every
DOF group's ids to rank-local 0..n_local-1 (`MPIDiagonalSOE.cpp:408-429,
"renumber DOFs 0 through size"`). A post-analysis dump therefore shows
rank-LOCAL ids: shared boundary nodes legitimately DISAGREE across ranks, and
per-rank ids are dense. This is by design, not a defect — verified by tracing
`numberDOF`'s own assignment (correct globals) against the final dump (globals
compacted per rank, exact -offset match). Consequently the oracle is TWO decks:

  Deck A — `system Mumps` (keeps global ids; no setID in Mumps*SOE):
      the STRICT numberer oracle.
      A1 shared boundary nodes carry IDENTICAL ids in every owning rank's view
      A2 the deduped union of ids is exactly {0..neqn-1}, each once
      A3 two identical runs -> byte-identical dumps (determinism)
  Deck B — `system MPIDiagonal` (the production explicit stack):
      the END-STATE oracle, post-localization.
      B1 each rank's free ids are a dense per-rank bijection 0..n_r-1
      B2 the fully-fixed node reads -1 on every dof
      B3 two identical runs -> byte-identical dumps
  MUT — planted faults make every comparator FAIL (an oracle that cannot
      fail is not an oracle); no binary needed.

For G1 (N2/N3): Deck A byte-compare `LadrunoParallelRCM` vs stock
`ParallelRCM` is the strict numberer-identity gate; Deck B byte-compare is the
production end-state gate. Both decks read the numberer verb from
$::env(ADR74_NUMBERER).

Model: a 2D truss chain, SPMD-partitioned by rank arithmetic (rank r owns
elements r*EPR+1..(r+1)*EPR; segment-boundary nodes are shared).

Theory: Ladruno_implementation/74_ladruno_parallel_numberer_adr.md (§Gates,
§Implementation plan N0).
"""
import os
import shutil
import subprocess
from pathlib import Path

import pytest

pytestmark = [pytest.mark.zone_a]

REPO = Path(__file__).resolve().parents[1]
EPR = 5  # elements per rank

# ---------------------------------------------------------------------------
# binary discovery — env override first, then the packaged dist tree
# ---------------------------------------------------------------------------


def _find_mp():
    exe = os.environ.get("LADRUNO_MP_EXE")
    mpi = os.environ.get("LADRUNO_MPIEXEC")
    if not exe:
        cand = REPO / "dist" / "bin" / ("OpenSeesMP.exe" if os.name == "nt" else "OpenSeesMP")
        exe = str(cand) if cand.exists() else None
    if not mpi:
        cand = REPO / "dist" / "openseesmp" / ("mpiexec.exe" if os.name == "nt" else "mpiexec")
        mpi = str(cand) if cand.exists() else shutil.which("mpiexec")
    return exe, mpi


MP_EXE, MPIEXEC = _find_mp()
needs_mp = pytest.mark.skipif(
    MP_EXE is None or MPIEXEC is None,
    reason="OpenSeesMP / mpiexec not found (set LADRUNO_MP_EXE / LADRUNO_MPIEXEC)",
)

# %(analysis)s is the solver-stack block; %(stride)d spreads node tags for the
# G1c sparse-tag-space case (stride 1 = the dense default). A gravity-style
# lateral load makes the 1-step solution a physics probe for G1b same-physics.
DECK = r"""
set pid [getPID]
set np  [getNP]
set EPR %(epr)d
set STR %(stride)d

model BasicBuilder -ndm 2 -ndf 2
uniaxialMaterial Elastic 1 1000.0

set n0 [expr $pid*$EPR]
set n1 [expr ($pid+1)*$EPR]
for {set i $n0} {$i <= $n1} {incr i} {
    node [expr $i*$STR+1] [expr double($i)] 0.0
    mass [expr $i*$STR+1] 1.0 1.0
}
if {$pid == 0} { fix 1 1 1 }
for {set e [expr $n0+1]} {$e <= $n1} {incr e} {
    element truss $e [expr ($e-1)*$STR+1] [expr $e*$STR+1] 1.0 1
}
pattern Plain 1 Linear {
    load [expr $n1*$STR+1] 1.0 0.0
}

constraints Transformation
if {[info exists ::env(ADR74_NUMBERER)]} {
    eval "numberer $::env(ADR74_NUMBERER)"
} else {
    numberer ParallelRCM
}
%(analysis)s

ladrunoNumbering numbering_[getPID].txt
set fd [open disp_[getPID].txt w]
puts $fd [format "%%.17g" [nodeDisp [expr $n1*$STR+1] 1]]
close $fd
"""

# Deck A: MUMPS keeps the numberer's GLOBAL ids (no setID in Mumps*SOE) —
# implicit Newmark, 1 step. The strict oracle.
ANALYSIS_MUMPS = r"""
system Mumps
test NormDispIncr 1e-8 20
algorithm Linear
integrator Newmark 0.5 0.25
analysis Transient
analyze 1 1.0e-3
"""

# Deck B: the production explicit stack. MPIDiagonalSOE::setSize LOCALIZES
# ids after numbering (see module docstring) — end-state oracle.
ANALYSIS_MPIDIAG = r"""
system MPIDiagonal
test NormDispIncr 1e-8 20
algorithm Linear
integrator CentralDifferenceLadruno
analysis Transient
analyze 1 1.0e-6
"""


def _run_mp(workdir: Path, np: int, analysis: str, numberer: str = "ParallelRCM",
            stride: int = 1) -> dict:
    """Run the chain deck at np ranks; return {rank: {nodeTag: (ids...)}}."""
    workdir.mkdir(parents=True, exist_ok=True)
    deck = workdir / "chain.tcl"
    deck.write_text(DECK % {"epr": EPR, "analysis": analysis, "stride": stride},
                    encoding="ascii")

    env = dict(os.environ)
    env["ADR74_NUMBERER"] = numberer
    env["LADRUNO_OPENSEES_QUIET"] = "1"
    # the fork ships its own MPI runtime next to mpiexec (impi.dll etc.)
    env["PATH"] = str(Path(MPIEXEC).parent) + os.pathsep + env.get("PATH", "")
    # a build-tree exe cannot find init.tcl on its own (the same gap the cluster
    # wrapper openseesmp.sh plugs). NB a failed Tcl init still exits rc=0 on
    # this stack — the dump-count assertion below is the real failure signal.
    if "TCL_LIBRARY" not in env:
        tcl_lib = REPO / "dist" / "lib" / "tcl8.6"
        if tcl_lib.exists():
            env["TCL_LIBRARY"] = str(tcl_lib)

    proc = subprocess.run(
        [MPIEXEC, "-n", str(np), MP_EXE, deck.name],
        stdin=subprocess.DEVNULL,   # load-bearing on Windows -- see quirks ledger
        cwd=workdir, env=env, capture_output=True, text=True, timeout=300,
    )
    (workdir / "run_output.txt").write_text(proc.stdout + proc.stderr,
                                            encoding="utf-8", errors="replace")
    dumps = sorted(workdir.glob("numbering_*.txt"))
    assert proc.returncode == 0, (
        f"OpenSeesMP np{np} failed rc={proc.returncode}\n"
        f"stdout:\n{proc.stdout[-2000:]}\nstderr:\n{proc.stderr[-2000:]}"
    )
    assert len(dumps) == np, (
        f"expected {np} dumps, got {[d.name for d in dumps]}\n"
        f"stdout:\n{proc.stdout[-2000:]}\nstderr:\n{proc.stderr[-2000:]}"
    )
    return {int(d.stem.split("_")[1]): parse_dump(d) for d in dumps}


# ---------------------------------------------------------------------------
# the comparators — module-level on purpose: G1 (N2/N3) imports these
# ---------------------------------------------------------------------------


def parse_dump(path: Path) -> dict:
    """`<nodeTag> <ndf> <id0> ...` lines -> {nodeTag: (ids...)}."""
    out = {}
    for line in path.read_text().splitlines():
        parts = line.split()
        if not parts:
            continue
        tag, ndf, ids = int(parts[0]), int(parts[1]), tuple(int(x) for x in parts[2:])
        assert len(ids) == ndf, f"{path.name}: node {tag} ndf={ndf} but {len(ids)} ids"
        assert tag not in out, f"{path.name}: duplicate node {tag}"
        out[tag] = ids
    return out


def assert_shared_consistent(dumps: dict) -> dict:
    """GLOBAL-id oracle only (Deck A): shared nodes identical across ranks.
    Invalid after MPIDiagonalSOE localization — see module docstring."""
    merged = {}
    for rank, view in sorted(dumps.items()):
        for tag, ids in view.items():
            if tag in merged:
                assert merged[tag] == ids, (
                    f"node {tag}: rank views disagree {merged[tag]} vs {ids} (rank {rank})"
                )
            merged[tag] = ids
    return merged


def assert_bijection(merged: dict, fixed_nodes=(1,)):
    """GLOBAL-id oracle only (Deck A): free ids = {0..neqn-1} exactly once."""
    free = []
    for tag, ids in merged.items():
        if tag in fixed_nodes:
            assert all(i == -1 for i in ids), f"fixed node {tag} ids {ids} != all -1"
            continue
        for i in ids:
            assert i >= 0, f"node {tag}: unnumbered/constrained marker {i} survived setup"
            free.append(i)
    assert sorted(free) == list(range(len(free))), (
        f"equation ids are not a bijection onto 0..{len(free)-1}"
    )
    return len(free)


def assert_rank_local_bijection(dumps: dict, fixed_nodes=(1,)):
    """LOCALIZED oracle (Deck B): each rank's free ids dense 0..n_r-1, once."""
    total = 0
    for rank, view in sorted(dumps.items()):
        free = []
        for tag, ids in view.items():
            if tag in fixed_nodes:
                assert all(i == -1 for i in ids), (
                    f"rank {rank}: fixed node {tag} ids {ids} != all -1"
                )
                continue
            for i in ids:
                assert i >= 0, f"rank {rank}: node {tag} carries marker {i}"
                free.append(i)
        assert sorted(free) == list(range(len(free))), (
            f"rank {rank}: local ids not dense 0..{len(free)-1}"
        )
        total += len(free)
    return total


def assert_dumps_identical(a: dict, b: dict):
    """The G1 bit-compare: same ranks, same nodes, same ids, exactly."""
    assert sorted(a) == sorted(b), f"rank sets differ: {sorted(a)} vs {sorted(b)}"
    for rank in a:
        assert a[rank] == b[rank], f"rank {rank}: numbering differs"


def _expected_neqn(np: int) -> int:
    n_nodes = np * EPR + 1
    return 2 * (n_nodes - 1)  # node 1 fully fixed, ndf=2 elsewhere


# ---------------------------------------------------------------------------
# Deck A — strict global oracle (system Mumps)
# ---------------------------------------------------------------------------


@needs_mp
@pytest.mark.parametrize("np", [2, 8])
def test_global_oracle_mumps(tmp_path, np):
    dumps = _run_mp(tmp_path, np, ANALYSIS_MUMPS)
    merged = assert_shared_consistent(dumps)
    assert len(merged) == np * EPR + 1
    assert assert_bijection(merged) == _expected_neqn(np)


@needs_mp
def test_global_oracle_determinism(tmp_path):
    a = _run_mp(tmp_path / "a", 2, ANALYSIS_MUMPS)
    b = _run_mp(tmp_path / "b", 2, ANALYSIS_MUMPS)
    assert_dumps_identical(a, b)
    for r in (0, 1):
        fa = (tmp_path / "a" / f"numbering_{r}.txt").read_bytes()
        fb = (tmp_path / "b" / f"numbering_{r}.txt").read_bytes()
        assert fa == fb, f"rank {r}: dump bytes differ between identical runs"


# ---------------------------------------------------------------------------
# Deck B — end-state oracle on the production stack (system MPIDiagonal)
# ---------------------------------------------------------------------------


@needs_mp
@pytest.mark.parametrize("np", [2, 8])
def test_localized_oracle_mpidiagonal(tmp_path, np):
    dumps = _run_mp(tmp_path, np, ANALYSIS_MPIDIAG)
    # per-rank dense localization (MPIDiagonalSOE.cpp:408-429)
    total = assert_rank_local_bijection(dumps)
    # shared nodes are double-counted across ranks in the localized view:
    # total >= global neqn, with equality only at np==1
    assert total >= _expected_neqn(np)
    assert sum(len(v) for v in dumps.values()) == np * (EPR + 1)


@needs_mp
def test_localized_oracle_determinism(tmp_path):
    a = _run_mp(tmp_path / "a", 2, ANALYSIS_MPIDIAG)
    b = _run_mp(tmp_path / "b", 2, ANALYSIS_MPIDIAG)
    assert_dumps_identical(a, b)


# ---------------------------------------------------------------------------
# N1 gate — LadrunoParallelNumberer delegate mode vs stock, bit-identical.
# Simultaneously a null-test of the oracle and of the subclass wiring
# (promotion, tag-forwarding ctors, verb registration). ADR-74 §plan N1.
# ---------------------------------------------------------------------------


@needs_mp
@pytest.mark.parametrize("np", [2, 8])
@pytest.mark.parametrize("deck,stock,ladruno", [
    (ANALYSIS_MUMPS, "ParallelRCM", "LadrunoParallelRCM"),      # strict oracle
    (ANALYSIS_MPIDIAG, "ParallelRCM", "LadrunoParallelRCM"),    # production end-state
], ids=["mumps-rcm", "mpidiag-rcm"])
def test_g1_rcm_identity(tmp_path, np, deck, stock, ladruno):
    """G1: the RCM verb is bit-identical to stock — the load-bearing gate.
    (Delegate-proven at N1, must SURVIVE the T0 engine swap.)"""
    a = _run_mp(tmp_path / "stock", np, deck, numberer=stock)
    b = _run_mp(tmp_path / "ladruno", np, deck, numberer=ladruno)
    assert_dumps_identical(a, b)


def _read_disp(workdir: Path, np: int) -> dict:
    return {r: float((workdir / f"disp_{r}.txt").read_text().strip())
            for r in range(np)}


@needs_mp
@pytest.mark.parametrize("np", [2, 8])
def test_g1b_plain_bijection_and_physics(tmp_path, np):
    """G1b: LadrunoParallelPlain fixes the upstream tag-0 quirk, so it is NOT
    bit-identical to stock ParallelPlain — the gate is a valid numbering
    (strict global oracle on the Mumps deck) + same physics as stock."""
    a = _run_mp(tmp_path / "stock", np, ANALYSIS_MUMPS, numberer="ParallelPlain")
    b = _run_mp(tmp_path / "ladruno", np, ANALYSIS_MUMPS, numberer="LadrunoParallelPlain")
    # both are valid global numberings...
    for dumps in (a, b):
        merged = assert_shared_consistent(dumps)
        assert assert_bijection(merged) == _expected_neqn(np)
    # ...and the physics agrees to solver precision
    da, db = _read_disp(tmp_path / "stock", np), _read_disp(tmp_path / "ladruno", np)
    for r in range(np):
        assert da[r] == pytest.approx(db[r], rel=1e-10, abs=1e-14), (
            f"rank {r}: plain-verb physics diverged: {da[r]} vs {db[r]}"
        )


@needs_mp
@pytest.mark.parametrize("np", [2, 8])
def test_g1c_sparse_tag_fallback(tmp_path, np):
    """G1c: a strided tag space (maxTag >> numVertex) crosses the kappa guard
    into the hash fallback — numbering must STILL be bit-identical to stock,
    and the fallback must OBSERVABLY engage (review: a silently-dense run
    would pass without testing the hash path — an oracle that cannot fail)."""
    stride = 1_000_000  # maxTag ~ 1e7 for tens of nodes => dense budget blown
    a = _run_mp(tmp_path / "stock", np, ANALYSIS_MUMPS,
                numberer="ParallelRCM", stride=stride)
    b = _run_mp(tmp_path / "ladruno", np, ANALYSIS_MUMPS,
                numberer="LadrunoParallelRCM", stride=stride)
    assert_dumps_identical(a, b)
    merged = assert_shared_consistent(b)
    assert assert_bijection(merged) == _expected_neqn(np)
    out = (tmp_path / "ladruno" / "run_output.txt").read_text(errors="replace")
    assert "RefIndex -> hash fallback" in out, (
        "kappa guard never engaged — the hash path went untested"
    )


# Plain-handler stack: equalDOF MPs survive to the numberer as -4 markers,
# which the -4 fixup (the rewritten "second quadratic") must resolve.
MP_TIES = r"""
# two MPs on ONE constrained node (different retained nodes, different dofs):
# exercises the -4 fixup multi-MP application-order preservation
equalDOF 2 4 1
equalDOF 3 4 2
"""


@needs_mp
def test_g1_mp_fixup_identity(tmp_path):
    """Review MAJOR: the rewritten -4 fixup had zero gate coverage. Under
    `constraints Plain`, equalDOF leaves -4 markers the numberer resolves
    itself. Stock-vs-Ladruno must stay bit-identical. (Bijection does NOT
    hold here: constrained dofs share eqn ids with their retained dofs.)

    np1 ON PURPOSE: at np>=2 this deck HANGS under the STOCK numberer too —
    the `constraints Plain` + equalDOF + parallel-MUMPS stack deadlocks in
    the solve (observed: both ranks spin in MPI busy-poll; reproduced with
    `numberer ParallelRCM`, i.e. upstream, not ADR-74 code). The -4 fixup
    itself is rank-local logic (build the constrainedNode->MPs index, resolve
    this rank's -4 markers), which np1 exercises completely. The upstream
    MP+Plain+Mumps hang is ledger material, not an ADR-74 gate."""
    analysis = "constraints Plain" + MP_TIES + ANALYSIS_MUMPS
    a = _run_mp(tmp_path / "stock", 1, analysis, numberer="ParallelRCM")
    b = _run_mp(tmp_path / "ladruno", 1, analysis, numberer="LadrunoParallelRCM")
    assert_dumps_identical(a, b)
    assert_shared_consistent(b)
    # the ties must actually be present and resolved: constrained node 4's
    # dof-1 id == node 2's dof-1 id, dof-2 id == node 3's dof-2 id
    view = b[0]
    assert view[4][0] == view[2][0] and view[4][1] == view[3][1], (
        f"equalDOF ties unresolved: n2={view[2]} n3={view[3]} n4={view[4]}"
    )


def _run_mp_nofix(workdir, np, analysis, numberer):
    """The tag-0 micro-case needs DOF-group tag 0 to be a FREE node — with the
    standard deck the fixed node 1 IS vertex 0 (zero dofs), so its ordering
    position is invisible and stock/Ladruno coincide. Strip the fix; use the
    explicit stack (no factorization => a floating chain is fine for 1 step)."""
    workdir.mkdir(parents=True, exist_ok=True)
    deck = workdir / "chain.tcl"
    body = DECK % {"epr": EPR, "analysis": analysis, "stride": 1}
    body = body.replace("if {$pid == 0} { fix 1 1 1 }", "")
    deck.write_text(body, encoding="ascii")
    env = dict(os.environ)
    env["ADR74_NUMBERER"] = numberer
    env["LADRUNO_OPENSEES_QUIET"] = "1"
    env["PATH"] = str(Path(MPIEXEC).parent) + os.pathsep + env.get("PATH", "")
    tcl_lib = REPO / "dist" / "lib" / "tcl8.6"
    if tcl_lib.exists():
        env.setdefault("TCL_LIBRARY", str(tcl_lib))
    proc = subprocess.run([MPIEXEC, "-n", str(np), MP_EXE, deck.name],
                          stdin=subprocess.DEVNULL,
                          cwd=workdir, env=env, capture_output=True, text=True,
                          timeout=300)
    dumps = sorted(workdir.glob("numbering_*.txt"))
    assert len(dumps) == np, f"expected {np} dumps\n{proc.stdout[-1500:]}"
    return {int(d.stem.split("_")[1]): parse_dump(d) for d in dumps}


@needs_mp
def test_g1b_tag0_divergence_documented(tmp_path):
    """The tag-0 micro-case (review): stock ParallelPlain's zero-filled
    theOrderedRefs treats vertex tag 0 as already-ordered and numbers it
    LAST; LadrunoParallelPlain orders it naturally. With vertex 0 FREE the
    dumps must DIFFER (the fix is real and visible) while both remain valid
    per-rank numberings (MPIDiagonal end-state: rank-local dense)."""
    a = _run_mp_nofix(tmp_path / "stock", 2, ANALYSIS_MPIDIAG, "ParallelPlain")
    b = _run_mp_nofix(tmp_path / "ladruno", 2, ANALYSIS_MPIDIAG,
                      "LadrunoParallelPlain")
    for dumps in (a, b):
        assert_rank_local_bijection(dumps, fixed_nodes=())
    with pytest.raises(AssertionError):
        assert_dumps_identical(a, b)   # identical would mean the fix vanished


@needs_mp
def test_lagrange_negative_ref_fails_loudly(tmp_path):
    """Review MAJOR-1 companion: a Lagrange handler under MP produces
    node-less DOF groups (ref=-1) that stock silently FUSES across ranks.
    The Ladruno verb must fail LOUDLY (FATAL + abort), never silently
    misnumber, and never hang past the harness timeout."""
    analysis = "constraints Lagrange" + ANALYSIS_MUMPS
    workdir = tmp_path
    workdir.mkdir(parents=True, exist_ok=True)
    deck = workdir / "chain.tcl"
    body = DECK % {"epr": EPR, "analysis": analysis, "stride": 1}
    body = body.replace("constraints Transformation", "")
    deck.write_text(body, encoding="ascii")

    env = dict(os.environ)
    env["ADR74_NUMBERER"] = "LadrunoParallelRCM"
    env["LADRUNO_OPENSEES_QUIET"] = "1"
    env["PATH"] = str(Path(MPIEXEC).parent) + os.pathsep + env.get("PATH", "")
    tcl_lib = REPO / "dist" / "lib" / "tcl8.6"
    if tcl_lib.exists():
        env.setdefault("TCL_LIBRARY", str(tcl_lib))

    # NB: on Windows Intel MPI, P0's exit(-1) does NOT tear down fellow ranks —
    # the non-root spins in MPI_Recv (empirically: orphaned OpenSeesMP at 100%
    # CPU). This is exactly the review-MAJOR failure mode the FATAL exists to
    # surface; the harness must reap the process TREE itself.
    proc = subprocess.Popen([MPIEXEC, "-n", "2", MP_EXE, deck.name],
                            cwd=workdir, env=env,
                            stdin=subprocess.DEVNULL,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            text=True)
    try:
        out, _ = proc.communicate(timeout=90)
    except subprocess.TimeoutExpired:
        out = ""
    finally:
        if os.name == "nt":
            subprocess.run(["taskkill", "/F", "/T", "/PID", str(proc.pid)],
                           stdin=subprocess.DEVNULL, capture_output=True)
        else:
            proc.kill()
        try:
            rest, _ = proc.communicate(timeout=15)
            out = (out or "") + (rest or "")
        except Exception:
            pass
    assert "FATAL LadrunoParallelNumberer" in (out or ""), (
        f"expected the loud ref<0 abort; got tail:\n{(out or '')[-1500:]}"
    )


# ---------------------------------------------------------------------------
# MUT — the comparators must be able to fail (no binary needed)
# ---------------------------------------------------------------------------


def test_oracle_mutation_must_fail():
    base = {0: {1: (-1, -1), 2: (0, 1), 3: (2, 3)},
            1: {3: (2, 3), 4: (4, 5)}}
    merged = assert_shared_consistent(base)
    assert_bijection(merged)

    # (a) duplicate id -> global bijection fails
    mutated = {0: {1: (-1, -1), 2: (0, 1), 3: (2, 3)},
               1: {3: (2, 3), 4: (4, 4)}}
    with pytest.raises(AssertionError):
        assert_bijection(assert_shared_consistent(mutated))

    # (b) shared-node disagreement -> Deck-A consistency fails
    split = {0: {3: (2, 3)}, 1: {3: (2, 99)}}
    with pytest.raises(AssertionError):
        assert_shared_consistent(split)

    # (c) identity comparator fails on any difference
    with pytest.raises(AssertionError):
        assert_dumps_identical(base, mutated)

    # (d) unnumbered placeholder rejected by both oracles
    unnumbered = {0: {1: (-1, -1), 2: (-2, 0)}}
    with pytest.raises(AssertionError):
        assert_bijection(assert_shared_consistent(unnumbered))
    with pytest.raises(AssertionError):
        assert_rank_local_bijection(unnumbered)

    # (e) non-dense local ids -> Deck-B bijection fails
    gappy = {0: {1: (-1, -1), 2: (0, 2)}}
    with pytest.raises(AssertionError):
        assert_rank_local_bijection(gappy)

    # (f) and the localized oracle ACCEPTS what Deck A must reject:
    # overlapping per-rank dense ids (the MPIDiagonal end state)
    localized = {0: {1: (-1, -1), 2: (0, 1)}, 1: {2: (0, 1), 3: (2, 3)}}
    assert_rank_local_bijection(localized)
    with pytest.raises(AssertionError):
        assert_bijection(assert_shared_consistent({0: {2: (0, 1)}, 1: {3: (0, 1)}}))
