"""ADR-76 R4 acceptance test 4 -- checker.

Reads adr76_factoronce.h5 (written by adr76_factoronce_model.tcl) and asserts
that `algorithm ModifiedNewton -initial -factoronce` honours BOTH options, in
either argument order.

    caseA  -initial               formTangent calls == NSTEPS   (no latch)
    caseB  -initial -factoronce   formTangent calls == 1        (latched)
    caseC  -factoronce -initial   identical to caseB            (order-independent)
    caseD  -factoronce            formTangent calls == 1, but the CURRENT tangent
                                  is latched, so it converges in far fewer
                                  iterations than caseB/caseC

Pre-fix (upstream `if (OPS_GetNumRemainingInputArgs() > 0)`) the parser read one
option and stopped, so caseB collapsed onto caseA and caseC onto caseD. Each
assertion below therefore fails loudly on a regression to the old parser.

Converged displacements must agree across all cases. NOTE the scope of that claim:
it holds here because this deck uses a FORCE-residual test (NormUnbalance). Under
NormDispIncr / RelativeNormDispIncr / EnergyIncr the acceptance criterion is
dU = K_f^-1 R, so a different K_f accepts at a different point and the "K_f only
sets the path" reassurance does NOT hold. See ADR-76 s3.

Scope limit of this smoke, stated so it is not over-read: the profiled window is
LINEAR (all ten bars sit in the same Steel01 hardening branch under monotonic
load), so Ki and Kt are each constant through it. That is what makes the parser
assertions crisp, but it also means this test validates the OPTION PARSER and
nothing else -- it cannot detect any of the tangent-staleness hazards that
-factoronce carries (domain change under a latch, re-arm on convergence failure).

    python3.12 adr76_factoronce_check.py <dir-containing-the-h5>
"""

from __future__ import annotations

import os
import sys

sys.path.insert(
    0,
    os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "..", "..", "Ladruno_tools", "profiler_viewer",
    ),
)

from profiler_results import ProfilerResults  # noqa: E402

NSTEPS = 5
CASES = ("caseA", "caseB", "caseC", "caseD", "caseE", "caseF")


def _find(node: dict, name: str):
    """Depth-first search for the first rollup node with the given name."""
    if node.get("name") == name:
        return node
    for ch in node.get("children", []):
        hit = _find(ch, name)
        if hit is not None:
            return hit
    return None


def _calls(pr: ProfilerResults, run: str, phase: str) -> int:
    node = _find(pr.rollup(run), phase)
    if node is None:
        raise AssertionError(f"{run}: no '{phase}' node in the rollup tree")
    return int(node["calls"])


def main() -> int:
    out = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(out, "adr76_factoronce.h5")
    problems = 0

    def check(ok: bool, msg: str) -> None:
        nonlocal problems
        print(("  PASS  " if ok else "  FAIL  ") + msg)
        if not ok:
            problems += 1

    with ProfilerResults(path) as pr:
        have = set(pr.runs())
        missing = set(CASES) - have
        if missing:
            print(f"FAIL: missing runs {sorted(missing)} (have {sorted(have)})")
            return 1

        tan = {c: _calls(pr, c, "formTangent") for c in CASES}
        # convTest fires once per iteration in ModifiedNewton, so it is the
        # iteration counter.
        its = {c: _calls(pr, c, "convTest") for c in CASES}

        print("\nformTangent calls:", tan)
        print("iterations       :", its, "\n")

        # -- `-factoronce` takes effect when it is NOT the first option --------
        check(tan["caseA"] == NSTEPS,
              f"caseA (-initial): formTangent == {NSTEPS} per analysis "
              f"(one per step, no latch) [got {tan['caseA']}]")
        check(tan["caseB"] == 1,
              f"caseB (-initial -factoronce): formTangent == 1 -- `-factoronce` "
              f"survived a preceding `-initial` [got {tan['caseB']}]")

        # -- `-initial` takes effect when it is NOT the first option ----------
        check(tan["caseC"] == 1,
              f"caseC (-factoronce -initial): formTangent == 1 [got {tan['caseC']}]")
        check(its["caseC"] == its["caseB"],
              f"caseC == caseB iteration-for-iteration -- option order is "
              f"irrelevant [got {its['caseC']} vs {its['caseB']}]")
        # -- the spelling variants the adversarial review caught -----------------
        for c, spelling in (("caseE", "-initial -factorOnce"),
                            ("caseF", "-Initial -factorOnce")):
            check(tan[c] == tan["caseB"] and its[c] == its["caseB"],
                  f"{c} (`{spelling}`) behaves exactly like caseB -- the camelCase "
                  f"/ capitalised spellings accepted by Linear, ExpressNewton and "
                  f"NewtonRaphson are honoured here too "
                  f"[got {tan[c]}/{its[c]} vs {tan['caseB']}/{its['caseB']}]")

        # Asserted, not merely printed: the review found both of these were only
        # claimed in prose. caseA vs caseB is the strongest discriminator available
        # (same iteration count, 5 assemblies vs 1).
        check(tan["caseD"] == 1,
              f"caseD (-factoronce): formTangent == 1 [got {tan['caseD']}]")
        check(its["caseA"] == its["caseB"],
              f"caseA and caseB run the IDENTICAL iteration count on 5 assemblies "
              f"vs 1 [got {its['caseA']} vs {its['caseB']}]")
        check(its["caseD"] == NSTEPS,
              f"caseD converges in exactly one iteration per step (exact tangent on "
              f"a linear window) [got {its['caseD']}, expected {NSTEPS}]")

        check(its["caseD"] < its["caseB"],
              f"caseD (-factoronce alone) latches the CURRENT tangent and so "
              f"converges in fewer iterations than caseB/caseC, which latch the "
              f"INITIAL one -- proves `-initial` is not silently dropped "
              f"[got {its['caseD']} vs {its['caseB']}]")

    # -- the equilibrium itself must not move ---------------------------------
    # K_f sets the path, not the answer. All cases converge to the same
    # tolerance, so agree to that tolerance rather than bit-for-bit.
    disp_path = os.path.join(out, "adr76_factoronce_disp.txt")
    disp = {}
    with open(disp_path) as fh:
        for line in fh:
            run, val = line.split()
            disp[run] = float(val)
    print("tip ux           :", disp, "\n")
    ref = disp["caseA"]
    spread = max(abs(v - ref) for v in disp.values())
    check(spread <= 1.0e-9 * max(abs(ref), 1.0e-30),
          f"all cases converge to the same tip displacement "
          f"[max deviation {spread:.3e} about {ref:.12e}]")

    # Absolute anchor (2026-07-26): mutual agreement alone would pass a
    # regression that drives every case to the same WRONG answer. The chain's
    # converged tip is analytically 7.0: the profiled window advances the
    # load from the lambda=150 pre-yield to lambda = 150 + 5*2.0 = 160, and
    # in a series chain every bar carries that full axial force. Steel01 with
    # Fy=100, E=1000, b=0.1, A=1: post-yield eps = 0.1 + (sigma-100)/(b*E) =
    # 0.1 + 60/100 = 0.7 per unit-length bar, x10 bars = 7.0. caseD lands on
    # it exactly; the -initial cases converge to within their 1e-8 force
    # tolerance, measured 9.4e-10 away. 1e-7 accepts that spread and still
    # fails on any answer that is wrong rather than merely iterated.
    ANCHOR = 7.0
    worst = max(abs(v - ANCHOR) for v in disp.values())
    check(worst <= 1.0e-7,
          f"every case lands on the analytic tip ux == {ANCHOR} "
          f"[max |ux - {ANCHOR}| = {worst:.3e}]")

    print()
    if problems:
        print(f"ADR-76 R4 smoke: {problems} FAILURE(S)")
        return 1
    print("ADR-76 R4 smoke: all checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
