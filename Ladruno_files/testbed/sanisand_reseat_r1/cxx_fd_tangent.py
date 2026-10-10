"""Tangent consistency of the WP-151 C++ (SAS-ME, IntScheme 129) where R1 acts.

For tiny increments d*u from a loaded state, the Richardson finite difference
R = 2 (s(d)-s0)/d - (s(2d)-s0)/(2d) must equal T u, with T the replay's tangent_ep:
the continuum tangent at the LOADED state, plastic loading assumed (what the first
SAS-ME stage integrates to first order).  Only plastic replays (no elastic update or
stage, rc 0) are compared.

  wall   the 5 wall states x 32 directions, prototypes 11 (DM04), 13 (R1), 14 (floor), 16 (f+h)
  ring   the 80 ring rows x 32 directions, prototype 16 (floor + hysteresis)
  cap    states where the CAP binds at the start: along the prototype-13 fan trials whose
         stages hit the cap, lambda = 1/40 .. 1 of the trial increment, the first state whose
         tiny replays bind the cap
  kink   the worst cap case with d from 1e-9 down to 1e-11

    py -3.12 -S cxx_fd_tangent.py <worktree> [<bin_dir>] [wall ring cap kink]
    (PYTHONPATH = the 3.12 site-packages; the default bin_dir is <worktree>/dist/bin.
     Snapshot dist/bin first if a build may run meanwhile: a loaded pyd blocks the build.)
Measured on 7b66acf64/8a7884fa2 (win32): wall median 2e-6 (DM04 1.6e-7); ring median 8e-7;
cap 18 states, median 9e-6, one 9e-4 at the cap's switch-on that falls to 1e-8 as d shrinks."""
import math
import os
import sys

W = os.path.abspath(sys.argv[1])
ARGS = sys.argv[2:]
BIN = os.path.abspath(ARGS.pop(0)) if ARGS and os.path.isdir(ARGS[0]) else os.path.join(W, "dist", "bin")
PARTS = ARGS or ["wall", "ring", "cap", "kink"]
os.add_dll_directory(BIN)
sys.path[:0] = [BIN, os.path.join(W, "tests"), os.path.join(W, "Ladruno_scripts")]
import opensees as ops  # noqa: E402
assert os.path.normcase(os.path.dirname(os.path.abspath(ops.__file__))) == os.path.normcase(BIN)
import sanisand_replay as sr  # noqa: E402
import wp151_reseat_tools as T  # noqa: E402

D0 = 1.0e-9


def nrm(v):
    return math.sqrt(sum(x * x for x in v))


def mv(M, u):
    return [sum(M[i][j] * u[j] for j in range(6)) for i in range(6)]


def fd(tag, st, u, D=D0):
    """(rel error, sas stats of the d replay) or None if not a plastic replay."""
    o1 = T.replay(ops, tag, st, [D * x for x in u])
    o2 = T.replay(ops, tag, st, [2 * D * x for x in u])
    if int(o1["rc"]) or int(o2["rc"]):
        return None
    s1 = o1["sas"] or {}
    if s1.get("elastic", 0) or s1.get("elasticStages", 0):
        return None
    s0 = st["sigma"]
    fd1 = [(a - b) / D for a, b in zip(o1["sigma"], s0)]
    fd2 = [(a - b) / (2 * D) for a, b in zip(o2["sigma"], s0)]
    rich = [2 * a - b for a, b in zip(fd1, fd2)]
    Tu = mv(o1["tangent_ep"], u)
    return nrm([a - b for a, b in zip(rich, Tu)]) / max(nrm(Tu), 1e-30), s1


def report(name, errs):
    errs = sorted(errs, key=lambda e: -e[0])
    if not errs:
        print(f"{name}: no plastic replays")
        return
    print(f"{name}: {len(errs)} plastic replays; rel |FD - T u| median {errs[len(errs) // 2][0]:.2e}, "
          f"max {errs[0][0]:.2e} ({errs[0][1]})")


ops.wipe()
T.define(ops, tags=[11, 13, 14, 16])
dirs = T.fib_sphere(32)

if "wall" in PARTS:
    for tag, name in ((11, "wall DM04"), (13, "wall R1 (f+h+cap)"), (14, "wall floor alone"),
                      (16, "wall floor+hyst")):
        errs = []
        for st in T.refusers():
            for u3 in dirs:
                r = fd(tag, st, T.deps_of(u3, 1.0))
                if r:
                    errs.append((r[0], f"{st['leg']}/{st['k']} floored {r[1].get('hFloored', 0):g}"))
        report(name, errs)

if "ring" in PARTS:
    errs = []
    for path in sr.RING_CSVS:
        for row in sr.read_ring_csv(path):
            for u3 in dirs:
                r = fd(16, row, T.deps_of(u3, 1.0))
                if r:
                    errs.append((r[0], f"{os.path.basename(path)} {row['element']}/{row['gp']}"))
    report("ring floor+hyst", errs)

worst = None
if "cap" in PARTS or "kink" in PARTS:
    errs = []
    for kind, key, st, de in T.jobs():
        if kind != "fan":
            continue
        o = T.replay(ops, 13, st, de)
        if int(o["rc"]) or not (o["sas"] or {}).get("hSoftCapped", 0):
            continue
        u = [x / nrm(de) for x in de]
        for k in range(1, 41):
            om = T.replay(ops, 13, st, [k / 40.0 * x for x in de])
            if int(om["rc"]):
                continue
            mid = dict(sigma=list(om["sigma"]), alpha=list(om["alpha"]),
                       alpha_in=list(om["alpha_in"]), z=list(om["z"]), e=om["e"])
            r = fd(13, mid, u)
            if r and r[1].get("hSoftCapped", 0):
                errs.append((r[0], f"{key} lambda={k / 40.0:g}", mid, u))
                break
    report("cap bound at the start", errs)
    if errs:
        worst = max(errs, key=lambda e: e[0])

if "kink" in PARTS and worst:
    print(f"kink: {worst[1]}")
    for D in (1e-9, 3e-10, 1e-10, 3e-11, 1e-11):
        r = fd(13, worst[2], worst[3], D)
        if r:
            print(f"   d {D:.0e}: rel {r[0]:.2e}, cap bound {r[1].get('hSoftCapped', 0):g}")
