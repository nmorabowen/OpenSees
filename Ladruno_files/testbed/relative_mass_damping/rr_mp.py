"""T6: the drag toy R_ff deck split over 2 and 4 OpenSeesMP ranks against the serial R_ff run.

Usage: py -3.12 rr_mp.py <OpenSeesMP.exe> <mpiexec> <serial R_ff folder>
Partitions (flat per rank, shared nodes carry their mass on one rank only, as the STKO decks):
  np 2: rank 0 soil column + base; rank 1 structure (node 43 shared, its mass on rank 0)
  np 4: rank 0 soil 1-22 + base; rank 1 soil 22-43; rank 2 node 43, foundation 200, node 201, spring 1001;
        rank 3 nodes 201 (mass), 202, 203, springs 1002-1003 (201 shared, its mass on rank 3)
"""
import re
import shutil
import subprocess
import sys
import pathlib

import numpy as np

EXE, MPIEXEC, SER = sys.argv[1], sys.argv[2], pathlib.Path(sys.argv[3])
MODE = sys.argv[4] if len(sys.argv) > 4 else "R"   # R: reference command; A: compensating loads (oracle)
import tempfile
WORK = pathlib.Path(tempfile.gettempdir()) / "relative_mass_damping" / ("rr_mp_work" + ("" if MODE == "R" else "_A"))
CMD = "relativeMassDamping"

src = (SER / "main.tcl").read_text().splitlines()
nodes, soil_eles = {}, {}
for l in src:
    m = re.match(r"node (\d+) (\S+)(?: -mass (\S+))?$", l)
    if m:
        nodes[int(m.group(1))] = (m.group(2), m.group(3))
    m = re.match(r"element zeroLength (\d+) (\d+) (\d+) -mat 1 -dir 1$", l)
    if m:
        soil_eles[int(m.group(1))] = (int(m.group(2)), int(m.group(3)))
mats = [l for l in src if l.startswith("uniaxialMaterial")]
lys = next(l for l in src if l.startswith("element zeroLength 99 "))
springs = {int(re.match(r"element zeroLength (\d+)", l).group(1)): l for l in src if l.startswith("element zeroLength 100")}
reg2 = next(l for l in src if l.startswith("region 2 "))
reg1 = next(l for l in src if l.startswith("region 1 "))
beta = reg2.split("-rayleigh")[1]
alpha = reg1.split("-rayleigh")[1]
base = [l for l in src if l.startswith("timeSeries Path 11") or l.startswith("pattern Plain 11")]
ser13 = next(l for l in src if l.startswith("timeSeries Path 13"))


def node_line(t, with_mass):
    x, m = nodes[t]
    return f"node {t} {x}" + (f" -mass {m}" if (m and with_mass) else "")


def rank_block(np_, r):
    L = []
    if np_ == 2:
        soil = range(1, 44) if r == 0 else []
        soil_mass = set(soil)
        se = [e for e in soil_eles if r == 0]
        struct = [] if r == 0 else [43, 200, 201, 202, 203]
        smass = {200, 201, 202, 203} if r == 1 else set()
        spr = [1001, 1002, 1003] if r == 1 else []
        eq = r == 1
        lysb = r == 0
    else:
        soil = {0: range(1, 23), 1: range(22, 44)}.get(r, [])
        soil_mass = set(range(1, 23)) if r == 0 else set(range(23, 44)) if r == 1 else set()
        se = [e for e, (i, j) in soil_eles.items() if (r == 0 and j <= 22) or (r == 1 and i >= 22)]
        struct = {2: [43, 200, 201], 3: [201, 202, 203]}.get(r, [])
        smass = {2: {200}, 3: {201, 202, 203}}.get(r, set())
        spr = {2: [1001], 3: [1002, 1003]}.get(r, [])
        eq = r == 2
        lysb = r == 0
    for t in soil:
        L.append(node_line(t, t in soil_mass))
    for t in struct:
        if t in soil:
            continue
        L.append(node_line(t, (t in smass) or (t == 43 and t in soil_mass)))
    if lysb:
        L += [node_line(100, False), "fix 100 1"]
    L += mats
    L += [f"element zeroLength {e} {soil_eles[e][0]} {soil_eles[e][1]} -mat 1 -dir 1" for e in se]
    if lysb:
        L.append(lys)
        L += base
    if eq:
        L.append("equalDOF 43 200 1")
    L += [springs[e] for e in spr]
    if spr:
        L.append(f"region 2 -ele {' '.join(map(str, spr))} -rayleigh{beta}")
    if smass:
        L.append(f"region 1 -nodeOnly {' '.join(map(str, sorted(smass)))} -rayleigh{alpha}")
        if MODE == "R":
            L += [ser13, f"set rr [{CMD} 1 -vel 1 13 -region 1]", 'puts "RR [getPID] $rr"']
        else:
            a = float(alpha.split()[0])
            L += [ser13, "pattern Plain 13 13 {"] + [f"  load {t} {a * float(nodes[t][1])!r}" for t in sorted(smass)] + ["}"]
    rec = [t for t in [43, 200, 201, 202, 203] if t in (set(soil) | set(struct)) and (t in smass or t in soil_mass or (t == 43 and t in soil_mass))]
    if rec:
        for q in ("disp", "vel"):
            L.append(f"recorder Node -file {q}.{r}.txt -time -precision 17 -node {' '.join(map(str, rec))} -dof 1 {q}")
        L.append(f"set recnodes {{{' '.join(map(str, rec))}}}")
    return L


for np_ in (2, 4):
    d = WORK / f"np{np_}"
    if d.exists():
        shutil.rmtree(d)
    d.mkdir(parents=True)
    for f in ("vrec.txt", "vff.txt"):
        shutil.copy(SER / f, d / f)
    T = ["wipe", "model basic -ndm 1 -ndf 1"]
    for r in range(np_):
        T.append(f"{'if' if r == 0 else 'elseif'} {{[getPID] == {r}}} {{" if r == 0 else f"}} elseif {{[getPID] == {r}}} {{")
        T += ["  " + l for l in rank_block(np_, r)]
    T.append("}")
    T += ["constraints Penalty 1e14 1e14", "numberer ParallelPlain", "system Mumps", "test NormDispIncr 1e-12 10",
          "algorithm Linear", "integrator Newmark 0.5 0.25", "analysis Transient",
          "set ok [analyze 28000 0.0025]", 'puts "done [getPID] $ok [getTime]"', "wipe"]
    (d / "main.tcl").write_text("\n".join(T) + "\n")
    r = subprocess.run([MPIEXEC, "-localonly", "-n", str(np_), EXE, "main.tcl"], cwd=d, capture_output=True, text=True, timeout=3600)
    txt = r.stdout + r.stderr
    (d / "log.txt").write_text(txt)
    print(f"np {np_}:", [l for l in txt.splitlines() if l.startswith(("RR", "done"))])
    # compare with serial
    sd = np.loadtxt(SER / "disp.txt"); sv = np.loadtxt(SER / "vel.txt")
    scol = {t: i + 1 for i, t in enumerate([43, 200, 201, 202, 203])}
    got = {}
    for rr in range(np_):
        f = d / f"disp.{rr}.txt"
        if not f.exists():
            continue
        tcl = (d / "main.tcl").read_text()
        blk = tcl.split(f"[getPID] == {rr}}}")[1].split("set recnodes {")[1].split("}")[0]
        tags = list(map(int, blk.split()))
        a = np.loadtxt(f); b = np.loadtxt(d / f"vel.{rr}.txt")
        for i, t in enumerate(tags):
            got[t] = (a[:, i + 1], b[:, i + 1])
    for t in sorted(got):
        ed = np.abs(got[t][0] - sd[:, scol[t]]).max() / np.abs(sd[:, scol[t]]).max()
        ev = np.abs(got[t][1] - sv[:, scol[t]]).max() / np.abs(sv[:, scol[t]]).max()
        print(f"  np {np_} node {t}: disp {ed:.2e} vel {ev:.2e} (relative to serial peak)")
    if 203 in got and 200 in got:
        roof = got[203][0] - got[200][0]; roof_s = sd[:, scol[203]] - sd[:, scol[200]]
        print(f"  np {np_} roof drift vs serial {np.abs(roof - roof_s).max() / np.abs(roof_s).max():.2e}")
