#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# =============================================================================
# ADR-68 P0.5 — LysmerTriangle + EnergyBalanceRecorder sign-pinning experiment
#
# PURPOSE
#   The fork's EnergyBalanceRecorder integrates IE = int F_resisting . v dt
#   from Element::getResistingForce() at commit. LysmerTriangle (stage 0)
#   returns its stale `internalForces` member there, which is filled by
#   getResistingForceIncInertia() with C . v_gnd (the compliant-base
#   injection tractions from a LysmerVelocityLoader elemental load, NOT the
#   physical damper force C . v_node — note the `0*v(i) + gnd_velocity(i)`
#   in LysmerTriangle.cpp:500-502). Suspected accounting leak:
#     (a) leak ~= -W_inject : RES ~= 0 (accidentally closes), IE strongly NEGATIVE
#     (b) leak ~= +W_inject : RES ~= -2 W_inject, IE positive
#
# DISCOVERED CONSTRAINTS (grep evidence, see RESULTS.md)
#   * LysmerVelocityLoader has NO interpreter hook: no Tcl or Python eleLoad
#     type creates it anywhere in SRC. It is C++-only dead code from the
#     interpreter's point of view. Run B below proves this at runtime.
#   * LysmerTriangle itself is Tcl-only (registered in TclElementCommands.cpp
#     but absent from the Python OpenSeesElementCommands.cpp map), so the
#     whole experiment drives the INSTALLED OpenSees.exe (Tcl) from Python.
#   * The identical accounting pathway (elemental-load work rides
#     getResistingForce into IE with resisting-force sign, while ULW only sees
#     Node::getUnbalancedLoad) IS reachable via stdBrick body forces
#     (eleLoad -type -BrickW). Run C uses that as a mechanism-identical
#     surrogate to pin the sign empirically.
#
# RUNS
#   A1 control-initvel : Lysmer base, initial vx on top nodes, no loads.
#                        Radiation-only. Lysmer contributes 0 to IE (no
#                        loader => gnd_velocity = 0 => stale force = 0).
#   A2 control-pulse   : Lysmer base, Ricker NODAL force on top nodes.
#                        Nodal loads are counted in ULW => closure test.
#   B  loader-attempt  : try to invoke the velocity loader from Tcl; capture
#                        the parser rejection verbatim.
#   C  surrogate       : Ricker ELEMENTAL drive (BrickW body force, x-dir)
#                        + Lysmer base. W_inject integrated independently
#                        from recorded nodal velocities. Pins the sign.
#
# Uses the INSTALLED build only — nothing is compiled here.
# =============================================================================

import math
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
EXE = r"C:\Program Files\Ladruno\OpenSees\bin\OpenSees.exe"

# ---- model constants --------------------------------------------------------
RHO = 2000.0
VS = 200.0
NU = 0.25
G = RHO * VS * VS                       # 8.0e7
E = 2.0 * G * (1.0 + NU)                # 2.0e8
VP = VS * math.sqrt((2.0 - 2.0 * NU) / (1.0 - 2.0 * NU))   # ~346.41
NL = 10                                  # brick layers, 1x1x10 m column
DT = 1.0e-3
NSTEP = 600
F0 = 10.0                                # Ricker center frequency [Hz]
T0 = 0.15                                # Ricker delay [s]


def ricker(t):
    a = math.pi * F0 * (t - T0)
    return (1.0 - 2.0 * a * a) * math.exp(-a * a)


def ricker_path_values(n=701):
    return [ricker(i * DT) for i in range(n)]


# ---- Tcl deck builder -------------------------------------------------------

def common_model(with_body_force):
    """Nodes 4k+c (k=0..NL, c=1..4), bricks 1..NL, Lysmer 101/102 on base."""
    bf = " 1.0 0.0 0.0" if with_body_force else ""
    L = []
    L.append("model basic -ndm 3 -ndf 3")
    L.append(f"nDMaterial ElasticIsotropic 1 {E:.8e} {NU} {RHO}")
    L.append(f"# soil column 1x1x{NL} m, Vs={VS}, Vp={VP:.4f}")
    for k in range(NL + 1):
        b = 4 * k
        L.append(f"node {b+1} 0.0 0.0 {float(k)}")
        L.append(f"node {b+2} 1.0 0.0 {float(k)}")
        L.append(f"node {b+3} 1.0 1.0 {float(k)}")
        L.append(f"node {b+4} 0.0 1.0 {float(k)}")
    # 1-D x-polarized shear: fix y,z everywhere, x free (base held by Lysmer)
    for n in range(1, 4 * (NL + 1) + 1):
        L.append(f"fix {n} 0 1 1")
    for k in range(NL):
        b, t = 4 * k, 4 * (k + 1)
        L.append(
            f"element stdBrick {k+1} {b+1} {b+2} {b+3} {b+4} "
            f"{t+1} {t+2} {t+3} {t+4} 1{bf}"
        )
    # base face z=0 covered by 2 LysmerTriangle, default stage 0 (pure damping)
    L.append(f"element LysmerTriangle 101 1 2 3 {RHO} {VP:.8f} {VS}")
    L.append(f"element LysmerTriangle 102 1 3 4 {RHO} {VP:.8f} {VS}")
    return L


def analysis_block(tag):
    return [
        f'recorder EnergyBalance -file "{tag}_energy.txt" -time',
        f'recorder Node -file "{tag}_vels.txt" -time '
        f"-nodeRange 1 {4*(NL+1)} -dof 1 vel",
        "constraints Plain",
        "numberer RCM",
        "if {[catch {system UmfPack} m]} { "
        'if {[catch {system ProfileSPD} m2]} { system FullGeneral } }',
        "integrator Newmark 0.5 0.25",
        "algorithm Linear",
        "analysis Transient",
        f"set ok [analyze {NSTEP} {DT}]",
        'puts "ANALYZE_RETURN $ok"',
        "wipe",
        'puts "DONE"',
    ]


def path_series(tsTag):
    vals = " ".join(f"{v:.10e}" for v in ricker_path_values())
    return [f"timeSeries Path {tsTag} -dt {DT} -values {{{vals}}}"]


def deck_A1():
    L = common_model(False)
    for n in range(4 * NL + 1, 4 * NL + 5):           # top layer nodes 41..44
        L.append(f"setNodeVel {n} 1 0.1 -commit")
    L += analysis_block("A1")
    return L


def deck_A2():
    L = common_model(False)
    L += path_series(1)
    L.append("pattern Plain 1 1 {")
    for n in range(4 * NL + 1, 4 * NL + 5):
        L.append(f"    load {n} 250.0 0.0 0.0")        # 1000 N total, Ricker-scaled
    L.append("}")
    L += analysis_block("A2")
    return L


def deck_B():
    L = common_model(False)
    L += path_series(1)
    L.append("pattern Plain 1 1 {")
    L.append('    if {[catch {eleLoad -ele 101 102 -type -lysmerVelocityLoader 1} msg]} '
             '{ puts "LOADER_TRY1_FAIL: $msg" } else { puts "LOADER_TRY1_OK" }')
    L.append('    if {[catch {eleLoad -ele 101 102 -type -LysmerVelocityLoader 1} msg]} '
             '{ puts "LOADER_TRY2_FAIL: $msg" } else { puts "LOADER_TRY2_OK" }')
    L.append('    if {[catch {eleLoad -ele 101 102 -type -lysmer 1} msg]} '
             '{ puts "LOADER_TRY3_FAIL: $msg" } else { puts "LOADER_TRY3_OK" }')
    L.append("}")
    L.append('puts "DONE"')
    L.append("wipe")
    return L


def deck_B2():
    """Decisive loader test: drive ONLY via the velocity-loader eleLoad.

    The Tcl eleLoad handler silently returns TCL_OK for unknown -type
    strings (TclModelBuilder.cpp tail: 'if get here we have successfully
    created the load'), so run B's OK is not proof. If the loader were
    actually created, C.v_gnd would shake the column; if it is a silent
    no-op the response stays identically zero.
    """
    L = common_model(False)
    L += path_series(1)
    L.append("pattern Plain 1 1 {")
    L.append("    eleLoad -ele 101 102 -type -lysmerVelocityLoader 1")
    L.append("    eleLoad -ele 101 102 -type -fooBarBazNotALoad 1")
    L.append("}")
    L += analysis_block("B2")
    return L


def deck_C():
    L = common_model(True)                             # bricks carry b=(1,0,0)
    L += path_series(1)
    ele_list = " ".join(str(i) for i in range(1, NL + 1))
    L.append("pattern Plain 1 1 {")
    L.append(f"    eleLoad -ele {ele_list} -type -BrickW")
    L.append("}")
    L += analysis_block("C")
    return L


# ---- run + parse ------------------------------------------------------------

def run_deck(tag, lines):
    tcl = os.path.join(HERE, f"{tag}.tcl")
    with open(tcl, "w") as f:
        f.write("\n".join(lines) + "\nexit\n")
    p = subprocess.run([EXE, tcl], cwd=HERE, capture_output=True, text=True,
                       timeout=300)
    out = (p.stdout or "") + (p.stderr or "")
    with open(os.path.join(HERE, f"{tag}_console.log"), "w") as f:
        f.write(out)
    return out


def parse_energy(tag):
    rows = []
    with open(os.path.join(HERE, f"{tag}_energy.txt")) as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 7:
                rows.append([float(x) for x in parts[:7]])
    return rows  # t KE IE DW ULW RES ERR%


def summarize(tag, rows):
    def col(i):
        return [r[i] for r in rows]
    names = ["KE", "IE", "DW", "ULW", "RES", "ERR%"]
    last = rows[-1]
    s = [f"[{tag}] rows={len(rows)}  t_end={last[0]:.3f}"]
    for i, nm in enumerate(names, start=1):
        c = col(i)
        s.append(f"  {nm:4s} final={last[i]: .6e}  min={min(c): .6e}  "
                 f"max={max(c): .6e}")
    return "\n".join(s)


def w_inject_C():
    """Independent W_inject = sum_nodes int f_i(t) vx_i(t) dt for run C.

    f_i(t) = loadFactor(t) * b_x * V_trib(i); regular trilinear hex =>
    V/8 per node per brick => 0.125 m^3 (single-layer nodes) / 0.25 (interior).
    """
    times, vels = [], []
    with open(os.path.join(HERE, "C_vels.txt")) as f:
        for line in f:
            parts = [float(x) for x in line.split()]
            if len(parts) == 4 * (NL + 1) + 1:
                times.append(parts[0])
                vels.append(parts[1:])
    trib = []
    for k in range(NL + 1):
        tv = 0.125 if (k == 0 or k == NL) else 0.25
        trib += [tv] * 4
    W = 0.0
    prev_p = None
    for it, t in enumerate(times):
        fac = ricker(t)
        power = sum(fac * 1.0 * trib[i] * vels[it][i] for i in range(len(trib)))
        if prev_p is None:
            W += power * (times[0])       # rectangle seed like the recorder
        else:
            W += 0.5 * (prev_p + power) * (t - times[it - 1])
        prev_p = power
    return W


def main():
    print("=" * 78)
    print("ADR-68 P0.5 sign-pinning | installed exe:", EXE)
    print("python:", sys.executable)
    if not os.path.isfile(EXE):
        print("FATAL: installed OpenSees.exe not found"); sys.exit(2)

    results = {}

    for tag, deck in (("A1", deck_A1()), ("A2", deck_A2())):
        print("-" * 78)
        out = run_deck(tag, deck)
        ok = re.search(r"ANALYZE_RETURN (-?\d+)", out)
        print(f"[{tag}] analyze return:", ok.group(1) if ok else "MISSING", flush=True)
        rows = parse_energy(tag)
        results[tag] = rows
        print(summarize(tag, rows))

    print("-" * 78)
    outB = run_deck("B", deck_B())
    for m in re.findall(r"LOADER_TRY\d+_(?:FAIL|OK).*", outB):
        print("[B]", m.strip())
    print("[B] console excerpt:",
          " | ".join(l.strip() for l in outB.splitlines()
                     if "eleLoad" in l or "WARNING" in l)[:400])

    print("-" * 78)
    outB2 = run_deck("B2", deck_B2())
    ok = re.search(r"ANALYZE_RETURN (-?\d+)", outB2)
    print("[B2] analyze return:", ok.group(1) if ok else "MISSING")
    vmax = 0.0
    with open(os.path.join(HERE, "B2_vels.txt")) as f:
        for line in f:
            parts = [float(x) for x in line.split()]
            if len(parts) > 1:
                vmax = max(vmax, max(abs(x) for x in parts[1:]))
    print(f"[B2] max |vx| over all nodes/steps under loader-only drive: {vmax:.3e}")
    print("[B2] =>", "loader ENGAGED (response nonzero)" if vmax > 0.0
          else "loader is a SILENT NO-OP (zero response; eleLoad accepted the "
               "unknown type without creating any load)")

    print("-" * 78)
    outC = run_deck("C", deck_C())
    ok = re.search(r"ANALYZE_RETURN (-?\d+)", outC)
    print("[C] analyze return:", ok.group(1) if ok else "MISSING")
    rowsC = parse_energy("C")
    results["C"] = rowsC
    print(summarize("C", rowsC))

    W = w_inject_C()
    last = rowsC[-1]
    KEe, IEe, DWe, ULWe, RESe = last[1], last[2], last[3], last[4], last[5]
    print(f"\n[C] independent W_inject (sum int f.v dt) = {W: .6e}")
    print(f"[C] IE_end/W = {IEe/W: .4f}   RES_end/W = {RESe/W: .4f}   "
          f"DW_end/W = {DWe/W: .4f}   ULW_end/W = {ULWe/W: .4f}")

    print("=" * 78)
    if IEe < -0.5 * abs(W):
        print("VERDICT: (a) -- injection work leaks into IE with NEGATIVE sign")
        print("(leak = -W_inject; IE_end/W above). ULW is blind to the elemental")
        print("load. NOTE: RES here reads ~ +0.5*W (not the textbook 0) because a")
        print("SEPARATE defect superposes: Lysmer stage-0 absorption under implicit")
        print("Newmark books only ~half the absorbed energy in DW (see control A2,")
        print("RES = +0.50*ULW). With energy-consistent absorption (DW -> W) the")
        print("-W_inject leak alone would drive RES -> 0, i.e. hypothesis (a):")
        print("balance accidentally closes while IE trends strongly negative.")
    elif IEe > 0.5 * abs(W) and RESe < -1.5 * abs(W):
        print("VERDICT: (b) -- leak ~ +W_inject; RES ~ -2 W_inject.")
    else:
        print("VERDICT: inconclusive by thresholds -- inspect numbers above.")
    print("=" * 78)


if __name__ == "__main__":
    main()
