"""relativeMassDamping: mass-proportional Rayleigh damping relative to a prescribed
reference velocity, f_D = alphaM * M * (v - v_ref(t)).

Classic-Tcl gates (OpenSees.exe; the SSI decks that use the feature are classic
Tcl / OpenSeesMP), plus one openseespy registration smoke.

  T1  SDOF m u'' + alpha m (u' - v_ref) + k u = p(t), Newmark average
      acceleration: matches the closed recursion, and equals the same model with
      the reference written as a nodal load +alpha m v_ref (the oracle).
  T2  `rayleighForces` (reaction flag 2) and `nodalRayleighForces` return the
      relative damping force.
  T3  EnergyBalance -v2 -relativeMassSplit: DW_rel >= 0 and non-decreasing,
      DW_rel + W_ref = DW_nod, RES identical to the load-oracle run.
  T4  A zero reference series is bitwise identical to no command (the guarded
      branch reproduces the original arithmetic).
  T5  Integrators that form the nodal damping force through
      DOF_Group::addD_Force (HHT_TP) are refused with a FATAL line.
  T6  Command audit return and refusals.
"""
import os
import sys
import math
import subprocess
from pathlib import Path

import pytest

pytestmark = [pytest.mark.zone_a]

REPO = Path(__file__).resolve().parents[1]
CMD = "relativeMassDamping"


def _find_exe():
    exe = os.environ.get("LADRUNO_TCL_EXE")
    if exe:
        return exe if Path(exe).exists() else None
    cand = REPO / "dist" / "bin" / ("OpenSees.exe" if os.name == "nt" else "OpenSees")
    return str(cand) if cand.exists() else None


TCL_EXE = _find_exe()
needs_exe = pytest.mark.skipif(TCL_EXE is None, reason="classic OpenSees exe not found (build it, or set LADRUNO_TCL_EXE)")

M, K, ALPHA, DT, N = 2.0, 800.0, 0.7, 0.005, 2000


def vref(t):
    return 0.3 * math.sin(2 * math.pi * 1.3 * t) + 0.1 * math.sin(2 * math.pi * 7.0 * t)


def pload(t):  # triangular pulse 0 -> 50 -> 0 over 0.4 s, as the Path series below
    return 50.0 * t / 0.2 if t <= 0.2 else (50.0 * (0.4 - t) / 0.2 if t <= 0.4 else 0.0)


def run(tmp, name, body, integ="Newmark 0.5 0.25", algo="Newton"):
    d = tmp / name
    d.mkdir()
    with open(d / "vr.txt", "w") as fh:
        fh.write("\n".join(repr(vref(i * DT)) for i in range(N + 1)) + "\n")
    with open(d / "vz.txt", "w") as fh:
        fh.write("\n".join("0.0" for _ in range(N + 1)) + "\n")
    deck = [
        "wipe", "model basic -ndm 1 -ndf 1", "node 1 0.0", f"node 2 0.0 -mass {M!r}", "fix 1 1",
        f"uniaxialMaterial Elastic 1 {K!r}", "element zeroLength 1 1 2 -mat 1 -dir 1",
        f"region 1 -nodeOnly 2 -rayleigh {ALPHA!r} 0.0 0.0 0.0",
        f"timeSeries Path 5 -dt {DT!r} -filePath vr.txt",
        f"timeSeries Path 7 -dt {DT!r} -filePath vz.txt",
        "timeSeries Path 6 -time {0.0 0.2 0.4 1000.0} -values {0.0 50.0 0.0 0.0}",
        "pattern Plain 6 6 { load 2 1.0 }",
        "recorder Node -file d.txt -time -precision 17 -node 2 -dof 1 disp",
        "recorder Node -file v.txt -time -precision 17 -node 2 -dof 1 vel",
        "recorder Node -file rf.txt -time -precision 17 -node 2 -dof 1 rayleighForces",
        "recorder Node -file nrf.txt -time -precision 17 -node 2 -dof 1 nodalRayleighForces",
        "recorder EnergyBalance -time -precision 17 -v2 -relativeMassSplit -file e2.txt",
    ] + body + [
        "constraints Plain", "numberer Plain", "system BandGeneral", "test NormDispIncr 1e-14 20",
        f"algorithm {algo}", f"integrator {integ}", "analysis Transient",
        f"set ok [analyze {N} {DT!r}]", 'puts "done $ok"', "wipe",
    ]
    (d / "main.tcl").write_text("\n".join(deck) + "\n")
    p = subprocess.run([TCL_EXE, "main.tcl"], cwd=str(d), env=dict(os.environ, LADRUNO_OPENSEES_QUIET="1"),
                       stdin=subprocess.DEVNULL, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                       encoding="utf-8", errors="replace", timeout=300)
    return d, p.stdout, p.returncode


def col(d, f, c=1):
    return [float(l.split()[c]) for l in (d / f).read_text().split("\n") if l.strip()]


def rows(d, f):
    return [[float(x) for x in l.split()] for l in (d / f).read_text().split("\n") if l.strip()]


REF = [f'puts "RR [{CMD} 1 -vel 1 5 -node 2]"']
LOAD = [f"pattern Plain 5 5 {{ load 2 {ALPHA * M!r} }}"]


def recursion():
    b, g = 0.25, 0.5
    u = [0.0] * (N + 1); v = [0.0] * (N + 1); a = [0.0] * (N + 1)
    keff = K + ALPHA * M * g / (b * DT) + M / (b * DT * DT)
    for i in range(N):
        t1 = (i + 1) * DT
        rhs = (pload(t1) + ALPHA * M * vref(t1)
               + M * (u[i] / (b * DT * DT) + v[i] / (b * DT) + (1 / (2 * b) - 1) * a[i])
               + ALPHA * M * (g / (b * DT) * u[i] + (g / b - 1) * v[i] + DT * (g / (2 * b) - 1) * a[i]))
        u[i + 1] = rhs / keff
        a[i + 1] = (u[i + 1] - u[i]) / (b * DT * DT) - v[i] / (b * DT) - (1 / (2 * b) - 1) * a[i]
        v[i + 1] = v[i] + DT * ((1 - g) * a[i] + g * a[i + 1])
    return u, v


def relerr(x, y):
    return max(abs(p - q) for p, q in zip(x, y)) / max(abs(q) for q in y)


@needs_exe
def test_sdof_recursion_and_load_oracle(tmp_path):
    dr, out, rc = run(tmp_path, "R", REF)
    assert rc == 0 and "done 0" in out, out
    da, out_a, rc_a = run(tmp_path, "A", LOAD)
    assert rc_a == 0 and "done 0" in out_a, out_a
    u, v = recursion()
    assert relerr(col(dr, "d.txt"), u[1:]) < 1e-9
    assert relerr(col(dr, "d.txt"), col(da, "d.txt")) < 1e-12


@needs_exe
def test_rayleigh_force_responses_are_relative(tmp_path):
    dr, out, rc = run(tmp_path, "R", REF)
    assert rc == 0, out
    v = col(dr, "v.txt")
    exp = [ALPHA * M * (vi - vref((i + 1) * DT)) for i, vi in enumerate(v)]
    assert relerr(col(dr, "rf.txt"), exp) < 1e-10              # reaction flag 2: +alphaM M (v - v_ref)
    assert relerr([-x for x in col(dr, "nrf.txt")], exp) < 1e-10  # getResponse: -alphaM M (v - v_ref)


@needs_exe
def test_energy_split(tmp_path):
    dr, out, rc = run(tmp_path, "R", REF)
    assert rc == 0, out
    assert "DW_rel W_ref" in out, "split columns not echoed:\n" + out
    da, _, _ = run(tmp_path, "A", LOAD)
    R = rows(dr, "e2.txt"); A = rows(da, "e2.txt")
    # R: time KE_ele KE_nod IE DW_ele DW_nod ULW DW_rel W_ref RES ERR ; A: same (split flag, columns zero)
    assert len(R[0]) == 11 and len(A[0]) == 11
    dw_rel = [r[7] for r in R]
    assert min(dw_rel) >= 0.0
    assert all(b >= a - 1e-15 for a, b in zip(dw_rel, dw_rel[1:]))
    scale = max(abs(x) for r in R for x in r[1:9])
    assert max(abs(r[7] + r[8] - r[5]) for r in R) < 1e-12 * scale
    assert all(r[7] == 0.0 and r[8] == 0.0 for r in A)            # no reference -> zero split
    assert max(abs(r[9] - a[9]) for r, a in zip(R, A)) < 1e-12 * scale   # RES identical to the oracle
    # DW_nod(A) - W_A = DW_nod(R), with W_A = ULW(A) - ULW(R)
    assert abs((A[-1][5] - (A[-1][6] - R[-1][6])) - R[-1][5]) < 1e-12 * scale


@needs_exe
def test_zero_reference_is_bitwise_absent(tmp_path):
    d0, out0, rc0 = run(tmp_path, "none", [])
    dz, outz, rcz = run(tmp_path, "zero", [f"{CMD} 1 -vel 1 7 -node 2"])
    assert rc0 == 0 and rcz == 0, out0 + outz
    for f in ("d.txt", "v.txt", "rf.txt", "nrf.txt"):
        assert (d0 / f).read_bytes() == (dz / f).read_bytes(), f
    e0 = rows(d0, "e2.txt"); ez = rows(dz, "e2.txt")
    assert [r[:7] + r[-2:] for r in e0] == [r[:7] + r[-2:] for r in ez]


@needs_exe
def test_addD_force_integrators_are_refused(tmp_path):
    d, out, rc = run(tmp_path, "tp", REF, integ="HHT_TP 0.9")
    assert rc != 0 and "FATAL" in out and CMD in out, out
    d, out, rc = run(tmp_path, "hht", REF, integ="HHT 0.9")
    assert rc == 0 and "done 0" in out, out


@needs_exe
def test_command_audit_and_refusals(tmp_path):
    d, out, rc = run(tmp_path, "R", REF)
    rr = [l for l in out.splitlines() if l.startswith("RR")][0].split()[1:]
    assert [float(x) for x in rr] == [1.0, 0.0, pytest.approx(ALPHA * M)]
    for bad in ([f"{CMD} 1 -vel 1 5"],                          # no selector
                [f"{CMD} 1 -vel 1 99 -node 2"],                 # missing series
                [f"{CMD} 1 -node 2"],                           # no -vel
                [f"{CMD} 1 -vel 1 5 -node 2 -allMassNodes"],    # two selectors
                [f"{CMD} 1 -vel 1 5 -node 2", f"{CMD} 1 -vel 1 5 -node 2"]):   # duplicate tag
        d, out, rc = run(tmp_path, "bad%d" % len(list(tmp_path.iterdir())), bad)
        assert "done 0" not in out, bad


def test_openseespy_registration():
    pyd = REPO / "dist" / "bin" / "opensees.pyd"
    if not pyd.exists():
        pytest.skip("openseespy not built in this worktree")
    code = ("import sys; sys.path.insert(0, r'%s'); import opensees as ops; ops.wipe(); ops.model('basic','-ndm',1,'-ndf',1);"
            "ops.node(1,0.0,'-mass',2.0); ops.rayleigh(0.5,0.,0.,0.); ops.timeSeries('Constant',1);"
            "print('RR', ops.%s(1,'-vel',1,1,'-allMassNodes'))" % (pyd.parent, CMD))
    p = subprocess.run([sys.executable, "-S", "-c", code], stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT, encoding="utf-8", errors="replace", timeout=120,
                       env=dict(os.environ, LADRUNO_OPENSEES_QUIET="1"))
    line = [l for l in p.stdout.splitlines() if l.startswith("RR")]
    assert line and "1.0" in line[0], p.stdout
