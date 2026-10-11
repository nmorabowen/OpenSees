"""Viscous zeroLength dashpot and EnergyBalance: SDOF m=1 on zeroLength (spring k=100 + dashpot c=2), load pulse.
Variants: V0 Elastic+Viscous materials -doRayleigh 0 (two materials, same dir); P0 Parallel(Elastic,Viscous) -factors
-doRayleigh 0 (the T2S form); P1 the same with -doRayleigh 1 (Rayleigh factors all zero); S separate dashpot element."""
import subprocess, sys, pathlib, numpy as np
EXE = sys.argv[1]
import tempfile
H = pathlib.Path(tempfile.gettempdir()) / "relative_mass_damping" / "dashpot"
H.mkdir(parents=True, exist_ok=True)
m, k, c, dt, n = 1.0, 100.0, 2.0, 0.001, 5000
base = ["wipe", "model basic -ndm 1 -ndf 1", "node 1 0.0", "node 2 0.0 -mass 1.0", "fix 1 1",
        f"uniaxialMaterial Elastic 1 {k}", f"uniaxialMaterial Viscous 2 {c} 1.0",
        "uniaxialMaterial Elastic 11 1.0", "uniaxialMaterial Viscous 12 1.0 1.0",
        f"uniaxialMaterial Parallel 3 11 12 -factors {k} {c}"]
var = {"V0": ["element zeroLength 1 1 2 -mat 1 2 -dir 1 1 -doRayleigh 0"],
       "P0": ["element zeroLength 1 1 2 -mat 3 -dir 1 -doRayleigh 0"],
       "P1": ["element zeroLength 1 1 2 -mat 3 -dir 1 -doRayleigh 1"],
       "S":  ["element zeroLength 1 1 2 -mat 1 -dir 1 -doRayleigh 0", "element zeroLength 2 1 2 -mat 2 -dir 1 -doRayleigh 0"]}
for name, ele in var.items():
    d = H / name; d.mkdir(parents=True, exist_ok=True)
    L = base + ele + ["timeSeries Path 1 -time {0.0 0.05 0.1 100.0} -values {0.0 10.0 0.0 0.0}",
         "pattern Plain 1 1 { load 2 1.0 }",
         "recorder Node -file v.txt -time -precision 17 -node 2 -dof 1 vel",
         "recorder EnergyBalance -time -precision 12 -file e.txt",
         "recorder EnergyBalance -time -precision 12 -v2 -file e2.txt",
         "constraints Plain", "numberer Plain", "system BandGeneral", "test NormDispIncr 1e-14 20",
         "algorithm Newton", "integrator Newmark 0.5 0.25", "analysis Transient",
         f"set ok [analyze {n} {dt}]", 'puts "done $ok"', "wipe"]
    (d / "main.tcl").write_text("\n".join(L) + "\n")
    r = subprocess.run([EXE, "main.tcl"], cwd=d, capture_output=True, text=True)
    assert "done 0" in r.stdout + r.stderr, (r.stdout + r.stderr)[-1500:]
    v = np.loadtxt(d / "v.txt"); vv = np.r_[0.0, v[:, 1]]
    Dvis = np.sum(0.5 * (c * vv[:-1] ** 2 + c * vv[1:] ** 2) * dt)        # trapezoidal dashpot dissipation
    e = np.loadtxt(d / "e.txt")[-1]   # time KE IE DW ULW RES ERR
    print(f"{name}: KE {e[1]:.6e} IE {e[2]:.6e} DW {e[3]:.6e} ULW {e[4]:.6e} RES {e[5]:.6e} | D_visc {Dvis:.6e}  "
          f"RES/D_visc {e[5] / Dvis:+.6f}  DW/D_visc {e[3] / Dvis:.6f}")
