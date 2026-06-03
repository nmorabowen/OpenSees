"""Mesh-convergence study for the notched-bar single-point over-prediction (README §4.6).

Pulls the DEN notched bar (J2 plasticity, NO damage -> stable at fine mesh) at h=4,2,1,0.5
with both the single-point `ssp` formulation (a.k.a. the deprecated `eas` name) and the
fully-integrated `bbar`, records load@u=0.4 and wall time per solve, and writes
figures/lemaitre_C_convergence.png. Confirms both converge to a common ~20 kN limit and the
single-point over-prediction (17.4% -> 3.7%) vanishes ~O(h).

Usage:  python Ladruno_scripts/lemaitre_convergence_study.py  <dist\\bin>
"""
import os, sys, time, math
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
DIST = sys.argv[1] if len(sys.argv) > 1 else \
    r"C:\Users\nmora\Github\OpenSees_Compile\OpenSees\.claude\worktrees\wonderful-nobel-1377f4\dist\bin"
os.add_dll_directory(DIST); sys.path.insert(0, DIST)
sys.path.insert(0, os.path.join(ROOT, "tests")); os.environ["LADRUNO_OPENSEES_QUIET"] = "1"
import opensees as ops
import test_lemaitre_notched_bar as C
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

ISO = ("-iso", "voce", 400.0, 120.0, 20.0, 200.0)        # J2 + hardening, no damage
FIG = os.path.join(ROOT, "Ladruno_implementation", "lemaitre_validation", "figures",
                   "lemaitre_C_convergence.png")


def _form_name():
    """Use 'ssp' if the binary supports it (post eas->ssp rename), else the deprecated 'eas'."""
    ops.wipe(); ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, (x, y, z) in {1:(0,0,0),2:(1,0,0),3:(1,1,0),4:(0,1,0),5:(0,0,1),6:(1,0,1),7:(1,1,1),8:(0,1,1)}.items():
        ops.node(t, float(x), float(y), float(z))
    ops.nDMaterial("ElasticIsotropic", 1, 1.0e5, 0.3)
    try:
        ops.element("LadrunoBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1, "-formulation", "ssp")
        return "ssp"
    except Exception:
        return "eas"


SP = _form_name()


def solve(h, form, umax=0.4, n=40):
    nodes, hexes = C._mesh(h)
    keep = [hx for hx in hexes if not C._in_notch(sum(nodes[z][0] for z in hx)/8, sum(nodes[z][1] for z in hx)/8)]
    ops.wipe(); ops.model("basic", "-ndm", 3, "-ndf", 3)
    used = set(z for hx in keep for z in hx)
    for t in used: ops.node(int(t), *nodes[t])
    ops.nDMaterial("LadrunoJ2", 1, C.K, C.G, *ISO)
    for i, hx in enumerate(keep, 1):
        ops.element("LadrunoBrick", i, *[int(z) for z in hx], 1, "-formulation", form)
    xmin = min(nodes[t][0] for t in used); xmax = max(nodes[t][0] for t in used)
    left = sorted(int(t) for t in used if abs(nodes[t][0]-xmin) < 1e-6)
    for t in used:
        if abs(nodes[t][2]) < 1e-6: ops.fix(int(t), 0, 0, 1)
    for t in left: ops.fix(int(t), 1, 1 if t == left[0] else 0, 0)
    pull = sorted(int(t) for t in used if abs(nodes[t][0]-xmax) < 1e-6); m = pull[0]
    for nd in pull:
        if nd != m: ops.equalDOF(m, nd, 1)
    ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1); ops.load(m, 1.0, 0, 0)
    ops.constraints("Transformation"); ops.numberer("RCM"); ops.system("UmfPack")
    ops.test("NormDispIncr", 1e-6, 40, 0); ops.analysis("Static")
    base = umax/n; u = P = 0.0
    for _ in range(n):
        ops.integrator("DisplacementControl", m, 1, base); ops.algorithm("Newton")
        if ops.analyze(1) != 0: break
        ops.reactions(); P = -sum(ops.nodeReaction(z)[0] for z in left); u = ops.nodeDisp(m)[0]
    return P/1000.0, len(keep)


meshes = [4.0, 2.0, 1.0, 0.5]
res = {SP: {}, "bbar": {}}; wall = {SP: {}, "bbar": {}}; ne = {}
for h in meshes:
    for form in (SP, "bbar"):
        t0 = time.time(); P, n = solve(h, form); res[form][h] = P; ne[h] = n; wall[form][h] = time.time()-t0
        print(f"  {form:4s} h={h:<4} P={P:7.3f} kN  ne={n:5d}  {wall[form][h]:.0f}s", flush=True)

hs = sorted(meshes)
fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
xb = [1.0/h for h in hs]
ax[0].plot(xb, [res["bbar"][h] for h in hs], "o-", color="C0", label="bbar (8 GP)")
ax[0].plot(xb, [res[SP][h] for h in hs], "s--", color="C3", label=f"{SP} (1 pt)")
ax[0].axhline(res["bbar"][min(hs)], ls=":", color="0.5", label=f"common limit ~{res['bbar'][min(hs)]:.0f} kN")
ax[0].set(xlabel="mesh refinement  1/h", ylabel="notch load @ u=0.4 [kN]",
          title="C  notch-load convergence (J2 plasticity)\nsingle-point & bbar -> common limit as h->0")
ax[0].legend(); ax[0].grid(alpha=0.3)
gap = [100*(res[SP][h]/res["bbar"][h]-1.0) for h in hs]
ax[1].loglog(xb, gap, "D-", color="C2")
for h, g in zip(hs, gap): ax[1].annotate(f"{g:.1f}%", (1.0/h, g), fontsize=8, ha="left", va="bottom")
ax[1].plot(xb, [gap[0]*(xb[0]/x) for x in xb], ":", color="0.6", label="O(h) slope")
ax[1].set(xlabel="mesh refinement  1/h", ylabel="single-point over-prediction vs bbar [%]",
          title="C  over-prediction vanishes under refinement\n(~O(h))"); ax[1].legend(); ax[1].grid(alpha=0.3, which="both")
fig.tight_layout(); fig.savefig(FIG); plt.close(fig)
print("wrote", FIG)
