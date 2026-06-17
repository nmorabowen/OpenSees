"""LadrunoRCConcrete — meshed squat RC wall, quasi-static EXPLICIT cyclic (Zone-B, gmsh).

The RC shell ADR (19, Phase 2b.2c) established that cyclic pinching is a PANEL effect:
a single material point / single element under cyclic shear is concrete-damage-dominated
and the aggregate-interlock terms are nearly inert (the crack sits on the principal plane,
so the crack-tangential slip g_nt ~ 0). The distinctive cyclic interlock behaviour only
emerges when the principal directions ROTATE across a meshed wall.

It also established that the cyclic-softening solve will NOT converge with implicit Newton
(the consistent tangent goes indefinite on the softening branch). The fix is to drive the
wall as a QUASI-STATIC EXPLICIT problem (CentralDifferenceLadruno): no tangent is formed,
so the indefinite-tangent stall simply does not exist and the reversals integrate through.

These two cases prove both points end-to-end on a gmsh-meshed squat wall (aspect H/L=0.75)
under a held axial load + a reversing top-drift history:

1. **Explicit completes the cyclic softening history** — the full ±drift schedule runs to
   completion (both the cyclic and the monotone-bound material), where implicit Newton
   diverges at the first reversal. The headline robustness result.

2. **The cyclic interlock degradation is load-bearing at the panel scale** — with the
   cyclic friction-slip + X-crack wear ON (`-cyclic -xcrack`), the wall dissipates
   measurably LESS hysteretic energy and reaches a LOWER peak shear than the monotone
   `-interlock` bound (which cannot degrade across reversals). At single-element scale this
   ablation is invisible (damage-dominated); at the mesh scale it is a clear margin.

Honest scope: this demonstrates the panel MECHANISM + the degradation ablation + explicit
robustness. A quantitative pinched-loop match to a named experiment (Tran-Wallace) is a
further calibration step (specimen geometry + reinforcement layout) and is not asserted here.

Zone-B: needs gmsh; explicit quasi-static run is heavier than Zone-A, runs nightly.
"""
import math

import pytest

from _testbed import ops

gmsh = pytest.importorskip("gmsh")
pytestmark = [pytest.mark.zone_b]

# --- materials (N, mm) ---
E, NU, KC = 30000.0, 0.2, 2.0 / 3.0
CE = [0.0, 0.0007, 0.0020, 0.0100]; CS = [0.0, 24.0, 30.0, 8.0]
CD = [0.0, 0.0, 0.25, 1.0 - 8.0 / (E * 0.0100)]
FT = 3.0
TE = [0.0, FT / E, 2.0e-3, 1.5e-2]; TS = [0.0, FT, 1.8, 0.15]
TD = [0.0, 0.0, 1.0 - 1.8 / (E * 2.0e-3), 1.0 - 0.15 / (E * 1.5e-2)]
Lw, Hw, T = 2000.0, 1500.0, 150.0      # squat wall, aspect H/L = 0.75
NLAY = 4
RHO = 2.4e-9                            # tonne/mm^3 -> element mass -> stable explicit dt


def _mesh(h):
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("wall")
        occ = gmsh.model.occ
        occ.addRectangle(0, 0, 0, Lw, Hw)
        occ.synchronize()
        nx, ny = max(2, round(Lw / h)), max(2, round(Hw / h))
        for c in gmsh.model.getEntities(1):
            bb = gmsh.model.getBoundingBox(1, c[1])
            horiz = (bb[3] - bb[0]) > (bb[4] - bb[1])
            gmsh.model.mesh.setTransfiniteCurve(c[1], (nx + 1) if horiz else (ny + 1))
        gmsh.model.mesh.setTransfiniteSurface(1)
        gmsh.model.mesh.setRecombine(2, 1)
        occ.synchronize()
        gmsh.model.mesh.generate(2)
        nt, nc, _ = gmsh.model.mesh.getNodes()
        nodes = {int(t): (nc[3 * i], nc[3 * i + 1]) for i, t in enumerate(nt)}
        quads = []
        for et, _, en in zip(*gmsh.model.mesh.getElements(2)):
            if et == 3:
                en = [int(x) for x in en]
                quads += [en[4 * k:4 * k + 4] for k in range(len(en) // 4)]
        return nodes, quads
    finally:
        gmsh.finalize()


def _mat(tag, cyclic, xcrack):
    a = ["LadrunoRCConcrete", tag, E, NU, "-Ce", *CE, "-Cs", *CS, "-Cd", *CD,
         "-Te", *TE, "-Ts", *TS, "-Td", *TD, "-Kc", KC, "-beta", "-rho", RHO,
         "-interlock", "-agg", 16.0]
    if cyclic:
        a += ["-cyclic"]
    if xcrack:
        a += ["-xcrack"]
    ops.nDMaterial(*a)


def _run(h, cyclic, xcrack, peaks, axial=3.0, steps_per_seg=2500, dt_frac=0.2):
    """Build the meshed wall and drive a reversing top-drift history under
    CentralDifferenceLadruno (quasi-static explicit). Returns
    (loop[(u_top, base_shear)], reached_drift, n_quads)."""
    nodes, quads = _mesh(h)
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    for t, (x, y) in nodes.items():
        ops.node(int(t), x, y, 0.0)
    _mat(1, cyclic, xcrack)
    layers = []
    for _ in range(NLAY):
        layers += [1, T / NLAY]
    ops.section("LayeredShell", 1, NLAY, *layers)
    for i, q in enumerate(quads, 1):
        ops.element("ASDShellQ4", i, *[int(n) for n in q], 1)
    ymin = min(y for _, y in nodes.values()); ymax = max(y for _, y in nodes.values())
    base = [int(t) for t, (x, y) in nodes.items() if abs(y - ymin) < 1e-6]
    top = [int(t) for t, (x, y) in nodes.items() if abs(y - ymax) < 1e-6]
    for t in base:
        ops.fix(t, 1, 1, 1, 1, 1, 1)
    for t in nodes:                            # planar membrane: lock z + all rotations
        if int(t) not in base:
            ops.fix(int(t), 0, 0, 1, 1, 1, 1)

    # ASDShellQ4 provides no per-element dt_cr -> manual wave-speed bound dt ~ frac*h/c
    c = math.sqrt(E / RHO)
    dt = dt_frac * (h / c)
    seg_time = steps_per_seg * dt

    Fax = axial * Lw * T                        # held axial compression on the top edge
    ops.timeSeries("Constant", 1); ops.pattern("Plain", 1, 1)
    for t in top:
        ops.load(t, 0.0, -Fax / len(top), 0, 0, 0, 0)

    NSEG = 40                                   # smooth cosine drift, dt-based time grid
    seg = []; prev = 0.0
    for pk in peaks:
        for k in range(1, NSEG + 1):
            seg.append(prev + (pk - prev) * (1 - math.cos(math.pi * k / NSEG)) / 2.0)
        prev = pk
    node_dt = seg_time / NSEG
    times = [i * node_dt for i in range(len(seg) + 1)]
    ops.timeSeries("Path", 2, "-time", *times, "-values", 0.0, *seg, "-useLast")
    ops.pattern("Plain", 2, 2)
    for t in top:                               # rigid top: same prescribed u_x on all top nodes
        ops.sp(t, 1, 1.0)

    ops.constraints("Transformation"); ops.numberer("RCM")
    ops.system("Diagonal"); ops.algorithm("Linear")
    ops.integrator("CentralDifferenceLadruno", "-cflAbort", "-lump", "diagonal")
    ops.analysis("Transient")

    target_end = times[-1]
    loop = []; last_u = 0.0; done = 0.0; step = 0
    sample_every = max(1, steps_per_seg // 30)
    maxstep = len(peaks) * steps_per_seg + 50
    while done < target_end - dt and step < maxstep:
        if ops.analyze(1, dt) != 0:
            break                               # graceful stop (explicit should not need this)
        done += dt; step += 1
        if step % sample_every == 0:
            ops.reactions()
            bs = sum(ops.nodeReaction(t)[0] for t in base)
            last_u = ops.nodeDisp(top[0], 1)
            loop.append((last_u, -bs))
            if abs(bs) > 1e10:                  # instability tell
                break
    reached = max((abs(u) for u, _ in loop), default=0.0)
    return loop, reached, len(quads)


def _hyst_energy(loop):
    fs = [f for _, f in loop]; us = [u for u, _ in loop]
    return abs(sum(0.5 * (fs[i] + fs[i - 1]) * (us[i] - us[i - 1]) for i in range(1, len(loop))))


# ===========================================================================
def test_squat_wall_explicit_completes_cyclic_softening():
    """The headline robustness result: the full reversing drift schedule runs to
    completion under quasi-static explicit on the meshed squat wall — the reversals
    that diverge implicit Newton integrate through cleanly."""
    peaks = [4.0, -4.0, 8.0, -8.0]
    loop, reached, nq = _run(500.0, cyclic=True, xcrack=True, peaks=peaks)
    assert nq >= 9, f"mesh too coarse ({nq} quads)"
    assert reached > 7.0, f"explicit did not complete the cyclic history (reached {reached:.2f} mm of 8)"
    fpk = max(abs(f) for _, f in loop)
    assert 1.0e5 < fpk < 1.0e7, f"base shear {fpk:.0f} N out of physical range (instability?)"


def test_squat_wall_cyclic_interlock_degradation_is_load_bearing():
    """Panel-scale ablation: with the cyclic friction-slip + X-crack wear ON the wall
    dissipates LESS hysteretic energy and reaches a LOWER peak shear than the monotone
    -interlock bound (which cannot degrade across reversals). Invisible at single-element
    scale (damage-dominated); a clear margin once principal rotation engages the interlock."""
    peaks = [4.0, -4.0, 8.0, -8.0]
    cyc, rc, _ = _run(500.0, cyclic=True, xcrack=True, peaks=peaks)
    mon, rm, _ = _run(500.0, cyclic=False, xcrack=False, peaks=peaks)
    assert rc > 7.0 and rm > 7.0, f"a run did not complete (cyclic {rc:.2f}, monotone {rm:.2f})"
    W_cyc, W_mon = _hyst_energy(cyc), _hyst_energy(mon)
    F_cyc = max(abs(f) for _, f in cyc); F_mon = max(abs(f) for _, f in mon)
    assert W_cyc < 0.92 * W_mon, (
        f"cyclic degradation not load-bearing: hysteretic energy {W_cyc:.0f} vs monotone "
        f"{W_mon:.0f} N*mm (expected a clear reduction)")
    assert F_cyc < F_mon, (
        f"cyclic peak shear {F_cyc:.0f} not below monotone bound {F_mon:.0f} N")
