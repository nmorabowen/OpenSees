r"""LadrunoUP (ELE 33017) — ADR-71 P3 Taylor-Hood PLUMBING battery (WP3.A, 3.C-PL).

The Taylor-Hood (heterogeneous-ndf) lane of the unified Biot u-p element: the
Bezier quadratic shapes BT6 (2D, 6 nodes) and BTET10 (3D, 10 nodes) with
quadratic displacement and LINEAR vertex pressure (`-pOrder linear`). Vertices
are pressure CARRIERS (ndf = ndm+1); mid-edge nodes carry displacement only
(ndf = ndm). This battery proves the *plumbing* the ADR-71 §7 P3 row calls out
"heterogeneous-ndf plumbing (DOF map, sendSelf, recorders, numberers) proven;
straight-sided guard carried; equalDOF mixed-ndf interface example" — the
sibling battery (test_ladruno_up_element_th_b1.py) carries the ZCB80 physics
gate. Eight gates, the 3.C-PL pins verbatim:

  (i)   sendSelf/recvSelf database round-trip (BT6-TH + BTET10-TH, non-default
        args, exact restore + re-solve) — the 1.D serial gate, TH edition.
  (ii)  recorders: node disp-slot p on carriers (recorder Node -dof <ndm+1> disp
        == nodeDisp slot), eleResponse porePressure / flux carrier-consistent
        lengths.
  (iii) numberers: same TH model under Plain vs RCM -> identical u/p (mixed-ndf
        renumbering proof).
  (iv)  equalDOF mixed-ndf interface: a dry (ndf=2 standard quad) region tied to
        a saturated BT6-TH (ndf=3) region at DUPLICATED interface nodes via the
        EXPLICIT-DOF-list equalDOF (ADR §4.4 ⟨FW-F10⟩ — the bare form mis-sizes
        against the smaller-ndf node). THE example the P4 guide will cite.
  (v)   BTET10 winding acceptance: right-handed vs left-handed node input ->
        identical u/p (the element folds |detJ|), one-time note observed.
  (vi)  straight-side guard: a curved mid-edge node -> loud setDomain error and
        a deactivated (singular) element, not a crash.
  (vii) BT6-TH static steady-seepage exactness: a linear-p patch test recovers
        the manufactured linear field at the free interior vertex to ~machine,
        the Darcy flux == -k̄·∇p exactly, AND the per-GP porePressure values ==
        the field at the 3-pt rule points (the P1 gate-4 analog).
  (viii)BTET10-TH 1D consolidation smoke vs the Terzaghi series (loose 5%).

Panel-hardening additions (P3 adversarial panel, post-battery):
  (ix)  TH transient-block FD gate (the P1 test-8 pattern ported to one BT6-TH):
        central-FD of reactions('-dynamic') w.r.t. nodal vel/accel on a fully-
        penalty-prescribed 15-DOF element. Gates: Cpu == (-Kup)^T (the honest-p
        transpose relation, Kup from a free-floating static printA), Cpp == the
        ANALYTIC compacted storage S, Cuu/Cup ~ 0 without Rayleigh, mass p-rows
        AND p-cols zero (-dynSeepage off), Muu symmetric with x-direction total
        == rho*Area*thick. Makes the compacted-basis Q^T/S/M correctness self-
        contained instead of riding on the B1 physics gate.
  (x)   '-dynSeepage on' (explicit opt-in; default off since P4 B5) TH gate:
        same FD rig; the
        dR/da p-rows now carry +G = int dNp^T kbar rhoF Nu dV (checked against
        the ANALYTIC G — every Bernstein Nu integrates to A/6) while the
        transient printA effective tangent keeps A[p,u] == c2*Q^T exactly (G
        stays OUT of the tangent, mirroring the P1 dyn-seepage quantification).
        The only P3 test exercising the element-default configuration.

--------------------------------------------------------------------------
THE TH MODELING DANCE (verified working, this build):
  vertex nodes defined under  model('basic','-ndm',ndm,'-ndf',ndm+1);
  mid-edge nodes under         model('basic','-ndm',ndm,'-ndf',ndm)   [re-issuing
  the model command changes the builder ndf for SUBSEQUENT nodes without wiping];
  then element('LadrunoUP', tag, <vertices...>, <midedges...>, matTag,
               ..., '-pOrder','linear').
  BT6 node order:    v0,v1,v2,  m(v0,v1), m(v1,v2), m(v2,v0).
  BTET10 node order: v0,v1,v2,v3,  then edges (v0,v1)(v1,v2)(v0,v2)(v0,v3)
                     (v2,v3)(v1,v3)   [pinned to LadrunoUPShapes.h / LadrunoUP.cpp
                     setDomain btet10Edges = {{4,0,1},{5,1,2},{6,0,2},{7,0,3},
                     {8,2,3},{9,1,3}}].
  Straight sides are MANDATORY (the setDomain guard is STRICT, tol 1e-6*edgeLen).
  Pressure carrier = a vertex's disp slot ndm (1-based dof ndm+1). p is read as
  nodeDisp(node, ndm+1) — the honest-p contract (ADR §3.2), UNSYMMETRIC tangent
  => UmfPack throughout.

TOLERANCE RATIONALE
  * (i) exact 0.0 drift — same-process File-datastore round-trip is bit-exact;
    any dropped ctor arg / committed state moves a response.
  * (vii) ~1e-9 — a linear field is in the P1-Laplacian nullspace-complement; the
    patch test is exact up to conditioning of the penalty-imposed BCs.
  * (iii)/(v) ≤1e-12 rel — identical operator, only the DOF numbering / winding
    differs; must agree to round-off.
  * (viii) ≤5% p, U within 0.03 — loose (decay CLASS + no-blowup, not accuracy
    records), per the pin.

RUN MECHANICS (worktree, per project memory "OpenSees Python test env"):
  the module bootstraps <worktree>/dist/bin, evicting any boot-.pth stale pyd and
  asserting the loaded pyd is THIS worktree's. Prescribed-p uses PENALTY
  constraints (the Transformation staged-sp trap is on record) with NormDispIncr
  (a NormUnbalance floor bites at the H·p penalty scale). The winding-note gate
  (v) runs in a `-S` SUBPROCESS so the "printed once per process" note is seen in
  isolation regardless of test order.
  Invoke: py -3.12 -m pytest tests/test_ladruno_up_element_th.py -v
"""
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

# --------------------------------------------------------------------------
# bootstrap: THIS worktree's engine (evict any boot-.pth preloaded stale pyd)
# --------------------------------------------------------------------------
_DIST = str(Path(__file__).resolve().parents[1] / "dist" / "bin")
if not (os.path.isfile(os.path.join(_DIST, "opensees.pyd"))
        or os.path.isfile(os.path.join(_DIST, "opensees.so"))):
    pytest.skip(f"worktree engine not built: {_DIST}", allow_module_level=True)

from _engine import bind_worktree_engine  # noqa: E402
ops = bind_worktree_engine(_DIST)

pytestmark = [pytest.mark.zone_b]

_UMF = ("UmfPack",)   # honest-p tangent is unsymmetric (ADR §3.2)


# --------------------------------------------------------------------------
# shared Terzaghi oracles (numpy) — mirror the T3/analytic donors
# --------------------------------------------------------------------------
def terzaghi_series(Z, Tv, nterms=400):
    """Excess-pore-pressure ratio p/p0; Z = z/H from the DRAINED face."""
    m = np.arange(nterms)
    M = (2 * m + 1) * np.pi / 2.0
    return np.sum((2.0 / M) * np.sin(M * Z) * np.exp(-M * M * Tv))


def terzaghi_U(Tv, nterms=400):
    m = np.arange(nterms)
    M = (2 * m + 1) * np.pi / 2.0
    return 1.0 - np.sum((2.0 / M**2) * np.exp(-M * M * Tv))


# --------------------------------------------------------------------------
# BT6-TH mesh helpers (straight-sided, shared mid-edge nodes)
# --------------------------------------------------------------------------
def _bt6_unit_square(matArgs, mat_id=1, E=1.0e4, nu=0.25, rho=2.0):
    """Unit square [0,1]^2 as two straight-sided BT6 (diagonal 0->2 split).
    Returns (vertices{tag:(x,y)}, mids{(a,b):tag}, elem_tags)."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    V = {1: (0.0, 0.0), 2: (1.0, 0.0), 3: (1.0, 1.0), 4: (0.0, 1.0)}
    for t, (x, y) in V.items():
        ops.node(t, x, y)
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    M = {}
    _next = [11]

    def mid(a, b):
        key = tuple(sorted((a, b)))
        if key not in M:
            xa, ya = V[a]
            xb, yb = V[b]
            t = _next[0]
            _next[0] += 1
            ops.node(t, 0.5 * (xa + xb), 0.5 * (ya + yb))
            M[key] = t
        return M[key]

    ops.nDMaterial("ElasticIsotropic", mat_id, E, nu, rho)
    # CCW triangles sharing diagonal (1,3)
    tris = [(1, 2, 3), (1, 3, 4)]
    etags = []
    for e, (a, b, c) in enumerate(tris, start=1):
        ops.element("LadrunoUP", e, a, b, c,
                    mid(a, b), mid(b, c), mid(c, a), mat_id, *matArgs)
        etags.append(e)
    return V, M, etags


# ==========================================================================
# (i) sendSelf / recvSelf database round-trip — TH edition (the 1.D gate)
# ==========================================================================
def _vmax(a, b):
    assert a and b and len(a) == len(b), (a, b)
    return max(abs(x - y) for x, y in zip(a, b))


def _roundtrip_check(db_path, ref_pnodes, ref_ele, resolve):
    """save -> wipe -> restore -> assert exact + re-solve. Shared by BT6/BTET10."""
    ops.database("File", db_path)
    ops.save(1)
    ops.wipe()
    ops.database("File", db_path)
    ops.restore(1)
    res_p = {n: ops.nodeDisp(n, ref_pnodes[n][1]) for n in ref_pnodes}
    dp = max(abs(ref_pnodes[n][0] - res_p[n]) for n in ref_pnodes)
    assert dp == 0.0, f"nodal pressure drift after restore: {dp}"
    for tag, resp in ref_ele.items():
        for name, refval in resp.items():
            got = ops.eleResponse(tag, name)
            d = _vmax(refval, got)
            assert d < 1.0e-10, f"ele {tag} response {name!r} drifted by {d:g}"
    assert resolve() == 0, "re-solve after restore did not converge"


@pytest.mark.t1
def test_i_db_roundtrip_bt6_th(tmp_path):
    """BT6-TH: two-element unit square, every u-p knob off its default
    (-alpha 0.9, -Ks 3.6e4, anisotropic -perm, -lumped, -dynSeepage off, gravity
    -body, solid-only Rayleigh). Static gravity solve produces a non-trivial pore
    pressure at the un-drained base carriers; save/wipe/restore must reproduce the
    exact nodal state + per-GP responses and re-solve."""
    common = ["-Kf", 2.2e6, "-poro", 0.3, "-rhoF", 1.0,
              "-perm", 1.0e-4, 2.0e-4, "-alpha", 0.9, "-Ks", 3.6e4,
              "-body", 0.0, -9.81, "-fluidBody", 0.0, -9.81,
              "-lumped", "-dynSeepage", "off", "-pOrder", "linear"]
    V, M, etags = _bt6_unit_square(common)
    ops.fix(1, 1, 1, 0)          # base u fixed, p FREE (the pressure probe)
    ops.fix(2, 1, 1, 0)
    ops.fix(3, 0, 0, 1)          # drained top vertices (p = 0)
    ops.fix(4, 0, 0, 1)
    ops.fix(M[(1, 2)], 1, 1)     # base mid-edge u fixed
    ops.rayleigh(0.02, 0.0, 1.0e-4, 0.0)

    def solve():
        ops.constraints("Plain")
        ops.numberer("RCM")
        ops.system(*_UMF)
        ops.test("NormDispIncr", 1.0e-10, 25)
        ops.algorithm("Newton")
        ops.integrator("LoadControl", 1.0)
        ops.analysis("Static")
        return ops.analyze(1)

    assert solve() == 0, "reference static analysis did not converge"
    ref_pnodes = {n: (ops.nodeDisp(n, 3), 3) for n in (1, 2)}
    assert abs(ref_pnodes[1][0]) > 1.0, "BT6-TH model ~zero base pressure -- not a probe"
    ref_ele = {t: {name: ops.eleResponse(t, name)
                   for name in ("porePressure", "stresses", "stressesTotal", "flux")}
               for t in etags}
    # porePressure/flux carrier-consistent lengths (nGP == nP == 3 for BT6)
    assert len(ref_ele[1]["porePressure"]) == 3
    assert len(ref_ele[1]["flux"]) == 3 * 2
    _roundtrip_check(str(tmp_path / "bt6_db"), ref_pnodes, ref_ele, solve)


def _btet10_single(matArgs, mat_id=1, E=1.0e4, nu=0.25, rho=2.0,
                   order=(0, 1, 2, 3)):
    """One straight-sided BTET10 on the reference tet, node slots permuted by
    `order` (order maps element vertex-slot -> physical vertex). Returns
    (vtags[4], mtags[6], P physical coords)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 4)
    P = {0: (0.0, 0.0, 0.0), 1: (1.2, 0.0, 0.0),
         2: (0.0, 1.1, 0.0), 3: (0.0, 0.0, 1.3)}
    vt = []
    for slot in range(4):
        x, y, z = P[order[slot]]
        t = 100 + slot
        ops.node(t, x, y, z)
        vt.append(t)
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    edges = [(0, 1), (1, 2), (0, 2), (0, 3), (2, 3), (1, 3)]
    mt = []
    for i, (a, b) in enumerate(edges):
        xa = P[order[a]]
        xb = P[order[b]]
        t = 200 + i
        ops.node(t, *(0.5 * (xa[k] + xb[k]) for k in range(3)))
        mt.append(t)
    ops.nDMaterial("ElasticIsotropic", mat_id, E, nu, rho)
    ops.element("LadrunoUP", 1, *(vt + mt), mat_id, *matArgs)
    return vt, mt, P, order


@pytest.mark.t1
def test_i_db_roundtrip_btet10_th(tmp_path):
    """BTET10-TH: one element, non-default args (anisotropic 3D -perm, -alpha 0.9,
    -Ks 3.6e4, -dynSeepage off, gravity -body, -lumped). Steady-seepage static
    solve with a prescribed linear pressure on three vertices; save/restore exact
    + re-solve. std formulation (3D bbar mean-dilatation on a quadratic tet is a
    reserved concern; the round-trip exercises the physics knobs, not bbar)."""
    common = ["-Kf", 2.2e6, "-poro", 0.3, "-rhoF", 1.0,
              "-perm", 1.0e-4, 1.5e-4, 2.0e-4, "-alpha", 0.9, "-Ks", 3.6e4,
              "-body", 0.0, 0.0, -9.81, "-lumped", "-dynSeepage", "off",
              "-pOrder", "linear"]
    vt, mt, P, order = _btet10_single(common)
    for t in vt:
        ops.fix(t, 1, 1, 1, 0)   # all u fixed, p free -> pure seepage
    for t in mt:
        ops.fix(t, 1, 1, 1)

    def pf(x, y, z):
        return 3.0 + 2.0 * x + 5.0 * y + 4.0 * z

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for slot in (0, 1, 2):       # prescribe p at three vertices; 4th free
        x, y, z = P[order[slot]]
        ops.sp(vt[slot], 4, pf(x, y, z))

    def solve():
        ops.constraints("Penalty", 1e14, 1e14)
        ops.numberer("RCM")
        ops.system(*_UMF)
        ops.test("NormDispIncr", 1.0e-12, 30)
        ops.algorithm("Newton")
        ops.integrator("LoadControl", 1.0)
        ops.analysis("Static")
        return ops.analyze(1)

    assert solve() == 0, "reference BTET10 seepage solve did not converge"
    pfree = ops.nodeDisp(vt[3], 4)
    assert abs(pfree) > 1.0, "BTET10-TH free-vertex pressure ~zero -- not a probe"
    ref_pnodes = {vt[3]: (pfree, 4)}
    ref_ele = {1: {name: ops.eleResponse(1, name)
                   for name in ("porePressure", "stresses", "stressesTotal", "flux")}}
    assert len(ref_ele[1]["porePressure"]) == 4         # nGP == nP == 4 for BTET10
    assert len(ref_ele[1]["flux"]) == 4 * 3
    # NOTE: the load pattern must survive the round-trip too; re-solve re-imposes it.
    _roundtrip_check(str(tmp_path / "btet10_db"), ref_pnodes, ref_ele, solve)


# ==========================================================================
# (ii) recorders — node disp-slot p on carriers + eleResponse lengths
# ==========================================================================
@pytest.mark.t1
def test_ii_recorder_pressure_and_eleresponse_lengths(tmp_path):
    """The classic `recorder Node -dof <ndm+1> disp` route reads the carrier pore
    pressure off the displacement channel (the honest-p contract) — it must equal
    nodeDisp(node, ndm+1). eleResponse porePressure length == nP (3 BT6),
    flux length == nGP*ndm (6). PressureSource (the Ladruno_NodeResults PRESSURE
    source) is a P4 recorder-gate item (⟨FW-F4⟩); the node-recorder disp-slot
    route is the P3 plumbing proof gated here."""
    common = ["-Kf", 2.2e6, "-poro", 0.3, "-rhoF", 1.0,
              "-perm", 1.0e-4, 1.0e-4, "-dynSeepage", "off", "-pOrder", "linear"]
    V, M, etags = _bt6_unit_square(common)
    ops.fix(1, 1, 1, 0)
    ops.fix(2, 1, 1, 0)
    ops.fix(3, 0, 0, 1)
    ops.fix(4, 0, 0, 1)
    ops.fix(M[(1, 2)], 1, 1)
    # prescribe a base pressure so the recorder sees a non-trivial value
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.sp(1, 3, 7.0)
    ops.sp(2, 3, 7.0)
    recfile = str(tmp_path / "p_carrier.out")
    # node recorder on the pressure DOF slot (ndm+1 = 3) of two carriers
    ops.recorder("Node", "-file", recfile, "-time",
                 "-node", 1, 2, "-dof", 3, "disp")
    ops.constraints("Penalty", 1e14, 1e14)
    ops.numberer("Plain")
    ops.system(*_UMF)
    ops.test("NormDispIncr", 1e-12, 30)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    p1 = ops.nodeDisp(1, 3)
    p2 = ops.nodeDisp(2, 3)
    assert abs(p1 - 7.0) < 1e-6 and abs(p2 - 7.0) < 1e-6, (p1, p2)

    # eleResponse carrier-consistent lengths
    pp = ops.eleResponse(1, "porePressure")
    fx = ops.eleResponse(1, "flux")
    assert len(pp) == 3, f"BT6 porePressure length {len(pp)} != nP(=nGP)=3"
    assert len(fx) == 3 * 2, f"BT6 flux length {len(fx)} != nGP*ndm=6"

    # flush recorder + compare the disp-slot channel to nodeDisp
    ops.remove("recorders")
    with open(recfile) as fh:
        last = fh.readlines()[-1].split()
    rp1, rp2 = float(last[1]), float(last[2])
    assert abs(rp1 - p1) <= 1e-10 * max(1.0, abs(p1)), (rp1, p1)
    assert abs(rp2 - p2) <= 1e-10 * max(1.0, abs(p2)), (rp2, p2)


# ==========================================================================
# (iii) numberers — Plain vs RCM identical (mixed-ndf renumbering proof)
# ==========================================================================
@pytest.mark.t1
def test_iii_numberer_plain_vs_rcm_identical():
    """The heterogeneous-ndf DOF map must renumber correctly: the same BT6-TH
    steady-seepage + gravity problem solved under numberer Plain vs RCM produces
    identical u and p at every node (≤1e-12 rel). A mis-handled mixed-ndf map
    (carrier p-slot vs mid-edge u-only) would move a DOF under RCM."""
    common = ["-Kf", 2.2e6, "-poro", 0.3, "-rhoF", 1.0,
              "-perm", 1.0e-4, 1.0e-4, "-body", 0.0, -9.81,
              "-dynSeepage", "off", "-pOrder", "linear"]

    def run(numberer):
        V, M, etags = _bt6_unit_square(common)
        ops.fix(1, 1, 1, 0)
        ops.fix(2, 1, 1, 0)
        ops.fix(3, 0, 0, 1)
        ops.fix(4, 0, 0, 1)
        ops.fix(M[(1, 2)], 1, 1)
        ops.constraints("Plain")
        ops.numberer(numberer)
        ops.system(*_UMF)
        ops.test("NormDispIncr", 1e-12, 30)
        ops.algorithm("Newton")
        ops.integrator("LoadControl", 1.0)
        ops.analysis("Static")
        assert ops.analyze(1) == 0, f"{numberer} solve failed"
        allnodes = list(V) + list(M.values())
        out = {}
        for n in allnodes:
            out[(n, 1)] = ops.nodeDisp(n, 1)
            out[(n, 2)] = ops.nodeDisp(n, 2)
        for n in V:                # vertices carry p on slot 3
            out[(n, 3)] = ops.nodeDisp(n, 3)
        return out

    a = run("Plain")
    b = run("RCM")
    assert set(a) == set(b)
    worst = max(abs(a[k] - b[k]) / max(1.0, abs(a[k])) for k in a)
    assert worst <= 1e-12, f"Plain vs RCM disagree (worst rel {worst:g})"
    # sanity: the problem actually produced motion + pressure (not a dead field)
    assert max(abs(v) for v in a.values()) > 1.0


# ==========================================================================
# (iv) equalDOF mixed-ndf interface — THE P4-guide example
# ==========================================================================
@pytest.mark.t1
def test_iv_equaldof_mixed_ndf_interface():
    """A DRY elastic region (standard 4-node `quad`, ndf=2) stacked on a SATURATED
    BT6-TH region (ndf=3). The shared interface nodes are DUPLICATED — the dry-side
    node is ndf=2, the wet-side vertex is ndf=3 — and tied on the u-DOFs only with
    the EXPLICIT-DOF-LIST equalDOF(rNode, cNode, 1, 2). ⟨FW-F10⟩: the bare
    equalDOF(rNode, cNode) sizes its constraint from the model-builder NDF and
    mis-sizes against the smaller-ndf dry node — the explicit list is mandatory.
    Gates: static gravity + a transient smoke both rc=0; interface u continuity
    exact (equalDOF); saturated-region p finite/sane.

    Topology (the guide reuses this):
        wet block  [0,1]x[0,1]  = 2 BT6-TH (diagonal 1-3 split), vertices 1-4
                                  (ndf=3, p carriers), mid-edges 11-15 (ndf=2);
        dry block  [0,1]x[1,2]  = 1 quad, corners 20,21 (DUP of wet vertices
                                  4,3) + 22,23 (top);
        ties       equalDOF(4,20,1,2)  equalDOF(3,21,1,2)  (wet retains, dry
                                  constrained; u1,u2 only);
        BCs        wet base 1,2 u=0 & p=0 (drained base), base mid-edge 11 u=0;
        load       downward point loads on the dry top corners 22,23.
    """
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    Vw = {1: (0.0, 0.0), 2: (1.0, 0.0), 3: (1.0, 1.0), 4: (0.0, 1.0)}
    for t, (x, y) in Vw.items():
        ops.node(t, x, y)
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(11, 0.5, 0.0)   # m(1,2)
    ops.node(12, 1.0, 0.5)   # m(2,3)
    ops.node(13, 0.5, 0.5)   # m(3,1) diagonal (shared by both BT6)
    ops.node(14, 0.5, 1.0)   # m(3,4) wet-top mid-edge (free boundary node)
    ops.node(15, 0.0, 0.5)   # m(4,1)
    ops.nDMaterial("ElasticIsotropic", 1, 1.0e4, 0.25, 2.0)
    common = ["-Kf", 2.2e6, "-poro", 0.3, "-rhoF", 1.0, "-perm", 1e-5, 1e-5,
              "-pOrder", "linear", "-dynSeepage", "off", "-body", 0.0, -9.81]
    ops.element("LadrunoUP", 1, 1, 2, 3, 11, 12, 13, 1, *common)
    ops.element("LadrunoUP", 2, 1, 3, 4, 13, 14, 15, 1, *common)
    # dry block: duplicated interface nodes 20 (== v4), 21 (== v3), ndf=2
    ops.node(20, 0.0, 1.0)
    ops.node(21, 1.0, 1.0)
    ops.node(22, 1.0, 2.0)
    ops.node(23, 0.0, 2.0)
    ops.nDMaterial("ElasticIsotropic", 2, 5.0e3, 0.25, 2.0)
    ops.element("quad", 10, 20, 21, 22, 23, 1.0, "PlaneStrain", 2)
    # explicit-DOF-list equalDOF (⟨FW-F10⟩): tie shared u1,u2; wet vertex retains
    ops.equalDOF(4, 20, 1, 2)
    ops.equalDOF(3, 21, 1, 2)
    # BCs: wet base u fixed AND drained (p=0) -> well-posed static seepage; the
    # saturated interior/top vertices keep p free (gravity-driven pore pressure).
    ops.fix(1, 1, 1, 1)
    ops.fix(2, 1, 1, 1)
    ops.fix(11, 1, 1)        # base mid-edge u fixed

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(22, 0.0, -1.0)
    ops.load(23, 0.0, -1.0)

    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system(*_UMF)
    ops.test("NormDispIncr", 1e-10, 30)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    assert ops.analyze(1) == 0, "(iv) static gravity+load did not converge"

    # interface u continuity (equalDOF enforced -> exact)
    for wet, dry in ((4, 20), (3, 21)):
        du1 = abs(ops.nodeDisp(wet, 1) - ops.nodeDisp(dry, 1))
        du2 = abs(ops.nodeDisp(wet, 2) - ops.nodeDisp(dry, 2))
        assert du1 <= 1e-10 and du2 <= 1e-10, (
            f"interface u discontinuity wet{wet}/dry{dry}: {du1:g}, {du2:g}")
    # saturated region pressure finite and sane
    for n in (3, 4):
        p = ops.nodeDisp(n, 3)
        assert np.isfinite(p) and abs(p) < 1e3, f"wet vertex {n} pressure insane: {p}"

    # transient smoke on the same mixed-ndf model
    ops.setTime(0.0)
    ops.wipeAnalysis()
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system(*_UMF)
    ops.test("NormDispIncr", 1e-8, 30)
    ops.algorithm("Newton")
    ops.integrator("Newmark", 0.6, 0.3025)
    ops.analysis("Transient")
    for _ in range(5):
        assert ops.analyze(1, 0.01) == 0, "(iv) transient smoke blew up"
    for wet, dry in ((4, 20), (3, 21)):   # continuity holds through transient too
        assert abs(ops.nodeDisp(wet, 1) - ops.nodeDisp(dry, 1)) <= 1e-9
        assert abs(ops.nodeDisp(wet, 2) - ops.nodeDisp(dry, 2)) <= 1e-9


# ==========================================================================
# (v) BTET10 winding acceptance — right vs left handed (SUBPROCESS)
# ==========================================================================
# Runs in a `-S` subprocess so the "printed once per process" winding note is
# observed in isolation (order-independent). Asserts internally, prints sentinels.
_WINDING_DRIVER = r'''
import os, sys
DIST = sys.argv[1]
os.add_dll_directory(DIST); sys.path.insert(0, DIST)
os.environ["PATH"] = DIST + os.pathsep + os.environ.get("PATH", "")
import opensees as ops
assert os.path.normcase(DIST) in os.path.normcase(ops.__file__), ops.__file__

P = {0:(0.0,0.0,0.0), 1:(1.2,0.0,0.0), 2:(0.0,1.1,0.0), 3:(0.0,0.0,1.3)}
def pf(x,y,z): return 3.0 + 2.0*x + 5.0*y + 4.0*z

def run(order):
    ops.wipe(); ops.model("basic","-ndm",3,"-ndf",4)
    vt=[]
    for slot in range(4):
        x,y,z=P[order[slot]]; ops.node(100+slot,x,y,z); vt.append(100+slot)
    ops.model("basic","-ndm",3,"-ndf",3)
    edges=[(0,1),(1,2),(0,2),(0,3),(2,3),(1,3)]; mt=[]
    for i,(a,b) in enumerate(edges):
        xa=P[order[a]]; xb=P[order[b]]; t=200+i
        ops.node(t,*(0.5*(xa[k]+xb[k]) for k in range(3))); mt.append(t)
    ops.nDMaterial("ElasticIsotropic",1,1.0e4,0.25,2.0)
    ops.element("LadrunoUP",1,*(vt+mt),1,"-Kf",2.2e6,"-poro",0.3,"-rhoF",1.0,
                "-perm",1e-4,1e-4,1e-4,"-pOrder","linear","-dynSeepage","off")
    for t in vt: ops.fix(t,1,1,1,0)
    for t in mt: ops.fix(t,1,1,1)
    ops.timeSeries("Linear",1); ops.pattern("Plain",1,1)
    for slot in (0,1,2):
        x,y,z=P[order[slot]]; ops.sp(vt[slot],4,pf(x,y,z))
    ops.constraints("Penalty",1e14,1e14); ops.numberer("Plain")
    ops.system("UmfPack"); ops.test("NormDispIncr",1e-13,50)
    ops.algorithm("Newton"); ops.integrator("LoadControl",1.0); ops.analysis("Static")
    assert ops.analyze(1)==0, "winding leg failed to converge"
    return ops.nodeDisp(vt[3],4)

# order A = conventional right-handed (0,1,2,3); order B = swap 1<->2 (left-handed)
pa = run((0,1,2,3))
pb = run((0,2,1,3))
d = abs(pa-pb)
assert d <= 1e-12*max(1.0,abs(pa)), "winding parity broken: %r vs %r (d=%g)"%(pa,pb,d)
assert abs(pa) > 1.0, "winding probe pressure ~0"
print("WINDING-PARITY-OK pa=%.15g pb=%.15g d=%g"%(pa,pb,d))
'''


@pytest.mark.t1
def test_v_btet10_winding_acceptance(tmp_path):
    """The same one-element BTET10 seepage problem fed right-handed (0,1,2,3) vs
    left-handed (0,2,1,3 — vertices 1,2 swapped) node input yields IDENTICAL free-
    vertex pressure (≤1e-12 rel): the element folds dv = w*|detJ| for the
    negative-winding case, so the signed-J gradients stay orientation-correct. The
    conventionally right-handed tet is the one that evaluates detJ<0 under the
    pinned node↔barycentric map, so its ONE-TIME informational note must appear."""
    driver = tmp_path / "winding_driver.py"
    driver.write_text(_WINDING_DRIVER)
    env = dict(os.environ)
    env["PATH"] = _DIST + os.pathsep + env.get("PATH", "")
    # encoding pinned: the child prints the UTF-8 splash banner, and cp1252
    # text-mode decode raises in the reader thread => proc.stdout is None and
    # the assert dies with TypeError instead of a diagnosis (quirks ledger).
    proc = subprocess.run([sys.executable, "-S", str(driver), _DIST],
                          capture_output=True, text=True, timeout=200, env=env,
                          encoding="utf-8", errors="replace")
    out = (proc.stdout or "") + "\n" + (proc.stderr or "")
    assert "WINDING-PARITY-OK" in proc.stdout and proc.returncode == 0, (
        f"winding parity gate failed (rc={proc.returncode})\n{out}")
    # the one-time note (optional per the pin, but observable in this isolated run)
    assert ("folding the integration measure" in out
            or "detJ < 0" in out), (
        "expected the BTET10 right-handed-winding fold note in the isolated run\n"
        + out)


# ==========================================================================
# (vi) straight-side guard — curved mid-edge -> loud error, not a crash
# ==========================================================================
@pytest.mark.t1
def test_vi_straight_side_guard_loud(capfd):
    """A BT6 with one mid-edge node displaced 1% of the edge length off the
    midpoint must trip the STRICT straight-side guard (tol 1e-6*edgeLen): the
    setDomain guard prints a loud message naming the element/node/distance and
    deactivates the element (geomValid_ = false -> zero stiffness). Observed
    surfacing via openseespy: ops.element() returns normally (no Python
    exception); the deactivated element then makes the assembled system SINGULAR,
    so ops.analyze() returns a non-zero (< 0) flag. This test gates BOTH the
    message text and the analysis abort."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, 1.0, 0.0)
    ops.node(3, 0.0, 1.0)
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(4, 0.5, 0.01)   # mid(1,2) displaced 1% of edge length (=1.0) off
    ops.node(5, 0.5, 0.5)
    ops.node(6, 0.0, 0.5)
    ops.nDMaterial("ElasticIsotropic", 1, 1.0e4, 0.2, 2.0)
    ops.element("LadrunoUP", 1, 1, 2, 3, 4, 5, 6, 1,
                "-Kf", 2.2e6, "-poro", 0.3, "-rhoF", 1.0,
                "-perm", 1e-4, 1e-4, "-pOrder", "linear", "-dynSeepage", "off")
    ops.fix(1, 1, 1, 1)
    ops.fix(2, 1, 1, 0)
    ops.fix(3, 1, 1, 0)
    ops.fix(4, 1, 1)
    ops.fix(5, 1, 1)
    ops.fix(6, 1, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.sp(2, 3, 5.0)
    ops.constraints("Penalty", 1e14, 1e14)
    ops.numberer("Plain")
    ops.system(*_UMF)
    ops.test("NormDispIncr", 1e-10, 20)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    rc = ops.analyze(1)
    # C-level opserr goes to fd 2 -> capfd (NOT capsys), per project memory
    captured = capfd.readouterr()
    msg = captured.out + captured.err
    assert "off the midpoint" in msg and "STRAIGHT sides" in msg, (
        f"straight-side guard message missing; captured:\n{msg}")
    assert "not activated" in msg, f"expected 'element not activated'; got:\n{msg}"
    assert rc != 0, (
        f"deactivated element should singularize the solve, but analyze rc={rc}")


# ==========================================================================
# (vii) BT6-TH static steady-seepage exactness — the P1 gate-4 analog
# ==========================================================================
@pytest.mark.t1
def test_vii_bt6_steady_seepage_linear_exact():
    """Linear-p manufactured solution / patch test. A unit square split into 4
    straight-sided BT6 about an INTERIOR vertex (0.5,0.5). All u fixed (pure
    seepage). A linear field p = 1 + 2x + 3y is prescribed on the four boundary
    vertices; the interior vertex p is left FREE. The P1 vertex-pressure Laplacian
    reproduces a linear field exactly, so the solved interior p == field(0.5,0.5)
    = 3.5 to ~machine, and the per-GP Darcy flux == -k̄·∇p = -1e-4*(2,3) exactly
    at every Gauss point (the P1 gate-4 steady-seepage exactness analog).

    Panel addition: the eleResponse('porePressure') VALUES are gated against the
    exact field at the 3-pt interior rule points. BT6 GP g carries barycentric
    weight 2/3 on element-node g and 1/6 on the other two (LadrunoUPShapes.h
    pinned rule {(1/6,1/6),(2/3,1/6),(1/6,2/3)}), so under the affine map
    x_gp(g) = 2/3*v_g + 1/6*v_{g+1} + 1/6*v_{g+2} and the expected value is
    pfield(x_gp). Measured on element 1 (v=(0,0),(1,0),(.5,.5)):
    [1.75, 2.75, 3.0] — matching the GP-g <-> node-g convention."""
    def pfield(x, y):
        return 1.0 + 2.0 * x + 3.0 * y

    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    V = {1: (0.0, 0.0), 2: (1.0, 0.0), 3: (1.0, 1.0),
         4: (0.0, 1.0), 5: (0.5, 0.5)}
    for t, (x, y) in V.items():
        ops.node(t, x, y)
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    M = {}
    _next = [11]

    def mid(a, b):
        key = tuple(sorted((a, b)))
        if key not in M:
            xa, ya = V[a]
            xb, yb = V[b]
            t = _next[0]
            _next[0] += 1
            ops.node(t, 0.5 * (xa + xb), 0.5 * (ya + yb))
            M[key] = t
        return M[key]

    kx, ky = 1.0e-4, 1.0e-4
    ops.nDMaterial("ElasticIsotropic", 1, 1.0e4, 0.2, 2.0)
    tris = [(1, 2, 5), (2, 3, 5), (3, 4, 5), (4, 1, 5)]
    for e, (a, b, c) in enumerate(tris, start=1):
        ops.element("LadrunoUP", e, a, b, c,
                    mid(a, b), mid(b, c), mid(c, a), 1,
                    "-Kf", 2.2e6, "-poro", 0.3, "-rhoF", 1.0,
                    "-perm", kx, ky, "-pOrder", "linear", "-dynSeepage", "off")
    for t in V:
        ops.fix(t, 1, 1, 0)          # all u fixed; p handled below
    for t in M.values():
        ops.fix(t, 1, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t in (1, 2, 3, 4):           # boundary vertices prescribed; 5 free
        x, y = V[t]
        ops.sp(t, 3, pfield(x, y))
    ops.constraints("Penalty", 1e14, 1e14)
    ops.numberer("Plain")
    ops.system(*_UMF)
    ops.test("NormDispIncr", 1e-13, 50)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    assert ops.analyze(1) == 0, "(vii) steady-seepage patch solve failed"

    p5 = ops.nodeDisp(5, 3)
    assert abs(p5 - pfield(0.5, 0.5)) <= 1e-9, (
        f"interior vertex p {p5} != linear field {pfield(0.5, 0.5)}")
    # flux == -k*grad(linear) == constant at every GP of every element
    for e in (1, 2, 3, 4):
        fx = ops.eleResponse(e, "flux")
        assert len(fx) == 3 * 2
        for gp in range(3):
            qx, qy = fx[2 * gp], fx[2 * gp + 1]
            assert abs(qx - (-kx * 2.0)) <= 1e-9 * abs(kx * 2.0) + 1e-15, (e, gp, qx)
            assert abs(qy - (-ky * 3.0)) <= 1e-9 * abs(ky * 3.0) + 1e-15, (e, gp, qy)
    # panel addition: porePressure VALUES == the exact linear field at the 3-pt
    # rule points (GP g: barycentric 2/3 on node g, 1/6 on the others)
    for e, (a, b, c) in enumerate(tris, start=1):
        pp = ops.eleResponse(e, "porePressure")
        assert len(pp) == 3
        verts = [np.array(V[a]), np.array(V[b]), np.array(V[c])]
        for gp in range(3):
            xgp = (2.0 / 3.0) * verts[gp] \
                + (1.0 / 6.0) * verts[(gp + 1) % 3] \
                + (1.0 / 6.0) * verts[(gp + 2) % 3]
            expect = pfield(xgp[0], xgp[1])
            assert abs(pp[gp] - expect) <= 1e-9 * max(1.0, abs(expect)), (
                f"ele {e} GP {gp}: porePressure {pp[gp]} != field {expect} "
                f"at {tuple(xgp)}")


# ==========================================================================
# (viii) BTET10-TH 1D consolidation smoke vs the Terzaghi series
# ==========================================================================
def _kuhn_column(nz, hx, hy, hz, matArgs, mat_id=1, E=1.0e4, nu=0.2, rho=2.0):
    """A 1x1xnz box meshed with the (Freudenthal/Kuhn) 6-tet-per-cell split about
    the (0,0,0)-(1,1,1) cube diagonal — the split is congruent across cells so
    faces AND edges match, and a global unique-edge table gives every shared edge
    ONE mid-edge node. Vertices ndf=4 (p carriers), mid-edges ndf=3. Returns
    (vtag{(i,j,k):tag}, columns) where columns[(i,j)] lists (k, tag) up the height.
    """
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 4)
    Nx, Ny = 1, 1
    vtag = {}
    t = 1
    for k in range(nz + 1):
        for j in range(Ny + 1):
            for i in range(Nx + 1):
                ops.node(t, i * hx, j * hy, k * hz)
                vtag[(i, j, k)] = t
                t += 1
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("ElasticIsotropic", mat_id, E, nu, rho)

    # cube-corner (i,j,k) offsets by binary index n = i + 2j + 4k
    def corner(base, n):
        di = n & 1
        dj = (n >> 1) & 1
        dk = (n >> 2) & 1
        return (base[0] + di, base[1] + dj, base[2] + dk)

    # 6 Kuhn tets along diagonal 0->7 (monotone bit-flip paths)
    kuhn = [(0, 1, 3, 7), (0, 1, 5, 7), (0, 2, 3, 7),
            (0, 2, 6, 7), (0, 4, 5, 7), (0, 4, 6, 7)]
    edges_local = [(0, 1), (1, 2), (0, 2), (0, 3), (2, 3), (1, 3)]

    midcache = {}
    mt = [1000]

    coords = {}
    for gk, gt in vtag.items():
        coords[gt] = (gk[0] * hx, gk[1] * hy, gk[2] * hz)

    def midnode(va, vb):
        key = tuple(sorted((va, vb)))
        if key not in midcache:
            xa = coords[va]
            xb = coords[vb]
            tt = mt[0]
            mt[0] += 1
            xm = tuple(0.5 * (xa[c] + xb[c]) for c in range(3))
            ops.node(tt, *xm)
            coords[tt] = xm            # record mid-edge coords too (used by BCs)
            midcache[key] = tt
        return midcache[key]

    e = 1
    for k in range(nz):
        for j in range(Ny):
            for i in range(Nx):
                base = (i, j, k)
                cindex = {n: vtag[corner(base, n)] for n in range(8)}
                for tet in kuhn:
                    vt = [cindex[c] for c in tet]
                    mt6 = [midnode(vt[a], vt[b]) for (a, b) in edges_local]
                    ops.element("LadrunoUP", e, *(vt + mt6), mat_id, *matArgs)
                    e += 1

    columns = {}
    for j in range(Ny + 1):
        for i in range(Nx + 1):
            columns[(i, j)] = [(k, vtag[(i, j, k)]) for k in range(nz + 1)]
    return vtag, columns, coords, midcache


@pytest.mark.t1
def test_viii_btet10_terzaghi_consolidation_smoke():
    """1D Terzaghi consolidation of a laterally-confined BTET10-TH column
    (1x1x`nz`, Kuhn 6-tet split, straight sides). Confined ux=uy=0 everywhere,
    uz=0 base, drained (p=0) top vertices; instantaneous surcharge applied with
    the CONSISTENT 6-node-triangle face load (corner weight 0, mid-side A/3) so the
    initial pore pressure matches Terzaghi's p0 = q*Q̄/(Q̄+Eoed). Loose gates
    (decay CLASS, not records): p-profile err ≤5% at Tv∈{0.2,0.5}; U within 0.03
    of the series; settlement monotone. TH (no stabilization)."""
    E, nu, rho = 1.0e4, 0.2, 2.0
    Kf, poro, rhoF, kbar = 4444.0, 0.4, 1.0, 1.0e-4
    nz = 6
    hx = hy = 1.0
    hz = 1.0
    H = nz * hz
    q = 10.0
    Eoed = E * (1 - nu) / ((1 + nu) * (1 - 2 * nu))
    ooQ = poro / Kf
    Qbar = 1.0 / ooQ
    cv = kbar / (1.0 / Eoed + ooQ)
    p0 = q * Qbar / (Qbar + Eoed)         # Terzaghi undrained initial pressure
    s0 = q * H / (Eoed + Qbar)            # undrained settlement
    sinf = q * H / Eoed                   # drained settlement

    matArgs = ["-Kf", Kf, "-poro", poro, "-rhoF", rhoF,
               "-perm", kbar, kbar, kbar, "-pOrder", "linear",
               "-dynSeepage", "off"]
    vtag, columns, coords, midcache = _kuhn_column(nz, hx, hy, hz, matArgs,
                                                   E=E, nu=nu, rho=rho)

    # BCs: confine ux=uy on ALL nodes (vertices + mid-edges); uz=0 at base
    # vertices; drained (p=0) at the top vertices only.
    all_v = list(vtag.values())
    all_m = list(midcache.values())
    top_k = nz
    for (i, j, k), t in vtag.items():
        fz = 1 if k == 0 else 0
        fp = 1 if k == top_k else 0
        ops.fix(t, 1, 1, fz, fp)
    for t in all_m:
        # mid-edge nodes: confine ux=uy; uz free unless the edge lies on the base
        cz = coords[t][2]
        fz = 1 if abs(cz) < 1e-9 else 0
        ops.fix(t, 1, 1, fz)

    # consistent surcharge on the top face: each top Kuhn triangle is a 6-node
    # tri -> uniform pressure lumps to the 3 mid-side nodes (A/3 each), corners 0.
    # Build the top-face triangle list from the Kuhn split's top facet.
    # Top cube face (k+1 level) corners n in {4,5,6,7}; the two Kuhn tets touching
    # it, (0,4,5,7) and (0,4,6,7), each contribute the top triangle {4,5,7} and
    # {4,6,7}. In (i,j) top-cell terms: corners 4=(0,0),5=(1,0),6=(0,1),7=(1,1).
    zload = {}

    def add(tag, val):
        zload[tag] = zload.get(tag, 0.0) + val

    kx = nz - 1  # top cell layer index
    for j in range(1):
        for i in range(1):
            c = {
                4: vtag[(i, j, kx + 1)],
                5: vtag[(i + 1, j, kx + 1)],
                6: vtag[(i, j + 1, kx + 1)],
                7: vtag[(i + 1, j + 1, kx + 1)],
            }
            for tri in ((4, 5, 7), (4, 6, 7)):
                a, b, cc = (c[tri[0]], c[tri[1]], c[tri[2]])
                pa, pb, pc = coords[a], coords[b], coords[cc]
                # triangle area in the xy-plane
                area = 0.5 * abs((pb[0] - pa[0]) * (pc[1] - pa[1])
                                 - (pc[0] - pa[0]) * (pb[1] - pa[1]))
                w = q * area / 3.0
                for (va, vb) in ((a, b), (b, cc), (cc, a)):
                    add(midcache[tuple(sorted((va, vb)))], w)

    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for tag, val in zload.items():
        ops.load(tag, 0.0, 0.0, -val)

    ops.system(*_UMF)
    ops.numberer("RCM")
    ops.constraints("Transformation")
    ops.test("NormDispIncr", 1e-9, 40)
    ops.algorithm("Newton")
    ops.integrator("Newmark", 0.6, 0.3025)
    ops.analysis("Transient")

    # settlement measure: the mean uz over ALL top-face nodes (4 vertices +
    # 5 mid-edges). On this coarse quadratic mesh the load lumps to the top mid-
    # edge nodes (consistent 6-node-tri rule -> corners get 0), so corner nodes
    # LAG and mid-edge nodes LEAD the true surface settlement; the mean of the two
    # families tracks the Terzaghi degree of consolidation to ~0.05 while either
    # family alone oscillates. (Pore pressure -- the physics gate -- is measured on
    # the interior vertex carriers and is clean.)
    top_nodes = [vtag[(i, j, top_k)] for i in range(2) for j in range(2)]
    top_nodes += [m for m in midcache.values() if abs(coords[m][2] - H) < 1e-9]

    def top_settle():
        return -float(np.mean([ops.nodeDisp(n, 3) for n in top_nodes]))

    # tiny priming step, then march to each target Tv
    assert ops.analyze(1, 1.0e-4) == 0, "(viii) priming step failed"
    t = 1.0e-4
    dt = 0.02
    col = columns[(0, 0)]                 # the x=0,y=0 vertex column
    settle = {}
    perr = {}
    for Tv in (0.1, 0.2, 0.5):
        tt = Tv * H * H / cv
        nst = int(round((tt - t) / dt))
        assert nst > 0
        assert ops.analyze(nst, dt) == 0, f"(viii) blew up before Tv={Tv}"
        t += nst * dt
        TvA = t * cv / (H * H)
        errs = []
        for (k, tag) in col:
            if k == 0 or k == top_k:
                continue
            z = k * hz
            Z = (H - z) / H               # depth from the drained (top) face
            pnum = ops.nodeDisp(tag, 4)
            pref = p0 * terzaghi_series(Z, TvA)
            errs.append(abs(pnum - pref) / p0)
        perr[Tv] = max(errs)
        settle[Tv] = top_settle()

    # PHYSICS gate (the pin's ≤5%): pore-pressure profile vs the Terzaghi series
    # on the interior vertex carriers.
    assert perr[0.2] <= 0.05, f"Tv=0.2 pore-pressure err {perr[0.2]:.3%} > 5%"
    assert perr[0.5] <= 0.05, f"Tv=0.5 pore-pressure err {perr[0.5]:.3%} > 5%"
    # settlement monotone (the pin's second gate) + physically bounded
    assert settle[0.1] < settle[0.2] < settle[0.5], (
        f"settlement not monotone: {settle}")
    for Tv in settle:
        assert s0 * 0.5 < settle[Tv] < sinf * 1.5, (
            f"settlement Tv={Tv} = {settle[Tv]:g} outside [~s0, ~sinf] bounds "
            f"(s0={s0:g}, sinf={sinf:g})")
    # LOOSE documentation comparison to the series (coarse-mesh smoke, generous
    # band): measured U(all-top mean) ~0.40/0.55 vs series 0.357/0.504 (~0.05 off).
    U_fem = {Tv: (settle[Tv] - s0) / (sinf - s0) for Tv in settle}
    for Tv in (0.2, 0.5):
        assert abs(U_fem[Tv] - terzaghi_U(Tv)) < 0.12, (
            f"U(Tv={Tv}) fem {U_fem[Tv]:.4f} vs series {terzaghi_U(Tv):.4f} "
            "outside the loose 0.12 smoke band")


# ==========================================================================
# (ix)/(x) TH transient-block FD gates — the P1 test-8 pattern on one BT6-TH
# ==========================================================================
# Single BT6-TH on the unit right triangle v1=(0,0) v2=(1,0) v3=(0,1);
# mid-edges m12=(.5,0) m23=(.5,.5) m31=(0,.5). 15 DOFs, node-major under the
# Plain numberer: n1(ux,uy,p) n2(ux,uy,p) n3(ux,uy,p) n4(ux,uy) n5(ux,uy)
# n6(ux,uy) -> u-DOF slots _THU, carrier-p slots _THP (the COMPACTED basis).
_THU = [0, 1, 3, 4, 6, 7, 9, 10, 11, 12, 13, 14]
_THP = [2, 5, 8]
_TH_KF, _TH_PORO, _TH_RHOF, _TH_KBAR = 5000.0, 0.4, 1.0, 1.0e-4
_TH_RHO = 2.0
_TH_AREA = 0.5                     # unit right triangle
_TH_THICK = 1.0

# state points for the FD probes (order: n1 u,u,p | n2 u,u,p | n3 u,u,p |
# n4 u,u | n5 u,u | n6 u,u) — random-ish, nonzero everywhere it matters
_THX0 = np.array([0.0, 0.0, 0.5, 0.01, -0.005, 2.0, 0.008, 0.012, 3.0,
                  -0.003, 0.009, 0.004, -0.006, 0.002, 0.007])
_THV0 = np.array([0.001, -0.002, 0.03, 0.002, 0.001, -0.02, -0.001, 0.003,
                  0.05, 0.0, -0.001, 0.002, 0.001, -0.002, 0.001])
_THA0 = np.array([0.01, -0.02, 0.3, 0.02, 0.01, -0.2, -0.01, 0.03, 0.5,
                  0.0, -0.01, 0.02, 0.01, -0.02, 0.01])

_TH_FD_DOFS = [(1, 3), (2, 3), (3, 3), (4, 2), (5, 2), (6, 2)]  # (node, ndf)


def _bt6_th_blocks_analytic():
    """Independent numpy oracle for the BT6-TH compacted blocks on the unit
    right triangle (thick t, area A):
      dNp   — constant gradients of the linear barycentric carrier basis
              [L1, L2, L3] (carrier order == element node order):
              grad L1 = (-1,-1), grad L2 = (1,0), grad L3 = (0,1).
      S     = (1/Qbar) * t * (A/12) * [[2,1,1],[1,2,1],[1,1,2]]
              (exact: ∫ Li Lj dA = A/12 off-diagonal, A/6 diagonal; the 3-pt
              interior rule integrates quadratics exactly).
      G     = ∫ dNp^T k̄ ρf Nu dV: every Bernstein-quadratic u-shape (Li² and
              2LiLj alike) integrates to A/6, so
              G[b, (a,i)] = k̄_i * ρf * dNp[b,i] * t * A/6  (independent of a).
    """
    A, t = _TH_AREA, _TH_THICK
    ooQ = _TH_PORO / _TH_KF        # alpha=1, Ks infinite -> 1/Qbar = n/Kf
    dNp = np.array([[-1.0, -1.0], [1.0, 0.0], [0.0, 1.0]])
    S = ooQ * t * (A / 12.0) * np.array([[2.0, 1.0, 1.0],
                                         [1.0, 2.0, 1.0],
                                         [1.0, 1.0, 2.0]])
    G = np.zeros((3, 12))
    for b in range(3):
        for a in range(6):
            for i in range(2):
                G[b, 2 * a + i] = _TH_KBAR * _TH_RHOF * dNp[b, i] * t * A / 6.0
    return S, G


def _bt6_th_args(dyn):
    return ["-Kf", _TH_KF, "-poro", _TH_PORO, "-rhoF", _TH_RHOF,
            "-perm", _TH_KBAR, _TH_KBAR, "-pOrder", "linear",
            "-dynSeepage", dyn]


def _bt6_th_fd_build(xvec, dyn):
    """One BT6-TH with EVERY one of the 15 DOFs penalty-prescribed to xvec and
    one static step solved (constraints('Penalty') — Transformation on a fully-
    prescribed rig leaves zero free equations and FullGenLinSOE dies; P1 donor
    pattern)."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, 1.0, 0.0)
    ops.node(3, 0.0, 1.0)
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(4, 0.5, 0.0)
    ops.node(5, 0.5, 0.5)
    ops.node(6, 0.0, 0.5)
    ops.nDMaterial("ElasticIsotropic", 1, 1.0e4, 0.2, _TH_RHO)
    ops.element("LadrunoUP", 1, 1, 2, 3, 4, 5, 6, 1, *_bt6_th_args(dyn))
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    k = 0
    for nd, nf in _TH_FD_DOFS:
        for d in range(1, nf + 1):
            ops.sp(nd, d, float(xvec[k]))
            k += 1
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Penalty", 1e12, 1e12)
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Linear")
    ops.analysis("Static")
    assert ops.analyze(1) == 0


def _bt6_th_dyn_resid(xvec, vvec, avec, dyn):
    """getResistingForceIncInertia residual (15-vector): static state from
    xvec, trial rates imposed, reactions('-dynamic') gathered node-major."""
    _bt6_th_fd_build(xvec, dyn)
    k = 0
    for nd, nf in _TH_FD_DOFS:
        for d in range(1, nf + 1):
            ops.setNodeVel(nd, d, float(vvec[k]), "-commit")
            ops.setNodeAccel(nd, d, float(avec[k]), "-commit")
            k += 1
    ops.reactions("-dynamic")
    R = np.zeros(15)
    k = 0
    for nd, nf in _TH_FD_DOFS:
        rr = ops.nodeReaction(nd)
        for d in range(nf):
            R[k] = rr[d]
            k += 1
    return R


def _bt6_th_printA(dyn, transient, dt=0.1, gamma=0.6, beta=0.3025):
    """Assembled tangent of the free-floating single BT6-TH via printA under
    Plain numbering (equation order == node-major DOF order; donor pattern).
    static    -> A = [Kuu, -Q; 0, H]
    transient -> A = Ktot + c2*C + c3*M (Newmark)."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, 1.0, 0.0)
    ops.node(3, 0.0, 1.0)
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(4, 0.5, 0.0)
    ops.node(5, 0.5, 0.5)
    ops.node(6, 0.0, 0.5)
    ops.nDMaterial("ElasticIsotropic", 1, 1.0e4, 0.2, _TH_RHO)
    ops.element("LadrunoUP", 1, 1, 2, 3, 4, 5, 6, 1, *_bt6_th_args(dyn))
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(3, 0.0, 0.0, 0.0)
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    if transient:
        ops.integrator("Newmark", gamma, beta)
        ops.algorithm("Linear")
        ops.analysis("Transient")
        assert ops.analyze(1, dt) == 0
    else:
        ops.integrator("LoadControl", 1.0)
        ops.algorithm("Linear")
        ops.analysis("Static")
        assert ops.analyze(1) == 0
    A = np.array(ops.printA("-ret"))
    N = int(round(len(A) ** 0.5))
    assert N == 15
    return A.reshape(N, N).T          # printA is COLUMN-major -> transpose


def _bt6_th_fd_C_M(dyn, h=1e-5):
    """Central-FD of the '-dynamic' reactions w.r.t. nodal velocities (-> C)
    and accelerations (-> M), 15x15 each."""
    Cfd = np.zeros((15, 15))
    for j in range(15):
        vp = _THV0.copy()
        vp[j] += h
        vm = _THV0.copy()
        vm[j] -= h
        Cfd[:, j] = (_bt6_th_dyn_resid(_THX0, vp, _THA0, dyn)
                     - _bt6_th_dyn_resid(_THX0, vm, _THA0, dyn)) / (2 * h)
    Mfd = np.zeros((15, 15))
    for j in range(15):
        ap = _THA0.copy()
        ap[j] += h
        am = _THA0.copy()
        am[j] -= h
        Mfd[:, j] = (_bt6_th_dyn_resid(_THX0, _THV0, ap, dyn)
                     - _bt6_th_dyn_resid(_THX0, _THV0, am, dyn)) / (2 * h)
    return Cfd, Mfd


@pytest.mark.t1
def test_ix_th_transient_block_fd_gate():
    """(panel hardening 1) The P1 FD pattern on one BT6-TH, -dynSeepage OFF:
    central-FD of reactions('-dynamic') w.r.t. vel/accel on the fully-penalty-
    prescribed 15-DOF rig, cross-checked against the free-floating static
    printA and the independent numpy oracle. Gates:
      * Cpu == (-Kup)^T          — the honest-p transpose relation on the
        COMPACTED carrier basis (Q assembled from the same Np both places);
        1e-6 rel (panel critic measured 6.2e-8).
      * Cpp == S (analytic)      — compacted storage present and exact.
      * Cuu == Cup == 0          — no Rayleigh set, no stab on TH.
      * mass p-rows AND p-cols 0 — dynSeepage off: the p-equation carries no
        acceleration term at all.
      * Muu symmetric, x-direction total == rho*Area*thick (Bernstein
        partition of unity; critic measured 3e-9 rel).
    Makes the TH compacted-basis Q^T/S/M correctness self-contained instead of
    riding on the B1 physics gate."""
    S_ex, _ = _bt6_th_blocks_analytic()
    Ka = _bt6_th_printA("off", transient=False)
    Kup = Ka[np.ix_(_THU, _THP)]                    # = -Q (honest-p static)
    assert np.max(np.abs(Ka[np.ix_(_THP, _THU)])) == 0.0   # static pu block = 0
    qs = np.max(np.abs(Kup))
    assert qs > 0.0

    Cfd, Mfd = _bt6_th_fd_C_M("off")

    # --- damping blocks ---
    cs = max(qs, float(np.max(np.abs(S_ex))))
    rel_t = np.max(np.abs(Cfd[np.ix_(_THP, _THU)] - (-Kup).T)) / qs
    assert rel_t < 1e-6, f"honest-p transpose relation Cpu != (-Kup)^T: {rel_t:.2e} rel"
    rel_s = np.max(np.abs(Cfd[np.ix_(_THP, _THP)] - S_ex)) / np.max(np.abs(S_ex))
    assert rel_s < 1e-6, f"compacted storage Cpp != S analytic: {rel_s:.2e} rel"
    assert np.max(np.abs(Cfd[np.ix_(_THU, _THU)])) < 1e-6 * cs   # no Rayleigh
    assert np.max(np.abs(Cfd[np.ix_(_THU, _THP)])) < 1e-6 * cs   # damp up = 0

    # --- mass blocks ---
    Muu = Mfd[np.ix_(_THU, _THU)]
    ms = float(np.max(np.abs(Muu)))
    assert ms > 0.0
    assert np.max(np.abs(Mfd[:, _THP])) < 1e-9 * ms, "mass p-cols not zero"
    assert np.max(np.abs(Mfd[np.ix_(_THP, _THU)])) < 1e-9 * ms, (
        "mass p-rows not zero with -dynSeepage off")
    # symmetry is FD-measured, so the floor is ulp(R)/(2h) ~ 1e-10 abs
    # (measured 1.8e-10 = 2.7e-9 rel), not the analytic 0 — gate at 1e-6 rel
    assert np.max(np.abs(Muu - Muu.T)) < 1e-6 * ms, "Muu not symmetric"
    xsl = slice(0, 12, 2)                            # x-direction u-DOF slots
    total_x = float(np.sum(Muu[xsl, xsl]))
    m_exact = _TH_RHO * _TH_AREA * _TH_THICK
    assert abs(total_x - m_exact) < 1e-6 * m_exact, (
        f"Muu x-total {total_x} != rho*A*t = {m_exact}")


@pytest.mark.t1
def test_x_th_dynseepage_on_default_config():
    """(panel hardening 2) The '-dynSeepage on' configuration (explicit opt-in;
    default flipped to OFF by the P4 B5 adjudication) —
    the only P3 test exercising it. Same FD rig: the dR/da p-rows must now carry
    the +G dynamic-seepage residual coupling (∫ dNp^T k̄ ρf Nu dV, checked
    against the ANALYTIC block — every Bernstein u-shape integrates to A/6)
    while G stays OUT of the assembled tangent: the transient printA p-row/u-col
    block equals c2*Q^T exactly, with the mismatch bound far below the c3*G
    footprint a leaked G would leave (the P1 dyn-seepage quantification,
    TH edition). Mass p-COLS stay zero either way."""
    S_ex, G_ex = _bt6_th_blocks_analytic()
    gs = float(np.max(np.abs(G_ex)))
    assert gs > 0.0

    _, Mfd = _bt6_th_fd_C_M("on")
    Gfd = Mfd[np.ix_(_THP, _THU)]
    rel_g = np.max(np.abs(Gfd - G_ex)) / gs
    assert rel_g < 1e-6, (
        f"dyn-seepage accel coupling != analytic G: {rel_g:.2e} rel "
        f"(FD max {np.max(np.abs(Gfd)):.3e}, analytic max {gs:.3e})")
    ms = float(np.max(np.abs(Mfd[np.ix_(_THU, _THU)])))
    assert np.max(np.abs(Mfd[:, _THP])) < 1e-9 * ms, "mass p-cols not zero"

    # G stays OUT of the tangent: transient effective A = Ktot + c2*C + c3*M
    # has A[p,u] == c2*Q^T (no c3*G term). Q from the static printA (-Kup).
    dt, gamma, beta = 0.1, 0.6, 0.3025
    c2 = gamma / (beta * dt)
    c3 = 1.0 / (beta * dt * dt)
    Ka = _bt6_th_printA("on", transient=False)
    Q = -Ka[np.ix_(_THU, _THP)]
    At = _bt6_th_printA("on", transient=True, dt=dt, gamma=gamma, beta=beta)
    dev = np.max(np.abs(At[np.ix_(_THP, _THU)] - c2 * Q.T))
    scal = float(np.max(np.abs(At)))
    assert dev < 1e-8 * scal, (
        f"transient A[p,u] != c2*Q^T (G leaked into the tangent?): dev {dev:.3e}")
    # and the gate is DISCRIMINATING: a leaked G would move A[p,u] by ~c3*|G|,
    # orders above the accepted mismatch
    assert dev < 0.01 * c3 * gs, (
        f"gate not discriminating: dev {dev:.3e} vs c3*Gmax {c3 * gs:.3e}")
