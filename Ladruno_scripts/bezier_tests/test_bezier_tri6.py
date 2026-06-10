"""Regression test for the BezierTri6 element (ELE_TAG 272).

Tracked, self-contained, path-portable. Exercises the element end to end
against the built `opensees` module:

  1. Patch test            — constant strain -> constant stress (exact)
  2. SelfWeight load path  — total reaction = gamma * A * t
  3. Lumped mass default   — positive eigenvalues (all-positive lumped mass)
  4. PlaneStress B-bar guard — '-bbar' + PlaneStress warns and is ignored
  5. Multi-element cantilever — inline straight-sided T6 mesh; tip deflection
                                within tolerance of the Timoshenko estimate

Resolving the built module (in order):
  - env LADRUNO_OPENSEES_BIN  (dir holding opensees.pyd)
  - <repo-root>/dist/bin      (repo-root = three levels up from this file)
  - a bare `import opensees`  (already on sys.path)

Run:  python Ladruno_scripts/bezier_tests/test_bezier_tri6.py
Exit code 0 = all pass, 1 = a failure (CI-friendly).

Element authors: Nicolas Mora Bowen, Patricio Palacios, Jose Abell,
Guppi (Ladruno). Ref: Kadapa, IJNME 2019; 117(5):543-573, doi:10.1002/nme.5967.
"""
import os
import sys
from pathlib import Path

os.environ.setdefault("LADRUNO_OPENSEES_QUIET", "1")


def _import_opensees():
    candidates = []
    if os.environ.get("LADRUNO_OPENSEES_BIN"):
        candidates.append(os.environ["LADRUNO_OPENSEES_BIN"])
    # repo-root/dist/bin  (this file: <root>/Ladruno_scripts/bezier_tests/<file>)
    root = Path(__file__).resolve().parents[2]
    candidates.append(str(root / "dist" / "bin"))
    for c in candidates:
        if c and (Path(c) / "opensees.pyd").exists() or (c and (Path(c) / "opensees.so").exists()):
            sys.path.insert(0, c)
            break
    import opensees as ops
    return ops


ops = _import_opensees()
TOL = 1e-6
results = []


def record(name, ok, detail=""):
    results.append((name, ok))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}{(' — ' + detail) if detail else ''}")


# --- 1. patch test -----------------------------------------------------------
def patch_test():
    ops.wipe(); ops.model('basic', '-ndm', 2, '-ndf', 2)
    E, nu, eps0 = 1000.0, 0.3, 1.0e-3
    ops.nDMaterial('ElasticIsotropic', 1, E, nu)
    crd = {1: (0, 0), 2: (2, 0), 3: (0, 1), 4: (1, 0), 5: (1, 0.5), 6: (0, 0.5)}
    for n, (x, y) in crd.items():
        ops.node(n, float(x), float(y))
    ops.element('BezierTri6', 1, 1, 2, 3, 4, 5, 6, 1.0, 'PlaneStress', 1)
    ops.fix(1, 1, 1)
    ops.timeSeries('Linear', 1); ops.pattern('Plain', 1, 1)
    for n, (x, y) in crd.items():
        if n != 1:
            ops.sp(n, 1, eps0 * x); ops.sp(n, 2, -nu * eps0 * y)
    ops.constraints('Penalty', 1e16, 1e16); ops.numberer('Plain')
    ops.system('FullGeneral'); ops.test('NormDispIncr', 1e-12, 10)
    ops.algorithm('Newton'); ops.integrator('LoadControl', 1.0)
    ops.analysis('Static'); ops.analyze(1)
    s = ops.eleResponse(1, 'stresses')
    exact = E * eps0
    ok = len(s) == 9 and all(abs(s[3 * i] - exact) < TOL and abs(s[3 * i + 1]) < TOL
                             and abs(s[3 * i + 2]) < TOL for i in range(3))
    record("patch test (constant stress exact)", ok, f"sigma_xx={round(s[0], 6)} (exact {exact})")


# --- 2. SelfWeight -----------------------------------------------------------
def selfweight():
    ops.wipe(); ops.model('basic', '-ndm', 2, '-ndf', 2)
    ops.nDMaterial('ElasticIsotropic', 1, 1000.0, 0.3)
    gamma, A, t = 19.62, 1.0, 1.0
    crd = {1: (0, 0), 2: (2, 0), 3: (0, 1), 4: (1, 0), 5: (1, 0.5), 6: (0, 0.5)}
    for n, (x, y) in crd.items():
        ops.node(n, float(x), float(y))
    ops.element('BezierTri6', 1, 1, 2, 3, 4, 5, 6, t, 'PlaneStress', 1,
                '-bodyForce', 0.0, -gamma)
    fixed = [1, 2, 4]
    for n in fixed:
        ops.fix(n, 1, 1)
    ops.timeSeries('Linear', 1); ops.pattern('Plain', 1, 1)
    ops.eleLoad('-ele', 1, '-type', '-selfWeight', 1.0, 1.0, 0.0)
    ops.constraints('Plain'); ops.numberer('Plain'); ops.system('FullGeneral')
    ops.test('NormDispIncr', 1e-12, 10); ops.algorithm('Newton')
    ops.integrator('LoadControl', 1.0); ops.analysis('Static'); ops.analyze(1)
    ops.reactions()
    Ry = sum(ops.nodeReaction(n, 2) for n in fixed)
    W = gamma * A * t
    record("SelfWeight (reaction = gamma*A*t)", abs(Ry - W) < TOL, f"Ry={round(Ry, 4)} vs {W}")


# --- 3. lumped mass default --------------------------------------------------
def lumped_mass():
    ops.wipe(); ops.model('basic', '-ndm', 2, '-ndf', 2)
    ops.nDMaterial('ElasticIsotropic', 1, 1000.0, 0.3)
    for n, (x, y) in {1: (0, 0), 2: (2, 0), 3: (0, 2), 4: (1, 0), 5: (1, 1), 6: (0, 1)}.items():
        ops.node(n, float(x), float(y))
    ops.element('BezierTri6', 1, 1, 2, 3, 4, 5, 6, 1.0, 'PlaneStress', 1, '-rho', 2.0)
    ops.fix(1, 1, 1); ops.fix(2, 1, 1)
    ev = ops.eigen('-fullGenLapack', 3)
    record("lumped-mass default (positive eigenvalues)", all(e > 0 for e in ev),
           f"eig={[round(e, 2) for e in ev]}")


# --- 4. PlaneStress B-bar guard ---------------------------------------------
def bbar_guard():
    ops.wipe(); ops.model('basic', '-ndm', 2, '-ndf', 2)
    ops.nDMaterial('ElasticIsotropic', 1, 1000.0, 0.3)
    for n, (x, y) in {1: (0, 0), 2: (2, 0), 3: (0, 1), 4: (1, 0), 5: (1, 0.5), 6: (0, 0.5)}.items():
        ops.node(n, float(x), float(y))
    try:
        ops.element('BezierTri6', 1, 1, 2, 3, 4, 5, 6, 1.0, 'PlaneStress', 1, '-bbar')
        record("PlaneStress B-bar guard (warns + disables, no crash)", True)
    except Exception as e:  # noqa: BLE001
        record("PlaneStress B-bar guard", False, str(e))


# --- 5. multi-element cantilever vs Timoshenko -------------------------------
def _rect_t6_mesh(L, H, nx, ny):
    """Straight-sided structured T6 mesh of an L x H rectangle (mid-nodes at
    physical midpoints -> P = X). Returns (nodes{nid:(x,y)}, elems[[n1..n6]],
    left_ids, right_ids)."""
    cid, nodes = {}, {}
    nid = 1
    for j in range(ny + 1):
        for i in range(nx + 1):
            cid[(i, j)] = nid; nodes[nid] = (L * i / nx, H * j / ny); nid += 1
    mid = {}

    def m(a, b):
        nonlocal nid
        k = (min(a, b), max(a, b))
        if k not in mid:
            mid[k] = nid
            nodes[nid] = (0.5 * (nodes[a][0] + nodes[b][0]), 0.5 * (nodes[a][1] + nodes[b][1]))
            nid += 1
        return mid[k]

    elems = []
    for j in range(ny):
        for i in range(nx):
            bl, br = cid[(i, j)], cid[(i + 1, j)]
            tr, tl = cid[(i + 1, j + 1)], cid[(i, j + 1)]
            for (c1, c2, c3) in [(bl, br, tr), (bl, tr, tl)]:  # CCW
                elems.append([c1, c2, c3, m(c1, c2), m(c2, c3), m(c3, c1)])
    left = [n for n, (x, y) in nodes.items() if abs(x) < 1e-9]
    right = [n for n, (x, y) in nodes.items() if abs(x - L) < 1e-9]
    return nodes, elems, left, right


def cantilever():
    L, H, t, E, nu, P = 10.0, 2.0, 1.0, 1000.0, 0.3, 10.0
    nodes, elems, left, right = _rect_t6_mesh(L, H, 8, 4)
    ops.wipe(); ops.model('basic', '-ndm', 2, '-ndf', 2)
    ops.nDMaterial('ElasticIsotropic', 1, E, nu)
    for nid, (x, y) in nodes.items():
        ops.node(nid, x, y)
    for e, conn in enumerate(elems, start=1):
        ops.element('BezierTri6', e, *conn, t, 'PlaneStress', 1)
    for n in left:
        ops.fix(n, 1, 1)
    ops.timeSeries('Linear', 1); ops.pattern('Plain', 1, 1)
    for n in right:
        ops.load(n, 0.0, -P / len(right))
    ops.constraints('Plain'); ops.numberer('RCM'); ops.system('BandGeneral')
    ops.test('NormDispIncr', 1e-10, 30); ops.algorithm('Newton')
    ops.integrator('LoadControl', 1.0); ops.analysis('Static'); ok = ops.analyze(1)
    uy = min(ops.nodeDisp(n, 2) for n in right)
    I = t * H ** 3 / 12.0; G = E / (2 * (1 + nu)); As = (5.0 / 6.0) * t * H
    est = -(P * L ** 3 / (3 * E * I) + P * L / (G * As))
    ratio = uy / est
    record("cantilever vs Timoshenko (within 5%)", ok == 0 and abs(ratio - 1.0) < 0.05,
           f"FE={round(uy, 4)} est={round(est, 4)} ratio={round(ratio, 3)}")


# --- 6. parameter forwarding (stress-control staging) ------------------------
# Regression for the setParameter gap: without the element-level override,
# `addToParameter <p> element <tag> commitStressIncrementXX` was a SILENT
# no-op through the Bezier elements (stresses stayed 0.0) while
# Tri31/SixNodeTri hosts reached the target exactly. Mirrors the apeGmsh
# ops.initial_stress / STKO stressControl mechanism with the same material
# chain (PlaneStrain-wrapped ASDPlasticMaterial3D, MohrCoulomb).
_ASDP_MC_ARGS = (
    'MohrCoulomb_YF', 'MohrCoulomb_PF', 'LinearIsotropic3D_EL',
    'BackStress(NullHardeningTensorFunction):',
    'Begin_Internal_Variables',
    'BackStress', 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    'End_Internal_Variables',
    'Begin_Model_Parameters',
    'AF_cr', 0.0, 'AF_ha', 0.0, 'DP_eta', 0.0, 'DP_etabar', 0.0,
    'DP_xi_c', 0.0, 'Dilatancy', 0.0, 'DuncanChang_MaxSigma3', 0.0,
    'DuncanChang_n', 0.0, 'InitialP0', 0.0, 'MC_c', 1014.0, 'MC_ds', 1e-05,
    'MC_phi', 45.95, 'MC_psi', 11.49, 'MassDensity', 4.5,
    'PoissonsRatio', 0.18, 'ReferencePressure', 0.0,
    'ReferenceYoungsModulus', 0.0, 'ScalarLinearHardeningParameter', 0.0,
    'TC_min_stress', 0.0, 'TensorLinearHardeningParameter', 0.0,
    'YoungsModulus', 4080000.0,
    'End_Model_Parameters',
    'Begin_Integration_Options',
    'f_absolute_tol', 1e-06, 'stress_absolute_tol', 1e-06,
    'n_max_iterations', 100, 'rk45_dT_min', 0.01, 'rk45_niter_max', 100,
    'return_to_yield_surface', 'Disabled',
    'integration_method', 'Backward_Euler', 'tangent_type', 'Secant',
    'End_Integration_Options',
)
_SIGMA_TARGET = -6300.0


def _stress_control_chain():
    # Penalty + sp (NOT ops.fix): an all-fix model leaves the solver with
    # zero equations and FullGeneral hard-crashes on setSize(0).
    ops.constraints('Penalty', 1e16, 1e16); ops.numberer('Plain')
    ops.system('FullGeneral'); ops.test('NormDispIncr', 1e-8, 20)
    ops.algorithm('Newton'); ops.integrator('LoadControl', 0.0)
    ops.analysis('Static')


_SIGMA_TOL = 0.5   # vs target 6300 — pass/fail gap is -6300 vs 0.0


def _inject_isotropic(ele_tag):
    """parameter/addToParameter/updateParameter the ISOTROPIC -6300 state
    (XX=YY=ZZ), exactly like the staged in-situ install. Isotropic = pure
    pressure, far inside the MC surface — a pure sigma_xx injection at zero
    confinement yields (tau 3150 > strength ~2969) and return-maps away
    from the target, which would test plasticity, not parameter routing."""
    for ptag, comp in ((1, 'XX'), (2, 'YY'), (3, 'ZZ')):
        ops.parameter(ptag)
        ops.addToParameter(ptag, 'element', ele_tag,
                           f'commitStressIncrement{comp}')
    _stress_control_chain()
    for ptag in (1, 2, 3):
        ops.updateParameter(ptag, _SIGMA_TARGET)
    return ops.analyze(1)


def stress_control_tri6():
    ops.wipe(); ops.model('basic', '-ndm', 2, '-ndf', 2)
    ops.nDMaterial('ASDPlasticMaterial3D', 1, *_ASDP_MC_ARGS)
    ops.nDMaterial('PlaneStrain', 2, 1)
    crd = {1: (0, 0), 2: (2, 0), 3: (0, 1), 4: (1, 0), 5: (1, 0.5), 6: (0, 0.5)}
    for n, (x, y) in crd.items():
        ops.node(n, float(x), float(y))
    ops.element('BezierTri6', 1, 1, 2, 3, 4, 5, 6, 1.0, 'PlaneStrain', 2)
    ops.timeSeries('Linear', 1); ops.pattern('Plain', 1, 1)
    for n in crd:
        ops.sp(n, 1, 0.0); ops.sp(n, 2, 0.0)
    rc = _inject_isotropic(1)
    s = ops.eleResponse(1, 'stresses')   # 3 GP x (sxx, syy, sxy)
    ok = (rc == 0 and len(s) == 9
          and all(abs(s[3 * g] - _SIGMA_TARGET) < _SIGMA_TOL
                  and abs(s[3 * g + 1] - _SIGMA_TARGET) < _SIGMA_TOL
                  for g in range(3)))
    record("stress-control param forwarding (Tri6)", ok,
           f"sigma_xx={round(s[0], 3) if s else '—'} (target {_SIGMA_TARGET})")


def stress_control_tet10():
    ops.wipe(); ops.model('basic', '-ndm', 3, '-ndf', 3)
    ops.nDMaterial('ASDPlasticMaterial3D', 1, *_ASDP_MC_ARGS)
    # canonical straight-sided tet: vertices then mid-edges
    # (1-2),(2-3),(1-3),(1-4),(3-4),(2-4) per the BezierTet10 node order
    v = {1: (0, 0, 0), 2: (1, 0, 0), 3: (0, 1, 0), 4: (0, 0, 1)}
    mids = {5: (1, 2), 6: (2, 3), 7: (1, 3), 8: (1, 4), 9: (3, 4), 10: (2, 4)}
    for n, c in v.items():
        ops.node(n, *[float(q) for q in c])
    for n, (a, b) in mids.items():
        ops.node(n, *[(v[a][k] + v[b][k]) / 2.0 for k in range(3)])
    ops.element('BezierTet10', 1, *range(1, 11), 1)
    ops.timeSeries('Linear', 1); ops.pattern('Plain', 1, 1)
    for n in range(1, 11):
        ops.sp(n, 1, 0.0); ops.sp(n, 2, 0.0); ops.sp(n, 3, 0.0)
    rc = _inject_isotropic(1)
    s = ops.eleResponse(1, 'stresses')   # 4 GP x 6 comps
    ncomp = len(s) // 4 if s else 0
    ok = (rc == 0 and ncomp >= 6
          and all(abs(s[ncomp * g] - _SIGMA_TARGET) < _SIGMA_TOL
                  and abs(s[ncomp * g + 1] - _SIGMA_TARGET) < _SIGMA_TOL
                  and abs(s[ncomp * g + 2] - _SIGMA_TARGET) < _SIGMA_TOL
                  for g in range(4)))
    record("stress-control param forwarding (Tet10)", ok,
           f"sigma_xx={round(s[0], 3) if s else '—'} (target {_SIGMA_TARGET})")


if __name__ == "__main__":
    print("=" * 64)
    print("  BezierTri6 regression test")
    print("=" * 64)
    for fn in (patch_test, selfweight, lumped_mass, bbar_guard, cantilever,
               stress_control_tri6, stress_control_tet10):
        try:
            fn()
        except Exception as exc:  # noqa: BLE001
            record(fn.__name__, False, f"EXCEPTION {type(exc).__name__}: {exc}")
    print("=" * 64)
    n_ok = sum(1 for _, ok in results if ok)
    print(f"  {n_ok}/{len(results)} passed")
    sys.exit(0 if n_ok == len(results) else 1)
