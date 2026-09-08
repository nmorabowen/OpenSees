"""ADR-97 P0 oracle 6/6 -- ELEMENT-LEVEL free-DOF finite-difference tangent check
for the fork itself.  Model-independent: it never consults a numpy oracle, only
the binary's own assembled residual.

WHAT IT MEASURES
----------------
Two runs of the SAME one-element, one-step path, from the same virgin state:

* **load-driven rig** -- top-face DOFs FREE under nodal loads.  After the step,
  ``printA('-sparse','-ret')`` is the tangent OpenSees assembled from
  ``NDMaterial::getTangent()``: ``K_asm``.
* **displacement-driven rig** -- every one of those same DOFs prescribed with
  ``sp`` to the converged value ``u*`` plus/minus ``h`` on one DOF at a time.
  With all of them prescribed the free-DOF count is zero, so the analysis is a
  pure state determination and ``nodeReaction`` returns the assembled internal
  force ``R(u)``.  Central differences give ``K_fd = dR/du``.

``K_fd`` is the true consistent tangent of the map the binary COMMITS.  A
material whose ``tangent_type`` is the algorithmic tangent of its own return map
lands at ``rel_err ~ 1e-9``; everything else exposes itself.

TRAPS OBEYED
------------
* ``printA('-sparse','-ret')`` -- the dense ``-ret`` path is empty except under
  ``FullGeneral``, and ``FullGeneral`` crashes the fully-prescribed rig.  Both
  rigs use ``UmfPack``.
* Equation numbers come from ``ops.nodeDOFs`` -- never assumed from the
  numberer.
* One step only.  The load-driven and displacement-driven rigs traverse the
  IDENTICAL material path only when the step is taken in one increment (the
  load path is not a straight line in u-space once plastic).
* ``ops.ladrunoBuild()`` is printed first -- a wrong or stale binary is the
  first thing to rule out.

NEGATIVE CONTROL (measured on 3622d6214, ``stdBrick`` + VonMises, dEps_zz =
-3.6e-3): ``Backward_Euler`` with ``tangent_type Continuum`` gives
``rel_err = 0.5734``.  ADR-94 M3 quotes 57.3% for the same state against a
numpy reference tangent; this driver reproduces it WITHOUT the oracle, which
independently confirms that the ADR-94 reference was the consistent tangent.

Use from pytest::

    from fd_tangent_driver import fd_check, mat_vm
    r = fd_check(lambda t: mat_vm(t, tangent="Algorithmic",
                                  method="Closest_Point"))
    assert r["rel_err"] < 1e-6

Run standalone::

    LADRUNO_OPENSEES_QUIET=1 PYTHONPATH=<...>/dist/bin \\
        python3.12 Ladruno_implementation/adr97_oracle/fd_tangent_driver.py
"""
import contextlib
import os
import sys

import numpy as np

try:                                   # inside tests/ (pytest puts it on path)
    from _testbed import ops
except ImportError:                    # standalone
    try:
        import opensees as ops         # local source build
    except ModuleNotFoundError:        # installed wheel
        import openseespy.opensees as ops

@contextlib.contextmanager
def quiet_ops():
    """Silence the C++ side's own stdout (ASDPlasticMaterial3D prints a page per
    construction).  Set ADR97_VERBOSE_OPS=1 to keep it."""
    if os.environ.get("ADR97_VERBOSE_OPS"):
        yield
        return
    sys.stdout.flush()
    saved = os.dup(1)
    null = os.open(os.devnull, os.O_WRONLY)
    os.dup2(null, 1)
    os.close(null)
    try:
        yield
    finally:
        sys.stdout.flush()
        os.dup2(saved, 1)
        os.close(saved)


NODES = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
         (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
TOP = (5, 6, 7, 8)

# ---------------------------------------------------------------------------
# the ADR-94 von Mises control model (same numbers as
# tests/test_adr94_hlist_numerics.py::mat_vm)
# ---------------------------------------------------------------------------
E_VM, NU_VM, SY_VM, H_VM = 70000.0, 0.3, 30.0, 7000.0
IV_VM = ("BackStress(TensorLinearHardeningFunction):"
         "YieldStress(ScalarLinearHardeningFunction):")


def mat_vm(tag, tangent="Continuum", method="Backward_Euler", hiso=H_VM,
           niter=100):
    ops.nDMaterial(
        "ASDPlasticMaterial3D", tag,
        "VonMises_YF", "VonMises_PF", "LinearIsotropic3D_EL", IV_VM,
        "Begin_Model_Parameters",
        "YoungsModulus", E_VM, "PoissonsRatio", NU_VM,
        "ScalarLinearHardeningParameter", hiso,
        "TensorLinearHardeningParameter", 0.0, "MassDensity", 0.0,
        "End_Model_Parameters",
        "Begin_Internal_Variables", "YieldStress", SY_VM,
        "BackStress", 0., 0., 0., 0., 0., 0., "End_Internal_Variables",
        "Begin_Integration_Options",
        "integration_method", method, "tangent_type", tangent,
        "n_max_iterations", int(niter),
        "End_Integration_Options",
    )


def mat_elastic(tag, E=E_VM, nu=NU_VM):
    ops.nDMaterial("ElasticIsotropic", tag, E, nu)


# ---------------------------------------------------------------------------
def _dof_list(rig):
    """(node, dof) pairs that are FREE in the load-driven rig."""
    if rig == "uniaxial":
        return [(k, 3) for k in TOP]
    return [(k, d) for k in TOP for d in (1, 2, 3)]


def _build(mat_fn, ele, rig):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for k, (x, y, z) in enumerate(NODES):
        ops.node(k + 1, float(x), float(y), float(z))
    for k in range(1, 5):
        ops.fix(k, 1, 1, 1)
    if rig == "uniaxial":
        for k in TOP:
            ops.fix(k, 1, 1, 0)
    mat_fn(1)
    ops.element(ele, 1, *range(1, 9), 1)


def _analysis(tol=1e-13, maxit=100):
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("UmfPack")            # NOT FullGeneral: it crashes the sp rig
    ops.test("NormDispIncr", tol, maxit, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")


def _run_load(mat_fn, ele, rig, load):
  with quiet_ops():
      _build(mat_fn, ele, rig)
      ops.timeSeries("Linear", 1)
      ops.pattern("Plain", 1, 1)
      for k in TOP:
          ops.load(k, *load)
      _analysis()
      rc = ops.analyze(1)
      dofs = _dof_list(rig)
      eqn = {}
      for k in TOP:
          for d, e in enumerate(ops.nodeDOFs(k)):
              eqn[(k, d + 1)] = e
      n = len(dofs)
      d = ops.printA("-sparse", "-ret")
      ndof = max(eqn.values()) + 1
      Kfull = np.zeros((ndof, ndof))
      for i, j, v in zip(d["rowIndices"], d["colIndices"], d["values"]):
          Kfull[i, j] += v
      idx = [eqn[p] for p in dofs]
      K = Kfull[np.ix_(idx, idx)]
      u = np.array([ops.nodeDisp(k, dd) for (k, dd) in dofs])
      eps = np.array(list(ops.eleResponse(1, "strains"))[0:6])
      sig = np.array(list(ops.eleResponse(1, "stresses"))[0:6])
      return rc, K, u, eps, sig


def _run_prescribed(mat_fn, ele, rig, u):
  with quiet_ops():
      _build(mat_fn, ele, rig)
      ops.timeSeries("Linear", 1)
      ops.pattern("Plain", 1, 1)
      for (k, dd), val in zip(_dof_list(rig), u):
          ops.sp(k, dd, float(val))
      _analysis()
      rc = ops.analyze(1)
      ops.reactions()
      R = np.array([ops.nodeReaction(k, dd) for (k, dd) in _dof_list(rig)])
      return rc, R


def fd_check(mat_fn, ele="stdBrick", rig="uniaxial", load=(0.0, 0.0, -60.0),
             h=1e-7, verbose=True, label=""):
    """Assembled free-DOF tangent vs a central difference of the binary's own
    internal force.  Returns ``dict(rel_err, rel_fro, K_fd, K_asm, u, eps, sig,
    rc)``."""
    rc, K_asm, u, eps, sig = _run_load(mat_fn, ele, rig, load)
    if rc != 0:
        if verbose:
            print(f"  {label:52s} LOAD-DRIVEN RIG DID NOT CONVERGE "
                  f"(analyze -> {rc}); no measurement")
        return dict(rel_err=float("nan"), rel_fro=float("nan"), K_fd=None,
                    K_asm=K_asm, u=u, eps=eps, sig=sig, rc=rc)
    n = len(u)
    K_fd = np.zeros((n, n))
    for j in range(n):
        s = h * max(1.0, abs(u[j]))
        up = u.copy(); up[j] += s
        um = u.copy(); um[j] -= s
        rcp, Rp = _run_prescribed(mat_fn, ele, rig, up)
        rcm, Rm = _run_prescribed(mat_fn, ele, rig, um)
        assert rcp == 0 and rcm == 0, "prescribed rig failed to converge"
        K_fd[:, j] = (Rp - Rm) / (2.0 * s)
    scale = max(float(np.max(np.abs(K_fd))), 1e-30)
    rel = float(np.max(np.abs(K_fd - K_asm))) / scale
    fro = float(np.linalg.norm(K_fd - K_asm) / max(np.linalg.norm(K_fd), 1e-30))
    if verbose:
        print(f"  {label:52s} rel_err(max) = {rel:.6e}   rel_fro = {fro:.6e}"
              f"   |K_fd|_F = {np.linalg.norm(K_fd):.3f}")
    return dict(rel_err=rel, rel_fro=fro, K_fd=K_fd, K_asm=K_asm, u=u,
                eps=eps, sig=sig, rc=0)


def try_material(**kw):
    """Return a mat_fn, or None if this build's parser refuses the options."""
    def fn(tag):
        mat_vm(tag, **kw)
    try:
        with quiet_ops():
            ops.wipe()
            ops.model("basic", "-ndm", 3, "-ndf", 3)
            fn(1)
            ops.wipe()
        return fn
    except Exception:
        with quiet_ops():
            ops.wipe()
        return None


# ---------------------------------------------------------------------------
def main():
    print("=" * 78)
    print("ADR-97 oracle 6/6 -- element-level free-DOF FD tangent check")
    print("=" * 78)
    print("build:", ops.ladrunoBuild())
    print()

    print("0. SELF-TEST -- ElasticIsotropic: K_asm IS the exact tangent, so the")
    print("   finite difference must reproduce it to round-off.")
    r0 = fd_check(mat_elastic, label="ElasticIsotropic (self-test)")
    assert r0["rel_err"] < 1e-8, r0["rel_err"]
    print("   -> the rig and the sign convention are verified.\n")

    print("1. NEGATIVE CONTROL -- ASDPlasticMaterial3D VonMises, "
          "Backward_Euler,\n   one step to sigma_zz = -240 (deep plastic), "
          "the five tangent_type options:")
    res = {}
    for tg in ("Continuum", "Secant", "Elastic",
               "Numerical_Algorithmic_FirstOrder",
               "Numerical_Algorithmic_SecondOrder"):
        fn = try_material(tangent=tg)
        if fn is None:
            print(f"  {tg:52s} NOT ACCEPTED by this build's parser")
            continue
        res[tg] = fd_check(fn, label=f"Backward_Euler / {tg}")
    print()
    print("   eps  =", np.array2string(res["Continuum"]["eps"], precision=10))
    print("   sig  =", np.array2string(res["Continuum"]["sig"], precision=8))
    print("   u*   =", np.array2string(res["Continuum"]["u"], precision=12))
    print("   K_asm (Continuum) =")
    for row in res["Continuum"]["K_asm"]:
        print("     " + " ".join("%14.6f" % v for v in row))
    print("   K_fd  (the consistent tangent of the committed map) =")
    for row in res["Continuum"]["K_fd"]:
        print("     " + " ".join("%14.6f" % v for v in row))

    print("\n2. THE ADR-97 TARGET -- integration_method Closest_Point with")
    print("   tangent_type Algorithmic (does not exist yet; P1 adds it):")
    fn = try_material(tangent="Algorithmic", method="Closest_Point")
    if fn is None:
        print("  Closest_Point / Algorithmic                          NOT "
              "AVAILABLE in this build (expected before P1)")
    else:
        r = fd_check(fn, label="Closest_Point / Algorithmic")
        print(f"  gate: rel_err <= 1e-6  ->  "
              f"{'PASS' if r['rel_err'] <= 1e-6 else 'FAIL'}")

    print("\n3. A SHEARED, 12-FREE-DOF RIG (top face fully free) -- the same")
    print("   measurement where the tangent's shear block is exercised:")
    for tg in ("Continuum", "Numerical_Algorithmic_FirstOrder"):
        fn = try_material(tangent=tg)
        if fn is not None:
            fd_check(fn, rig="full", load=(-3.0, 0.0, -13.0),
                     label=f"12-DOF rig, Backward_Euler / {tg}")

    print("\n" + "=" * 78)
    print("ELEMENT-LEVEL FD REFERENCE BLOCK (values C++ tests pin)")
    print("=" * 78)
    print("  elastic self-test rel_err          = %.3e" % r0["rel_err"])
    for tg, r in res.items():
        print(f"  Backward_Euler / {tg:34s} rel_err = {r['rel_err']:.6f}")
    print("  (ADR-94 M3 quoted 57.3 / 79.9 / 102.5 / 4.6 / 4.6 % against a "
          "numpy\n   reference tangent; this driver reproduces them from the "
          "binary alone.)")


if __name__ == "__main__":
    main()
