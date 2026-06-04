"""T0 / T0m element & material verification primitives, with the locked
tolerances from the canonical test bed (Axis 2). See
Ladruno_implementation/testbed/00_canonical_testbed.md §2.

All tolerances are RELATIVE and conditioning-aware — never machine-eps, which
false-fails because zero-energy eigenvalues come back at ~1e-9*lambda_max.
"""
from ._ops import ops

# locked constants
ZERO_ENERGY_REL = 1e-8   # lambda_i < ZERO_ENERGY_REL * lambda_max  => zero mode
PATCH_RTOL = 1e-9        # gauss-point stress vs closed-form constant stress
RIGID_FORCE_REL = 1e-8   # ||F_resisting|| < RIGID_FORCE_REL * E * A
TANGENT_FD_RTOL = 1e-6   # consistent tangent vs finite difference


def n_rigid_body_modes(ndm: int) -> int:
    """Expected zero-energy mode count for a free element: 1 (1D), 3 (2D), 6 (3D)."""
    return {1: 1, 2: 3, 3: 6}[ndm]


def zero_energy_mode_count(node_tags, ndf, tol_rel=ZERO_ENERGY_REL):
    """Count zero-energy modes of a SINGLE FREE element already defined in the
    current domain. Sets unit mass on every DOF so the generalized eigenproblem
    K.phi = lambda.M.phi reduces to the eigenvalues of K, then counts
    lambda < tol_rel * lambda_max.

    Returns (count, eigenvalues).
    """
    node_tags = [int(n) for n in node_tags]
    for n in node_tags:
        ops.mass(n, *([1.0] * ndf))
    n_modes = len(node_tags) * ndf
    eigs = ops.eigen("-fullGenLapack", n_modes)
    lam_max = max(abs(e) for e in eigs)
    cut = tol_rel * lam_max
    count = sum(1 for e in eigs if abs(e) < cut)
    return count, eigs


def assert_zero_energy(node_tags, ndf, ndm, tol_rel=ZERO_ENERGY_REL):
    """Assert the free element has EXACTLY n_rigid_body_modes(ndm) zero modes.
    n_rbm+1 means a spurious hourglass mode / rank deficiency."""
    count, eigs = zero_energy_mode_count(node_tags, ndf, tol_rel)
    expected = n_rigid_body_modes(ndm)
    assert count == expected, (
        f"zero-energy mode count {count} != {expected} rigid-body modes "
        f"(extra zeros = spurious mechanism/hourglass). eigenvalues={eigs}"
    )
    return count, eigs


def check_constant_stress(ele_tag, n_gauss, target_stress, E, strain_ref,
                          rtol=PATCH_RTOL):
    """Patch-test core: every Gauss point of a constant-strain patch must report
    the same closed-form constant stress. atol floor scales with E*strain so the
    check is meaningful for any material stiffness.

    target_stress : list of stress components expected at every Gauss point.
    """
    atol = 1e-9 * E * abs(strain_ref)
    for gp in range(1, n_gauss + 1):
        s = ops.eleResponse(int(ele_tag), "material", gp, "stress")
        assert s, (
            f"element {ele_tag} returned no material stress at gp {gp} — "
            "missing setResponse('material', ...) branch?"
        )
        for comp, (si, ti) in enumerate(zip(s, target_stress)):
            assert abs(si - ti) <= rtol * abs(ti) + atol, (
                f"patch test failed: ele {ele_tag} gp {gp} comp {comp}: "
                f"{si} != {ti} (rtol {rtol}, atol {atol})"
            )


def uniaxial_tangent_fd(mat_tag, strain0, dstrain=1e-7, rtol=TANGENT_FD_RTOL):
    """T0m: central finite-difference the consistent tangent of a uniaxial
    material against getTangent() at strain0. Returns (k_analytic, k_fd).

    WARNING — valid only for PATH-INDEPENDENT (nonlinear-elastic) materials.
    OpenSees ``setStrain`` calls ``setTrialStrain`` + ``commitState`` (see
    OpenSeesCommandsPython.cpp), so for a plasticity material the ``strain0-d``
    probe UNLOADS elastically and the central difference returns ~(E+E_alg)/2,
    not the consistent tangent. For path-dependent materials probe each point as
    an independent one-step return from a fresh material instead (see the V6 test
    in tests/test_ladrunoUniaxialJ2_material.py)."""
    ops.testUniaxialMaterial(int(mat_tag))
    ops.setStrain(strain0)
    k = ops.getTangent()
    ops.setStrain(strain0 - dstrain)
    s_m = ops.getStress()
    ops.setStrain(strain0 + dstrain)
    s_p = ops.getStress()
    k_fd = (s_p - s_m) / (2.0 * dstrain)
    assert abs(k - k_fd) <= rtol * abs(k) + 1e-9, (
        f"tangent FD mismatch: analytic {k} vs FD {k_fd} (rtol {rtol})"
    )
    return k, k_fd
