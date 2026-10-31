"""Wave-propagation verification for LadrunoCPML (Anchors 3a, 3b, and 4).

Covers:
  1. Anchor 3a — 2D Full-Space Cylindrical Wave Radiation (Green's Function)
     across broad incidence angles (0 deg, 45 deg corner, ~68 deg high-oblique)
     on both structured Q4 (quad + LadrunoCPML) and unstructured/skewed
     T3 (LadrunoCST + LadrunoCPML) meshes against an extended-mesh reference.
  2. Anchor 3b — 3D Point-Source Spherical Wave Radiation across 3D coordinate
     faces, edges, and the 3D (x, y, z) corner (where the 3-stage
     phi_2 o g_1 o g_0 cascade is active) on both H8 (stdBrick + LadrunoCPML)
     and TET4 (FourNodeTetrahedron + LadrunoCPML) meshes.
  3. Anchor 4 — Lamb's Problem on a 2D Elastic Half-Space (surface point load
     exciting Rayleigh surface waves and grazing P/SV waves): head-to-head
     comparison of LadrunoCPML vs. ASDAbsorbingBoundary2D against an
     extended-mesh reference.
"""

import math
import numpy as np
import pytest
import opensees as ops

pytestmark = [pytest.mark.zone_a]


def _ricker_series(dt: float, n_steps: int, f0: float, t0: float, amp: float = 1.0) -> list[float]:
    """Generate discrete Ricker wavelet time series."""
    vals = []
    for k in range(n_steps + 2):
        t = k * dt
        u = math.pi * f0 * (t - t0)
        u2 = u * u
        vals.append(amp * (1.0 - 2.0 * u2) * math.exp(-u2))
    return vals


def _rel_l2_and_peak_error(u_test: np.ndarray, u_ref: np.ndarray) -> tuple[float, float]:
    """Compute relative L2 and relative peak reflection errors over vector waveform."""
    diff = u_test - u_ref
    l2_ref = float(np.linalg.norm(u_ref))
    l2_err = float(np.linalg.norm(diff)) / max(l2_ref, 1.0e-30)

    mag_ref = np.linalg.norm(u_ref, axis=1) if u_ref.ndim > 1 else np.abs(u_ref)
    mag_diff = np.linalg.norm(diff, axis=1) if diff.ndim > 1 else np.abs(diff)
    peak_ref = float(np.max(mag_ref))
    peak_err = float(np.max(mag_diff)) / max(peak_ref, 1.0e-30)
    return l2_err, peak_err


# ==============================================================================
# TEST GROUP 1 — ANCHOR 3a: 2D CYLINDRICAL WAVE RADIATION (Q4 & UNSTRUCTURED T3)
# ==============================================================================

def _run_2d_cylindrical_q4(mode: str) -> dict[str, np.ndarray]:
    """Run 2D quadrant cylindrical radiation model on structured Q4 mesh.

    mode:
      - 'reference': extended domain [0, 22m]^2 of FourNodeQuad ('quad') elements
      - 'cpml': [0, 6m]^2 'quad' interior + [6, 9m] LadrunoCPML layer
      - 'truncated_fixed': [0, 6m]^2 'quad' interior with fixed outer wall
    """
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)

    dx = 0.5
    L_int = 6.0
    L_pml = 3.0
    if mode == "reference":
        L_tot = 22.0
    elif mode == "cpml":
        L_tot = L_int + L_pml
    elif mode == "truncated_fixed":
        L_tot = L_int
    else:
        raise ValueError(f"Unknown mode: {mode}")

    nx = int(round(L_tot / dx))
    ny = nx

    # Isotropic elastic properties: Vp = 1000 m/s, Vs = 577.35 m/s
    rho = 1000.0
    nu = 0.25
    E = 8.333333333333334e8
    thick = 1.0

    ops.nDMaterial("ElasticIsotropic", 1, E, nu, rho)

    def nid(i: int, j: int) -> int:
        return 1 + j * (nx + 1) + i

    for j in range(ny + 1):
        y = j * dx
        for i in range(nx + 1):
            x = i * dx
            n = nid(i, j)
            ops.node(n, x, y)
            fix_x = 1 if (i == 0 or (mode in ("cpml", "truncated_fixed") and i == nx)) else 0
            fix_y = 1 if (j == 0 or (mode in ("cpml", "truncated_fixed") and j == ny)) else 0
            if fix_x or fix_y:
                ops.fix(n, fix_x, fix_y)

    etag = 1
    for j in range(ny):
        yc = (j + 0.5) * dx
        for i in range(nx):
            xc = (i + 0.5) * dx
            n1 = nid(i, j)
            n2 = nid(i + 1, j)
            n3 = nid(i + 1, j + 1)
            n4 = nid(i, j + 1)
            in_pml = (mode == "cpml") and (xc > L_int or yc > L_int)
            if in_pml:
                ops.element(
                    "LadrunoCPML", etag, n1, n2, n3, n4,
                    "-E", E, "-nu", nu, "-rho", rho,
                    "-thick", thick,
                    "-box", 0.0, L_int, 0.0, L_int,
                    "-L", L_pml,
                    "-m", 2.0, "-R0", 1.0e-4, "-kappaMax", 3.0, "-f0", 55.0,
                    "-lumped",
                )
            else:
                ops.element(
                    "quad", etag, n1, n2, n3, n4,
                    thick, "PlaneStrain", 1, 0.0, rho, 0.0, 0.0,
                )
            etag += 1

    # Transient radial Ricker source distributed over origin cluster
    dt = 0.0002
    n_steps = 165
    f0 = 55.0
    t0 = 0.014
    ricker = _ricker_series(dt, n_steps, f0, t0, amp=1.0e5)
    ops.timeSeries("Path", 1, "-dt", dt, "-values", *ricker)
    ops.pattern("Plain", 1, 1)

    r_src = 1.2
    for j in range(4):
        for i in range(4):
            x = i * dx
            y = j * dx
            r = math.hypot(x, y)
            if 1.0e-9 < r <= r_src:
                w = math.exp(-2.5 * (r / r_src) ** 2)
                ops.load(nid(i, j), w * (x / r), w * (y / r))

    rec_nodes = {
        "deg0": nid(10, 1),    # (5.0, 0.5)
        "deg45": nid(10, 10),  # (5.0, 5.0)
        "deg68": nid(11, 4),   # (5.5, 2.0) -> ~70 deg / high-oblique
    }

    ops.system("UmfPack")
    ops.numberer("RCM")
    ops.constraints("Plain")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.algorithm("Linear")
    ops.analysis("Transient")

    traces = {k: np.zeros((n_steps, 2), dtype=float) for k in rec_nodes}
    for step in range(n_steps):
        ok = ops.analyze(1, dt)
        assert ok == 0, f"2D Q4 analysis failed at step {step} in mode={mode}"
        for k, node_tag in rec_nodes.items():
            traces[k][step, 0] = ops.nodeDisp(node_tag, 1)
            traces[k][step, 1] = ops.nodeDisp(node_tag, 2)

    return traces


def _run_2d_cylindrical_t3(mode: str) -> dict[str, np.ndarray]:
    """Run 2D quadrant cylindrical radiation model on unstructured/skewed T3 mesh.

    Uses LadrunoCST in the interior and 3-node LadrunoCPML in the exterior, with
    skewed interior/exterior node perturbations and alternating diagonal splits.
    """
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)

    dx = 0.5
    L_int = 6.0
    L_pml = 3.0
    if mode == "reference":
        L_tot = 22.0
    elif mode == "cpml":
        L_tot = L_int + L_pml
    elif mode == "truncated_fixed":
        L_tot = L_int
    else:
        raise ValueError(f"Unknown mode: {mode}")

    nx = int(round(L_tot / dx))
    ny = nx
    i_int = int(round(L_int / dx))

    rho = 1000.0
    nu = 0.25
    E = 8.333333333333334e8
    thick = 1.0

    ops.nDMaterial("ElasticIsotropic", 1, E, nu, rho)

    def nid(i: int, j: int) -> int:
        return 1 + j * (nx + 1) + i

    protected = {(10, 1), (10, 10), (11, 4)}

    for j in range(ny + 1):
        for i in range(nx + 1):
            x = i * dx
            y = j * dx
            if (i, j) == (11, 4):
                x, y = 5.5, 2.2
            elif (
                (i, j) not in protected
                and i not in (0, i_int, nx)
                and j not in (0, i_int, ny)
                and (i > 3 or j > 3)
            ):
                x += 0.06 * dx * math.sin(1.7 * i + 2.3 * j)
                y += 0.06 * dx * math.cos(2.1 * i - 1.3 * j)

            n = nid(i, j)
            ops.node(n, x, y)
            fix_x = 1 if (i == 0 or (mode in ("cpml", "truncated_fixed") and i == nx)) else 0
            fix_y = 1 if (j == 0 or (mode in ("cpml", "truncated_fixed") and j == ny)) else 0
            if fix_x or fix_y:
                ops.fix(n, fix_x, fix_y)

    etag = 1
    for j in range(ny):
        for i in range(nx):
            n1 = nid(i, j)
            n2 = nid(i + 1, j)
            n3 = nid(i + 1, j + 1)
            n4 = nid(i, j + 1)
            in_pml = (mode == "cpml") and (i >= i_int or j >= i_int)

            if (i + j) % 2 == 0:
                tris = [(n1, n2, n3), (n1, n3, n4)]
            else:
                tris = [(n1, n2, n4), (n2, n3, n4)]

            for ta, tb, tc in tris:
                if in_pml:
                    ops.element(
                        "LadrunoCPML", etag, ta, tb, tc,
                        "-E", E, "-nu", nu, "-rho", rho,
                        "-thick", thick,
                        "-box", 0.0, L_int, 0.0, L_int,
                        "-L", L_pml,
                        "-m", 2.0, "-R0", 1.0e-4, "-kappaMax", 3.0, "-f0", 55.0,
                        "-lumped",
                    )
                else:
                    ops.element(
                        "LadrunoCST", etag, ta, tb, tc, 1,
                        "-type", "PlaneStrain", "-thick", thick, "-rho", rho,
                    )
                etag += 1

    dt = 0.0002
    n_steps = 165
    f0 = 55.0
    t0 = 0.014
    ricker = _ricker_series(dt, n_steps, f0, t0, amp=1.0e5)
    ops.timeSeries("Path", 1, "-dt", dt, "-values", *ricker)
    ops.pattern("Plain", 1, 1)

    r_src = 1.2
    for j in range(4):
        for i in range(4):
            x = i * dx
            y = j * dx
            r = math.hypot(x, y)
            if 1.0e-9 < r <= r_src:
                w = math.exp(-2.5 * (r / r_src) ** 2)
                ops.load(nid(i, j), w * (x / r), w * (y / r))

    rec_nodes = {
        "deg0": nid(10, 1),    # (5.0, 0.5)
        "deg45": nid(10, 10),  # (5.0, 5.0)
        "deg68": nid(11, 4),   # (5.5, 2.2)
    }

    ops.system("UmfPack")
    ops.numberer("RCM")
    ops.constraints("Plain")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.algorithm("Linear")
    ops.analysis("Transient")

    traces = {k: np.zeros((n_steps, 2), dtype=float) for k in rec_nodes}
    for step in range(n_steps):
        ok = ops.analyze(1, dt)
        assert ok == 0, f"2D T3 analysis failed at step {step} in mode={mode}"
        for k, node_tag in rec_nodes.items():
            traces[k][step, 0] = ops.nodeDisp(node_tag, 1)
            traces[k][step, 1] = ops.nodeDisp(node_tag, 2)

    return traces


def test_anchor_3a_2d_cylindrical_radiation_q4():
    """Anchor 3a (Q4): 2D cylindrical radiation across 0 deg, 45 deg corner, 68 deg."""
    ref = _run_2d_cylindrical_q4("reference")
    cpml = _run_2d_cylindrical_q4("cpml")
    trunc = _run_2d_cylindrical_q4("truncated_fixed")

    for angle_key in ("deg0", "deg45", "deg68"):
        l2_cpml, peak_cpml = _rel_l2_and_peak_error(cpml[angle_key], ref[angle_key])
        l2_trunc, peak_trunc = _rel_l2_and_peak_error(trunc[angle_key], ref[angle_key])
        print(
            f"[Q4 {angle_key}] CPML L2={l2_cpml:.4f}, peak={peak_cpml:.4f} | "
            f"Trunc L2={l2_trunc:.4f}, peak={peak_trunc:.4f}"
        )
        assert l2_cpml <= 0.03, f"Q4 {angle_key} L2 error {l2_cpml:.4f} > 3.0%"
        assert peak_cpml <= 0.03, f"Q4 {angle_key} peak error {peak_cpml:.4f} > 3.0%"
        assert max(l2_trunc, peak_trunc) > 0.80, (
            f"Q4 {angle_key} truncated boundary error {max(l2_trunc, peak_trunc):.4f} <= 80%"
        )


def test_anchor_3a_2d_cylindrical_radiation_unstructured_t3():
    """Anchor 3a (T3): 2D cylindrical radiation on unstructured/skewed T3 mesh."""
    ref = _run_2d_cylindrical_t3("reference")
    cpml = _run_2d_cylindrical_t3("cpml")
    trunc = _run_2d_cylindrical_t3("truncated_fixed")

    for angle_key in ("deg0", "deg45", "deg68"):
        l2_cpml, peak_cpml = _rel_l2_and_peak_error(cpml[angle_key], ref[angle_key])
        l2_trunc, peak_trunc = _rel_l2_and_peak_error(trunc[angle_key], ref[angle_key])
        print(
            f"[T3 {angle_key}] CPML L2={l2_cpml:.4f}, peak={peak_cpml:.4f} | "
            f"Trunc L2={l2_trunc:.4f}, peak={peak_trunc:.4f}"
        )
        assert l2_cpml <= 0.03, f"T3 {angle_key} L2 error {l2_cpml:.4f} > 3.0%"
        assert peak_cpml <= 0.03, f"T3 {angle_key} peak error {peak_cpml:.4f} > 3.0%"
        assert max(l2_trunc, peak_trunc) > 0.80, (
            f"T3 {angle_key} truncated boundary error {max(l2_trunc, peak_trunc):.4f} <= 80%"
        )


# ==============================================================================
# TEST GROUP 2 — ANCHOR 3b: 3D SPHERICAL WAVE RADIATION (H8 & TET4)
# ==============================================================================

def _run_3d_spherical(topo: str, mode: str) -> dict[str, np.ndarray]:
    """Run 3D octant spherical wave radiation model for H8 or TET4 meshes.

    topo: 'H8' (stdBrick + LadrunoCPML) or 'TET4' (FourNodeTetrahedron + LadrunoCPML)
    mode: 'reference' ([0, 4.5m]^3), 'cpml' ([0, 1.5m]^3 + 2.5m PML shell),
          'truncated_fixed' ([0, 1.5m]^3 with fixed outer faces)
    """
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)

    dx = 0.5
    L_int = 1.5
    L_pml = 2.5
    if mode == "reference":
        L_tot = 4.5
    elif mode == "cpml":
        L_tot = L_int + L_pml
    elif mode == "truncated_fixed":
        L_tot = L_int
    else:
        raise ValueError(f"Unknown mode: {mode}")

    n_side = int(round(L_tot / dx))
    i_int = int(round(L_int / dx))

    rho = 1000.0
    nu = 0.25
    E = 8.333333333333334e8

    ops.nDMaterial("ElasticIsotropic", 1, E, nu, rho)

    def nid(i: int, j: int, k: int) -> int:
        return 1 + (k * (n_side + 1) + j) * (n_side + 1) + i

    protected_3d = {(2, 1, 1), (2, 2, 2)}
    for k in range(n_side + 1):
        for j in range(n_side + 1):
            for i in range(n_side + 1):
                x = i * dx
                y = j * dx
                z = k * dx
                if (
                    topo == "TET4"
                    and (i, j, k) not in protected_3d
                    and i not in (0, i_int, n_side)
                    and j not in (0, i_int, n_side)
                    and k not in (0, i_int, n_side)
                    and (i > 1 or j > 1 or k > 1)
                ):
                    # Deterministic unstructured 3D node skew (including inside 3D PML faces, edges & corner)
                    x += 0.05 * dx * math.sin(1.3 * i + 1.9 * j + 2.7 * k)
                    y += 0.05 * dx * math.cos(2.1 * i - 1.5 * j + 1.1 * k)
                    z += 0.05 * dx * math.sin(1.7 * i + 2.3 * j - 1.9 * k)
                n = nid(i, j, k)
                ops.node(n, x, y, z)
                fx = 1 if (i == 0 or (mode in ("cpml", "truncated_fixed") and i == n_side)) else 0
                fy = 1 if (j == 0 or (mode in ("cpml", "truncated_fixed") and j == n_side)) else 0
                fz = 1 if (k == 0 or (mode in ("cpml", "truncated_fixed") and k == n_side)) else 0
                if fx or fy or fz:
                    ops.fix(n, fx, fy, fz)

    etag = 1
    for k in range(n_side):
        for j in range(n_side):
            for i in range(n_side):
                in_pml = (mode == "cpml") and (i >= i_int or j >= i_int or k >= i_int)
                p000 = nid(i, j, k)
                p100 = nid(i + 1, j, k)
                p110 = nid(i + 1, j + 1, k)
                p010 = nid(i, j + 1, k)
                p001 = nid(i, j, k + 1)
                p101 = nid(i + 1, j, k + 1)
                p111 = nid(i + 1, j + 1, k + 1)
                p011 = nid(i, j + 1, k + 1)

                if topo == "H8":
                    if in_pml:
                        ops.element(
                            "LadrunoCPML", etag,
                            p000, p100, p110, p010, p001, p101, p111, p011,
                            "-E", E, "-nu", nu, "-rho", rho,
                            "-box", 0.0, L_int, 0.0, L_int, 0.0, L_int,
                            "-L", L_pml,
                            "-m", 2.5, "-R0", 1.5e-3, "-kappaMax", 1.0, "-alphaMax", 15.0,
                        )
                    else:
                        ops.element(
                            "stdBrick", etag,
                            p000, p100, p110, p010, p001, p101, p111, p011,
                            1, 0.0, 0.0, 0.0, rho,
                        )
                    etag += 1
                elif topo == "TET4":
                    tets = [
                        (p000, p100, p110, p111),
                        (p000, p100, p111, p101),
                        (p000, p010, p111, p110),
                        (p000, p010, p011, p111),
                        (p000, p001, p101, p111),
                        (p000, p001, p111, p011),
                    ]
                    for n1, n2, n3, n4 in tets:
                        if in_pml:
                            ops.element(
                                "LadrunoCPML", etag, n1, n2, n3, n4,
                                "-E", E, "-nu", nu, "-rho", rho,
                                "-box", 0.0, L_int, 0.0, L_int, 0.0, L_int,
                                "-L", L_pml,
                                "-m", 3.0, "-R0", 5.0e-4, "-kappaMax", 1.0, "-alphaMax", 15.0,
                            )
                        else:
                            ops.element("FourNodeTetrahedron", etag, n1, n2, n3, n4, 1)
                        etag += 1
                else:
                    raise ValueError(f"Unknown topo: {topo}")

    dt = 0.00015
    n_steps = 32
    f0 = 115.0
    t0 = 0.0022
    ricker = _ricker_series(dt, n_steps, f0, t0, amp=1.0e5)
    ops.timeSeries("Path", 1, "-dt", dt, "-values", *ricker)
    ops.pattern("Plain", 1, 1)

    for k in range(2):
        for j in range(2):
            for i in range(2):
                x, y, z = i * dx, j * dx, k * dx
                r = math.sqrt(x * x + y * y + z * z)
                if r > 1.0e-9:
                    ops.load(nid(i, j, k), x / r, y / r, z / r)

    rec_nodes = {
        "axis": nid(2, 1, 1),    # (1.0, 0.5, 0.5) near coordinate face
        "corner": nid(2, 2, 2),  # (1.0, 1.0, 1.0) along (1,1,1) space diagonal into 3D corner
    }

    ops.system("UmfPack")
    ops.numberer("RCM")
    ops.constraints("Plain")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.algorithm("Linear")
    ops.analysis("Transient")

    traces = {k: np.zeros((n_steps, 3), dtype=float) for k in rec_nodes}
    for step in range(n_steps):
        ok = ops.analyze(1, dt)
        assert ok == 0, f"3D {topo} analysis failed at step {step} in mode={mode}"
        for k, node_tag in rec_nodes.items():
            traces[k][step, 0] = ops.nodeDisp(node_tag, 1)
            traces[k][step, 1] = ops.nodeDisp(node_tag, 2)
            traces[k][step, 2] = ops.nodeDisp(node_tag, 3)

    return traces


@pytest.mark.parametrize("topo", ["H8", "TET4"])
def test_anchor_3b_3d_spherical_radiation_faces_edges_corners(topo: str):
    """Anchor 3b: 3D spherical wave radiation across faces, edges, and 3D corner."""
    ref = _run_3d_spherical(topo, "reference")
    cpml = _run_3d_spherical(topo, "cpml")
    trunc = _run_3d_spherical(topo, "truncated_fixed")

    for rec_key in ("axis", "corner"):
        l2_cpml, peak_cpml = _rel_l2_and_peak_error(cpml[rec_key], ref[rec_key])
        l2_trunc, peak_trunc = _rel_l2_and_peak_error(trunc[rec_key], ref[rec_key])
        reduction_l2 = l2_trunc / max(l2_cpml, 1.0e-12)
        print(
            f"[3D {topo} {rec_key}] CPML L2={l2_cpml:.4f}, peak={peak_cpml:.4f} | "
            f"Trunc L2={l2_trunc:.4f}, peak={peak_trunc:.4f} (reduction={reduction_l2:.1f}x)"
        )
        assert l2_cpml <= 0.04, f"3D {topo} {rec_key} L2 error {l2_cpml:.4f} > 4.0%"
        assert peak_cpml <= 0.04, f"3D {topo} {rec_key} peak error {peak_cpml:.4f} > 4.0%"
        assert reduction_l2 >= 10.0, (
            f"3D {topo} {rec_key} error reduction {reduction_l2:.2f}x < 10x vs unabsorbed boundary"
        )


# ==============================================================================
# TEST GROUP 3 — ANCHOR 4: LAMB'S PROBLEM (RAYLEIGH WAVE VS. ASD ABSORBING 2D)
# ==============================================================================

def _run_lamb_halfspace(mode: str) -> np.ndarray:
    """Run Lamb's problem on a 2D half-space with half-symmetry at x = 0 and free surface y = 0.

    Interior domain: [0, 8.0m] x [-6.0m, 0.0m] with dx = dy = 0.5m ('quad' elements).
    mode:
      - 'reference': extended domain [0, 28.0m] x [-20.0m, 0.0m]
      - 'asd': truncated at x = 8.0m, y = -6.0m with ASDAbsorbingBoundary2D (R, B, BR)
      - 'cpml': truncated at x = 8.0m, y = -6.0m with 3.0m LadrunoCPML layer
                (-box -100.0 8.0 -6.0 100.0 -L 3.0)
    Returns surface displacement history (n_steps, 2) at receiver (6.5m, 0.0m).
    """
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)

    dx = 0.5
    Lx = 8.0
    Ly = 6.0
    L_pml = 3.0

    if mode == "reference":
        Lx_tot = 28.0
        Ly_tot = 20.0
    elif mode == "cpml":
        Lx_tot = Lx + L_pml
        Ly_tot = Ly + L_pml
    elif mode == "asd":
        Lx_tot = Lx
        Ly_tot = Ly
    else:
        raise ValueError(f"Unknown mode: {mode}")

    nx = int(round(Lx_tot / dx))
    ny = int(round(Ly_tot / dx))
    nx_int = int(round(Lx / dx))
    ny_int = int(round(Ly / dx))

    rho = 1000.0
    nu = 0.25
    E = 8.333333333333334e8
    G = E / (2.0 * (1.0 + nu))
    thick = 1.0

    ops.nDMaterial("ElasticIsotropic", 1, E, nu, rho)

    def nid(i: int, j: int) -> int:
        return 1 + j * (nx + 2) + i

    for j in range(ny + 1):
        y = -j * dx
        for i in range(nx + 1):
            x = i * dx
            n = nid(i, j)
            ops.node(n, x, y)
            fix_x = 1 if (i == 0 or (mode in ("reference", "cpml") and (i == nx or j == ny))) else 0
            fix_y = 1 if (mode in ("reference", "cpml") and (i == nx or j == ny)) else 0
            if fix_x or fix_y:
                ops.fix(n, fix_x, fix_y)

    etag = 1
    for j in range(ny):
        for i in range(nx):
            n1 = nid(i, j + 1)
            n2 = nid(i + 1, j + 1)
            n3 = nid(i + 1, j)
            n4 = nid(i, j)
            in_pml = (mode == "cpml") and (i >= nx_int or j >= ny_int)
            if in_pml:
                ops.element(
                    "LadrunoCPML", etag, n1, n2, n3, n4,
                    "-E", E, "-nu", nu, "-rho", rho,
                    "-thick", thick,
                    "-box", -100.0, Lx, -Ly, 100.0,
                    "-L", L_pml,
                    "-m", 2.0, "-R0", 1.0e-3, "-kappaMax", 2.0, "-f0", 45.0,
                    "-lumped",
                )
            else:
                ops.element(
                    "quad", etag, n1, n2, n3, n4,
                    thick, "PlaneStrain", 1, 0.0, rho, 0.0, 0.0,
                )
            etag += 1

    if mode == "asd":
        for j in range(ny + 1):
            ops.node(nid(nx + 1, j), (nx + 1) * dx, -j * dx)
        for i in range(nx + 1):
            n_bot = nid(i, ny + 1)
            ops.node(n_bot, i * dx, -(ny + 1) * dx)
            if i == 0:
                ops.fix(n_bot, 1, 0)
        ops.node(nid(nx + 1, ny + 1), (nx + 1) * dx, -(ny + 1) * dx)

        asd_tags = []
        for j in range(ny):
            n1 = nid(nx, j + 1)
            n2 = nid(nx + 1, j + 1)
            n3 = nid(nx + 1, j)
            n4 = nid(nx, j)
            ops.element("ASDAbsorbingBoundary2D", etag, n1, n2, n3, n4, G, nu, rho, thick, "R")
            asd_tags.append(etag)
            etag += 1

        for i in range(nx):
            n1 = nid(i, ny + 1)
            n2 = nid(i + 1, ny + 1)
            n3 = nid(i + 1, ny)
            n4 = nid(i, ny)
            ops.element("ASDAbsorbingBoundary2D", etag, n1, n2, n3, n4, G, nu, rho, thick, "B")
            asd_tags.append(etag)
            etag += 1

        n1 = nid(nx, ny + 1)
        n2 = nid(nx + 1, ny + 1)
        n3 = nid(nx + 1, ny)
        n4 = nid(nx, ny)
        ops.element("ASDAbsorbingBoundary2D", etag, n1, n2, n3, n4, G, nu, rho, thick, "BR")
        asd_tags.append(etag)
        etag += 1

        for atag in asd_tags:
            ops.setParameter("-val", 1, "-ele", atag, "stage")

    # Vertical Ricker surface load near (0, 0) on the free surface y = 0
    dt = 0.0002
    n_steps = 190
    f0 = 45.0
    t0 = 0.016
    ricker = _ricker_series(dt, n_steps, f0, t0, amp=-1.0e5)
    ops.timeSeries("Path", 1, "-dt", dt, "-values", *ricker)
    ops.pattern("Plain", 1, 1)
    # Gaussian-weighted surface traction around origin (0, 0) to avoid grid-cutoff dispersion
    ops.load(nid(0, 0), 0.0, 0.5)
    ops.load(nid(1, 0), 0.0, 0.8)
    ops.load(nid(2, 0), 0.0, 0.4)
    ops.load(nid(3, 0), 0.0, 0.1)

    rec_node = nid(13, 0)  # (6.5m, 0.0m)

    ops.system("UmfPack")
    ops.numberer("RCM")
    ops.constraints("Plain")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.algorithm("Linear")
    ops.analysis("Transient")

    u_rec = np.zeros((n_steps, 2), dtype=float)
    for step in range(n_steps):
        ok = ops.analyze(1, dt)
        assert ok == 0, f"Lamb half-space analysis failed at step {step} in mode={mode}"
        u_rec[step, 0] = ops.nodeDisp(rec_node, 1)
        u_rec[step, 1] = ops.nodeDisp(rec_node, 2)

    return u_rec


def test_anchor_4_lambs_problem_rayleigh_wave_cpml_vs_asd():
    """Anchor 4: Lamb's problem head-to-head (LadrunoCPML vs. ASDAbsorbingBoundary2D)."""
    u_ref = _run_lamb_halfspace("reference")
    u_asd = _run_lamb_halfspace("asd")
    u_cpml = _run_lamb_halfspace("cpml")

    l2_cpml, peak_cpml = _rel_l2_and_peak_error(u_cpml, u_ref)
    l2_asd, peak_asd = _rel_l2_and_peak_error(u_asd, u_ref)

    ratio_l2 = l2_asd / max(l2_cpml, 1.0e-12)
    ratio_peak = peak_asd / max(peak_cpml, 1.0e-12)
    print(
        f"[Lamb (6.5, 0)] CPML L2={l2_cpml:.4f}, peak={peak_cpml:.4f} | "
        f"ASD L2={l2_asd:.4f}, peak={peak_asd:.4f} | "
        f"Improvement: L2={ratio_l2:.2f}x, peak={ratio_peak:.2f}x"
    )

    assert l2_cpml <= 0.03, f"Lamb CPML L2 error {l2_cpml:.4f} > 3.0%"
    assert peak_cpml <= 0.03, f"Lamb CPML peak error {peak_cpml:.4f} > 3.0%"
    assert ratio_l2 >= 4.0, (
        f"Lamb CPML L2 improvement over ASDAbsorbingBoundary2D is {ratio_l2:.2f}x < 4.0x"
    )
    assert ratio_peak >= 4.0, (
        f"Lamb CPML peak improvement over ASDAbsorbingBoundary2D is {ratio_peak:.2f}x < 4.0x"
    )
