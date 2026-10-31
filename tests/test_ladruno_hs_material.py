# ------------------------------------------------------------------------------
# Ladruno OpenSees fork — file added by the Ladruno project
# Author:    Nicolas Mora Bowen — APE Ingenieria (ladruno.com)
# Provenance: Original fork-authored source; not present in upstream OpenSees.
# License:   See root LICENSE.md (upstream academic/non-commercial terms apply
#            to upstream-derived portions; fork additions are (c) the author).
# ------------------------------------------------------------------------------
"""Zone-A OpenSeesPy + Tcl verification battery for LadrunoHS / LadrunoHSSmall (WP-210, classTag 33029).

Validates the full 4-tier Hardening Soil (Schanz, Vermeer & Bonnier 1999) + HSsmall
(Benz 2007) + Matsuoka-Nakai smooth Lode cross-section (Panteghini & Lagioia 2014)
campaign through the compiled OpenSees binaries (dist/bin/opensees.pyd and dist/bin/OpenSees.exe):

  * Tier 1 (T0m): Parser fail-closed guards, wild-increment refusal & WP-99 commit latch,
                  all response channels, and 3D prescribed-strain parity vs Python oracle.
  * Tier 2 (T1):  Hostun sand virgin oedometer (Eoed_ref, K0nc = 1 - sin(phi)) + unload,
                  CIDC drained triaxial compression (mc == mn parity in TXC, Rowe dilatancy +
                  emax cutoff), HSsmall cyclic simple shear (G0 reset on reversal), and
                  Python vs OpenSees.exe Tcl parity.
  * Tier 3 (T2):  Plane-strain biaxial & strip-footing BVP demonstrating that -lodeShape mn
                  strictly bounds -lodeShape mc from above (while coinciding in TXC/TXE),
                  stage 0 -> 1 (updateMaterialStage) gravity-to-plastic transition, -implex
                  extrapolation, all 5 dimensional views (ThreeDimensional, PlaneStrain,
                  AxiSymmetric, PlaneStress, PlateFiber), and 2-round FileDatastore
                  sendSelf/recvSelf checkpoint/restore.
"""

from __future__ import annotations

import math
import os
import subprocess
import sys
import tempfile

import numpy as np
import pytest

from _testbed import ops
from _testbed.roundtrip import database_roundtrip, setresponse_smoke

sys.path.insert(0, os.path.dirname(__file__))
import ladruno_hs_reference as ref  # noqa: E402

pytestmark = [pytest.mark.zone_a]

_CUBE = {
    1: (0.0, 0.0, 0.0),
    2: (1.0, 0.0, 0.0),
    3: (1.0, 1.0, 0.0),
    4: (0.0, 1.0, 0.0),
    5: (0.0, 0.0, 1.0),
    6: (1.0, 0.0, 1.0),
    7: (1.0, 1.0, 1.0),
    8: (0.0, 1.0, 1.0),
}

_REGIME_NAME_TO_CODE = {
    "elastic": 0,
    "shear": 1,
    "failure": 2,
    "cap": 3,
    "shear_cap": 4,
    "failure_cap": 5,
    "tension": 6,
}


def _hostun_args(tag: int = 1, **kw) -> list:
    """Build argument list for nDMaterial LadrunoHS / LadrunoHSSmall."""
    cmd = kw.pop("cmd", "LadrunoHS")
    E50 = kw.pop("E50_ref", 35000.0)
    Eoed = kw.pop("Eoed_ref", 35000.0)
    Eur = kw.pop("Eur_ref", 105000.0)
    m = kw.pop("m", 0.5)
    c = kw.pop("c", 0.0)
    phi = kw.pop("phi", 42.0)
    psi = kw.pop("psi", 12.0)
    args = [cmd, tag, E50, Eoed, Eur, m, c, phi, psi]
    if cmd == "LadrunoHSSmall":
        g0, gam07 = kw.pop("smallStrain")
        args += [g0, gam07]
    elif "smallStrain" in kw:
        g0, gam07 = kw.pop("smallStrain")
        args += ["-smallStrain", g0, gam07]
    if "nu" in kw:
        args += ["-nu", kw.pop("nu")]
    if "pref" in kw:
        args += ["-pref", kw.pop("pref")]
    if "Rf" in kw:
        args += ["-Rf", kw.pop("Rf")]
    if "K0nc" in kw:
        args += ["-K0nc", kw.pop("K0nc")]
    if "sigmat" in kw:
        args += ["-sigmat", kw.pop("sigmat")]
    if "OCR" in kw:
        args += ["-OCR", kw.pop("OCR")]
    if "pop" in kw:
        args += ["-pop", kw.pop("pop")]
    if "e0" in kw:
        args += ["-e0", kw.pop("e0")]
    if "emax" in kw:
        args += ["-emax", kw.pop("emax")]
    if "rho" in kw:
        args += ["-rho", kw.pop("rho")]
    if "lodeShape" in kw:
        args += ["-lodeShape", kw.pop("lodeShape")]
    if "stiffStress" in kw:
        args += ["-stiffStress", kw.pop("stiffStress")]
    if "alphaCap" in kw:
        args += ["-alphaCap", kw.pop("alphaCap")]
    if "betaCap" in kw:
        args += ["-betaCap", kw.pop("betaCap")]
    if "sig0" in kw:
        args += ["-sig0", kw.pop("sig0")]
    if "sig03" in kw:
        s11, s22, s33 = kw.pop("sig03")
        args += ["-sig03", s11, s22, s33]
    if "stage" in kw:
        args += ["-stage", kw.pop("stage")]
    if kw.pop("implex", False):
        args += ["-implex"]
    assert not kw, f"Unconsumed keyword args: {kw}"
    return args


def _build_prescribed_cube(mat_args: list):
    """Single stdBrick unit cube with all 8 nodes prescribed (homogeneous strain probe)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for n, xyz in _CUBE.items():
        ops.node(n, *xyz)
        ops.fix(n, 1, 1, 1)
    ops.nDMaterial(*mat_args)
    ops.element("stdBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, mat_args[1])
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.test("NormDispIncr", 1.0e-10, 20, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")


def _impose_cube_strain_open(eps6_open: np.ndarray) -> int:
    """Prescribe nodal displacements on the unit cube for OpenSees (tension-positive) Voigt strain."""
    exx, eyy, ezz, gxy, gyz, gxz = [float(v) for v in eps6_open]
    for n, (x, y, z) in _CUBE.items():
        ux = exx * x + 0.5 * gxy * y + 0.5 * gxz * z
        uy = 0.5 * gxy * x + eyy * y + 0.5 * gyz * z
        uz = 0.5 * gxz * x + 0.5 * gyz * y + ezz * z
        ops.setNodeDisp(n, 1, ux, "-commit")
        ops.setNodeDisp(n, 2, uy, "-commit")
        ops.setNodeDisp(n, 3, uz, "-commit")
    return ops.analyze(1)


# ==============================================================================
# Tier 1 (T0m): Parser Guards, Wild-Increment Refusal, Responses & Oracle Parity
# ==============================================================================


def test_parser_fail_closed_guards():
    """Verify parser rejects invalid parameters, unknown flags, and malformed options."""
    bad_cases = [
        # Missing required args
        ["LadrunoHS", 1, 35000.0, 35000.0],
        # Non-positive E50_ref
        _hostun_args(1, E50_ref=0.0),
        # Eur_ref < E50_ref
        _hostun_args(1, E50_ref=35000.0, Eur_ref=20000.0),
        # m outside [0, 1]
        _hostun_args(1, m=1.2),
        _hostun_args(1, m=-0.1),
        # phi out of (0, 89)
        _hostun_args(1, phi=0.0),
        _hostun_args(1, phi=89.5),
        # nu out of (-1, 0.5)
        _hostun_args(1, nu=0.5),
        # Rf out of (0, 1]
        _hostun_args(1, Rf=0.0),
        _hostun_args(1, Rf=1.05),
        # OCR < 1
        _hostun_args(1, OCR=0.8),
        # Unknown flag
        _hostun_args(1) + ["-bogusFlag", 1.0],
        # Unknown lodeShape
        _hostun_args(1, lodeShape="vonMises"),
        # Unknown stiffStress
        _hostun_args(1, stiffStress="octahedral"),
        # LadrunoHSSmall without G0_ref / gamma_07
        ["LadrunoHSSmall", 1, 35000.0, 35000.0, 105000.0, 0.5, 0.0, 42.0, 12.0],
        # -smallStrain with non-positive G0_ref
        _hostun_args(1, smallStrain=(0.0, 1.5e-4)),
        # Negative -sig0
        _hostun_args(1, sig0=-100.0),
    ]
    for case in bad_cases:
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 3)
        with pytest.raises(Exception):
            ops.nDMaterial(*case)


def test_prescribed_cube_matches_python_oracle_mc_and_mn():
    """Drive a 3D stdBrick through general 3D compression + shear and compare GP1 to Python oracle."""
    for shape in ("mc", "mn"):
        p = ref.HSParams(
            E50_ref=35000.0,
            Eoed_ref=35000.0,
            Eur_ref=105000.0,
            m=0.5,
            c=5.0,
            phi_deg=38.0,
            psi_deg=8.0,
            nu_ur=0.2,
            p_ref=100.0,
            Rf=0.9,
            OCR=1.2,
            e0=0.65,
            emax=0.85,
            lode_shape=shape,
            use_small_strain=True,
            G0_ref=140000.0,
            gamma_07=2.0e-4,
        )
        sig0_comp = np.array([100.0, 80.0, 120.0, 0.0, 0.0, 0.0])
        st = ref.initialize_state_from_stress(p, sig0_comp)

        _build_prescribed_cube(
            _hostun_args(
                1,
                E50_ref=35000.0,
                Eoed_ref=35000.0,
                Eur_ref=105000.0,
                m=0.5,
                c=5.0,
                phi=38.0,
                psi=8.0,
                nu=0.2,
                pref=100.0,
                Rf=0.9,
                OCR=1.2,
                e0=0.65,
                emax=0.85,
                lodeShape=shape,
                smallStrain=(140000.0, 2.0e-4),
                sig03=(100.0, 80.0, 120.0),
            )
        )

        # Check all response channels before stepping
        setresponse_smoke(
            1,
            [
                ("material", 1, "stress"),
                ("material", 1, "strain"),
                ("material", 1, "tangent"),
                ("material", 1, "state"),
                ("material", 1, "gamma_p"),
                ("material", 1, "pp"),
                ("material", 1, "voidRatio"),
                ("material", 1, "regime"),
                ("material", 1, "capParams"),
            ],
        )

        cap_ops = np.array(ops.eleResponse(1, "material", 1, "capParams"))
        assert abs(cap_ops[0] - p.alpha_cap) < 1e-12
        assert abs(cap_ops[1] - p.beta_cap) < 1e-12

        # Multi-step general 3D geotechnical (compression-positive) strain path
        deps_comp = np.array([5.0e-4, -1.5e-4, -1.0e-4, 6.0e-4, -2.0e-4, 1.5e-4])
        eps_comp = np.zeros(6)
        for _ in range(8):
            eps_comp = eps_comp + deps_comp
            st, _, regime_name = ref.return_map_3d(p, st, eps_comp)
            assert _impose_cube_strain_open(-eps_comp) == 0

            # OpenSees stress is tension-positive (-st.sig)
            sig_ops = np.array(ops.eleResponse(1, "material", 1, "stress"))
            state_ops = np.array(ops.eleResponse(1, "material", 1, "state"))
            np.testing.assert_allclose(sig_ops, -st.sig, rtol=1e-11, atol=1e-10)
            assert abs(state_ops[3] - st.gamma_p) < 1e-12
            assert abs(state_ops[4] - st.p_p) < 1e-10
            assert abs(state_ops[5] - st.void_ratio) < 1e-12
            assert int(round(state_ops[6])) == _REGIME_NAME_TO_CODE[regime_name]


def test_wild_strain_increment_refusal_and_commit_latch():
    """Verify a wild strain increment (||deps|| > 5.0) refuses setTrialStrain and commitState."""
    _build_prescribed_cube(_hostun_args(1, sig0=100.0))
    # First valid step (eps_xx_open = -1.0e-3 -> compressive +1.0e-3)
    assert _impose_cube_strain_open(np.array([-1.0e-3, 0.0, 0.0, 0.0, 0.0, 0.0])) == 0
    sig_committed = np.array(ops.eleResponse(1, "material", 1, "stress"))

    # Impose a wild strain jump (deps_xx = -10.0 > 5.0 limit) -> analyze(1) must fail (< 0)
    rc = _impose_cube_strain_open(np.array([-10.0, 0.0, 0.0, 0.0, 0.0, 0.0]))
    assert rc < 0

    # Restore valid displacement and re-solve -> matches committed state
    assert _impose_cube_strain_open(np.array([-1.0e-3, 0.0, 0.0, 0.0, 0.0, 0.0])) == 0
    sig_after = np.array(ops.eleResponse(1, "material", 1, "stress"))
    np.testing.assert_allclose(sig_after, sig_committed, rtol=1e-12, atol=1e-12)


# ==============================================================================
# Tier 2 (T1): Oedometer, Drained Triaxial (mc == mn), HSsmall Cyclic & Tcl Parity
# ==============================================================================


def test_oedometer_virgin_compression_k0nc_and_unload():
    """Single stdBrick in 1D oedometer compression: verifies K0nc = 1 - sin(phi) and Eoed_ref."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for n, (x, y, z) in _CUBE.items():
        ops.node(n, x, y, z)
        # Fix x and y everywhere (zero lateral strain), fix z=0 at bottom, free z at top (5,6,7,8)
        ops.fix(n, 1, 1, 1 if z == 0.0 else 0)

    # Start at K0nc state: compressive sig_v0 = 50 kPa, sig_h0 = K0nc * 50 kPa
    phi_deg = 35.0
    k0nc = 1.0 - math.sin(math.radians(phi_deg))
    sig_v0 = 50.0
    sig_h0 = k0nc * sig_v0

    ops.nDMaterial(
        *_hostun_args(
            1,
            E50_ref=30000.0,
            Eoed_ref=30000.0,
            Eur_ref=90000.0,
            m=0.5,
            c=0.0,
            phi=phi_deg,
            psi=5.0,
            nu=0.2,
            pref=100.0,
            K0nc=k0nc,
            Rf=0.9,
            OCR=1.0,
            sig03=(sig_h0, sig_h0, sig_v0),
        )
    )
    ops.element("stdBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)

    for n in (6, 7, 8):
        ops.equalDOF(5, n, 3)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(5, 0.0, 0.0, -1.0)

    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.test("NormDispIncr", 1.0e-10, 30, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")

    # Compress axially by eps_zz = -0.008 in 40 steps
    ops.integrator("DisplacementControl", 5, 3, -0.008 / 40)
    for _ in range(40):
        assert ops.analyze(1) == 0

    sig = np.array(ops.eleResponse(1, "material", 1, "stress"))
    sig_h, sig_v = -sig[0], -sig[2]
    # Along virgin oedometer compression, K0 = sig_h / sig_v stays close to 1 - sin(phi)
    assert abs((sig_h / sig_v) - k0nc) / k0nc < 0.04
    pp_after_load = ops.eleResponse(1, "material", 1, "pp")[0]
    assert pp_after_load > 180.0

    # Unload by +0.001 in 10 steps -> purely elastic unload (pp unchanged)
    ops.integrator("DisplacementControl", 5, 3, +0.001 / 10)
    for _ in range(10):
        assert ops.analyze(1) == 0
    pp_after_unload = ops.eleResponse(1, "material", 1, "pp")[0]
    assert abs(pp_after_unload - pp_after_load) < 1e-9


def test_drained_triaxial_compression_mc_equals_mn_and_tcl_parity():
    """Drained triaxial compression (CIDC) on a 3D stdBrick: verifies mc == mn parity and Tcl parity."""
    def run_cidc_py(lode_shape: str):
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 3)
        for n, (x, y, z) in _CUBE.items():
            ops.node(n, x, y, z)
        ops.fix(1, 1, 1, 1)
        ops.fix(2, 0, 1, 1)
        ops.fix(3, 0, 0, 1)
        ops.fix(4, 1, 0, 1)
        ops.fix(5, 1, 1, 0)
        ops.fix(6, 0, 1, 0)
        ops.fix(8, 1, 0, 0)
        for n in (3, 6, 7):
            ops.equalDOF(2, n, 1)
        for n in (3, 7, 8):
            ops.equalDOF(4, n, 2)
        for n in (6, 7, 8):
            ops.equalDOF(5, n, 3)

        ops.nDMaterial(
            *_hostun_args(
                1,
                E50_ref=35000.0,
                Eoed_ref=35000.0,
                Eur_ref=105000.0,
                m=0.5,
                c=0.0,
                phi=42.0,
                psi=12.0,
                nu=0.2,
                pref=100.0,
                Rf=0.9,
                OCR=1.5,
                e0=0.65,
                emax=0.653,
                lodeShape=lode_shape,
                sig0=100.0,
            )
        )
        ops.element("stdBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)

        # Constant lateral confining load (-100 kPa on +x and +y faces) held by Constant timeSeries
        ops.timeSeries("Constant", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(2, -100.0, 0.0, 0.0)
        ops.load(4, 0.0, -100.0, 0.0)

        # Axial reference load on node 5 z-DOF for DisplacementControl
        ops.timeSeries("Linear", 2)
        ops.pattern("Plain", 2, 2)
        ops.load(5, 0.0, 0.0, -1.0)

        ops.system("FullGeneral")
        ops.numberer("Plain")
        ops.constraints("Plain")
        ops.test("NormDispIncr", 1.0e-9, 40, 0)
        ops.algorithm("Newton")
        ops.integrator("DisplacementControl", 5, 3, -0.04 / 40)
        ops.analysis("Static")

        q_hist = []
        ev_hist = []
        e_hist = []
        for _ in range(40):
            assert ops.analyze(1) == 0
            sig = np.array(ops.eleResponse(1, "material", 1, "stress"))
            eps = np.array(ops.eleResponse(1, "material", 1, "strain"))
            vr = ops.eleResponse(1, "material", 1, "voidRatio")[0]
            q_hist.append(sig[0] - sig[2])  # (-sig3) - (-sig1) = sig1_comp - sig3_comp > 0
            ev_hist.append(-(eps[0] + eps[1] + eps[2]))  # geotechnical volumetric strain
            e_hist.append(vr)
        return np.array(q_hist), np.array(ev_hist), np.array(e_hist)

    q_mc, ev_mc, e_mc = run_cidc_py("mc")
    q_mn, ev_mn, e_mn = run_cidc_py("mn")

    # 1. Exact TXC parity between Mohr-Coulomb (-lodeShape mc) and Matsuoka-Nakai (-lodeShape mn)
    np.testing.assert_allclose(q_mc, q_mn, rtol=1e-12, atol=1e-10)
    np.testing.assert_allclose(ev_mc, ev_mn, rtol=1e-12, atol=1e-12)

    # 2. Hyperbolic mobilization toward q_f = 100 * (2 sin 42) / (1 - sin 42) = 404.4 kPa
    assert 330.0 < q_mc[-1] < 405.0
    # 3. Rowe stress-dilatancy: initially contractive (ev > 0), then dilative, then capped at emax=0.653
    assert np.max(ev_mc[:5]) > 0.0
    assert e_mc[-1] >= 0.653 - 1e-4
    assert abs(ev_mc[-1] - ev_mc[-2]) < 1.5e-4

    # 4. Tcl parity check via dist/bin/OpenSees.exe subprocess
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    exe_path = os.path.join(repo_root, "dist", "bin", "OpenSees.exe")
    if os.path.isfile(exe_path):
        tcl_script = """
model basic -ndm 3 -ndf 3
node 1 0.0 0.0 0.0; node 2 1.0 0.0 0.0; node 3 1.0 1.0 0.0; node 4 0.0 1.0 0.0
node 5 0.0 0.0 1.0; node 6 1.0 0.0 1.0; node 7 1.0 1.0 1.0; node 8 0.0 1.0 1.0
fix 1 1 1 1; fix 2 0 1 1; fix 3 0 0 1; fix 4 1 0 1
fix 5 1 1 0; fix 6 0 1 0; fix 8 1 0 0
equalDOF 2 3 1; equalDOF 2 6 1; equalDOF 2 7 1
equalDOF 4 3 2; equalDOF 4 7 2; equalDOF 4 8 2
equalDOF 5 6 3; equalDOF 5 7 3; equalDOF 5 8 3
nDMaterial LadrunoHS 1 35000.0 35000.0 105000.0 0.5 0.0 42.0 12.0 -nu 0.2 -pref 100.0 -Rf 0.9 -OCR 1.5 -e0 0.65 -emax 0.653 -lodeShape mn -sig0 100.0
element stdBrick 1 1 2 3 4 5 6 7 8 1
timeSeries Constant 1
pattern Plain 1 1 {
    load 2 -100.0 0.0 0.0
    load 4 0.0 -100.0 0.0
}
timeSeries Linear 2
pattern Plain 2 2 {
    load 5 0.0 0.0 -1.0
}
system FullGeneral
numberer Plain
constraints Plain
test NormDispIncr 1.0e-9 40 0
algorithm Newton
integrator DisplacementControl 5 3 -0.001
analysis Static
analyze 40
set s [eleResponse 1 material 1 stress]
puts "TCL_Q_FINAL [expr {[lindex $s 0] - [lindex $s 2]}]"
exit
"""
        with tempfile.TemporaryDirectory() as td:
            tcl_file = os.path.join(td, "cidc.tcl")
            with open(tcl_file, "w", encoding="utf-8") as f:
                f.write(tcl_script)
            proc = subprocess.run(
                [exe_path, tcl_file],
                stdin=subprocess.DEVNULL,
                capture_output=True,
                text=True,
                timeout=30,
                check=True,
            )
            out = proc.stdout + proc.stderr
            q_tcl = None
            for line in out.splitlines():
                if line.startswith("TCL_Q_FINAL"):
                    q_tcl = float(line.split()[1])
            assert q_tcl is not None
            assert abs(q_tcl - q_mn[-1]) < 1e-9


def test_hssmall_alias_cyclic_simple_shear_reversal():
    """Verify LadrunoHSSmall degrades shear tangent inside the elastic domain and resets to G0 on reversal."""
    _build_prescribed_cube(
        _hostun_args(
            1,
            cmd="LadrunoHSSmall",
            OCR=2.5,
            smallStrain=(140000.0, 1.5e-4),
            stiffStress="p",
            sig0=100.0,
        )
    )

    # First pre-expand the shear yield cone to gamma_xy = +5.0e-3 so subsequent unload/reload is inside the cone
    assert _impose_cube_strain_open(np.array([0.0, 0.0, 0.0, 5.0e-3, 0.0, 0.0])) == 0

    # Unload from +5.0e-3 down to +3.5e-3 in 15 steps (inside the shear cone -> Regime 0 elastic HSsmall)
    g_tan_unload = []
    for k in range(1, 16):
        eps_open = np.array([0.0, 0.0, 0.0, 5.0e-3 - 1.0e-4 * k, 0.0, 0.0])
        assert _impose_cube_strain_open(eps_open) == 0
        tan = np.array(ops.eleResponse(1, "material", 1, "tangent")).reshape((6, 6))
        g_tan_unload.append(tan[3, 3])

    # Immediately upon reversal (k=1, dgamma = 1.0e-4), tangent shear modulus resets near G0,
    # then degrades monotonically down to Gur(p) = 43,750 * (p / 100)^0.5 kPa
    sig_end = np.array(ops.eleResponse(1, "material", 1, "stress"))
    p_end = -float(np.mean(sig_end[:3]))
    gur_expected = 43750.0 * (p_end / 100.0) ** 0.5
    assert g_tan_unload[0] > 70000.0
    assert abs(g_tan_unload[-1] - gur_expected) < 1e-6
    assert all(g_tan_unload[i] >= g_tan_unload[i + 1] - 1e-9 for i in range(len(g_tan_unload) - 1))

    # Second reversal (+3.5e-3 -> +3.5e-3 + 1.0e-5) -> Benz (2007) history tensor resets tangent near G0 again!
    eps_rev2 = np.array([0.0, 0.0, 0.0, 3.5e-3 + 1.0e-5, 0.0, 0.0])
    assert _impose_cube_strain_open(eps_rev2) == 0
    tan_rev2 = np.array(ops.eleResponse(1, "material", 1, "tangent")).reshape((6, 6))
    assert tan_rev2[3, 3] > 120000.0


# ==============================================================================
# Tier 3 (T2): Plane-Strain Bounding (mn > mc), Multi-Dim Views, Stage & DB
# ==============================================================================


def test_plane_strain_biaxial_and_footing_mn_strictly_bounds_mc():
    """Verify User Decision 1: in PlaneStrain, -lodeShape mn strictly bounds -lodeShape mc from above."""
    def run_plane_strain_cell(lode_shape: str) -> float:
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 2)
        ops.node(1, 0.0, 0.0)
        ops.node(2, 1.0, 0.0)
        ops.node(3, 1.0, 1.0)
        ops.node(4, 0.0, 1.0)
        ops.fix(1, 1, 1)
        ops.fix(2, 0, 1)
        ops.fix(4, 1, 0)
        ops.equalDOF(2, 3, 1)
        ops.equalDOF(4, 3, 2)

        ops.nDMaterial(
            *_hostun_args(
                1,
                E50_ref=35000.0,
                Eoed_ref=35000.0,
                Eur_ref=105000.0,
                m=0.5,
                c=5.0,
                phi=38.0,
                psi=8.0,
                nu=0.2,
                pref=100.0,
                Rf=0.9,
                OCR=1.5,
                lodeShape=lode_shape,
                sig0=100.0,
            )
        )
        ops.element("quad", 1, 1, 2, 3, 4, 1.0, "PlaneStrain", 1)

        ops.timeSeries("Constant", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(2, -100.0, 0.0)

        ops.timeSeries("Linear", 2)
        ops.pattern("Plain", 2, 2)
        ops.load(4, 0.0, -1.0)

        ops.system("FullGeneral")
        ops.numberer("Plain")
        ops.constraints("Plain")
        ops.test("NormDispIncr", 1.0e-9, 40, 0)
        ops.algorithm("Newton")
        ops.integrator("DisplacementControl", 4, 2, -0.025 / 25)
        ops.analysis("Static")

        for _ in range(25):
            assert ops.analyze(1) == 0
        sig = np.array(ops.eleResponse(1, "material", 1, "stress"))
        return float(-sig[1])

    sig_v_mc = run_plane_strain_cell("mc")
    sig_v_mn = run_plane_strain_cell("mn")
    # In plane strain (intermediate principal stress sigma_2 between sigma_1 and sigma_3),
    # Matsuoka-Nakai (-lodeShape mn) strictly bounds Mohr-Coulomb (-lodeShape mc) from above
    assert sig_v_mn > sig_v_mc * 1.05
    assert sig_v_mn < sig_v_mc * 1.25

    # Multi-element plane-strain strip-footing BVP (4x2 quad mesh, footing on left half x in [0, 1]):
    def run_strip_footing_bvp(lode_shape: str) -> float:
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 2)
        nx, ny = 4, 2
        xs = np.linspace(0.0, 2.0, nx + 1)
        ys = np.linspace(0.0, 1.0, ny + 1)
        for j, y in enumerate(ys):
            for i, x in enumerate(xs):
                nid = j * (nx + 1) + i + 1
                ops.node(nid, float(x), float(y))
                if j == 0:
                    ops.fix(nid, 1, 1)
                elif i == 0 or i == nx:
                    ops.fix(nid, 1, 0)
        ops.nDMaterial(
            *_hostun_args(
                1,
                E50_ref=35000.0,
                Eoed_ref=35000.0,
                Eur_ref=105000.0,
                m=0.5,
                c=10.0,
                phi=38.0,
                psi=8.0,
                nu=0.2,
                pref=100.0,
                Rf=0.9,
                OCR=1.5,
                lodeShape=lode_shape,
                sig0=100.0,
            )
        )
        eid = 1
        for j in range(ny):
            for i in range(nx):
                n1 = j * (nx + 1) + i + 1
                n2 = n1 + 1
                n3 = n2 + (nx + 1)
                n4 = n1 + (nx + 1)
                ops.element("quad", eid, n1, n2, n3, n4, 1.0, "PlaneStrain", 1)
                eid += 1
        top_left = ny * (nx + 1) + 1
        top_mid1 = top_left + 1
        top_mid2 = top_left + 2
        ops.equalDOF(top_left, top_mid1, 2)
        ops.equalDOF(top_left, top_mid2, 2)
        # Constant top-surface overburden traction balancing -sig0 100.0 kPa on the 2.0 m wide mesh (dx = 0.5 m):
        ops.timeSeries("Constant", 10)
        ops.pattern("Plain", 10, 10)
        for i in range(nx + 1):
            nid = ny * (nx + 1) + i + 1
            trib_w = 0.25 if (i == 0 or i == nx) else 0.50
            ops.load(nid, 0.0, -100.0 * trib_w)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(top_left, 0.0, -1.0)
        ops.system("FullGeneral")
        ops.numberer("Plain")
        ops.constraints("Plain")
        ops.test("NormDispIncr", 1.0e-8, 40, 0)
        ops.algorithm("Newton")
        ops.integrator("DisplacementControl", top_left, 2, -0.015 / 30)
        ops.analysis("Static")
        for _ in range(30):
            assert ops.analyze(1) == 0
        return float(ops.getLoadFactor(1))

    r_mc = run_strip_footing_bvp("mc")
    r_mn = run_strip_footing_bvp("mn")
    assert r_mn > r_mc * 1.02


def test_stage_transition_and_implex_mode():
    """Verify -stage 0 (elastic gravity) -> updateMaterialStage -stage 1 (plastic), -implex, and user cap preservation across multiple elements."""
    for use_implex in (False, True):
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 2)
        node_coords = ((1, (0.0, 0.0)), (2, (1.0, 0.0)), (3, (1.0, 1.0)), (4, (0.0, 1.0)), (5, (2.0, 0.0)), (6, (2.0, 1.0)))
        for n, (x, y) in node_coords:
            ops.node(n, x, y)
            ops.fix(n, 1, 1)

        ops.nDMaterial(
            *_hostun_args(
                1,
                cmd="LadrunoHSSmall",
                smallStrain=(140000.0, 1.5e-4),
                alphaCap=0.85,
                betaCap=1.2e-4,
                stage=0,
                implex=use_implex,
                sig0=100.0,
            )
        )
        ops.element("quad", 1, 1, 2, 3, 4, 1.0, "PlaneStrain", 1)
        ops.element("quad", 2, 2, 5, 6, 3, 1.0, "PlaneStrain", 1)
        ops.system("FullGeneral")
        ops.numberer("Plain")
        ops.constraints("Plain")
        ops.test("NormDispIncr", 1.0e-10, 20, 0)
        ops.algorithm("Newton")
        ops.integrator("LoadControl", 1.0)
        ops.analysis("Static")

        # Step 1 in stage 0 (elastic): gamma_p must remain 0.0 under elastic shear (tau = 21.875 kPa < tau_f)
        for n, (x, y) in node_coords:
            ops.setNodeDisp(n, 1, 0.0005 * y, "-commit")
        assert ops.analyze(1) == 0
        assert ops.eleResponse(1, "material", 1, "gamma_p")[0] == 0.0
        assert ops.eleResponse(2, "material", 1, "gamma_p")[0] == 0.0
        sig_stage0 = np.array(ops.eleResponse(1, "material", 1, "stress"))

        # Switch to stage 1 (plastic) via updateMaterialStage (must flip both element 1 and element 2!)
        ops.updateMaterialStage("-material", 1, "-stage", 1)
        # First small step in stage 1 (dgamma = 2e-5) must increment from end-of-stage-0 strain (no jump!)
        for n, (x, y) in node_coords:
            ops.setNodeDisp(n, 1, (0.0005 + 2.0e-5) * y, "-commit")
        assert ops.analyze(1) == 0
        sig_step1 = np.array(ops.eleResponse(1, "material", 1, "stress"))
        assert abs(sig_step1[2] - sig_stage0[2]) < 5.0

        for k in range(2, 6):
            for n, (x, y) in node_coords:
                ops.setNodeDisp(n, 1, 0.002 * k * y, "-commit")
            assert ops.analyze(1) == 0
        assert ops.eleResponse(1, "material", 1, "gamma_p")[0] > 1.0e-3
        assert ops.eleResponse(2, "material", 1, "gamma_p")[0] > 1.0e-3

        # Verify adding a new element (eid=3) in Stage 1 via prototype getCopy() adopts stage 1
        # and does NOT clobber existing elements back to stage 0 ([M-2]):
        ops.node(7, 3.0, 0.0)
        ops.node(8, 3.0, 1.0)
        ops.fix(7, 1, 1)
        ops.fix(8, 1, 1)
        ops.element("quad", 3, 5, 7, 8, 6, 1.0, "PlaneStrain", 1)
        for n, (x, y) in (*node_coords, (7, (3.0, 0.0)), (8, (3.0, 1.0))):
            ops.setNodeDisp(n, 1, 0.012 * y, "-commit")
        assert ops.analyze(1) == 0
        assert ops.eleResponse(3, "material", 1, "gamma_p")[0] > 1.0e-3

        # Verify parameterID == 13 ("stage" / "materialState") also updates ladrunoHSTagStageMap ([M-1]):
        ops.parameter(20, "element", 1, "material", 1, "stage")
        ops.updateParameter(20, 0)
        ops.updateParameter(20, 1)
        for n, (x, y) in (*node_coords, (7, (3.0, 0.0)), (8, (3.0, 1.0))):
            ops.setNodeDisp(n, 1, 0.015 * y, "-commit")
        assert ops.analyze(1) == 0
        assert ops.eleResponse(1, "material", 1, "gamma_p")[0] > 1.5e-3
        assert ops.eleResponse(2, "material", 1, "gamma_p")[0] > 1.5e-3

        # Verify parameter update preserves explicit user -alphaCap / -betaCap:
        ops.parameter(10, "element", 1, "material", 1, "E50")
        ops.updateParameter(10, 40000.0)
        cap_after = np.array(ops.eleResponse(1, "material", 1, "capParams"))
        assert abs(cap_after[0] - 0.85) < 1e-12
        assert abs(cap_after[1] - 1.2e-4) < 1e-15


def test_reduced_dimensional_views_planestress_and_platefiber():
    """Verify PlaneStress and PlateFiber dimensional projections with -sig0 and non-trivial plastic response."""
    # 1. PlaneStress (order 3: [exx, eyy, gxy] with condensed sig_zz = 0)
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(1, 0.0, 0.0)
    ops.node(2, 1.0, 0.0)
    ops.node(3, 1.0, 1.0)
    ops.node(4, 0.0, 1.0)
    for n in (1, 2, 3, 4):
        ops.fix(n, 1, 1)
    ops.nDMaterial(*_hostun_args(1, c=25.0, sigmat=20.0, OCR=1.2, sig0=50.0))
    ops.element("quad", 1, 1, 2, 3, 4, 1.0, "PlaneStress", 1)
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.test("NormDispIncr", 1e-10, 20, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    for k in range(1, 5):
        for n, (x, y) in ((1, (0, 0)), (2, (1, 0)), (3, (1, 1)), (4, (0, 1))):
            ops.setNodeDisp(n, 1, -0.001 * k * x + 0.0005 * k * y, "-commit")
        assert ops.analyze(1) == 0
    sig_ps = np.array(ops.eleResponse(1, "material", 1, "stress"))
    tan_ps = np.array(ops.eleResponse(1, "material", 1, "tangent")).reshape((3, 3))
    gam_ps = ops.eleResponse(1, "material", 1, "gamma_p")[0]
    assert len(sig_ps) == 3
    assert -sig_ps[0] > 60.0  # non-trivial compressive stress growth above initial 50 kPa
    assert gam_ps > 1.0e-4    # active plastic shear mobilization in PlaneStress
    assert np.all(np.diag(tan_ps) > 1000.0)

    # 2. PlateFiber via LayeredShellFiberSection + ASDShellQ4 (order 5)
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, 1.0, 0.0, 0.0)
    ops.node(3, 1.0, 1.0, 0.0)
    ops.node(4, 0.0, 1.0, 0.0)
    for n in (1, 2, 3, 4):
        ops.fix(n, 1, 1, 1, 1, 1, 1)
    ops.nDMaterial(*_hostun_args(1, c=25.0, sigmat=20.0, OCR=1.2, sig0=50.0))
    ops.section("LayeredShell", 1, 3, 1, 0.1, 1, 0.1, 1, 0.1)
    ops.element("ASDShellQ4", 1, 1, 2, 3, 4, 1)
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.test("NormDispIncr", 1e-10, 20, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    for k in range(1, 5):
        for n, (x, y) in ((1, (0, 0)), (2, (1, 0)), (3, (1, 1)), (4, (0, 1))):
            ops.setNodeDisp(n, 1, -0.001 * k * x + 0.0005 * k * y, "-commit")
        assert ops.analyze(1) == 0
    ops.reactions()
    r_shell_x = sum(abs(ops.nodeReaction(n, 1)) for n in (1, 2, 3, 4))
    assert r_shell_x > 10.0  # non-trivial shell membrane force from PlateFiber condensation


def test_database_roundtrip_two_rounds_mid_plastic_hssmall():
    """Verify 2-round FileDatastore sendSelf/recvSelf checkpoint/restore mid-plastic path."""
    def build_and_step():
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 3)
        for n, (x, y, z) in _CUBE.items():
            ops.node(n, x, y, z)
        ops.fix(1, 1, 1, 1)
        ops.fix(2, 0, 1, 1)
        ops.fix(3, 0, 0, 1)
        ops.fix(4, 1, 0, 1)
        ops.fix(5, 1, 1, 0)
        ops.fix(6, 0, 1, 0)
        ops.fix(8, 1, 0, 0)
        for n in (3, 6, 7):
            ops.equalDOF(2, n, 1)
        for n in (3, 7, 8):
            ops.equalDOF(4, n, 2)
        for n in (6, 7, 8):
            ops.equalDOF(5, n, 3)

        ops.nDMaterial(
            *_hostun_args(
                1,
                cmd="LadrunoHSSmall",
                OCR=1.2,
                lodeShape="mn",
                smallStrain=(140000.0, 1.5e-4),
                sig0=100.0,
            )
        )
        ops.element("stdBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
        ops.timeSeries("Constant", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(2, -100.0, 0.0, 0.0)
        ops.load(4, 0.0, -100.0, 0.0)
        ops.timeSeries("Linear", 2)
        ops.pattern("Plain", 2, 2)
        ops.load(5, 0.0, 0.0, -1.0)
        ops.system("FullGeneral")
        ops.numberer("Plain")
        ops.constraints("Plain")
        ops.test("NormDispIncr", 1.0e-9, 40, 0)
        ops.algorithm("Newton")
        ops.integrator("DisplacementControl", 5, 3, -0.002)
        ops.analysis("Static")
        assert ops.analyze(5) == 0

    def probe_state():
        return (
            list(ops.eleResponse(1, "material", 1, "stress"))
            + list(ops.eleResponse(1, "material", 1, "state"))
            + list(ops.eleResponse(1, "material", 1, "capParams"))
        )

    database_roundtrip(
        build_and_step,
        probe_nodes=[2, 4, 5, 7],
        ndf=3,
        dbname="ladruno_hs_rt",
        probe_fn=probe_state,
        rounds=2,
    )


def test_drained_triaxial_extension_and_non_auxetic_hssmall():
    """Verify:
      1. Drained Triaxial Extension (TXE, constant minor principal stress sig3 = 100 kPa on z-face,
         equal radial compression on x and y faces) gives exact mc == mn parity and matches the
         analytical constant-sig3 hyperbolic curve.
      2. Non-auxetic LadrunoHSSmall with G0_ref / Gur_ref = 5.0 preserves nu_eff == nu_ur = 0.20 > 0.
    """
    def run_cide_py(lode_shape: str):
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 3)
        for n, (x, y, z) in _CUBE.items():
            ops.node(n, x, y, z)
        ops.fix(1, 1, 1, 1)
        ops.fix(2, 0, 1, 1)
        ops.fix(3, 0, 0, 1)
        ops.fix(4, 1, 0, 1)
        ops.fix(5, 1, 1, 0)
        ops.fix(6, 0, 1, 0)
        ops.fix(8, 1, 0, 0)
        for n in (3, 6, 7):
            ops.equalDOF(2, n, 1)
        for n in (3, 7, 8):
            ops.equalDOF(4, n, 2)
        for n in (6, 7, 8):
            ops.equalDOF(5, n, 3)

        ops.nDMaterial(
            *_hostun_args(
                1,
                E50_ref=35000.0,
                Eoed_ref=35000.0,
                Eur_ref=105000.0,
                m=0.5,
                c=0.0,
                phi=42.0,
                psi=12.0,
                nu=0.2,
                pref=100.0,
                Rf=0.9,
                OCR=10.0,
                lodeShape=lode_shape,
                sig0=100.0,
            )
        )
        ops.element("stdBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)

        # Constant axial stress sig_3 = 100 kPa on +z face; equal radial compression on +x and +y faces:
        ops.timeSeries("Constant", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(5, 0.0, 0.0, -100.0)

        ops.timeSeries("Linear", 2)
        ops.pattern("Plain", 2, 2)
        ops.load(2, -1.0, 0.0, 0.0)
        ops.load(4, 0.0, -1.0, 0.0)

        ops.system("FullGeneral")
        ops.numberer("Plain")
        ops.constraints("Plain")
        ops.test("NormDispIncr", 1.0e-10, 40, 0)
        ops.algorithm("Newton")
        ops.integrator("DisplacementControl", 2, 1, -0.015 / 30)
        ops.analysis("Static")

        q_hist = []
        gam_hist = []
        for _ in range(30):
            assert ops.analyze(1) == 0
            sig = np.array(ops.eleResponse(1, "material", 1, "stress"))
            gam_p = ops.eleResponse(1, "material", 1, "gamma_p")[0]
            # In TXE: sig1_comp = sig2_comp = -sig[0], sig3_comp = -sig[2] = 100 kPa
            assert abs(-sig[2] - 100.0) < 1e-6
            assert abs(sig[0] - sig[1]) < 1e-8
            q_hist.append((-sig[0]) - (-sig[2]))
            gam_hist.append(gam_p)
        return np.array(q_hist), np.array(gam_hist)

    q_mc, gam_mc = run_cide_py("mc")
    q_mn, gam_mn = run_cide_py("mn")
    np.testing.assert_allclose(q_mc, q_mn, rtol=1e-11, atol=1e-10)
    np.testing.assert_allclose(gam_mc, gam_mn, rtol=1e-11, atol=1e-12)

    # Check exact constant-sig3 hyperbolic formula at step 15:
    sin_phi = math.sin(math.radians(42.0))
    qf = 100.0 * (2.0 * sin_phi / (1.0 - sin_phi))
    qa = qf / 0.9
    ei_ref = 2.0 * 35000.0 / (2.0 - 0.9)
    q_mid = q_mc[14]
    gam_ana = (2.0 / ei_ref) * (q_mid / (1.0 - q_mid / qa)) - (2.0 * q_mid / 105000.0)
    assert abs(gam_mc[14] - gam_ana) < 1e-5

    # Check non-auxetic LadrunoHSSmall tangent with G0_ref / Gur_ref = 5.0:
    # First mobilize the shear yield locus to eps_x = -2e-4, then apply a small elastic
    # unloading reversal (eps_x = -1.99e-4) where regime == 0 and Ge resets to G0_ref:
    _build_prescribed_cube(
        _hostun_args(
            1,
            cmd="LadrunoHSSmall",
            E50_ref=25000.0,
            Eoed_ref=25000.0,
            Eur_ref=72000.0,
            nu=0.20,
            OCR=20.0,
            smallStrain=(150000.0, 2.0e-4),
            sig0=100.0,
        )
    )
    assert _impose_cube_strain_open(np.array([-2.0e-4, 0.0, 0.0, 0.0, 0.0, 0.0])) == 0
    sig_pre = np.array(ops.eleResponse(1, "material", 1, "stress"))
    assert _impose_cube_strain_open(np.array([-1.999e-4, 0.0, 0.0, 0.0, 0.0, 0.0])) == 0
    assert int(ops.eleResponse(1, "material", 1, "regime")[0]) == 0
    sig_hss = np.array(ops.eleResponse(1, "material", 1, "stress"))
    dsig_hss = sig_hss - sig_pre
    nu_eff = dsig_hss[1] / (dsig_hss[0] + dsig_hss[1])
    assert abs(nu_eff - 0.20) < 1e-6
    ctan = np.array(ops.eleResponse(1, "material", 1, "tangent")).reshape(6, 6)
    assert ctan[0, 1] > 0.0
    assert abs(ctan[0, 1] / (ctan[0, 0] + ctan[0, 1]) - 0.20) < 1e-4


def test_implex_hssmall_and_condensed_planestress_and_responses():
    """Verify:
      1. LadrunoHSSmall -implex uses the non-auxetic small-strain secant modulus Ge and
         Ke = (Kur_ref/Gur_ref)*Ge (preserving nu_ur = 0.20 and SPD tangent), with zero
         spurious dsig_p drift across purely elastic steps.
      2. Recordable response aliases plasticStrain, void_ratio, G_sec, G_tan, and
         smallStrainModulus return the exact internal variables.
      3. Condensed PlaneStress mode enforces sig_zz == 0 and produces a consistent 3x3 tangent.
    """
    _build_prescribed_cube(
        _hostun_args(
            1,
            cmd="LadrunoHSSmall",
            E50_ref=25000.0,
            Eoed_ref=25000.0,
            Eur_ref=72000.0,
            nu=0.20,
            OCR=20.0,
            smallStrain=(150000.0, 2.0e-4),
            sig0=100.0,
            implex=True,
        )
    )
    # Mobilize shear cone first, then take two consecutive purely elastic unloading steps:
    assert _impose_cube_strain_open(np.array([-2.0e-4, 0.0, 0.0, 0.0, 0.0, 0.0])) == 0
    assert _impose_cube_strain_open(np.array([-1.999e-4, 0.0, 0.0, 0.0, 0.0, 0.0])) == 0
    assert int(ops.eleResponse(1, "material", 1, "regime")[0]) == 0

    ctan_imp = np.array(ops.eleResponse(1, "material", 1, "tangent")).reshape(6, 6)
    g_sec = ops.eleResponse(1, "material", 1, "G_sec")[0]
    g_tan = ops.eleResponse(1, "material", 1, "G_tan")[0]
    ss_mod = ops.eleResponse(1, "material", 1, "smallStrainModulus")
    assert abs(ss_mod[0] - g_sec) < 1e-12
    assert abs(ss_mod[1] - g_tan) < 1e-12
    assert g_sec > 140000.0  # Near G0_ref = 150000 kPa on elastic reversal (vs Gur ~ 31500 kPa)
    assert abs(ctan_imp[3, 3] - g_sec) / g_sec < 1e-4
    assert abs(ctan_imp[0, 1] / (ctan_imp[0, 0] + ctan_imp[0, 1]) - 0.20) < 1e-12

    # Second elastic unloading step: verify IMPL-EX stress matches implicit elastic stress
    # (confirming dsig_p was identically zero on the preceding elastic step!)
    assert _impose_cube_strain_open(np.array([-1.998e-4, 0.0, 0.0, 0.0, 0.0, 0.0])) == 0
    sig_imp2 = np.array(ops.eleResponse(1, "material", 1, "stress"))

    _build_prescribed_cube(
        _hostun_args(
            1,
            cmd="LadrunoHSSmall",
            E50_ref=25000.0,
            Eoed_ref=25000.0,
            Eur_ref=72000.0,
            nu=0.20,
            OCR=20.0,
            smallStrain=(150000.0, 2.0e-4),
            sig0=100.0,
            implex=False,
        )
    )
    assert _impose_cube_strain_open(np.array([-2.0e-4, 0.0, 0.0, 0.0, 0.0, 0.0])) == 0
    assert _impose_cube_strain_open(np.array([-1.999e-4, 0.0, 0.0, 0.0, 0.0, 0.0])) == 0
    assert _impose_cube_strain_open(np.array([-1.998e-4, 0.0, 0.0, 0.0, 0.0, 0.0])) == 0
    sig_std2 = np.array(ops.eleResponse(1, "material", 1, "stress"))
    np.testing.assert_allclose(sig_imp2, sig_std2, rtol=1e-12, atol=1e-10)

    # Check response aliases plasticStrain and void_ratio:
    gam_p1 = ops.eleResponse(1, "material", 1, "gamma_p")[0]
    gam_p2 = ops.eleResponse(1, "material", 1, "plasticStrain")[0]
    assert abs(gam_p1 - gam_p2) < 1e-15
    vr1 = ops.eleResponse(1, "material", 1, "voidRatio")[0]
    vr2 = ops.eleResponse(1, "material", 1, "void_ratio")[0]
    assert abs(vr1 - vr2) < 1e-15

    # Check 2D PlaneStress condensed mode central-FD 3x3 tangent consistency:
    def eval_planestress(eps_inplane: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 2)
        ops.node(1, 0.0, 0.0)
        ops.node(2, 1.0, 0.0)
        ops.node(3, 1.0, 1.0)
        ops.node(4, 0.0, 1.0)
        for n in (1, 2, 3, 4):
            ops.fix(n, 1, 1)
        ops.nDMaterial(*_hostun_args(1, c=25.0, sigmat=20.0, OCR=1.2, sig0=50.0))
        ops.element("quad", 1, 1, 2, 3, 4, 1.0, "PlaneStress", 1)
        ops.system("FullGeneral")
        ops.numberer("Plain")
        ops.constraints("Plain")
        ops.test("NormDispIncr", 1e-10, 20, 0)
        ops.algorithm("Newton")
        ops.integrator("LoadControl", 1.0)
        ops.analysis("Static")
        exx, eyy, gxy = float(eps_inplane[0]), float(eps_inplane[1]), float(eps_inplane[2])
        for n, (x, y) in ((1, (0.0, 0.0)), (2, (1.0, 0.0)), (3, (1.0, 1.0)), (4, (0.0, 1.0))):
            ops.setNodeDisp(n, 1, exx * x + 0.5 * gxy * y, "-commit")
            ops.setNodeDisp(n, 2, eyy * y + 0.5 * gxy * x, "-commit")
        assert ops.analyze(1) == 0
        sig3 = np.array(ops.eleResponse(1, "material", 1, "stress"))
        tan3 = np.array(ops.eleResponse(1, "material", 1, "tangent")).reshape(3, 3)
        return sig3, tan3

    eps_ps0 = np.array([-1.2e-3, -0.4e-3, 0.5e-3])
    _, tan_ps = eval_planestress(eps_ps0)
    delta_ps = 1.0e-7
    tan_ps_fd = np.zeros((3, 3))
    for j in range(3):
        ep = eps_ps0.copy()
        em = eps_ps0.copy()
        ep[j] += delta_ps
        em[j] -= delta_ps
        sp, _ = eval_planestress(ep)
        sm, _ = eval_planestress(em)
        tan_ps_fd[:, j] = (sp - sm) / (2.0 * delta_ps)
    rel_err_ps = np.linalg.norm(tan_ps - tan_ps_fd) / max(1.0, np.linalg.norm(tan_ps_fd))
    assert rel_err_ps < 1.0e-5

