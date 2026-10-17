"""Standalone C++17 verification suite for SRC/element/ladrunoMembrane/LadrunoMembraneKernel.h
cross-checking C++ kernel outputs against tests/_testbed/membrane_load_ref.py across
T1.1-T1.3 and T2.1-T3.5 to < 1e-11.
"""

from __future__ import annotations

import json
import math
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
_TESTBED = os.path.join(_HERE, "_testbed")
if _TESTBED not in sys.path:
    sys.path.insert(0, _TESTBED)

import membrane_load_ref as ref

_KERNEL_H = Path(_ROOT) / "SRC" / "element" / "ladrunoMembrane" / "LadrunoMembraneKernel.h"


_CPP_DRIVER_SOURCE = r"""
#include "LadrunoMembraneKernel.h"
#include <iomanip>
#include <iostream>
#include <sstream>

using namespace LadrunoMembraneKernel;

static void print_vec(const std::vector<double>& v) {
    std::cout << "[";
    for (std::size_t i = 0; i < v.size(); ++i) {
        if (i > 0) std::cout << ",";
        std::cout << std::setprecision(17) << v[i];
    }
    std::cout << "]";
}

template <std::size_t N>
static void print_arr(const std::array<double, N>& a) {
    std::cout << "[";
    for (std::size_t i = 0; i < N; ++i) {
        if (i > 0) std::cout << ",";
        std::cout << std::setprecision(17) << a[i];
    }
    std::cout << "]";
}

static void print_dist_result(const DistributionResult& r) {
    std::cout << "{";
    std::cout << "\"mode\":\"" << r.mode << "\",";
    std::cout << "\"total_applied_load\":" << std::setprecision(17) << r.total_applied_load << ",";
    std::cout << "\"total_applied_mx\":" << std::setprecision(17) << r.total_applied_mx << ",";
    std::cout << "\"total_applied_my\":" << std::setprecision(17) << r.total_applied_my << ",";
    std::cout << "\"plate_support_rank\":" << r.plate_support_rank << ",";
    std::cout << "\"support_node_reactions\":";
    print_vec(r.support_node_reactions);
    std::cout << ",\"frame_loads\":[";
    for (std::size_t k = 0; k < r.frame_loads.size(); ++k) {
        if (k > 0) std::cout << ",";
        const auto& fl = r.frame_loads[k];
        std::cout << "{";
        std::cout << "\"frame_id\":" << fl.frame_id << ",";
        std::cout << "\"aL\":" << std::setprecision(17) << fl.aL << ",";
        std::cout << "\"bL\":" << std::setprecision(17) << fl.bL << ",";
        std::cout << "\"wa_eq\":" << std::setprecision(17) << fl.wa_eq << ",";
        std::cout << "\"wb_eq\":" << std::setprecision(17) << fl.wb_eq << ",";
        std::cout << "\"F_total\":" << std::setprecision(17) << fl.F_total << ",";
        std::cout << "\"M_origin\":" << std::setprecision(17) << fl.M_origin << ",";
        std::cout << "\"overhang_torque\":" << std::setprecision(17) << fl.overhang_torque << ",";
        std::cout << "\"p0_2d\":";
        print_arr(fl.p0_2d);
        std::cout << ",\"q0_2d\":";
        print_arr(fl.q0_2d);
        std::cout << "}";
    }
    std::cout << "]}";
}

int main() {
    std::cout << "{\n";

    // 1. T1.1-T1.3 FEF 2D/3D + section sp(x)
    auto fef2d = partial_trapezoid_fef_2d(6.0, -12.0, -30.0, 4.0, 10.0, 0.2, 0.75, 1, 0.45);
    auto fef3d = partial_trapezoid_fef_3d(6.0, -12.0, -30.0, 8.0, -16.0, 4.0, 10.0, 0.2, 0.75, 1, 2, 0.2, 0.35);
    auto sp3d = force_beam_section_sp_3d(3.2, 8.0, -10.0, -35.0, 15.0, -25.0, 5.0, 12.0, 0.2, 0.8);

    std::cout << "\"fef2d_p0\":"; print_arr(fef2d.p0); std::cout << ",\n";
    std::cout << "\"fef2d_q0\":"; print_arr(fef2d.q0); std::cout << ",\n";
    std::cout << "\"fef2d_ql0\":"; print_arr(fef2d.ql0); std::cout << ",\n";
    std::cout << "\"fef3d_p0\":"; print_arr(fef3d.p0); std::cout << ",\n";
    std::cout << "\"fef3d_q0\":"; print_arr(fef3d.q0); std::cout << ",\n";
    std::cout << "\"fef3d_ql0\":"; print_arr(fef3d.ql0); std::cout << ",\n";
    std::cout << "\"sp3d\":[" << std::setprecision(17)
              << sp3d.P << "," << sp3d.Mz << "," << sp3d.Vy << "," << sp3d.My << "," << sp3d.Vz << "],\n";

    // 2. T2.1-T2.3: Acute 30-deg rhombus + partial beam + degeneracy check
    double deg30 = 30.0 * 3.14159265358979323846 / 180.0;
    double Ls = 5.0;
    std::vector<Point2D> skew_poly = {
        {0.0, 0.0},
        {Ls, 0.0},
        {Ls + Ls * std::cos(deg30), Ls * std::sin(deg30)},
        {Ls * std::cos(deg30), Ls * std::sin(deg30)}
    };
    auto m_skew = build_pslg_cdt_mesh(skew_poly);
    std::cout << "\"skew_area\":" << std::setprecision(17) << m_skew.total_area << ",\n";
    std::cout << "\"skew_nverts\":" << m_skew.vertices.size() << ",\n";
    std::cout << "\"skew_acute0\":" << (m_skew.is_acute[0] ? 1 : 0) << ",\n";
    std::cout << "\"skew_acute2\":" << (m_skew.is_acute[2] ? 1 : 0) << ",\n";

    bool caught_degen = false;
    try {
        build_pslg_cdt_mesh({{0.0, 0.0}, {2.0, 0.0}, {4.0, 0.0}});
    } catch (const std::invalid_argument&) {
        caught_degen = true;
    }
    std::cout << "\"caught_degen\":" << (caught_degen ? 1 : 0) << ",\n";

    // Partial edge + T-junction
    std::vector<Point2D> sq6 = {{0.0, 0.0}, {6.0, 0.0}, {6.0, 6.0}, {0.0, 6.0}};
    std::vector<Segment2D> frames_tj = {
        {{0.0, 0.0}, {6.0, 0.0}},
        {{0.0, 6.0}, {6.0, 6.0}},
        {{0.0, 3.0}, {6.0, 3.0}},
        {{3.0, 0.0}, {3.0, 3.0}}
    };
    auto res_tj = distribute_surface_load(sq6, 10.0, "tributary", frames_tj);
    std::cout << "\"res_tj\":"; print_dist_result(res_tj); std::cout << ",\n";

    std::vector<Segment2D> frames_part = {{{1.5, 0.0}, {4.5, 0.0}}};
    auto res_part = distribute_surface_load(sq6, 10.0, "tributary", {{{0.0, 0.0}, {6.0, 0.0}}});
    // Also test slab overhanging a beam vs beam partially along slab edge:
    std::vector<Segment2D> frames_sub = {{{0.0, -2.0}, {0.0, 6.0}}};
    std::vector<Point2D> slab_sub = {{0.0, 0.0}, {4.0, 0.0}, {4.0, 2.0}, {0.0, 2.0}};
    auto res_sub = distribute_surface_load(slab_sub, 10.0, "tributary", frames_sub);
    std::cout << "\"res_sub\":"; print_dist_result(res_sub); std::cout << ",\n";

    // 3. T3.1 Rectangular bay tributary (Lx=4, Ly=6)
    std::vector<Point2D> rect46 = {{0.0, 0.0}, {4.0, 0.0}, {4.0, 6.0}, {0.0, 6.0}};
    std::vector<Segment2D> frames46 = {
        {{0.0, 0.0}, {4.0, 0.0}},
        {{4.0, 0.0}, {4.0, 6.0}},
        {{4.0, 6.0}, {0.0, 6.0}},
        {{0.0, 6.0}, {0.0, 0.0}}
    };
    auto res_t31 = distribute_surface_load(rect46, 10.0, "tributary", frames46);
    std::cout << "\"res_t31\":"; print_dist_result(res_t31); std::cout << ",\n";

    // 4. T3.2 Oblique oneway (theta = 37 deg)
    double th = 37.0 * 3.14159265358979323846 / 180.0;
    double c = std::cos(th), s = std::sin(th);
    auto rot = [&](double u, double v) -> Point2D {
        return {u * c - v * s, u * s + v * c};
    };
    double L_u = 5.0, L_v = 3.5;
    std::vector<Point2D> poly_ob = {rot(0, 0), rot(L_u, 0), rot(L_u, L_v), rot(0, L_v)};
    std::vector<Segment2D> frames_ob = {
        {rot(0, 0), rot(0, L_v)},
        {rot(L_u, 0), rot(L_u, L_v)}
    };
    auto res_t32 = distribute_surface_load(poly_ob, 12.0, "oneway", frames_ob, {}, {}, {c, s});
    std::cout << "\"res_t32\":"; print_dist_result(res_t32); std::cout << ",\n";

    // 5. T3.3 fea_poisson isotropic & orthotropic
    std::vector<Point2D> sq4 = {{0.0, 0.0}, {4.0, 0.0}, {4.0, 4.0}, {0.0, 4.0}};
    std::vector<Segment2D> frames4 = {
        {{0.0, 0.0}, {4.0, 0.0}},
        {{4.0, 0.0}, {4.0, 4.0}},
        {{4.0, 4.0}, {0.0, 4.0}},
        {{0.0, 4.0}, {0.0, 0.0}}
    };
    auto res_t33_iso = distribute_surface_load(sq4, 8.0, "fea_poisson", frames4, {}, {}, {1.0, 0.0}, 1.0);
    auto res_t33_ortho = distribute_surface_load(sq4, 8.0, "fea_poisson", frames4, {}, {}, {1.0, 0.0}, 1e-4);
    std::cout << "\"res_t33_iso\":"; print_dist_result(res_t33_iso); std::cout << ",\n";
    std::cout << "\"res_t33_ortho\":"; print_dist_result(res_t33_ortho); std::cout << ",\n";

    // 6. T3.4 fea_plate DKT 3-edge cantilever & 1-edge single girder
    std::vector<Segment2D> frames_3edge = {
        {{0.0, 0.0}, {4.0, 0.0}},
        {{0.0, 4.0}, {4.0, 4.0}},
        {{0.0, 0.0}, {0.0, 4.0}}
    };
    std::vector<Segment2D> frames_1edge = {
        {{0.0, 0.0}, {0.0, 4.0}}
    };
    auto res_t34_3 = distribute_surface_load(sq4, 10.0, "fea_plate", frames_3edge);
    auto res_t34_1 = distribute_surface_load(sq4, 10.0, "fea_plate", frames_1edge);
    std::cout << "\"res_t34_3\":"; print_dist_result(res_t34_3); std::cout << ",\n";
    std::cout << "\"res_t34_1\":"; print_dist_result(res_t34_1); std::cout << ",\n";

    // 7. T3.5 Guyan static condensation of in-plane CST diaphragm
    std::vector<Point2D> poly_c = {{0.0, 0.0}, {4.0, 0.0}, {4.5, 3.0}, {0.5, 3.5}};
    auto mesh_c = build_pslg_cdt_mesh(poly_c, {}, {{{0.0, 0.0}, {4.0, 0.0}}}, {}, 1.2);
    auto cond = condense_inplane_cst(
        mesh_c.vertices,
        mesh_c.triangles,
        mesh_c.boundary_vertex_indices,
        30000.0, 0.2, 0.2, 2.4);
    std::cout << "\"cond_nb_dof\":" << cond.nb_dof << ",\n";
    std::cout << "\"cond_ni_dof\":" << cond.ni_dof << ",\n";
    std::cout << "\"cond_K\":"; print_vec(cond.K_cond); std::cout << ",\n";
    std::cout << "\"cond_M\":"; print_vec(cond.M_cond); std::cout << ",\n";
    std::cout << "\"cond_G\":"; print_vec(cond.G); std::cout << "\n";

    std::cout << "}\n";
    return 0;
}
"""


def _compile_and_run_cpp_driver(tmp_path: Path) -> dict:
    cpp_file = tmp_path / "driver.cpp"
    exe_file = tmp_path / ("driver.exe" if os.name == "nt" else "driver")
    cpp_file.write_text(_CPP_DRIVER_SOURCE, encoding="utf-8")

    env = os.environ.copy()
    include_dir = str(_KERNEL_H.parent)

    gpp = shutil.which("g++", path=env.get("PATH"))
    cl = shutil.which("cl", path=env.get("PATH"))
    if gpp:
        cmd = [gpp, "-O2", "-std=c++17", f"-I{include_dir}", str(cpp_file), "-o", str(exe_file)]
    elif cl:
        cmd = [cl, "/nologo", "/EHsc", "/O2", "/std:c++17", f"/I{include_dir}", str(cpp_file), f"/Fe:{exe_file}"]
    elif os.name == "nt":  # ci-coverage: portable fall back to MSVC detection on Windows when g++ is absent
        vswhere = Path(r"C:\Program Files (x86)\Microsoft Visual Studio\Installer\vswhere.exe")
        vcvars = None
        if vswhere.exists():
            res_vs = subprocess.run(
                [
                    str(vswhere),
                    "-latest",
                    "-products",
                    "*",
                    "-requires",
                    "Microsoft.VisualStudio.Component.VC.Tools.x86.x64",
                    "-property",
                    "installationPath",
                ],
                env=env,
                capture_output=True,
                text=True,
                check=False,
            )
            install_dir = res_vs.stdout.strip().splitlines()
            if install_dir:
                cand = Path(install_dir[0].strip()) / "VC" / "Auxiliary" / "Build" / "vcvars64.bat"
                if cand.exists():
                    vcvars = cand
        if vcvars is None:
            fallback = Path(
                r"C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat"
            )
            if fallback.exists():
                vcvars = fallback
        if vcvars is None:
            pytest.skip("Neither g++, cl.exe, nor vcvars64.bat found")
        bat_file = tmp_path / "build_driver.bat"
        bat_file.write_text(
            f'@echo off\ncall "{vcvars}" >nul\n'
            f'cl /nologo /EHsc /O2 /std:c++17 /I"{include_dir}" "{cpp_file}" /Fe:"{exe_file}"\n',
            encoding="utf-8",
        )
        cmd = ["cmd.exe", "/c", str(bat_file)]
    else:
        pytest.skip("Neither g++ nor cl.exe found on PATH")

    proc_compile = subprocess.run(
        cmd,
        cwd=str(tmp_path),
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc_compile.returncode == 0, (
        f"C++ compilation failed:\nSTDOUT:\n{proc_compile.stdout}\nSTDERR:\n{proc_compile.stderr}"
    )

    proc_run = subprocess.run(
        [str(exe_file)],
        cwd=str(tmp_path),
        env=env,
        capture_output=True,
        text=True,
        check=True,
    )
    return json.loads(proc_run.stdout)


def _assert_dist_matches(cpp_res: dict, py_res: ref.DistributionResult, atol: float = 1e-11) -> None:
    assert cpp_res["mode"] == py_res.mode
    np.testing.assert_allclose(cpp_res["total_applied_load"], py_res.total_applied_load, rtol=atol, atol=atol)
    np.testing.assert_allclose(cpp_res["total_applied_mx"], py_res.total_applied_mx, rtol=atol, atol=atol)
    np.testing.assert_allclose(cpp_res["total_applied_my"], py_res.total_applied_my, rtol=atol, atol=atol)
    assert cpp_res["plate_support_rank"] == py_res.plate_support_rank
    np.testing.assert_allclose(cpp_res["support_node_reactions"], py_res.support_node_reactions, rtol=atol, atol=atol)
    assert len(cpp_res["frame_loads"]) == len(py_res.frame_loads)
    for c_fl, p_fl in zip(cpp_res["frame_loads"], py_res.frame_loads):
        assert c_fl["frame_id"] == p_fl.frame_id
        np.testing.assert_allclose(c_fl["aL"], p_fl.aL, rtol=atol, atol=atol)
        np.testing.assert_allclose(c_fl["bL"], p_fl.bL, rtol=atol, atol=atol)
        np.testing.assert_allclose(c_fl["wa_eq"], p_fl.wa_eq, rtol=atol, atol=atol)
        np.testing.assert_allclose(c_fl["wb_eq"], p_fl.wb_eq, rtol=atol, atol=atol)
        np.testing.assert_allclose(c_fl["F_total"], p_fl.F_total, rtol=atol, atol=atol)
        np.testing.assert_allclose(c_fl["M_origin"], p_fl.M_origin, rtol=atol, atol=atol)
        np.testing.assert_allclose(c_fl["overhang_torque"], p_fl.overhang_torque, rtol=atol, atol=atol)
        np.testing.assert_allclose(c_fl["p0_2d"], p_fl.p0_2d, rtol=atol, atol=atol)
        np.testing.assert_allclose(c_fl["q0_2d"], p_fl.q0_2d, rtol=atol, atol=atol)


@pytest.mark.zone_a
@pytest.mark.t1
def test_cpp_kernel_matches_python_oracle_across_t1_to_t3(tmp_path: Path):
    """Compile LadrunoMembraneKernel.h with C++17 and verify < 1e-11 agreement
    against membrane_load_ref.py across T1.1-T1.3 and T2.1-T3.5.
    """
    data = _compile_and_run_cpp_driver(tmp_path)

    # 1. T1.1-T1.3
    py_fef2d = ref.partial_trapezoid_fef_2d(6.0, -12.0, -30.0, 4.0, 10.0, 0.2, 0.75, release=1, phi=0.45)
    py_fef3d = ref.partial_trapezoid_fef_3d(
        6.0, -12.0, -30.0, 8.0, -16.0, 4.0, 10.0, 0.2, 0.75, releasez=1, releasey=2, phiz=0.2, phiy=0.35
    )
    py_sp3d = ref.force_beam_section_sp_3d(3.2, 8.0, -10.0, -35.0, 15.0, -25.0, 5.0, 12.0, 0.2, 0.8)

    np.testing.assert_allclose(data["fef2d_p0"], py_fef2d["p0"], atol=1e-12)
    np.testing.assert_allclose(data["fef2d_q0"], py_fef2d["q0"], atol=1e-12)
    np.testing.assert_allclose(data["fef2d_ql0"], py_fef2d["ql0"], atol=1e-12)
    np.testing.assert_allclose(data["fef3d_p0"], py_fef3d["p0"], atol=1e-12)
    np.testing.assert_allclose(data["fef3d_q0"], py_fef3d["q0"], atol=1e-12)
    np.testing.assert_allclose(data["fef3d_ql0"], py_fef3d["ql0"], atol=1e-12)
    np.testing.assert_allclose(
        data["sp3d"],
        [py_sp3d["P"], py_sp3d["Mz"], py_sp3d["Vy"], py_sp3d["My"], py_sp3d["Vz"]],
        atol=1e-12,
    )

    # 2. T2.1-T2.3
    deg30 = math.radians(30.0)
    Ls = 5.0
    skew_poly = [
        (0.0, 0.0),
        (Ls, 0.0),
        (Ls + Ls * math.cos(deg30), Ls * math.sin(deg30)),
        (Ls * math.cos(deg30), Ls * math.sin(deg30)),
    ]
    py_skew = ref.build_pslg_cdt_mesh(skew_poly)
    assert abs(data["skew_area"] - py_skew.total_area) < 1e-12
    assert data["skew_nverts"] == len(py_skew.vertices)
    assert data["skew_acute0"] == 1
    assert data["skew_acute2"] == 1
    assert data["caught_degen"] == 1

    sq6 = [(0.0, 0.0), (6.0, 0.0), (6.0, 6.0), (0.0, 6.0)]
    frames_tj = [
        ((0.0, 0.0), (6.0, 0.0)),
        ((0.0, 6.0), (6.0, 6.0)),
        ((0.0, 3.0), (6.0, 3.0)),
        ((3.0, 0.0), (3.0, 3.0)),
    ]
    py_tj = ref.distribute_surface_load(sq6, 10.0, mode="tributary", frames=frames_tj)
    _assert_dist_matches(data["res_tj"], py_tj, atol=1e-11)

    slab_sub = [(0.0, 0.0), (4.0, 0.0), (4.0, 2.0), (0.0, 2.0)]
    frames_sub = [((0.0, -2.0), (0.0, 6.0))]
    py_sub = ref.distribute_surface_load(slab_sub, 10.0, mode="tributary", frames=frames_sub)
    _assert_dist_matches(data["res_sub"], py_sub, atol=1e-11)
    assert abs(data["res_sub"]["frame_loads"][0]["aL"] - 0.25) < 1e-12
    assert abs(data["res_sub"]["frame_loads"][0]["bL"] - 0.50) < 1e-12

    # 3. T3.1 Rectangular bay tributary
    rect46 = [(0.0, 0.0), (4.0, 0.0), (4.0, 6.0), (0.0, 6.0)]
    frames46 = [
        ((0.0, 0.0), (4.0, 0.0)),
        ((4.0, 0.0), (4.0, 6.0)),
        ((4.0, 6.0), (0.0, 6.0)),
        ((0.0, 6.0), (0.0, 0.0)),
    ]
    py_t31 = ref.distribute_surface_load(rect46, 10.0, mode="tributary", frames=frames46)
    _assert_dist_matches(data["res_t31"], py_t31, atol=1e-11)

    # 4. T3.2 Oblique oneway
    th = math.radians(37.0)
    c, s = math.cos(th), math.sin(th)
    rot = lambda u, v: (u * c - v * s, u * s + v * c)
    L_u, L_v = 5.0, 3.5
    poly_ob = [rot(0, 0), rot(L_u, 0), rot(L_u, L_v), rot(0, L_v)]
    frames_ob = [(rot(0, 0), rot(0, L_v)), (rot(L_u, 0), rot(L_u, L_v))]
    py_t32 = ref.distribute_surface_load(poly_ob, 12.0, mode="oneway", frames=frames_ob, span_dir=(c, s))
    _assert_dist_matches(data["res_t32"], py_t32, atol=1e-11)

    # 5. T3.3 fea_poisson
    sq4 = [(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)]
    frames4 = [
        ((0.0, 0.0), (4.0, 0.0)),
        ((4.0, 0.0), (4.0, 4.0)),
        ((4.0, 4.0), (0.0, 4.0)),
        ((0.0, 4.0), (0.0, 0.0)),
    ]
    py_t33_iso = ref.distribute_surface_load(sq4, 8.0, mode="fea_poisson", frames=frames4, orthotropy_eta=1.0)
    py_t33_ortho = ref.distribute_surface_load(
        sq4, 8.0, mode="fea_poisson", frames=frames4, span_dir=(1.0, 0.0), orthotropy_eta=1e-4
    )
    _assert_dist_matches(data["res_t33_iso"], py_t33_iso, atol=1e-11)
    _assert_dist_matches(data["res_t33_ortho"], py_t33_ortho, atol=1e-11)

    # 6. T3.4 fea_plate
    frames_3edge = [((0.0, 0.0), (4.0, 0.0)), ((0.0, 4.0), (4.0, 4.0)), ((0.0, 0.0), (0.0, 4.0))]
    frames_1edge = [((0.0, 0.0), (0.0, 4.0))]
    py_t34_3 = ref.distribute_surface_load(sq4, 10.0, mode="fea_plate", frames=frames_3edge)
    py_t34_1 = ref.distribute_surface_load(sq4, 10.0, mode="fea_plate", frames=frames_1edge)
    _assert_dist_matches(data["res_t34_3"], py_t34_3, atol=1e-11)
    _assert_dist_matches(data["res_t34_1"], py_t34_1, atol=1e-11)

    # 7. T3.5 Guyan static condensation of in-plane CST
    poly_c = [(0.0, 0.0), (4.0, 0.0), (4.5, 3.0), (0.5, 3.5)]
    mesh_c = ref.build_pslg_cdt_mesh(poly_c, frames=[((0.0, 0.0), (4.0, 0.0))], h_target=1.2)
    py_cond = ref.condense_inplane_cst(
        mesh_c.vertices,
        mesh_c.triangles,
        mesh_c.boundary_vertex_indices,
        E=30000.0,
        nu=0.2,
        t=0.2,
        rho=2.4,
    )
    nb_dof = data["cond_nb_dof"]
    ni_dof = data["cond_ni_dof"]
    assert nb_dof == py_cond["K_cond"].shape[0]
    assert ni_dof == py_cond["G"].shape[0]
    cpp_K = np.array(data["cond_K"], dtype=float).reshape((nb_dof, nb_dof))
    cpp_M = np.array(data["cond_M"], dtype=float).reshape((nb_dof, nb_dof))
    cpp_G = np.array(data["cond_G"], dtype=float).reshape((ni_dof, nb_dof))
    np.testing.assert_allclose(cpp_K, py_cond["K_cond"], atol=1e-11, rtol=1e-11)
    np.testing.assert_allclose(cpp_M, py_cond["M_cond"], atol=1e-11, rtol=1e-11)
    np.testing.assert_allclose(cpp_G, py_cond["G"], atol=1e-11, rtol=1e-11)
