"""ADR-87 D2 Mutation Gate for LadrunoMembrane (-DLADRUNO_MUTATE_MEMBRANE).

Verifies:
  1. Source-level contract: `#ifdef LADRUNO_MUTATE_MEMBRANE` is wired into
     `SRC/element/ladrunoMembrane/LadrunoMembrane.cpp` (both `applyCurrentLoads`
     and `getTangentStiff` / `getInitialStiff`) and
     `SRC/element/ladrunoMembrane/LadrunoMembraneKernel.h`.
  2. Compiled C++17 mutation kill score: compiles the C++ membrane kernel in
     both baseline (unmutated) and mutant (`-DLADRUNO_MUTATE_MEMBRANE`) modes
     and runs 6 physics probes (`tributary`, `oneway`, `fea_poisson` isotropic,
     `fea_poisson` orthotropic, `fea_plate` DKT, and Guyan-condensed CST
     diaphragm stiffness `K_cond`). Asserts 6/6 pass on the baseline binary and
     6/6 are killed (score = 1.00 >= 0.60 floor) on the mutant binary.
"""
from __future__ import annotations

import json
import math
import os
import shutil
import subprocess
from pathlib import Path

import pytest

_HERE = Path(__file__).resolve().parent
_ROOT = _HERE.parent
_ELE_DIR = _ROOT / "SRC" / "element" / "ladrunoMembrane"
_KERNEL_H = _ELE_DIR / "LadrunoMembraneKernel.h"
_ELE_CPP = _ELE_DIR / "LadrunoMembrane.cpp"

_PROBE_CPP = r"""
#include "LadrunoMembraneKernel.h"
#include <cmath>
#include <iomanip>
#include <iostream>

using namespace LadrunoMembraneKernel;

static double sum_frame_forces(const DistributionResult& r) {
    double s = 0.0;
    for (const auto& fl : r.frame_loads) {
        s += fl.F_total;
    }
    return s;
}

int main() {
    std::vector<Point2D> rect = {{0.0, 0.0}, {6.0, 0.0}, {6.0, 4.0}, {0.0, 4.0}};
    std::vector<Segment2D> frames = {
        {{0.0, 0.0}, {6.0, 0.0}},
        {{6.0, 0.0}, {6.0, 4.0}},
        {{6.0, 4.0}, {0.0, 4.0}},
        {{0.0, 4.0}, {0.0, 0.0}}
    };
    double q = 10.0;

    auto r_trib = distribute_surface_load(rect, q, "tributary", frames);
    auto r_oneway = distribute_surface_load(rect, q, "oneway", frames, {}, {}, {0.0, 1.0});
    auto r_pois_iso = distribute_surface_load(rect, q, "fea_poisson", frames);
    auto r_pois_ortho = distribute_surface_load(rect, q, "fea_poisson", frames, {}, {}, {0.0, 1.0}, 1e-4);
    auto r_plate = distribute_surface_load(rect, q, "fea_plate", frames);

    auto mesh = build_pslg_cdt_mesh(rect, {}, frames, {}, 1.5);
    auto cond = condense_inplane_cst(
        mesh.vertices, mesh.triangles, mesh.boundary_vertex_indices, 30000.0, 0.2, 0.2, 2.4);
    double k_trace = 0.0;
    for (int i = 0; i < cond.nb_dof; ++i) {
        k_trace += cond.K_cond[i * cond.nb_dof + i];
    }

    std::cout << std::setprecision(17) << "{"
              << "\"tributary_sum_f\":" << sum_frame_forces(r_trib) << ","
              << "\"oneway_sum_f\":" << sum_frame_forces(r_oneway) << ","
              << "\"poisson_iso_sum_f\":" << sum_frame_forces(r_pois_iso) << ","
              << "\"poisson_ortho_sum_f\":" << sum_frame_forces(r_pois_ortho) << ","
              << "\"plate_sum_f\":" << sum_frame_forces(r_plate) << ","
              << "\"cst_k_trace\":" << k_trace
              << "}\n";
    return 0;
}
"""


def _compile_and_run_probe(tmp_path: Path, mutate: bool) -> dict[str, float]:
    tag = "mut" if mutate else "base"
    cpp_file = tmp_path / f"probe_{tag}.cpp"
    exe_file = tmp_path / (f"probe_{tag}.exe" if os.name == "nt" else f"probe_{tag}")
    cpp_file.write_text(_PROBE_CPP, encoding="utf-8")

    env = os.environ.copy()
    include_dir = str(_ELE_DIR)

    gpp = shutil.which("g++", path=env.get("PATH"))
    cl = shutil.which("cl", path=env.get("PATH"))
    define_gpp = ["-DLADRUNO_MUTATE_MEMBRANE=1"] if mutate else []
    define_cl = ["/DLADRUNO_MUTATE_MEMBRANE=1"] if mutate else []

    if gpp:
        cmd = [gpp, "-O2", "-std=c++17", *define_gpp, f"-I{include_dir}", str(cpp_file), "-o", str(exe_file)]
        res = subprocess.run(cmd, cwd=str(tmp_path), env=env, capture_output=True, text=True, check=False)
    elif cl:
        cmd = [cl, "/nologo", "/EHsc", "/O2", "/std:c++17", *define_cl, f"/I{include_dir}", str(cpp_file), f"/Fe:{exe_file}"]
        res = subprocess.run(cmd, cwd=str(tmp_path), env=env, capture_output=True, text=True, check=False)
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
            pytest.skip("No C++17 compiler (g++ or cl) found in PATH.")
        def_str = " /DLADRUNO_MUTATE_MEMBRANE=1" if mutate else ""
        bat_file = tmp_path / f"build_{tag}.bat"
        bat_file.write_text(
            f'@echo off\ncall "{vcvars}" >nul\n'
            f'cl /nologo /EHsc /O2 /std:c++17{def_str} /I"{include_dir}" "{cpp_file}" /Fe:"{exe_file}"\n',
            encoding="utf-8",
        )
        res = subprocess.run(["cmd.exe", "/c", str(bat_file)], cwd=str(tmp_path), env=env, capture_output=True, text=True, check=False)
    else:
        pytest.skip("No C++17 compiler found in PATH.")

    assert res.returncode == 0, f"Compile failed ({tag}):\nSTDOUT:\n{res.stdout}\nSTDERR:\n{res.stderr}"
    run_res = subprocess.run([str(exe_file)], cwd=str(tmp_path), env=env, capture_output=True, text=True, check=True)
    return json.loads(run_res.stdout)


def test_ladruno_mutate_membrane_source_wiring():
    """Ensure LADRUNO_MUTATE_MEMBRANE is wired in both LadrunoMembrane.cpp and LadrunoMembraneKernel.h."""
    cpp_text = _ELE_CPP.read_text(encoding="utf-8")
    kernel_text = _KERNEL_H.read_text(encoding="utf-8")
    assert cpp_text.count("LADRUNO_MUTATE_MEMBRANE") >= 2
    assert kernel_text.count("LADRUNO_MUTATE_MEMBRANE") >= 2


def test_ladruno_mutate_membrane_kills_all_physics_probes(tmp_path: Path):
    """Compile baseline vs -DLADRUNO_MUTATE_MEMBRANE and assert 6/6 physics probes are killed."""
    base = _compile_and_run_probe(tmp_path, mutate=False)
    mut = _compile_and_run_probe(tmp_path, mutate=True)

    expected_load = 10.0 * 6.0 * 4.0  # 240.0

    def _evaluate_probes(data: dict[str, float]) -> dict[str, bool]:
        return {
            "tributary": math.isclose(data["tributary_sum_f"], expected_load, rel_tol=1e-6),
            "oneway": math.isclose(data["oneway_sum_f"], expected_load, rel_tol=1e-6),
            "poisson_iso": math.isclose(data["poisson_iso_sum_f"], expected_load, rel_tol=1e-6),
            "poisson_ortho": math.isclose(data["poisson_ortho_sum_f"], expected_load, rel_tol=1e-6),
            "plate": math.isclose(data["plate_sum_f"], expected_load, rel_tol=1e-6),
            "cst_k_trace": data["cst_k_trace"] > 1.0e3,
        }

    base_results = _evaluate_probes(base)
    mut_results = _evaluate_probes(mut)

    # All 6 probes must pass on the unmutated baseline
    assert all(base_results.values()), f"Baseline probes failed: {base_results}"

    # All 6 probes must FAIL (be killed) on the -DLADRUNO_MUTATE_MEMBRANE binary
    killed = [name for name, passed in mut_results.items() if not passed]
    kill_score = len(killed) / len(mut_results)
    assert kill_score == 1.0, f"Expected 100% kill score under LADRUNO_MUTATE_MEMBRANE, got {kill_score} (mut={mut})"
