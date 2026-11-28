"""WP-260b BEFORE vs AFTER bit-identity (G-BYTE) MD5 + signbit harness.

Runs a battery of explicit-wave, T7-skip, -lnvd Transformation+equalDOF,
formUri(0) damage-path, and lumped Rayleigh-alphaM decks on two opensees.pyd
directories (`--before-bin` and `--after-bin`) inside isolated `python -S`
subprocesses, verifying IEEE-754 byte MD5, signbit identity, and max rel diff.

Refuses trivial self-comparisons (`before_sha == after_sha`) and supports
optional `--expect-before <sha>` / `--expect-after <sha>` prefix assertions
(minimum 9 characters).

Usage:
    py -3.12 Ladruno_scripts/wp260b_ab_md5.py \\
        --before-bin <base_dist_bin> \\
        --after-bin  <head_dist_bin> \\
        [--expect-before <sha9+>] [--expect-after <sha9+>]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import struct
import subprocess
import sys


_PROBE_SHA = r'''
import os, sys

BIN_DIR = sys.argv[1]
if hasattr(os, "add_dll_directory") and os.path.isdir(BIN_DIR):
    os.add_dll_directory(BIN_DIR)
sys.path.insert(0, BIN_DIR)
os.environ["PATH"] = BIN_DIR + os.pathsep + os.environ.get("PATH", "")
os.environ["LADRUNO_OPENSEES_QUIET"] = "1"

import opensees as ops
print("SHA:" + str(ops.ladrunoBuild()).strip())
'''


_CHILD = r'''
import hashlib, json, math, os, struct, sys

BIN_DIR = sys.argv[1]
CASE = sys.argv[2]

if hasattr(os, "add_dll_directory") and os.path.isdir(BIN_DIR):
    os.add_dll_directory(BIN_DIR)
sys.path.insert(0, BIN_DIR)
os.environ["PATH"] = BIN_DIR + os.pathsep + os.environ.get("PATH", "")
os.environ["LADRUNO_OPENSEES_QUIET"] = "1"

import opensees as ops
ops.ladrunoThreads(1)
build_sha = str(ops.ladrunoBuild()).strip()


def _hex8_nodes(nx, ny, nz, L=1.0):
    def nid(i, j, k):
        return 1 + i + (nx + 1) * (j + (ny + 1) * k)
    for k in range(nz + 1):
        for j in range(ny + 1):
            for i in range(nx + 1):
                ops.node(nid(i, j, k), i * L, j * L, k * L)
    return nid


def _run_d3_cpml(integrator_kind):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    E, nu, rho = 1.0e8, 0.25, 2000.0
    ops.nDMaterial("ElasticIsotropic", 1, E, nu, rho)
    nx, ny, nz = 6, 6, 6
    nid = _hex8_nodes(nx, ny, nz, 1.0)
    for k in range(nz + 1):
        for j in range(ny + 1):
            for i in range(nx + 1):
                if i == 0 or i == nx or j == 0 or j == ny or k == 0:
                    ops.fix(nid(i, j, k), 1, 1, 1)
    etag = 1
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                ns = (
                    nid(i, j, k), nid(i + 1, j, k), nid(i + 1, j + 1, k), nid(i, j + 1, k),
                    nid(i, j, k + 1), nid(i + 1, j, k + 1), nid(i + 1, j + 1, k + 1), nid(i, j + 1, k + 1),
                )
                in_core = (1 <= i <= nx - 2) and (1 <= j <= ny - 2) and (k >= 1)
                if in_core:
                    ops.element("LadrunoBrick", etag, *ns, 1, "-formulation", "uri", "-lumped")
                else:
                    ops.element(
                        "LadrunoCPML", etag, *ns,
                        "-E", E, "-nu", nu, "-rho", rho,
                        "-box", 1.0, float(nx - 1), 1.0, float(ny - 1), 1.0, float(nz),
                        1.0, 1.0, 1.0,
                        "-m", 2.0, "-R0", 1e-4, "-kappaMax", 1.0, "-alphaMax", 0.0,
                        "-lumped",
                    )
                etag += 1
    ops.timeSeries("Trig", 1, 0.0, 0.02, 0.01)
    ops.pattern("Plain", 1, 1)
    ops.load(nid(nx // 2, ny // 2, nz), 1.0e5, -2.0e4, 3.0e4)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    if integrator_kind == "cdl":
        ops.integrator("CentralDifferenceLadruno")
    else:
        ops.integrator("ExplicitBathe", 0.54)
    ops.analysis("Transient")
    for _ in range(40):
        assert ops.analyze(1, 4.0e-4) == 0
    vals = []
    for k in range(nz + 1):
        for j in range(ny + 1):
            for i in range(nx + 1):
                vals.extend(ops.nodeDisp(nid(i, j, k)))
    return vals


def _run_eb_lnvd_transf_equaldof():
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.uniaxialMaterial("Elastic", 1, 100.0)
    ops.uniaxialMaterial("Elastic", 2, 250.0)
    ops.node(1, 0.0, 0.0)
    ops.node(2, 1.0, 0.0)
    ops.node(3, 1.0, 1.0)
    ops.node(4, 0.0, 1.0)
    ops.node(5, 2.0, 0.0)
    ops.fix(1, 1, 1)
    ops.fix(4, 1, 1)
    ops.equalDOF(2, 3, 1, 2)
    ops.mass(2, 1.5, 1.5)
    ops.mass(3, 2.0, 2.0)
    ops.mass(5, 1.2, 1.2)
    ops.element("Truss", 1, 1, 2, 1.0, 1)
    ops.element("Truss", 2, 4, 3, 1.0, 2)
    ops.element("Truss", 3, 2, 5, 1.0, 1)
    ops.element("Truss", 4, 3, 5, 1.0, 2)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 15.0, -5.0)
    ops.load(5, 10.0, 8.0)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    ops.integrator("ExplicitBathe", 0.54, "-lnvd", 0.8)
    ops.analysis("Transient")
    vals = []
    for _ in range(120):
        assert ops.analyze(1, 0.02) == 0
        for nd in (2, 3, 5):
            vals.extend(ops.nodeDisp(nd))
    return vals


def _run_t7_brick20():
    vals = []
    corners = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
               (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
    mids = [(0.5, 0, 0), (1, 0.5, 0), (0.5, 1, 0), (0, 0.5, 0),
            (0.5, 0, 1), (1, 0.5, 1), (0.5, 1, 1), (0, 0.5, 1),
            (0, 0, 0.5), (1, 0, 0.5), (1, 1, 0.5), (0, 1, 0.5)]
    for form in ("std", "uri"):
        for lumped in (True, False):
            ops.wipe()
            ops.model("basic", "-ndm", 3, "-ndf", 3)
            ops.nDMaterial("ElasticIsotropic", 1, 1.0e7, 0.25, 2.0)
            for idx, (x, y, z) in enumerate(corners + mids, start=1):
                ops.node(idx, float(x), float(y), float(z))
                if z == 0.0:
                    ops.fix(idx, 1, 1, 1)
            extra = ["-lumped"] if lumped else []
            ops.element("LadrunoBrick20", 1, *range(1, 21), 1, "-formulation", form, *extra)
            ops.timeSeries("Constant", 1)
            ops.pattern("Plain", 1, 1)
            ops.load(7, 100.0, -50.0, 25.0)
            ops.constraints("Plain")
            ops.numberer("Plain")
            ops.system("Diagonal" if lumped else "BandGeneral")
            ops.test("NormDispIncr", 1e-12, 1)
            ops.algorithm("Linear")
            ops.integrator("CentralDifferenceLadruno")
            ops.analysis("Transient")
            for _ in range(30):
                assert ops.analyze(1, 5.0e-5) == 0
            for idx in range(5, 21):
                vals.extend(ops.nodeDisp(idx))
    return vals


def _run_t7_quad_all():
    vals = []
    # 2x2 LadrunoQuad mesh (4 elements, 6 free nodes = 12 DOFs) across all 6
    # (formulation, geom) combinations, recorded every 10 steps over 40 steps
    # -> 6 * 4 * 12 = 288 values.
    def qnid(i, j):
        return 1 + i + 3 * j

    for form, geom in (("std", "linear"), ("bbar", "linear"), ("ssp", "linear"),
                       ("eas", "linear"), ("std", "finite"), ("bbar", "finite")):
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 2)
        ops.nDMaterial("ElasticIsotropic", 1, 1.0e7, 0.25, 2.0)
        mat_tag = 1
        if geom == "finite":
            ops.nDMaterial("LogStrain2D", 2, 1)
            mat_tag = 2
        for j in range(3):
            for i in range(3):
                nd = qnid(i, j)
                ops.node(nd, float(i), float(j))
                if j == 0:
                    ops.fix(nd, 1, 1)
        etag = 1
        for j in range(2):
            for i in range(2):
                ops.element(
                    "LadrunoQuad", etag,
                    qnid(i, j), qnid(i + 1, j), qnid(i + 1, j + 1), qnid(i, j + 1),
                    mat_tag, "-type", "PlaneStrain",
                    "-thick", 1.0, "-formulation", form, "-geom", geom,
                )
                etag += 1
        ops.timeSeries("Constant", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(qnid(0, 2), 300.0, 100.0)
        ops.load(qnid(1, 2), 400.0, -150.0)
        ops.load(qnid(2, 2), 500.0, -250.0)
        ops.constraints("Plain")
        ops.numberer("Plain")
        ops.system("Diagonal")
        ops.test("NormDispIncr", 1e-12, 1)
        ops.algorithm("Linear")
        ops.integrator("ExplicitBathe", 0.54)
        ops.analysis("Transient")
        for step in range(1, 41):
            assert ops.analyze(1, 5.0e-5) == 0
            if step % 10 == 0:
                for j in (1, 2):
                    for i in range(3):
                        vals.extend(ops.nodeDisp(qnid(i, j)))
    assert len(vals) >= 200, f"t7_quad_all emitted only {len(vals)} values"
    return vals


def _run_t7_shell_callers():
    vals = []

    def g3nid(i, j):
        return 1 + i + 3 * j

    # 1) LadrunoCST: 2x2 grid split into 8 CST triangles, 6 free nodes (12 DOFs),
    # recorded every 5 steps over 40 steps -> 8 * 12 = 96 values.
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.nDMaterial("ElasticIsotropic", 1, 1.0e7, 0.25, 2.0)
    for j in range(3):
        for i in range(3):
            nd = g3nid(i, j)
            ops.node(nd, float(i), float(j))
            if j == 0:
                ops.fix(nd, 1, 1)
    etag = 1
    for j in range(2):
        for i in range(2):
            n00, n10 = g3nid(i, j), g3nid(i + 1, j)
            n01, n11 = g3nid(i, j + 1), g3nid(i + 1, j + 1)
            ops.element("LadrunoCST", etag, n00, n10, n11, 1, "-type", "PlaneStrain", "-thick", 1.0)
            etag += 1
            ops.element("LadrunoCST", etag, n00, n11, n01, 1, "-type", "PlaneStrain", "-thick", 1.0)
            etag += 1
    ops.timeSeries("Constant", 1); ops.pattern("Plain", 1, 1)
    ops.load(g3nid(0, 2), 150.0, 50.0)
    ops.load(g3nid(1, 2), 200.0, -100.0)
    ops.load(g3nid(2, 2), 250.0, -120.0)
    ops.constraints("Plain"); ops.numberer("Plain"); ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1); ops.algorithm("Linear")
    ops.integrator("CentralDifferenceLadruno"); ops.analysis("Transient")
    for step in range(1, 41):
        assert ops.analyze(1, 5e-5) == 0
        if step % 5 == 0:
            for j in (1, 2):
                for i in range(3):
                    vals.extend(ops.nodeDisp(g3nid(i, j)))

    # 2) LadrunoLST: 2 six-node triangles on [0,1]x[0,1], 6 free nodes (12 DOFs),
    # recorded every 5 steps over 35 steps -> 7 * 12 = 84 values.
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.nDMaterial("ElasticIsotropic", 1, 1.0e7, 0.25, 2.0)
    for j in range(3):
        for i in range(3):
            nd = g3nid(i, j)
            ops.node(nd, 0.5 * i, 0.5 * j)
            if j == 0:
                ops.fix(nd, 1, 1)
    ops.element(
        "LadrunoLST", 1,
        g3nid(0, 0), g3nid(2, 0), g3nid(2, 2),
        g3nid(1, 0), g3nid(2, 1), g3nid(1, 1),
        1, "-type", "PlaneStrain", "-thick", 1.0,
    )
    ops.element(
        "LadrunoLST", 2,
        g3nid(0, 0), g3nid(2, 2), g3nid(0, 2),
        g3nid(1, 1), g3nid(1, 2), g3nid(0, 1),
        1, "-type", "PlaneStrain", "-thick", 1.0,
    )
    ops.timeSeries("Constant", 1); ops.pattern("Plain", 1, 1)
    ops.load(g3nid(0, 2), 200.0, -100.0)
    ops.load(g3nid(2, 2), 150.0, -80.0)
    ops.constraints("Plain"); ops.numberer("Plain"); ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1); ops.algorithm("Linear")
    ops.integrator("CentralDifferenceLadruno"); ops.analysis("Transient")
    for step in range(1, 36):
        assert ops.analyze(1, 5e-5) == 0
        if step % 5 == 0:
            for j in (1, 2):
                for i in range(3):
                    vals.extend(ops.nodeDisp(g3nid(i, j)))

    # 3) BezierTri6 (consistent mass via BandGeneral + CDL): 2 six-node triangles
    # on [0,1]x[0,1], 6 free nodes (12 DOFs), recorded every 5 steps over 35 steps -> 84 values.
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.nDMaterial("ElasticIsotropic", 1, 1.0e7, 0.25, 2.0)
    for j in range(3):
        for i in range(3):
            nd = g3nid(i, j)
            ops.node(nd, 0.5 * i, 0.5 * j)
            if j == 0:
                ops.fix(nd, 1, 1)
    ops.element(
        "BezierTri6", 1,
        g3nid(0, 0), g3nid(2, 0), g3nid(2, 2),
        g3nid(1, 0), g3nid(2, 1), g3nid(1, 1),
        1.0, "PlaneStrain", 1,
    )
    ops.element(
        "BezierTri6", 2,
        g3nid(0, 0), g3nid(2, 2), g3nid(0, 2),
        g3nid(1, 1), g3nid(1, 2), g3nid(0, 1),
        1.0, "PlaneStrain", 1,
    )
    ops.timeSeries("Constant", 1); ops.pattern("Plain", 1, 1)
    ops.load(g3nid(0, 2), 200.0, -100.0)
    ops.load(g3nid(2, 2), 150.0, -80.0)
    ops.constraints("Plain"); ops.numberer("Plain"); ops.system("BandGeneral")
    ops.test("NormDispIncr", 1e-12, 1); ops.algorithm("Linear")
    ops.integrator("CentralDifferenceLadruno"); ops.analysis("Transient")
    for step in range(1, 36):
        assert ops.analyze(1, 5e-5) == 0
        if step % 5 == 0:
            for j in (1, 2):
                for i in range(3):
                    vals.extend(ops.nodeDisp(g3nid(i, j)))

    assert len(vals) >= 200, f"t7_shell_callers emitted only {len(vals)} values"
    return vals


def _run_brick_formuri_damage():
    vals = []
    # Concrete3D active-cracking damage path on a 2x1x1 LadrunoBrick -formulation uri
    # patch (triggers damageResponse != 0 and damageScale() < 1 in LadrunoBrick::formUri(0),
    # verifying skipped getTangent() + stiffness hourglass control).
    # 6 free top nodes (18 DOFs) recorded every 2 steps across 40 steps -> 20 * 18 = 360 values.
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("LadrunoConcrete3D", 1, 3.0e10, 0.2, 3.0e7, 3.0e6, 100.0, 25000.0, "-rho", 2400.0)
    nid = _hex8_nodes(2, 1, 1, 1.0)
    for j in range(2):
        for i in range(3):
            ops.fix(nid(i, j, 0), 1, 1, 1)
    for i in range(2):
        ns = (
            nid(i, 0, 0), nid(i + 1, 0, 0), nid(i + 1, 1, 0), nid(i, 1, 0),
            nid(i, 0, 1), nid(i + 1, 0, 1), nid(i + 1, 1, 1), nid(i, 1, 1),
        )
        ops.element("LadrunoBrick", i + 1, *ns, 1,
                    "-formulation", "uri", "-lumped", "-hourglass", "stiffness", 0.05)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    # Strong tension + asymmetric shear so Concrete3D actively cracks across both elements
    for j in range(2):
        for i in range(3):
            ops.load(nid(i, j, 1), 1.5e6, 0.3e6 * (i - 1), 2.0e6)
    ops.constraints("Plain"); ops.numberer("Plain"); ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1); ops.algorithm("Linear")
    ops.integrator("CentralDifferenceLadruno"); ops.analysis("Transient")
    for step in range(1, 41):
        assert ops.analyze(1, 5.0e-5) == 0
        if step % 2 == 0:
            for j in range(2):
                for i in range(3):
                    vals.extend(ops.nodeDisp(nid(i, j, 1)))
    for etag in (1, 2):
        dmg = ops.eleResponse(etag, "material", "1", "damage")
        assert dmg and max(dmg) > 0.5, f"ele {etag}: expected active Concrete3D damage > 0.5, got {dmg}"
    assert len(vals) >= 200, f"brick_formuri_damage emitted only {len(vals)} values"
    return vals


def _run_quad_bv_cd_cache():
    vals = []
    # Allowlisted ElasticIsotropic on LadrunoQuad with -bulkViscosity (Finding 7b)
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.nDMaterial("ElasticIsotropic", 1, 1.0e7, 0.25, 2.0)
    for i, (x, y) in enumerate([(0, 0), (1, 0), (1, 1), (0, 1)], 1):
        ops.node(i, float(x), float(y))
        if y == 0:
            ops.fix(i, 1, 1)
    ops.element("LadrunoQuad", 1, 1, 2, 3, 4, 1, "-type", "PlaneStrain", "-thick", 1.0,
                "-formulation", "std", "-bulkViscosity", 0.5, 1.2)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(3, 200.0, -500.0)
    ops.load(4, -100.0, -500.0)
    ops.constraints("Plain"); ops.numberer("Plain"); ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1); ops.algorithm("Linear")
    ops.integrator("ExplicitBathe", 0.54); ops.analysis("Transient")
    for _ in range(50):
        assert ops.analyze(1, 5.0e-5) == 0
        vals.extend(ops.nodeDisp(3))
        vals.extend(ops.nodeDisp(4))
    return vals


def _run_rayleigh_alphaM_all():
    vals = []
    # 1) LadrunoBrick (uri, lumped)
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("ElasticIsotropic", 1, 1.0e7, 0.25, 2.0)
    coords = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
              (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
    for i, (x, y, z) in enumerate(coords, 1):
        ops.node(i, float(x), float(y), float(z))
        if z == 0: ops.fix(i, 1, 1, 1)
    ops.element("LadrunoBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1, "-formulation", "uri", "-lumped")
    ops.rayleigh(25.0, 0.0, 0.0, 0.0)
    ops.timeSeries("Constant", 1); ops.pattern("Plain", 1, 1); ops.load(7, 1000.0, -300.0, 200.0)
    ops.constraints("Plain"); ops.numberer("Plain"); ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1); ops.algorithm("Linear")
    ops.integrator("CentralDifferenceLadruno"); ops.analysis("Transient")
    for _ in range(40): assert ops.analyze(1, 1e-4) == 0
    for n in (5, 6, 7, 8): vals.extend(ops.nodeDisp(n))

    # 2) LadrunoBrick20 (lumped)
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("ElasticIsotropic", 1, 1.0e7, 0.25, 2.0)
    mids = [(0.5, 0, 0), (1, 0.5, 0), (0.5, 1, 0), (0, 0.5, 0),
            (0.5, 0, 1), (1, 0.5, 1), (0.5, 1, 1), (0, 0.5, 1),
            (0, 0, 0.5), (1, 0, 0.5), (1, 1, 0.5), (0, 1, 0.5)]
    for idx, (x, y, z) in enumerate(coords + mids, start=1):
        ops.node(idx, float(x), float(y), float(z))
        if z == 0.0: ops.fix(idx, 1, 1, 1)
    ops.element("LadrunoBrick20", 1, *range(1, 21), 1, "-formulation", "std", "-lumped")
    ops.rayleigh(25.0, 0.0, 0.0, 0.0)
    ops.timeSeries("Constant", 1); ops.pattern("Plain", 1, 1); ops.load(7, 100.0, -50.0, 20.0)
    ops.constraints("Plain"); ops.numberer("Plain"); ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1); ops.algorithm("Linear")
    ops.integrator("CentralDifferenceLadruno"); ops.analysis("Transient")
    for _ in range(40): assert ops.analyze(1, 5e-5) == 0
    for idx in range(5, 21): vals.extend(ops.nodeDisp(idx))

    # 3) LadrunoQuad
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.nDMaterial("ElasticIsotropic", 1, 1.0e7, 0.25, 2.0)
    for i, (x, y) in enumerate([(0, 0), (1, 0), (1, 1), (0, 1)], 1):
        ops.node(i, float(x), float(y))
        if y == 0: ops.fix(i, 1, 1)
    ops.element("LadrunoQuad", 1, 1, 2, 3, 4, 1, "-type", "PlaneStrain", "-thick", 1.0)
    ops.rayleigh(25.0, 0.0, 0.0, 0.0)
    ops.timeSeries("Constant", 1); ops.pattern("Plain", 1, 1); ops.load(3, 100.0, -50.0)
    ops.constraints("Plain"); ops.numberer("Plain"); ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1); ops.algorithm("Linear")
    ops.integrator("CentralDifferenceLadruno"); ops.analysis("Transient")
    for _ in range(40): assert ops.analyze(1, 1e-4) == 0
    vals.extend(ops.nodeDisp(3)); vals.extend(ops.nodeDisp(4))

    # 4) LadrunoCPML (lumped)
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for i, (x, y, z) in enumerate(coords, 1):
        ops.node(i, float(x), float(y), float(z))
        if z == 0: ops.fix(i, 1, 1, 1)
    ops.element(
        "LadrunoCPML", 1, 1, 2, 3, 4, 5, 6, 7, 8,
        "-E", 1.0e8, "-nu", 0.25, "-rho", 2000.0,
        "-box", 0.0, 0.5, 0.0, 0.5, 0.0, 0.5, 1.0, 1.0, 1.0,
        "-m", 2.0, "-R0", 1e-4, "-kappaMax", 1.0, "-alphaMax", 0.0,
        "-lumped",
    )
    ops.rayleigh(25.0, 0.0, 0.0, 0.0)
    ops.timeSeries("Constant", 1); ops.pattern("Plain", 1, 1); ops.load(7, 1000.0, -400.0, 200.0)
    ops.constraints("Plain"); ops.numberer("Plain"); ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1); ops.algorithm("Linear")
    ops.integrator("CentralDifferenceLadruno"); ops.analysis("Transient")
    for _ in range(40): assert ops.analyze(1, 5e-4) == 0
    for n in (5, 6, 7, 8): vals.extend(ops.nodeDisp(n))
    return vals


CASES = {
    "d3a_cdl_cpml": lambda: _run_d3_cpml("cdl"),
    "d3b_eb_cpml": lambda: _run_d3_cpml("eb"),
    "eb_lnvd_transf_equaldof": _run_eb_lnvd_transf_equaldof,
    "t7_brick20": _run_t7_brick20,
    "t7_quad_all": _run_t7_quad_all,
    "t7_shell_callers": _run_t7_shell_callers,
    "brick_formuri_damage": _run_brick_formuri_damage,
    "quad_bv_cd_cache": _run_quad_bv_cd_cache,
    "rayleigh_alphaM_all": _run_rayleigh_alphaM_all,
}

vals = CASES[CASE]()
raw = struct.pack(f"<{len(vals)}d", *vals)
md5 = hashlib.md5(raw).hexdigest()
print("JSON:" + json.dumps({"sha": build_sha, "md5": md5, "vals": vals}))
'''


CASE_NAMES = [
    "d3a_cdl_cpml",
    "d3b_eb_cpml",
    "eb_lnvd_transf_equaldof",
    "t7_brick20",
    "t7_quad_all",
    "t7_shell_callers",
    "brick_formuri_damage",
    "quad_bv_cd_cache",
    "rayleigh_alphaM_all",
]


def _probe_sha(bin_dir: str) -> str:
    proc = subprocess.run(
        [sys.executable, "-S", "-c", _PROBE_SHA, os.path.abspath(bin_dir)],
        stdin=subprocess.DEVNULL,
        capture_output=True,
        text=True,
        timeout=60,
        encoding="utf-8",
        errors="replace",
    )
    if proc.returncode != 0:
        raise RuntimeError(
            f"failed to probe ladrunoBuild() in {bin_dir} (rc={proc.returncode})\n"
            f"STDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}"
        )
    for line in proc.stdout.splitlines():
        if line.startswith("SHA:"):
            return line[4:].strip()
    raise RuntimeError(f"missing SHA output from {bin_dir}:\n{proc.stdout}\n{proc.stderr}")


def _run_one(bin_dir: str, case: str) -> dict:
    proc = subprocess.run(
        [sys.executable, "-S", "-c", _CHILD, os.path.abspath(bin_dir), case],
        stdin=subprocess.DEVNULL,
        capture_output=True,
        text=True,
        timeout=180,
        encoding="utf-8",
        errors="replace",
    )
    if proc.returncode != 0:
        raise RuntimeError(
            f"case {case} failed on {bin_dir} (rc={proc.returncode})\n"
            f"STDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}"
        )
    for line in proc.stdout.splitlines():
        if line.startswith("JSON:"):
            return json.loads(line[5:])
    raise RuntimeError(f"missing JSON output for {case} on {bin_dir}:\n{proc.stdout}\n{proc.stderr}")


def _check_expect_prefix(flag: str, expected: str, actual: str, label: str) -> str | None:
    exp = expected.strip()
    if not exp:
        return None
    if len(exp) < 9:
        return f"ERROR: {flag} requires at least 9 characters (got {len(exp)}: {expected!r})"
    if not actual.lower().startswith(exp.lower()):
        return f"ERROR: {label} SHA {actual!r} does not match {flag} prefix {exp!r}"
    return None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--before-bin", "--before", dest="before_bin", required=True,
                    help="dist/bin directory of BEFORE (base) build")
    ap.add_argument("--after-bin", "--after", dest="after_bin", required=True,
                    help="dist/bin directory of AFTER (head) build")
    ap.add_argument("--expect-before", default="",
                    help="expected BEFORE ladrunoBuild SHA prefix (at least 9 characters)")
    ap.add_argument("--expect-after", default="",
                    help="expected AFTER ladrunoBuild SHA prefix (at least 9 characters)")
    ap.add_argument("--out", default="", help="optional path to write markdown table")
    args = ap.parse_args()

    for flag, val in (("--expect-before", args.expect_before), ("--expect-after", args.expect_after)):
        s = val.strip()
        if s and len(s) < 9:
            print(f"ERROR: {flag} requires at least 9 characters (got {len(s)}: {val!r})", file=sys.stderr)
            return 2

    before_sha = _probe_sha(args.before_bin)
    after_sha = _probe_sha(args.after_bin)

    if not before_sha or before_sha.lower() == "unknown":
        print(f"ERROR: invalid BEFORE ladrunoBuild SHA: {before_sha!r}", file=sys.stderr)
        return 2
    if not after_sha or after_sha.lower() == "unknown":
        print(f"ERROR: invalid AFTER ladrunoBuild SHA: {after_sha!r}", file=sys.stderr)
        return 2

    for flag, exp, act, lbl in (
        ("--expect-before", args.expect_before, before_sha, "BEFORE"),
        ("--expect-after", args.expect_after, after_sha, "AFTER"),
    ):
        err = _check_expect_prefix(flag, exp, act, lbl)
        if err is not None:
            print(err, file=sys.stderr)
            return 2

    if before_sha == after_sha:
        print(
            f"ERROR: BEFORE and AFTER builds have the same ladrunoBuild SHA ({before_sha}); "
            f"refusing trivial self-comparison.",
            file=sys.stderr,
        )
        return 2

    rows = []
    any_fail = False
    for case in CASE_NAMES:
        b = _run_one(args.before_bin, case)
        a = _run_one(args.after_bin, case)
        if b["sha"] != before_sha or a["sha"] != after_sha:
            print(
                f"ERROR: SHA changed during run on case {case}: "
                f"before={b['sha']!r} (expected {before_sha!r}), "
                f"after={a['sha']!r} (expected {after_sha!r})",
                file=sys.stderr,
            )
            return 2
        bv, av = b["vals"], a["vals"]
        assert len(bv) == len(av) and len(bv) > 0
        sign_diffs = sum(
            1 for x, y in zip(bv, av)
            if (x != y) or (math.copysign(1.0, x) != math.copysign(1.0, y))
        )
        peak = max(max(abs(x) for x in bv), 1e-30)
        max_rel = max(abs(x - y) for x, y in zip(bv, av)) / peak
        ok = (b["md5"] == a["md5"]) and (sign_diffs == 0) and (peak > 0.0)
        if not ok:
            any_fail = True
        rows.append((case, len(bv), b["md5"], a["md5"], sign_diffs, max_rel, ok))

    counts_summary = ", ".join(f"{case}={n}" for case, n, *_ in rows)
    lines = [
        f"BEFORE SHA: `{before_sha}`",
        f"AFTER  SHA: `{after_sha}`",
        f"Value counts per case: {counts_summary}",
        "",
        "| Case | Values | MD5 BEFORE | MD5 AFTER | Signbit diffs | Max rel diff | Status |",
        "|---|---:|---|---|---:|---:|---|",
    ]
    for case, n, mb, ma, sd, mr, ok in rows:
        lines.append(
            f"| `{case}` | {n} | `{mb}` | `{ma}` | `{sd}` | `{mr:.3e}` | {'PASS (G-BYTE)' if ok else 'FAIL'} |"
        )
    text = "\n".join(lines) + "\n"
    print(text)
    if args.out:
        with open(args.out, "w", encoding="utf-8") as fh:
            fh.write(text)
    return 1 if any_fail else 0


if __name__ == "__main__":
    raise SystemExit(main())
