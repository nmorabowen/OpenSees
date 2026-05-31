"""MPCO_Ladruno (.ladruno) reader + schema-v1 validator + canonical normalizers.

This module is the **executable form of the schema spec** (mpco_ladruno_schema_v1.md).
Three jobs:

1. ``LadrunoReader`` -- walk a .ladruno HDF5 file per schema v1.
2. ``validate`` -- assert every structural invariant the spec promises (L2 of the test
   plan). Returns a list of problems; empty == conformant.
3. ``normalize_nodal`` / ``normalize_element`` -- flatten results to canonical
   ``{(locator...): value}`` dicts so the L1 parity gate can diff a .ladruno against a
   legacy .mpco *by information, not layout*.

Pure h5py + numpy. No OpenSees, no build. When the C++ recorder emits its first file,
this validates it immediately.
"""

from __future__ import annotations

import sys
from typing import Any

import h5py
import numpy as np

import ladruno_basis

GENERATOR = "MPCO_Ladruno"
FORMAT_VERSION = 1

# COLUMN_MAP/LEVELS path codes (schema §7.2)
LEVEL_ELEMENT, LEVEL_GAUSS, LEVEL_SECTION, LEVEL_FIBER, LEVEL_MATERIAL = 0, 1, 2, 3, 4


def _attr(obj: h5py.HLObject, name: str) -> Any:
    v = obj.attrs[name]
    if isinstance(v, bytes):
        return v.decode("utf-8")
    if isinstance(v, np.ndarray):
        if v.dtype.kind == "S":
            decoded = [x.decode("utf-8") for x in v.reshape(-1)]
            return decoded[0] if v.size == 1 else decoded
        # The MPCO-style h5 wrapper writes "scalar" int/double attributes as
        # 1-element 1-D arrays. Squeeze those back to a Python scalar so callers
        # can int()/float() them; keep genuine vectors (size>1) as arrays.
        if v.size == 1:
            return v.reshape(-1)[0].item()
    return v


def _decode_lines(raw: Any) -> list[str]:
    """COMP_NAMES is newline-separated, one line per COLUMN_MAP block."""
    if isinstance(raw, bytes):
        raw = raw.decode("utf-8")
    return [ln for ln in str(raw).split("\n")]


# ---------------------------------------------------------------------------
# Reader
# ---------------------------------------------------------------------------


class LadrunoReader:
    """Minimal schema-v1 reader. Open with a path; use as a context manager."""

    def __init__(self, path: str):
        self.path = path
        self.f = h5py.File(path, "r")

    def __enter__(self) -> "LadrunoReader":
        return self

    def __exit__(self, *exc) -> None:
        self.f.close()

    # -- INFO --
    def info(self) -> dict:
        g = self.f["INFO"]
        out = {k: _attr(g, k) for k in g.attrs}
        return out

    # -- stages --
    def stages(self) -> list[str]:
        return [k for k in self.f if k.startswith("MODEL_STAGE")]

    def stage_attrs(self, stage: str) -> dict:
        g = self.f[stage]
        return {k: _attr(g, k) for k in g.attrs}

    # -- model --
    def nodes(self, stage: str) -> tuple[np.ndarray, np.ndarray]:
        g = self.f[f"{stage}/MODEL/NODES"]
        return g["ID"][...], g["COORDINATES"][...]

    def element_groups(self, stage: str) -> dict[str, dict]:
        """name -> {connectivity, topology, family, order, param_domain, rational,
        num_ctrl, ndir, num_gp, gp_param, gp_weight, ctrl_weight?, global_gp_coords?}.

        QUADRATURE-TOLERANT: a group whose integration rule is not tabulated emits
        no QUADRATURE child -- then gp_param/gp_weight/num_gp/ndir are None (the
        caller warns, it is not an error). NDIR (schema-v1 D4) is authoritative for
        the parametric dimension, decoupled from len(ORDER); it falls back to the
        GP_PARAM column count for older files that predate the attribute.
        GLOBAL_GP_COORDS (D2) is optional, stored 2-D [nElem x (NUM_GP*ndim)]."""
        out: dict[str, dict] = {}
        base = self.f[f"{stage}/MODEL/ELEMENTS"]
        for name in base:
            grp = base[name]
            d = {
                "connectivity": grp["CONNECTIVITY"][...],
                "topology": _attr(grp, "TOPOLOGY"),
                "family": _attr(grp, "FAMILY"),
                "order": tuple(int(x) for x in np.atleast_1d(_attr(grp, "ORDER")))
                if "ORDER" in grp.attrs else (),
                "param_domain": _attr(grp, "PARAM_DOMAIN"),
                "rational": int(_attr(grp, "RATIONAL")),
                "num_ctrl": int(_attr(grp, "NUM_CTRL")),
            }
            if "QUADRATURE/GP_PARAM" in grp:
                gp_param = grp["QUADRATURE/GP_PARAM"][...]
                d["gp_param"] = gp_param
                d["num_gp"] = (int(_attr(grp, "NUM_GP")) if "NUM_GP" in grp.attrs
                               else gp_param.shape[0])
                d["ndir"] = (int(_attr(grp, "NDIR")) if "NDIR" in grp.attrs
                             else gp_param.shape[1])
                d["gp_weight"] = (grp["QUADRATURE/GP_WEIGHT"][...]
                                  if "QUADRATURE/GP_WEIGHT" in grp else None)
            else:
                d["gp_param"] = None
                d["num_gp"] = None
                d["ndir"] = None
                d["gp_weight"] = None
            if "CTRL_WEIGHT" in grp:
                d["ctrl_weight"] = grp["CTRL_WEIGHT"][...]
            if "GLOBAL_GP_COORDS" in grp:
                d["global_gp_coords"] = grp["GLOBAL_GP_COORDS"][...]
            out[name] = d
        return out

    def local_axes(self, stage: str) -> dict[str, dict]:
        out: dict[str, dict] = {}
        path = f"{stage}/MODEL/LOCAL_AXES"
        if path not in self.f:
            return out
        for name in self.f[path]:
            grp = self.f[path][name]
            out[name] = {"id": grp["ID"][...], "frame": grp["FRAME"][...]}
        return out

    def section_assignments(self, stage: str) -> dict[str, dict]:
        out: dict[str, dict] = {}
        path = f"{stage}/MODEL/SECTION_ASSIGNMENTS"
        if path not in self.f:
            return out
        for name in self.f[path]:
            grp = self.f[path][name]
            d = {"kind": _attr(grp, "KIND"), "assignment": grp["ASSIGNMENT"][...]}
            if "FIBER_DATA" in grp:
                d["fiber_data"] = grp["FIBER_DATA"][...]
                d["fiber_materials"] = grp["FIBER_MATERIALS"][...]
            out[name] = d
        return out


# ---------------------------------------------------------------------------
# COLUMN_MAP decode (schema §7.2)
# ---------------------------------------------------------------------------


def decode_columns(cmap_grp: h5py.Group) -> list[dict]:
    """Expand COLUMN_MAP into a flat per-column list:
    [{col, gauss_id, section_tag, fiber_id, comp}].

    A block i covers MULTIPLICITY[i] consecutive items (incrementing FIBER_ID, or
    GAUSS_ID if no fiber level) each with NUM_COMP[i] components named in COMP_NAMES[i].
    """
    levels = cmap_grp["LEVELS"][...]
    gauss = cmap_grp["GAUSS_ID"][...]
    sect = cmap_grp["SECTION_TAG"][...]
    fiber = cmap_grp["FIBER_ID"][...]
    ncomp = cmap_grp["NUM_COMP"][...]
    mult = cmap_grp["MULTIPLICITY"][...]
    names = _decode_lines(_attr(cmap_grp, "COMP_NAMES"))

    cols: list[dict] = []
    col = 0
    for i in range(len(ncomp)):
        comp_names = [c for c in names[i].split(",") if c]
        has_fiber = LEVEL_FIBER in np.atleast_1d(levels[i]) if levels.ndim > 1 else (
            int(levels[i]) == LEVEL_FIBER
        )
        for m in range(int(mult[i])):
            fib = int(fiber[i]) + m if int(fiber[i]) >= 0 else -1
            gp = int(gauss[i]) + (m if (not _has_fiber(levels, i) and int(gauss[i]) >= 0) else 0)
            for c in range(int(ncomp[i])):
                cols.append(
                    {
                        "col": col,
                        "gauss_id": gp,
                        "section_tag": int(sect[i]),
                        "fiber_id": fib,
                        "comp": comp_names[c] if c < len(comp_names) else f"c{c}",
                    }
                )
                col += 1
    return cols


def _has_fiber(levels: np.ndarray, i: int) -> bool:
    row = levels[i]
    return LEVEL_FIBER in np.atleast_1d(row)


def iter_step_slices(grp):
    """Yield (step_id:int, arr2d[nIds×nComp]) for a result/bucket group, handling
    BOTH layouts transparently:

    * chunked (schema D3, .ladruno): DATA is a single [T×nIds×nComp] dataset with
      a sibling STEP[T] axis giving each slab's commit/step id;
    * legacy (frozen .mpco): DATA is a group of per-step STEP_<k> datasets, each
      carrying a scalar STEP attribute.

    This lets every reader/checker diff a chunked .ladruno against a per-step
    .mpco with one code path."""
    data = grp["DATA"]
    if isinstance(data, h5py.Dataset):
        full = np.asarray(data[...])
        steps = (np.asarray(grp["STEP"][...]).reshape(-1) if "STEP" in grp
                 else np.arange(full.shape[0]))
        for t in range(full.shape[0]):
            yield int(steps[t]), np.atleast_2d(full[t])
    else:
        for sname in data:
            ds = data[sname]
            yield int(_attr(ds, "STEP")), np.atleast_2d(ds[...])


def _check_data_shape(grp, n_rows: int, n_cols: int, err, ctx: str) -> None:
    """Validate the DATA payload shape for both layouts."""
    data = grp["DATA"]
    if isinstance(data, h5py.Dataset):
        if data.ndim != 3 or tuple(data.shape[1:]) != (n_rows, n_cols):
            err(f"{ctx}: DATA {data.shape} != (T, {n_rows}, {n_cols})")
        T = data.shape[0]
        for ax in ("STEP", "TIME"):
            if ax in grp and np.asarray(grp[ax][...]).reshape(-1).shape[0] != T:
                err(f"{ctx}: {ax} length != T={T}")
    else:
        for sn in data:
            arr = data[sn]
            if arr.shape != (n_rows, n_cols):
                err(f"{ctx}/{sn}: {arr.shape} != ({n_rows}, {n_cols})")


# ---------------------------------------------------------------------------
# Canonical normalizers (for the L1 parity diff)
# ---------------------------------------------------------------------------


def normalize_nodal(reader: LadrunoReader, stage: str) -> dict[tuple, float]:
    """{(node_tag, component_name, step): value} over all ON_NODES results."""
    out: dict[tuple, float] = {}
    base_path = f"{stage}/RESULTS/ON_NODES"
    if base_path not in reader.f:
        return out
    base = reader.f[base_path]
    for rname in base:
        grp = base[rname]
        ids = np.asarray(grp["ID"][...]).reshape(-1)
        comps = [c for c in _attr(grp, "COMPONENTS").split(",") if c]
        for k, arr in iter_step_slices(grp):
            for r, tag in enumerate(ids):
                for c, cn in enumerate(comps):
                    out[(int(tag), f"{rname}:{cn}", k)] = float(arr[r, c])
    return out


def normalize_element(reader: LadrunoReader, stage: str) -> dict[tuple, float]:
    """{(elem_tag, gauss_id, section_tag, fiber_id, component, step): value}."""
    out: dict[tuple, float] = {}
    base_path = f"{stage}/RESULTS/ON_ELEMENTS"
    if base_path not in reader.f:
        return out
    base = reader.f[base_path]
    for rname in base:
        rgrp = base[rname]
        for bname in rgrp:
            bucket = rgrp[bname]
            if "COLUMN_MAP" not in bucket:  # e.g. LOCAL_AXES time series -> skip here
                continue
            cols = decode_columns(bucket["COLUMN_MAP"])
            ids = np.asarray(bucket["ID"][...]).reshape(-1)
            for k, arr in iter_step_slices(bucket):
                for r, tag in enumerate(ids):
                    for cd in cols:
                        out[
                            (
                                int(tag),
                                cd["gauss_id"],
                                cd["section_tag"],
                                cd["fiber_id"],
                                f"{rname}:{cd['comp']}",
                                k,
                            )
                        ] = float(arr[r, cd["col"]])
    return out


# ---------------------------------------------------------------------------
# Validator (L2) -- the executable spec
# ---------------------------------------------------------------------------


def validate(path: str) -> list[str]:
    """Return a list of schema-v1 violations. Empty list == conformant file."""
    problems: list[str] = []

    def err(msg: str) -> None:
        problems.append(msg)

    with LadrunoReader(path) as r:
        info = r.info()
        if info.get("GENERATOR") != GENERATOR:
            err(f"INFO/GENERATOR must be {GENERATOR!r}, got {info.get('GENERATOR')!r}")
        if int(info.get("FORMAT_VERSION", -1)) != FORMAT_VERSION:
            err(f"INFO/FORMAT_VERSION must be {FORMAT_VERSION}, got {info.get('FORMAT_VERSION')}")
        ndim = int(info.get("SPATIAL_DIM", 0))
        if ndim not in (1, 2, 3):
            err(f"INFO/SPATIAL_DIM must be 1|2|3, got {ndim}")

        # Partition manifest (parallel output). Only present/meaningful when this is
        # one ".part-N.ladruno" of a partitioned set; PARTITION_ID must index into
        # the 0-based-contiguous NUM_PARTITIONS the apeGmsh stitcher globs for.
        if int(info.get("PARTITIONED", 0)):
            pid = int(info.get("PARTITION_ID", -1))
            nparts = int(info.get("NUM_PARTITIONS", 0))
            if not (0 <= pid < nparts):
                err(f"INFO: PARTITION_ID {pid} not in [0, NUM_PARTITIONS={nparts})")

        stages = r.stages()
        if not stages:
            err("no MODEL_STAGE[*] group present")

        for stage in stages:
            sa = r.stage_attrs(stage)
            if sa.get("KIND") not in ("transient", "static", "eigen"):
                err(f"{stage}: KIND must be transient|static|eigen, got {sa.get('KIND')!r}")

            nid, ncoord = r.nodes(stage)
            if ncoord.shape[0] != nid.shape[0]:
                err(f"{stage}: NODES COORDINATES rows {ncoord.shape[0]} != ID {nid.shape[0]}")
            if ncoord.ndim != 2 or ncoord.shape[1] != ndim:
                err(f"{stage}: NODES COORDINATES must be [n x {ndim}], got {ncoord.shape}")

            # element groups: BASIS + QUADRATURE invariants
            groups = r.element_groups(stage)
            for name, g in groups.items():
                if g["connectivity"].shape[1] != 1 + g["num_ctrl"]:
                    err(f"{stage}/{name}: CONNECTIVITY cols {g['connectivity'].shape[1]} "
                        f"!= 1+NUM_CTRL ({1 + g['num_ctrl']})")
                if g["gp_param"] is None:
                    # QUADRATURE-tolerant (schema-v1 decision): an untabulated rule
                    # legitimately emits no QUADRATURE. Warn (stderr), do not fail.
                    print(f"WARN {stage}/{name}: no QUADRATURE group (untabulated "
                          f"rule); skipping GP checks", file=sys.stderr)
                else:
                    # NDIR is authoritative for the parametric dim (D4), NOT len(ORDER)
                    # -- simplices have ndir 2/3 while ORDER stays total-degree.
                    ndir = g["ndir"]
                    if g["gp_param"].shape != (g["num_gp"], ndir):
                        err(f"{stage}/{name}: GP_PARAM {g['gp_param'].shape} != "
                            f"(NUM_GP={g['num_gp']}, NDIR={ndir})")
                    if g["param_domain"] == "bary" and ndir not in (2, 3):
                        err(f"{stage}/{name}: bary PARAM_DOMAIN needs NDIR in (2,3), got {ndir}")
                    if g["gp_weight"] is not None and g["gp_weight"].shape[0] != g["num_gp"]:
                        err(f"{stage}/{name}: GP_WEIGHT len {g['gp_weight'].shape[0]} != NUM_GP")
                    # GLOBAL_GP_COORDS (D2) optional; if present, [nElem x NUM_GP*ndim].
                    if "global_gp_coords" in g:
                        nelem = g["connectivity"].shape[0]
                        want = (nelem, g["num_gp"] * ndim)
                        if g["global_gp_coords"].shape != want:
                            err(f"{stage}/{name}: GLOBAL_GP_COORDS {g['global_gp_coords'].shape}"
                                f" != (nElem, NUM_GP*ndim)={want}")
                if g["topology"] not in (
                    "line", "tri", "quad", "tet", "hex", "wedge", "pyramid", "custom"
                ):
                    err(f"{stage}/{name}: unknown TOPOLOGY {g['topology']!r}")
                if g["rational"]:
                    if "ctrl_weight" not in g:
                        err(f"{stage}/{name}: RATIONAL=1 but no CTRL_WEIGHT")
                    elif g["ctrl_weight"].shape != (
                        g["connectivity"].shape[0], g["num_ctrl"]
                    ):
                        err(f"{stage}/{name}: CTRL_WEIGHT {g['ctrl_weight'].shape} != "
                            f"(nElem, NUM_CTRL)")

            # local axes alignment
            for name, la in r.local_axes(stage).items():
                if la["frame"].shape != (la["id"].shape[0], 4):
                    err(f"{stage}/LOCAL_AXES/{name}: FRAME {la['frame'].shape} != (nE, 4)")

            # section assignments
            for name, s in r.section_assignments(stage).items():
                if s["kind"] not in ("fiber", "resultant"):
                    err(f"{stage}/SECTION_ASSIGNMENTS/{name}: KIND must be fiber|resultant")
                if s["kind"] == "fiber":
                    if "fiber_data" not in s:
                        err(f"{stage}/SECTION_ASSIGNMENTS/{name}: fiber kind needs FIBER_DATA")
                    elif s["fiber_data"].shape[0] != s["fiber_materials"].shape[0]:
                        err(f"{stage}/SECTION_ASSIGNMENTS/{name}: FIBER_DATA rows != "
                            f"FIBER_MATERIALS")

            _validate_element_results(r, stage, err)
            _validate_nodal_results(r, stage, err)

    return problems


# ---------------------------------------------------------------------------
# Round-trip ORACLE (D5) -- the write-time geometry self-check
# ---------------------------------------------------------------------------


def round_trip_oracle(path: str, tol: float = 1e-12) -> list[str]:
    """Mandatory write->read->reconstruct geometry check (schema-v1 D5).

    For every element group that carries BOTH a GP_PARAM table and the static
    ``GLOBAL_GP_COORDS`` (the belt-and-suspenders D2 dataset), independently
    reconstruct ``x(GP_PARAM[k]) = sum_i N_i(xi_k) X_i`` from the file's own
    MODEL/NODES + CONNECTIVITY using ``ladruno_basis.shape_functions`` and compare
    to ``GLOBAL_GP_COORDS[k]`` (<= ``tol``). This catches GP-ordering and basis
    bugs *at write time*: if the recorder's C++ ``computeGlobalGP`` and this
    independent Python basis disagree, one is wrong.

    Groups with no GLOBAL_GP_COORDS or no QUADRATURE are skipped (nothing to
    cross-check); a topology the oracle basis does not cover (e.g. a custom-rule
    bernstein/IGA element whose own gate validates it) is warned and skipped, not
    failed. Returns a list of problems; empty == oracle-clean.
    """
    problems: list[str] = []
    with LadrunoReader(path) as r:
        ndim = int(r.info().get("SPATIAL_DIM", 0))
        for stage in r.stages():
            nid, ncoord = r.nodes(stage)
            nid = np.asarray(nid).reshape(-1)
            coord_of = {int(t): np.asarray(ncoord[i], dtype=float) for i, t in enumerate(nid)}
            for name, g in r.element_groups(stage).items():
                if g.get("global_gp_coords") is None or g["gp_param"] is None:
                    continue
                topo = g["topology"]
                conn = g["connectivity"]
                num_gp = g["num_gp"]
                nelem = conn.shape[0]
                num_nodes = conn.shape[1] - 1
                ggp = np.asarray(g["global_gp_coords"], dtype=float)
                if ggp.size != nelem * num_gp * ndim:
                    problems.append(
                        f"{stage}/{name}: GLOBAL_GP_COORDS size {ggp.size} != "
                        f"nElem*NUM_GP*ndim ({nelem}*{num_gp}*{ndim})")
                    continue
                ggp3 = ggp.reshape(nelem, num_gp, ndim)
                max_err = 0.0
                try:
                    for e in range(nelem):
                        X = np.array([coord_of[int(t)][:ndim] for t in conn[e, 1:]])
                        for k in range(num_gp):
                            N = ladruno_basis.shape_functions(topo, num_nodes, g["gp_param"][k])
                            x_rec = N @ X
                            max_err = max(max_err, float(np.max(np.abs(x_rec - ggp3[e, k]))))
                except ValueError as exc:
                    print(f"WARN {stage}/{name}: oracle skipped ({exc})", file=sys.stderr)
                    continue
                if max_err > tol:
                    problems.append(
                        f"{stage}/{name}: round-trip x(GP_PARAM) vs GLOBAL_GP_COORDS "
                        f"max_err={max_err:.3e} > tol {tol:.0e}")
    return problems


def _validate_nodal_results(r: LadrunoReader, stage: str, err) -> None:
    path = f"{stage}/RESULTS/ON_NODES"
    if path not in r.f:
        return
    for rname in r.f[path]:
        grp = r.f[path][rname]
        ids = grp["ID"][...]
        ncomp = len([c for c in _attr(grp, "COMPONENTS").split(",") if c])
        _check_data_shape(grp, ids.shape[0], ncomp, err, f"{stage}/ON_NODES/{rname}")


def _validate_element_results(r: LadrunoReader, stage: str, err) -> None:
    path = f"{stage}/RESULTS/ON_ELEMENTS"
    if path not in r.f:
        return
    for rname in r.f[path]:
        rgrp = r.f[path][rname]
        for bname in rgrp:
            bucket = rgrp[bname]
            if "COLUMN_MAP" not in bucket:
                continue
            num_cols = int(_attr(bucket, "NUM_COLUMNS"))
            cm = bucket["COLUMN_MAP"]
            mult = cm["MULTIPLICITY"][...]
            ncomp = cm["NUM_COMP"][...]
            k = len(ncomp)
            for arr_name in ("LEVELS", "GAUSS_ID", "SECTION_TAG", "FIBER_ID", "MULTIPLICITY"):
                if len(cm[arr_name][...]) != k:
                    err(f"{stage}/ON_ELEMENTS/{rname}/{bname}: COLUMN_MAP/{arr_name} len "
                        f"!= {k}")
            if len(_decode_lines(_attr(cm, "COMP_NAMES"))) != k:
                err(f"{stage}/ON_ELEMENTS/{rname}/{bname}: COMP_NAMES lines != {k} blocks")
            expect = int(np.sum(mult * ncomp))
            if expect != num_cols:
                err(f"{stage}/ON_ELEMENTS/{rname}/{bname}: NUM_COLUMNS {num_cols} != "
                    f"sum(MULT*NUM_COMP) {expect}")
            ids = bucket["ID"][...]
            _check_data_shape(bucket, ids.shape[0], num_cols, err,
                              f"{stage}/ON_ELEMENTS/{rname}/{bname}")
            if "SECTION_MAP" in bucket:
                # SECTION_MAP is [nElem x NUM_GP]; resolve NUM_GP from the model group
                sm = bucket["SECTION_MAP"][...]
                if sm.shape[0] != ids.shape[0]:
                    err(f"{stage}/ON_ELEMENTS/{rname}/{bname}: SECTION_MAP rows "
                        f"{sm.shape[0]} != nElem {ids.shape[0]}")
