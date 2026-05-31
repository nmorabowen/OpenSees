"""Build a tiny, spec-conformant .ladruno file in pure Python.

Two purposes:
  1. Self-test fixture -- lets the validator / reader / reconstruction be tested today,
     before the C++ recorder exists.
  2. Reference artifact -- the C++ implementers target *this exact layout*. If the real
     recorder's output diverges from this file, one of them is wrong.

Model (2-D): a unit-square FourNodeQuad (stress) + a 2-node DispBeamColumn2d with a
3-fiber section (section.fiber.stress), exercising BASIS/QUADRATURE, LOCAL_AXES,
SECTION_ASSIGNMENTS(fiber), SECTION_MAP, and COLUMN_MAP with a multiplicity>1 block.
"""

from __future__ import annotations

import h5py
import numpy as np

INV_SQRT3 = 1.0 / np.sqrt(3.0)


def _s(g: h5py.Group, name: str, value: str) -> None:
    g.attrs[name] = np.bytes_(value.encode("utf-8"))


def _chunked(g: h5py.Group, slabs, steps, times) -> None:
    """Write the canonical chunked time-series layout (schema D3): a single
    DATA[T×nIds×nComp] extensible dataset + STEP[T] + TIME[T] axes. `slabs` is a
    list of [nIds×nComp] arrays, one per step."""
    data = np.stack([np.atleast_2d(a) for a in slabs], axis=0)
    g.create_dataset("DATA", data=data, maxshape=(None,) + data.shape[1:],
                     chunks=True, shuffle=True, compression="gzip", compression_opts=4)
    g.create_dataset("STEP", data=np.asarray(steps, dtype=np.int32), maxshape=(None,))
    g.create_dataset("TIME", data=np.asarray(times, dtype=np.float64), maxshape=(None,))


def build(path: str) -> None:
    with h5py.File(path, "w") as f:
        info = f.create_group("INFO")
        _s(info, "GENERATOR", "MPCO_Ladruno")
        info.attrs["FORMAT_VERSION"] = 1
        _s(info, "SOLVER_NAME", "OpenSees")
        info.attrs["SOLVER_VERSION"] = np.array([3, 5, 1], dtype=np.int32)
        info.attrs["SPATIAL_DIM"] = 2
        prov = info.create_group("PROVENANCE")
        _s(prov, "MODEL_HASH", "")
        _s(prov, "RESULTS_HASH", "")

        st = f.create_group("MODEL_STAGE[0]")
        st.attrs["STEP"] = 0
        st.attrs["TIME"] = 0.0
        _s(st, "KIND", "static")

        # ---- MODEL/NODES (quad 1..4 unit square, beam 4-5 vertical) ----
        model = st.create_group("MODEL")
        nodes = model.create_group("NODES")
        node_ids = np.array([1, 2, 3, 4, 5], dtype=np.int64)
        coords = np.array(
            [[0, 0], [1, 0], [1, 1], [0, 1], [0, 2]], dtype=np.float64
        )
        nodes.create_dataset("ID", data=node_ids)
        nodes.create_dataset("COORDINATES", data=coords)

        # ---- MODEL/ELEMENTS ----
        elems = model.create_group("ELEMENTS")
        # quad: lagrange order [1,1], 4 GP at +/- 1/sqrt3
        quad = elems.create_group("156-FourNodeQuad[0]")
        quad.create_dataset("CONNECTIVITY", data=np.array([[10, 1, 2, 3, 4]], dtype=np.int64))
        _s(quad, "TOPOLOGY", "quad")
        _s(quad, "FAMILY", "lagrange")
        quad.attrs["ORDER"] = np.array([1, 1], dtype=np.int32)
        _s(quad, "PARAM_DOMAIN", "[-1,1]")
        quad.attrs["RATIONAL"] = 0
        quad.attrs["NUM_CTRL"] = 4
        quad.attrs["NUM_GP"] = 4
        gpq = quad.create_group("QUADRATURE")
        gpq.create_dataset(
            "GP_PARAM",
            data=np.array(
                [[-INV_SQRT3, -INV_SQRT3], [INV_SQRT3, -INV_SQRT3],
                 [INV_SQRT3, INV_SQRT3], [-INV_SQRT3, INV_SQRT3]],
                dtype=np.float64,
            ),
        )
        gpq.create_dataset("GP_WEIGHT", data=np.ones(4, dtype=np.float64))

        # beam: line lagrange order [1], 2 GP
        beam = elems.create_group("28-DispBeamColumn2d[0]")
        beam.create_dataset("CONNECTIVITY", data=np.array([[20, 4, 5]], dtype=np.int64))
        _s(beam, "TOPOLOGY", "line")
        _s(beam, "FAMILY", "lagrange")
        beam.attrs["ORDER"] = np.array([1], dtype=np.int32)
        _s(beam, "PARAM_DOMAIN", "[-1,1]")
        beam.attrs["RATIONAL"] = 0
        beam.attrs["NUM_CTRL"] = 2
        beam.attrs["NUM_GP"] = 2
        gpb = beam.create_group("QUADRATURE")
        gpb.create_dataset(
            "GP_PARAM", data=np.array([[-INV_SQRT3], [INV_SQRT3]], dtype=np.float64)
        )
        gpb.create_dataset("GP_WEIGHT", data=np.ones(2, dtype=np.float64))

        # ---- MODEL/LOCAL_AXES (beam only; vertical -> local x = global y) ----
        la = model.create_group("LOCAL_AXES")
        la_beam = la.create_group("28-DispBeamColumn2d[0]")
        la_beam.create_dataset("ID", data=np.array([20], dtype=np.int64))
        # quaternion for the rotation taking global->local (placeholder valid unit quat)
        la_beam.create_dataset(
            "FRAME", data=np.array([[0.70710678, 0.0, 0.0, 0.70710678]], dtype=np.float64)
        )

        # ---- MODEL/SECTION_ASSIGNMENTS (fiber section on the beam) ----
        secs = model.create_group("SECTION_ASSIGNMENTS")
        sec1 = secs.create_group("SECTION_1[FiberSection2d]")
        sec1.attrs["ID"] = 1
        _s(sec1, "NAME", "FiberSection2d")
        _s(sec1, "KIND", "fiber")
        sec1.create_dataset("ASSIGNMENT", data=np.array([[20, 0], [20, 1]], dtype=np.int64))
        sec1.create_dataset(
            "FIBER_DATA",
            data=np.array([[-0.1, 0, 0.01], [0.0, 0, 0.01], [0.1, 0, 0.01]], dtype=np.float64),
        )
        sec1.create_dataset("FIBER_MATERIALS", data=np.array([1, 1, 1], dtype=np.int64))

        # ---- RESULTS ----
        res = st.create_group("RESULTS")

        # ON_NODES/DISPLACEMENT, 2 steps
        on_n = res.create_group("ON_NODES")
        disp = on_n.create_group("DISPLACEMENT")
        _s(disp, "DISPLAY_NAME", "Displacement")
        _s(disp, "COMPONENTS", "Ux,Uy")
        _s(disp, "DIMENSION", "L")
        _s(disp, "DESCRIPTION", "")
        disp.attrs["TYPE"] = 0
        disp.attrs["DATA_TYPE"] = 1
        disp.create_dataset("ID", data=node_ids)
        _chunked(
            disp,
            [(k + 1) * 1e-3 * np.arange(1, 11, dtype=np.float64).reshape(5, 2) for k in (0, 1)],
            steps=[0, 1], times=[0.0, 1.0])

        # ON_ELEMENTS/stress (quad): 4 GP x 3 comps = 12 cols, multiplicity 1
        on_e = res.create_group("ON_ELEMENTS")
        sb = on_e.create_group("stress").create_group("156-FourNodeQuad[0:0]")
        sb.attrs["NUM_COLUMNS"] = 12
        cm = sb.create_group("COLUMN_MAP")
        cm.create_dataset("LEVELS", data=np.array([[0, 1]] * 4, dtype=np.int32))
        cm.create_dataset("GAUSS_ID", data=np.array([0, 1, 2, 3], dtype=np.int64))
        cm.create_dataset("SECTION_TAG", data=np.array([-1, -1, -1, -1], dtype=np.int64))
        cm.create_dataset("FIBER_ID", data=np.array([-1, -1, -1, -1], dtype=np.int64))
        cm.create_dataset("NUM_COMP", data=np.array([3, 3, 3, 3], dtype=np.int64))
        cm.create_dataset("MULTIPLICITY", data=np.array([1, 1, 1, 1], dtype=np.int64))
        _s(cm, "COMP_NAMES", "\n".join(["sigma_xx,sigma_yy,sigma_xy"] * 4))
        sb.create_dataset("SECTION_MAP", data=np.array([[-1, -1, -1, -1]], dtype=np.int64))
        sb.create_dataset("ID", data=np.array([10], dtype=np.int64))
        _chunked(sb, [np.arange(12, dtype=np.float64).reshape(1, 12)], steps=[0], times=[0.0])

        # ON_ELEMENTS/section.fiber.stress (beam): 2 GP x (1 sec x 3 fibers) x 1 comp = 6
        fb = on_e.create_group("section.fiber.stress").create_group("28-DispBeamColumn2d[1000:1:0]")
        fb.attrs["NUM_COLUMNS"] = 6
        cmf = fb.create_group("COLUMN_MAP")
        # one block per GP; multiplicity 3 compresses the 3 identical fibers
        cmf.create_dataset("LEVELS", data=np.array([[0, 1, 2, 3], [0, 1, 2, 3]], dtype=np.int32))
        cmf.create_dataset("GAUSS_ID", data=np.array([0, 1], dtype=np.int64))
        cmf.create_dataset("SECTION_TAG", data=np.array([1, 1], dtype=np.int64))
        cmf.create_dataset("FIBER_ID", data=np.array([0, 0], dtype=np.int64))  # start fiber
        cmf.create_dataset("NUM_COMP", data=np.array([1, 1], dtype=np.int64))
        cmf.create_dataset("MULTIPLICITY", data=np.array([3, 3], dtype=np.int64))
        _s(cmf, "COMP_NAMES", "sigma\nsigma")
        fb.create_dataset("SECTION_MAP", data=np.array([[1, 1]], dtype=np.int64))
        fb.create_dataset("ID", data=np.array([20], dtype=np.int64))
        _chunked(fb, [np.arange(6, dtype=np.float64).reshape(1, 6) * 10.0], steps=[0], times=[0.0])


if __name__ == "__main__":
    import sys

    out = sys.argv[1] if len(sys.argv) > 1 else "synthetic.ladruno"
    build(out)
    print(f"wrote {out}")
