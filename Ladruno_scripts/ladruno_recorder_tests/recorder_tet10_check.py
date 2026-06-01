"""Recorder round-trip — BezierTet10 .ladruno checker (venv python + h5py).

Verifies the inverted-dispatch wiring produced a correctly self-described
element block: FAMILY=bernstein, topology=tet, NUM_CTRL=10, NUM_GP=4, the
exact 4×3 GP_PARAM and 4 GP_WEIGHT (=1/24).

    python recorder_tet10_check.py <out_dir>
"""
import sys
import os
import numpy as np
import h5py

OUT = sys.argv[1] if len(sys.argv) > 1 else "."
path = os.path.join(OUT, "test_tet10.ladruno")


def squeeze(v):
    """Ladruno h5 writes scalar attrs as 1-elem arrays; bytes for strings."""
    if isinstance(v, bytes):
        return v.decode("utf-8", "replace")
    a = np.asarray(v)
    if a.ndim >= 1 and a.size == 1:
        v = a.reshape(-1)[0]
        return v.decode("utf-8", "replace") if isinstance(v, bytes) else v
    if a.dtype.kind == "S":
        return [x.decode("utf-8", "replace") for x in a.reshape(-1)]
    return a


passed, failed = [], []


def check(name, ok, detail=""):
    (passed if ok else failed).append(name)
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}{('  — ' + detail) if detail else ''}")


with h5py.File(path, "r") as f:
    # GENERATOR provenance
    gen = None
    if "INFO" in f:
        for k in f["INFO"].attrs:
            if k.upper() == "GENERATOR":
                gen = squeeze(f["INFO"].attrs[k])
    check("INFO/GENERATOR == 'Ladruno'", gen == "Ladruno", f"got {gen!r}")

    # ladruno schema (PR #18): MODEL/ELEMENTS/<name> is a GROUP holding a
    # CONNECTIVITY dataset + a QUADRATURE child group (GP_PARAM, GP_WEIGHT);
    # BASIS attrs (FAMILY/TOPOLOGY/ORDER/...) live on the group.
    blk = {"grp": None}

    def visit(name, obj):
        base = name.split("/")[-1]
        if isinstance(obj, h5py.Group) and "BezierTet10" in base and blk["grp"] is None:
            blk["grp"] = obj

    f.visititems(visit)

    grp = blk["grp"]
    check("BezierTet10 element block present (GROUP)", grp is not None,
          grp.name if grp is not None else "")
    if grp is not None:
        a = {k: squeeze(v) for k, v in grp.attrs.items()}
        print(f"    attrs: {a}")
        check("TOPOLOGY == 'tet'", a.get("TOPOLOGY") == "tet", str(a.get("TOPOLOGY")))
        check("FAMILY == 'bernstein' (inverted dispatch)",
              a.get("FAMILY") == "bernstein", str(a.get("FAMILY")))
        check("RATIONAL == 0", int(a.get("RATIONAL", -1)) == 0, str(a.get("RATIONAL")))
        check("NUM_CTRL == 10", int(a.get("NUM_CTRL", -1)) == 10, str(a.get("NUM_CTRL")))
        check("NUM_GP == 4", int(a.get("NUM_GP", -1)) == 4, str(a.get("NUM_GP")))
        order = np.asarray(a.get("ORDER", [])).reshape(-1)
        check("ORDER == [2,2,2] (one per GP_PARAM column)",
              order.size == 3 and np.all(order == 2), str(order))
        check("PARAM_DOMAIN == 'bary'", a.get("PARAM_DOMAIN") == "bary",
              str(a.get("PARAM_DOMAIN")))

        # CONNECTIVITY child dataset: 1 elem × (tag + 10 nodes) = (1, 11)
        conn = grp.get("CONNECTIVITY")
        check("CONNECTIVITY child dataset (1, 11)",
              conn is not None and tuple(conn.shape) == (1, 11),
              str(None if conn is None else conn.shape))

        # QUADRATURE child group with GP_PARAM [4×3] + GP_WEIGHT [4]
        quad = grp.get("QUADRATURE")
        check("QUADRATURE child group present", isinstance(quad, h5py.Group),
              str(None if quad is None else quad.name))
        if isinstance(quad, h5py.Group):
            gpp = quad.get("GP_PARAM")
            check("GP_PARAM dataset (4, 3)",
                  gpp is not None and tuple(np.asarray(gpp).shape) == (4, 3),
                  str(None if gpp is None else np.asarray(gpp).shape))
            if gpp is not None and np.asarray(gpp).shape == (4, 3):
                g = np.asarray(gpp)
                a_, b_ = 0.585410196624968, 0.138196601125011
                ok_vals = True
                for r in g:
                    bary = np.append(r, 1.0 - r.sum())
                    na = np.sum(np.abs(bary - a_) < 1e-9)
                    nb = np.sum(np.abs(bary - b_) < 1e-9)
                    if not (na == 1 and nb == 3):
                        ok_vals = False
                check("GP_PARAM values are the 4-pt tet rule", bool(ok_vals),
                      f"rows~{np.round(g,4).tolist()}")
            gpw = quad.get("GP_WEIGHT")
            if gpw is not None:
                gpw = np.asarray(gpw).reshape(-1)
                check("GP_WEIGHT == [1/24]*4",
                      gpw.size == 4 and np.allclose(gpw, 1/24.0),
                      str(np.round(gpw, 6).tolist()))
            else:
                check("GP_WEIGHT dataset present", False)

print("=" * 60)
print(f"  {len(passed)} passed, {len(failed)} failed")
sys.exit(0 if not failed else 1)
