"""Parse a profiler h5 -> phase rollup + elem_by_type. Usage: python parse_prof.py <file.h5>"""
import os, sys
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
import h5py

path = sys.argv[1]
with h5py.File(path, "r") as f:
    run = f["runs"][list(f["runs"])[0]]
    print("META:", {k: run.attrs[k] for k in run.attrs})
    def walk(g, p="", depth=0):
        for k in g:
            it = g[k]
            if isinstance(it, h5py.Group):
                a = dict(it.attrs)
                pp = p + "/" + k if p else k
                if "wall_ns" in a:
                    print(f"{'  '*depth}{pp}: {int(a['wall_ns'])/1e6:.2f} ms / {int(a.get('calls',0))} calls")
                walk(it, pp, depth + 1)
            elif k == "elem_by_type":
                for row in it[()]:
                    print(f"{'  '*depth}[{p}] classTag={int(row['classTag'])} count={int(row['count'])} "
                          f"wall={int(row['wall_ns'])/1e6:.2f} ms  {float(row['wall_ns_per_elem'])/1e3:.2f} us/ele")
    if "rollup" in run and "root" in run["rollup"]:
        walk(run["rollup"]["root"], "root")
