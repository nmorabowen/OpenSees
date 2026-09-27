"""WP-126 check: what does each layer report for the reaction at the shared support node 1?

  1. raw part files (h5py): the per-partition values the recorder STORED
  2. apeGmsh LadrunoReader on the serial (np=1) file: the reference
  3. apeGmsh LadrunoMultiPartitionReader over the np=2 part files: what a user gets
  4. the -envelope files (MAX / ABSMAX per partition)

Run with a Python that has apeGmsh + h5py (here: Python 3.11, editable apeGmsh):
    python repro_check.py <out_dir>
"""
from __future__ import annotations

import sys
from pathlib import Path

import h5py
import numpy as np
from apeGmsh.results.readers._ladruno import LadrunoReader
from apeGmsh.results.readers._ladruno_multi import LadrunoMultiPartitionReader
from apeGmsh.results.readers._protocol import ResultLevel

OUT = Path(sys.argv[1])
NODE = 1


def dump_raw(path):
    print(f"-- raw {path.name}")
    with h5py.File(path, "r") as f:
        def visit(name, obj):
            if isinstance(obj, h5py.Dataset) and ("REACTION" in name.upper() or "ENVELOPE" in name.upper()) \
                    and obj.size < 200:
                print(f"   {name}: {np.array2string(np.asarray(obj[()]), precision=6)}")
        f.visititems(visit)


def via_reader(reader, label):
    st = reader.stages()[0]
    sid = getattr(st, "id", getattr(st, "stage_id", st))
    comps = [c for c in reader.available_components(sid, ResultLevel.NODES) if "reaction" in str(c).lower()]
    out = {}
    for c in comps:
        sl = reader.read_nodes(sid, c, node_ids=np.array([NODE]))
        v = np.asarray(sl.values)
        out[c] = v[:, 0] if v.size else v
    print(f"-- {label}: " + "; ".join(f"{c} = {np.array2string(v, precision=6)}" for c, v in out.items()))
    return out


def main():
    for p in sorted(OUT.glob("*.ladruno")):
        dump_raw(p)
    ref = via_reader(LadrunoReader(OUT / "mp_np1_stream.ladruno"), "serial (np=1) reader")
    parts = sorted(OUT.glob("mp_np2_stream.part-*.ladruno"))
    got = via_reader(LadrunoMultiPartitionReader(parts), "stitched np=2 reader")
    print()
    bad = 0
    for c, rv in ref.items():
        gv = got.get(c)
        same = gv is not None and gv.shape == rv.shape and np.allclose(gv, rv, rtol=0, atol=1e-9)
        bad += 0 if same else 1
        print(f"{'OK  ' if same else 'DIFF'} {c}: serial {rv}  vs  stitched {gv}")
    print("\nREPRO:", "stitched reactions MATCH serial" if bad == 0 else f"{bad} component(s) DIFFER")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
