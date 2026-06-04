"""Bare element-verb robustness — check (L1/L2).

Asserts the bare-keyword guard:
  1. both bare-verb files were written (the model ran analyze() to completion,
     i.e. no segfault — the model runner failing would already fail the gate);
  2. each file is schema-conformant;
  3. neither file carries an ON_ELEMENTS bucket (a bare keyword with no sub-verb
     legitimately produces no element result);
  4. the valid nodal request DID record (so the recorder otherwise worked).

Pure h5py+numpy (runs under the venv check python).
    python shell_bare_verb_check.py <bare_section.ladruno> <bare_material.ladruno>
"""
from __future__ import annotations

import os
import sys

os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import ladruno_format as lf  # noqa: E402


def _element_buckets(reader, stage):
    base_path = f"{stage}/RESULTS/ON_ELEMENTS"
    if base_path not in reader.f:
        return 0
    n = 0
    base = reader.f[base_path]
    for rname in base:
        for bname in base[rname]:
            if "COLUMN_MAP" in base[rname][bname]:
                n += 1
    return n


def main(*paths: str) -> int:
    fails = []
    for p in paths:
        label = os.path.basename(p)
        if not os.path.exists(p):
            fails.append(f"{label}: file not written (model likely crashed)")
            continue
        problems = lf.validate(p)
        if problems:
            fails.append(f"{label}: schema problems: {problems[:3]}")
        with lf.LadrunoReader(p) as r:
            stage = r.stages()[0]
            nb = _element_buckets(r, stage)
            if nb != 0:
                fails.append(f"{label}: expected 0 element buckets for a bare verb, got {nb}")
            nodal = lf.normalize_nodal(r, stage)
            if len(nodal) == 0:
                fails.append(f"{label}: the valid -N displacement request recorded nothing")

    if fails:
        for f in fails:
            print(f"[FAIL] SHELL BARE VERB: {f}")
        return 1
    print("SHELL_BARE_VERB_CHECK: ALL PASS (bare -E section / -E material: no crash, "
          "no element bucket, nodal result intact)")
    return 0


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("usage: shell_bare_verb_check.py <bare_section.ladruno> <bare_material.ladruno>")
        sys.exit(2)
    sys.exit(main(sys.argv[1], sys.argv[2]))
