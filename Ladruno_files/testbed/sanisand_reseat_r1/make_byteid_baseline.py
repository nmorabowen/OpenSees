"""WP-151: record the SAS-ME byte-identity baseline on a build that does NOT carry
the WP-151 change (its SANISAND sources == origin/ladruno before this branch).

    py -3.12 -S make_byteid_baseline.py <bin_dir_without_wp151> [out.json]

Writes tests/data/wp151_sasme_byteid_baseline.json: the rows of
tests/wp151_reseat_tools.jobs() for prototype 11 (flags absent).  The test
compares the WP-151 build's prototypes 11 (absent) and 12 (given as 0) to it."""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
bin_dir = os.path.abspath(sys.argv[1])
out = sys.argv[2] if len(sys.argv) > 2 else os.path.join(ROOT, "tests", "data",
                                                         "wp151_sasme_byteid_baseline.json")
os.add_dll_directory(bin_dir)
sys.path.insert(0, bin_dir)
sys.path.insert(0, os.path.join(ROOT, "tests"))
import opensees as ops  # noqa: E402
assert os.path.normcase(os.path.dirname(os.path.abspath(ops.__file__))) == os.path.normcase(bin_dir), ops.__file__
import wp151_reseat_tools as T  # noqa: E402

ops.wipe()
T.define(ops, tags=[11])
rows = T.run_rows(ops, 11)
json.dump(dict(note="SAS-ME (IntScheme 129) replays, prototype 11 (the WP-138 E_B material, "
                    "no WP-151 flags), recorded on a build WITHOUT WP-151 (SANISAND sources == "
                    "origin/ladruno 6bd6905b3; SAS-ME last changed at beb6d8333)",
               platform=sys.platform, rows=rows), open(out, "w"), indent=0)
print("wrote", out, len(rows), "rows; refused:", sum(1 for r in rows.values() if r[0] != 0))
