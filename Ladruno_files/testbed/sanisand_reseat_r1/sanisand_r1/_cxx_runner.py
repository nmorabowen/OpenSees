"""Subprocess side of the C++ cross-check (runs under the fork's CPython 3.12
with -S; imports NOTHING from sanisand_reference).

    python -S _cxx_runner.py <bin_dir> <site_packages> <jobs.json> <out.json>

jobs.json = {"params": [18 OpenSees material args],
             "protos": {"name": [IntScheme, TolR, honorTolR, maxSubsteps], ...},
             "pmin": 0.0101, "presidual": 0.0,
             "jobs": [{"proto": name, "sigma": [6], "alpha": [6], "alpha_in": [6],
                       "z": [6], "e": e, "deps": [6]}, ...]}
All states COMPRESSION POSITIVE (the model's internal convention), alpha/z
tensor components, deps Voigt with engineering shear.  Each job is ONE
`ladrunoSANISANDReplay` call (WP-127, F21) on a private copy of the prototype.
"""
import json
import os
import sys


def main():
    bin_dir, site, jobs_path, out_path = sys.argv[1:5]
    os.add_dll_directory(bin_dir)
    sys.path.insert(0, bin_dir)
    if site not in sys.path:
        sys.path.append(site)
    import opensees as ops
    here = os.path.normcase(os.path.dirname(os.path.abspath(ops.__file__)))
    assert here == os.path.normcase(os.path.abspath(bin_dir)), ops.__file__

    spec = json.load(open(jobs_path))
    P = spec["params"]
    ops.wipe()
    tags = {}
    for k, (name, (scheme, tolr, honor, maxsub)) in enumerate(sorted(spec["protos"].items())):
        tag = k + 1
        ops.nDMaterial("LadrunoSANISAND", tag, *P, int(scheme), 0, 1, 1.0e-7, float(tolr),
                       "-Presidual", float(spec.get("presidual", 0.0)),
                       "-Pmin", float(spec.get("pmin", 0.0101)),
                       "-maxSubsteps", int(maxsub), "-honorTolR", int(honor),
                       "-flipAlphaIn", "init")
        tags[name] = tag
    out = []
    for j in spec["jobs"]:
        args = [tags[j["proto"]], "-convention", "compressionPositive",
                "-sigma", *map(float, j["sigma"]), "-alpha", *map(float, j["alpha"]),
                "-alphaIn", *map(float, j["alpha_in"]), "-fabric", *map(float, j["z"]),
                "-voidRatio", float(j["e"]), "-dStrain", *map(float, j["deps"]),
                "-type", j.get("type", "3D"), "-trace", 0, "-dt", 1.0,
                "-primed", 1, "-prevIncrNorm", 0.0]
        r = ops.ladrunoSANISANDReplay(*args)
        if r is None or len(r) < 6 or int(r[0]) != 1:
            out.append(dict(ok=False, raw=list(r) if r else None))
            continue
        r = list(r)
        rc, nst = int(r[1]), int(r[2])
        stats = r[6:6 + nst]
        st = r[6 + nst:6 + nst + 34]
        out.append(dict(ok=True, rc=rc, substeps=stats[2], forced=stats[5],
                        abandoned=stats[8], cap=stats[9], pn_resets=stats[11],
                        entry_pmin=stats[10],
                        sigma=st[0:6], alpha=st[6:12], alpha_in=st[12:18], z=st[18:24],
                        e=st[24], p=st[25], q=st[26], f_before=st[27], f_after=st[28],
                        path=st[29]))
    json.dump(out, open(out_path, "w"))


if __name__ == "__main__":
    main()
