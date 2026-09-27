"""WP-134 attribution probe: the C++ ModifiedEuler per-substep trace (T, dT, err,
code) on the cross-check outliers.  Shows the one-substep, err = 0 acceptance.

Runs under the fork's CPython 3.12 (-S), NOT the reference interpreter:
    <py312> -S Ladruno_files/testbed/sanisand_reference/probe_cpp_trace.py
Cases are written as literal states (from sanisand_reference.crosscheck's
benign_states: on the yield surface, alpha_in = 0, z = 0, e = 0.72)."""
import json
import math
import os
import sys

BIN = os.environ.get("SANISAND_REF_OPENSEES_BIN",
                     r"C:\Users\nmora\Github\OpenSees_Compile\OpenSees\.claude\worktrees"
                     r"\tims-implementation-review-3733c6\dist\bin")
os.add_dll_directory(BIN)
sys.path.insert(0, BIN)
sys.path.append(r"C:\Users\nmora\AppData\Local\Python\pythoncore-3.12-64\Lib\site-packages")
import opensees as ops  # noqa: E402

assert os.path.normcase(os.path.dirname(os.path.abspath(ops.__file__))) == os.path.normcase(BIN)

P = [264.32, 0.312885, 0.6944, 1.3309, 0.71, 0.027, 0.83, 0.45, 101.0, 0.005,
     1.3, 0.968, 3.5, 0.05, 5.75, 12.5, 1100.0, 2.0]
M, SQ23 = 0.005, math.sqrt(2.0 / 3.0)


def on_yield(p, eta, d):
    # d: deviatoric direction as (xx, yy, zz, xy); unit-normalised here
    tr = (d[0] + d[1] + d[2]) / 3.0
    v = [d[0] - tr, d[1] - tr, d[2] - tr, d[3]]
    nn = math.sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2 + 2 * v[3] ** 2)
    n = [x / nn for x in v]
    s = [eta * p * SQ23 * x for x in n]
    sig = [p + s[0], p + s[1], p + s[2], s[3], 0.0, 0.0]
    al = [s[i] / p - SQ23 * M * n[i] for i in range(4)] + [0.0, 0.0]
    return sig, al


CASES = {
    "p20_TE txLoad 1e-4": (20.0, 0.4 * 1.3309 * 0.71, (-2, 1, 1, 0), [1e-4, 0, 0, 0, 0, 0]),
    "p20_TC isoComp 1e-5": (20.0, 0.5 * 1.3309, (2, -1, -1, 0), [1e-5, 1e-5, 0, 0, 0, 0]),
    "p100_TCshear txLoad 1e-5": (100.0, 0.8 * 1.3309, (1, -0.5, -0.5, 0.6), [1e-5, 0, 0, 0, 0, 0]),
}


def main():
    ops.wipe()
    ops.nDMaterial("LadrunoSANISAND", 1, *P, 1, 0, 1, 1.0e-7, 1.0e-8, "-Presidual", 0.0,
                   "-Pmin", 0.0101, "-maxSubsteps", 200000, "-honorTolR", 1,
                   "-flipAlphaIn", "init")
    out = {}
    for name, (p, eta, d, de) in CASES.items():
        sig, al = on_yield(p, eta, d)
        r = list(ops.ladrunoSANISANDReplay(1, "-convention", "compressionPositive", "-sigma", *sig,
                                           "-alpha", *al, "-alphaIn", *[0.0] * 6, "-fabric", *[0.0] * 6,
                                           "-voidRatio", 0.72, "-dStrain", *de, "-trace", 50))
        nst, nrec, w = int(r[2]), int(r[3]), int(r[4])
        base = 6 + nst + 34
        recs = [r[base + k * w: base + (k + 1) * w] for k in range(nrec)]
        path = r[6 + nst + 29]
        out[name] = dict(path=path, substeps=r[6 + 2], trace=[dict(T=x[0], dT=x[1], err=x[2], code=int(x[3])) for x in recs])
        print(name, "path", path, "substeps", r[8], "trace", [(round(x[0], 6), x[1], x[2], int(x[3])) for x in recs[:8]])
    here = os.path.dirname(os.path.abspath(__file__))
    json.dump(out, open(os.path.join(here, "out", "cpp_trace_probe.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
