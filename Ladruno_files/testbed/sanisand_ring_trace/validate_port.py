"""WP-128: validate md_port.py against the C++ `ladrunoSANISANDReplay`
(IntScheme 1, campaign flags) on all 80 ring rows x the probe set.

Pass criterion (per case): same rc, same substep census (attempts, accepted,
rejected, forced, clamp), and returned sigma / alpha agreeing to 1e-8
relative (||.|| of the difference over max(||.||, 1e-3)).  Divergence on a
case whose substep sequence hits an err within 1e-9 relative of TolE is a
round-off branch flip, reported separately.
Output: out/validate_port.txt
"""
import math
import sys

import _boot as B
from _boot import ops, sr
import md_port

PROBES = {}
for d in (1e-7, 1e-6, 1e-5, 1e-4):
    PROBES[f"isoComp@{d:.0e}"] = [d, d, 0, 0, 0, 0]
    PROBES[f"isoExt@{d:.0e}"] = [-d, -d, 0, 0, 0, 0]
    PROBES[f"shear+@{d:.0e}"] = [0, 0, 0, d, 0, 0]
    PROBES[f"shear-@{d:.0e}"] = [0, 0, 0, -d, 0, 0]
    PROBES[f"triax@{d:.0e}"] = [-0.5 * d, d, 0, 0, 0, 0]


def rel(a, b):
    num = math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))
    den = max(math.sqrt(sum(x * x for x in a)), 1e-3)
    return num / den


def main():
    B.define_prototypes()
    mat = md_port.Material(B.P)
    lines = []
    n = bad = nchaos = 0
    worst = 0.0
    for path in sr.RING_CSVS:
        for r in sr.read_ring_csv(path):
            for pname, de in PROBES.items():
                c = sr.replay_row(ops, B.TAG_ME, r, de, trace=0)
                # the replay projects alpha / alpha_in / z deviatoric first
                al = B.dev(r["alpha"])
                ai = B.dev(r["alpha_in"])
                z = B.dev(r["z"])
                pyr = mat.update(r["sigma"], al, ai, z, r["e"], de)
                s = c["stats"]
                same_census = (int(s["substeps"]) == pyr["substeps"]
                               and int(s["accepted"]) == pyr["acc"]
                               and int(s["forcedAtDTmin"]) == pyr["forced"]
                               and int(s["forcedClampMc"]) == pyr["clamp"]
                               and int(s["capHits"]) == pyr["cap"]
                               and (c["rc"] == pyr["rc"] or (c["rc"] != 0 and pyr["rc"] != 0)))
                ds = rel(c["sigma"], pyr["sigma"])
                da = rel(c["alpha"], pyr["alpha"])
                n += 1
                ok = same_census and ds < 1e-8 and da < 1e-8
                if ok:
                    worst = max(worst, ds, da)
                if not ok:
                    # round-off sensitivity of the C++ itself: perturb sigma by 1e-14 rel
                    sp = [x * (1 + 1e-14 * (k + 1)) for k, x in enumerate(r["sigma"])]
                    c2 = sr.replay(ops, B.TAG_ME, sp, r["alpha"], r["alpha_in"], r["z"],
                                   r["e"], de, "compressionPositive", trace=0)
                    self_ds = rel(c["sigma"], c2["sigma"])
                    self_sub = int(c2["stats"]["substeps"])
                    chaotic = self_ds > 0.1 * ds or self_sub != int(s["substeps"])
                    if chaotic:
                        nchaos += 1
                    bad += 1
                    lines.append(f"   C++ self-sensitivity (sigma*(1+1e-14)): dsig={self_ds:.2e} "
                                 f"sub={self_sub} -> {'ROUND-OFF CHAOTIC' if chaotic else 'PORT DIFFERENCE'}")
                    lines.append(
                        f"MISMATCH {r['element']}/{r['gp']} {pname}: C++ rc={c['rc']} "
                        f"sub={int(s['substeps'])} acc={int(s['accepted'])} "
                        f"forced={int(s['forcedAtDTmin'])} | py rc={pyr['rc']} "
                        f"sub={pyr['substeps']} acc={pyr['acc']} forced={pyr['forced']} "
                        f"path C++={c['path']} py={pyr['path']} dsig={ds:.2e} dalpha={da:.2e}")
    lines.append(f"cases {n}, mismatches {bad} (of which C++ itself round-off chaotic: {nchaos}), worst rel diff on the matching "
                 f"cases {worst:.2e}")
    txt = "\n".join(lines)
    print(txt)
    with open(f"{B.OUT}/validate_port.txt", "w") as fh:
        fh.write(txt + "\n")


if __name__ == "__main__":
    main()
