"""ADR-94 wp/94c byte-identity comparison: pre-94c vs post-94c committed histories."""
import json
import sys

a = json.load(open(sys.argv[1]))
b = json.load(open(sys.argv[2]))
keys = sorted(set(a) | set(b))
for k in keys:
    va, vb = a.get(k, {}), b.get(k, {})
    if "stress" not in va or "stress" not in vb:
        print("%-52s MISSING" % k)
        continue
    same_codes = va["codes"] == vb["codes"]
    n = min(len(va["stress"]), len(vb["stress"]))
    worst = 0.0
    scale = 0.0
    nbit = 0
    for i in range(n):
        for x, y in zip(va["stress"][i], vb["stress"][i]):
            if x != y:
                nbit += 1
            worst = max(worst, abs(x - y))
            scale = max(scale, abs(x), abs(y))
    rel = worst / scale if scale else 0.0
    tag = "IDENTICAL" if (nbit == 0 and same_codes and
                          len(va["stress"]) == len(vb["stress"])) else "DIFFERS"
    print("%-52s %-10s codes %s->%s  rows %d->%d  max|d|=%.3e rel=%.3e" % (
        k, tag, len(va["codes"]), len(vb["codes"]),
        len(va["stress"]), len(vb["stress"]), worst, rel))
