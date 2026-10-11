"""F3-d diagnostic: replay the exact openseespy call stream apeGmsh emitted
for three legs (MD, LS default/init, LS -flipAlphaIn vanilla) and compare."""
import json, sys
import opensees as ops
print("build:", ops.ladrunoBuild())
L = json.load(open(sys.argv[1]))
legs, cur = [], None
for name, args, kw, rec in L:
    if name == "wipe":
        cur = {"q": [], "rec": []}; legs.append(cur)
    r = getattr(ops, name)(*args, **kw)
    if name in ("nodeDisp", "eleResponse"):
        cur["q"].append(r if isinstance(r, list) else [r]); cur["rec"].append(rec if isinstance(rec, list) else [rec])
names = ["MD", "LS_init_default", "LS_vanilla"]
for n, leg in zip(names, legs):
    same = leg["q"] == leg["rec"]
    print(f"{n}: {len(leg['q'])} queries; identical to Windows-recorded values: {same}")
md = legs[0]["q"]
for n, leg in zip(names[1:], legs[1:]):
    diff = [i for i, (a, b) in enumerate(zip(md, leg["q"])) if a != b]
    print(f"{n} vs MD: bit-identical over all {len(md)} queries: {not diff}; first differing query idx: {diff[0] if diff else None}")
# final settlement = mean of last 4 nodeDisp of each leg
def uz(leg):
    nd = [c for c in leg["q"] if len(c) == 1]
    return sum(x[0] for x in nd[-4:]) / 4
for n, leg in zip(names, legs):
    print(f"final uz {n}: {uz(leg)!r}")
# Per-GP alpha_in at the first plastic step vs alpha at the flip (snapshot
# layout from trace.py: 4 nodeDisp + GP1 alpha/alpha_in/stress + 8x(alpha, alpha_in)).
for n, leg in zip(names, legs):
    s0, s1 = leg["q"][0:23], leg["q"][23:46]
    same = [g + 1 for g in range(8) if s1[8 + 2 * g] == s0[7 + 2 * g]]
    print(f"{n}: GPs with alpha_in(step1) == alpha(flip): {same}")
