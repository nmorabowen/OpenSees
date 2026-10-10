import json, random, sys
exec(open(r"C:\Users\nmora\AppData\Local\Temp\claude\C--Users-nmora-Github-OpenSees-Compile-OpenSees--claude-worktrees-tims-implementation-review-3733c6\e034e494-142a-4e05-862e-9d267344b7ad\scratchpad\r1\h.py").read())

base = (2, 2, 1, 1e-7, 1e-7, "-Presidual", 0.0, "-Pmin", 0.0101, "-flipAlphaIn", "init")
MATS = {
    1: base,
    2: base + ("-cppmStart", "explicit", "-cppmHalvings", 0, "-cppmOnFail", "refuse"),
    3: base + ("-cppmStart", "explicit"),
    4: base + ("-cppmLineSearch", "on", "-cppmHalvings", 0, "-cppmOnFail", "refuse"),
    5: base + ("-cppmOnFail", "refuse", "-cppmHalvings", 0),
}
ops.wipe()
ops.model("basic", "-ndm", 3, "-ndf", 3)
for t, o in MATS.items():
    ops.nDMaterial("LadrunoSANISAND", t, *CAMPAIGN, *o)

rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
Mc = 1.3309
out = []
for it in range(int(sys.argv[2]) if len(sys.argv) > 2 else 300):
    p0 = rng.choice([2.0, 10.0, 50.0, 150.0])
    e = rng.choice([0.62, 0.72, 0.85])
    # deviatoric alpha of random direction, magnitude fraction of Mc
    d = [rng.gauss(0, 1) for _ in range(3)]
    m = sum(d) / 3; d = [x - m for x in d]
    sh = [rng.gauss(0, 0.3) for _ in range(3)]
    nrm = math.sqrt(sum(x * x for x in d) + 2 * sum(x * x for x in sh))
    frac = rng.choice([0.2, 0.6, 0.9, 1.1])
    a = [frac * math.sqrt(2 / 3) * Mc * x / nrm for x in d + sh]
    sig = [p0 * (1 + a[i]) if i < 3 else p0 * a[i] for i in range(6)]  # on the yield cone axis
    ain = list(a) if rng.random() < 0.5 else [0.0] * 6
    z = [0.0] * 6
    mag = rng.choice([1e-4, 5e-4, 2e-3, 1e-2])
    dd = [rng.gauss(0, 1) for _ in range(6)]
    vol = rng.choice([0.0, 0.3, -0.3])
    tr = sum(dd[:3]) / 3
    dd = [dd[i] - tr + vol if i < 3 else dd[i] for i in range(6)]
    nd = math.sqrt(sum(x * x for x in dd))
    de = [mag * x / nd for x in dd]
    row = dict(sig=sig, alpha=a, ain=ain, z=z, e=e, de=de, res={})
    for t in MATS:
        r = sr.replay(ops, t, sig, a, ain, z, e, de, "compressionPositive")
        row["res"][t] = dict(rc=r["rc"], sigma=r["sigma"], alpha=r["alpha"], z=r["z"], p=r["p"],
                             f_after=r["f_after"], stats=r["stats"], path=r["path"])
    out.append(row)
json.dump(out, open(R + r"\p2_out_%s.json" % (sys.argv[1] if len(sys.argv) > 1 else "1"), "w"))
# summary
def dist(x, y):
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(x, y)))
n_guess = n_div = 0
for i, row in enumerate(out):
    r1, r2, r3 = row["res"][1], row["res"][2], row["res"][3]
    if r2["stats"]["cppmGuessOk"] > 0:
        n_guess += 1
        s = max(1.0, max(abs(x) for x in r1["sigma"]))
        d12 = dist(r1["sigma"], r2["sigma"]) / s
        if d12 > 0.05:
            n_div += 1
            print(i, "guess-ok result differs from default ladder by", round(d12, 3),
                  "p0", row["sig"][0], "mag", max(map(abs, row["de"])),
                  "r1 path", r1["path"], dict((k, r1["stats"][k]) for k in ("cppmNewtonFail", "cppmHalvings", "cppmExplicitFail")),
                  "p1", round(r1["p"], 3), "p2", round(r2["p"], 3), "f2", r2["f_after"])
print("rows", len(out), "guessOk rows", n_guess, "diverging >5%", n_div)
for t in MATS:
    rcs = [row["res"][t]["rc"] for row in out]
    print(t, "rc counts", {x: rcs.count(x) for x in set(rcs)},
          "refusals", sum(row["res"][t]["stats"]["cppmRefusals"] for row in out),
          "negp", sum(1 for row in out if row["res"][t]["p"] < 0))
