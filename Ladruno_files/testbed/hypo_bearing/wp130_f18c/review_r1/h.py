exec(open(r"C:\Users\nmora\AppData\Local\Temp\claude\C--Users-nmora-Github-OpenSees-Compile-OpenSees--claude-worktrees-tims-implementation-review-3733c6\e034e494-142a-4e05-862e-9d267344b7ad\scratchpad\r1\boot.py").read())
import math
import sanisand_replay as sr

PARAMS = [264.32, 0.3129, 0.6944, 1.33090, 0.71, 0.027, 0.83, 0.45, 101.0, 0.005,
          1.3, 0.968, 3.5, 0.05, 5.75, 12.5, 1100.0, 2.0]
CAMPAIGN = list(PARAMS); CAMPAIGN[1] = 0.312885
XY = [(0., 0.), (1., 0.), (1., 1.), (0., 1.)]


def analysis(tol=1e-13):
    ops.constraints("Transformation"); ops.numberer("Plain"); ops.system("FullGeneral")
    ops.test("NormDispIncr", tol, 25, 0); ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0); ops.analysis("Static")


def series(lat_path, ax_path):
    lat = list(lat_path) + [lat_path[-1]]
    ax = list(ax_path) + [ax_path[-1]]
    ops.timeSeries("Path", 1, "-dt", 1.0, "-values", *lat)
    ops.timeSeries("Path", 2, "-dt", 1.0, "-values", *ax)


def paths(n_conf, e_conf, incs):
    lat = [i / n_conf for i in range(n_conf + 1)]
    ax = list(lat)
    for d_lat, d_ax in incs:
        lat.append(lat[-1] + d_lat / e_conf)
        ax.append(ax[-1] + d_ax / e_conf)
    return lat, ax


def build_3d(matcmd, params, opts, n_conf, e_conf, incs, ele="stdBrick"):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for k in range(2):
        for j, (x, y) in enumerate(XY):
            ops.node(4 * k + j + 1, x, y, float(k))
    ops.nDMaterial(matcmd, 1, *params, *opts)
    ops.element(ele, 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
    for k in range(2):
        for j, (x, y) in enumerate(XY):
            ops.fix(4 * k + j + 1, 1 if x == 0. else 0, 1 if y == 0. else 0, 1 if k == 0 else 0)
    series(*paths(n_conf, e_conf, incs))
    ops.pattern("Plain", 1, 1)
    for k in range(2):
        for j, (x, y) in enumerate(XY):
            n = 4 * k + j + 1
            if x == 1.: ops.sp(n, 1, -e_conf)
            if y == 1.: ops.sp(n, 2, -e_conf)
    ops.pattern("Plain", 2, 2)
    for j in range(4):
        ops.sp(4 + j + 1, 3, -e_conf)
    analysis()


def build_ps(params, opts, n_conf, e_conf, incs, ele="quad"):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for j, (x, y) in enumerate(XY):
        ops.node(j + 1, x, y)
    ops.nDMaterial("LadrunoSANISAND", 1, *params, *opts)
    if ele == "quad":
        ops.element("quad", 1, 1, 2, 3, 4, 1.0, "PlaneStrain", 1)
    else:
        ops.element("SSPquad", 1, 1, 2, 3, 4, 1, "PlaneStrain", 1.0)
    for j, (x, y) in enumerate(XY):
        ops.fix(j + 1, 1 if x == 0. else 0, 1 if y == 0. else 0)
    series(*paths(n_conf, e_conf, incs))
    ops.pattern("Plain", 1, 1)
    for j, (x, y) in enumerate(XY):
        if x == 1.: ops.sp(j + 1, 1, -e_conf)
    ops.pattern("Plain", 2, 2)
    for j, (x, y) in enumerate(XY):
        if y == 1.: ops.sp(j + 1, 2, -e_conf)
    analysis()


def grab(ele=1, gp=1):
    g = lambda name: list(ops.eleResponse(ele, "material", gp, name))
    return dict(sig=g("stress"), eps=g("strain"), alpha=g("alpha"), ain=g("alpha_in"),
                z=g("fabric"), e=g("state")[24], tan=g("tangent"), st=g("substepStats"))


def fd_from(k, k1, prev, h=1e-8, dirs=None):
    """FD of the return map from committed state k over k1's increment."""
    de = [a - b for a, b in zip(k1["eps"], k["eps"])]
    dn = math.sqrt(sum((a - b) ** 2 for a, b in zip(k["eps"][:3], prev["eps"][:3]))
                   + 0.5 * sum((a - b) ** 2 for a, b in zip(k["eps"][3:], prev["eps"][3:])))
    run = lambda d: sr.replay(ops, 1, k["sig"], k["alpha"], k["ain"], k["z"], k["e"], d,
                              "tensionPositive", prev_incr_norm=dn)
    base = run(de)
    rep = max(abs(a - b) for a, b in zip(base["sigma"], k1["sig"])) / max(map(abs, k1["sig"]))
    if dirs is None:
        dirs = [[1.0 if i == j else 0.0 for i in range(6)] for j in range(6)]
    cols = []
    for v in dirs:
        dp = [a + h * b for a, b in zip(de, v)]
        dm = [a - h * b for a, b in zip(de, v)]
        rp, rm = run(dp), run(dm)
        cols.append([(rp["sigma"][i] - rm["sigma"][i]) / (2 * h) for i in range(6)])
    T = [[k1["tan"][6 * i + j] for j in range(6)] for i in range(6)]
    Tcols = [[sum(T[i][j] * v[j] for j in range(6)) for i in range(6)] for v in dirs]
    return cols, Tcols, rep, base


def relerr(A, B):
    n = math.sqrt(sum(x * x for c in B for x in c))
    return math.sqrt(sum((a - b) ** 2 for ca, cb in zip(A, B) for a, b in zip(ca, cb))) / n
