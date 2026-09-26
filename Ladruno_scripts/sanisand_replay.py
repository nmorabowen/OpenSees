"""LadrunoSANISAND material-point state replay (WP-127 / TIMs F21).

Puts ONE LadrunoSANISAND material point into a given state (sigma, alpha,
alpha_in, fabric z, void ratio e), applies ONE prescribed strain increment,
and returns what the integrator did: return code, substep census, the
per-substep ModifiedEuler trace (T, dT, err, outcome, at-dT_min flag) and the
returned state.  A thin wrapper over the `ladrunoSANISANDReplay` command
(SRC/material/nD/LadrunoSANISAND.cpp); see the SANISAND section of the
opensees-ladruno guide for the command itself.

The command works on a PRIVATE getCopy() of the nDMaterial prototype `tag`:
no element, no domain, no analysis, and nothing it does survives the call
(the process-wide ManzariDafalias stage flag and `ops_Dt` are set for the
call and restored).

SIGN CONVENTION -- ALWAYS EXPLICIT.  `convention` is required:
  "compressionPositive"  sigma is the model's INTERNAL stress (mSigma,
                         compression positive), dStrain compression positive.
                         The TIMs ring-point CSVs are in THIS convention even
                         though their README says otherwise: on all 80 rows
                         p_kPa == +tr(sigma)/3 (WP-127 plan, finding A).
  "tensionPositive"      sigma as `eleResponse ... stress` returns it and
                         dStrain as an element strain (OpenSees convention).
alpha, alpha_in and z are ALWAYS the raw internal tensors (as the `alpha` /
`state` responses report them) -- they are stress RATIOS and carry no sign
flip.  Shear strains are ENGINEERING (gamma) in both conventions.

usage (CLI): replay both attached CSVs with the documented probes

    python -S <bootstrap> Ladruno_scripts/sanisand_replay.py [--delta 1e-5] [--csv PATH ...]

(the fork's test bootstrap: CPython 3.12 with -S and this worktree's
dist/bin on sys.path -- see Ladruno_internal/BUILD_GOTCHAS.md).
"""
import csv
import math
import os

# --- the TIMs 2d-model campaign material (README of the attachments) --------
# G0, nu, e_init, Mc, c, lambda_c, e0, ksi, P_atm, m, h0, ch, nb, A0, nd,
# z_max, cz, Den.  nu = 0.312885 is the Jaky K0 substitution.
CAMPAIGN_PARAMS = [264.32, 0.312885, 0.6944, 1.3309, 0.71, 0.027, 0.83, 0.45,
                   101.0, 0.005, 1.3, 0.968, 3.5, 0.05, 5.75, 12.5, 1100.0,
                   2.0]
# IntScheme 1, TanType 0, JacoType 1, TolF, TolR, then the fork flags.
CAMPAIGN_OPTS = (1, 0, 1, 1.0e-7, 1.0e-7,
                 "-Presidual", 0.0, "-Pmin", 0.0101, "-maxSubsteps", 20000,
                 "-flipAlphaIn", "init")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
ATTACH = os.path.join(ROOT, "Ladruno_implementation",
                      "_tims_2d_model_requests_2026-09-25")
RING_CSVS = [os.path.join(ATTACH, "ring_points_b8.csv"),
             os.path.join(ATTACH, "ring_points_b16.csv")]

# --- output layout of ladrunoSANISANDReplay (format version 1) ------------
STAT_NAMES = [
    # cumulative since revertToStart (on the private copy: this call only)
    "updates", "meCalls", "substeps", "accepted", "rejectedErr",
    "forcedAtDTmin", "forcedClampMc", "rejectedLowP", "abandonedLowP",
    "capHits", "entryPminClamps", "pnResets", "maxSubstepsOneUpdate",
    # the last update that entered ModifiedEuler
    "lastSubsteps", "lastForcedAtDTmin", "lastAbandonedLowP", "lastCapHit",
]
STATE_NAMES = (["sigma"] * 6 + ["alpha"] * 6 + ["alpha_in"] * 6 + ["z"] * 6
               + ["e", "p", "q", "fBefore", "fAfter", "path", "elasticRatio",
                  "trAlpha0", "trAlphaIn0", "trZ0"])
TRACE_FIELDS = ["T", "dT", "err", "code", "atDTmin"]
TRACE_CODES = {
    0: "accept",
    1: "rejectErr",          # err > TolE, dT > dT_min: retried at q*dT
    2: "forcedAtDTmin",      # err > TolE at dT == dT_min: ACCEPTED (finding C)
    3: "forcedAtDTminClampMc",  # ... and the radial eta -> Mc clamp fired
    4: "rejectLowPStage1",   # p < p_r after the first stage: dT cut
    5: "rejectLowPEnd",      # p < p_r at the end of the substep: dT cut
    6: "abandonLowP",        # p < p_r at dT == dT_min: RETURNS with T < 1
    7: "capHit",             # -maxSubsteps fired: update refused
}
PATH_CODES = {
    -1: "notExplicit", 0: "elastic", 1: "startOutsideYield",
    2: "elasticToPlastic", 3: "plastic", 4: "unloadThenPlastic",
    5: "pnBelowPresidualReset",
}


def define_campaign_material(ops, tag=1):
    """`nDMaterial LadrunoSANISAND` with the TIMs campaign set."""
    ops.nDMaterial("LadrunoSANISAND", tag, *CAMPAIGN_PARAMS, *CAMPAIGN_OPTS)


def _vec6(v, name):
    v = [float(x) for x in v]
    if len(v) != 6:
        raise ValueError(f"{name}: need 6 components, got {len(v)}")
    return v


def replay(ops, tag, sigma, alpha, alpha_in, z, e, dstrain, convention,
           mat_type="3D", trace=10000, dt=1.0, primed=True,
           prev_incr_norm=0.0):
    """Run one replay; returns a dict (see module docstring for semantics)."""
    if convention not in ("compressionPositive", "tensionPositive"):
        raise ValueError("convention must be 'compressionPositive' or "
                         "'tensionPositive' (no default, on purpose)")
    args = [tag, "-convention", convention,
            "-sigma", *_vec6(sigma, "sigma"),
            "-alpha", *_vec6(alpha, "alpha"),
            "-alphaIn", *_vec6(alpha_in, "alpha_in"),
            "-fabric", *_vec6(z, "z"),
            "-voidRatio", float(e),
            "-dStrain", *_vec6(dstrain, "dstrain"),
            "-type", mat_type, "-trace", int(trace), "-dt", float(dt),
            "-primed", 1 if primed else 0, "-prevIncrNorm", float(prev_incr_norm)]
    out = ops.ladrunoSANISANDReplay(*args)
    if out is None or len(out) < 6:
        raise RuntimeError(f"ladrunoSANISANDReplay failed (returned {out!r}); "
                           "see stderr")
    out = list(out)
    version, rc, n_stats, n_rec, width, dropped = out[:6]
    if int(version) != 1:
        raise RuntimeError(f"unknown replay output format {version}")
    n_stats, n_rec, width = int(n_stats), int(n_rec), int(width)
    i = 6
    stats = dict(zip(STAT_NAMES, out[i:i + n_stats]))
    i += n_stats
    st = out[i:i + 34]
    i += 34
    recs = []
    for k in range(n_rec):
        r = out[i + k * width:i + (k + 1) * width]
        recs.append(dict(T=r[0], dT=r[1], err=r[2], code=int(r[3]),
                         outcome=TRACE_CODES.get(int(r[3]), "?"),
                         atDTmin=bool(r[4])))
    return dict(rc=int(rc), stats=stats, trace=recs,
                trace_dropped=int(dropped),
                sigma=st[0:6], alpha=st[6:12], alpha_in=st[12:18], z=st[18:24],
                e=st[24], p=st[25], q=st[26], f_before=st[27], f_after=st[28],
                path=PATH_CODES.get(int(st[29]), str(st[29])),
                elastic_ratio=st[30], tr_alpha0=st[31], tr_alpha_in0=st[32],
                tr_z0=st[33])


def read_ring_csv(path):
    """The TIMs ring-point CSV -> list of dicts (numbers as floats).

    The sigma columns are INTERNAL, compression-positive (finding A) --
    replay them with convention="compressionPositive"."""
    rows = []
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh):
            d = {k: float(v) for k, v in r.items()}
            d["element"] = int(d["element"])
            d["gp"] = int(d["gp"])
            for key in ("sigma", "alpha", "alpha_in", "z"):
                d[key] = [d[f"{key}_{i}"] for i in range(6)]
            rows.append(d)
    return rows


def check_sign_convention(row):
    """Finding A, row by row: returns (p_csv, +tr/3, -tr/3)."""
    tr3 = sum(row["sigma"][:3]) / 3.0
    return row["p_kPa"], tr3, -tr3


def probes(delta):
    """The documented WP-127 probes, COMPRESSION-POSITIVE (internal) strain,
    engineering shear, plane-strain-admissible (zz = yz = zx = 0).

    The act's local strain direction at failure is unknown, so two canonical
    directions are used: an in-plane isotropic compression and a simple
    shear.  Both at magnitude `delta`."""
    return {
        "isoComp": [delta, delta, 0.0, 0.0, 0.0, 0.0],
        "shear": [0.0, 0.0, 0.0, delta, 0.0, 0.0],
    }


def replay_row(ops, tag, row, dstrain, **kw):
    return replay(ops, tag, row["sigma"], row["alpha"], row["alpha_in"],
                  row["z"], row["e"], dstrain, "compressionPositive", **kw)


def _summary(res):
    s = res["stats"]
    codes = {}
    for r in res["trace"]:
        codes[r["outcome"]] = codes.get(r["outcome"], 0) + 1
    errs = [r["err"] for r in res["trace"] if math.isfinite(r["err"])]
    dts = [r["dT"] for r in res["trace"]]
    return (f"rc={res['rc']:>3} path={res['path']:<18} substeps={int(s['substeps']):>6} "
            f"acc={int(s['accepted']):>6} rej={int(s['rejectedErr']):>6} "
            f"forced@dTmin={int(s['forcedAtDTmin']):>5} abandonLowP={int(s['abandonedLowP'])} "
            f"cap={int(s['capHits'])} minDT={min(dts) if dts else float('nan'):.1e} "
            f"maxErr={max(errs) if errs else float('nan'):.2e} p_out={res['p']:.4g}")


def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--csv", nargs="*", default=RING_CSVS)
    ap.add_argument("--delta", type=float, nargs="*", default=[1.0e-5])
    ap.add_argument("--rows", type=int, default=0,
                    help="replay only the first N rows of each CSV (0 = all)")
    args = ap.parse_args(argv)
    try:
        import opensees as ops
    except ModuleNotFoundError:
        import openseespy.opensees as ops
    ops.wipe()
    define_campaign_material(ops, 1)
    for path in args.csv:
        rows = read_ring_csv(path)
        if args.rows:
            rows = rows[:args.rows]
        print(f"== {os.path.basename(path)}: {len(rows)} rows")
        bad = [r for r in rows
               if abs(check_sign_convention(r)[0] - check_sign_convention(r)[1])
               > 1e-4 * max(1.0, abs(r["p_kPa"]))]
        print(f"   finding A: p_kPa == +tr(sigma)/3 on {len(rows) - len(bad)}/{len(rows)} rows")
        for delta in args.delta:
            for name, de in probes(delta).items():
                for r in rows:
                    res = replay_row(ops, 1, r, de)
                    print(f"   el {r['element']:>5} gp {r['gp']} eta/Mb {r['eta_over_Mb_compression']:6.3f} "
                          f"{name}@{delta:.0e}: {_summary(res)}")


if __name__ == "__main__":
    main()
