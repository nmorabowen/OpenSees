import shutil, sys
p = sys.argv[1]
shutil.copy(p, p + ".bak_pre_ladders")
s = open(p).read()

def rep(old, new):
    global s
    assert s.count(old) == 1, old
    s = s.replace(old, new)

rep("""def build(mat, scheme, extra, maxsub, dp_phi, dp_g, deterministic, tolr=1.0e-7,
          tantype=0):""",
"""def build(mat, scheme, extra, maxsub, dp_phi, dp_g, deterministic, tolr=1.0e-7,
          tantype=0, presidual=0.0, einit=None):""")
rep("""        args = ["LadrunoSANISAND", MAT, *SAN, int(scheme)""",
"""        san = list(SAN)
        if einit is not None:
            san[2] = float(einit)
        args = ["LadrunoSANISAND", MAT, *san, int(scheme)""")
rep(""""-maxSubsteps", int(maxsub), "-Presidual", 0.0,""",
    """"-maxSubsteps", int(maxsub), "-Presidual", float(presidual),""")
rep("""    ap.add_argument("--out", required=True)""",
"""    ap.add_argument("--presidual", type=float, default=0.0,
                    help="LadrunoSANISAND -Presidual (kPa); 0 = the campaign's (cohesionless)")
    ap.add_argument("--einit", type=float, default=E_INIT,
                    help="initial void ratio e_init (positional 3rd SANISAND param); "
                         "default = the campaign's 0.6944")
    ap.add_argument("--out", required=True)""")
rep("""                 bool(args.deterministic), args.tolr, args.tantype)""",
"""                 bool(args.deterministic), args.tolr, args.tantype,
                 presidual=args.presidual, einit=args.einit)
    log(f"ladder knobs: Presidual = {args.presidual} kPa, e_init = {args.einit} "
        f"(campaign 0.6944). Initial stress is e_init-INDEPENDENT: gamma' = {GAMMA} "
        f"kN/m3 fixed, K0 = nu*/(1-nu*) = {K0:.4f} fixed by nu* (stage-0 elastic); "
        f"e_init changes only the elastic moduli magnitude G(e,p) and psi")""")
rep("""                   dead_end_first_step={int(k): int(v) for k, v in DEADEND.items()})
""",
"""                   dead_end_first_step={int(k): int(v) for k, v in DEADEND.items()})
    if san and CENSUS is CENSUS_SAS:
        try:
            Fz = read_field(deck, want_stats=True, want_f=False)
            tot = Fz["stats"].sum(axis=0)
            summary["sas_totals"] = {CENSUS["names"][i]: float(tot[i])
                                     for i in range(min(len(tot), len(CENSUS["names"])))}
            summary["refusals_by_code"] = {SAS_REF_CODES[i]: int(tot[i])
                                           for i in SAS_REF_CODES if i < len(tot)}
            log(f"cumulative sasStats refusals by code: {summary['refusals_by_code']}")
        except Exception as exc:
            log(f"NOTE sas totals failed: {exc}")
""")
open(p, "w").write(s)
print("patched")
