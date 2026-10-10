import shutil, sys
p = sys.argv[1]
shutil.copy(p, p + ".bak_pre_ablation")
s = open(p).read()

def rep(old, new):
    global s
    assert s.count(old) == 1, old
    s = s.replace(old, new)

rep("""          tantype=0, presidual=0.0, einit=None):""",
"""          tantype=0, presidual=0.0, einit=None, over=None):""")
rep("""        if einit is not None:
            san[2] = float(einit)
""",
"""        if einit is not None:
            san[2] = float(einit)
        # ablation overrides: SAN index -> value (h0 10, nb 12, A0 13, nd 14, zmax 15)
        for i, v in (over or {}).items():
            san[i] = float(v)
""")
rep("""    ap.add_argument("--out", required=True)""",
"""    ap.add_argument("--zmax", type=float, default=SAN[15], help="fabric z_max (campaign 12.5)")
    ap.add_argument("--nb", type=float, default=SAN[12], help="n^b (campaign 3.5)")
    ap.add_argument("--nd", type=float, default=SAN[14], help="n^d (campaign 5.75)")
    ap.add_argument("--A0", type=float, default=SAN[13], help="A0 dilatancy (campaign 0.05)")
    ap.add_argument("--h0", type=float, default=SAN[10], help="h0 hardening (campaign 1.3)")
    ap.add_argument("--out", required=True)""")
rep("""                 presidual=args.presidual, einit=args.einit)
""",
"""                 presidual=args.presidual, einit=args.einit,
                 over={15: args.zmax, 12: args.nb, 14: args.nd, 13: args.A0, 10: args.h0})
    global NB, ND
    NB, ND = args.nb, args.nd          # derived() diagnostics follow the ablation
    log(f"ablation knobs: zmax = {args.zmax}, nb = {args.nb}, nd = {args.nd}, "
        f"A0 = {args.A0}, h0 = {args.h0} (campaign 12.5 / 3.5 / 5.75 / 0.05 / 1.3)")
""")
open(p, "w").write(s)
print("patched")
