import csv,json,os,glob,sys
pat=sys.argv[1:] or ["K5p_*"]
for d in sorted(sum((glob.glob(p) for p in pat),[])):
    f=d+"/steps.csv"
    if not os.path.exists(f): print(d,"no steps"); continue
    R=list(csv.DictReader(open(f)))
    if not R: print(d,"empty"); continue
    S=[float(r["s_over_B"]) for r in R]; Q=[float(r["q_kPa"]) for r in R]
    at=lambda x: Q[min(range(len(S)),key=lambda i:abs(S[i]-x))] if S[-1]>=x-1e-3 else float("nan")
    im=max(range(len(Q)),key=Q.__getitem__)
    s=json.load(open(d+"/summary.json")) if os.path.exists(d+"/summary.json") else {}
    pl=s.get("plateau") or {}
    print("%-30s end %.4f qmax %7.1f @ %.4f | .02 .05 .10 .15 .20 .25 .30 .35 .40: %s | mode %s fails %s wall %s plateau %s k_norm %s refused %s"%(d,S[-1],Q[im],S[im]," ".join("%.0f"%at(x) for x in (.02,.05,.1,.15,.2,.25,.3,.35,.4)),s.get("mode"),s.get("failed_attempts"),int(s.get("total_wall_s") or 0),pl.get("verdict"),pl.get("k_norm"),(s.get("norsand") or {}).get("refused_trials")))
