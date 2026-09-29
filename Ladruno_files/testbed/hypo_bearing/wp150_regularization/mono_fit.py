"""WP-150: a physically bounded MONOTONIC parameter set for the campaign sand (owner GO 2026-09-28).

FIXED (TIMs' anchors / unsourced): G0 264.32, nu 0.312885, e_init 0.6944, Mc 1.3309, c 0.71, lambda_c 0.027, e0 0.83,
xi 0.45, P_atm 101, m 0.005, ch 0.968, zmax 12.5, cz 1100.  FREE: nb, A0, nd, h0 (the deck's --nb --A0 --nd --h0).
TARGETS (Bolton 1986) at the working D_r, with I_R = D_r (10 - ln p'_peak) - 1, clipped to [0, 4]:
  TX:  phi_peak - phi_cs,TX = 3 I_R,  (-de_v/de_1)_max = 0.3 I_R;   PS:  phi_peak - phi_cs,PS = 5 I_R;
  peak axial strain within 1-5 % (TX and PS).
phi_cs,TX = asin(3 Mc/(6+Mc)) = 33.0 deg; phi_cs,PS = 39.5 deg ESTIMATED (TX + 6.5 deg, the PS-TX critical-state offset
the DM04 Toyoura run shows; memo section 10 footnote 1).
Integrated EXACTLY (WP-134 oracle, uw_model options = the C++'s equations), drained PS and TX from isotropic
p0 = 10, 50, 150, 500 kPa to 25 % axial strain.
    python mono_fit.py fit Dr=0.47 [start=nb,A0,nd,h0] [fix=nd:3.5]  -> Nelder-Mead in log-parameters (free ones), parallel
    python mono_fit.py eval Dr=0.47 nb=.. A0=.. nd=.. h0=..  -> the fit report for one set
"""
import json, math, os, sys
from concurrent.futures import ProcessPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "Ladruno_scripts"))

P0S = (10.0, 50.0, 150.0, 500.0)
PHI_CS_TX = math.degrees(math.asin(3 * 1.3309 / (6 + 1.3309)))
PHI_CS_PS = PHI_CS_TX + 6.5
E0 = 0.6944


def _run(args):
    kind, p0, nb, A0, nd, h0 = args
    import numpy as np
    from dataclasses import replace
    from sanisand_reference import CAMPAIGN
    from sanisand_reference.driver import isotropic_state
    from sanisand_reference.integrator import Control, integrate, path_table, _unpack
    from sanisand_reference.ring import ring_variants
    P = replace(CAMPAIGN, e_init=E0, nb=nb, A0=A0, nd=nd, h0=h0)
    O = ring_variants()["uw_model"]
    mask = (True, False, True, True, True, True) if kind == "PS" else (True, False, False, True, True, True)
    res = integrate(isotropic_state(p0, E0), Control(mask, (0.25, 0, 0, 0, 0, 0)), P, O, rtol=1e-7)
    tab = path_table(res, P, O)
    ea = np.array([t["eps"][0] for t in tab]); eps = np.array([t["eps"][:3] for t in tab]); ev = eps.sum(1)
    sig = np.array([t["sigma"] for t in tab])
    s1, s3 = sig[:, 0], np.minimum(sig[:, 1], sig[:, 2])
    phi = np.degrees(np.arcsin(np.clip((s1 - s3) / (s1 + s3), -1, 1)))
    ip = int(np.argmax(phi))
    de1, dev = np.gradient(ea), np.gradient(ev)
    with np.errstate(invalid="ignore", divide="ignore"):
        if kind == "PS":
            dil = np.degrees(np.arcsin(np.clip(-dev / (de1 - np.gradient(eps[:, 1])), -1, 1)))
        else:
            dil = -dev / de1
    # fabric activity under monotonic loading: max <z:n> along the path
    zn = 0.0
    Y = res.path["y"]
    for k in range(0, len(Y), max(1, len(Y) // 60)):
        sg, al, z, e, _ = _unpack(Y[k])
        p = np.trace(sg) / 3.0
        nn = (sg - p * np.eye(3)) - p * al
        nrm = np.sqrt(np.sum(nn * nn))
        if nrm > 0:
            zn = max(zn, float(np.sum(z * nn) / nrm))
    return dict(kind=kind, p0=p0, status=res.status, phi_peak=float(phi[ip]), eps_peak=float(ea[ip]),
                p_peak=float(tab[ip]["p"]), phi_end=float(phi[-1]), dil_max=float(np.nanmax(dil)),
                psi0=float(tab[0]["psi"]), zn_max=zn,
                path=dict(eps_a=ea[::4].tolist(), q=[t["q"] for t in tab][::4], ev=ev[::4].tolist()))


def IR(Dr, p):
    return min(max(Dr * (10.0 - math.log(p)) - 1.0, 0.0), 4.0)


def evaluate(pars, Dr, pool):
    nb, A0, nd, h0 = pars
    jobs = [(k, p0, nb, A0, nd, h0) for k in ("TX", "PS") for p0 in P0S]
    out = list(pool.map(_run, jobs))
    err, rows = 0.0, []
    for r in out:
        ir = IR(Dr, r["p_peak"])
        if r["kind"] == "TX":
            t_dphi, t_dil = 3 * ir, 0.3 * ir
            dphi = r["phi_peak"] - PHI_CS_TX
            e = ((dphi - t_dphi) / 1.0) ** 2 + ((r["dil_max"] - t_dil) / 0.05) ** 2
        else:
            t_dphi, t_dil = 5 * ir, float("nan")
            dphi = r["phi_peak"] - PHI_CS_PS
            e = ((dphi - t_dphi) / 1.0) ** 2
        ep = r["eps_peak"]
        e += (max(0.0, ep - 0.05) / 0.005) ** 2 + (max(0.0, 0.01 - ep) / 0.002) ** 2
        if r["status"] != "ok":
            e += 1e4
        err += e
        rows.append(dict(r, IR=ir, target_dphi=t_dphi, dphi=dphi, target_dil=t_dil))
    return err, rows


def report(rows, pars, Dr, label):
    L = [f"### {label}: nb {pars[0]:.3f}, A0 {pars[1]:.3f}, nd {pars[2]:.3f}, h0 {pars[3]:.3f} (D_r {Dr})", "",
         "| test | p0 | ψ0 | p′ at peak | I_R | φ′_peak ° | Δφ model / Bolton ° | ε at peak | max dilatancy model / Bolton | max ⟨z:n⟩ |",
         "|---|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        dil = (f"{r['dil_max']:.3f} / {r['target_dil']:.3f}" if r["kind"] == "TX" else f"ψ_max {r['dil_max']:.1f}° / —")
        L.append(f"| {r['kind']} | {r['p0']:.0f} | {r['psi0']:+.3f} | {r['p_peak']:.0f} | {r['IR']:.2f} | {r['phi_peak']:.1f} | "
                 f"{r['dphi']:.1f} / {r['target_dphi']:.1f} | {100*r['eps_peak']:.1f} % | {dil} | {r['zn_max']:+.3f} |")
    return "\n".join(L) + "\n"


def main():
    mode = sys.argv[1]
    kw = dict(a.split("=") for a in sys.argv[2:])
    Dr = float(kw.get("Dr", 0.47))
    with ProcessPoolExecutor(max_workers=8) as pool:
        if mode == "eval":
            pars = [float(kw[k]) for k in ("nb", "A0", "nd", "h0")]
            err, rows = evaluate(pars, Dr, pool)
            txt = report(rows, pars, Dr, kw.get("label", "set")) + f"\nobjective {err:.2f}\n"
            print(txt)
            tag = kw.get("label", "set")
            open(os.path.join(HERE, f"out_mono_{tag}.md"), "w", encoding="utf-8").write(txt)
            json.dump(dict(pars=pars, Dr=Dr, rows=rows), open(os.path.join(HERE, f"out_mono_{tag}.json"), "w"))
            return
        import numpy as np
        from scipy.optimize import minimize
        start = [float(v) for v in kw.get("start", "1.3,1.0,3.0,6.0").split(",")]
        names = ["nb", "A0", "nd", "h0"]
        fixed = dict((k, float(v)) for k, v in (t.split(":") for t in kw.get("fix", "").split(",") if t))
        free = [i for i, n in enumerate(names) if n not in fixed]
        for n, v in fixed.items():
            start[names.index(n)] = v
        hist = []

        def full(xf):
            pars = list(start)
            for j, i in enumerate(free):
                pars[i] = float(np.exp(xf[j]))
            return np.array(pars)

        def f(x):
            pars = full(x)
            err, _ = evaluate(pars, Dr, pool)
            hist.append((err, pars.tolist()))
            print(f"  eval {len(hist):3d}: obj {err:9.3f}  nb {pars[0]:.3f} A0 {pars[1]:.3f} nd {pars[2]:.3f} h0 {pars[3]:.3f}", flush=True)
            return err
        r = minimize(f, np.log([start[i] for i in free]), method="Nelder-Mead",
                     options=dict(maxfev=int(kw.get("maxfev", 160)), xatol=0.01, fatol=0.05, initial_simplex=None))
        best = full(r.x).tolist()
        err, rows = evaluate(best, Dr, pool)
        tag = kw.get("label", f"fit_Dr{Dr}")
        txt = report(rows, best, Dr, tag) + f"\nobjective {err:.2f} ({len(hist)} evaluations)\n"
        print(txt)
        open(os.path.join(HERE, f"out_mono_{tag}.md"), "w", encoding="utf-8").write(txt)
        json.dump(dict(pars=best, Dr=Dr, rows=rows, hist=hist), open(os.path.join(HERE, f"out_mono_{tag}.json"), "w"))


if __name__ == "__main__":
    main()
