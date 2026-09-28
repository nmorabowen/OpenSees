"""WP-150 T5: element-level physics of the TIMs campaign SANISAND set (flag OFF baseline).

Drained plane-strain and triaxial compression from isotropic states at p0 = 10, 50, 150, 500 kPa,
e0 = 0.6944 (the footing's e_init), integrated EXACTLY by the WP-134 reference integrator
(`Ladruno_scripts/sanisand_reference`, Radau IIA, uw_model options = the UW constitutive additions
U1-U5 with the paper's alpha_in rule: the oracle the C++ SAS-ME reproduces).

Physics checks that need NO relative density (Bolton 1986, Geotechnique 36(1)):
  plane strain : phi_max - phi_cs = 0.8 * psi_max          (psi_max = peak dilation angle)
  triaxial     : phi_max - phi_cs = 3 I_R,  (-de_v/de_1)_max = 0.3 I_R  ->  phi_max - phi_cs ~ 10 (-de_v/de_1)_max [deg]
and the pressure dependence: the D_r implied by Bolton's I_R = D_r (10 - ln p') - 1 from each p0
should be the same at every p0 if the model's pressure dependence is sand-like.
Compression positive; axial = x.  Usage:  python t5_element_physics.py [eps_axial_max] [A0=... etc]
"""
import json, math, os, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "Ladruno_scripts"))

import numpy as np  # noqa: E402
from sanisand_reference import CAMPAIGN, TOYOURA, e_critical  # noqa: E402
from sanisand_reference.driver import isotropic_state  # noqa: E402
from sanisand_reference.integrator import Control, integrate, path_table  # noqa: E402
from sanisand_reference.ring import ring_variants  # noqa: E402

EPS_MAX = float(sys.argv[1]) if len(sys.argv) > 1 and "=" not in sys.argv[1] else 0.25
over = dict(a.split("=") for a in sys.argv[1:] if "=" in a)
SET = over.pop("SET", "campaign")
E0 = float(over.pop("E0", 0.6944))
P = {"campaign": CAMPAIGN, "toyoura": TOYOURA}[SET]
if over:
    from dataclasses import replace
    P = replace(P, **{k: float(v) for k, v in over.items()})
O = ring_variants()["uw_model"]
P0S = (10.0, 50.0, 150.0, 500.0)
TESTS = {
    # mask: True = strain-controlled component, value = its total increment; False = d(sigma) = 0
    "PS": ((True, False, True, True, True, True), (EPS_MAX, 0.0, 0.0, 0.0, 0.0, 0.0)),
    "TX": ((True, False, False, True, True, True), (EPS_MAX, 0.0, 0.0, 0.0, 0.0, 0.0)),
}


def mobilized(sig):
    s1 = sig[0]
    s3 = min(sig[1], sig[2])
    return math.degrees(math.asin(max(-1.0, min(1.0, (s1 - s3) / (s1 + s3)))))


def analyse(tab, kind):
    ea = np.array([t["eps"][0] for t in tab])
    eps = np.array([t["eps"][:3] for t in tab])
    ev = eps.sum(1)
    sig = [t["sigma"] for t in tab]
    phi = np.array([mobilized(s) for s in sig])
    ip = int(np.argmax(phi))
    # dilatancy by finite differences on the recorded path (compression positive -> dilation = dev < 0)
    de1 = np.gradient(ea); dev = np.gradient(ev)
    if kind == "PS":
        de3 = np.gradient(eps[:, 1])
        with np.errstate(invalid="ignore", divide="ignore"):
            sinpsi = -dev / (de1 - de3)
        dil = np.degrees(np.arcsin(np.clip(sinpsi, -1, 1)))
    else:
        with np.errstate(invalid="ignore", divide="ignore"):
            dil = -dev / de1
    dmax = float(np.nanmax(dil))
    return dict(phi_peak=float(phi[ip]), eps_a_peak=float(ea[ip]), phi_end=float(phi[-1]),
                eta_peak=float(tab[ip]["eta"]), psi0=float(tab[0]["psi"]), psi_end=float(tab[-1]["psi"]),
                ev_end=float(ev[-1]), dil_max=dmax, dil_at_peak=float(dil[ip]),
                eps_a_end=float(ea[-1]))


def main():
    out, rows = {}, []
    phi_cs_tx = math.degrees(math.asin(3 * P.Mc / (6 + P.Mc)))
    for kind, (mask, vals) in TESTS.items():
        for p0 in P0S:
            t0 = time.time()
            st = isotropic_state(p0, E0)
            res = integrate(st, Control(mask, vals), P, O, rtol=1e-8)
            tab = path_table(res, P, O)
            a = analyse(tab, kind)
            a.update(kind=kind, p0=p0, status=res.status, wall=time.time() - t0)
            rows.append(a)
            out[f"{kind}_{int(p0)}"] = [dict(eps_a=t["eps"][0], ev=sum(t["eps"][:3]), p=t["p"], q=t["q"],
                                             phi=mobilized(t["sigma"]), psi=t["psi"]) for t in tab]
            print(f"{kind} p0={p0:5.0f} status={res.status:>14s} phi_peak={a['phi_peak']:.2f} @ eps_a={a['eps_a_peak']:.4f} "
                  f"phi_end={a['phi_end']:.2f} psi0={a['psi0']:+.3f} psi_end={a['psi_end']:+.3f} "
                  f"dil_max={a['dil_max']:.3f} ev_end={a['ev_end']:+.4f} ({a['wall']:.0f}s)", flush=True)
    # Bolton checks
    L = ["| test | p0 kPa | ψ0 | φ′_peak ° | ε_a at peak | φ′ at ε_a end ° | ψ end | max dilatancy | Bolton check |",
         "|---|---|---|---|---|---|---|---|---|"]
    for a in rows:
        if a["kind"] == "PS":
            phi_cs = a["phi_end"] if abs(a["psi_end"]) < 0.01 else float("nan")
            pred = 0.8 * a["dil_max"]
            chk = f"Δφ = {a['phi_peak'] - a['phi_end']:.1f}° vs 0.8·ψ_max = {pred:.1f}°"
            dil = f"ψ_max = {a['dil_max']:.1f}°"
        else:
            pred = 10.0 * a["dil_max"]
            chk = f"Δφ = {a['phi_peak'] - phi_cs_tx:.1f}° vs 10·(−dε_v/dε_1)max = {pred:.1f}°"
            dil = f"(−dε_v/dε_1)max = {a['dil_max']:.3f}"
            # D_r implied by Bolton's triaxial 3 I_R, I_R = D_r (10 - ln p) - 1
            IR = (a["phi_peak"] - phi_cs_tx) / 3.0
            Dr = (IR + 1.0) / (10.0 - math.log(a["p0"]))
            chk += f"; implied D_r (Δφ = 3 I_R) = {Dr:.2f}"
        L.append(f"| {a['kind']} | {a['p0']:.0f} | {a['psi0']:+.3f} | {a['phi_peak']:.1f} | {a['eps_a_peak']:.4f} | "
                 f"{a['phi_end']:.1f} | {a['psi_end']:+.3f} | {dil} | {chk} |")
    txt = (f"T5 {SET} set{(' with ' + str(over)) if over else ''}, e0 = {E0}, ε_a to {EPS_MAX}; "
           f"φ′_cs (triaxial, from Mc) = {phi_cs_tx:.2f}°\n\n" + "\n".join(L) + "\n")
    print(txt)
    tag = "_".join([SET] * (SET != "campaign") + [f"e{E0}"] * (E0 != 0.6944) + [f"{k}{v}" for k, v in over.items()])
    open(os.path.join(HERE, f"out_t5{('_' + tag) if tag else ''}.md"), "w", encoding="utf-8").write(txt)
    json.dump(dict(rows=rows, paths=out), open(os.path.join(HERE, f"out_t5{('_' + tag) if tag else ''}.json"), "w"))


if __name__ == "__main__":
    main()
