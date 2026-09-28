"""WP-150 GATE 0: the fork's C++ LadrunoSANISAND (SAS-ME, IntScheme 129) driven through element tests by a MIXED-control
loop of `ladrunoSANISANDReplay` (WP-127) increments, at the material point (no element, no global Newton).

  python -S gate0_cxx_driver.py <dist/bin> <site-packages> <spec.json> <out.json>

spec.json = {"params": [18 LadrunoSANISAND positionals, e_init is replaced per test], "tolr": 1e-4,
             "variants": {"off": [], "on": ["-sasHFloor", 1, "-sasReseatHyst", 1, "-sasSoftCap", 0.5]},
             "tests": [[label, kind, e0, p0, eps_max], ...]}   kind: "TXu" | "TXd" | "PSd"
Compression positive (the model's convention); axial = xx. TXu: d_eps = (da, -da/2, -da/2) (isochoric).
TXd: lateral yy = zz solved each increment so that sigma_yy stays at p0. PSd: zz = 0 (plane strain), yy solved so
that sigma_yy stays at p0. Newton with a finite-difference slope, |d sigma| < 1e-9 p0 + 1e-9.
"""
import json, os, sys


def main():
    bin_dir, site, spec_path, out_path = sys.argv[1:5]
    os.add_dll_directory(bin_dir)
    sys.path.insert(0, bin_dir)
    if site not in sys.path:
        sys.path.append(site)
    import opensees as ops
    assert os.path.normcase(os.path.dirname(os.path.abspath(ops.__file__))) == os.path.normcase(os.path.abspath(bin_dir))
    spec = json.load(open(spec_path))
    P0list = list(spec["params"])
    NSTEP = int(spec.get("nstep", 1000))
    out = {"build": ops.ladrunoBuild().strip().splitlines()[0]}
    for vname, flags in spec["variants"].items():
        out[vname] = {}
        for label, kind, e0, p0, emax in spec["tests"]:
            ops.wipe()
            P = list(P0list); P[2] = float(e0)
            ops.nDMaterial("LadrunoSANISAND", 1, *P, 129, 0, 1, 1.0e-7, float(spec.get("tolr", 1e-4)),
                           "-Pmin", 0.0101, "-Presidual", 0.0, "-maxSubsteps", 20000, "-flipAlphaIn", "init",
                           *flags)
            st = dict(sigma=[p0, p0, p0, 0, 0, 0], alpha=[0.0] * 6, alpha_in=[0.0] * 6, z=[0.0] * 6, e=float(e0))

            def step(d):
                r = ops.ladrunoSANISANDReplay(1, "-convention", "compressionPositive", "-sigma", *st["sigma"],
                                              "-alpha", *st["alpha"], "-alphaIn", *st["alpha_in"],
                                              "-fabric", *st["z"], "-voidRatio", st["e"], "-dStrain", *d,
                                              "-type", "3D", "-trace", 0, "-dt", 1.0, "-primed", 1,
                                              "-prevIncrNorm", 0.0)
                r = list(r)
                if int(r[0]) != 1:
                    return None
                rc, nst = int(r[1]), int(r[2])
                s = r[6 + nst:6 + nst + 34]
                return dict(rc=rc, sigma=s[0:6], alpha=s[6:12], alpha_in=s[12:18], z=s[18:24], e=s[24], p=s[25], q=s[26])

            da = emax / NSTEP
            rec = dict(eps_a=[0.0], p=[p0], q=[0.0], e=[e0], ev=[0.0], rc=[])
            ea = ev = 0.0
            status = "ok"
            for k in range(NSTEP):
                if kind == "TXu":
                    d = [da, -0.5 * da, -0.5 * da, 0, 0, 0]
                    res = step(d)
                else:
                    lat = -0.3 * da if not rec["rc"] else lat_prev
                    for it in range(30):
                        d = [da, lat, lat if kind == "TXd" else 0.0, 0, 0, 0]
                        res = step(d)
                        if res is None or res["rc"] != 0:
                            break
                        f = res["sigma"][1] - p0
                        if abs(f) < 1e-9 * p0 + 1e-9:
                            break
                        h = 1e-6 * max(abs(da), 1e-9)
                        d2 = [da, lat + h, (lat + h) if kind == "TXd" else 0.0, 0, 0, 0]
                        r2 = step(d2)
                        if r2 is None or r2["rc"] != 0:
                            res = r2; break
                        slope = (r2["sigma"][1] - res["sigma"][1]) / h
                        lat -= f / slope
                    lat_prev = lat
                if res is None or res["rc"] != 0:
                    status = f"refused at step {k} rc {None if res is None else res['rc']}"
                    break
                st.update(sigma=res["sigma"], alpha=res["alpha"], alpha_in=res["alpha_in"], z=res["z"], e=res["e"])
                ea += da; ev += d[0] + d[1] + d[2]
                rec["eps_a"].append(ea); rec["p"].append(res["p"]); rec["q"].append(res["q"])
                rec["e"].append(res["e"]); rec["ev"].append(ev); rec["rc"].append(res["rc"])
            rec["status"] = status
            rec.pop("rc")
            out[vname][label] = rec
            print(f"{vname:4s} {label:22s} {status:28s} q_peak {max(rec['q']):8.1f} q_end {rec['q'][-1]:8.1f} "
                  f"p_end {rec['p'][-1]:8.1f} e_end {rec['e'][-1]:.4f}", flush=True)
    json.dump(out, open(out_path, "w"))


if __name__ == "__main__":
    main()
