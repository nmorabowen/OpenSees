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

            if kind == "CTXu":
                # cyclic undrained triaxial: here `emax` is the CSR; q_signed = s_xx - s_yy cycles between
                # +-2 CSR p0 under isochoric strain steps; stop at 5 % double-amplitude axial strain, p' < 1 kPa,
                # or `ncyc` cycles. Reports the cycle count N at 5 % DA.
                csr, qamp = float(emax), 2.0 * float(emax) * p0
                da, sgn, half, ncyc = 2e-5, 1.0, 0, int(spec.get("ncyc", 40))
                ea, emn, emx = 0.0, 0.0, 0.0
                rec = dict(N=None, status="ok", p_min=p0, cycles=0, trace=[])
                for k in range(2000000):
                    d = [sgn * da, -0.5 * sgn * da, -0.5 * sgn * da, 0, 0, 0]
                    res = step(d)
                    if res is None or res["rc"] != 0:
                        rec["status"] = f"refused at half-cycle {half} rc {None if res is None else res['rc']}"
                        break
                    st.update(sigma=res["sigma"], alpha=res["alpha"], alpha_in=res["alpha_in"], z=res["z"], e=res["e"])
                    ea += d[0]
                    emn, emx = min(emn, ea), max(emx, ea)
                    qs = res["sigma"][0] - res["sigma"][1]
                    rec["p_min"] = min(rec["p_min"], res["p"])
                    if k % 200 == 0:
                        rec["trace"].append((half, ea, res["p"], qs))
                    if (sgn > 0 and qs >= qamp) or (sgn < 0 and qs <= -qamp):
                        sgn, half = -sgn, half + 1
                    if emx - emn >= 0.05:
                        rec["N"] = half / 2.0; break
                    if res["p"] < 1.0:
                        rec["N"] = half / 2.0; rec["status"] = "p' < 1 kPa"; break
                    if half >= 2 * ncyc:
                        break
                rec["cycles"] = half / 2.0
                out[vname][label] = rec
                print(f"{vname:4s} {label:22s} CSR {csr}: N(5% DA) {rec['N']}  cycles run {rec['cycles']}  "
                      f"p'_min {rec['p_min']:.2f}  {rec['status']}", flush=True)
                continue
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
                    tol = 1e-9 * p0 + 1e-9
                    conv = False

                    def lat_res(x):
                        dd = [da, x, x if kind == "TXd" else 0.0, 0, 0, 0]
                        rr = step(dd)
                        return dd, rr, (None if rr is None or rr["rc"] != 0 else rr["sigma"][1] - p0)
                    for it in range(30):
                        d, res, f = lat_res(lat)
                        if f is None:
                            break
                        if abs(f) < tol:
                            conv = True; break
                        h = 1e-6 * max(abs(da), 1e-9)
                        _, r2, f2 = lat_res(lat + h)
                        if f2 is None:
                            res = r2; break
                        slope = (f2 - f) / h
                        lat -= f / slope
                    if not conv and res is not None and res["rc"] == 0:
                        # WP-150 fix: Newton can stall where the response is non-smooth (e.g. an alpha_in re-seat at
                        # phase transformation); fall back to a bracketed bisection on the lateral strain
                        lo, hi = lat - abs(da), lat + abs(da)
                        _, _, flo = lat_res(lo); _, _, fhi = lat_res(hi)
                        for _ in range(40):
                            if flo is None or fhi is None or flo * fhi <= 0:
                                break
                            lo, hi = lo - 2 * abs(da), hi + 2 * abs(da)
                            _, _, flo = lat_res(lo); _, _, fhi = lat_res(hi)
                        if flo is not None and fhi is not None and flo * fhi <= 0:
                            for _ in range(80):
                                mid = 0.5 * (lo + hi)
                                d, res, fm = lat_res(mid)
                                if fm is None:
                                    break
                                if abs(fm) < tol:
                                    break
                                if fm * flo < 0:
                                    hi = mid
                                else:
                                    lo, flo = mid, fm
                            lat = 0.5 * (lo + hi)
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
