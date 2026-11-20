"""WP-229 (TIMs F26) gate, C++ side: the 17 DM04 Toyoura triaxial simulations on LadrunoSANISAND
(SAS-ME, IntScheme 129) with `-sasPlasticExp nh`, at the material point through `ladrunoSANISANDReplay`.

  py -3.12 -S gate_nh_cxx.py <dist/bin> <site-packages>      -> out_gate_nh_cxx.json, out_gate_nh_cxx.md

The matrix, the constants and the 500 equal axial increments are those of gate_nh_oracle.py. Undrained:
d_eps = (da, -da/2, -da/2). Drained: the lateral strain is solved each increment so that sigma_yy stays at
p0 (a bracketed Illinois root on the lateral strain).

Variants: R1 off, and R1 on (-sasHFloor 1 -sasReseatHyst 1 -sasSoftCap 0.5), each at nh 0.5 / 0.2 / 0.25 /
0.3 (and 0.35 .. 0.45 for the trend). "nh 0.5" is run twice: flag absent and `-sasPlasticExp 0.5`; the two
must agree bit for bit. Compared with out_gate_nh_oracle.json when it is present.
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
NH = (0.5, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45)
NSTEP = 500
TOYOURA = [125.0, 0.05, 0.8, 1.25, 0.712, 0.019, 0.934, 0.7, 100.0, 0.01, 7.05, 0.968, 1.1, 0.704, 3.5,
           4.0, 600.0, 2.0]
R1 = ("-sasHFloor", 1.0, "-sasReseatHyst", 1.0, "-sasSoftCap", 0.5)
TESTS = ([(f"F5_u_e0.735_p{p}", False, 0.735, float(p)) for p in (100, 1000, 2000, 3000)]
         + [(f"F6_u_e0.833_p{p}", False, 0.833, float(p)) for p in (100, 1000, 2000, 3000)]
         + [(f"F7_u_e0.907_p{p}", False, 0.907, float(p)) for p in (100, 1000, 2000)]
         + [(f"F8_d_p500_e{e}", True, e, 500.0) for e in (0.960, 0.886, 0.810)]
         + [(f"F9_d_p100_e{e}", True, e, 100.0) for e in (0.996, 0.917, 0.831)])


def run(ops, flags, drained, e0, p0):
    ops.wipe()
    P = list(TOYOURA)
    P[2] = float(e0)
    ops.nDMaterial("LadrunoSANISAND", 1, *P, 129, 0, 1, 1.0e-7, 1.0e-4, "-Pmin", 0.0101, "-Presidual", 0.0,
                   "-maxSubsteps", 20000, "-flipAlphaIn", "init", *flags)
    st = dict(sigma=[p0, p0, p0, 0, 0, 0], alpha=[0.0] * 6, alpha_in=[0.0] * 6, z=[0.0] * 6, e=float(e0))
    cen = dict(hFloored=0.0, hSoftCapped=0.0, stages=0.0, max_lat_res=0.0)

    def step(d):
        r = list(ops.ladrunoSANISANDReplay(
            1, "-convention", "compressionPositive", "-sigma", *st["sigma"], "-alpha", *st["alpha"],
            "-alphaIn", *st["alpha_in"], "-fabric", *st["z"], "-voidRatio", st["e"], "-dStrain", *d,
            "-type", "3D", "-trace", 0, "-dt", 1.0, "-primed", 1, "-prevIncrNorm", 0.0))
        if int(r[0]) != 1:
            return None
        rc, nst, nrec, width = int(r[1]), int(r[2]), int(r[3]), int(r[4])
        s = r[6 + nst:6 + nst + 34]
        i = 6 + nst + 34 + nrec * width
        sas = r[i + 2:i + 2 + int(r[i + 1])] if len(r) > i + 1 and int(r[i]) == 129 else []
        return dict(rc=rc, sigma=s[0:6], alpha=s[6:12], alpha_in=s[12:18], z=s[18:24], e=s[24], p=s[25],
                    q=s[26], sas=sas)

    da = 0.25 / NSTEP
    q, p, status = [0.0], [p0], "ok"
    for k in range(NSTEP):
        if not drained:
            res = step([da, -0.5 * da, -0.5 * da, 0, 0, 0])
        else:
            tol = 1e-9 * p0 + 1e-9

            def g(x):
                rr = step([da, x, x, 0, 0, 0])
                return rr, (None if rr is None or rr["rc"] != 0 else rr["sigma"][1] - p0)
            # Bracketed root (Illinois) on the lateral strain in [-2 da, da/2]: sigma_yy is monotonic in it.
            # No secant on a finite-difference slope: at SAS-ME's TolR the residual is not smooth at the
            # 1e-6 da scale, the slope is noise and the iterate leaves the admissible range.
            lo, hi = -2.0 * da, 0.5 * da
            res, flo = g(lo)
            res, fhi = g(hi)
            f = None
            if flo is not None and fhi is not None and flo * fhi < 0:
                # stop on the residual OR on the bracket width: SAS-ME's result has jumps of the order of
                # TolR * sigma in the increment (a different substep sequence), so a residual below them is
                # not always reachable. The residual left is recorded (max_lat_res, kPa).
                for _ in range(400):
                    mid = (lo * fhi - hi * flo) / (fhi - flo)
                    if not (min(lo, hi) < mid < max(lo, hi)):
                        mid = 0.5 * (lo + hi)
                    res, f = g(mid)
                    if f is None or abs(f) < tol or abs(hi - lo) < 1e-10 * da:
                        break
                    if f * fhi < 0:
                        lo, flo = hi, fhi
                    else:
                        flo *= 0.5
                    hi, fhi = mid, f
            if f is None:
                status = f"lateral control failed at step {k}"
                break
            cen["max_lat_res"] = max(cen["max_lat_res"], abs(f))
        if res is None or res["rc"] != 0:
            status = f"refused at step {k} rc {None if res is None else res['rc']}"
            break
        st.update(sigma=res["sigma"], alpha=res["alpha"], alpha_in=res["alpha_in"], z=res["z"], e=res["e"])
        q.append(res["q"])
        p.append(res["p"])
        if len(res["sas"]) > 34:      # sasStats columns 33 / 34: hFloored / hSoftCapped (per replay)
            cen["hFloored"] += res["sas"][33]
            cen["hSoftCapped"] += res["sas"][34]
            cen["stages"] += 2.0 * res["sas"][2]
    return dict(status=status, q_peak=max(q), q_end=q[-1], p_end=p[-1], e_end=st["e"], q=q, p=p, census=cen,
                state_hex=[float(x).hex() for x in st["sigma"] + st["alpha"] + st["alpha_in"] + st["z"]])


def main():
    bin_dir, site = sys.argv[1:3]
    os.add_dll_directory(bin_dir)
    sys.path.insert(0, bin_dir)
    if site not in sys.path:
        sys.path.append(site)
    import opensees as ops
    assert os.path.normcase(os.path.dirname(os.path.abspath(ops.__file__))) == os.path.normcase(
        os.path.abspath(bin_dir))
    out = dict(build=ops.ladrunoBuild().strip().splitlines()[0], runs={})
    for r1name, r1 in (("R1off", ()), ("R1on", R1)):
        for nh in ("absent",) + NH:
            flags = tuple(r1) + (() if nh == "absent" else ("-sasPlasticExp", nh))
            for label, drained, e0, p0 in TESTS:
                rec = run(ops, flags, drained, e0, p0)
                out["runs"][f"{r1name}|{nh}|{label}"] = rec
                print(f"{r1name} nh {nh} {label:18s} {rec['status']:24s} q_peak {rec['q_peak']:8.1f} "
                      f"q_end {rec['q_end']:8.1f}", flush=True)
    json.dump(out, open(os.path.join(HERE, "out_gate_nh_cxx.json"), "w"))
    report(out)


def report(out):
    R = out["runs"]
    orc = {}
    op = os.path.join(HERE, "out_gate_nh_oracle.json")
    if os.path.isfile(op):
        orc = {(r["nh"], r["label"]): r for r in json.load(open(op))}
    L = [f"build {out['build']}", ""]
    ident = all(R[f"{v}|absent|{t[0]}"]["state_hex"] == R[f"{v}|0.5|{t[0]}"]["state_hex"]
                and R[f"{v}|absent|{t[0]}"]["q"] == R[f"{v}|0.5|{t[0]}"]["q"]
                for v in ("R1off", "R1on") for t in TESTS)
    L.append(f"`-sasPlasticExp 0.5` against the flag absent, 34 runs x 500 increments: "
             f"{'bit-identical (every q, and the end state by float.hex)' if ident else 'DIFFERENT'}")
    for v in ("R1off", "R1on"):
        L += ["", f"### {v}: C++ ratios to the C++ nh = 0.5 run (q_peak / q at 25 %)", "",
              "| test | nh 0.5 q_peak / q_25% | " + " | ".join(f"nh {n}" for n in NH[1:]) + " |",
              "|---|---|" + "---|" * (len(NH) - 1)]
        bad = []
        for label, *_ in TESTS:
            ref = R[f"{v}|0.5|{label}"]
            row = f"| {label} | {ref['q_peak']:.1f} / {ref['q_end']:.1f} ({ref['status']}) |"
            for n in NH[1:]:
                r = R[f"{v}|{n}|{label}"]
                a, b = r["q_peak"] / ref["q_peak"], r["q_end"] / ref["q_end"]
                okr = r["status"] == "ok" and ref["status"] == "ok"
                if n in (0.2, 0.25, 0.3) and (not okr or abs(a - 1) > 0.1 or abs(b - 1) > 0.1):
                    bad.append((n, label, a, b, r["status"]))
                row += f" {a:.4f} / {b:.4f}{'' if okr else ' (' + r['status'] + ')'} |"
            L.append(row)
        L += ["", f"outside +-10 % (or not ok) at nh 0.2 / 0.25 / 0.3: {len(bad)} of 51"]
        L += [f"- nh {b[0]} {b[1]}: q_peak x{b[2]:.4f}, q_25% x{b[3]:.4f}, {b[4]}" for b in bad]
    if orc:
        L += ["", "### C++ (R1 off) against the oracle at the same nh: max |q_cxx - q_oracle| / q_max over the curve",
              "", "| test | " + " | ".join(f"nh {n}" for n in NH) + " |", "|---|" + "---|" * len(NH)]
        worst = 0.0
        for label, *_ in TESTS:
            row = f"| {label} |"
            for n in NH:
                r, o = R[f"R1off|{n}|{label}"], orc.get((n, label))
                if o is None or r["status"] != "ok" or o["status"] != "ok" or len(r["q"]) != len(o["q"]):
                    row += " n/a |"
                    continue
                d = max(abs(a - b) for a, b in zip(r["q"], o["q"])) / max(o["q"])
                worst = max(worst, d)
                row += f" {d:.2e} |"
            L.append(row)
        L += ["", f"worst: {worst:.2e}"]
    L += ["", "largest lateral-stress residual accepted in a drained run (kPa; the bracket closed on a jump "
          f"of the integrator's output): {max(r['census']['max_lat_res'] for r in R.values()):.2e}"]
    L += ["", "### R1 on against R1 off at the same nh (max |dq| / q_max), and the R1 census of the R1-on run",
          "", "| test | " + " | ".join(f"nh {n}: dq; floored / capped stages" for n in (0.5, 0.25)) + " |",
          "|---|---|---|"]
    for label, *_ in TESTS:
        row = f"| {label} |"
        for n in (0.5, 0.25):
            a, b = R[f"R1on|{n}|{label}"], R[f"R1off|{n}|{label}"]
            if a["status"] != "ok" or b["status"] != "ok":
                row += f" {a['status']} / {b['status']} |"
                continue
            d = max(abs(x - y) for x, y in zip(a["q"], b["q"])) / max(b["q"])
            c = a["census"]
            row += f" {d:.2e}; {c['hFloored']:.0f} / {c['hSoftCapped']:.0f} of {c['stages']:.0f} |"
        L.append(row)
    txt = "\n".join(L)
    open(os.path.join(HERE, "out_gate_nh_cxx.md"), "w").write(txt + "\n")
    print(txt)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--report":
        report(json.load(open(os.path.join(HERE, "out_gate_nh_cxx.json"))))
    else:
        main()
