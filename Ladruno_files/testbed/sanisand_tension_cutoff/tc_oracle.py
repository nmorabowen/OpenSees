"""WP-152 tension cutoff (separation): the state machine as a DEFINITION, driven around
the UNMODIFIED WP-134 oracle (the WP-151 testbed copy `sanisand_r1`, which carries the R1
toggles).  The oracle's `p_floor` stop -- the exact trajectory reaching p -> 0 inside an
increment -- is the counterpart of SAS-ME's tension refusal (code 6), i.e. trigger E1.
The oracle integrates exactly, so it has no accuracy/cost failures: trigger E2 (codes 4/9
at p0 < p_sep) has no oracle counterpart and is exercised on the C++ side only.

State machine (Ladruno_implementation/152_sanisand_tension_cutoff.md):
  NORMAL: integrate the increment.  status p_floor -> SEPARATED at the END of the increment:
          sigma = (p_min - p_r) I (model p = p_min), alpha = alpha_in = 0, fabric kept,
          e follows the strain, tr(eps_entry) := tr(eps) at the end of this increment.
  SEPARATED: strain absorbed; g = tr(eps) - tr(eps_entry) (compression positive).
          g >= g_c = (p_contact - p_min) / K(p_contact) -> NORMAL at the end of this
          increment with sigma = (p_re - p_r) I, p_re = p_min + K(p_contact) g, alpha =
          alpha_in = 0.  Else sigma stays (p_min - p_r) I.
Strain increments are compression positive, Voigt engineering shear (the oracle's)."""
from __future__ import annotations

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "sanisand_reseat_r1"))
import r1common as C  # noqa: E402
from sanisand_r1 import State, integrate  # noqa: E402
from sanisand_r1.model import elastic_moduli, t2v  # noqa: E402

I3 = np.eye(3)
Z3 = np.zeros((3, 3))


def _e_after(e, dv, P, O):
    ref = P.e_init if O.void_ratio_law == "initial" else e
    return e - (1.0 + ref) * dv


def model_p(st, O):
    return float(np.trace(st.sigma)) / 3.0 + O.p_residual


def drive(state0, deps_list, O, p_sep, p_contact, P=None):
    """Returns one record per increment: mode after it ('N'/'S'), event, sigma (Voigt,
    compression positive), model p, g, the oracle status."""
    P = C.P if P is None else P
    st, mode, g, tr_entry = state0, "N", 0.0, None
    out = []
    for k, de in enumerate(deps_list):
        de = [float(x) for x in de]
        dv = de[0] + de[1] + de[2]
        ev, status = None, None
        if mode == "N":
            r = integrate(st, de, P, O, record=False)
            status = r.status
            if status == "ok":
                st = r.state
            elif status == "p_floor":                                   # E1 (tension)
                e_end = _e_after(st.e, dv, P, O)
                st = State((O.p_min - O.p_residual) * I3, Z3.copy(), st.z.copy(), e_end, Z3.copy())
                mode, g, ev = "S", 0.0, "enter_tension"
            else:
                out.append(dict(k=k, mode=mode, event="refused", status=status))
                break
        else:
            e_end = _e_after(st.e, dv, P, O)
            g += dv
            _, K = elastic_moduli(p_contact, e_end, P, O)
            gc = (p_contact - O.p_min) / K
            if g >= gc:
                p_re = O.p_min + K * g
                st = State((p_re - O.p_residual) * I3, Z3.copy(), st.z.copy(), e_end, Z3.copy())
                mode, ev = "N", "recontact"
            else:
                st = State((O.p_min - O.p_residual) * I3, Z3.copy(), st.z.copy(), e_end, Z3.copy())
        out.append(dict(k=k, mode=mode, event=ev, status=status, sigma=t2v(st.sigma).tolist(),
                        p=model_p(st, O), g=g, e=st.e))
    return out


def work(records, deps_list, sigma0):
    """Net work sum sigma_mid : deps (trapezoid over increments; Voigt engineering shear, so
    the shear products count once)."""
    s_prev = np.array(sigma0, dtype=float)
    W = 0.0
    for r, de in zip(records, deps_list):
        if "sigma" not in r:
            break
        s = np.array(r["sigma"])
        W += float(np.dot(0.5 * (s + s_prev), np.array(de, dtype=float)))
        s_prev = s
    return W


# ---------------------------------------------------------------- the element paths
def iso_state(p0, e0, O):
    return State((p0 - O.p_residual) * I3, Z3.copy(), Z3.copy(), e0, Z3.copy())


def path_iso(n_out=40, n_hold=10, n_back=60, d=1.0e-5):
    """Isotropic unloading past separation, then reloading past re-contact."""
    return [[-d, -d, -d, 0, 0, 0]] * n_out + [[-d, -d, -d, 0, 0, 0]] * n_hold + [[d, d, d, 0, 0, 0]] * n_back


def path_te(n_out=60, n_back=90, d=1.0e-5):
    """Uniaxial extension along x (lateral strains fixed) past separation, then back."""
    return [[-d, 0, 0, 0, 0, 0]] * n_out + [[d, 0, 0, 0, 0, 0]] * n_back


def path_cycles(n_cyc=4, n_leg=60, d=1.0e-5):
    """Isotropic open/close cycles across separation and re-contact."""
    out = []
    for _ in range(n_cyc):
        out += [[-d, -d, -d, 0, 0, 0]] * n_leg + [[d, d, d, 0, 0, 0]] * n_leg
    return out


def fixture(worktree):
    """tests/data/wp152_oracle_paths.json: the three element paths, from the post-flip
    state the C++ driver recorded (tests/data/wp152_flip_state.json), R1 (the full set) +
    the cutoff (p_sep 0.5, p_contact 1.0).  The increments are those of
    tests/wp152_cutoff_tools.paths(), stored with the result so a drift is caught."""
    import json
    data = os.path.join(worktree, "tests", "data")
    flip = json.load(open(os.path.join(data, "wp152_flip_state.json")))
    O = C.variants()["T1B1S"]
    p_sep, p_contact = 0.5, 1.0
    s = flip["sigma"]
    sig = np.array([[s[0], s[3], s[5]], [s[3], s[1], s[4]], [s[5], s[4], s[2]]])
    st0 = State(sig, Z3.copy(), Z3.copy(), flip["e"], Z3.copy())
    out = dict(note="WP-152 oracle expectations: tc_oracle.py --fixture (unmodified WP-134 oracle, "
                    "R1 T1B1S, cutoff p_sep 0.5 p_contact 1.0, p_min 0.0101, p_r 0)",
               flip=flip, p_sep=p_sep, p_contact=p_contact, paths={})
    for name, path in (("iso", path_iso()), ("te", path_te()), ("cyc", path_cycles())):
        rec = drive(st0, path, O, p_sep, p_contact)
        out["paths"][name] = dict(increments=path, records=rec,
                                  events=[(r["k"], r["event"]) for r in rec if r.get("event")],
                                  work=work(rec, path, t2v(st0.sigma)))
        print(name, out["paths"][name]["events"], f"work {out['paths'][name]['work']:.4g}")
    json.dump(out, open(os.path.join(data, "wp152_oracle_paths.json"), "w"))


if __name__ == "__main__":
    if len(sys.argv) > 2 and sys.argv[1] == "--fixture":
        fixture(sys.argv[2])
        sys.exit(0)
    O = C.variants()["T1B1S"]                   # R1, the full set (the campaign setting)
    p_sep, p_contact = 0.5, 1.0
    e0 = C.P.e_init
    for name, path, p0 in (("iso", path_iso(), 2.0), ("te", path_te(), 2.0), ("cyc", path_cycles(), 2.0)):
        st0 = iso_state(p0, e0, O)
        rec = drive(st0, path, O, p_sep, p_contact)
        evs = [(r["k"], r["event"]) for r in rec if r.get("event")]
        print(f"{name}: {len(rec)} steps; events {evs}; last mode {rec[-1].get('mode')}, p_end "
              f"{rec[-1].get('p', float('nan')):.4g}; net work {work(rec, path, t2v(st0.sigma)):.4g}")
