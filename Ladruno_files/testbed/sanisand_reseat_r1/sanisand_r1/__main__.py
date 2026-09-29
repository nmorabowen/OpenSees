"""CLI for the SANISAND reference integrator (WP-134).

    cd Ladruno_scripts
    python -m sanisand_reference increment --sigma 20 20 20 0 0 0 --alpha 0 0 0 0 0 0 \\
        --e 0.7 --deps 0 1e-4 0 0 0 0 [--options uw|uw_me|paper] [--alpha-in ...] [--z ...]
    python -m sanisand_reference reproducer
    python -m sanisand_reference ring --mesh b8 --row 0 --probe shear --delta 1e-5
    python -m sanisand_reference triaxial --set toyoura --p0 100 --e0 0.833 --strain 0.3 --undrained

All stresses COMPRESSION POSITIVE (kPa), tensor components xx yy zz xy yz zx;
strains Voigt with ENGINEERING shear.  Prints a JSON summary."""
import argparse
import json
import sys

from .driver import triaxial
from .integrator import Control, integrate
from .model import CAMPAIGN, TOYOURA, Options, State
from .ring import REPRODUCER_DEPS, load_ring_csv, probes, reproducer_state, row_state

PRESETS = {"paper": Options, "uw": Options.uw, "uw_me": Options.uw_me}


def _opts(a):
    O = PRESETS[a.options]()
    if a.alpha_in_rule:
        O = O.with_(alpha_in_rule=a.alpha_in_rule)
    return O


def main(argv=None):
    ap = argparse.ArgumentParser(prog="sanisand_reference", description=__doc__.split("\n")[0])
    sp = ap.add_subparsers(dest="cmd", required=True)
    for name in ("increment", "reproducer", "ring", "triaxial"):
        s = sp.add_parser(name)
        s.add_argument("--options", choices=list(PRESETS), default="paper")
        s.add_argument("--alpha-in-rule", choices=["paper", "uw"], default=None)
        s.add_argument("--rtol", type=float, default=1e-10)
        if name == "increment":
            s.add_argument("--sigma", type=float, nargs=6, required=True)
            s.add_argument("--alpha", type=float, nargs=6, required=True)
            s.add_argument("--alpha-in", type=float, nargs=6, default=[0.0] * 6)
            s.add_argument("--z", type=float, nargs=6, default=[0.0] * 6)
            s.add_argument("--e", type=float, required=True)
            s.add_argument("--deps", type=float, nargs=6, required=True)
        if name == "ring":
            s.add_argument("--mesh", default="b8")
            s.add_argument("--row", type=int, default=0)
            s.add_argument("--probe", choices=["isoComp", "shear"], default="shear")
            s.add_argument("--delta", type=float, default=1e-5)
        if name == "triaxial":
            s.add_argument("--set", choices=["toyoura", "campaign"], default="toyoura")
            s.add_argument("--p0", type=float, default=100.0)
            s.add_argument("--e0", type=float, default=0.833)
            s.add_argument("--strain", type=float, default=0.3)
            s.add_argument("--undrained", action="store_true")
    a = ap.parse_args(argv)
    O = _opts(a)
    if a.cmd == "increment":
        st = State.from_voigt(a.sigma, a.alpha, a.z, a.e, a.alpha_in)
        r = integrate(st, Control.strain(a.deps), CAMPAIGN, O, rtol=a.rtol)
    elif a.cmd == "reproducer":
        r = integrate(reproducer_state(), REPRODUCER_DEPS, CAMPAIGN, O, rtol=a.rtol)
    elif a.cmd == "ring":
        row = load_ring_csv(a.mesh)[a.row]
        r = integrate(row_state(row), probes(a.delta)[a.probe], CAMPAIGN, O, rtol=a.rtol)
    else:
        P = TOYOURA if a.set == "toyoura" else CAMPAIGN
        r, tab = triaxial(P, a.p0, a.e0, a.strain, drained=not a.undrained, O=O,
                          rtol=max(a.rtol, 1e-8))
    out = r.summary()
    out["start"] = r.start
    out["end"] = r.end
    out["segments"] = r.segments
    out["reseats"] = r.reseats
    json.dump(out, sys.stdout, indent=1, default=float)
    print()


if __name__ == "__main__":
    main()
