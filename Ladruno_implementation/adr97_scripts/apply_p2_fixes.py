"""ADR-97 P2 (wp/97c) -- follow-up C++ edits found by running the gates.

1. The admissibility self-check of the principal-space return was written against
   `yf_tolerance()` alone, which defaults to the ABSOLUTE `f_absolute_tol = 1e-6`
   (`f_relative_tol` defaults to 0).  On the ADR-84 MCTC deck -- kPa, |sigma| ~
   5.4e3, strength scale c cos(phi) = 94 -- the header's own `f` recomputed from
   the reassembled `Q diag(y) Q^T` came back at 3.4e-6, i.e. 6e-10 RELATIVE, and
   the step was refused.  That is round-off, not a bad return, and it is
   amplified because an EDGE return lands exactly on a corner where the Lode
   angle is ill conditioned (`dA/dtheta` is not stationary there and
   `dtheta/dJ3 ~ 1/cos(3 theta)` diverges).  Refusing on it is ADR-94 M5 all over
   again: the same problem in Pa and in kPa would disagree.

   The check is replaced by TWO checks:
   * an EXACT one in principal space -- all three surface functions evaluated at
     the returned principal point, which is where the return was computed and is
     perfectly conditioned.  This is the one that catches a coding error (a wrong
     region, an edge return that overshoots past another surface);
   * the header's COMPOSITE `f` -- the real gate for MohrCoulombTensionCutoff,
     where the plain-MC return must still respect the cutoff plane -- at a
     tolerance that is RELATIVE to the returned stress magnitude.  A genuine
     fall-through error there is O(|sigma|); round-off is not.

Run from the worktree root, after apply_p2_core.py:
    python3.12 Ladruno_implementation/adr97_scripts/apply_p2_fixes.py
"""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SRC = os.path.join(ROOT, "SRC", "material", "nD", "ASDPlasticMaterial3D")

EDITS = []


def edit(relpath, anchor, new, tag):
    EDITS.append((relpath, anchor, new, tag))


edit(
    "ASDPlasticMaterial3D.h",
    """        // Admissibility against the header's OWN composite f.  For plain
        // Mohr-Coulomb this is a self-check (the return is exact by
        // construction); for MohrCoulombTensionCutoff it is the real gate: this
        // path is only reached after `special_return` declined, and the MC return
        // it performs must still respect the cutoff plane.
        {
            const double f_ret = yf(sigma_ret, iv_storage, parameters_storage);
            if (!(f_ret <= tol_f))
            {
                opserr << "ASDPlasticMaterial3D::Closest_Point (tag " << ASDP_TAG
                       << ") - the principal-space return landed OUTSIDE this yield"
                       << " function's own surface (f = " << f_ret << " > tol = "
                       << tol_f << ", region " << region
                       << ") -- rejecting step (ADR-97 P2)" << endln;
                return LADRUNO_MATERIAL_REFUSED;
            }
        }
""",
    """        // Admissibility, in TWO places, for two different reasons.
        //
        // (a) EXACTLY, in principal space, against ALL THREE surfaces of the
        //     sextant.  This is where the return was computed and it is perfectly
        //     conditioned, so it is the check that catches a coding error -- a
        //     misclassified region, or an edge return that overshoots past a
        //     third surface.
        //
        // (b) LOOSELY, against the header's own COMPOSITE f.  For plain
        //     Mohr-Coulomb (b) is redundant with (a); for
        //     MohrCoulombTensionCutoff it is the real gate, because this path is
        //     only reached after `special_return` declined and the plain-MC
        //     return it then performs must still respect the cutoff plane.  Its
        //     tolerance is RELATIVE to the returned stress magnitude, NOT the
        //     bare `yf_tolerance()`: that accessor defaults to the ABSOLUTE
        //     `f_absolute_tol = 1e-6` (`f_relative_tol` defaults to 0), and on
        //     the ADR-84 MCTC deck (kPa, |sigma| ~ 5.4e3, strength scale 94) the
        //     header's f recomputed from the reassembled `Q diag(y) Q^T` is
        //     3.4e-6 -- 6e-10 relative, i.e. round-off, amplified because an EDGE
        //     return lands exactly on a corner where the Lode angle is ill
        //     conditioned (dtheta/dJ3 ~ 1/cos(3 theta)).  Refusing on that would
        //     be ADR-94 M5 again: the same model in Pa and in kPa disagreeing.
        //     A genuine fall-through error is O(|sigma|), which 1e-8 relative
        //     still catches by eight orders of magnitude.
        double sig_max = 0.0;
        for (int i = 0; i < 6; ++i)
        {
            const double a = (sigma_ret(i) < 0) ? -sigma_ret(i) : sigma_ret(i);
            if (a > sig_max) sig_max = a;
        }
        double scale_ref = yf.strength_scale(iv_storage, parameters_storage);
        if (scale_ref < 0) scale_ref = -scale_ref;
        const double ref_mag = (sig_max > scale_ref) ? sig_max : scale_ref;
        {
            double f_princ = -1e300;
            for (int t = 0; t < 3; ++t)
            {
                const double v = av[t].dot(yv) - k_coh;
                if (v > f_princ) f_princ = v;
            }
            if (!(f_princ <= 1e-10 * ref_mag))
            {
                opserr << "ASDPlasticMaterial3D::Closest_Point (tag " << ASDP_TAG
                       << ") - the principal-space return landed OUTSIDE the"
                       << " Mohr-Coulomb cone (max_k a_k . y - k = " << f_princ
                       << ", region " << region
                       << ") -- rejecting step (ADR-97 P2)" << endln;
                return LADRUNO_MATERIAL_REFUSED;
            }
        }
        {
            const double f_ret = yf(sigma_ret, iv_storage, parameters_storage);
            const double tol_adm = (tol_f > 1e-8 * ref_mag) ? tol_f : 1e-8 * ref_mag;
            if (!(f_ret <= tol_adm))
            {
                opserr << "ASDPlasticMaterial3D::Closest_Point (tag " << ASDP_TAG
                       << ") - the principal-space return landed OUTSIDE this yield"
                       << " function's own surface (f = " << f_ret << " > tol = "
                       << tol_adm << ", region " << region
                       << ") -- rejecting step (ADR-97 P2)" << endln;
                return LADRUNO_MATERIAL_REFUSED;
            }
        }
""",
    "relative admissibility tolerance + exact principal check",
)

# `scale_ref` is now computed above; drop the local duplicate in the tangent block.
edit(
    "ASDPlasticMaterial3D.h",
    """            double sref = yf.strength_scale(iv_storage, parameters_storage);
            if (sref < 0) sref = -sref;
            const double eps_deg = 1e-9 * ((xmax > sref) ? xmax : sref);""",
    """            const double eps_deg = 1e-9 * ((xmax > scale_ref) ? xmax : scale_ref);""",
    "reuse scale_ref for the degeneracy threshold",
)


def main():
    applied, skipped = [], []
    for relpath, anchor, new, tag in EDITS:
        path = os.path.join(SRC, relpath)
        if not os.path.isfile(path):
            print("MISSING FILE: %s" % path)
            return 2
        with open(path, "r", encoding="utf-8", newline="") as fh:
            txt = fh.read()
        crlf = "\r\n" in txt
        an = anchor.replace("\n", "\r\n") if crlf else anchor
        nw = new.replace("\n", "\r\n") if crlf else new
        if nw in txt:
            skipped.append(tag)
            continue
        n = txt.count(an)
        if n != 1:
            print("ANCHOR MISS (%d matches) in %s for %r" % (n, relpath, tag))
            return 3
        txt = txt.replace(an, nw, 1)
        with open(path, "w", encoding="utf-8", newline="") as fh:
            fh.write(txt)
        applied.append(tag)
    print("applied  : %d" % len(applied))
    for t in applied:
        print("   + %s" % t)
    print("already  : %d" % len(skipped))
    for t in skipped:
        print("   = %s" % t)
    return 0


if __name__ == "__main__":
    sys.exit(main())
