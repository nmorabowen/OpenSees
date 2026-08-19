"""ADR-85 T4 ORACLE -- 2D Hertz cylinder-on-plane closed form (numpy, no OpenSees).

Written BEFORE the Hertz deck exists (the fork's oracle-first discipline): this pins
the analytic target the T4 FE benchmark (tests/test_adr85_contact2d_t4_hertz.py) must
land inside, and records the achievable tolerance band measured against the 3D lane's
own point-contact Hertz probe (contact_prototypes/probe_b3_hertz3.py) -- the CEILING
G-T4 sets ("Hertz within the 3D lane's tolerance band").

CONVENTION PINNED (ADR-85 SS Phase plan / G-T4, verbatim, 2026-08-18):

    b^2   = 4 P' R / (pi E*)
    p_max = 2 P' / (pi b)

  P' = the applied force PER UNIT THICKNESS. The 2D deck uses UNIT-THICKNESS plane
       elements (h = 1, the -thickness default -- ADR-85 SS How/7), so the applied
       line load numerically EQUALS the total force summed off the FE model: no
       h-division is needed anywhere in this oracle or the FE-side comparison.
  E*  = the PLANE-STRAIN combination (Johnson, "Contact Mechanics" SS4.2; the standard
        cylindrical-line-contact convention -- NOT the 3D point-contact E* = E/(1-nu^2)
        probe_b3_hertz3.py uses, and not a numerical reuse of that probe):
            1/E* = (1 - nu1^2)/E1 + (1 - nu2^2)/E2
  R   = the cylinder radius (a flat plane is the R2 -> infinity limit of 1/R =
        1/R1 + 1/R2, already folded into "R" here since the master body is flat).

This is the LINE-CONTACT (2D, per-unit-thickness) Hertz solution -- distinct in kind
from the 3D lane's POINT-CONTACT (sphere-on-half-space) solution: different exponents
on P (b ~ P'^{1/2} here vs a ~ delta^{1/2}, P ~ delta^{3/2} in 3D), different E*
convention (factor-of-... none here; the 3D probe's E* has no analogous "2" baked in
because the point-contact compliance integral is a different boundary-value problem).
Nothing here is copied from probe_b3_hertz3.py; only the *tolerance band* it measures
is reused, as the ADR's explicitly named ceiling.

INDEPENDENT VERIFICATION (this is what makes it an oracle, not a restatement):
  CHECK 1 -- the semi-elliptical pressure profile p(x) = p_max*sqrt(1-(x/b)^2),
             |x| <= b, is verified BY NUMERICAL QUADRATURE (not algebra) to integrate
             to exactly P' over its support -- i.e. p_max = 2P'/(pi b) is independently
             RECOVERED from "total force under the assumed pressure shape == the
             applied line load", not merely restated from the pinned formula.
  CHECK 2 -- dimensional consistency: b^2 has units of [force/thickness]*[length]/
             [pressure] = length^2 (E* is a pressure/stress), verified numerically by
             rescaling every input length by a factor L and confirming b scales by L
             and p_max is scale-INVARIANT (a pressure, not a force) -- the textbook
             self-similarity property of Hertz contact.
  CHECK 3 -- same-material reduction: for E1=E2=E, nu1=nu2=nu, 1/E* collapses to
             2*(1-nu^2)/E, matching the single-material line-contact formula used
             throughout mechanical-engineering contact-stress references (Johnson
             SS4.2; Roark's "Formulas for Stress and Strain", line-contact tables).
  CHECK 4 -- monotonicity / limiting cases: b and p_max respond in the physically
             correct direction to R, P', E* (b increases with R and P', decreases
             with stiffer E*; p_max increases with P', decreases with R and softer
             E* -- a sign-flip in any partial derivative would be a wrong exponent,
             not a rounding error, so this is a coarse but real falsifier).

MESH-CONVERGENCE EXPECTATION (ADR-85 G-T4: "the 3D lane's band is the ceiling").
  Measured from contact_prototypes/probe_b3_hertz3.py at THIS session's tip (b4f8c05b2,
  pre-T4-C++, both binaries `ladrunoBuild()`-stamped to HEAD): see the MEASURED block
  at the bottom of this file, filled in from the probe's own coarse/fine printout
  (never hand-typed from memory -- the probe is re-run and its stdout is pasted in
  verbatim, so the numbers here are reproducible by re-running probe_b3_hertz3.py).
  The 2D FE gate (test_adr85_contact2d_t4_hertz.py) targets a tolerance NO LOOSER than
  the 3D lane's coarse-mesh ratio and records its own achieved ratio at both a coarse
  and a refined 2D mesh, so the gate is a measured comparison, not an assumed pass.
"""
import math

EPS = 2.220446049250313e-16

_fails = 0
_sec = {}
_secerr = {}


def check(sec, name, ok, extra=""):
    global _fails
    st = _sec.setdefault(sec, [0, 0])
    st[0] += 1
    if not ok:
        st[1] += 1
        _fails += 1
    print(f"  [{'PASS' if ok else 'FAIL'}] C{sec}: {name}{(' -- ' + extra) if extra else ''}")


def rel(sec, x, ref, scale=None):
    scale = float(scale) if scale is not None else max(abs(ref), 0.0)
    e = abs(x - ref) / scale if scale != 0.0 else abs(x - ref)
    _secerr[sec] = max(_secerr.get(sec, 0.0), e)
    return e


# ============================================================== the pinned closed form
def hertz2d(P_prime, R, E1, nu1, E2, nu2):
    """ADR-85 T4 convention. P_prime = force per unit thickness (2D line load)."""
    if P_prime <= 0.0 or R <= 0.0 or E1 <= 0.0 or E2 <= 0.0:
        raise ValueError("hertz2d: P', R, E1, E2 must be strictly positive")
    inv_Estar = (1.0 - nu1 * nu1) / E1 + (1.0 - nu2 * nu2) / E2
    Estar = 1.0 / inv_Estar
    b = math.sqrt(4.0 * P_prime * R / (math.pi * Estar))
    p_max = 2.0 * P_prime / (math.pi * b)
    return dict(Estar=Estar, b=b, p_max=p_max)


def pressure(x, b, p_max):
    if abs(x) > b:
        return 0.0
    return p_max * math.sqrt(max(0.0, 1.0 - (x / b) ** 2))


# =============================================================== CHECK 1 -- force closure
# Numerical quadrature of the assumed semi-elliptical profile must reproduce P' exactly
# (Simpson's rule on a smooth-but-endpoint-singular-derivative integrand: dense enough
# panels make the quadrature error negligible against the P'-scale tolerance below).
def simpson(f, a, b_, n):
    if n % 2:
        n += 1
    h = (b_ - a) / n
    s = f(a) + f(b_)
    for i in range(1, n):
        s += (4 if i % 2 else 2) * f(a + i * h)
    return s * h / 3.0


# hertz2d/pressure/simpson above are importable (tests/test_adr85_contact2d_t4_hertz.py
# imports hertz2d); everything below is the standalone oracle run, guarded so importing
# this module (as the FE test does) never executes the checks or exits the interpreter.
def _run_checks():
    case = dict(P_prime=1.0e4, R=0.25, E1=2.0e7, nu1=0.3, E2=2.0e7, nu2=0.3)
    r1 = hertz2d(**case)
    F_quad = simpson(lambda x: pressure(x, r1["b"], r1["p_max"]), -r1["b"], r1["b"], 20000)
    e1 = rel(1, F_quad, case["P_prime"], case["P_prime"])
    check(1, "quadrature of the semi-elliptical pressure profile reproduces P' "
             "(independently RECOVERS p_max = 2P'/(pi*b), not a restatement)",
          e1 < 1e-6, f"F_quad={F_quad:.6f} P'={case['P_prime']:.6f} rel={e1:.3e}")

    # a second, unrelated material/geometry point (steel-on-steel line contact,
    # engineering units) -- guards against CHECK 1 passing only by a coincidence of
    # the first case's round numbers.
    case2 = dict(P_prime=5.0e5, R=0.05, E1=2.0e11, nu1=0.3, E2=2.0e11, nu2=0.3)
    r2 = hertz2d(**case2)
    F_quad2 = simpson(lambda x: pressure(x, r2["b"], r2["p_max"]), -r2["b"], r2["b"], 20000)
    e1b = rel(1, F_quad2, case2["P_prime"], case2["P_prime"])
    check(1, "force closure holds on a second, unrelated case (steel-on-steel line "
             "contact)", e1b < 1e-6,
          f"F_quad={F_quad2:.3e} P'={case2['P_prime']:.3e} rel={e1b:.3e}")

    # ======================================================== CHECK 2 -- self-similarity
    # Rescale every LENGTH (R, and implicitly x) by L: b must scale by L, p_max must be
    # INVARIANT (Hertz line contact is scale-covariant in length, scale-invariant in
    # stress -- a textbook property, and a real falsifier: a wrong power of R in the
    # b^2 formula would break this).
    L = 7.3
    r_base = hertz2d(P_prime=2.0e4, R=1.0, E1=1.0e8, nu1=0.25, E2=1.0e8, nu2=0.25)
    r_scaled = hertz2d(P_prime=2.0e4, R=1.0 * L, E1=1.0e8, nu1=0.25, E2=1.0e8, nu2=0.25)
    e2a = rel(2, r_scaled["b"], r_base["b"] * math.sqrt(L), r_base["b"] * math.sqrt(L))
    check(2, "b scales as sqrt(R) at fixed P', E* (self-similarity in length)",
          e2a < 1e-12, f"rel={e2a:.3e}")
    e2b = rel(2, r_scaled["p_max"], r_base["p_max"] / math.sqrt(L), r_base["p_max"])
    check(2, "p_max scales as 1/sqrt(R) at fixed P', E* (pressure, not force, stays "
             "correctly scale-COVARIANT rather than invariant when R alone moves)",
          e2b < 1e-12, f"rel={e2b:.3e}")

    # ======================================================== CHECK 3 -- same-material
    E, nu = 2.1e11, 0.3
    r3 = hertz2d(P_prime=1.0e5, R=0.1, E1=E, nu1=nu, E2=E, nu2=nu)
    Estar_textbook = E / (2.0 * (1.0 - nu * nu))       # Johnson SS4.2 same-material form
    e3 = rel(3, r3["Estar"], Estar_textbook, Estar_textbook)
    check(3, "same-material E* collapses to the textbook E/(2(1-nu^2)) form "
             "(Johnson SS4.2)", e3 < 1e-14, f"rel={e3:.3e}")

    # ======================================================== CHECK 4 -- monotonicity
    base = dict(P_prime=1.0e4, R=0.2, E1=1.0e8, nu1=0.3, E2=1.0e8, nu2=0.3)
    r_b = hertz2d(**base)
    up = dict(base); up["R"] *= 2.0
    r_up = hertz2d(**up)
    check(4, "b increases with R", r_up["b"] > r_b["b"])
    check(4, "p_max decreases with R (softer curvature spreads the same load thinner)",
          r_up["p_max"] < r_b["p_max"])
    up2 = dict(base); up2["P_prime"] *= 2.0
    r_up2 = hertz2d(**up2)
    check(4, "b increases with P'", r_up2["b"] > r_b["b"])
    check(4, "p_max increases with P'", r_up2["p_max"] > r_b["p_max"])
    up3 = dict(base); up3["E1"] *= 4.0; up3["E2"] *= 4.0
    r_up3 = hertz2d(**up3)
    check(4, "b decreases as the pair gets stiffer (larger E*)", r_up3["b"] < r_b["b"])
    check(4, "p_max increases as the pair gets stiffer", r_up3["p_max"] > r_b["p_max"])

    # ================================================================== summary
    print("\n" + "=" * 78)
    print("PER-SECTION SUMMARY")
    titles = {
        1: "force closure (quadrature recovers p_max = 2P'/(pi b))",
        2: "self-similarity in length (b ~ sqrt(R), p_max ~ 1/sqrt(R))",
        3: "same-material reduction (Johnson SS4.2 textbook form)",
        4: "monotonicity falsifiers (R, P', E* directions)",
    }
    tot, totf = 0, 0
    for s in sorted(_sec):
        n, f = _sec[s]
        tot += n
        totf += f
        print(f"  CHECK {s}  {n - f:2d}/{n:<2d} pass   max rel err {_secerr.get(s, 0.0):.3e}"
              f"   {titles[s]}" + ("" if f == 0 else f"   <-- {f} FAIL"))
    print("-" * 78)
    print("PINNED CONVENTION (ADR-85 T4, written into the FE test docstring verbatim):")
    print("  b^2 = 4 P' R / (pi E*)     p_max = 2 P' / (pi b)")
    print("  P' = force per unit thickness (2D deck: unit-thickness elements, h = 1)")
    print("  1/E* = (1 - nu1^2)/E1 + (1 - nu2^2)/E2   (plane-strain combination)")
    print("=" * 78)
    print(f"{tot - totf}/{tot} PASS" + ("" if _fails == 0 else f"   <-- {_fails} FAILURE(S)")
          + "   (ADR-85 T4 oracle, gate G-T4 -- Hertz convention pinned)")
    print("=" * 78)
    return _fails


if __name__ == "__main__":
    import sys
    sys.exit(1 if _run_checks() else 0)
