"""Cross-examination: LadrunoUniaxialJ2 (+Lemaitre) vs ASDSteel1D (Zone-A).

ASDSteel1D (Casalucci, Petracca, Camata — ASDEA) is a peer-reviewed plastic-damage
steel-bar model: J2 + a 2-term Armstrong-Frederick (Chaboche) kinematic backbone,
ductile fracture, buckling, bond-slip, IMPL-EX, char-length regularization. It shares
the *backbone family* with LadrunoUniaxialJ2 and is an independent second
implementation of the same ductile-steel physics — an ideal internal cross-check for
the Lemaitre damage mode.

Two legs, matched to what is and isn't shared:

  D1  BACKBONE (quantitative).  ASDSteel1D calibrates its two backstresses internally
      from (E, sy, su, eu):  H1 = E/(400 eu)*0.9, g1 = E/(400 eu (su-sy)),
      H2 = E/(400 eu)*5, g2 = g1*50, saturation H_k/g_k.  LadrunoUniaxialJ2 uses the
      SAME (C,g) Chaboche convention, so we replicate those exact backstresses in
      `-kin 2` and compare the undamaged tension backbone.  The only difference is the
      AF time-integration (ASDSteel1D integrates the AF ODE exactly-exponentially;
      Ladruno uses backward-Euler), so the two agree to ~1% and Ladruno CONVERGES to
      ASDSteel1D under step refinement.

  D2  DUCTILE FRACTURE (qualitative).  The damage *laws* differ (energy/triaxiality
      Lemaitre vs the ASDEA strain/energy fracture), so this is a shape check, not an
      overlay: on the same coupon BOTH reach a hardening peak, then SOFTEN to a small
      fraction of the peak with a MONOTONE damage variable approaching rupture.

  D3  IMPL-EX completion.  Both ship an IMPL-EX path for the softening-robustness
      escape; both must complete a full push-to-rupture with `-implex` and damage.

Pure OpenSees (both materials are built in) => Zone-A.
"""
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

# steel coupon (N, mm): mild structural steel
_E, _SY, _SU, _EU = 200000.0, 400.0, 550.0, 0.10


def _asd_chaboche(E, sy, su, eu):
    """The exact internal Chaboche calibration ASDSteel1D derives from (E,sy,su,eu)
    (ASDSteel1DMaterial.cpp OPS_ASDSteel1DMaterial)."""
    H1b = E / (400.0 * eu)
    g1 = E / (400.0 * eu * (su - sy))
    return (H1b * 0.9, g1, H1b * 5.0, g1 * 50.0)   # (C1,g1,C2,g2)


def _truss(mat_fn, hist, dmg_key=None):
    """Drive a strain history; return list of (eps, sigma, D|None). dmg_key reads the
    material's damage response under that name (e.g. 'damage' or 'Damage')."""
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(1, 0.0)
    ops.node(2, 1.0)
    mat_fn(1)
    ops.fix(1, 1)
    ops.element("Truss", 1, 1, 2, 1.0, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 1.0)
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.test("NormDispIncr", 1.0e-8, 200, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")
    out = []
    prev = 0.0
    for eps in hist:
        ops.integrator("DisplacementControl", 2, 1, eps - prev)
        if ops.analyze(1) != 0:
            break
        prev = eps
        s = ops.eleResponse(1, "material", 1, "stress")
        D = None
        if dmg_key:
            d = ops.eleResponse(1, "material", 1, dmg_key)
            D = d[0] if d else None
        out.append((eps, s[0] if s else 0.0, D))
    return out


def _ramp(eps_t, n):
    return [eps_t * i / n for i in range(1, n + 1)]


# ==========================================================================
# D1 — undamaged backbone overlay + convergence under refinement
# ==========================================================================
def _backbone_err(nsteps):
    C1, g1, C2, g2 = _asd_chaboche(_E, _SY, _SU, _EU)

    def asd(t):
        ops.uniaxialMaterial("ASDSteel1D", t, _E, _SY, _SU, _EU)   # no fracture => pure plasticity

    def lad(t):
        ops.uniaxialMaterial("LadrunoUniaxialJ2", t, _E,
                             "-iso", "voce", _SY, 0.0, 0.0, 0.0,    # constant yield = sy
                             "-kin", 2, C1, g1, C2, g2)             # ASDSteel1D backstresses

    hist = _ramp(0.08, nsteps)                                      # stay below eu (no fracture)
    a = _truss(asd, hist)
    l = _truss(lad, hist)
    worst = 0.0
    for (ea, sa, _), (el, sl, _) in zip(a, l):
        if abs(sa) < 1.0:                                           # skip the tiny elastic foot
            continue
        worst = max(worst, abs(sa - sl) / abs(sa))
    return worst


def test_D1_backbone_overlay_and_convergence():
    """Replicating ASDSteel1D's internal Chaboche in LadrunoUniaxialJ2 reproduces its
    undamaged steel backbone to ~1%, and the residual (the AF backward-Euler vs
    exact-exponential difference) SHRINKS under step refinement."""
    err_coarse = _backbone_err(20)
    err_fine = _backbone_err(160)
    assert err_coarse < 0.02, f"coarse backbone error {err_coarse:.2%} (> 2%)"
    assert err_fine < err_coarse, \
        f"not converging: fine {err_fine:.3%} !< coarse {err_coarse:.3%}"
    assert err_fine < 0.005, f"refined backbone error {err_fine:.3%} (> 0.5%)"


# ==========================================================================
# D2 — ductile-fracture signature parity (qualitative: peak -> soften, monotone D)
# ==========================================================================
def _fracture_signature(mat_fn, dmg_key, eps_t=0.20, n=160):
    res = _truss(mat_fn, _ramp(eps_t, n), dmg_key=dmg_key)
    sig = [s for (_e, s, _D) in res]
    Ds = [D for (_e, _s, D) in res if D is not None]
    peak = max(sig)
    ipk = sig.index(peak)
    tail = min(sig[ipk:])                                  # lowest stress after the peak
    mono = all(b >= a - 1e-9 for a, b in zip(Ds, Ds[1:])) if Ds else False
    return {"peak": peak, "tail": tail, "Dfinal": (Ds[-1] if Ds else 0.0),
            "mono": mono, "nsteps": len(res), "ipk": ipk}


def test_D2_ductile_fracture_signature_parity():
    """Both models, driven to rupture on the same coupon, must show the ductile-
    fracture signature: a hardening peak, post-peak softening to a small fraction of
    the peak, and a monotone damage variable approaching rupture. (The curves are NOT
    expected to overlay — the damage laws differ.)"""
    C1, g1, C2, g2 = _asd_chaboche(_E, _SY, _SU, _EU)

    def asd(t):
        ops.uniaxialMaterial("ASDSteel1D", t, _E, _SY, _SU, _EU, "-fracture")

    def lad(t):
        ops.uniaxialMaterial("LadrunoUniaxialJ2", t, _E,
                             "-iso", "voce", _SY, 0.0, 0.0, 0.0,
                             "-kin", 2, C1, g1, C2, g2,
                             "-damage", "lemaitre", 0.12, 1.0, 0.03, 0.95)

    A = _fracture_signature(asd, "Damage")
    L = _fracture_signature(lad, "damage")

    for tag, r in (("ASDSteel1D", A), ("Ladruno+Lemaitre", L)):
        assert r["peak"] > _SY, f"{tag}: never hardened past yield (peak {r['peak']:.0f})"
        assert r["ipk"] > 0, f"{tag}: peak at step 0 (no hardening branch)"
        assert r["tail"] < 0.3 * r["peak"], \
            f"{tag}: did not soften (tail {r['tail']:.0f} vs peak {r['peak']:.0f})"
        assert r["mono"], f"{tag}: damage not monotone"
        assert r["Dfinal"] > 0.8, f"{tag}: final damage {r['Dfinal']:.2f} (< 0.8)"


# ==========================================================================
# D2b — rough CALIBRATION of Lemaitre (p_D, r) to ASDSteel1D -fracture
# ==========================================================================
def _energy(res):
    """Area under the σ–ε curve (list of (eps, sigma, D|None))."""
    W = 0.0
    for (e0, s0, _), (e1, s1, _) in zip(res, res[1:]):
        W += 0.5 * (s0 + s1) * (e1 - e0)
    return W


def calibrate_lemaitre_to_asd(E, sy, su, eu, C1, g1, C2, g2, eps_t=0.20, n=400):
    """Fit Lemaitre (r,s,p_D,D_c) to ASDSteel1D `-fracture` on a monotone tension coupon:
    p_D := eupl = eu - sy/E (shared necking ONSET); r bisected to match dissipated ENERGY.
    The softening SHAPE still differs — that is the irreducible damage-law difference."""
    hist = _ramp(eps_t, n)
    asd = _truss(lambda t: ops.uniaxialMaterial("ASDSteel1D", t, E, sy, su, eu, "-fracture"),
                 hist, dmg_key="Damage")
    W_target = _energy(asd)
    pD = round(eu - sy / E, 3)
    s, Dc = 1.0, 0.95

    def W_of(r):
        l = _truss(lambda t: ops.uniaxialMaterial(
            "LadrunoUniaxialJ2", t, E, "-iso", "voce", sy, 0, 0, 0,
            "-kin", 2, C1, g1, C2, g2, "-damage", "lemaitre", r, s, pD, Dc), hist, dmg_key="damage")
        return _energy(l)

    lo, hi = 0.02, 0.6                        # larger r -> slower damage -> more energy
    for _ in range(22):
        r = 0.5 * (lo + hi)
        if W_of(r) < W_target:
            lo = r
        else:
            hi = r
    return (round(0.5 * (lo + hi), 4), s, pD, Dc)


def test_D2b_calibration_matches_onset_and_energy():
    """A two-parameter rough calibration (p_D=eupl, r energy-matched) makes the Lemaitre
    fracture coupon reproduce ASDSteel1D's PEAK and dissipated ENERGY — the residual is only
    the softening-shape difference of the two damage laws (README §5). Quantitative backing
    for the 'hysteresis aligns once calibrated' claim."""
    C1, g1, C2, g2 = _asd_chaboche(_E, _SY, _SU, _EU)
    r, s, pD, Dc = calibrate_lemaitre_to_asd(_E, _SY, _SU, _EU, C1, g1, C2, g2)
    assert abs(pD - (_EU - _SY / _E)) < 1e-6, f"p_D should equal eupl, got {pD}"

    hist = _ramp(0.20, 400)
    asd = _truss(lambda t: ops.uniaxialMaterial("ASDSteel1D", t, _E, _SY, _SU, _EU, "-fracture"),
                 hist, dmg_key="Damage")
    cal = _truss(lambda t: ops.uniaxialMaterial(
        "LadrunoUniaxialJ2", t, _E, "-iso", "voce", _SY, 0, 0, 0,
        "-kin", 2, C1, g1, C2, g2, "-damage", "lemaitre", r, s, pD, Dc), hist, dmg_key="damage")

    pk_a = max(x[1] for x in asd); pk_c = max(x[1] for x in cal)
    W_a, W_c = _energy(asd), _energy(cal)
    assert abs(pk_c - pk_a) / pk_a < 0.03, f"calibrated peak {pk_c:.0f} vs ASD {pk_a:.0f} (>3%)"
    assert abs(W_c - W_a) / W_a < 0.05, f"calibrated energy {W_c:.1f} vs ASD {W_a:.1f} (>5%)"
    # the uncalibrated default is visibly worse on the peak (it damages early) — proves the
    # calibration is doing real work, not a no-op.
    unc = _truss(lambda t: ops.uniaxialMaterial(
        "LadrunoUniaxialJ2", t, _E, "-iso", "voce", _SY, 0, 0, 0,
        "-kin", 2, C1, g1, C2, g2, "-damage", "lemaitre", 0.12, 1.0, 0.03, 0.95), hist, dmg_key="damage")
    pk_u = max(x[1] for x in unc)
    assert abs(pk_c - pk_a) < abs(pk_u - pk_a), \
        f"calibration did not improve the peak match (cal {pk_c:.0f}, unc {pk_u:.0f}, ASD {pk_a:.0f})"


# ==========================================================================
# D3 — IMPL-EX completes a softening push-to-rupture for both models
# ==========================================================================
def test_D3_implex_completion_parity():
    """Both models expose an IMPL-EX path for the indefinite softening tangent; with
    `-implex` both complete the full push and accumulate damage (no Newton stall)."""
    C1, g1, C2, g2 = _asd_chaboche(_E, _SY, _SU, _EU)
    n = 200

    def asd(t):
        ops.uniaxialMaterial("ASDSteel1D", t, _E, _SY, _SU, _EU, "-fracture", "-implex")

    def lad(t):
        ops.uniaxialMaterial("LadrunoUniaxialJ2", t, _E,
                             "-iso", "voce", _SY, 0.0, 0.0, 0.0,
                             "-kin", 2, C1, g1, C2, g2,
                             "-damage", "lemaitre", 0.12, 1.0, 0.03, 0.95, "-implex")

    a = _truss(asd, _ramp(0.20, n), dmg_key="Damage")
    l = _truss(lad, _ramp(0.20, n), dmg_key="damage")
    assert len(a) == n, f"ASDSteel1D IMPL-EX completed only {len(a)}/{n} steps"
    assert len(l) == n, f"Ladruno IMPL-EX completed only {len(l)}/{n} steps"
    assert a[-1][2] is not None and a[-1][2] > 0.5, "ASDSteel1D IMPL-EX: low final damage"
    assert l[-1][2] is not None and l[-1][2] > 0.5, "Ladruno IMPL-EX: low final damage"
