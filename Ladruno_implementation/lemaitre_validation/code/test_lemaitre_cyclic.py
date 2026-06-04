"""Lemaitre coupled ductile damage — CYCLIC validation (Zone-A).

Monotonic tests (test_lemaitre_validation.py) barely exercise the Armstrong-Frederick
(Chaboche) kinematic hardening, and Lemaitre's signature regime is low-cycle fatigue
(LCF) — damage accruing with the accumulated plastic strain cycle after cycle. This file
adds the cyclic battery:

  E1  multi-cycle hysteresis oracle: LadrunoUniaxialJ2 -damage lemaitre vs the independent
      from-the-equations reference integrator over several symmetric strain cycles
      (Bauschinger + back-stress saturation + damage), to ~1e-8.
  E2  LCF damage accumulation: under fixed-amplitude symmetric cycles D grows monotonically
      cycle-by-cycle; the per-cycle increment is steady once the loop stabilises; and a
      larger strain amplitude reaches rupture in FEWER cycles (the Coffin-Manson trend).
  E3  ASDSteel1D cyclic cross-exam: matched J2+Chaboche backstresses ⇒ the undamaged
      hysteresis loops overlay and CONVERGE under step refinement; with damage on both
      models show degrading (softening) loops with monotone damage.

All uniaxial (Truss + uniaxial material) ⇒ Zone-A.
"""
import pytest

from _testbed import ops
from _testbed.lemaitre_ref import Lemaitre1D
import test_lemaitre_validation as V
import test_lemaitre_vs_asdsteel1d as X

pytestmark = [pytest.mark.zone_a]


def _cyclic(amp, ncyc, per_branch):
    """Symmetric triangular strain history 0→+amp→−amp→+amp… for ncyc full cycles."""
    seq = [0.0, amp]
    for _ in range(ncyc):
        seq += [-amp, amp]
    hist = []
    for a, b in zip(seq, seq[1:]):
        for i in range(1, per_branch + 1):
            hist.append(a + (b - a) * i / per_branch)
    return hist


# ==========================================================================
# E1 — multi-cycle hysteresis oracle (code vs reference integrator)
# ==========================================================================
def test_E1_cyclic_oracle():
    hist = _cyclic(0.03, 4, 40)
    code = V._run_truss(lambda t: V._mat1d(t, dmg=V._DMG), hist)
    ref = V._ref(V._DMG)
    maxes = maxd = 0.0
    grew = False
    for (eps, sig, D, _p) in code:
        rs, rD, _rp, _rse = ref.step(eps)
        maxes = max(maxes, abs(sig - rs) / (abs(rs) + 1.0))
        maxd = max(maxd, abs(D - rD))
        grew = grew or D > 1e-3
    assert grew, "damage never accumulated over the cyclic history"
    assert maxes < 1e-8, f"cyclic stress vs reference: max rel err {maxes:.2e}"
    assert maxd < 1e-8, f"cyclic damage vs reference: max abs err {maxd:.2e}"


# ==========================================================================
# E2 — LCF: damage accumulates per cycle, steady increment, amplitude → life
# ==========================================================================
def _run_to_damage(amp, Dstop=0.6, maxcyc=400, per_branch=20):
    """Cycle at fixed amplitude (0→+amp, then ∓amp peaks) until D >= Dstop. Returns
    (cycles_to_stop, [per-cycle ΔD]); reads D at the end of each full +amp→−amp→+amp cycle."""
    V._build_truss(lambda t: V._mat1d(t, dmg=[1.5, 1.0, 0.005, 0.999]))
    state = {"prev": 0.0, "ok": True}

    def ramp(a, b):
        for i in range(1, per_branch + 1):
            e = a + (b - a) * i / per_branch
            ops.integrator("DisplacementControl", 2, 1, e - state["prev"])
            if ops.analyze(1) != 0:
                state["ok"] = False
                return False
            state["prev"] = e
        return True

    per, Dprev = [], 0.0
    if not ramp(0.0, amp):                        # initial quarter
        return 0, per
    for c in range(1, maxcyc + 1):
        if not (ramp(amp, -amp) and ramp(-amp, amp)):   # one full cycle
            return c - 1, per
        D = ops.eleResponse(1, "material", 1, "damage")[0]
        per.append(D - Dprev); Dprev = D
        if D >= Dstop:
            return c, per
    return maxcyc, per


def test_E2_lcf_damage_per_cycle_monotone_and_steady():
    # amplitude well above elastic shakedown => sustained plasticity every cycle
    cyc, per = _run_to_damage(0.05, Dstop=0.5)
    assert cyc >= 5, f"reached D-stop in too few cycles to assess ({cyc})"
    # damage increment each cycle is positive (monotone D)
    assert all(d > 0 for d in per[1:]), f"per-cycle damage not all positive: {per}"
    # once the loop stabilises the per-cycle increment is ~steady (last few within 40%)
    tail = per[-4:]
    assert max(tail) <= 1.4 * min(tail), f"per-cycle ΔD not steady in tail: {tail}"


def test_E2_amplitude_shortens_life():
    # both amplitudes above shakedown (finite life); Coffin-Manson: larger amp => fewer cycles
    n_small, _ = _run_to_damage(0.04, Dstop=0.5)
    n_large, _ = _run_to_damage(0.08, Dstop=0.5)
    assert n_large < n_small, \
        f"larger amplitude should rupture sooner: N(0.08)={n_large} !< N(0.04)={n_small}"


# ==========================================================================
# E3 — ASDSteel1D cyclic cross-exam
# ==========================================================================
def test_E3_cyclic_backbone_overlay_and_convergence():
    """Matched J2+Chaboche backstresses ⇒ the undamaged cyclic hysteresis loops overlay,
    and the residual (AF backward-Euler vs ASDEA exact-exponential) shrinks under
    step refinement."""
    C1, g1, C2, g2 = X._asd_chaboche(X._E, X._SY, X._SU, X._EU)

    def asd(t):
        ops.uniaxialMaterial("ASDSteel1D", t, X._E, X._SY, X._SU, X._EU)

    def lad(t):
        ops.uniaxialMaterial("LadrunoUniaxialJ2", t, X._E,
                             "-iso", "voce", X._SY, 0.0, 0.0, 0.0, "-kin", 2, C1, g1, C2, g2)

    def err(per_branch):
        hist = _cyclic(0.02, 3, per_branch)
        a = X._truss(asd, hist)
        l = X._truss(lad, hist)
        return max(abs(sa - sl) / (abs(sa) + 1.0) for (_e, sa, _), (_e2, sl, _) in zip(a, l))

    e_coarse = err(15)
    e_fine = err(120)
    assert e_coarse < 0.08, f"cyclic backbone error {e_coarse:.2%} (>8%)"
    assert e_fine < e_coarse, f"not converging: fine {e_fine:.3%} !< coarse {e_coarse:.3%}"
    assert e_fine < 0.025, f"refined cyclic backbone error {e_fine:.3%} (>2.5%)"


def test_E3_cyclic_damage_parity():
    """With damage/fracture on, both models show degrading (softening) cyclic loops and a
    monotone damage variable — same LCF signature (laws differ ⇒ qualitative)."""
    C1, g1, C2, g2 = X._asd_chaboche(X._E, X._SY, X._SU, X._EU)

    def asd(t):
        ops.uniaxialMaterial("ASDSteel1D", t, X._E, X._SY, X._SU, X._EU, "-fracture")

    def lad(t):
        ops.uniaxialMaterial("LadrunoUniaxialJ2", t, X._E, "-iso", "voce", X._SY, 0.0, 0.0, 0.0,
                             "-kin", 2, C1, g1, C2, g2, "-damage", "lemaitre", 0.12, 1.0, 0.03, 0.95)
    hist = _cyclic(0.05, 6, 30)
    for name, mat, key in (("ASDSteel1D", asd, "Damage"), ("Ladruno", lad, "damage")):
        res = X._truss(mat, hist, dmg_key=key)
        Ds = [D for (_e, _s, D) in res if D is not None]
        assert Ds and Ds[-1] > 0.3, f"{name}: low final cyclic damage {Ds[-1] if Ds else None}"
        assert all(b >= a - 1e-9 for a, b in zip(Ds, Ds[1:])), f"{name}: cyclic damage not monotone"
        # peak stress in the last cycle is below the first (loops degrade)
        n = len(res)
        peak_first = max(abs(s) for (_e, s, _D) in res[:n // 6])
        peak_last = max(abs(s) for (_e, s, _D) in res[-n // 6:])
        assert peak_last < peak_first, f"{name}: loops did not degrade ({peak_last:.0f} !< {peak_first:.0f})"
