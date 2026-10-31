"""WP-202: responseSpectrumAnalysis -scale was parsed and silently ignored.

Stock ResponseSpectrumAnalysis stored the -scale factor in m_scale but never
applied it, so a script that scaled a unit spectrum (Sa in g times 9.81, an
R/Ie reduction, an SLE scale) got the UNSCALED displacements with no warning.

Oracle: the same model run without -scale. With -scale s every modal
displacement must be exactly s times the unscaled one, on every path that
reads the spectrum: per-mode (-mode k), the full per-mode loop (last mode
left committed), -combine, and the timeSeries form.
"""
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

TN = [0.02, 0.1, 0.3, 0.6, 1.0, 3.0]
SA = [0.40, 1.00, 1.20, 0.90, 0.60, 0.20]
NODES = (2, 3, 4)


def _setup():
    ops.wipe()
    ops.model('basic', '-ndm', 1, '-ndf', 1)
    ops.node(1, 0.0)
    ops.fix(1, 1)
    for i, mm in enumerate((2.0, 1.5, 1.0), start=2):
        ops.node(i, float(i - 1))
        ops.mass(i, mm)
    for e, kk in enumerate((1200.0, 800.0, 500.0), start=1):
        ops.uniaxialMaterial('Elastic', e, kk)
        ops.element('Truss', e, e, e + 1, 1.0, e)
    ops.eigen('-fullGenLapack', 3)
    ops.modalProperties()


def _field():
    return [ops.nodeDisp(n, 1) for n in NODES]


def _assert_scaled(base, scaled, s):
    assert any(abs(b) > 1e-12 for b in base), "reference field is zero"
    for b, v in zip(base, scaled):
        assert v == pytest.approx(s * b, rel=1e-12, abs=1e-15), (base, scaled, s)


@pytest.mark.parametrize("mode", [1, 2, 3])
def test_scale_per_mode(mode):
    _setup()
    ops.responseSpectrumAnalysis(1, '-Tn', *TN, '-Sa', *SA, '-mode', mode)
    base = _field()
    ops.responseSpectrumAnalysis(1, '-Tn', *TN, '-Sa', *SA, '-mode', mode, '-scale', 9.81)
    _assert_scaled(base, _field(), 9.81)


def test_scale_full_loop_and_order_independent():
    # -scale before the lists must parse the same as after them
    _setup()
    ops.responseSpectrumAnalysis(1, '-Tn', *TN, '-Sa', *SA)
    base = _field()   # last mode committed
    ops.responseSpectrumAnalysis(1, '-scale', 0.25, '-Tn', *TN, '-Sa', *SA)
    _assert_scaled(base, _field(), 0.25)


@pytest.mark.parametrize("rule", ["SRSS", "CQC"])
def test_scale_combined(rule):
    _setup()
    extra = ['-damp', 0.05] if rule == "CQC" else []
    ops.responseSpectrumAnalysis(1, '-Tn', *TN, '-Sa', *SA, '-combine', rule, *extra)
    base = _field()
    ops.responseSpectrumAnalysis(1, '-Tn', *TN, '-Sa', *SA, '-combine', rule, *extra,
                                 '-scale', 2.5)
    _assert_scaled(base, _field(), 2.5)


def test_scale_timeseries_form():
    _setup()
    ops.timeSeries('Path', 7, '-time', *TN, '-values', *SA)
    ops.responseSpectrumAnalysis(7, 1, '-mode', 1)
    base = _field()
    ops.responseSpectrumAnalysis(7, 1, '-mode', 1, '-scale', 3.0)
    _assert_scaled(base, _field(), 3.0)


def test_scale_one_is_default():
    _setup()
    ops.responseSpectrumAnalysis(1, '-Tn', *TN, '-Sa', *SA, '-mode', 2)
    base = _field()
    ops.responseSpectrumAnalysis(1, '-Tn', *TN, '-Sa', *SA, '-mode', 2, '-scale', 1.0)
    assert _field() == base
