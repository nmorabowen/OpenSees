"""WP-214: mode-kind tags in `modalProperties -print`, and the tighter Penalty
state-check margin in `residualModes`.

ADR: Ladruno_implementation/213_ladruno_modal_basis_completeness_adr.md, D6.

Tags. When the WP-213 mode-kind registry holds a valid record for the domain,
every mode row of the report (sections 2 and 6-10, screen and `-file`) ends in
a trailing comment `  # R` (residual pseudo-mode) or `  # Z` (Ritz vector);
eigen rows get nothing. Section 2 adds one `# NOTE:` line per kind present.
The comment keeps the rows numeric for a `#`-aware reader (numpy.loadtxt).
`-return` gains a `modeKind` key (0 / 1 / 2 per mode) only in that case.
Without a record (plain `eigen`, or a new `eigen` after `residualModes`) the
report is byte-identical to bac089465: pinned against
tests/data/wp214_modalprops_baseline_v3.txt, captured from that build
(`python tests/test_wp214_mode_tags.py --capture`, run on a bac089465 build).
Line endings are compared after universal-newline translation (the `-file`
stream is text mode, so it writes CRLF on Windows on both builds).

Penalty margin (Opus review of #1023). The `residualModes` "model changed
after eigen" check now refuses above max(1e-4, 10 f) and warns above
max(1e-7, f), f = eps |phi|^T |K| |phi| / phi^T K phi. On the tied-floor
Penalty model (kip-in, 10 storeys, f1 = 4e-3, measured round-off 2.5e-4) a 20 %
stiffness loss after `eigen 1` is now refused (it only warned at 100 f).
"""
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from _testbed import ops

import test_wp213_residual_modes as t213

pytestmark = [pytest.mark.zone_a]

HERE = Path(__file__).resolve().parent
BASELINE = HERE / 'data' / 'wp214_modalprops_baseline_v3.txt'
MARK = '# MODAL ANALYSIS REPORT'
MODE_SECTIONS = ('* 2.', '* 6.', '* 7.', '* 8.', '* 9.', '* 10.')


def _report(tmp_path, steps):
    """Run `steps` in a child (opserr is the C stderr), then
    `modalProperties -print -file`; return (screen report, file report)."""
    fname = tmp_path / 'mp.txt'
    script = tmp_path / 'child.py'
    script.write_text('\n'.join([
        "import sys",
        f"sys.path.insert(0, {str(HERE)!r})",
        "import test_wp213_residual_modes as t",
        "from _testbed import ops",
        steps,
        "sys.stderr.flush()",
        f"ops.modalProperties('-print', '-file', {str(fname)!r})",
    ]) + '\n')
    env = {k: v for k, v in os.environ.items() if k != 'LADRUNO_OPENSEES_QUIET'}
    r = subprocess.run([sys.executable, str(script)], env=env, stdin=subprocess.DEVNULL,
                       capture_output=True, text=True, timeout=120)
    assert r.returncode == 0, r.stderr[-2000:]
    err = r.stderr
    i = err.index(MARK)
    file_txt = fname.read_text()
    return err[i:], file_txt


def _mode_rows(report):
    """{section: [(mode, line), ...]} for the mode rows of sections 2, 6-10."""
    rows, sec = {}, None
    for line in report.split('\n'):
        if line.startswith('* '):
            sec = next((s for s in MODE_SECTIONS if line.startswith(s + ' ')), None)
            if sec:
                rows[sec] = []
            continue
        if sec and line.startswith(' ') and line.split():
            tok = line.split()[0]
            if tok.isdigit():
                rows[sec].append((int(tok), line))
    return rows


V3_EIGEN = "t._v3_build(); ops.eigen('-fullGenLapack', 3)"


def test_plain_eigen_report_is_byte_identical(tmp_path):
    screen, file_txt = _report(tmp_path, V3_EIGEN)
    base = BASELINE.read_text()
    assert file_txt == base
    assert screen[:len(base)] == base
    assert '# NOTE' not in screen and '  # R' not in screen


def test_new_eigen_clears_the_tags_byte_identical(tmp_path):
    steps = V3_EIGEN + "; ops.residualModes('-dir', 1); ops.eigen('-fullGenLapack', 3)"
    screen, file_txt = _report(tmp_path, steps)
    base = BASELINE.read_text()
    assert file_txt == base
    assert screen[:len(base)] == base


def test_residual_rows_tagged_R(tmp_path):
    screen, file_txt = _report(tmp_path, V3_EIGEN + "; ops.residualModes('-dir', 1)")
    for rep in (screen, file_txt):
        rows = _mode_rows(rep)
        assert set(rows) == set(MODE_SECTIONS), rows.keys()
        for sec, lines in rows.items():
            assert [m for m, _ in lines] == [1, 2, 3, 4, 5], (sec, lines)
            for m, line in lines:
                if m >= 4:
                    assert line.endswith('  # R'), (sec, line)
                    # the row stays numeric for a '#'-aware reader
                    vals = line.split('#')[0].split()
                    assert all(np.isfinite(float(v)) for v in vals)
                else:
                    assert '#' not in line, (sec, line)
        assert rep.count('  # R\n') == 12
        assert ('# NOTE: modes 4-5 (R) are residual pseudo-modes (static correction); '
                'their periods are Rayleigh, not natural.\n') in rep
    # the eigen rows are those of the plain report, character for character
    base_rows = _mode_rows(BASELINE.read_text())
    tagged = _mode_rows(file_txt)
    assert tagged['* 2.'][:3] == base_rows['* 2.']
    assert screen[:len(file_txt)] == file_txt


def test_period_override_named_in_the_note(tmp_path):
    _, file_txt = _report(tmp_path, V3_EIGEN + "; ops.residualModes('-dir', 1, '-period', 0.05)")
    assert file_txt.count('  # R\n') == 12
    assert 'responseSpectrumAnalysis reads their Sa at T = 0.05' in file_txt


def test_single_residual_mode_note_is_singular(tmp_path):
    # 3 storeys, 1 eigenmode: span{r, K^-1 M r} adds two; with 2 eigenmodes
    # the third storey's residual adds one
    steps = ("t.MODEL_MK = ([1.0, 1.0, 1.0], [300.0, 300.0, 300.0]); t.N3 = 3; "
             "t._v3_build(); ops.eigen('-fullGenLapack', 2); ops.residualModes('-dir', 1)")
    _, file_txt = _report(tmp_path, steps)
    assert '# NOTE: mode 3 (R) is a residual pseudo-mode (static correction); ' \
           'its period is Rayleigh, not natural.\n' in file_txt
    assert file_txt.count('  # R\n') == 6


def test_return_dict_gains_modeKind_only_with_a_record():
    t213._v3_build()
    ops.eigen('-fullGenLapack', 3)
    plain = ops.modalProperties('-return')
    assert 'modeKind' not in plain
    ops.residualModes('-dir', 1)
    tagged = ops.modalProperties('-return')
    assert list(tagged['modeKind']) == [0, 0, 0, 1, 1]
    assert set(tagged) - set(plain) == {'modeKind'}
    ops.eigen('-fullGenLapack', 3)
    again = ops.modalProperties('-return')
    assert set(again) == set(plain)
    for key in plain:
        np.testing.assert_array_equal(again[key], plain[key])


# ---------------------------------------------------------------------------
# Penalty state-check margin
# ---------------------------------------------------------------------------
def _soften(factor, nf=10, k=500.0):
    eles = list(range(1, 2 * nf + 1))
    ops.setParameter('-val', factor * k, '-ele', *eles, 'E')


@pytest.mark.parametrize('solver', [('-fullGenLapack', 1), (1,)])
def test_penalty_changed_model_is_refused(solver):
    # 10 storeys kip-in, alpha = 1e14: f1 = 4e-3. A 20 % stiffness loss after
    # eigen moves mode 1's Rayleigh quotient by 0.2: below 100 f1 = 0.4 (old
    # margin: warned only), above 10 f1 = 0.04 (new margin: refused)
    t213._penalty_tied_floors(10, 1.0, 500.0, solver)
    _soften(0.8)
    with pytest.raises(Exception):
        ops.residualModes('-dir', 1)
    # re-run eigen on the current state: accepted (3 modes: with 1 mode the
    # kernel's fixed 1e-4 K-solve residual check is marginal under Penalty,
    # issue #1025)
    ops.eigen(*solver[:-1], 3)
    assert ops.residualModes('-dir', 1)


def test_penalty_small_change_warns(tmp_path):
    # a 2 % loss after eigen 1: rel 0.02, between f1 and 10 f1 (warned; it
    # passed silently at 10 f1). With more modes the higher ones (smaller f)
    # refuse it outright.
    script = tmp_path / 'child.py'
    script.write_text('\n'.join([
        "import sys",
        f"sys.path.insert(0, {str(HERE)!r})",
        "import test_wp213_residual_modes as t",
        "import test_wp214_mode_tags as t4",
        "from _testbed import ops",
        "t._penalty_tied_floors(10, 1.0, 500.0, ('-fullGenLapack', 1))",
        "t4._soften(0.98)",
        "try:",
        "    ops.residualModes('-dir', 1)",
        "except Exception:",
        "    pass   # the eigen-1 K-solve residual check is marginal, #1025",
    ]) + '\n')
    r = subprocess.run([sys.executable, str(script)], env=dict(os.environ),
                       stdin=subprocess.DEVNULL, capture_output=True, text=True, timeout=120)
    assert r.returncode == 0, r.stderr[-2000:]
    out = r.stdout + r.stderr
    assert 'residualModes: WARNING mode 1 deviates from the current tangent' in out, out[-2000:]
    assert 'model changed after eigen' not in out, out[-2000:]


@pytest.mark.parametrize('nf,m,k', [(10, 1.0, 500.0), (3, 0.5, 100.0), (10, 500.0, 1.0e5)])
@pytest.mark.parametrize('solver', [('-fullGenLapack', 3), (3,), ('-fullGenLapack', 2), (2,)])
def test_penalty_unchanged_model_still_accepted(nf, m, k, solver):
    t213._penalty_tied_floors(nf, m, k, solver)
    ops.residualModes('-dir', 1)
    assert ops.residualModes('-kinds')


if __name__ == '__main__' and '--capture' in sys.argv:
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        _, file_txt = _report(Path(d), V3_EIGEN)
    BASELINE.write_text(file_txt, newline='\n')
    print('wrote', BASELINE, len(file_txt), 'chars')
