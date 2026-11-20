"""WP-238: the initial void ratio of LadrunoSANISAND set PER ELEMENT, so a deck can apply a
void-ratio field (Ladruno_scripts/sanisand_void_ratio_field.py) to a mesh that uses ONE
material tag:

    setParameter -val $e -ele $eleTag voidRatio

What is pinned (Ladruno_implementation/238_sanisand_void_ratio_field.md):
  (a)  a point given the void ratio e before the first step IS a point constructed with e:
       the stress and alpha of every step, elastic stage and plastic stage, bit for bit --
       by the untagged element form, the tagged element form and the `parameter` object;
  (b)  the untagged form was a SILENT NO-OP before this WP: the run it leaves is now the
       constructed one and NOT the default one (the gate fails on a build without the WP);
  (c)  only the elements named change: in one model with one material tag, three
       disconnected bricks -- one left alone, two set by `vf.apply` -- each reproduce,
       bit for bit, the single-brick run constructed with its own void ratio;
  (d)  the initial elastic operator follows: a column with FREE equations solved in one
       Linear step moves exactly as the column constructed with that void ratio;
  (e)  refused, and the state left as it was: a value that is not a void ratio (0,
       negative, NaN), and a point with strain history in the plastic stage;
  (f)  allowed in the elastic stage with strain history (vanilla's meaning: the CURRENT
       void ratio becomes the value), read back through `psi`;
  (g)  a tagged call naming ANOTHER material's tag does not touch the point;
  (h)  plane strain: `quad` forwards the parameter and the point is the constructed one.

Runtime: a few seconds.
"""
import math

import pytest

from _testbed import ops
import sanisand_void_ratio_field as vf
import wp152_cutoff_tools as T

E_DEF = T.PARAMS[2]
E_NEW = 0.70
OPTS = T.BASE + T.R1
INCS = [[1.0e-5, -0.5e-5, -0.5e-5, 0, 0, 0]] * 30 + [[-1.0e-5, 0.5e-5, 0.5e-5, 0, 0, 0]] * 10
XY = [(0., 0.), (1., 0.), (1., 1.), (0., 1.)]


def _hex(v):
    return [float(x).hex() for x in v]


def _build(n_bricks, e_ctor=E_DEF, incs=INCS):
    """n disconnected unit bricks (2 m apart), ONE material tag, the same prescribed
    normal-strain path on each (T.build's recipe: zero free equations)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    P = list(T.PARAMS)
    P[2] = e_ctor
    ops.nDMaterial("LadrunoSANISAND", 1, *P, *OPTS)
    for b in range(n_bricks):
        n0 = 8 * b
        for k in range(2):
            for j, (x, y) in enumerate(XY):
                ops.node(n0 + 4 * k + j + 1, x + 2.0 * b, y, float(k))
                ops.fix(n0 + 4 * k + j + 1, 1 if x == 0. else 0, 1 if y == 0. else 0, 1 if k == 0 else 0)
        ops.element("stdBrick", b + 1, *[n0 + i for i in range(1, 9)], 1)
    T._series(1, incs)
    for i, ts in enumerate((1, 2, 3)):
        ops.pattern("Plain", 10 + i, ts)
        for b in range(n_bricks):
            n0 = 8 * b
            for k in range(2):
                for j, (x, y) in enumerate(XY):
                    n = n0 + 4 * k + j + 1
                    if i == 0 and x == 1.:
                        ops.sp(n, 1, 1.0)
                    if i == 1 and y == 1.:
                        ops.sp(n, 2, 1.0)
                    if i == 2 and k == 1:
                        ops.sp(n, 3, 1.0)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-12, 10, 0)
    ops.algorithm("Linear")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")


def _row(ele):
    return _hex(ops.eleResponse(ele, "material", 1, "stress")) + _hex(ops.eleResponse(ele, "material", 1, "alpha"))


def _drive(n_bricks, incs=INCS, mid=None):
    """Stage 0, the flip, the plastic path. Returns per element the rows of every step.
    `mid(k)` runs before plastic step k."""
    rows = {e: [] for e in range(1, n_bricks + 1)}
    ops.updateMaterialStage("-material", 1, "-stage", 0)
    for _ in range(T.N0):
        assert ops.analyze(1) == 0
        for e in rows:
            rows[e].append(_row(e))
    ops.updateMaterialStage("-material", 1, "-stage", 1)
    for k in range(len(incs)):
        if mid:
            mid(k)
        assert ops.analyze(1) == 0, k
        for e in rows:
            rows[e].append(_row(e))
    return rows


def _ctor_run(e):
    _build(1, e)
    r = _drive(1)[1]
    ops.wipe()
    return r


@pytest.fixture(scope="module")
def refs():
    return {e: _ctor_run(e) for e in (E_DEF, E_NEW, 0.66, 0.74)}


def _set(form, v, ele=1):
    if form == "ele":
        return ops.setParameter("-val", v, "-ele", ele, "voidRatio")
    if form == "ele_tag":
        return ops.setParameter("-val", v, "-ele", ele, "voidRatio", "1")
    ops.parameter(900 + ele, "element", ele, "voidRatio", "1")
    return ops.updateParameter(900 + ele, v)


# ----------------------------------------------------------------------- (a) (b)
def test_the_references_differ(refs):
    assert refs[E_NEW][-1] != refs[E_DEF][-1] and refs[E_NEW][T.N0 - 1] != refs[E_DEF][T.N0 - 1]


@pytest.mark.parametrize("form", ["ele", "ele_tag", "parameter"])
def test_set_before_the_first_step_is_a_point_constructed_with_it(refs, form):
    _build(1)
    _set(form, E_NEW)
    got = _drive(1)[1]
    ops.wipe()
    assert got != refs[E_DEF], "the parameter did nothing (the silent no-op of the untagged form)"
    assert got == refs[E_NEW]


# ----------------------------------------------------------------------- (c)
def test_only_the_named_elements_change_and_apply_sets_a_field(refs):
    _build(3)
    assert vf.apply(ops, [2, 3], [0.66, 0.74]) == 2
    got = _drive(3)
    ops.wipe()
    assert got[1] == refs[E_DEF]
    assert got[2] == refs[0.66]
    assert got[3] == refs[0.74]


def test_every_gauss_point_of_the_element_is_set():
    _build(1)
    _set("ele", E_NEW)
    ops.updateMaterialStage("-material", 1, "-stage", 0)
    assert ops.analyze(1) == 0
    psi = [ops.eleResponse(1, "material", gp, "psi")[0] for gp in range(1, 9)]
    ops.wipe()
    _build(1, E_NEW)
    ops.updateMaterialStage("-material", 1, "-stage", 0)
    assert ops.analyze(1) == 0
    want = [ops.eleResponse(1, "material", gp, "psi")[0] for gp in range(1, 9)]
    ops.wipe()
    assert _hex(psi) == _hex(want)


# ----------------------------------------------------------------------- (d)
def _column_disp(e_ctor, e_set=None):
    top = [-1.0e-6 * (k + 1) for k in range(3)]
    P0 = list(T.PARAMS)
    T.PARAMS[2] = e_ctor
    try:
        T.build_column(top, OPTS)
    finally:
        T.PARAMS[:] = P0
    if e_set is not None:
        vf.apply(ops, [1, 2], [e_set, e_set])
    ops.algorithm("Linear")
    ops.updateMaterialStage("-material", 1, "-stage", 0)
    assert ops.analyze(1) == 0
    u = _hex([ops.nodeDisp(n, 3) for n in (5, 6, 7, 8)])
    ops.wipe()
    return u


def test_the_initial_elastic_operator_follows():
    """Free equations, ONE Linear step: the middle nodes move by K^-1 f with K assembled
    from the initial operator mCe, which depends on the void ratio through G."""
    ctor, setp, default = _column_disp(E_NEW), _column_disp(E_DEF, E_NEW), _column_disp(E_DEF)
    assert ctor != default
    assert setp == ctor


# ----------------------------------------------------------------------- (e)
@pytest.mark.parametrize("bad", [0.0, -0.3, float("nan"), float("inf")])
def test_a_value_that_is_not_a_void_ratio_is_refused(refs, bad):
    _build(1)
    try:
        _set("ele", bad)
    except Exception:                                    # noqa: BLE001  (the interpreter may raise)
        pass
    got = _drive(1)[1]
    ops.wipe()
    assert got == refs[E_DEF]


def test_refused_in_the_plastic_stage_once_the_point_has_strain(refs, capfd):
    def mid(k):
        if k == 12:
            try:
                _set("ele", E_NEW)
            except Exception:                            # noqa: BLE001
                pass
    _build(1)
    capfd.readouterr()
    got = _drive(1, mid=mid)[1]
    err = capfd.readouterr().err
    ops.wipe()
    assert got == refs[E_DEF]
    assert "voidRatio refused" in err


def test_allowed_in_the_plastic_stage_while_the_point_is_strain_free():
    """A deck that sets the field with the material already at stage 1 (no strain yet).
    The reference runs the SAME stage sequence on a point constructed with the value."""
    def run(e_ctor, e_set):
        _build(1, e_ctor)
        ops.updateMaterialStage("-material", 1, "-stage", 1)
        if e_set is not None:
            _set("ele", e_set)
        rows = _drive(1)[1]
        ops.wipe()
        return rows
    ctor, setp, default = run(E_NEW, None), run(E_DEF, E_NEW), run(E_DEF, None)
    assert ctor != default
    assert setp == ctor


# ----------------------------------------------------------------------- (f)
def test_elastic_stage_with_strain_sets_the_current_void_ratio():
    _build(1)
    ops.updateMaterialStage("-material", 1, "-stage", 0)
    for _ in range(3):
        assert ops.analyze(1) == 0
    psi0 = ops.eleResponse(1, "material", 1, "psi")[0]
    ev = T.EV0 * 3.0 / T.N0                              # compression positive
    e_now = E_DEF - (1.0 + E_DEF) * ev
    _set("ele", E_NEW)
    psi1 = ops.eleResponse(1, "material", 1, "psi")[0]   # the stress has not moved
    ops.wipe()
    assert psi1 - psi0 == pytest.approx(E_NEW - e_now, abs=1e-12)
    assert math.isfinite(psi1)


# ----------------------------------------------------------------------- (g)
def test_a_tagged_call_for_another_material_does_not_touch_the_point(refs):
    _build(1)
    try:
        ops.setParameter("-val", E_NEW, "-ele", 1, "voidRatio", "2")
    except Exception:                                    # noqa: BLE001
        pass
    got = _drive(1)[1]
    ops.wipe()
    assert got == refs[E_DEF]


# ----------------------------------------------------------------------- (h)
def test_plane_strain_quad_forwards_the_parameter():
    """The footing decks are plane strain: `quad` with LadrunoSANISANDPlaneStrain must hand
    the untagged parameter to its Gauss points, and the point must be the constructed one."""
    hist, cur = [], [0.0, 0.0]
    for k in range(T.N0):
        hist.append([-T.EV0 / 2.0 * (k + 1) / T.N0] * 2)
    cur = list(hist[-1])
    for de in [[1.0e-5, -1.0e-5]] * 25:
        cur = [cur[i] - de[i] for i in range(2)]
        hist.append(list(cur))

    def run(e_ctor, e_set):
        P0 = list(T.PARAMS)
        T.PARAMS[2] = e_ctor
        try:
            T.build_plane(hist, OPTS)
        finally:
            T.PARAMS[:] = P0
        if e_set is not None:
            vf.apply(ops, [1], [e_set])
        rows = []
        ops.updateMaterialStage("-material", 1, "-stage", 0)
        for k in range(len(hist)):
            if k == T.N0:
                ops.updateMaterialStage("-material", 1, "-stage", 1)
            assert ops.analyze(1) == 0, k
            rows.append([_hex(ops.eleResponse(1, "material", gp, "stress")) for gp in range(1, 5)])
        ops.wipe()
        return rows
    ctor, setp, default = run(E_NEW, None), run(E_DEF, E_NEW), run(E_DEF, None)
    assert ctor != default
    assert setp == ctor
