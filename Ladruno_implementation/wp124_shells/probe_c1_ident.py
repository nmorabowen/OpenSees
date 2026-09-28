"""WP-124 C1 evidence probe -- meaningful ONLY in a CONTINUUM=IDENT mutant build (the
mutation_rows.py rows C1a/C1b set it). Run by mutation_rows.py through run_pytest.py.

Under IDENT the ADR-87 gate replaces the cached *Ki of LadrunoBrick with the identity. If
getInitialStiff hands the caller *Ki, ModifiedNewton -initial iterates against I and cannot
converge in a few iterations. If it hands out the class-static scratch (gap C1, std/bbar
branch before the fix), the caller gets the REAL K0 and converges at once: the mutation never
reaches the solver.
"""
import opensees as ops

HEX = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0),
       (0.0, 0.0, 1.0), (1.0, 0.0, 1.0), (1.1, 1.05, 1.1), (0.0, 1.0, 1.0)]


def test_ident_mutation_reaches_the_caller_of_getInitialStiff():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for i, xyz in enumerate(HEX):
        ops.node(i + 1, *xyz)
    for n in (1, 2, 3, 4):
        ops.fix(n, 1, 1, 1)
    ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.3)
    ops.element("LadrunoBrick", 1, *range(1, 9), 1)          # std: the C1 branch
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(7, 1.0, 0.3, 0.0)
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-12, 5, 0)
    ops.algorithm("ModifiedNewton", "-initial")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static", "-noWarnings")
    rc = ops.analyze(1)
    assert rc != 0, ("ModifiedNewton -initial converged in <= 5 iterations: the caller got the "
                     f"REAL K0, not the IDENT-mutated *Ki (iterations {ops.testIter()})")
