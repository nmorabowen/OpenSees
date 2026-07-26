# Regression: a first-pivot singularity must NOT report success.
#
# THE BUG (fixed by the accompanying change). BandGen/FullGen/BandSPD LAPACK
# solvers mapped a positive LAPACK `info` with `return -info+1;`, which C parses
# as (-info)+1. LAPACK sets info = i when U(i,i) is exactly zero, so info == 1 --
# a singularity at the FIRST pivot, which is exactly what an all-zero assembled A
# produces -- returned 0, i.e. SUCCESS.
#
# Why that is a silent wrong answer rather than a cosmetic return-code nit:
#   * the solvers copy B into X BEFORE calling LAPACK
#     (BandGenLinLapackSolver.cpp:116-118);
#   * DGBSV/DGESV/DPBSV do not compute X when info > 0, so X is left equal to B;
#   * the algorithm gates on `theSOE->solve() < 0` and therefore accepts the
#     LOAD VECTOR as a displacement increment;
#   * `theSOE->factored = true` is never reached on the error path, so the flag
#     stays false and every later iteration silently repeats the same non-solve.
#
# MODEL. Three trusses in a line with an Elastic material of E = 0, so the whole
# assembled A is zero and the failing pivot is guaranteed to be the first one.
# `algorithm Linear` does a single solve with no iteration that could mask it.
#
# ⚠ n >= 2 IS LOAD-BEARING. BandSPDLinLapackSolver.cpp:98 and
# ProfileSPDLinDirectSolver.cpp:150 carry 1x1 special-case guards that intercept
# a single-DOF singular system before the LAPACK call. An earlier 1-DOF version
# of this probe reported BandSPD as CLEAN for exactly that reason. A minimal
# repro can be too minimal.
#
# ProfileSPD is included as the control: it was always correct (hard -2 return),
# so it must be unaffected by the fix, and it proves the deck itself is sound.
#
# TWO ROUTES TO info == 1, and the distinction is the whole reason this survived
# so long. The defect fires ONLY at info == 1, because (-info)+1 is already
# negative for every info >= 2. Measured:
#   * free node numbered LAST  -> info = 3 -> old code returned -2 -> DETECTED
#   * free node numbered FIRST -> info = 1 -> old code returned  0 -> BUG
# So the "obvious" singular model (a forgotten free node) reproduces the bug only
# when its zero row lands on the first pivot. Case 1 below (E = 0 everywhere)
# pins info == 1 by construction; case 2 is the realistic shape that happens to
# hit it. Both are kept: the synthetic one can never stop reproducing, the
# realistic one is what a user would actually build.
#
#   PRE-FIX : BandGeneral/FullGeneral/BandSPD -> rc = 0, nodeDisp = 7.0 (the load)
#   POST-FIX: every solver -> rc != 0, nodeDisp = 0.0
#
# Run: OpenSees.exe lapack_singular_model.tcl   (exit status is the verdict)

set LOAD 7.0
set FAILURES 0

# model 1 -- synthetic: E = 0 everywhere, so ALL of A is zero and the failing
# pivot is the first one by construction. info == 1 guaranteed.
proc build_zeroE {} {
    global LOAD
    model basic -ndm 1 -ndf 1
    node 1 0.0 ; node 2 1.0 ; node 3 2.0 ; node 4 3.0
    fix 1 1
    uniaxialMaterial Elastic 1 0.0     ;# zero stiffness => A is identically zero
    element truss 1 1 2 1.0 1
    element truss 2 2 3 1.0 1
    element truss 3 3 4 1.0 1
    pattern Plain 1 Linear { load 4 $LOAD }
}

# model 2 -- realistic: node 1 is unfixed and carries no element, so its zero row
# is numbered FIRST under `numberer Plain` => info == 1. Measured: moving the
# free node to the end instead gives info == 3, which the OLD code detected
# correctly (-3+1 = -2). Position, not singularity, is what triggered the defect.
proc build_freenode {} {
    global LOAD
    model basic -ndm 1 -ndf 1
    node 1 0.0 ; node 2 1.0 ; node 3 2.0 ; node 4 3.0
    fix 4 1
    uniaxialMaterial Elastic 1 1000.0
    element truss 2 2 3 1.0 1
    element truss 3 3 4 1.0 1
    pattern Plain 1 Linear { load 2 $LOAD }
}

proc run_case {sysCmd builder} {
    global LOAD
    wipe
    $builder

    constraints Plain
    numberer Plain
    eval system $sysCmd
    test NormUnbalance 1.0e-8 10 0
    algorithm Linear
    integrator LoadControl 1.0
    analysis Static

    set rc [analyze 1]
    # The tell-tale of the bug is that the "solution" is the un-solved RHS, so
    # read the loaded DOF -- node 4 for model 1, node 2 for model 2.
    if {$builder eq "build_zeroE"} { set d [nodeDisp 4 1] } else { set d [nodeDisp 2 1] }
    return [list $rc $d]
}

puts "=== LAPACK singular-matrix regression: 2 models x 4 solvers ==="
foreach model {build_zeroE build_freenode} {
foreach sys {BandGeneral FullGeneral BandSPD ProfileSPD} {
    set res [run_case $sys $model]
    set rc  [lindex $res 0]
    set d   [lindex $res 1]

    set ok 1
    if {$rc == 0} { set ok 0 }
    # The tell-tale of the bug: the "solution" is the un-solved load vector.
    if {[expr {abs($d - $LOAD)}] < 1.0e-12} { set ok 0 }

    if {$ok} {
        puts [format "  PASS  %-15s system %-12s rc = %3d (nonzero), nodeDisp = %g" \
                  $model $sys $rc $d]
    } else {
        puts [format "  FAIL  %-15s system %-12s rc = %3d, nodeDisp = %g  <-- singular system reported SUCCESS" \
                  $model $sys $rc $d]
        incr FAILURES
    }
}
}

wipe
if {$FAILURES > 0} {
    puts "=== LAPACK singular regression: $FAILURES FAILURE(S) ==="
    exit 1
}
puts "=== LAPACK singular regression: all solvers reject the singular system ==="
exit 0
