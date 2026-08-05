# ADR-76 R4 acceptance test 4 -- "algorithm ModifiedNewton -initial -factoronce:
# both options take effect", in BOTH argument orders.
#
# Why a pre-yielded model. On a virgin model Ki == Kt, so `-factoronce` and
# `-initial -factoronce` latch the SAME matrix and the two spellings are
# indistinguishable. We therefore push the chain well past yield with a plain
# Newton phase (unprofiled), and only then start the profiled case. At that
# point Kt = b*Ki = 0.1*Ki, so the two spellings latch visibly different
# matrices and separate by iteration count.
#
# The four cases and what each isolates:
#
#   caseA  -initial                 formTangent == nSteps   (baseline: no latch)
#   caseB  -initial -factoronce     formTangent == 1, SLOW  (Ki latched)
#   caseC  -factoronce -initial     must equal caseB exactly (order-independence)
#   caseD  -factoronce              formTangent == 1, FAST  (Kt latched)
#
# Pre-fix, OPS_ModifiedNewton read exactly one option, so caseB collapsed onto
# caseA (`-factoronce` dropped) and caseC collapsed onto caseD (`-initial`
# dropped). Post-fix, B == C and both differ from A and from D.
#
# Run:  OpenSees.exe adr76_factoronce_model.tcl
# Then: python3.12 adr76_factoronce_check.py <dir>

set OUT [file dirname [info script]]
set H5  [file join $OUT adr76_factoronce.h5]
set TXT [file join $OUT adr76_factoronce_disp.txt]
file delete -force $H5
file delete -force $TXT

set N       10        ;# trusses in series
set FY      100.0
set E0      1000.0
set BHARD   0.1       ;# Kt/Ki past yield -- keeps initial-stiffness iteration
                       # convergent in a sane number of passes
set PRESTEP 15        ;# unprofiled Newton steps, dlambda 10 -> lambda 150 >> Fy
set NSTEPS  5         ;# profiled steps

proc build_model {} {
    global N FY E0 BHARD
    wipe
    model basic -ndm 2 -ndf 2

    for {set i 1} {$i <= [expr {$N + 1}]} {incr i} {
        node $i [expr {($i - 1) * 1.0}] 0.0
    }
    fix 1 1 1
    for {set i 2} {$i <= [expr {$N + 1}]} {incr i} { fix $i 0 1 }

    uniaxialMaterial Steel01 1 $FY $E0 $BHARD
    for {set i 1} {$i <= $N} {incr i} {
        element truss $i $i [expr {$i + 1}] 1.0 1
    }

    pattern Plain 1 Linear { load [expr {$N + 1}] 1.0 0.0 }
}

# Push past yield with full Newton, profiler OFF. Leaves the domain in the
# yielded state every case starts from.
proc preyield {} {
    global PRESTEP
    constraints Plain
    numberer RCM
    system BandGeneral
    test NormUnbalance 1.0e-10 100 0
    algorithm Newton
    integrator LoadControl 10.0
    analysis Static
    if {[analyze $PRESTEP] != 0} { error "pre-yield phase failed" }
}

# One profiled case. `args` is the option list handed to `algorithm ModifiedNewton`.
proc run_case {runid args} {
    global N NSTEPS H5 TXT
    build_model
    preyield

    wipeAnalysis
    constraints Plain
    numberer RCM
    system BandGeneral
    # Initial-stiffness iteration past yield contracts at (1 - Kt/Ki) = 0.9 per
    # pass, so it needs room; caseD (current tangent) will finish in a handful.
    test NormUnbalance 1.0e-8 5000 0
    eval algorithm ModifiedNewton $args
    integrator LoadControl 2.0
    analysis Static

    profiler start
    set ok [analyze $NSTEPS]
    profiler stop
    if {$ok != 0} { error "case $runid ($args) failed to converge" }

    profiler report $H5 -run $runid
    profiler reset

    # The converged answer must agree across all cases. That holds under this
    # deck's FORCE-residual test (NormUnbalance); it would NOT hold under
    # NormDispIncr/EnergyIncr, which accept on dU = K_f^-1 R. Emitted so the
    # checker can assert it.
    set ux [nodeDisp [expr {$N + 1}] 1]
    set fh [open $TXT a]
    puts $fh [format "%s %.17e" $runid $ux]
    close $fh

    puts [format "  %-6s algorithm ModifiedNewton %-24s ux(tip) = %.12e" \
              $runid $args $ux]
}

puts "=== ADR-76 R4 smoke: -initial / -factoronce option-loop ==="
run_case caseA -initial
run_case caseB -initial -factoronce
run_case caseC -factoronce -initial
run_case caseD -factoronce
# caseE/caseF: the spelling variants the adversarial review caught. `-factorOnce`
# (camelCase) is what Linear, ExpressNewton and both commands.cpp sites accept and
# what LEDGER_quirks writes, so it is the likelier carry-over spelling -- and it
# was STILL being silently dropped by the first cut of the R4 fix. `-Initial` is
# accepted by NewtonRaphson and was not by ModifiedNewton. Both must now behave
# exactly like caseB.
run_case caseE -initial -factorOnce
run_case caseF -Initial -factorOnce
puts "=== wrote $H5 ==="
wipe
