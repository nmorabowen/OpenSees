# ADR-80 P3 acceptance gate — LadrunoLoadControl -tangentPredictor
# (candidate C: the Kratos-style `b -= K*du_D` route, D6).
#
# SAME model, sizing and marches as g123 / g5 / s1 so the numbers are directly
# comparable to the S1 verdict in 80c: 1D stdBrick chain, nu = 0, converged
# sigma = 300 MPa < fy = 379.5, so NOTHING yields physically and the entire
# penalty is spurious.
#
# WHAT S1 MEASURED (80c), and therefore what this deck must beat:
#
#   march      variant                 inc   cutbacks   iters
#   fixed      elastic control          10       0        20
#   fixed      stock / -extrapolate 0   10       0        60
#   fixed      -extrapolate 1.0         10       0        35
#   adaptive   elastic control           6       0        12
#   adaptive   stock / -extrapolate 0   43      23       224
#   adaptive   -extrapolate 1.0 (persist) 42     23       201   <- gate FAILED
#
# The adaptive row is the gate that matters: linear extrapolation moved
# iterations by -10 % and moved CUTBACKS BY NOTHING. -tangentPredictor exists
# because it eliminates the overstrain instead of reducing it.
#
# It is also STATELESS, which is a claim this deck tests directly: the adaptive
# ladder below RE-ISSUES `integrator LadrunoLoadControl ...` every step, which
# constructs a new object. That idiom made -extrapolate completely inert (80c
# §2, armed = 0). If the tangent route works here, it works under the fork's
# real robust-solve driver.
#
#   OpenSees.exe p3_tangent_predictor_acceptance.tcl
#       -> ./p3_tangent_predictor_acceptance.json

if {![info exists ::P3_OUT]} { set ::P3_OUT "p3_tangent_predictor_acceptance.json" }

set L      100.0
set AX      10.0
set EMOD 200000.0
set FY     379.5
set HISO  2000.0
set DELTA    0.15

proc build {N mat} {
    global L AX EMOD FY HISO
    wipe
    model basic -ndm 3 -ndf 3
    set h [expr $L/double($N)]
    for {set k 0} {$k <= $N} {incr k} {
        set z [expr $k*$h]; set b [expr 4*$k]
        node [expr $b+1] 0.0 0.0 $z
        node [expr $b+2] $AX 0.0 $z
        node [expr $b+3] $AX $AX $z
        node [expr $b+4] 0.0 $AX $z
    }
    if {$mat eq "j2"} {
        nDMaterial LadrunoJ2 1 [expr $EMOD/3.0] [expr $EMOD/2.0] \
            -iso voce $FY 0.0 0.0 $HISO
    } else {
        nDMaterial ElasticIsotropic 1 $EMOD 0.0
    }
    for {set k 0} {$k < $N} {incr k} {
        set b [expr 4*$k]
        element stdBrick [expr $k+1] [expr $b+1] [expr $b+2] [expr $b+3] [expr $b+4] \
                                     [expr $b+5] [expr $b+6] [expr $b+7] [expr $b+8] 1
    }
    for {set k 0} {$k <= $N} {incr k} {
        set b [expr 4*$k]
        for {set j 1} {$j <= 4} {incr j} {
            if {$k == 0} { fix [expr $b+$j] 1 1 1 } else { fix [expr $b+$j] 1 1 0 }
        }
    }
}

proc setup {N drive mat maxIter} {
    global EMOD AX DELTA L
    build $N $mat
    set top [expr 4*$N]
    timeSeries Linear 1
    pattern Plain 1 1 {
        if {$drive eq "disp"} {
            for {set j 1} {$j <= 4} {incr j} { sp [expr $top+$j] 3 $DELTA }
        } else {
            set F [expr $EMOD*$AX*$AX*$DELTA/$L/4.0]
            for {set j 1} {$j <= 4} {incr j} { load [expr $top+$j] 0.0 0.0 $F }
        }
    }
    constraints Transformation
    numberer RCM
    system BandGeneral
    test NormDispIncr 1.0e-10 $maxIter 0
    algorithm KrylovNewton
    analysis Static
    return $top
}

# THE DRIVEN NODE IS NOT A DISCRIMINATOR. Reporting u at the sp'd face would
# report the prescribed value back to itself, and would have passed even for the
# S2 silent-wrong-answer (only the boundary layer moved). Probe the MIDPOINT,
# whose exact answer under the linear elastic chain is delta/2 = 0.075 and which
# is 0.0 exactly in the failure mode.
proc umid {N} { return [nodeDisp [expr 4*($N/2)+1] 3] }

# `variant` selects the integrator command issued for a step of size $dl:
#   stock      -> integrator LoadControl $dl
#   ext0       -> integrator LadrunoLoadControl $dl                 (safety row)
#   ext1       -> integrator LadrunoLoadControl $dl -extrapolate 1.0
#   tangent    -> integrator LadrunoLoadControl $dl -tangentPredictor
proc issue {variant dl} {
    switch -- $variant {
        stock   { integrator LoadControl $dl }
        ext0    { integrator LadrunoLoadControl $dl }
        ext1    { integrator LadrunoLoadControl $dl -extrapolate 1.0 }
        tangent { integrator LadrunoLoadControl $dl -tangentPredictor }
        default { error "unknown variant $variant" }
    }
}

# ---------------------------------------------------------------- fixed march
# nSteps equal increments, generous iteration budget so no step can fail. This
# isolates iterations-within-successful-steps; it CANNOT see cutbacks, which is
# why the adaptive march below exists (see g5's header).
#
# TWO IDIOMS, because they are not equivalent and 80c is the reason:
#   reissue -- `integrator ...` re-issued every step (the fork's robust-solve
#              driver, and g5's ladder). CONSTRUCTS A NEW OBJECT each step, so
#              any predictor state on it dies. This is what made -extrapolate
#              completely inert (80c 2, armed = 0).
#   persist -- constructed once, `analyze 1` repeatedly. The only idiom under
#              which -extrapolate can fire at all.
# -tangentPredictor is stateless and must score the same under both.
proc fixed_march {N drive mat variant {idiom reissue} {nSteps 10} {maxIter 40}} {
    set top [setup $N $drive $mat $maxIter]
    set dl [expr 1.0/double($nSteps)]
    set iters 0
    set inc 0
    if {$idiom eq "persist"} { issue $variant $dl }
    for {set s 0} {$s < $nSteps} {incr s} {
        if {$idiom eq "reissue"} { issue $variant $dl }
        if {[analyze 1] != 0} { return [list $inc 0 $iters -1 -1] }
        set iters [expr $iters + [testIter]]
        incr inc
    }
    return [list $inc 0 $iters [nodeDisp [expr $top+1] 3] [umid $N]]
}

# ------------------------------------------------------------- adaptive march
# Caller-driven ladder identical to g5: halve on failure (a CUTBACK), double
# after 2 consecutive successes, cap at dlmax. maxIter is deliberately TIGHT and
# SWEPT below, so the penalty converts into failed steps rather than merely
# longer ones; 5 is the budget S1 measured 43/23/224 at.
# NOTE the deliberate RE-ISSUE of the integrator every step (see header).
proc adaptive_march {N drive mat variant {maxIter 5}} {
    set top [setup $N $drive $mat $maxIter]
    set lam   0.0
    set dl    0.1
    set dlmin 1.0e-4
    set dlmax 0.2
    set inc   0
    set cuts  0
    set iters 0
    set good  0
    set guard 0
    while {$lam < 1.0 - 1.0e-12 && $guard < 4000} {
        incr guard
        if {$lam + $dl > 1.0} { set dl [expr 1.0 - $lam] }
        issue $variant $dl
        set rc [analyze 1]
        set iters [expr $iters + [testIter]]
        if {$rc == 0} {
            set lam [expr $lam + $dl]
            incr inc
            incr good
            if {$good >= 2} { set dl [expr min($dl*2.0, $dlmax)]; set good 0 }
        } else {
            incr cuts
            set good 0
            set dl [expr $dl/2.0]
            if {$dl < $dlmin} { return [list $inc $cuts $iters -1 -1] }
        }
    }
    return [list $inc $cuts $iters [nodeDisp [expr $top+1] 3] [umid $N]]
}

set N 20
set rows {}
set SPECS {
    {disp elastic stock}
    {disp j2      stock}
    {disp j2      ext0}
    {disp j2      ext1}
    {disp j2      tangent}
    {disp elastic tangent}
    {trac j2      stock}
    {trac j2      tangent}
}

puts "COLS march idiom/maxIter drive mat variant increments cutbacks iters u_top u_mid"

# --- fixed march, both caller idioms ------------------------------------------
foreach idiom {reissue persist} {
    foreach spec $SPECS {
        lassign $spec drive mat variant
        if {[catch {fixed_march $N $drive $mat $variant $idiom} res]} {
            puts "ROW fixed $idiom $drive $mat $variant ERROR: $res"
            lappend rows "  {\"march\": \"fixed\", \"idiom\": \"$idiom\", \"drive\": \"$drive\", \"mat\": \"$mat\", \"variant\": \"$variant\", \"status\": \"error\"}"
            continue
        }
        lassign $res inc cuts iters utop umid
        set conv [expr {$utop < 0 ? "false" : "true"}]
        puts "ROW fixed    [format %-8s $idiom] [format %-5s $drive] [format %-7s $mat] [format %-8s $variant] inc=[format %-3s $inc] cutbacks=[format %-3s $cuts] iters=[format %-4s $iters] utop=[format %.6f $utop] umid=[format %.6f $umid]"
        lappend rows "  {\"march\": \"fixed\", \"idiom\": \"$idiom\", \"drive\": \"$drive\", \"mat\": \"$mat\", \"variant\": \"$variant\", \"increments\": $inc, \"cutbacks\": $cuts, \"iters\": $iters, \"u_top\": [format %.10g $utop], \"u_mid\": [format %.10g $umid], \"completed\": $conv}"
    }
}

# --- adaptive march, SWEEPING the iteration budget ----------------------------
# Cutbacks are not intrinsic to the model: they appear only when maxIter falls
# BETWEEN the clean cost and the penalized cost. g5 established that, and
# reporting a single hand-picked budget would be choosing the answer. maxIter 5
# is the budget S1 measured 43/23/224 at (s1_extrapolate_acceptance.tcl:96).
foreach mi {4 5 6 8} {
    foreach spec $SPECS {
        lassign $spec drive mat variant
        if {[catch {adaptive_march $N $drive $mat $variant $mi} res]} {
            puts "ROW adaptive mi=$mi $drive $mat $variant ERROR: $res"
            lappend rows "  {\"march\": \"adaptive\", \"maxIter\": $mi, \"drive\": \"$drive\", \"mat\": \"$mat\", \"variant\": \"$variant\", \"status\": \"error\"}"
            continue
        }
        lassign $res inc cuts iters utop umid
        set conv [expr {$utop < 0 ? "false" : "true"}]
        puts "ROW adaptive mi=[format %-6s $mi] [format %-5s $drive] [format %-7s $mat] [format %-8s $variant] inc=[format %-3s $inc] cutbacks=[format %-3s $cuts] iters=[format %-4s $iters] utop=[format %.6f $utop] umid=[format %.6f $umid]"
        lappend rows "  {\"march\": \"adaptive\", \"maxIter\": $mi, \"drive\": \"$drive\", \"mat\": \"$mat\", \"variant\": \"$variant\", \"increments\": $inc, \"cutbacks\": $cuts, \"iters\": $iters, \"u_top\": [format %.10g $utop], \"u_mid\": [format %.10g $umid], \"completed\": $conv}"
    }
}

set fh [open $::P3_OUT w]
puts $fh "{"
puts $fh "  \"gate\": \"ADR-80 P3 — LadrunoLoadControl -tangentPredictor (candidate C, b -= K*du_D)\","
puts $fh "  \"N\": $N, \"delta\": $DELTA, \"fy\": $FY, \"converged_sigma\": [expr $EMOD*$DELTA/$L],"
puts $fh "  \"fixed_march\": {\"steps\": 10, \"maxIter\": 40, \"idioms\": [\"reissue\", \"persist\"]},"
puts $fh "  \"adaptive_march\": {\"maxIter_sweep\": [4, 5, 6, 8], \"idiom\": \"re-issue the integrator every step (kills -extrapolate; the tangent route is stateless)\"},"
puts $fh "  \"rows\": \["
puts $fh [join $rows ",\n"]
puts $fh "  \]"
puts $fh "}"
close $fh
puts "WROTE $::P3_OUT"
