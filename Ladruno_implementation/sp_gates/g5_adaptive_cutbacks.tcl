# ADR-80 gate G5 — ADAPTIVE march with CUTBACKS.
#
# WHY THIS EXISTS. The G1/G2/G3 harness (g123_predictor_gates.tcl) marches a
# FIXED number of steps that always succeed, so it can only ever measure
# iterations-within-successful-steps. But the Cerro Lindo x28.6 was dominated by
# **52 cutbacks** -- a failure mode a fixed march cannot exhibit at all. Using
# the fixed harness as S1's only acceptance gate would therefore be blind to the
# very thing that made the field case expensive. This deck closes that gap.
#
# Same physics and sizing as g123 (see that file): 1D stdBrick chain, nu = 0,
# converged sigma = 300 MPa < fy = 379.5 so NOTHING yields physically and the
# whole penalty is spurious. The only change is the MARCH.
#
# Caller-driven adaptive ladder, mimicking the fork's robust-solve driver:
#   try dlambda; on failure halve it and retry (a CUTBACK); on `grow`
#   consecutive successes, enlarge by `up` (capped). Report increments,
#   cutbacks, and total Newton iterations -- the three numbers the findings
#   report for Cerro Lindo.
#
#   OpenSees.exe g5_adaptive_cutbacks.tcl      -> ./g5_adaptive_cutbacks.json

if {![info exists ::G5_OUT]} { set ::G5_OUT "g5_adaptive_cutbacks.json" }

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

# Adaptive ladder. maxIter is deliberately TIGHT (8): a generous cap would let
# every step limp to convergence and no cutback would ever fire, which would
# reproduce the blindness this deck exists to remove.
proc march {N drive mat {maxIter 8}} {
    global L AX EMOD DELTA
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

    set lam    0.0
    set dl     0.1        ;# start at 10 % of the load factor
    set dlmin  1.0e-4
    set dlmax  0.2
    set inc    0
    set cuts   0
    set iters  0
    set good   0
    set guard  0
    while {$lam < 1.0 - 1.0e-12 && $guard < 4000} {
        incr guard
        if {$lam + $dl > 1.0} { set dl [expr 1.0 - $lam] }
        integrator LoadControl $dl
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
            if {$dl < $dlmin} { return [list $inc $cuts $iters -1] }
        }
    }
    return [list $inc $cuts $iters [nodeDisp [expr $top+1] 3]]
}

# SWEEP the iteration budget rather than pick one. Cutbacks are not an intrinsic
# property of the model -- they appear only when maxIter falls BETWEEN the clean
# cost (2 iters/step) and the penalized cost (6 iters/step). Reporting one
# hand-picked budget would be choosing the answer; reporting the sweep shows
# exactly where the penalty converts from "more iterations" into "failed steps",
# which is the regime the Cerro Lindo run was living in.
set rows {}
puts "COLS maxIter N drive mat increments cutbacks iters u_top"
foreach mi {3 4 5 6 8 12} {
    foreach N {20} {
        foreach drive {disp trac} {
            foreach mat {elastic j2} {
                if {[catch {march $N $drive $mat $mi} res]} {
                    puts "ROW mi=$mi N=$N $drive $mat ERROR"
                    lappend rows "  {\"maxIter\": $mi, \"N\": $N, \"drive\": \"$drive\", \"mat\": \"$mat\", \"status\": \"error\"}"
                    continue
                }
                lassign $res inc cuts iters utop
                set conv [expr {$utop < 0 ? "false" : "true"}]
                puts "ROW mi=[format %-3s $mi] N=$N [format %-5s $drive] [format %-7s $mat] inc=$inc cutbacks=$cuts iters=$iters utop=[format %.6f $utop]"
                lappend rows "  {\"maxIter\": $mi, \"N\": $N, \"drive\": \"$drive\", \"mat\": \"$mat\", \"increments\": $inc, \"cutbacks\": $cuts, \"iters\": $iters, \"u_top\": [format %.10g $utop], \"completed\": $conv}"
            }
        }
    }
}

set fh [open $::G5_OUT w]
puts $fh "{"
puts $fh "  \"gate\": \"ADR-80 G5 — adaptive march with cutbacks (closes the fixed-march blind spot)\","
puts $fh "  \"delta\": $DELTA, \"fy\": $FY, \"converged_sigma\": [expr $EMOD*$DELTA/$L],"
puts $fh "  \"maxIter\": 8,"
puts $fh "  \"rows\": \["
puts $fh [join $rows ",\n"]
puts $fh "  \]"
puts $fh "}"
close $fh
puts "WROTE $::G5_OUT"
