# ADR-80 S1 acceptance gate — LadrunoLoadControl -extrapolate.
#
# Same model/sizing as g123 and g5 (1D stdBrick chain, nu = 0, converged
# sigma = 300 MPa < fy = 379.5 so NOTHING yields physically and the whole
# penalty is spurious).
#
# Two things are measured, and the second one is a DESIGN QUESTION, not a
# formality:
#
#   IDIOM A "reissue"  -- the caller re-issues `integrator ... $dl` every step,
#                         which is what the G5 adaptive ladder and the fork's
#                         robust-solve driver do. Each re-issue CONSTRUCTS A NEW
#                         INTEGRATOR OBJECT, so any predictor state carried on
#                         that object is destroyed. If the predictor is dead
#                         under this idiom it is dead for the fork's real driver.
#   IDIOM B "persist"  -- the integrator is created ONCE and `analyze 1` is
#                         called repeatedly, so the object (and its predictor
#                         state) survives across steps.
#
# ADR-80 §6 anticipated "does the predictor interact badly with a caller-driven
# adaptive march?" -- this deck answers it with numbers.
#
#   OpenSees.exe s1_extrapolate_acceptance.tcl  -> ./s1_extrapolate_acceptance.json

if {![info exists ::S1_OUT]} { set ::S1_OUT "s1_extrapolate_acceptance.json" }

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
        nDMaterial LadrunoJ2 1 [expr $EMOD/3.0] [expr $EMOD/2.0] -iso voce $FY 0.0 0.0 $HISO
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
    return $top
}

# integ: "stock" | "ladruno"; frac only used when integ eq ladruno
proc integ_cmd {integ dl frac} {
    if {$integ eq "stock"} {
        integrator LoadControl $dl
    } else {
        integrator LadrunoLoadControl $dl 1 $dl $dl -extrapolate $frac
    }
}

# IDIOM A: re-issue the integrator every step (destroys per-object state)
proc march_reissue {N drive mat integ frac {maxIter 5}} {
    set top [setup $N $drive $mat $maxIter]
    set lam 0.0; set dl 0.1; set dlmin 1.0e-4; set dlmax 0.2
    set inc 0; set cuts 0; set iters 0; set good 0; set guard 0
    while {$lam < 1.0 - 1.0e-12 && $guard < 4000} {
        incr guard
        if {$lam + $dl > 1.0} { set dl [expr 1.0 - $lam] }
        integ_cmd $integ $dl $frac
        set rc [analyze 1]
        set iters [expr $iters + [testIter]]
        if {$rc == 0} {
            set lam [expr $lam + $dl]; incr inc; incr good
            if {$good >= 2} { set dl [expr min($dl*2.0,$dlmax)]; set good 0 }
        } else {
            incr cuts; set good 0; set dl [expr $dl/2.0]
            if {$dl < $dlmin} { return [list $inc $cuts $iters -1] }
        }
    }
    return [list $inc $cuts $iters [nodeDisp [expr $top+1] 3]]
}

# IDIOM B: create the integrator ONCE; fixed dl, repeated `analyze 1`
proc march_persist {N drive mat integ frac nsteps {maxIter 5}} {
    set top [setup $N $drive $mat $maxIter]
    integ_cmd $integ [expr 1.0/$nsteps] $frac
    set inc 0; set cuts 0; set iters 0
    for {set s 0} {$s < $nsteps} {incr s} {
        set rc [analyze 1]
        set iters [expr $iters + [testIter]]
        if {$rc == 0} { incr inc } else { incr cuts }
    }
    return [list $inc $cuts $iters [nodeDisp [expr $top+1] 3]]
}

# IDIOM C: the FIX. Same adaptive ladder as A, but the step is resized through
# `ladrunoLoadControl setDeltaLambda` instead of re-issuing `integrator ...`,
# so the object -- and its predictor state -- survives across steps AND cutbacks.
proc march_cmd {N drive mat frac {maxIter 5}} {
    set top [setup $N $drive $mat $maxIter]
    # wide min/max so the integrator's own clamp never interferes; setDeltaLambda
    # also pins numIncrLastStep = specNumIncrStep so the built-in dlambda
    # adaptation factor is exactly 1 and the caller stays in control.
    integrator LadrunoLoadControl 0.1 1 1.0e-8 1.0 -extrapolate $frac
    set lam 0.0; set dl 0.1; set dlmin 1.0e-4; set dlmax 0.2
    set inc 0; set cuts 0; set iters 0; set good 0; set guard 0
    while {$lam < 1.0 - 1.0e-12 && $guard < 4000} {
        incr guard
        if {$lam + $dl > 1.0} { set dl [expr 1.0 - $lam] }
        ladrunoLoadControl setDeltaLambda $dl
        set rc [analyze 1]
        set iters [expr $iters + [testIter]]
        if {$rc == 0} {
            set lam [expr $lam + $dl]; incr inc; incr good
            if {$good >= 2} { set dl [expr min($dl*2.0,$dlmax)]; set good 0 }
        } else {
            incr cuts; set good 0; set dl [expr $dl/2.0]
            if {$dl < $dlmin} { return [list $inc $cuts $iters -1] }
        }
    }
    return [list $inc $cuts $iters [nodeDisp [expr $top+1] 3]]
}

set rows {}
proc emit {idiom integ frac label} {
    global rows
    foreach mat {elastic j2} {
        if {$idiom eq "reissue"} {
            set res [march_reissue 20 disp $mat $integ $frac]
        } elseif {$idiom eq "cmd"} {
            set res [march_cmd 20 disp $mat $frac]
        } else {
            set res [march_persist 20 disp $mat $integ $frac 10 8]
        }
        lassign $res inc cuts iters utop
        puts "ROW [format %-8s $idiom] [format %-22s $label] [format %-7s $mat] inc=$inc cutbacks=$cuts iters=$iters utop=[format %.6f $utop]"
        lappend rows "  {\"idiom\": \"$idiom\", \"config\": \"$label\", \"mat\": \"$mat\", \"increments\": $inc, \"cutbacks\": $cuts, \"iters\": $iters, \"u_top\": [format %.10g $utop]}"
    }
}

puts "COLS idiom config mat increments cutbacks iters u_top"
# --- IDIOM A (the fork's real driver shape): does the predictor survive? ------
emit reissue stock   0.0 "stock LoadControl"
emit reissue ladruno 0.0 "Ladruno frac=0"
emit reissue ladruno 1.0 "Ladruno frac=1"
# --- IDIOM B (integrator persists): the predictor's best case -----------------
emit persist stock   0.0 "stock LoadControl"
emit persist ladruno 0.0 "Ladruno frac=0"
emit persist ladruno 0.5 "Ladruno frac=0.5"
emit persist ladruno 1.0 "Ladruno frac=1"
# --- IDIOM C (the fix): adaptive ladder via `ladrunoLoadControl setDeltaLambda` -
emit cmd ladruno 0.0 "Ladruno frac=0"
emit cmd ladruno 1.0 "Ladruno frac=1"

set fh [open $::S1_OUT w]
puts $fh "{"
puts $fh "  \"gate\": \"ADR-80 S1 acceptance — LadrunoLoadControl -extrapolate\","
puts $fh "  \"rows\": \["
puts $fh [join $rows ",\n"]
puts $fh "  \]"
puts $fh "}"
close $fh
puts "WROTE $::S1_OUT"
