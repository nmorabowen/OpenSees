# ADR-80 gate G4 — does `constraints Auto` + a non-homogeneous `sp` commit a
# wrong answer? Runs the same A/B on any binary (pre- or post-fix) and writes a
# JSON artifact, so the measurement stays citable.
#
#   OpenSees.exe g4_auto_handler.tcl            -> ./g4_auto_handler.json
#   set ::G4_OUT <path> before sourcing to redirect.
#
# Reference handler is Transformation. The gate is "Auto must match it".
# Pre-fix signature: Auto + NormDispIncr/EnergyIncr -> interior 0.0 in 1 iteration.

if {![info exists ::G4_OUT]} { set ::G4_OUT "g4_auto_handler.json" }

set E 200000.0
set A 100.0
set L 1000.0
set DELTA 3.0
set EXPECTED [expr 0.5*$DELTA]

proc run_truss {handler testtype} {
    global E A L DELTA
    wipe
    model basic -ndm 2 -ndf 2
    node 1 0.0 0.0
    node 2 $L 0.0
    node 3 [expr 2.0*$L] 0.0
    fix 1 1 1
    fix 2 0 1
    fix 3 0 1
    uniaxialMaterial Elastic 1 $E
    element Truss 1 1 2 $A 1
    element Truss 2 2 3 $A 1
    timeSeries Linear 1
    pattern Plain 1 1 { sp 3 1 $DELTA }
    constraints $handler
    numberer Plain
    system BandGeneral
    test $testtype 1e-10 30 0
    algorithm Newton
    integrator LoadControl 1.0
    analysis Static
    set ok [analyze 1]
    return [list $ok [nodeDisp 2 1] [testIter]]
}

# ux,uy fixed everywhere => a 1D chain of two identical brick springs;
# the mid plane must reach DELTA/2 by symmetry.
proc run_brick {handler testtype} {
    global DELTA
    wipe
    model basic -ndm 3 -ndf 3
    set n 0
    foreach z {0.0 1.0 2.0} {
        foreach {x y} {0.0 0.0 1.0 0.0 1.0 1.0 0.0 1.0} {
            incr n
            node $n $x $y $z
        }
    }
    nDMaterial ElasticIsotropic 1 200000.0 0.25
    element stdBrick 1 1 2 3 4 5 6 7 8 1
    element stdBrick 2 5 6 7 8 9 10 11 12 1
    for {set i 1} {$i <= 12} {incr i} {
        if {$i <= 4} { fix $i 1 1 1 } else { fix $i 1 1 0 }
    }
    timeSeries Linear 1
    pattern Plain 1 1 {
        for {set i 9} {$i <= 12} {incr i} { sp $i 3 $DELTA }
    }
    constraints $handler
    numberer Plain
    system BandGeneral
    test $testtype 1e-10 30 0
    algorithm Newton
    integrator LoadControl 1.0
    analysis Static
    set ok [analyze 1]
    return [list $ok [nodeDisp 5 3] [testIter]]
}

set rows {}
foreach {label proc_} {truss run_truss brick run_brick} {
    foreach testtype {NormDispIncr EnergyIncr NormUnbalance} {
        foreach handler {Transformation Auto Plain} {
            if {[catch {$proc_ $handler $testtype} res]} {
                lappend rows "  {\"model\": \"$label\", \"handler\": \"$handler\", \"test\": \"$testtype\", \"status\": \"error\"}"
            } else {
                lassign $res ok u it
                set match [expr {abs($u - $EXPECTED) < 1.0e-9 ? "true" : "false"}]
                lappend rows "  {\"model\": \"$label\", \"handler\": \"$handler\", \"test\": \"$testtype\", \"analyze_rc\": $ok, \"interior\": [format %.12g $u], \"iters\": $it, \"matches_expected\": $match}"
                puts "ROW $label $handler $testtype rc=$ok u=[format %.6f $u] iters=$it ok=$match"
            }
        }
    }
}

set fh [open $::G4_OUT w]
puts $fh "{"
puts $fh "  \"gate\": \"ADR-80 G4 — AutoConstraintHandler missing updateElement()\","
puts $fh "  \"expected_interior\": $EXPECTED,"
puts $fh "  \"delta\": $DELTA,"
puts $fh "  \"reference_handler\": \"Transformation\","
puts $fh "  \"rows\": \["
puts $fh [join $rows ",\n"]
puts $fh "  \]"
puts $fh "}"
close $fh
puts "WROTE $::G4_OUT"
