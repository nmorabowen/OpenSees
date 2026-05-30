# T1 Tcl twin of tests/test_forceBeamColumn_cantilever.py.
# Reads the SAME shared params file (D1). Appends PASSED/FAILED to results.out
# AND exits nonzero on failure so it gates CI even outside the suite (G1).
wipe

set here [file dirname [info script]]
source [file join $here _params.tcl]
read_params [file join $here .. params forceBeamColumn_cantilever.params] P
set L $P(L); set E $P(E); set A $P(A); set I $P(I); set Pload $P(P)

model basic -ndm 2 -ndf 3
node 1 0.0 0.0;  fix 1 1 1 1
node 2 $L  0.0
section Elastic 1 $E $A $I
beamIntegration Lobatto 1 1 3
geomTransf Linear 1
element forceBeamColumn 1 1 2 1 1
timeSeries Constant 1
pattern Plain 1 1 { load 2 0.0 $Pload 0.0 }
analysis Static
analyze 1

set tip   [nodeDisp 2 2]
set exact [expr {$Pload*$L*$L*$L/(3.0*$E*$I)}]
set rel   [expr {abs($tip-$exact)/abs($exact)}]

set results [open [file join $here .. .. EXAMPLES verification results.out] a]
if {$rel < 1e-9} {
    puts $results "PASSED : forceBeamColumn_cantilever tip (rel $rel)"
    close $results
    puts "PASSED forceBeamColumn_cantilever"
} else {
    puts $results "FAILED : forceBeamColumn_cantilever tip $tip vs $exact (rel $rel)"
    close $results
    puts "FAILED forceBeamColumn_cantilever tip $tip vs $exact"
    exit 2
}
