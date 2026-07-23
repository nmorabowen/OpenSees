# End-to-end replicated-MPI smoke test for the independent two-level CMS command.

wipe
model BasicBuilder -ndm 1 -ndf 1
set freeDOFs 16
for {set tag 1} {$tag <= $freeDOFs + 1} {incr tag} {
    node $tag [expr {double($tag - 1)}]
}
fix 1 1
uniaxialMaterial Elastic 1 1.0
for {set tag 1} {$tag <= $freeDOFs} {incr tag} {
    element truss $tag $tag [expr {$tag + 1}] 1.0 1
}
for {set tag 2} {$tag <= $freeDOFs + 1} {incr tag} {
    mass $tag 1.0
}

proc solveCMS {} {
    return [eigen -ladrunoCMS \
        -hierarchy logical -level1 2 -level2 2 \
        -modesL2 4 -modesL1 8 -tol 1.0e-8 -maxEnrich 2 \
        -denseMax 128 -verifyAssembly full \
        -verifyFullMaxBytes 1048576 4]
}

set first [solveCMS]
set second [solveCMS]
foreach label {first second} values [list $first $second] {
    if {[llength $values] != 4} {
        error "$label solve returned [llength $values] modes"
    }
    for {set index 0} {$index < 4} {incr index} {
        set mode [expr {$index + 1}]
        set reference [expr {2.0 - 2.0*cos((2*$mode - 1)*acos(-1.0)/(2*$freeDOFs + 1))}]
        set actual [lindex $values $index]
        set relativeError [expr {abs($actual - $reference)/max(abs($reference), 1.0e-30)}]
        if {$relativeError > 1.0e-8} {
            error "$label mode $mode relative eigenvalue error $relativeError"
        }
    }
}

if {[getPID] == 0} {
    puts "ladruno CMS OpenSeesMP two-level smoke passed"
}
wipe
