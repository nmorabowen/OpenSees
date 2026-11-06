# A later-constructed ManzariDafalias must not reset the class-wide stage (apeGmsh #2444).
#
# Laterally confined stdBrick column. Stage 0 settles cube 1 under a top load
# (elastic); the material is flipped to stage 1; cube 2 is activated on top
# (new top nodes ride on the settled interface, so it starts stress-free) and
# the load is increased. The interface settlement uz(5) must equal the same
# flip WITHOUT the activation (the new cube only carries the extra load down),
# and must be far from the elastic answer. Before the fix the activation's
# constructor reset mElastFlag = 0 and uz(5) came back elastic.
# Stock 3.8.0 on this deck: -3.4067686468e-3 (plastic), -2.8970896097e-3 (elastic).
# Terminal marker: "MD_STAGE_2444: PASS" or "MD_STAGE_2444: FAIL ...".

set nFail 0
set MD {125 0.05 0.8 1.25 0.712 0.019 0.934 0.7 100 0.01 7.05 0.968 1.1 0.704 3.5 4 600 1.42}
set XY {{0. 0.} {1. 0.} {1. 1.} {0. 1.}}

proc run {flip activate} {
    global MD XY
    wipe
    model BasicBuilder -ndm 3 -ndf 3
    for {set k 0} {$k < 2} {incr k} {
        for {set j 0} {$j < 4} {incr j} {
            set xy [lindex $XY $j]
            node [expr {4*$k + $j + 1}] [lindex $xy 0] [lindex $xy 1] [expr {double($k)}]
        }
    }
    nDMaterial ManzariDafalias 1 {*}$MD
    element stdBrick 1 1 2 3 4 5 6 7 8 1
    for {set n 1} {$n <= 8} {incr n} { fix $n 1 1 [expr {$n <= 4 ? 1 : 0}] }
    updateMaterialStage -material 1 -stage 0
    timeSeries Linear 1
    pattern Plain 1 1 { for {set n 5} {$n <= 8} {incr n} { load $n 0. 0. -25. } }
    constraints Plain
    numberer RCM
    system BandGeneral
    test NormDispIncr 1e-8 50 0
    algorithm Newton
    integrator LoadControl 0.25
    analysis Static
    if {[analyze 4] != 0} { return "FAILED_STAGE0" }
    loadConst -time 0.0
    if {$flip} { updateMaterialStage -material 1 -stage 1 }
    set top {5 6 7 8}
    if {$activate} {
        for {set j 0} {$j < 4} {incr j} {
            set xy [lindex $XY $j]
            node [expr {9 + $j}] [lindex $xy 0] [lindex $xy 1] 2.0
            fix [expr {9 + $j}] 1 1 0
            setNodeDisp [expr {9 + $j}] 3 [nodeDisp [expr {5 + $j}] 3] -commit
        }
        element stdBrick 2 5 6 7 8 9 10 11 12 1
        set top {9 10 11 12}
    }
    timeSeries Linear 2
    pattern Plain 2 2 { foreach n $top { load $n 0. 0. -25. } }
    integrator LoadControl 0.05
    if {[analyze 20] != 0} { return "FAILED_STAGE1" }
    return [nodeDisp 5 3]
}

set plasticNoAdd [run 1 0]
set plasticAdd [run 1 1]
set elastic [run 0 0]
puts "uz5: flip+activate=$plasticAdd flip=$plasticNoAdd elastic=$elastic"
if {![string is double -strict $plasticAdd] || ![string is double -strict $plasticNoAdd] ||
    ![string is double -strict $elastic]} {
    puts "FAIL a leg did not converge"
    incr nFail
} else {
    if {abs($plasticAdd - $plasticNoAdd) <= 1e-4 * abs($plasticNoAdd)} {
        puts "PASS activation keeps the plastic stage"
    } else {
        puts "FAIL activation after the flip changed uz5 (stage reset, #2444)"
        incr nFail
    }
    if {abs($plasticNoAdd - $elastic) > 0.1 * abs($elastic)} {
        puts "PASS plastic and elastic legs are distinguishable"
    } else {
        puts "FAIL plastic and elastic legs coincide"
        incr nFail
    }
}

if {$nFail == 0} {
    puts "MD_STAGE_2444: PASS"
} else {
    puts "MD_STAGE_2444: FAIL ($nFail checks)"
}
